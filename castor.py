#!/usr/bin/python3

"""
Main module of Castor error assessment and correction. Gathers parameters
and runs all modules.

- TODO: split all functions into proper modules
"""


import sys
import argparse
from operator import itemgetter
from collections import Counter
from castor_class import *
from itertools import groupby


def _parse_arguments():
    """ Determine Castor parameters

    Returns
    -------
    args
        Software parameters

    """
    parser = argparse.ArgumentParser(
            usage=("ref_correct.py [OPTIONS] <draft_genome.fa> "
                   "<ref_genome.pileup> <mapped_reads.pileup>"),
            description=("Correct a draft genome using reference genomes "
                         "and reads mapped to the draft. version 0.2.8"))
    parser.add_argument(
            "draft",
            nargs="?",
            type=str,
            help="Draft genome fastA input file")
    parser.add_argument(
            "ref",
            nargs="?",
            type=str,
            help=("Reference genomes mapped to the draft genome in samtools "
                  "mpileup format file"))
    parser.add_argument(
            "reads",
            nargs="+",
            type=str,
            help=("Mpileup alignment information used to adjust found "
                  "errors. Recommended to use read information for "
                  "adjustment. More than one mpileup file can be used for "
                  "multiple adjustment attempts. One caveat: substitution "
                  "adjustments only use the last alignment file so "
                  "substitution data must be contained in that file "
                  "otherwise all substitution would be considered "
                  "unsupported and deleted. Ex. Fastq reads mapped to the "
                  "draft genome in samtools mpileup format file"))
    parser.add_argument(
            "--low", "-l",
            nargs="?",
            type=str,
            default="",
            help=("Optional: File path to mpileup file containing greater "
                  "and more distant reference genome used to correct for low"
                  " depth regions. --close mpileup set will be used to "
                  "adjust threshold. Note: imbalanced data will bias "
                  "correction towards the phylogroup containing the greatest"
                  " number of genomes in the dataset."))
    parser.add_argument(
            "--sub", "-s",
            nargs="?",
            type=float,
            default=1.0,
            help=("Threshold to correct substitutions. 0 to 1. Use >1 to "
                  "turn off substitution correction. [Default = 1.0]"))
    parser.add_argument(
            "--indel", "-i",
            nargs="?",
            type=float,
            default=0.99,
            help=("Threshold to correct indel errors. 0 to 1. Use > 1 to "
                  "turn off indel correction. [Default = 0.99]"))
    parser.add_argument(
            "--adjust", "-a",
            nargs="?",
            type=float,
            default=0.25,
            help=("Threshold to used to consider a candidate for error "
                  "adjustment. Use >2 to turn off error adjustments "
                  "[Default = 0.25]"))
    parser.add_argument(
            "--low_threshold",
            nargs="?",
            type=float,
            default=0.99,
            help=("Threshold used for optional low depth mpileup file. "
                  "[Default = 0.99]"))
    parser.add_argument(
            "--depth", "-d",
            nargs="?",
            type=int,
            default=5,
            help=("Minimum depth required to consider the mpileup "
                  "information. [Default = 5]"))
    parser.add_argument(
            "--passes", "-p",
            nargs="?",
            type=int,
            default=4,
            help=("Number of passes through error rates to determine "
                  "putative errors [Default = 4]"))
    parser.add_argument(
            "--out", "-o",
            nargs="?",
            type=str,
            default="out",
            help="Prefix for output files [Default = 'out']")
    parser.add_argument(
            "--extraInfo", "-e",
            action="store_true",
            help=("Print out extra information. Includes homopolymer "
                  "information at error site, regions skipped because of "
                  "low depth, initial error calls, and errors called at "
                  "regions of low depth."))
    parser.add_argument(
            "--module", "-m",
            nargs="?",
            type=str,
            default="All",
            choices=["All", "Detect", "Adjust", "Correct", "Errors-only"],
            help=("Temporary: Specify a particular module to run. "
                  "Options: All, Detect, Adjust, Correct, Errors-only."))
    parser.add_argument(
            "--verbal", "-v",
            action="store_true",
            help="Turn on informational/debugging text")
    parser.add_argument(
            "--errors",
            nargs="?",
            type=str,
            default="",
            help=("Bypass calculation of errors and correct with input "
                  "error file. Enter a previous error output file. (.err)"))

    return parser.parse_args()


def _read_draft(draft):
    """ Reads in a fastA file and returns a dictionary of the sequences.

    Taken from biostar (https://www.biostars.org/p/710/)
    Credit: brentp

    Parameters
    ----------
    draft : str
        A string of path leading to draft assembly to be corrected

    Returns
    -------
    Dict
        header is the key and sequence is the value
    """

    draft = open(draft)
    draft_iter = (x[1] for x in groupby(draft, lambda line: line[0] == ">"))

    d_gen = {}
    for header in draft_iter:
        header = header.__next__()[1:].strip()
        seq = "".join(s.strip() for s in draft_iter.__next__())
        print("Reading contig: {}".format(header))
        d_gen[header] = seq

    draft.close()

    return d_gen


def _load_prev_error(err_file):
    """ Loads a pre-computed error output in memory.

    Parameters
    ----------
    err_file : str
        Error file path as a string. Can be absolute or relative to program
        directory

    Returns
    -------
    list
        contains the errors formatted as a list

    """
    err_loaded = []
    with open(err_file) as f:
        for line in f:
            curr_line = line.split("\t")
            err_can = curr_line[0]
            contig = curr_line[1]
            pos = int(curr_line[2])
            depth = int(curr_line[3])
            err_len = int(curr_line[4])
            ref_nt = curr_line[5]
            replace_nt = curr_line[6]
            err_depth = curr_line[7]
            try:
                err_index = int(curr_line[8])
                err_loaded.append([err_can, contig, pos, depth, err_len, ref_nt,
                                   replace_nt, err_depth, err_index])
            except IndexError:
                err_loaded.append([err_can, contig, pos, depth, err_len, ref_nt,
                                   replace_nt, err_depth])

    return err_loaded


def _load_full_mpile_info(mpile_file):
    """ Load alignment files into memory

    Parameters
    ----------
    mpile_file : str
        Path pointing to the alignment file

    Returns
    -------
    list
        alignment with each line as a separate element. Empty list if there
        is no alignment file

    """
    mpile = []

    try:
        with open(mpile_file, "r") as f:
            mpile = f.readlines()
    except FileNotFoundError:
        pass

    return mpile


def _print_info(out_file, write_type, info):
    """ Outputs data structures passed as a tab-delimited file

    Parameters
    ----------
    out_file : str
        Output filename as a string
    write_type : char
        Output type, either as write or append
    info : dict or list
        Dictionary or list to output

    Returns
    -------
    None

    """
    with open(out_file, write_type) as out:
        if len(info) < 1:
            return

        # Assumes dictionary values are lists
        if type(info).__name__ == "dict":
                for key, value in info.items():
                    out.write("%s\n" % "\t".join(map(str, value)))
        elif type(info).__name__ == "list":
            if type(info[0]).__name__ == "list":
                for value in info:
                    out.write("%s\n" % "\t".join(map(str, value)))
            else:
                for value in info:
                    out.write(str(value) + "\n")


def _print_erates(out_file, write_type, err_rates, start, end):
    """ Outputs error rates for diagnosis within specified region.

    Parameters
    ----------
    out_file : str
        Output filename as a string
    write_type : char
        Output type, either write or append
    err_rates : list
        Calculated error rates from alignment files as Error classes
    start : int
        Index 1 of error rate list. Start of region
    end : int
        Index 2 of error rate list. End of region

    Returns
    -------
    None

    """

    with open(out_file, write_type) as out:
        if len(err_rates) < 1:
            return

        out.write("%s\n" % "\t".join(["contig", "pos", "ref", "depth", "ferate",
                                      "sub_ferate", "del_ferate", "ins_ferate",
                                      "sub_serate", "del_serate1", "del_serate2",
                                      "ins_serate1", "ins_serate2"]))
        try:
            for pos in range(start, end):
                line = [err_rates[pos].contig,
                        err_rates[pos].pos,
                        err_rates[pos].ref,
                        err_rates[pos].depth,
                        "{0:.3f}".format(err_rates[pos].erate),
                        "{0:.3f}".format(err_rates[pos].fullSub),
                        "{0:.3f}".format(err_rates[pos].fullDel),
                        "{0:.3f}".format(err_rates[pos].fullIns)]

                # substitutions
                try:
                    line += ["{0:.3f}".format(err_rates[pos].sub.erate)]
                except AttributeError:
                    line += ["{0:.3f}".format(0)]

                # Deletions
                try:
                    line += ["{0:.3f}".format(err_rates[pos].del1.erate)]
                except AttributeError:
                    line += ["{0:.3f}".format(0)]
                try:
                    line += ["{0:.3f}".format(err_rates[pos].del2.erate)]
                except AttributeError:
                    line += ["{0:.3f}".format(0)]

                # Insertions
                try:
                    line += ["{0:.3f}".format(err_rates[pos].ins1.erate)]
                except AttributeError:
                    line += ["{0:.3f}".format(0)]
                try:
                    line += ["{0:.3f}".format(err_rates[pos].ins2.erate)]
                except AttributeError:
                    line += ["{0:.3f}".format(0)]

                out.write("%s\n" % "\t".join(map(str, line)))
        except IndexError:
            print(("WARNING: error rates were not printed due to incorrect "
                   "data structure passed."), file=sys.stderr)


def check_indel_len(align, start):
    """ Check if an indel length is greater than 9 nucleotides

    Determines the indel length. Parses characters in the alignment
    sequentially until next character is not a number. If greater than 9
    then indel length takes two characters rather than one; a length >99
    will take three characters. The function also spaces the parsing index
    accordingly to the length.

    Parameters
    ----------
    align : str
        Alignment line
    start : int
        Character index to start parsing
    Returns
    -------
    int
        length of indel
    int
        position of the next read (index) in the alignment line

    """
    index = start + 1
    while index < len(align):
        try:
            int(align[index])
        except ValueError:
            indel_len = "".join(align[start:index])
            return int(indel_len), index
        index += 1


def get_most_common_nt(lst, ref_nt="-", next_nt="-"):
    """Determines the best common nucleotide for correction.

    Nucleotides are ranked based on both frequency and matching to flanking
    nucleotides. Ranking is mostly if there are multiple possibilities. A
    rank of 0 is the best with it being closest to the reference nucleotide
    while a rank of 3 indicates the nucleotide is "random". One issue is
    that there is no tiebreaker. Most common nucleotide is chosen by random
    during ties.

    Parameters
    ----------
    lst : list
        A list of all nucleotides suggested as a possible correction
    ref_nt : char
        The reference nucleotide at the position
    next_nt : char
        The downstream reference nucleotide

    Returns
    -------
    char
        the most common nucleotide in the list

    """
    if len(lst) < 1:
        return 0

    # list the most common element based on frequency
    poss_nt = Counter(lst).most_common()
    max_nt = 0
    best_rank = 99

    # determine best by comparing rank of the most frequent nucleotides
    if ref_nt == "-":
        return poss_nt[0][0]
    else:
        best_nt = ""
        for nt in Counter(lst).most_common():
            # set max count
            if nt[1] < max_nt:
                break
            else:
                max_nt = nt[1]

            # multiple max
            if nt[0][0] == ref_nt:
                if nt[0][-1] == ref_nt and 0 < best_rank:
                    best_nt = nt[0]
                    best_rank = 0
                elif nt[0][-1] == next_nt and 1 < best_rank:
                    best_nt = nt[0]
                    best_rank = 1
                else:
                    if 2 < best_rank:
                        best_nt = nt[0]
                        best_rank = 2
            elif nt[0][-1] == next_nt and 2 < best_rank:
                # only match to next 
                best_nt = nt[0]
                best_rank = 2
            elif 3 < best_rank:
                # everything else
                best_nt = nt[0]
                best_rank = 3
            else:
                pass

        return best_nt


def pull_nlength_indels_info(indel_list, n, alignment):
    """Subsets out string element with a length of n from a list

    Pulls out indel corrections of the appropriate length for error
    determination. Original list is not changed in case of multiple subsets

    Parameters
    ----------
    indel_list : list
        List of possible indels
    n : int
        Indel length to consider
    alignment : int
        An int that stores respective alignment positions of the indels in
        the list

    Returns
    -------
    list
        indel list subset
    int
        alignment positions of the elements in the subset

    """
    # initialization
    sub_align = 1
    subset_nt = []

    # alignment check
    if (bin(alignment).count("1") - 1) != len(indel_list):
        print(("BUG: Number of indels does not match bit alignments in "
               "function pull_nlength_indels_info"), file=sys.stderr)

    align_counter = alignment.bit_length() - 1

    # subset
    for index in range(len(indel_list)):
        for offset in range(align_counter, 0, -1):
            mask = 1 << (offset - 1)
            sub_align <<= 1
            if ((alignment & mask) >> (offset - 1)) == 1:
                alignment ^= mask
                align_counter = offset - 1
                break

        err_length = len(indel_list[index])
        if err_length == n:
            subset_nt.append(indel_list[index])
            sub_align ^= 1

    # alignment adjustment for depth
    if alignment.bit_length() != sub_align.bit_length():
        sub_align <<= (alignment.bit_length() - sub_align.bit_length())

    return subset_nt, sub_align


def get_error_rate(err_type, pos_erate, err_info=None, indel_alignment=0):
    """Calculate the position's error rate of a specific type

    Parameters
    ----------
    err_type : str
        The error type to consider for calculation
    pos_erate : Error class
        The class to store error rate information
    err_info : list
        A list of containing the errors and suggested corrections
    indel_alignment : int
        Alignment information for the respective error

    Returns
    -------
    None
        results are stored in pos_erate

    """
    depth = pos_erate.depth

    # full error rate
    if err_info is None:
        f_erate = 1
    else:
        f_erate = (len(err_info)/depth)

    if err_type == "total":
        pos_erate.erate = 1 - f_erate
    elif err_type == "del":
        # adjusts bit length depending on depth
        bit_diff = pos_erate.readEnd.bit_length() - indel_alignment.bit_length()
        if bit_diff > 0:
            indel_alignment <<= bit_diff
        elif bit_diff < 0:
            print(("BUG: Difference in alignment bit length in function "
                   "get_error_rate"), file=sys.stderr)
        else:
            pass

        # deletion error rate
        pos_erate.fullDel = f_erate
        sub_del, sub_align = pull_nlength_indels_info(err_info, 1,
                                                      indel_alignment)
        if sub_del:
            pos_erate.del1 = ErrType(sub_del[0], len(sub_del)/depth,
                                     sub_align)
        sub_del, sub_align = pull_nlength_indels_info(err_info, 2,
                                                      indel_alignment)
        if sub_del:
            pos_erate.del2 = ErrType(sub_del[0], len(sub_del)/depth,
                                     sub_align)
    elif err_type == "ins":
        # insertion error rate
        pos_erate.fullIns = f_erate
        sub_ins, sub_align = pull_nlength_indels_info(err_info, 1,
                                                      indel_alignment)
        if sub_ins:
            pos_erate.ins1 = ErrType(get_most_common_nt(sub_ins,
                                                        pos_erate.ref,
                                                        pos_erate.next),
                                     len(sub_ins)/depth, sub_align)
        sub_ins, sub_align = pull_nlength_indels_info(err_info, 2,
                                                      indel_alignment)
        if sub_ins:
            pos_erate.ins2 = ErrType(get_most_common_nt(sub_ins,
                                                        pos_erate.ref,
                                                        pos_erate.next),
                                     len(sub_ins)/depth, sub_align)
    elif err_type.startswith("sub"):
        # substitution error rate
        pos_erate.fullSub = f_erate
        sub_nt = set(err_info)
        pos_erate.sub = ErrType(sub_nt, f_erate, 0)
        if err_type == "sub_adjust":
            pos_erate.subAllErr = {"A": err_info.count("A")/depth,
                                   "C": err_info.count("C")/depth,
                                   "G": err_info.count("G")/depth,
                                   "T": err_info.count("T")/depth}
    else:
        print(("WARNING: Could not determine error type for error rate "
               "storage."), file=sys.stderr)

    return


def parse_pile(erate, mapped, pile_info, low=False, get_only=None):
    """Parse an ASCII alignment line

    Parameters
    ----------
    erate : Error class
        Stores the calculated error rate and properties
    mapped : str
        Actual alignment as a string represented by ASCII symbols
    pile_info : MpileInfo class
        Stores important information for next mpileup line
    low : boolean
        Specifies whether alignment is from primary or supplementary genome
        set
    get_only : str
        Retrieve error rate of this type

    Returns
    -------
    None
        parsed information is directly stored in the Error class

    """
    # keep track of insertions alignment position
    ins_pos = 1

    for i, char in enumerate(mapped):
        # keep track of reads in alignment line (by character), skip if
        # needed
        if pile_info.skip > 0:
            pile_info.skip = pile_info.skip - 1
            continue
        elif char == "^":
            pile_info.skip = 1              # skip mapping quality
            continue
        elif char == "$":
            erate.record_read_end()
            continue

        if char == "." or char == ",":
            # store alignment information
            erate.add_read_alignment()
            ins_pos <<= 1

            # store matching information
            pile_info.match.append(char)
        elif char == "-":
            # check length of deletion
            indel_len, next_i = check_indel_len(mapped, (i+1))
            pile_info.skip = indel_len + (next_i - (i+1))
        elif char == "+":
            # alignment information
            ins_pos ^= 1

            # check length of insertion and store insertion error rate
            indel_len, next_i = check_indel_len(mapped, (i+1))
            (pile_info
             .insStr
             .append("".join(mapped[next_i:(next_i+indel_len)])
                     .upper()))
            pile_info.skip = indel_len + (next_i - (i + 1))
        elif char == "*":
            # alignment information
            erate.add_read_alignment()
            ins_pos <<= 1
        else:
            # alignment information
            erate.add_read_alignment()
            ins_pos <<= 1

            # store substitutions error rate
            pile_info.subStr.append(char.upper())

    # pull error rates of each type of error
    if get_only is None:
        if pile_info.match:
            if len(pile_info.insStr) > 0:
                del pile_info.match[0: len(pile_info.insStr)]
            get_error_rate("total", erate, pile_info.match)

        if not low and pile_info.prevDel:
            get_error_rate("del", erate, pile_info.prevDel,
                           pile_info.delReads)
        elif low and pile_info.subPrevDel:
            get_error_rate("del", erate, pile_info.subPrevDel,
                           pile_info.subDelReads)
        else:
            pass

        if pile_info.insStr:
            get_error_rate("ins", erate, pile_info.insStr, ins_pos)

        if pile_info.subStr:
            get_error_rate("sub", erate, pile_info.subStr)
    else:
        if get_only == "del" and not low and pile_info.prevDel:
            get_error_rate("del", erate, pile_info.prevDel,
                           pile_info.delReads)
        elif get_only == "del" and low and pile_info.subPrevDel:
            get_error_rate("del", erate, pile_info.subPrevDel,
                           pile_info.subDelReads)
        elif get_only == "ins" and pile_info.insStr:
            get_error_rate("ins", erate, pile_info.insStr, ins_pos)
        elif get_only == "sub" and pile_info.subStr:
            get_error_rate("sub_adjust", erate, pile_info.subStr)
        else:
            return


def parse_del_only(pile_info, mapped, low=False):
    """Determine deletions in alignments that start next position

    Deletions and the respective lengths refer to deletions that start at
    the next position in the assembly

    Parameters
    ----------
    pile_info : MpileInfo class
        Class needed to store parsed information for the next position
    mapped : str
        ASCII string representing the alignment at the position
    low : boolean
        Specifies if mapped information is from the primary or supplementary
        reference genome alignments

    Returns
    -------
    None
        Deletion information is stored within MpileInfo class

    """
    pile_info.skip = 0

    for i, char in enumerate(mapped):
        if pile_info.skip > 0:
            pile_info.skip -= 1
            continue
        elif char == "^":               # Mapped sequence begins
            pile_info.skip = 1
            continue                    # Mapped sequence ends
        elif char == "$":
            pile_info.remove_prev_del_read(low)
            continue

        if char == "-":
            pile_info.adjust_prev_read_as_del(low)
            indel_len, next_i = check_indel_len(mapped, (i+1))
            pile_info.append_prev_del_info(
                    "".join(mapped[next_i:(next_i+indel_len)]).upper(), low)
            pile_info.skip = indel_len + (next_i - (i+1))
        elif char == "+":
            indel_len, next_i = check_indel_len(mapped, (i+1))
            pile_info.skip = indel_len + (next_i - (i+1))
        else:
            pile_info.record_prev_del_read(low)
            continue

    return


def parse_read_end_only(mapped):
    """Determines sequences ending this position

    Stripped derivative of parse_del_only function that focuses on mapped
    sequences that are ending in the current position

    Parameters
    ----------
    mapped : str
        ASCII string representing the alignment at the position

    Returns
    -------
    int
        Binary alignment of positions of sequences that end

    """
    skip = 0
    read_align = 1

    for i, char in enumerate(mapped):
        if skip > 0:
            skip -= 1
            continue
        elif char == "^":       # Mapped sequence begins
            skip = 1
            continue
        elif char == "$":       # Mapped sequence ends
            read_align ^= 1
            continue

        if char == "-":
            indel_len, next_i = check_indel_len(mapped, (i+1))
            skip = indel_len + (next_i - (i+1))
        elif char == "+":
            indel_len, next_i = check_indel_len(mapped, (i+1))
            skip = indel_len + (next_i - (i+1))
        else:
            # record mapped sequence
            read_align <<= 1
            continue

    return read_align


def parse_mpileup_line(pile_info, mpile, low):
    """Parse the mpileup line and calculate error rates

    Parameters
    ----------
    pile_info : MpileInfo class
        Stored error information from the immediate previous position
    mpile : str
        Line to be parsed
    low : boolean
        Specifies if alignment file is the primary or supplementary

    Returns
    -------
    Error class
        Class containing the error rate and other properties

    """
    # Diagnostic counter
    global _low_depth_count

    # Note: python is 0-based but samtools mpile is 1-based
    mpile_line = mpile[0].split()
    depth = int(mpile_line[3])

    # initialize an Error object with contig, ref nucleotide and depth info
    pos_erate = Error(mpile_line[0], int(mpile_line[1]) - 1, mpile_line[2],
                      depth)

    if depth < _depth_thres:
        if args.low and low:
            _low_depth_count += 1
        elif not args.low:
            _low_depth_count += 1
        else:
            pass
        if not low:
            pos_erate.set_low_depth()

        # low Depth. Only get read ends and set everything to 0
        pos_erate.readEnd = parse_read_end_only(mpile_line[4])
        get_error_rate("total", pos_erate)

        # look for deletions that start next position
        if depth < 1:
            pile_info.set_next_del_to_zero(low)
            return pos_erate
        else:
            pile_info.set_next_del_to_zero(low)
            parse_del_only(pile_info, mpile_line[4], low)
            return pos_erate
    else:
        pass

    # Temporary information
    try:
        pos_erate.next = mpile[1].split()[2]
    except AttributeError:
        pos_erate.next = "-"

    # get error rates
    parse_pile(pos_erate, mpile_line[4], pile_info, low)

    # reset storage class for next mpileup line
    pile_info.reset_class()
    pile_info.set_next_del_to_zero(low)

    # delete store next nt (unnecessary for future calculations)
    del pos_erate.next

    # look for deletions that start next position
    parse_del_only(pile_info, mpile_line[4], low)

    return pos_erate


def populate_erate_list(pile_info, erate_list, main_mpile, sub_mpile):
    """Determine order to populate list

    Both calculated main and sub alignment errors rates are ordered based on
    depth. More specifically, positions where depth is below threshold in
    the primary alignment, the supplementary error information is set as the
    main error information. The primary alignment is not discarded but is
    placed as an Error subclass within the main Error class. The inverse is
    also true when the primary alignment is the main error information.
    Although very messy, this structure is to prevent collision between the
    two alignments at the edges of low depth regions.

    Sub alignment information is only kept within n (max search window size)
    positions from a low depth region to save memory. When the primary
    alignments are the main source of error information, supplementary
    error information is stored in as a list (in MpileInfo class) until n
    element is exceeded, in which the oldest element is discarded, or when
    reaching a low depth region and information is added to previous
    positions.

    Visually (n = 5):

    P - Primary; S - Supplementary; H - High depth; * - Low depth

     Sub Class  ----SSSSSPPPPPPPPPPPPPPPPPPPPPPPPPPPSSSSS--------
                         |                         |
    Main Class  PPPPPPPPPSSSSSSSSSSSSSSSSSSSSSSSSSSSPPPPPPPPPPPPP
                         |                         |
         Depth  HHHHHHHHH***************************HHHHHHHHHHHHH

    Parameters
    ----------
    pile_info : MpileInfo class
        Stored information from the immediate upstream position
    erate_list : list
        Calculated error rates of previous position
    main_mpile : Error class
        Calculated error rate from the primary genome alignment file
    sub_mpile : Error class
        Calculated error rate from the supplementary genome alignment

    Returns
    -------
    None

    """
    if main_mpile.depth < 1 and sub_mpile.depth < 1:
        # both low depth. Reset everything
        main_mpile.lowErates = sub_mpile
        erate_list.append(main_mpile)
        if pile_info.lowCounter and pile_info.lowErates:
            pile_info.add_erate_to_err_rates(erate_list)
            pile_info.reset_low_erate()
    elif main_mpile.low and sub_mpile.depth > 0:
        # low depth region
        if pile_info.lowErates:
            pile_info.add_erate_to_err_rates(erate_list)
            pile_info.reset_low_erate()
            pile_info.set_low_counter()

        sub_mpile.set_low_depth()
        sub_mpile.lowErate = main_mpile
        erate_list.append(sub_mpile)
    else:
        # add to right end
        if pile_info.lowCounter \
                and len(pile_info.lowErates) > ((args.passes + 2) * 10 + 1):
            pile_info.add_erate_to_err_rates(erate_list)
            pile_info.reset_low_erate()

        pile_info.save_low_erate((args.passes + 2) * 10 + 1, sub_mpile)
        erate_list.append(main_mpile)


def fill_hp_info(err_rates, last=False):
    """Fill homopolymer and length information

    This function evaluates the nucleotide and homopolymer information from
    the previous position and determines the homopolymer information for the
    latest position. If a new homopolymer starts, hp_start attribute is set
    to last position and length is set to 1. If continuing a homopolymer,
    hp_start is set to hp_start of previous position and hp_length is
    incremented at the hp start position.

    Parameters
    ----------
    err_rates : list
        Calculated errors as Error classes
    last : boolean
        Specifies if the current position is the last one in a contig

    Returns
    -------
    None
        Error classes is directly modified

    """
    index = len(err_rates) - 1
    ref_nt = err_rates[index].ref

    if index == 0:
        err_rates[0].hp_start = 0
        return
    else:
        if ref_nt == err_rates[index-1].ref:
            err_rates[index].hp_start = err_rates[index-1].hp_start
        else:
            prev_hp_start = err_rates[index-1].hp_start
            err_rates[prev_hp_start].hp_len = index - prev_hp_start
            err_rates[index].hp_start = index

    if last:
        hp_length = 1
        while index != 0:
            if ref_nt != err_rates[index-hp_length].ref:
                start = err_rates[index].hp_start
                err_rates[start].hp_len = hp_length
                break
            else:
                hp_length += 1
    else:
        return

    return


def get_mpile_error_rates(main_mpile, sub_mpile=None):
    """Calculate error rates per position

    If available, both main and sub alignment files are parsed together,
    otherwise the main alignment file is parsed alone. This function
    cycles through the files line by line calculating error rates. Parsing
    of the true position needs to lag behind one line to ensure that error
    starting position are

    Parameters
    ----------
    main_mpile : str
        Path to the primary reference genome alignment file
    sub_mpile : str
        Path to the low depth supplementary alignment file

    Returns
    -------
    list of Error classes
        Calculated error rates with each element as a separate position

    """
    # Initialization
    err_rates = []
    container = MpileInfo()

    # iterate over mpiles to get error rates
    if sub_mpile != "":
        # parse at same time
        with open(main_mpile) as main_file:
            with open(sub_mpile) as sub_file:
                # get first lines
                main_to_parse, sub_to_parse = next(main_file), next(sub_file)
                for mline, sline in zip(main_file, sub_file):
                    main_erate = parse_mpileup_line(container,
                                                    [main_to_parse, mline],
                                                    False)
                    sub_erate = parse_mpileup_line(container,
                                                   [sub_to_parse, sline],
                                                   True)
                    populate_erate_list(container, err_rates, main_erate,
                                        sub_erate)
                    fill_hp_info(err_rates)
                    main_to_parse, sub_to_parse = mline, sline

                # parse the last line
                main_erate = parse_mpileup_line(container,
                                                [main_to_parse, None],
                                                False)
                sub_erate = parse_mpileup_line(container, [sub_to_parse, None],
                                               True)
                populate_erate_list(container, err_rates, main_erate, sub_erate)
                fill_hp_info(err_rates, True)
    else:
        # parse only main mpile
        with open(main_mpile) as main_file:
            main_to_parse = next(main_file)
            for mline in main_file:
                main_erate = parse_mpileup_line(container,
                                                [main_to_parse, mline],
                                                False)
                err_rates.append(main_erate)
                fill_hp_info(err_rates)
                main_to_parse = mline

            # parse last line
            main_erate = parse_mpileup_line(container,
                                            [main_to_parse, None],
                                            False)
            err_rates.append(main_erate)
            fill_hp_info(err_rates, True)

    # finished getting all info
    return err_rates


def retrieve_pos_erate(err_rates, pos, adjust, low):
    """Retrieve either error rate from primary or supplementary dataset

    Mostly important for combination-based error detection

    Parameters
    ----------
    err_rates : list
        All calculated error rates and its respective
        properties
    pos : int
        Position to pull
    adjust : boolean
        Indicate there is no supplementary dataset when adjusting errors
        detected
    low : boolean
        Determine which dataset to pull error information

    Returns
    -------
    Error class
        Error information belonging to the specified position

    """
    if args.low and not adjust and err_rates[pos].low != low:
        try:
            return err_rates[pos].lowErate
        except AttributeError:
            return err_rates[pos]
    else:
        return err_rates[pos]


def adjust_read_alignment(read_align, read_ends):
    """Reduces binary sequence alignment if an sequence ends

    To do so, alignments are converted to a binary string  and split into a
    list to be directed adjusted at specific indices. Once adjusted for
    sequencing ending this position, the list is joined back into a binary
    string to be converted into an integer.

    Companion function to reconstruct_prev_alignment

    Example
    -------

    Seq_ending  0010010   (Two sequence ended)
     Alignment  1110111
                -------
     Fin_align  11011


    Parameters
    ----------
    read_align : int
        Binary alignment to be adjusted
    read_ends : int
        Binary alignment of sequences ending this position

    Returns
    -------
    int
        Shortened binary alignment

    """
    b_read_ends = bin(read_ends)
    next_align = list(bin(read_align))
    if len(b_read_ends) != len(next_align):
        print("BUG: Mismatch in alignment and depth in adjust_read_alignment",
              file=sys.stderr)

    for bit in range(3, len(bin(read_ends))):
        if b_read_ends[bit] == "1":
            # read at this bit ends. Remove from alignment
            next_align[bit] = ""
        else:
            pass

    return int("".join(next_align), 2)


def get_uniq_err_rates(pos_erate, err_type, err_len, alignment,
                       uniq_pos_only=False):
    """Retrieve specified error rates and error alignments

    This function only returns errors found from sequences not containing an
     error (unique error)

    Parameters
    ----------
    pos_erate : Error class
        A singlar position's error rate and its respective properties
    err_type : str
        The error type to consider during retrieval
    err_len : int
        Error length to consider during retrieval
    alignment : int
        Binary record of aligned sequences and whether an error has been
        found
    uniq_pos_only : boolean
        Return unique error alignments only

    Returns
    -------
    int
        Error rate of the unique errors
    int
        Binary record of aligned sequences with a unique error

    """
    try:
        if err_type == "del":
            if err_len == 1:
                if uniq_pos_only:
                    return pos_erate.del1.reads_w_uniq_err(alignment)
                return (pos_erate.del1.uniq_err_erate(alignment),
                        pos_erate.del1.subset_out_uniq_err(alignment))
            elif err_len == 2:
                if uniq_pos_only:
                    return pos_erate.del2.reads_w_uniq_err(alignment)
                return (pos_erate.del2.uniq_err_erate(alignment),
                        pos_erate.del2.subset_out_uniq_err(alignment))
            else:
                if uniq_pos_only:
                    return 0
                return pos_erate.fullDel, alignment
        elif err_type == "ins":
            if err_len == 1:
                if uniq_pos_only:
                    return pos_erate.ins1.reads_w_uniq_err(alignment)
                return (pos_erate.ins1.uniq_err_erate(alignment),
                        pos_erate.ins1.subset_out_uniq_err(alignment))
            elif err_len == 2:
                if uniq_pos_only:
                    return pos_erate.ins2.reads_w_uniq_err(alignment)
                return (pos_erate.ins2.uniq_err_erate(alignment),
                        pos_erate.ins2.subset_out_uniq_err(alignment))
            else:
                if uniq_pos_only:
                    return 0
                return pos_erate.fullIns, alignment
        else:
            print(("BUG: Could not determine error type in get_uniq_err_"
                   "rates. DEBUG."), file=sys.stderr)
    except AttributeError:
        if uniq_pos_only:
            return 0
        else:
            return 0, alignment


def reconstruct_prev_alignment(read_align, read_ends):
    """Reconstruct an upstream bit alignment

    This function recreates an upstream bit alignment based on reads
    that have previously ended. To do so, alignments are converted to binary
    and split into a list to be directed adjusted at specific indices. Once
    reconstructed, list is joined back into a binary string to be converted
    into an integer.

    The bit shifting conversion method is slower when less than 1 million
    conversions are made and thus was not implemented

    Companion function to adjust_read_alignment

    Examples
    --------

    Ex. 1

    Prev_ended  0010010     (Two sequence ended)
    Curr_align  11111       (No new sequences added)
                -------
     Fin_align  1101101   (Old sequence restored)


    Ex. 2

    Prev_ended  0000010     (One sequence ended)
    Curr_align  00111111    (Two new sequences added)
                ---------
     Fin_align  001110111   (Old sequence restored; new sequences remain)

    Parameters
    ----------
    read_align : int
        Binary alignment to be adjusted
    read_ends : int
        Binary alignment of sequences that have ended from the previous
        position (0x0 indicates no reads ended)

    Returns
    -------
    int
        Reconstructed binary alignment

    """
    b_read_ends = bin(read_ends)
    next_align = list(bin(read_align))
    align_index = 3
    while align_index < len(b_read_ends):
        if b_read_ends[align_index] == "1":
            next_align[align_index:align_index] = ["0"]
        align_index += 1

    return int("".join(next_align), 2)


def join_single_errs(err_rates, pos_list, low):
    """Form temporary errors of length 2 if on the same aligned sequence

    This function should only be involved when different length searching is
    implemented and when the original error length considered is 2 nt. This
    function merges two errors 1 nt errors if are from the same mapped
    sequence. To ensure that errors are on the same sequence, recorded
    error alignment positions are compared after adjustments to account for
    depth changes between the two positons.

        Align 1 -------------------------------> Adj_align 1
          Pos 1                                  Pos 2
    Adj_Align 2 <------------------------------- Align 2

    Alignment 1 is compared to adjusted alignment 2
    Alignment 2 is compared to adjusted alignment 1

    If both comparisons agree that error position is on the same sequence,
    the errors are joined into error length of 2.

    Parameters
    ----------
    err_rates : list
        Calculated error rates and respective properties
    pos_list : list
        The two positions with single errors
    low : boolean
        Indicates if the error rate are from the primary or supplementary
        alignment file

    Returns
    -------
    int
        error rate of merged errors
    int
        Alignment positions were errors were joined

    """
    pos_candidate = 0
    uniq_errs = 0

    for first in range(0, len(pos_list) - 1):
        if type(pos_list[first][0]).__name__ == "list":
            continue
        else:
            for second in range(first + 1, len(pos_list)):
                if type(pos_list[second][0]).__name__ == "list":
                    continue
                else:
                    pos1 = pos_list[first][0]
                    to_compare_v1 = pos_list[first][1]
                    pos2 = pos_list[second][0]
                    to_compare_v2 = pos_list[second][1]

                    # first step only affects forward 
                    pos1_read_end = retrieve_pos_erate(err_rates, pos1,
                                                       False, low).readEnd
                    if (pos1_read_end
                            ^ (1 << (pos1_read_end.bit_length() - 1))) > 0:
                        to_compare_v1 = adjust_read_alignment(to_compare_v1,
                                                              pos1_read_end)

                    for step in range(1, (pos2 - pos1 + 1)):
                        pos1_read_end = retrieve_pos_erate(err_rates, pos1+step,
                                                           False, low).readEnd

                        # pos 1: adjust for new reads
                        bit_diff = (pos1_read_end.bit_length()
                                    - to_compare_v1.bit_length())
                        if bit_diff > 0:
                            to_compare_v1 <<= bit_diff
                        elif bit_diff < 0:
                            print(("BUG: Alignment bit difference in function "
                                   "join_single_errs"), file=sys.stderr)
                        else:
                            pass

                        # pos 1: adjust for reads ending
                        if (pos1 + step < pos2) \
                                and (pos1_read_end
                                     ^ (1 << (pos1_read_end.bit_length()
                                              - 1))) > 0:
                            to_compare_v1 = adjust_read_alignment(to_compare_v1,
                                                                  pos1_read_end)

                        # Pos 2: step backwards for adjustments
                        pos2_read_end = retrieve_pos_erate(err_rates, pos2-step,
                                                           False, low).readEnd
                        if (pos2_read_end
                                ^ (1 << (pos2_read_end.bit_length() - 1))) > 0:
                            to_compare_v2 = reconstruct_prev_alignment(
                                to_compare_v2, pos2_read_end)

                        # Pos 2: remove any new reads that were incorporated
                        bit_diff = (to_compare_v2.bit_length()
                                    - pos2_read_end.bit_length())
                        if bit_diff > 0:
                            to_compare_v2 >>= bit_diff
                        elif bit_diff < 0:
                            print(("BUG: Alignment bit difference in function "
                                   "join_single_errs"), file=sys.stderr)
                        else:
                            pass

                    # make sure that corrected alignment are the same length
                    if (pos_list[first][1]).bit_length() \
                            != to_compare_v2.bit_length():
                        print(("BUG: Alignment bit difference in forward "
                               "alignment in function join_single_errs"),
                              file=sys.stderr)
                    if (pos_list[second][1]).bit_length() \
                            != to_compare_v1.bit_length():
                        print(("BUG: Alignment bit difference in reverse "
                               "alignment in function join_single_errs"),
                              file=sys.stderr)

                    uniq1 = (pos_list[first][1]
                             & (to_compare_v2
                                ^ (1 << (to_compare_v2.bit_length() - 1))))
                    uniq2 = (pos_list[second][1]
                             & (to_compare_v1
                                ^ (1 << (to_compare_v1.bit_length() - 1))))

                    if bin(uniq1).count("1") != bin(uniq2).count("1"):
                        print(("BUG: Mismatch between forward and reverse bit "
                               "alignments in function join_single_errs"),
                              file=sys.stderr)

                    # nothing to merge
                    if uniq1 == 0:
                        continue

                    # merge and correct pos 1 and 2 for future corrections
                    pos_list.append([[pos1, pos2], [uniq1, uniq2]])
                    pos_list[first][1] ^= uniq1
                    pos_list[second][1] ^= uniq2

                    # record to return 
                    pos2_erate = retrieve_pos_erate(err_rates, pos2, False, low)
                    pos_candidate += bin(uniq2).count("1")/pos2_erate.depth
                    uniq_errs ^= uniq2

    # remove elements that doesn't have any errors to merge
    for pos in range(len(pos_list)-1, -1, -1):
        if type(pos_list[pos][0]).__name__ == "list":
            continue
        else:
            align = (1 << (pos_list[pos][1].bit_length() - 1))
            if (pos_list[pos][1] ^ align) == 0:
                del pos_list[pos]

    return pos_candidate, uniq_errs


def add_single_length_errors(err_rates, current_pos, current_align, pos_list,
                             low):
    """Calculates an error rate of supplementary shorter errors

    This function should only be involved when different length searching is
    implemented and when the original error length considered is 2 nt. In
    the event where errors of the same length are not found, shorter errors
    may be considered. All potential single length errors in sequences
    not containing a error are combined. Using depth, the error rate and
    alignment position of these shorter errors are returned

    Parameters
    ----------
    err_rates : list
        Calculated error rate and respective error properties
    current_pos : int
        Current position
    current_align : int
        Record of sequences that contain or do not contain an error as
        binary
    pos_list : list
        A list of positions that contain a single length error
    low : boolean
        Determines if the error rate is from primary or the supplementary
        genome set

    Returns
    -------
    int
        error rate of all single errors from sequences without an error of
        length 2
    int
        Alignment position(s) of those single errors
    """
    single_err = 0
    for single in pos_list:
        single_pos = single[0]
        to_compare = single[1]

        # Alignment adjustment to account for depth changes
        if type(single_pos).__name__ == "list":
            continue
        else:
            for step in range(single_pos, current_pos+1):
                pos1_read_end = retrieve_pos_erate(err_rates, step, False,
                                                   low).readEnd
                # New aligned sequences
                bit_diff = pos1_read_end.bit_length() - to_compare.bit_length()
                if bit_diff > 0:
                    to_compare <<= bit_diff
                elif bit_diff < 0:
                    print(("BUG: Alignment bit difference in function "
                           "add_single_length_errors"), file=sys.stderr)
                else:
                    pass

                # Ending sequences
                if (step < current_pos) \
                        and (pos1_read_end
                             ^ (1 << (pos1_read_end.bit_length() - 1))) > 0:
                    to_compare = adjust_read_alignment(to_compare,
                                                       pos1_read_end)

            # compare adjusted, check, and subset out unique single errors
            uniq_err = current_align & to_compare
            if uniq_err != 0:
                uniq_err ^= (1 << (uniq_err.bit_length() - 1))
            current_align ^= uniq_err
            single_err += bin(uniq_err).count("1") / err_rates[single_pos].depth

    return single_err, current_align


def winshift_nlen_set_erate_zero(err_rates, err_type, err_len, pos, alignment,
                                 low):
    """Adjust error rates of the positions used different length searching

    This function is a specialized version of set_erates_to_zero for the
    different length searching. Errors of length 2 can be split as two errors
    of length 1. A split is treated as a new specialized case of length 1
    (del3 and ins3) while error of length 2 is set as zero to indicate 'use'
    Split errors can not be used as a potential position for error
    positioning
    
    Parameters
    ----------
    err_rates : list
        Calculated error rates and the respective error properties
    err_type : str
        Error type to consider when setting error rate to zero
    err_len : int
        Error length to adjust
    pos : int
        Position to adjust info
    alignment : int
        Current record of sequences containing error that was 'used' to
        calculate the final error rate
    low : boolean
        Specifies if error rate adjustment comes from the primary or
        supplementary dataset

    Returns
    -------
    None
        Error classes are modified directly

    """
    pos_erate = retrieve_pos_erate(err_rates, pos, False, low)
    try:
        if err_len == 1:
            if err_type == "del":
                if ((pos_erate.del1.reads ^ alignment)
                        ^ (1 << (pos_erate.del1.reads.bit_length() - 1))) == 0:
                    del pos_erate.del1
                else:
                    pos_erate.del1.erate -= bin(alignment).count("1")\
                                            / pos_erate.depth
                    pos_erate.del1.reads ^= alignment
            elif err_type == "ins":
                if ((pos_erate.ins1.reads ^ alignment)
                        ^ (1 << (pos_erate.ins1.reads.bit_length() - 1))) == 0:
                    del pos_erate.ins1
                else:
                    pos_erate.ins1.erate -= bin(alignment).count("1")\
                                            / pos_erate.depth
                    pos_erate.ins1.reads ^= alignment
        elif err_len == 2:
            if err_type == "del":
                old_reads = pos_erate.del2.reads
                if ((pos_erate.del2.reads ^ alignment)
                        ^ (1 << (pos_erate.del2.reads.bit_length() - 1))) == 0:
                    del pos_erate.del2
                else:
                    pos_erate.del2.erate -= bin(alignment).count("1")\
                                            / pos_erate.depth
                    pos_erate.del2.reads ^= alignment

                # Split a 2 nt error to a specialized single error
                try:
                    if ((pos_erate.del3.reads ^ alignment)
                            ^ (1 << (pos_erate.del3.reads.bit_length()-1)) == 0):
                        del pos_erate.del3
                    else:
                        pos_erate.del3.erate -= bin(alignment).count("1")\
                                                / pos_erate.depth
                        pos_erate.del3.reads ^= alignment
                except AttributeError:
                    pos_erate.del3 = ErrType("N",
                                             (bin(alignment).count("1")
                                              / pos_erate.depth),
                                             (1 << (old_reads.bit_length() - 1))
                                             | alignment)
            elif err_type == "ins":
                old_reads = pos_erate.ins2.reads
                if ((pos_erate.ins2.reads ^ alignment)
                        ^ (1 << (pos_erate.ins2.reads.bit_length() - 1)) == 0):
                    del pos_erate.ins2
                else:
                    pos_erate.ins2.erate -= bin(alignment).count("1")\
                                            / pos_erate.depth
                    pos_erate.ins2.reads ^= alignment
                try:
                    if ((pos_erate.ins3.reads ^ alignment)
                            ^ (1 << (pos_erate.ins3.reads.bit_length()-1)) == 0):
                        del pos_erate.ins3
                    else:
                        pos_erate.ins3.erate -= bin(alignment).count("1")\
                                                / pos_erate.depth
                        pos_erate.ins3.reads ^= alignment
                except AttributeError:
                    pos_erate.ins3 = ErrType("N",
                                             (bin(alignment).count("1")
                                              / pos_erate.depth),
                                             (1 << (old_reads.bit_length() - 1))
                                             | alignment)
            else:
                print(("BUG: Issue with determining error type in "
                       "winshift_nlen_set_erate_zero for position {}")
                      .format(pos))
        else:
            print(("BUG: Issue with determining error length in "
                   "winshift_nlen_set_erate_zero for position {}")
                  .format(pos))
    except AttributeError:
        print(("BUG: No indel attribute found in position. May have been "
               "prematurely removed"), file=sys.stderr)
    return


def get_singlar_indel_erate_in_range(pos_erate, err_type, err_len, pos, align,
                                     winshift_difflen, pos_info_storage):
    """Retrieves a singular indel error rate of specified error type and length

    In the event of multiple lengths are considered, the error rate of split
    errors are included (del3 and ins3) as part of the respective error
    length (1 nt) rather than the original error length (2 nt)

    Parameters
    ----------
    pos_erate : Error class
        A position's calculated error rate and respective
        properties
    err_type : str
        The error type to consider
    err_len : int
        The error length to consider
    pos : int
        The position of error candidate if the error rate is greater than 0
    align : int
        Current record of aligned reads with and without errors
    winshift_difflen : int
        Used as a boolean to indicate if different error lengths search was
        implemented during this round of error detection
    pos_info_storage : list
        All positions with any errors of the specified type and length

    Returns
    -------
    float
        The position's error rate of type err_type and length err_len
    int
        Updated alignment information of 'reads' with errors
    int
        Alignment of 'reads' containing errors of type err_type and length
        err_len

    """
    uniq_align = get_uniq_err_rates(pos_erate, err_type, err_len, align,
                                    uniq_pos_only=True)
    pos_candidate, read_align = get_uniq_err_rates(pos_erate, err_type, err_len,
                                                   align)
    if pos_candidate > 0:
        pos_info_storage.append([pos, uniq_align])

    if winshift_difflen:
        if err_len == 1:
            if err_type == "del":
                try:
                    uniq_align_2 = pos_erate.del3.reads_w_uniq_err(read_align)
                    pos_candidate_2 = pos_erate.del3.uniq_err_erate(read_align)
                    read_align = pos_erate.del3.subset_out_uniq_err(read_align)
                except AttributeError:
                    pos_candidate_2 = 0
                    uniq_align_2 = 0
            elif err_type == "ins":
                try:
                    uniq_align_2 = pos_erate.ins3.reads_w_uniq_err(read_align)
                    pos_candidate_2 = pos_erate.ins3.uniq_err_erate(read_align)
                    read_align = pos_erate.ins3.subset_out_uniq_err(read_align)
                except AttributeError:
                    pos_candidate_2 = 0
                    uniq_align_2 = 0
            else:
                pos_candidate_2 = 0
                uniq_align_2 = 0

            if pos_candidate_2 > 0:
                pos_info_storage.append([pos * -1, uniq_align_2])
                uniq_align_2 ^= (1 << (uniq_align_2.bit_length()-1))
                pos_candidate += pos_candidate_2
        else:
            uniq_align_2 = 0
    else:
        uniq_align_2 = 0

    return pos_candidate, read_align, uniq_align ^ uniq_align_2


# TODO: Split the messy long function into smaller functions
def get_best_indel_erate_in_range(err_rates, err_type, err_len, start, end,
                                  hp=False, winshift_difflen_thres=0,
                                  uniq_erates=None, find=0, adjust=False):
    """Determines the maximum indel error rate of a search window

    This function is large due to separate handling of deletion and
    insertions. Also, error finding and adjustments are treated slightly
    differently. There are also many try exceptions as not all positions
    have all potential errors; reduced ErrType classes reduces memory load.

    Parameters
    ----------
    err_rates : list
        Calculated error rates as Error classes within the list
    err_type : str
        The error type to consider when calculating the final indel
        error rate
    err_len : int
        The error length to consider when calculating the final indel
        error rate
    start : int
        Start position of search window
    end : int
        End position of search window
    hp : boolean
        Specifies when you want to consider error rates as singlar position
        or as homopolymers
    winshift_difflen_thres : int
        Threshold for window shifting and different error length parsing. If
        0, parsing is not implemented. Should be used only when finding
        errors
    uniq_erates : int
        Alignment positions, stored as bits, of error rates used in final
        indel error rate
    find : int
        Original position of the error rate candidate. If 0, error rate is
        usually adjusting unless an edge case
    adjust : boolean
        Specifies if indel error rates is to be used for adjustment or not.
        Very Slight changes in implementation

    Returns
    -------
    int
        Calculated error rate for the specified search window

    """
    # Search window check and adjust if needed
    if start < 0:
        start = 0
    if end > len(err_rates):
        end = len(err_rates)
    if hp and err_type == "ins" and (start - 1) >= 0:
        orig_start = start - 2
    else:
        orig_start = start - 1

    # Initialzation of needed variables
    orig_pos_nlen = []
    err_candidate = 0
    if adjust:
        top_non_hp = []
    if winshift_difflen_thres:
        pos_diff_nlen = []

    # Extract original candidate information before starting search
    if find > 0:
        orig_pos = find
        find_low = err_rates[orig_pos].low
        orig_align = (1 << err_rates[orig_pos].readEnd.bit_length()) - 1
        (pos_candidate, orig_align,
         uniq_align) = get_singlar_indel_erate_in_range(err_rates[orig_pos],
                                                        err_type, err_len,
                                                        orig_pos, orig_align,
                                                        winshift_difflen_thres,
                                                        orig_pos_nlen)
        err_candidate += pos_candidate

        # Go backwards from candidate to start of search window to determine
        # alignment needed to start searching
        for reconstruct_pos in range(orig_pos-1, orig_start, -1):
            reconstruct_erate = retrieve_pos_erate(err_rates, reconstruct_pos,
                                                   adjust, find_low)
            reconstruct_read_end = reconstruct_erate.readEnd

            if (reconstruct_read_end
                    ^ (1 << (reconstruct_read_end.bit_length() - 1)) > 0):
                orig_align = reconstruct_prev_alignment(orig_align,
                                                        reconstruct_read_end)

            bit_diff = (orig_align.bit_length()
                        - reconstruct_read_end.bit_length())
            if bit_diff > 0:
                orig_align >>= bit_diff
            elif bit_diff < 0:
                print(("BUG: Alignment bit difference in function get_best_"
                       "indel_erate_in_range"), file=sys.stderr)
            else:
                pass

        read_align = orig_align                                 # two copies
    else:
        find_low = False
        read_align = (1 << err_rates[orig_start + 1].readEnd.bit_length()) - 1

    # Determine if an insertion in the upstream position is a part of the
    # first homopoylmer in the search window (incorporate into error rate)
    if hp and err_type == "ins":
        if (start - 1) >= 0:
            current_pos_erate = retrieve_pos_erate(err_rates, start - 1, adjust,
                                                   find_low)

            if err_len == 1:
                try:
                    nt = current_pos_erate.ins1.nt
                    current = (nt == err_rates[start].ref)
                except AttributeError:
                    current = False
            else:
                try:
                    nt = current_pos_erate.ins2.nt
                    current = (nt[-1] == err_rates[start].ref
                               and nt[0] != current_pos_erate.ref)
                except AttributeError:
                    current = False

            if current:
                (pos_candidate, read_align,
                 uniq_align) = get_singlar_indel_erate_in_range(
                    current_pos_erate, err_type, err_len, start-1, read_align,
                    0, orig_pos_nlen)
                err_candidate += pos_candidate
            else:
                pass

            # Alignment check
            if (current_pos_erate.readEnd
                    ^ (1 << current_pos_erate.readEnd.bit_length() - 1) > 0):
                if find > 0:
                    adjusted_read_end = \
                        current_pos_erate.readEnd \
                        >> (current_pos_erate.readEnd.bit_length()
                            - orig_align.bit_length())
                    orig_align = adjust_read_alignment(orig_align,
                                                       adjusted_read_end)
                read_align = adjust_read_alignment(read_align,
                                                   current_pos_erate.readEnd)
            else:
                pass
        else:
            pass
    else:
        pass

    # Go through search window
    for pos in range(start, end):
        current_pos_erate = retrieve_pos_erate(err_rates, pos, adjust, find_low)

        # Alignment adjustment
        bit_diff = (current_pos_erate.readEnd.bit_length()
                    - read_align.bit_length())
        if bit_diff > 0:
            # adjust for new read additions.
            read_align <<= bit_diff
            read_align ^= ((1 << bit_diff) - 1)

            if find > 0 and pos <= orig_pos:
                orig_uniq_align = orig_pos_nlen[0][1]
                for reconstruct_pos in range(orig_pos, pos - 1, -1):
                    reconstruct_erate = retrieve_pos_erate(err_rates,
                                                           reconstruct_pos,
                                                           adjust, find_low)
                    reconstruct_read_end = reconstruct_erate.readEnd
                    orig_uniq_align = reconstruct_prev_alignment(
                            orig_uniq_align, reconstruct_read_end)

                    # remove any new reads that were incorporated this round
                    bit_diff_2 = (orig_uniq_align.bit_length()
                                  - reconstruct_read_end.bit_length())
                    if bit_diff_2 > 0:
                        orig_uniq_align >>= bit_diff_2
                    elif bit_diff_2 < 0:
                        print(("BUG: Alignment bit difference in correcting "
                               "alignment during finding errors"),
                              file=sys.stderr)
                    else:
                        pass

                # alignment check
                if read_align.bit_length() != orig_uniq_align.bit_length():
                    print(("BUG: Alignment bit difference when determing "
                           "remaining intact genomes"), file=sys.stderr)

                # keep track of original alignment
                orig_uniq_align &= ((1 << bit_diff) - 1)
                read_align ^= orig_uniq_align
            else:
                pass
        elif bit_diff < 0:
            print(("BUG: Unique errors alignment length greater than read "
                   "depth on contig {} at position {}")
                  .format(err_rates[pos].contig, pos), file=sys.stderr)
        else:
            pass

        # low depth : conserve alignment only
        if current_pos_erate.depth < _depth_thres:
            if (current_pos_erate.readEnd
                    ^ (1 << current_pos_erate.readEnd.bit_length() - 1) > 0):
                if find > 0:
                    adjusted_read_end = (current_pos_erate.readEnd >>
                                         (current_pos_erate.readEnd.bit_length()
                                          - orig_align.bit_length()))
                    orig_align = adjust_read_alignment(orig_align,
                                                       adjusted_read_end)
                read_align = adjust_read_alignment(read_align,
                                                   current_pos_erate.readEnd)
            continue

        # Calculate indel
        if err_type == "del":
            if not adjust and abs(find) == pos:
                pass
            elif not adjust:
                (pos_candidate,
                 read_align, uniq_align) = get_singlar_indel_erate_in_range(
                    current_pos_erate, err_type, err_len, pos, read_align,
                    winshift_difflen_thres, orig_pos_nlen)
                err_candidate += pos_candidate
            else:
                if err_len == 1:
                    pos_candidate = 0
                else:
                    (pos_candidate, placeholder,
                     uniq_align) = get_singlar_indel_erate_in_range(
                        current_pos_erate, err_type, 1, pos, read_align,
                        winshift_difflen_thres, orig_pos_nlen)
                err_candidate += current_pos_erate.fullDel - pos_candidate
        elif hp and err_type == "ins":
            # TODO: redo this area since it is a mess
            # consider error rates as homopolymers
            uniq_align = get_uniq_err_rates(current_pos_erate, err_type,
                                            err_len, read_align,
                                            uniq_pos_only=True)
            hp_candidate, sub_align = get_uniq_err_rates(current_pos_erate,
                                                         err_type, err_len,
                                                         read_align)
            if not adjust and abs(find) == pos:
                hp_candidate = 0
            if hp_candidate > 0:
                if err_len == 1:
                    try:
                        err_nt = current_pos_erate.ins1.nt
                    except AttributeError:
                        err_nt = ""
                    if adjust:
                        hp_candidate = current_pos_erate.fullIns
                elif err_len == 2:
                    try:
                        err_nt = current_pos_erate.ins2.nt
                    except AttributeError:
                        err_nt = ""
                    if adjust:
                        hp_candidate, placeholder = get_uniq_err_rates(
                            current_pos_erate, err_type, 1, read_align)
                        hp_candidate = (current_pos_erate.fullIns
                                        - hp_candidate)
                else:
                    if (err_len - 2) == 1:
                        try:
                            err_nt = current_pos_erate.ins1.nt
                            sub_align = current_pos_erate\
                                .ins1.subset_out_uniq_err(sub_align)
                        except AttributeError:
                            err_nt = ""
                            sub_align = read_align
                    else:
                        try:
                            err_nt = current_pos_erate.ins2.nt
                            sub_align = current_pos_erate.ins2\
                                .subset_out_uniq_err(sub_align)
                        except AttributeError:
                            err_nt = ""
                            sub_align = read_align

                # Check for hp edges
                if (pos + 1) == end:
                    if err_nt == "":
                        err_candidate += 0
                    elif err_nt[0] == current_pos_erate.ref:
                        err_candidate += hp_candidate
                        read_align = sub_align
                        orig_pos_nlen.append([pos, uniq_align])
                    else:
                        next_hp = start + err_rates[start].hp_len
                        try:
                            next_hp_nt = err_rates[next_hp].ref
                        except IndexError:
                            next_hp_nt = ""

                        if err_nt[-1] == next_hp_nt:
                            err_candidate += 0
                        else:
                            # matches nothing so treat it as unique
                            if adjust:
                                if len(top_non_hp) > 0 \
                                        and top_non_hp[1] >= hp_candidate:
                                    err_candidate += 0
                                elif len(top_non_hp) > 0 \
                                        and top_non_hp[1] < hp_candidate:
                                    top_non_hp[0] = pos
                                    top_non_hp[1] = hp_candidate
                                else:
                                    top_non_hp.extend([pos, hp_candidate])
                            else:
                                err_candidate += hp_candidate
                                read_align = sub_align
                                orig_pos_nlen.append([pos, uniq_align])
                else:
                    if err_nt != "":
                        if (adjust
                                and err_nt[0] == current_pos_erate.ref
                                and err_nt[0] == err_nt[-1]):
                            err_candidate += hp_candidate
                        elif (adjust
                              and (err_nt[0] != current_pos_erate.ref
                                   or err_nt[-1] != current_pos_erate.ref)):
                            if len(top_non_hp) > 0 \
                                    and top_non_hp[1] >= hp_candidate:
                                err_candidate += 0
                            elif len(top_non_hp) > 0 \
                                    and top_non_hp[1] < hp_candidate:
                                top_non_hp[0] = pos
                                top_non_hp[1] = hp_candidate
                            else:
                                top_non_hp.extend([pos, hp_candidate])
                        else:
                            err_candidate += hp_candidate
                            read_align = sub_align
                            orig_pos_nlen.append([pos, uniq_align])
                    else:
                        err_candidate += 0
            else:
                err_candidate += 0
        elif err_type == "ins":
            if not adjust and abs(find) == pos:
                pass
            else:
                (pos_candidate, read_align,
                 uniq_align) = get_singlar_indel_erate_in_range(
                    current_pos_erate, err_type, err_len, pos, read_align,
                    winshift_difflen_thres, orig_pos_nlen)
                err_candidate += pos_candidate
        else:
            print("BUG: could not parse error type. Debug", file=sys.stderr)

        # check for alternative error length
        if winshift_difflen_thres:
            err_len_2 = err_len - (-1)**err_len
            pos_candidate, sub_align = get_uniq_err_rates(current_pos_erate,
                                                          err_type, err_len_2,
                                                          read_align)

            if pos_candidate > 0:
                # found potential corroborating errors
                uniq_align = get_uniq_err_rates(current_pos_erate, err_type,
                                                err_len_2, read_align,
                                                uniq_pos_only=True)
                pos_diff_nlen.append([pos, uniq_align])

                if err_len == 1:
                    err_candidate += pos_candidate
                    read_align = sub_align

                    # to match up with uniq_alignment of err_length 2
                    pos_diff_nlen[-1][1] ^= (1 << (uniq_align.bit_length()-1))
                elif err_len == 2:
                    if len(pos_diff_nlen) > 1:
                        # add potential 1 nt error candidates together
                        pos_candidate, sub_align = join_single_errs(
                            err_rates, pos_diff_nlen, low=find_low)
                        err_candidate += pos_candidate
                        read_align ^= sub_align
                else:
                    print("BUG: cannot parse error length. DEBUG",
                          file=sys.stderr)
            else:
                pass
        else:
            pass

        # Alignment check and adjust
        if (current_pos_erate.readEnd
                ^ (1 << current_pos_erate.readEnd.bit_length() - 1)) > 0:
            if find > 0:
                adjusted_read_end = current_pos_erate.readEnd \
                                    >> (current_pos_erate.readEnd.bit_length()
                                        - orig_align.bit_length())
                orig_align = adjust_read_alignment(orig_align,
                                                   adjusted_read_end)
            read_align = adjust_read_alignment(read_align,
                                               current_pos_erate.readEnd)
        else:
            pass

    # Check for intact genomes and adjust indel error
    if find > 0 and err_candidate > _indel_thres:
        orig_align = read_align \
                     >> (read_align.bit_length() - orig_align.bit_length())
        intact_genomes = (bin(orig_align).count("1") - 1)
        if intact_genomes > 0:
            alt_err_candidate = 1 - (intact_genomes/err_rates[orig_pos].depth)
            err_candidate = min(err_candidate, alt_err_candidate)

    # set erate to zero prematurely
    # Prevents consideration when determining error position
    # TODO: clean up. This is a mess
    if winshift_difflen_thres and pos_diff_nlen:
        if err_candidate > winshift_difflen_thres:
            if err_len == 1:
                for diff_err in pos_diff_nlen:
                    winshift_nlen_set_erate_zero(err_rates, err_type, 2,
                                                 diff_err[0], diff_err[1],
                                                 low=find_low)
                else:
                    pass
            elif err_len == 2:
                for two_pos in pos_diff_nlen:
                    if type(two_pos[0]).__name__ != "list":
                        continue
                    for n in range(2):
                        diff_err = [two_pos[0][n], two_pos[1][n]]
                        winshift_nlen_set_erate_zero(err_rates, err_type, 1,
                                                     diff_err[0], diff_err[1],
                                                     low=find_low)
                    else:
                        pass
                pass
            else:
                print("BUG: cannot parse error length. DEBUG", file=sys.stderr)
        else:
            # Check for shorter error length. Only for error length 2
            if err_len == 2:
                pos_candidate, sub_align = add_single_length_errors(
                    err_rates, pos, read_align, pos_diff_nlen, low=find_low)
                err_candidate += pos_candidate

                # Intact genome check
                if find > 0 and err_candidate > _indel_thres:
                    orig_align = sub_align >> (sub_align.bit_length()
                                               - orig_align.bit_length())
                    intact_genomes = (bin(orig_align).count("1") - 1)
                    if intact_genomes > 0:
                        alt_err_candidate = 1 - (intact_genomes
                                                 / err_rates[orig_pos].depth)
                        err_candidate = min(err_candidate, alt_err_candidate)

                # Re-check final error rate
                if err_candidate > winshift_difflen_thres:
                    for index in range(len(pos_diff_nlen)-1, -1, -1):
                        if type(pos_diff_nlen[index][0]).__name__ != "list":
                            continue
                        for n in range(2):
                            diff_err = [pos_diff_nlen[index][0][n],
                                        pos_diff_nlen[index][1][n]]
                            winshift_nlen_set_erate_zero(err_rates, err_type,
                                                         1, diff_err[0],
                                                         diff_err[1],
                                                         low=find_low)
                        # delete when done adjusting
                        del pos_diff_nlen[index]

                    # readjust error info of length 2
                    for diff_err in orig_pos_nlen:
                        diff_err[1] ^= (1 << (diff_err[1].bit_length()-1))
                        winshift_nlen_set_erate_zero(err_rates, err_type, 2,
                                                     diff_err[0], diff_err[1],
                                                     low=find_low)

                    # negative error rate indicates a length change
                    err_candidate *= -1
                    del orig_pos_nlen[:]
                else:
                    pass
            else:
                pass
    else:
        pass

    # pass unique errors to parent function through list
    if not adjust and uniq_erates is not None:
        # TODO: use copy.deepcopy
        if err_candidate < 0 and err_len == 2:
            # Alignment check and adjustment for consistency
            for e in pos_diff_nlen:
                if type(e[0]).__name__ == "list":
                    continue
                adjusted_read_align = read_align
                for reconstruct_pos in range(end - 1, e[0] - 1, -1):
                    reconstruct_erate = retrieve_pos_erate(err_rates,
                                                           reconstruct_pos,
                                                           adjust, find_low)
                    reconstruct_read_end = reconstruct_erate.readEnd
                    if (reconstruct_read_end
                            ^ (1 << (reconstruct_read_end.bit_length()-1)) > 0):
                        adjusted_read_align = reconstruct_prev_alignment(
                            adjusted_read_align, reconstruct_read_end)

                    bit_diff = (adjusted_read_align.bit_length()
                                - reconstruct_read_end.bit_length())
                    if bit_diff > 0:
                        adjusted_read_align >>= bit_diff
                    elif bit_diff < 0:
                        print(("BUG: Alignment bit difference during unique "
                               "error determination of positions used in error "
                               "finding"), file=sys.stderr)
                    else:
                        pass

                current_pos_erate = retrieve_pos_erate(err_rates, e[0], adjust,
                                                       find_low)
                e[1] = get_uniq_err_rates(current_pos_erate, err_type,
                                          err_len_2, adjusted_read_align,
                                          uniq_pos_only=True)
                uniq_erates.append(e)
        else:
            for e in orig_pos_nlen:
                uniq_erates.append(e)
    elif adjust and uniq_erates is not None and top_non_hp:
        uniq_erates.extend(top_non_hp)
    else:
        pass

    return err_candidate


def determine_range_error_info(err_rates, err_type, err_len, start, end, find=0,
                               hp=False, adjust=False, existing_errs=None,
                               check_diff_nlen=False):
    """Determines error properties

    This function determines the best error position, depth, erroneous
    nucleotide and suggested correction. All errors in the search window are
    considered and scored accordingly. Greater weight is given to errors
    with suggested corrections associated long homopolymers. Scores are also
    penalized for breaking homopolymers

    Parameters
    ----------
    err_rates : list
        Calculated error rates as Error classes within the list
    err_type : str
        Limit extracted error properties to this error type
    err_len : int
        Limit extracted error properties to this error length
    start : int
        Starting position of search window
    end : int
        Ending position of search window
    find : int
        Original error candidate position. TODO: change to low depth boolean
    hp : boolean
        Consider error rates as groups of homopolymer rather than individual
        position (similar to RLE)
    adjust : boolean
        Minor adjustments to scoring based on adjustment rather than
        detection
    existing_errs : list
        A list of positions that currently have an error. Avoids collisions
        with multiple errors of the same error type and length are within a
        small search window. Avoids choosing the best candidate twice
    check_diff_nlen : boolean
        Removes the error length limitation. Only lengths of 1 or 2
        nucleotides are considered

    Returns
    -------
    int
        The position of the final error. A negative integer indicates a
        repeat in final positioning
    int
        Alignment depth at the final position
    Char
        Erroneous nucleotide. A nucleotide indicates a suggested deletion.
        '-' indicates a suggested insertion
    Char
        Suggested correction. A nucleotide indicates a suggested insertion.
        '-' indicates a suggested deletion

    """
    # Range check
    if start < 0:
        start = 0
    if end > len(err_rates):
        end = len(err_rates)

    if hp and err_type == "ins":
        if (start - 1) >= 0:
            start -= 1

    if args.low:
        find_low = err_rates[find].low
    else:
        find_low = False

    # initialization
    score_info = {}
    hp_err_pos = {}
    if existing_errs is None:
        existing_errs = []

    for pos in range(start, end):
        # retrieve homopolymer information for current position
        curr_hp = err_rates[pos].hp_start
        curr_hp_nt = err_rates[pos].ref
        curr_hp_len = err_rates[curr_hp].hp_len

        if (pos + 1) < (curr_hp + curr_hp_len):
            next_hp = curr_hp
            next_hp_nt = curr_hp_nt
            next_hp_len = curr_hp_len
        else:
            try:
                next_hp = curr_hp + curr_hp_len
                next_hp_nt = err_rates[next_hp].ref
                next_hp_len = err_rates[next_hp].hp_len
            except IndexError:
                next_hp = 0
                next_hp_nt = ""
                next_hp_len = 0

        hp_length = 0
        hp_match = "none"

        # Position information
        current_pos_erate = retrieve_pos_erate(err_rates, pos, adjust, find_low)

        if err_type == "del":
            if err_len == 1:
                try:
                    err_nt = current_pos_erate.del1.nt
                except AttributeError:
                    if check_diff_nlen:
                        try:
                            err_nt = current_pos_erate.del3.nt
                        except AttributeError:
                            continue
                    else:
                        continue
                if not adjust:
                    try:
                        erate = current_pos_erate.del1.erate
                    except AttributeError:
                        erate = 0
                    try:
                        if check_diff_nlen:
                            erate += current_pos_erate.del3.erate
                    except AttributeError:
                        pass
                else:
                    erate = current_pos_erate.fullDel
            else:
                try:
                    err_nt = current_pos_erate.del2.nt
                except AttributeError:
                    continue
                if not adjust:
                    erate = current_pos_erate.del2.erate
                else:
                    try:
                        erate = current_pos_erate.fullDel \
                                - current_pos_erate.del1.erate
                    except AttributeError:
                        erate = current_pos_erate.fullDel
        elif err_type == "ins":
            if err_len == 1:
                try:
                    err_nt = current_pos_erate.ins1.nt
                except AttributeError:
                    if check_diff_nlen:
                        try:
                            err_nt = current_pos_erate.ins3.nt
                        except AttributeError:
                            continue
                    else:
                        continue

                if not adjust:
                    try:
                        erate = current_pos_erate.ins1.erate
                    except AttributeError:
                        erate = 0
                    try:
                        if check_diff_nlen:
                            erate += current_pos_erate.ins3.erate
                    except AttributeError:
                        pass
                else:
                    erate = current_pos_erate.fullIns
            else:
                try:
                    err_nt = current_pos_erate.ins2.nt
                except AttributeError:
                    continue

                if not adjust:
                    erate = current_pos_erate.ins2.erate
                else:
                    try:
                        erate = current_pos_erate.fullIns \
                                - current_pos_erate.ins1.erate
                    except AttributeError:
                        erate = current_pos_erate.fullIns
        else:
            print("BUG: Could not access error nucleotide at position {}"
                  .format(pos), file=sys.stderr)

        # Determine if error matches current, next or no homopolymers
        if err_len == 1:
            if err_nt == curr_hp_nt:
                hp_length = curr_hp_len
                hp_match = "current"
            elif err_nt == next_hp_nt:
                hp_length = next_hp_len
                hp_match = "next"
            else:
                hp_length = 0.5
        elif err_len > 1:
            if pos == (next_hp - 1):
                # between junction of current and next hp
                # Ex. suggested insertion of "AC"
                # NCCC -> NACCCC vs NAAA -> NACAAA: former is scored higher
                if curr_hp_nt == err_nt[0] and next_hp_nt == err_nt[-1]:
                    if curr_hp_len >= next_hp_len:
                        hp_length = curr_hp_len
                        hp_match = "current"
                    else:
                        if adjust:
                            # Adjustment edge case. Both hp are valid
                            hp_length = curr_hp_len
                            hp_match = "current"
                        else:
                            hp_length = next_hp_len
                            hp_match = "next"
                elif curr_hp_nt == err_nt[0] and next_hp_nt != err_nt[-1]:
                    hp_length = curr_hp_len
                    hp_match = "current"
                elif curr_hp_nt != err_nt[0] and next_hp_nt == err_nt[-1]:
                    hp_length = next_hp_len
                    hp_match = "next"
                else:
                    hp_length = 0.5
            else:
                if err_nt[0] == err_nt[-1]:
                    if curr_hp_nt == err_nt[0]:
                        hp_length = curr_hp_len
                        hp_match = "current"
                    else:
                        hp_length = 0.5
                else:
                    if curr_hp_nt == err_nt[0] or curr_hp_nt == err_nt[-1]:
                        hp_length = curr_hp_len/2
                        hp_match = "current"
                    else:
                        hp_length = 0.25
        else:
                print(("BUG: Can't tell error length to determine error "
                       "candidate score with length {} at position {}.")
                      .format(err_len, pos), file=sys.stderr)

        # calculate score
        if adjust:
            if err_type == "ins" and (pos + 1) == end and hp_match == "next":
                # skip if edge of hp and matches next hp
                continue
            else:
                if err_len > 1 and err_nt[0] != err_nt[-1]:
                    hp_match = "junction"

            score = erate
        else:
            score = hp_length * erate
            if pos in existing_errs:
                score = 0

        if hp_match == "current":
            if curr_hp in score_info:
                score_info[curr_hp] += score
                hp_err_pos[curr_hp][pos] = score
            else:
                score_info[curr_hp] = score
                hp_err_pos[curr_hp] = {pos: score}
        elif hp_match == "next":
            if next_hp in score_info:
                score_info[next_hp] += score
                hp_err_pos[next_hp][pos] = score
            else:
                score_info[next_hp] = score
                hp_err_pos[next_hp] = {pos: score}
        else:
            # Negative key = unique score
            key = 0 - pos
            score_info[key] = score

    # get the max score position
    try:
        hp_pos = max(score_info, key=score_info.get)
    except ValueError:
        if adjust:
            return -1, 0, "-", "-"
        else:
            print("BUG: Could not determine best position of error",
                  file=sys.stderr)

    if hp_pos > 0:
        final_pos = max(hp_err_pos[hp_pos], key=hp_err_pos[hp_pos].get)
    elif hp_pos == 0:
        try:
            final_pos = max(hp_err_pos[hp_pos], key=hp_err_pos[hp_pos].get)
        except KeyError:
            final_pos = abs(hp_pos)
    else:
        final_pos = abs(hp_pos)

    final_pos_erate = retrieve_pos_erate(err_rates, final_pos, adjust, find_low)

    # Retrieve the rest of error properties
    if err_type == "del":
        rep_nt = "-"
        if err_len == 1:
            ref_nt = final_pos_erate.ref
        else:
            ref_nt = final_pos_erate.del2.nt
    elif err_type == "ins":
        ref_nt = "-"
        if err_len == 1:
            try:
                rep_nt = final_pos_erate.ins1.nt
            except AttributeError:
                try:
                    rep_nt = final_pos_erate.ins3.nt[0]
                except AttributeError:
                    print(("BUG: Could not determine replacement nt at "
                           "position {}. DEBUG.").format(pos), file=sys.stderr)
                    rep_nt = "-"
        else:
            rep_nt = final_pos_erate.ins2.nt
    else:
        print(("BUG: cannot determine error type to pull relevant error info "
               "at function determine_range_error_info. Debug"),
              file=sys.stderr)

    if final_pos_erate.pos in existing_errs:
        return -1, final_pos_erate.depth, ref_nt, rep_nt
    else:
        return final_pos_erate.pos, final_pos_erate.depth, ref_nt, rep_nt


def set_erate_to_zero(err_rates, err_type, err_len, uniq_errs, orig_index):
    """Remove error rates that were used to find errors

    This function substracts the unique error rates used in determining an
    error. The alignments associated with these error rates are also "turned"
    off by converting to a bit of 0 indicating no errors present in those
    aligned sequences.

    Parameters
    ----------
    err_rates : list
        Calculated error rates as a list of Error classes
    err_type : str
        The error type to be adjusted
    err_len : int
        The error length to be adjusted
    uniq_errs : int
        Binary of alignment positions of associated errors
    orig_index : int
        Index of the original position candidate for determination if
        adjustments are made to the primary or supplementary error rates

    Returns
    -------
    None
        Error rates are directly modified

    """
    low = err_rates[orig_index].low
    for erate in uniq_errs:
        index = erate[0]
        erate_to_remove = ((bin(erate[1]).count("1") - 1)
                           / err_rates[abs(index)].depth)
        pos_erate = retrieve_pos_erate(err_rates, abs(index), False, low)

        if err_type == "del":
            if err_len == 1:
                if index >= 0:
                    pos_erate.del1.remove_uniq_errs(erate[1])
                    pos_erate.del1.erate -= erate_to_remove
                else:
                    pos_erate.del3.remove_uniq_errs(erate[1])
                    pos_erate.del3.erate -= erate_to_remove

            else:
                pos_erate.del2.remove_uniq_errs(erate[1])
                pos_erate.del2.erate -= erate_to_remove
        elif err_type == "ins":
            if err_len == 1:
                if index >= 0:
                    pos_erate.ins1.remove_uniq_errs(erate[1])
                    pos_erate.ins1.erate -= erate_to_remove
                else:
                    try:
                        pos_erate.ins3.remove_uniq_errs(erate[1])
                        pos_erate.ins3.erate -= erate_to_remove
                    except AttributeError:
                        print(("BUG: Could not find previously split indel "
                               "of length 2 in function set_erate_to_zero"),
                              file=sys.stderr)
            else:
                pos_erate.ins2.remove_uniq_errs(erate[1])
                pos_erate.ins2.erate -= erate_to_remove
        else:
            print(("BUG: error type could not be determined in function "
                   "set_erate_to_zero"), file=sys.stderr)

    del uniq_errs[:]
    return


def find_errors(err_rates, err_list, nt_range, sub_mpile, winshift_nlen=0,
                existing_errs=None):
    """Finds putative errors from error rate info within the search window

    This function cycles through every position in the error rate file to
    determine if a position or region contains an error. Substitution are
    considered independently of indels. At each position, the most frequent
    indel error type and length is evaluated. The alignment status, which
    alignment contains which indel error, is monitored as a int with an
    active error set as a bit of 1 and no error as a 0 bit. Positions with
    an error rate within 0.3 to error threshold are considered potential
    errors and the function will expand the search window to evaluate
    whether the region's error rate is above the error threshold. Error
    rates "used" the detection of an error are removed for future passes.

    Errors are stored with the following information (in order): Error type,
    Contig, Position, Depth, Error Length, Erroneous nucleotide, Suggested
    Correction, Error rate, Index for adjustment alignment files

    Parameters
    ----------
    err_rates : list
        Calculated error rates stored in an Error class
    err_list : list
        List of previously detected errors
    nt_range : int
        Current pass number. 0 = First pass through error rate list
    sub_mpile : str
        Supplementary alignment file used for low depth
    winshift_nlen : boolean
        Indicates whether search window shifting and evaluation of different
        error lengths should be implemented
    existing_errs : list
        Positions already found with errors to avoid collisions of errors
        with overlapping regions

    Returns
    -------
    None
        Pre-existng error list (err_list)is modified instead

    """
    # Initialization
    erates_used = []
    if existing_errs is None:
        existing_errs = []

    for i in range(0, len(err_rates)):
        # determine error threshold to use
        if sub_mpile == "" and err_rates[i].low:
            continue
        elif err_rates[i].low and err_rates[i].depth < _depth_thres:
            continue
        elif err_rates[i].low:
            thres = args.low_threshold
        else:
            thres = _indel_thres
        depth = err_rates[i].depth

        if nt_range < 0:
            # Independent. Substiutions
            try:
                sub = err_rates[i].fullSub
            except AttributeError:
                sub = 0

            if sub >= _sub_thres:
                err_list.append(["sub",
                                err_rates[i].contig,
                                err_rates[i].pos,
                                err_rates[i].depth,
                                1, err_rates[i].ref,
                                err_rates[i].sub.nt,
                                err_rates[i].sub.erate,
                                i])
        else:
            pass

        # Parse for indels
        try:
            del1 = err_rates[i].del1.erate
        except AttributeError:
            del1 = 0

        try:
            del2 = err_rates[i].del2.erate
        except AttributeError:
            del2 = 0

        try:
            ins1 = err_rates[i].ins1.erate
        except AttributeError:
            ins1 = 0

        try:
            ins2 = err_rates[i].ins2.erate
        except AttributeError:
            ins2 = 0

        # Determine error type and length to evaluate
        # Both an insertion and deletion should not occur at the same time
        if del1 > del2:
            err_del = del1
        elif del1 < del2:
            err_del = del2
        else:
            err_del = del1

        if ins1 > ins2:
            err_ins = ins1
        elif ins1 < ins2:
            err_ins = ins2
        else:
            err_ins = ins1

        if err_del > err_ins:
            indel_err = err_del
            err_type = "del"
            if err_del == del1:
                err_len = 1
            else:
                err_len = 2
        else:
            indel_err = err_ins
            err_type = "ins"
            if err_ins == ins1:
                err_len = 1
            else:
                err_len = 2

        # Evaluate position's error rate
        if indel_err >= thres:
            if nt_range < 0:
                read_align = (1 << err_rates[i].readEnd.bit_length()) - 1
                pos = err_rates[i].pos
                if err_type == "del":
                    rep_nt = "-"
                    if err_len == 1:
                        ref_nt = err_rates[i].ref
                    else:
                        ref_nt = err_rates[i].del2.nt
                else:
                    ref_nt = "-"
                    if err_len == 1:
                        rep_nt = err_rates[i].ins1.nt
                    else:
                        rep_nt = err_rates[i].ins2.nt

                # keep track to set to zero
                uniq_align = get_uniq_err_rates(err_rates[i], err_type, err_len,
                                                read_align, uniq_pos_only=True)
                erates_used.append([i, uniq_align])
                set_erate_to_zero(err_rates, err_type, err_len, erates_used, i)
            else:
                print(("BUG: Did not catch an indel error above threshold on "
                       "first pass at position {}. DEBUG").format(i),
                      file=sys.stderr)
                continue
        elif 0.25 <= indel_err < thres and (indel_err * depth) >= 2:
            # Error candidates
            if nt_range < 0:
                continue
            elif nt_range == 0:
                # set search window to hp length
                start = err_rates[i].hp_start
                end = start + err_rates[start].hp_len
            elif nt_range > 0:
                start = i - nt_range
                end = i + nt_range + 1
            else:
                pass

            # TODO: condense by multiplying this by boolean
            if winshift_nlen:
                winshift_nlen_thres = thres
            else:
                winshift_nlen_thres = 0

            indel_err = get_best_indel_erate_in_range(err_rates, err_type,
                                                      err_len, start, end,
                                                      (not nt_range),
                                                      winshift_nlen_thres,
                                                      uniq_erates=erates_used,
                                                      find=i)

            if indel_err > 1.2:
                continue
            elif indel_err < 0:
                # Found evidence for other error length.
                # Flip length and indel error value.
                err_len -= (-1)**err_len
                indel_err = abs(indel_err)
            elif (not winshift_nlen or winshift_nlen < 0) \
                    and nt_range > 0 \
                    and indel_err < thres:
                # Window shifting. Evaluate upstream and downstream windows
                eused2, eused3 = [], []

                # determine degree of window shifting
                left_start = start - int(nt_range/2)
                left_end = end - int(nt_range/2)
                right_start = start + int(nt_range / 2)
                right_end = end + int(nt_range/2)

                for right in range(left_end, i, -1):
                    if err_type == "del":
                        right_flank_erate = err_rates[right].fullDel
                    elif err_type == "ins":
                        right_flank_erate = err_rates[right].fullIns
                    else:
                        pass

                    if right_flank_erate <= 0:
                        left_start -= 1
                    else:
                        break

                for left in range(right_start, i):
                    if err_type == "del":
                        left_flank_erate = err_rates[left].fullDel
                    elif err_type == "ins":
                        left_flank_erate = err_rates[left].fullIns
                    else:
                        pass

                    if left_flank_erate <= 0:
                        right_end += 1
                    else:
                        break

                left_end = right + 1
                right_start = left
                left_start += 1
                right_end -= 1

                # get the different window erates
                left_window_erate = get_best_indel_erate_in_range(
                        err_rates, err_type, err_len, left_start, left_end,
                        (not nt_range), winshift_nlen_thres, eused2, find=i)

                if left_window_erate >= 0:
                    right_window_erate = get_best_indel_erate_in_range(
                        err_rates, err_type, err_len, right_start, right_end,
                        (not nt_range), winshift_nlen_thres, eused3, find=i)
                else:
                    right_window_erate = 0

                if abs(left_window_erate) > indel_err \
                        or abs(right_window_erate) > indel_err:
                    if abs(left_window_erate) > abs(right_window_erate):
                        # set left window as window used
                        indel_err = left_window_erate
                        start = left_start
                        end = left_end
                        del erates_used[:]
                        for e in eused2:
                            erates_used.append(e)
                    elif abs(right_window_erate) > abs(left_window_erate):
                        # set right window as window used
                        indel_err = right_window_erate
                        start = right_start
                        end = right_end
                        del erates_used[:]
                        for e in eused3:
                            erates_used.append(e)
                    else:
                        # if equal take the leftmost one
                        indel_err = left_window_erate
                        start = left_start
                        end = left_end
                        del erates_used[:]
                        for e in eused2:
                            erates_used.append(e)

                    # for e length changes
                    if indel_err < 0:
                        err_len -= (-1)**err_len
                        indel_err = abs(indel_err)
                    del eused2[:]
                    del eused3[:]
                else:
                    pass
            else:
                pass

            if indel_err > 1.2:
                del erates_used[:]
                continue
            elif indel_err > thres:
                # determine error properties
                pos, depth, ref_nt, rep_nt = determine_range_error_info(
                    err_rates, err_type, err_len, start, end, i, (not nt_range),
                    adjust=False, existing_errs=existing_errs)

                if pos == -1:
                    # Ideal position is already considered as an error
                    pos, depth, ref_nt, rep_nt = determine_range_error_info(
                        err_rates, err_type, err_len, start, end, i,
                        (not nt_range), adjust=False,
                        existing_errs=existing_errs, check_diff_nlen=True)
                    if pos == -1:
                        print(("BUG: error rate determination is finding an "
                               "error but no suitable candidate found near "
                               "position {} within contig {}. DEBUG.")
                              .format(err_rates[i].pos, err_rates[i].contig),
                              file=sys.stderr)
                        continue

                set_erate_to_zero(err_rates, err_type, err_len, erates_used, i)
            else:
                del erates_used[:]
                continue
        else:
            continue

        adjusted_i = i - (err_rates[i].pos-pos)
        if pos in existing_errs:
            print(("BUG: Collision. Two errors were determined to be on contig "
                   "{} at position {}").format(err_rates[pos].contig, pos),
                  file=sys.stderr)
        else:
            existing_errs.append(pos)

        err_list.append([err_type,
                        err_rates[i].contig,
                        pos, depth,
                        err_len,
                        ref_nt,
                        rep_nt,
                        indel_err,
                        adjusted_i])


def get_initial_errors(ref_mpile, low_mpile):
    """Detect errors from reference genome alignments

    This function uses calculated error rates to determine all possible
    reference based error through multiple passes through the list of error
    rates. The first two passes are mandatory with set window sizes of 1 and
    homopolymer length. The lsat n passes start at a window size of 11 and
    increases by 10 each successive pass. All later passes consider window
    shifting and different length search parameters set.

    Parameters
    ----------
    ref_mpile : str
        Path to the primary reference genome alignments
    low_mpile : str
        Path to the supplementary reference genome alignments

    Returns
    -------
    List of lists.
        All detected errors found with its information in a nested list

    """

    err_rates = get_mpile_error_rates(ref_mpile, low_mpile)

    # initialization
    found_err = []
    found_err_pos = []
    num_of_passes = 2 + args.passes

    for n in range(-1, num_of_passes):
        nt_range = n * 5
        if n > 2:
            find_errors(err_rates, found_err, nt_range, low_mpile,
                        winshift_nlen=-1, existing_errs=found_err_pos)
        elif n > 0:
            find_errors(err_rates, found_err, nt_range, low_mpile,
                        winshift_nlen=((n - 1) % 2),
                        existing_errs=found_err_pos)
        else:
            find_errors(err_rates, found_err, nt_range, low_mpile,
                        winshift_nlen=0, existing_errs=found_err_pos)

    found_err.sort(key=itemgetter(8, 0))
    return found_err


def get_region_error_rates(mpile, error, pile_info, nt_range):
    """Calculate error rates for each position in the search window

    Parameters
    ----------
    mpile : open file
        Adjustment alignment file
    error : list
        Single error to be adjusted
    pile_info : class
        For information storage between positions
    nt_range : int
        Search window size

    Returns
    -------
    list
        Calculated error rates as Error classes

    """

    # single error information for readability
    err_type = error[0]
    err_contig = error[1]
    err_pos = error[2]

    err_rates = []

    # adjust start point depending on error type
    if err_type == "del" and len(mpile) >= ((nt_range * 2) + 1):
        pile_info.set_next_del_to_zero()
        parse_del_only(pile_info, mpile[0].split()[4])
        start = 1
    else:
        start = 0
    end = len(mpile)

    for index in range(start, end):
        mpile_line = mpile[index].split()
        depth = int(mpile_line[3])

        # Error check
        contig = mpile_line[0]
        pos = int(mpile_line[1]) - 1
        if contig != err_contig \
                or not ((err_pos - end) <= pos <= (err_pos + end)):
            print(("WARNING: Either contig does not match or position {} is "
                   "not within the range ({} to {}) being assessed. Double "
                   "check that all mpileup files are the same length and order")
                  .format(pos, start, end), file=sys.stderr)

        # initialize error rate information
        pos_erate = Error(contig, pos, mpile_line[2], depth)

        if depth < _depth_thres:
            pos_erate.set_low_depth()
            pos_erate.readEnd = parse_read_end_only(mpile_line[4])
            err_rates.append(pos_erate)

            pile_info.set_next_del_to_zero()
            parse_del_only(pile_info, mpile_line[4])
            if index == (end - 1):
                fill_hp_info(err_rates, True)
            else:
                fill_hp_info(err_rates)
            continue
        else:
            pass

        try:
            pos_erate.next = mpile[index + 1].split()[2]
        except IndexError:
            pos_erate.next = "-"

        # calculate error rate
        parse_pile(pos_erate, mpile_line[4], pile_info, get_only=error[0])

        # reset
        pile_info.reset_class()
        pile_info.set_next_del_to_zero()
        del pos_erate.next

        err_rates.append(pos_erate)

        # NOTE: hp info on the edge of the region will likely be incorrect,
        # but unlikely to use this edge info.
        # WARNING: No adjustments are made to account for this edge info
        if end == 1:
            return err_rates
        elif index == (end - 1):
            fill_hp_info(err_rates, True)
        else:
            fill_hp_info(err_rates)

        # set up deletion information for next position
        parse_del_only(pile_info, mpile_line[4])

    return err_rates


def get_adjustment_candidates(error, err_rates, last_pos):
    """Finds all adjustment candidates for the specified errors

    Using previously calculated error rates for the search window, this
    function combines error rates per homopolymer to determine the error
    rate for each homopolymer. Errors in the middle of a homopolymer are
    considered as uniq error rate positions. Both 1 and 2 nt errors are
    considered. If there is no a candidate of the same error length,
    multiple candidates of a shorter error length are returned.

    Parameters
    ----------
    error :  list
        A single error and its information
    err_rates : list
        Calculated error rates
    last_pos : int
        Last position so errors outsides of the search are not considered

    Returns
    -------
    List of lists.
        All adjustment candidates of n length

    """
    # initialization
    err_type = error[0]
    err_len = error[4]
    err_candidates = [[], []]
    start = -1

    for index in range(len(err_rates)):
        if start == err_rates[index].hp_start:
            continue
        start = err_rates[index].hp_start
        end = start + err_rates[start].hp_len
        top_non_hp_err1 = []
        top_non_hp_err2 = []

        # get best indel error rate for both lengths
        err1 = get_best_indel_erate_in_range(err_rates, err_type, 1, start, end,
                                             True, uniq_erates=top_non_hp_err1,
                                             adjust=True)
        err2 = get_best_indel_erate_in_range(err_rates, err_type, 2, start, end,
                                             True, uniq_erates=top_non_hp_err2,
                                             adjust=True)

        # Error length 1
        if top_non_hp_err1 \
                and top_non_hp_err1[1] >= err1 \
                and top_non_hp_err1[1] >= _adjust_thres:
            err_pos, depth, err_ref, err_rep = determine_range_error_info(
                err_rates, err_type, 1, top_non_hp_err1[0],
                top_non_hp_err1[0] + 1, hp=True, adjust=True)
            if err_pos == -1:
                # found no position candidate
                err_pos, depth, err_ref, err_rep = determine_range_error_info(
                    err_rates, err_type, 2, top_non_hp_err1[0],
                    top_non_hp_err1[0] + 1, hp=True, adjust=True)

            # save as a candidate if it wasn't previous used as a candidate.
            if err_pos >= 0 not in last_pos:
                err_candidates[0].append([err_type, error[1], err_pos, depth, 1,
                                          err_ref[0], err_rep[0],
                                          top_non_hp_err1[1]])
        if err1 >= _adjust_thres:
            err_pos, depth, err_ref, err_rep = determine_range_error_info(
                err_rates, err_type, 1, start, end, hp=True, adjust=True)
            if err_pos == -1:
                err_pos, depth, err_ref, err_rep = determine_range_error_info(
                    err_rates, err_type, 2, start, end, hp=True, adjust=True)
            if err_pos >= 0 and err_pos not in last_pos:
                err_candidates[0].append([err_type, error[1], err_pos, depth, 1,
                                          err_ref[0], err_rep[0], err1])

        # Error length 2
        if top_non_hp_err2 and top_non_hp_err2[1] >= err2 \
                and top_non_hp_err2[1] >= (_adjust_thres - 0.1):
            err_pos, depth, err_ref, err_rep = determine_range_error_info(
                err_rates, err_type, 2, top_non_hp_err2[0],
                top_non_hp_err2[0] + 1, hp=True, adjust=True)
            if err_pos >= 0 and err_pos not in last_pos:
                err_candidates[1].append([err_type, error[1], err_pos, depth, 2,
                                          err_ref, err_rep, top_non_hp_err2[1]])
        if err2 >= (_adjust_thres - 0.1):
            err_pos, depth, err_ref, err_rep = determine_range_error_info(
                err_rates, err_type, 2, start, end, hp=True, adjust=True)
            if err_pos >= 0 and err_pos not in last_pos:
                err_candidates[1].append([err_type, error[1], err_pos, depth, 2,
                                          err_ref, err_rep, err2])

    # return only the useful candidates
    if err_len == 1 and err_candidates[0]:
        return err_candidates[0]
    elif err_len == 2 and err_candidates[1]:
        return err_candidates[1]
    elif err_len == 2 and len(err_candidates[0]) >= 2:
        # Splitting error into multiple shorter candidates
        return err_candidates[0]
    else:
        return []


def adjust_error(error, can, len_change=False):
    """Replaces the current error with the adjustment candidate

    Parameters
    ----------
    error : list
        A single error and its information
    can : list
        Adjustment candidate and its information
    len_change : boolean
        Indicates if an error is being spllit. Ex. 2 nt error into two 1 nt
        errors

    Returns
    -------
    None
        Adjustments are directly modified in the single error

    """

    err_type = error[0]
    if len_change:
        error[4] = can[4]
    if error[2] == can[2]:
        # same position
        if error[6] != can[6]:
            error[6] = can[6]
    else:
        # adjust position
        error[2] = can[2]
        if err_type == "del":
            error[5] = can[5]
        elif err_type == "ins":
            error[6] = can[6]
        else:
            print(("BUG: Could not determine error type in function adjust_"
                   "error. DEBUG."), file=sys.stderr)

    # add support number
    error[7] = ";".join([str(error[7]), str(can[7])])

    return


def _load_read_file_offsets(file):
    """Determine and store bit offset for each line in a file

    This function the exact bit where each line starts in a file. These bits
    allows for jumping around in a file since line lengths are not equal.

    Credit: Adam Rosenfield from StackOverflow

    Parameters
    ----------
    file : open file
        File to read to determine the bit offsets

    Returns
    -------
    List
        Bit offsets as integers

    """

    line_offset = []
    offset = 0
    try:
        for line in file:
            line_offset.append(offset)
            offset += len(line)
    except UnicodeDecodeError:
        sys.exit("ERROR: Mpileup file is corrupted. Re-try with another "
                 "mpileup file.")

    return line_offset


def _retrieve_mpile_chunk(file, offset, region_size, chunk):
    """Load a subset of alignments from an mpileup file

    Takes a file's bit offset to retrieve n (region_chunk) lines beyond the
    initial line set by the offset.

    Parameters
    ----------
    file : open file
        Alignment file.
    offset : int
        Bit offset to indicate the position of the first line within
        the file
    region_size : int
        Number of lines to retrieve after the first line
    chunk : list
        Old chunk to overwrite

    Returns
    -------
    List
        A list of alignment information. Each element is a single line from
        the alignment file

    """
    file.seek(offset)
    # reset chunk
    del chunk[:]

    # receive new chunk       
    for n in range(region_size):
        chunk.append(file.readline())
        if chunk[-1] == "":
            del chunk[-1]
            break
    return chunk


def adjust_initial_error(reads_mpile, err_list, adjust_sub=True):
    """Adjusts the errors detected from reference data using genome info

    Imports the files to be used for adjustments. For each file, the code
    loops through error list pulling the appropriate lines based on file
    index. Files with multiple contigs need to ensure that contig and
    position indices matches the reference alignment file (same order).
    Error rates are calculated within the search window. Positions greater
    than the adjustment threshold are considered candidates with the best
    candidate replacing the respective error.

    Parameters
    ----------
    reads_mpile : str
        Path to adjustment files
    err_list : list
        Detected errors from reference data
    adjust_sub : boolean
        Consider substitutions for current adjustment file.

    Returns
    -------
    None
        Adjustments are directly modified in passed error list

    """

    # TODO: change from global with a settings class
    global _adjustment_count
    global _supported_sub

    # Retrieve bit offsets for each line in adjustment file
    reads = open(reads_mpile)
    offsets = _load_read_file_offsets(reads)

    # Initialization
    mpile_chunk = []
    container = MpileInfo()
    extra = []                                      # To store split errors
    last_adjust_pos = [-1]

    for error in err_list:
        # SUBSTITUTIONS
        if error[0] == "sub":
            if adjust_sub:
                _retrieve_mpile_chunk(reads, offsets[error[8]], 1,
                                      mpile_chunk)
                err_rates = get_region_error_rates(mpile_chunk, error,
                                                   container, 1)

                # retrieve nt to consider for substitution
                try:
                    poss_subs = list(err_rates[0].sub.nt.intersection(error[6]))
                except AttributeError:
                    poss_subs = []

                # decide best substitution error rate
                best_sub_erate = 0
                for nt in poss_subs:
                    try:
                        sub_erate = err_rates[0].subAllErr[nt]
                    except AttributeError:
                        sub_erate = ({"A": 0, "C": 0, "G": 0, "T": 0})[nt]
                    if sub_erate > best_sub_erate:
                        best_sub_erate = sub_erate
                        best_nt = nt
                else:
                    pass

                # evaluate adjustment
                if best_sub_erate > 0.1:
                    _supported_sub += 1
                    error[6] = best_nt
                    error[7] = ";".join([str(error[7]), str(best_sub_erate)])
                    del error[-1]
                else:
                    # delete substitution error due to little evidence
                    del error[:]
            continue

        # INDELS

        nt_range = 10           # window size

        # Calculate error rates in the window and find candidates
        try:
            if error[0] == "del" and (error[8] - 1) > 0:
                _retrieve_mpile_chunk(reads, offsets[(error[8] - nt_range-1)],
                                      (nt_range*2 + 1), mpile_chunk)
            else:
                _retrieve_mpile_chunk(reads, offsets[(error[8] - nt_range)],
                                      (nt_range*2), mpile_chunk)
        except IndexError:
            print(("BUG: Could not retrieve relevant mpileup chunk between "
                   "indices {} and {}.").format((error[8] - nt_range - 1),
                                                (nt_range * 2 + 1)),
                  file=sys.stderr)
        err_rates = get_region_error_rates(mpile_chunk, error, container,
                                           nt_range)
        adj_can = get_adjustment_candidates(error, err_rates, last_adjust_pos)

        # Evaluate all possible adjustment candidates
        if len(adj_can) == 1:
            # only one candidate
            _adjustment_count += 1
            adjust_error(error, adj_can[0])
            last_adjust_pos = [adj_can[0][2]]
        elif len(adj_can) > 1:
            _adjustment_count += 1
            # more than one candidate. Select the best one
            # TODO: undo mess
            hp_start = err_rates[nt_range].hp_start
            for index in range(len(adj_can)):
                if (adj_can[index][2] == error[2]
                    or err_rates[hp_start].pos == adj_can[index][2]) \
                        and (adj_can[index][6][0] == error[6][0]
                             or adj_can[index][6][-1] == error[6][-1]):
                    best_can = adj_can[index]
                    break
                else:
                    best_can = []

            if not best_can:
                best_can = max(adj_can, key=lambda k: k[7])
            last_adjust_pos = [best_can[2]]

            # select a second best candidate if splitting original error
            if error[4] != adj_can[0][4]:
                # first adjustment
                adjust_error(error, best_can, True)

                # Determine next candidate
                to_del = next(index for index in range(len(adj_can))
                              if adj_can[index] == best_can)
                del adj_can[to_del]
                best_can = max(adj_can, key=lambda k: k[7])
                last_adjust_pos.append(best_can[2])

                # add index for further adjustment if needed
                if not adjust_sub:
                    pos_diff = error[2] - best_can[2]
                    best_can.append((error[8] - pos_diff))

                # add second candidate to extras. Resolve last
                extra.append(best_can)
            else:
                adjust_error(error, best_can)
        else:
            error[7] = ";".join([str(error[7]), "0"])
            pass

        # Delete indel's index information if not needed
        if adjust_sub:
            del error[-1]

    err_list += extra
    del extra[:]

    # Remove unsupported substitutions
    if adjust_sub:
        del_count = 0
        for i in range(len(err_list) - 1, -1, -1):
            if err_list[i]:
                continue
            else:
                del_count += 1
                del err_list[i]
        print("{} substitutions had little read support and were deleted"
              .format(del_count))

    print("{} errors have been adjusted with {} file"
          .format(_adjustment_count, reads_mpile))
    _adjustment_count = 0
    reads.close()

    return


def get_flank_nt_pos(gen, pos):
    """Retrieve the upstream and downstream positions flanking a homopolymer

    This function determines the edges of homopolymer of the specified
    position and returns the positions of the nucleotides flanking it,

    Parameters
    ----------
    gen : dict
        Pre-corrected draft assembly
    pos : int
        Specified position

    Returns
    -------
    int
        Upstream position
    int
        Downstream position

    """
    index = pos
    pre_pos, suc_pos = pos, pos
    cur_nt = gen[pos]
    direction = 1

    # moves in positive or negative direction until reaching non-homopolymer
    while (len(gen) + 1) > index > -1:
        if gen[index] == cur_nt:
            index += direction
        else:
            if index > pos:
                # Downstream position was found. Store position.
                # Reset position to search upstream.
                suc_pos = index
                index = pos - 1
                direction = -1
            elif index < pos:
                pre_pos = index
                break         # get out of loop
            else:
                print("Did not find any position. Check your results")
                break

    return pre_pos, suc_pos


# TODO: adjust for multi-contig assemblies
def get_hp(gen, pos, direction):
    """Retrieve the sequence and length of homopolymer

    Parameters
    ----------
    gen : dict
        Pre-corrected draft assembly
    pos : int
        Position in the draft assembly
    direction : int
        Determine to look upstream or downstream of the position

    Returns
    -------
    char
        Homopolymer nt
    int
        Homopolymer length

    """
    index = pos
    cur_nt = gen[pos]

    while (len(gen) + 1) > index > -1:
        if gen[index] == cur_nt:
            index += direction
        else:
            index += (direction * -1)
            break

    # subset homopolymer from draft sequence
    if index < pos:
        hp = gen[index:(pos+1)]
    elif index > pos:
        hp = gen[pos:index+1]
    else:
        hp = gen[pos]

    return hp[0], len(hp)


def reverse_comp(seq):
    """Reverse complements a sequence

    Parameters
    ----------
    seq : str
        String of the sequence to be reverse complemented

    Returns
    -------
    String
        reverse complement of the input sequence

    """
    return seq.translate(seq.maketrans("ACGTacgt", "TGCATGCA"))


def get_diagnostic_hp_data(out, draft, errors):
    """Retrieve information on each error

    Information retrieve include upstream, current, and

    Parameters
    ----------
    out : str
        Output filename as a string
    draft : dict
        Pre-corrected draft assembly
    errors : list
        Detected error list

    Returns
    -------
    None

    """

    dia_data = open("{}.diagnostic.txt".format(out), "w")
    dia_data.write(("Contig\tPosition\tCurrentHp\tCurrentHpLen\tPreHp\t"
                    "PreHpLen\tPostHp\tPostHpLen\tFor5mer\tRev5mer\t"
                    "CreateHp\tAdjustHp\tBreakHp\n"))
    for contig, err_list in errors.items():
        contig_seq = draft[contig]

        for err in err_list:
            err_type = err[0]
            err_pos = int(err[2])
            orig_pos = err_pos
            err_len = err[4]

            # change the central position for insertion where inserting nt
            # is the same as the next homopolymer
            if err_type == "ins":
                err_nt = err[6]
                try:
                    if err_nt[0] != contig_seq[err_pos] \
                            and err_nt[-1] == contig_seq[(err_pos + 1)]:
                        err_pos += 1
                except IndexError:
                    pass
            elif err_type == "del":
                err_nt = err[5]
            else:
                err_nt = err[6]

            pre_nt_pos, suc_nt_pos = get_flank_nt_pos(contig_seq, err_pos)

            cur_nt = contig_seq[err_pos]
            cur_nt_len = suc_nt_pos - pre_nt_pos - 1
            pre_nt, pre_nt_len = get_hp(contig_seq, pre_nt_pos, -1)
            suc_nt, suc_nt_len = get_hp(contig_seq, suc_nt_pos, 1)

            # adjust the error position to last deletion if deletion greater
            # than one
            if err_type == "del" and err_len > 1:
                for_pos = err_pos
                rev_pos = err_pos + err_len
            elif err_type == "ins":
                for_pos = err_pos
                rev_pos = for_pos
            else:
                for_pos = err_pos
                rev_pos = err_pos + 1

            # flanking forward and reverse 5-mer
            try:
                forward_5mer = contig_seq[(err_pos - 5):for_pos]
            except IndexError:
                forward_5mer = contig_seq[0:for_pos]

            try:
                reverse_5mer = reverse_comp(contig_seq[rev_pos:(rev_pos + 5)])
            except IndexError:
                reverse_5mer = reverse_comp(contig_seq[rev_pos:len(contig_seq)])

            # break, adjust, create homopolymer information            
            if err_type == "ins":
                pre = contig_seq[orig_pos]
                post = contig_seq[orig_pos + 1]
            else:
                pre = contig_seq[err_pos - 1]
                post = contig_seq[rev_pos]

            create_hp = False
            adjust_hp = False
            break_hp = False                # when no flanking nt is similar
            if err_type == "ins":
                if err_nt[0] == pre:
                    if pre == cur_nt:
                        if cur_nt_len == 1:
                            create_hp = True
                        else:
                            adjust_hp = True
                    elif pre == pre_nt:
                        if pre_nt_len == 1:
                            create_hp = True
                        else:
                            adjust_hp = True
                    else:
                        pass

                if err_nt[-1] == post:
                    if post == cur_nt:
                        if cur_nt_len == 1:
                            create_hp = True
                        else:
                            adjust_hp = True
                    elif post == suc_nt:
                        if suc_nt_len == 1:
                            create_hp = True
                        else:
                            adjust_hp = True
                    else:
                        pass

                if pre == post:
                    break_hp = True
                else:
                    pass
            elif err_type == "del":
                if pre == post:
                    create_hp = True
                    if pre == cur_nt and cur_nt_len > 1:
                        adjust_hp = True
                    elif pre == pre_nt and pre_nt_len > 1:
                        adjust_hp = True
                    elif post == suc_nt and suc_nt_len > 1:
                        adjust_hp = True
                    else:
                        pass
                else:
                    if err_len == 1:
                        if cur_nt_len == 2:
                            break_hp = True
                        elif cur_nt_len == 1:
                            pass
                        else:
                            adjust_hp = True
                    else:
                        if err_nt[-1] == cur_nt and (cur_nt_len <= 3):
                            break_hp = True
                        elif err_nt[-1] == cur_nt and cur_nt_len > 3:
                            adjust_hp = True
                        elif err_nt[-1] == suc_nt:
                            if suc_nt_len == 2 or cur_nt_len == 2:
                                break_hp = True
                            if suc_nt_len > 2 or cur_nt_len > 2:
                                adjust_hp = True
                        else:
                            # end deletion is next hp
                            if cur_nt_len == 2:
                                break_hp = True
                            elif cur_nt_len == 1:
                                pass
                            else:
                                adjust_hp = True

                            if suc_nt_len == 2:
                                break_hp = True
                            elif suc_nt_len == 1:
                                pass
                            else:
                                adjust_hp = True
            else:
                # assume sub
                if err_nt == pre_nt \
                        and err_nt == suc_nt \
                        and err_nt == pre \
                        and err_nt == post:
                    create_hp = True
                    if pre_nt_len > 1 or suc_nt_len > 1:
                        adjust_hp = True
                elif pre == pre_nt and err_nt == pre_nt:
                    if pre_nt_len == 1:
                        create_hp = True
                    else:
                        adjust_hp = True

                    if cur_nt_len == 2:
                        break_hp = True
                    elif cur_nt_len > 2:
                        adjust_hp = True

                elif post == suc_nt and err_nt == suc_nt:
                    if suc_nt_len == 1:
                        create_hp = True
                    else:
                        adjust_hp = True

                    if cur_nt_len == 2:
                        break_hp = True
                    elif cur_nt_len > 2:
                        adjust_hp = True
                elif pre == cur_nt and post == cur_nt:
                    break_hp = True
                elif pre == cur_nt or post == cur_nt:
                    if cur_nt_len == 2:
                        break_hp = True
                    else:
                        adjust_hp = True
                else:
                    pass

            lst_info = [contig, err[2], cur_nt, cur_nt_len, pre_nt[0],
                        pre_nt_len, suc_nt[0], suc_nt_len, forward_5mer,
                        reverse_5mer, create_hp, adjust_hp, break_hp]
            dia_data.write("%s\n" % "\t".join(map(str, lst_info)))

    return


def split_errors_by_contig(errors, keep=False):
    """Split the error list by contig

    Contigs are split by name at position one. Error list must be ordered
    to group all errors for individual contigs together.

    Parameters
    ----------
    errors : list
        Detected error list
    keep : boolean
        Decide whether to keep original error list

    Returns
    -------
    Dictionary
        Split error list with contig name as the key

    """
    split_err = {}

    if keep:
        err_list = errors[:]
    else:
        err_list = errors

    # Splitting
    while len(err_list) != 0:
        current_contig = err_list[0][1]
        try:
            next_contig_index = next(index for index in range(len(err_list))
                                     if err_list[index][1] != current_contig)
            split_err[current_contig] = err_list[0:next_contig_index]
            del err_list[0:next_contig_index]
        except StopIteration:
            split_err[current_contig] = err_list[0:len(err_list)]
            del err_list[0:len(err_list)]
            break

    return split_err


def correct_genome(draft, errors):
    """Corrects the loaded draft assembly based on passed error list

    Parameters
    ----------
    draft : dict
        Draft assembly; dictionary of sequences
    errors : dict
        A dictionary of error lists for individual contigs. Key
        is contig names; these names should match sequence headers in the draft
        assembly

    Returns
    -------
    None

    """
    correction_count = 0

    # for multiple contigs
    for contig, err_list in errors.items():
        contig_seq = list(draft[contig])

        for correction in err_list:
            pos = correction[2]
            err_type = correction[0]
            err_len = correction[4]

            # actual correction
            if err_type == "sub":
                if len(contig_seq[pos]) > 1:
                    temp_pos_nt = list(contig_seq[pos])
                    temp_pos_nt[0] = correction[6]
                    contig_seq[pos] = "".join(temp_pos_nt)
                else:
                    contig_seq[pos] = correction[6]
            elif err_type == "del":
                if err_len == 1:
                    contig_seq[pos] = ""
                else:
                    contig_seq[pos] = ""
                    contig_seq[pos+1] = ""
            elif err_type == "ins":
                contig_seq[pos] += correction[6]
            else:
                print(("BUG: Could not determine error type at position {} "
                       "for sequence correction").format(pos), file=sys.stderr)

            correction_count += 1

        # done for this contig
        draft[contig] = "".join(contig_seq)

    print("Made {} corrections!".format(correction_count))


def _write_genome_fasta(fasta_name, write_type, gen):
    """Output corrected assembly as a fasta file

    Parameters
    ----------
    fasta_name : str
        Name of output file
    write_type : char
        Output type; write or append
    gen : dict
        Corrected assembly as a dictionary

    Returns
    -------
    None

    """
    with open(fasta_name, write_type) as fasta:
        if type(gen).__name__ == "dict":
            for contig, seq in gen.items():
                fasta.write("%s\n" % ">{}".format(contig))
                if type(seq).__name__ == "list":
                    fasta.write("%s\n" % "".join(seq))
                elif type(seq).__name__ == "str":
                    fasta.write("%s\n" % seq)
        else:
            print(("BUG: Input type is not a dictionary. Could not output "
                   "genome as fasta"), file=sys.stderr)


def main():
    """Runs Castor modules: Error detection, Error

    Returns
    -------
    None

    """

    # Error detection
    if args.errors != "" or _indel_thres >= 2:
        errors = _load_prev_error(args.errors)
    elif args.errors == "" and _indel_thres >= 2:
        sys.exit("Please enter an error file.")
    else:
        errors = get_initial_errors(args.ref, args.low)
        _print_info("{}.initial.err".format(args.out), "w", errors)

    # Error adjustment
    if args.adjust < 2:
        for file_num in range(len(args.reads)):
            errors.sort(key=itemgetter(1, 2))
            if file_num == (len(args.reads) - 1):
                adjust_initial_error(args.reads[file_num], errors,
                                     adjust_sub=True)
            else:
                adjust_initial_error(args.reads[file_num], errors,
                                     adjust_sub=False)
        _print_info("{}.adjusted.02.err".format(args.out), "w", errors)

    # Resort the errors before correction
    errors.sort(key=itemgetter(1, 2))

    if args.module in ("Correct", "correct", "All", "all"):
        draft = _read_draft(args.draft)
        split_err = split_errors_by_contig(errors)

        # TODO: set as extraInfo
        # write some diagnostic data for hp info of the reads
        # get_diagnostic_hp_data(args.out, draft, split_err)

        correct_genome(draft, split_err)
        _write_genome_fasta("{}.fa".format(args.out), "w", draft)


if __name__ == "__main__":
    args = _parse_arguments()

    # Global thresholds to maintain
    # TODO: update with proper settings
    if args.module in ("Detect", "detect"):
        _indel_thres = args.indel
        _sub_thres = args.sub
        _adjust_thres = 2
    elif args.module in ("Adjust", "adjust"):
        _indel_thres = 2
        _sub_thres = 2
        _adjust_thres = args.adjust
    elif args.module in ("Correct", "correct"):
        _indel_thres = 2
        _sub_thres = 2
        _adjust_thres = 2
    else:
        _indel_thres = args.indel
        _sub_thres = args.sub
        _adjust_thres = args.adjust

    _depth_thres = args.depth
    _low_depth_count = 0
    _adjustment_count = 0
    _supported_sub = 0

    main()

    # Output last information
    print("{} positions were considered low depth and not evaluated"
          .format(_low_depth_count))
    print("{} substitutions have supporting read data".format(_supported_sub))
    print("DONE!")




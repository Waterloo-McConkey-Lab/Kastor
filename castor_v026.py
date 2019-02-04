#!/usr/bin/python3

# Python script to use reference genomes and reads aligned to draft genome for correction. Input are sam pileup files.
# Note: key 
import sys
import argparse 
import re
from Bio import SeqIO
from Bio.Seq import Seq
from operator import itemgetter
from collections import Counter
from rf_class import *

#REMOVE
import pdb
import gc

# read input arguments
def _parse_Arguments():
    # Input arguments
    parser = argparse.ArgumentParser(
            usage = "ref_correct.py [OPTIONS] <draft_genome.fa> <ref_genome.pileup> <mapped_reads.pileup>",
            description = "Correct a draft genome using reference genomes and reads mapped to the draft. version 0.2")
    parser.add_argument(
            "draft",
            nargs = "?",
            type = str,
            help = "Draft genome FASTA input file")
    parser.add_argument(
            "ref",
            nargs = "?",
            type = str,
            help = "Reference genomes mapped to the draft genome in samtools mpileup format file")
    parser.add_argument(
            "reads",
            nargs = "+",
            type = str,
            help = """Mpileup alignment information used to adjust found errors. Recommended to use
                      read information for adjustment. More than one mpileup file can be used for 
                      multiple adjustment attempts. One caveat: substitution adjustments only use the
                      last alignment file so substitution data must be contained in that file otherwise
                      all substitution would be considered unsupported and deleted.
                      Ex. Fastq reads mapped to the draft genome in samtools mpileup format file""")
    parser.add_argument(
            "--low", "-l",
            nargs = "?",
            type = str,
            default = "",
            help = """Optional: File path to mpileup file containing greater and more dstant reference
                      genome used to correct for low depth regions. --close mpileup set will be used to 
                      adjust threshold. Note: imbalanced data will bias correction towards the 
                      phylogroup containing the greatest number of genomes in the dataset.""")
    parser.add_argument(
            "--sub", "-s",
            nargs = "?",
            type = float, 
            default = 1.0,
            help = """Threshold to correct substititions. 0 to 1. Use >1 to turn off substition 
                      correction. [Default = 0.8]""")
    parser.add_argument(
            "--indel", "-i",
            nargs = "?",
            type = float,
            default = 0.85,
            help = """Threshold to correct indel errors. 0 to 1. Use > 1 to turn off indel correction. 
                      [Default = 0.85]""")
    parser.add_argument(
            "--adjust", "-a",
            nargs = "?",
            type = float,
            default = 0.25,
            help = """Threshold to used to consider a candidate for error adjustment. Use >2 to turn off 
                      error adjustments [Default = 0.25]""")
    parser.add_argument(
            "--low_threshold", 
            nargs = "?",
            type = float,
            default = 0.99,
            help = """Threshold used for optional low depth mpileup file. [Default = 0.99]""")
    parser.add_argument(
            "--nRef", "-n",
            nargs = "?",
            type = float,
            default = 11,
            help = """Expected number of genomes to match each position. Used in conjunction with --relax.
                      Relax indel threshold by r amount when depth is above expected genome number 
                      [Default = 11]""")
    parser.add_argument(
            "--depth", "-d",
            nargs = "?",
            type = int,
            default = 5,
            help = "Minimum depth required to consider the mpileup information. [Default = 5]")
    parser.add_argument(
            "--passes", "-p",
            nargs = "?",
            type = int,
            default = 6,
            help = "Number of passes through error rates to determine putative errors [Default = 6]")    
    parser.add_argument(
            "--out", "-o",
            nargs = "?",
            type = str,
            default = "out",
            help = "Prefix for output files [Default: out]")
    parser.add_argument(
            "--extraInfo", "-e",
            action = "store_true",
            help = """Print out extra information. Includes homopolymer information at error site, 
                      regions skipped because of low depth, initial error calls, and errors called at 
                      regions of low depth.""")
    parser.add_argument(
            "--messy", "-m",
            action = "store_true",
            help = """Turn on to check if surrounding the area surrounding a putative error have high
                      incidences of indels and thus 'messy' via alignment. Putative errors in these regions 
                      will be ignored""")    
    parser.add_argument(
            "--verbal", "-v",
            action = "store_true",
            help = "Turn on informational/debugging text")
    parser.add_argument(
            "--errors",
            nargs = "?",
            type = str,
            default = "",
            help = """Bypass calculation of errors and correct with input error file. Enter a input.err 
                      file. Overflow file (.overflow) must be in same directory.""")    

    return parser.parse_args()

# read the draft genome
def _read_draft():
    draft = SeqIO.parse(open(args.draft), "fasta")
    d_gen = {}
    for contig in draft:
        head, seq = contig.id, str(contig.seq)
        print("Reading contig: {}".format(head))
        d_gen[head] = seq
        #d_gen.append(head)
        #d_gen.append(seq)
    return d_gen

# takes a pre-computed error file and loads it into memory
def _load_prev_error(err_file):
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
                err_loaded.append([err_can, contig, pos, depth, err_len, ref_nt, replace_nt, err_depth, err_index])
            except IndexError:
                err_loaded.append([err_can, contig, pos, depth, err_len, ref_nt, replace_nt, err_depth])

    return err_loaded

# load mpileup file as list of mpileup lines
def _load_full_mpile_info(mpile_file):
    mpile = []

    if mpile_file == "":
        return []           # empty path
    else:        
        # get mpile info from open file
        with open(mpile_file, "r") as f:
            mpile = f.readlines()

        return mpile

# print dictionary or list info given that there are no custom data structures
def _print_info(out_file, write_type, info):
    with open(out_file, write_type) as out:
        if len(info) < 1: return

        #assumes dictionary values are lists
        if type(info).__name__ == "dict":
                for key, value in info.items():
                    out.write("%s\n" % "\t".join(map(str, value)))
        elif type(info).__name__ == "list":
            if  type(info[0]).__name__ == "list":
                for value in info:
                    out.write("%s\n" % "\t".join(map(str, value)))
            else:
                for value in info:
                    out.write(str(value) + "\n")

# print err_rates. Assumes a list of erates
def _print_erates(out_file, write_type, err_rates, start, end):
    # bin_file = open("bin.{}".format(out_file), write_type)
    with open(out_file, write_type) as out:
        if len(err_rates) < 1: return

        out.write("%s\n" % "\t".join(["contig", "pos", "ref", "depth", "ferate", 
                                      "sub_ferate", "del_ferate", "ins_ferate",
                                      "sub_serate", "del_serate1", "del_serate2",
                                      "ins_serate1", "ins_serate2"]))
        if type(err_rates).__name__ == "list":
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
                             
        else:
            print("WARNING: error rates were not printed due to incorrect data structure passed.", file=sys.stderr)

    # bin_file.close()

# for indels greater than 9 since 10+ would be separated as two characters
# should return int of length size and end position of digit in list 
# (i.e., if digit is 15, should return int(15) and position )
def check_indel_len(char, start):
    index = start + 1
    while index < len(char):
        try: 
            int(char[index])
        except ValueError:
            indel_len = "".join(char[start:(index)])
            return int(indel_len), index
        index += 1
    # should exit while loop with ValueError

# mostly for ins 
# return the most common element in a list
def get_most_common_nt(lst, ref="-", next="-"):
    if len(lst) < 1: return 0

    # no longer use max(set(lst), key=lst.count) since ties are chosen at random
    # Limit: will add computation time
    poss_nt = Counter(lst).most_common()
    max_nt = 0
    best_rank = 99                  # rank the nt matches. 0 = best rank(closest to ref hp). 3 = random
    if ref == "-":
        # likely sub or not sure about ref so return first random best hit
        return poss_nt[0][0]
    else:
        # some info to decide between multiple top hits
        for nt in Counter(lst).most_common():
            # set max count
            if nt[1] < max_nt:
                break
            else:
                max_nt = nt[1]

            # multiple max
            if nt[0][0] == ref:
                if nt[0][-1] == ref and 0 < best_rank:
                    best_nt = nt[0]; best_rank = 0
                elif nt[0][-1] == next and 1 < best_rank:
                    best_nt = nt[0]; best_rank = 1
                else:
                    if 2 < best_rank:
                        best_nt = nt[0]; best_rank = 2
            elif nt[0][-1] == next and 2 < best_rank:
                # only match to next 
                best_nt = nt[0]; best_rank = 2
            elif 3 < best_rank:
                # everything else
                best_nt = nt[0]; best_rank = 3
            else:
                pass

        return best_nt

# subsets indels of n length from a list of all indel lengths
def pull_nlength_indels_info(indel_list, n, alignment):
    # return subsetted list and alignment positions of n length indels. 
    # Don't modify current list as it may be used in another subset
    sub_align = 1
    subset_nt = []

    # make sure that the number of indels matches the amount of alignments
    if (bin(alignment).count("1") - 1) != len(indel_list): 
        print("BUG: Number of indels does not match bit alignments in function pull_nlength_indels_info", file=sys.stderr)

    align_counter = alignment.bit_length() - 1
    # pull only errors that are n long
    for index in range(len(indel_list)):
        for offset in range(align_counter, 0, -1):
            mask = 1 << (offset - 1)
            sub_align <<=1
            if ((alignment & mask) >> (offset - 1)) == 1:
                alignment ^= mask
                align_counter = offset - 1 
                break

        err_length = len(indel_list[index])
        if err_length == n:
            subset_nt.append(indel_list[index])
            sub_align ^= 1

    # if needed, adjust bit length based on full alignment
    if alignment.bit_length() != sub_align.bit_length():
        sub_align <<= (alignment.bit_length() - sub_align.bit_length())

    return subset_nt, sub_align

# store error rates for each position depending on error type
# Note: nt reported is the most common error nt at that position. If tie, the first nt is taken. 
def get_error_rate(err_type, pos_erate, err_info=None, indel_alignment=0):
    depth = pos_erate.depth
    
    # total error rate
    if err_info is None:
        # nothing to parse
        f_erate = 1
    else:
        # parse the passed list
        f_erate = (len(err_info)/depth)

    if err_type == "total":
        # already recorded a blank or total error rate
        pos_erate.erate = 1 - f_erate
    elif err_type == "del":
        # if new reads were added this position, deletion calculated from last position 
        # is out of alignment. Adjustment of alignment length (depth) is needed
        bit_diff = pos_erate.readEnd.bit_length() - indel_alignment.bit_length()
        if bit_diff > 0:
            indel_alignment <<= bit_diff
        elif bit_diff < 0:
            print("BUG: Difference in alignment bit length in function get_error_rate", file=sys.stderr)
        else:
            pass

        # deletion error rate
        pos_erate.fullDel = f_erate
        sub_del, sub_align = pull_nlength_indels_info(err_info, 1, indel_alignment)
        if sub_del:
            pos_erate.del1 = ErrType(sub_del[0], len(sub_del)/depth, sub_align)
        sub_del, sub_align = pull_nlength_indels_info(err_info, 2, indel_alignment)
        if sub_del:
            pos_erate.del2 = ErrType(sub_del[0], len(sub_del)/depth, sub_align)
    elif err_type == "ins":
        # insertion error rate
        # Limitation: the most common insertion type may not be the best insertion to use
        pos_erate.fullIns = f_erate
        sub_ins, sub_align = pull_nlength_indels_info(err_info, 1, indel_alignment)
        if sub_ins:
            pos_erate.ins1 = ErrType(get_most_common_nt(sub_ins, pos_erate.ref, pos_erate.next), len(sub_ins)/depth, sub_align)
        sub_ins, sub_align = pull_nlength_indels_info(err_info, 2, indel_alignment)
        if sub_ins:
            pos_erate.ins2 = ErrType(get_most_common_nt(sub_ins, pos_erate.ref, pos_erate.next), len(sub_ins)/depth, sub_align)
    elif err_type.startswith("sub"):
        # substitution error rate
        pos_erate.fullSub = f_erate
        sub_nt = set(err_info)
        # pos_erate.sub = ErrType(sub_nt, err_info.count(get_most_common_nt(err_info))/depth, 0)
        pos_erate.sub = ErrType(sub_nt, f_erate, 0)

        # if error rates are for substitution adjustment
        if err_type == "sub_adjust":
            pos_erate.subAllErr = {"A": err_info.count("A")/depth,
                                   "C": err_info.count("C")/depth,
                                   "G": err_info.count("G")/depth,
                                   "T": err_info.count("T")/depth}
    else:
        print("WARNING: Could not determine error type for error rate storage.", file=sys.stderr)

    return

# assign erate information based on mpile information
def parse_pile(erate, mapped, pile_info, low=False, get_only=None):
    # initialize variable for bit recording of read positions    

    # keep track of insertions alignment
    ins_pos = 1

    for i,char in enumerate(mapped):
        # skip if needed
        if pile_info.skip > 0:
            pile_info.skip = pile_info.skip - 1
            continue
        elif char == "^":       # signify the mapped sequence begins at this position
            pile_info.skip = 1    # Following ASCII character signify mapping quality. Didn't use this information. 
            continue            # Limitation is not able to use read info in determining error rates
        elif char == "$":
            erate.record_read_end()
            continue

        indel_len = 0
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

            if args.messy: erate.coverage += 1 

            # check length of insertion
            indel_len, next_i = check_indel_len(mapped, (i+1))

            # store insertion information error rate
            pile_info.insStr.append("".join(mapped[next_i:(next_i+indel_len)]).upper())
            pile_info.skip = indel_len + (next_i - (i+1))
        elif char == "*":
            # store alignment information
            erate.add_read_alignment()
            ins_pos <<= 1

            # actual deletion of current position
            if args.messy: erate.coverage -= 1
        else:
            # store alignment information
            erate.add_read_alignment()
            ins_pos <<= 1            

            # all remaining possibilities should be only substitutions
            pile_info.subStr.append(char.upper())

    # pull error rates of each type of error
    if get_only is None:
        if pile_info.match:
            if len(pile_info.insStr) > 0:
                del pile_info.match[0: len(pile_info.insStr)]
            get_error_rate("total", erate, pile_info.match)

        if not low and pile_info.prevDel:
            get_error_rate("del", erate, pile_info.prevDel, pile_info.delReads)
        elif low and pile_info.subPrevDel:
            get_error_rate("del", erate, pile_info.subPrevDel, pile_info.subDelReads)
        else:
            pass

        if pile_info.insStr:
            get_error_rate("ins", erate, pile_info.insStr, ins_pos)      # insertions do not have quality

        if pile_info.subStr:
            get_error_rate("sub", erate, pile_info.subStr)
    else:
        if get_only == "del" and not low and pile_info.prevDel:
            get_error_rate("del", erate, pile_info.prevDel, pile_info.delReads)
        elif get_only == "del" and low and pile_info.subPrevDel:
            get_error_rate("del", erate, pile_info.subPrevDel, pile_info.subDelReads)
        elif get_only == "ins" and pile_info.insStr:
            get_error_rate("ins", erate, pile_info.insStr, ins_pos)      # insertions do not have quality
        elif get_only == "sub" and pile_info.subStr:
            get_error_rate("sub_adjust", erate, pile_info.subStr)
        else:
            return


# derivative of parse_pile but aimed at retriving next_del error only
def parse_del_only(pile_info, mapped, low=False):
    pile_info.skip = 0            # Counter

    for i,char in enumerate(mapped):
        # skip_counter if needed
        if pile_info.skip > 0:
            pile_info.skip -= 1
            continue
        elif char == "^":       # signify the mapped sequence begins at this position
            pile_info.skip = 1    # Following ASCII character signify mapping quality. Didn't use this information. 
            continue            # Limitation is not able to use read info in determining error rates
        elif char == "$":
            # last read ends so get rid of last bit in deletion alignment
            pile_info.remove_prev_del_read(low)
            continue
        
        indel_len = 0
        if char == "-":                                             # parse deletion
            pile_info.adjust_prev_read_as_del(low)                     # set the last read entered as a del
            indel_len, next_i = check_indel_len(mapped, (i+1))
            pile_info.append_prev_del_info("".join(mapped[next_i:(next_i+indel_len)]).upper(), low)
            pile_info.skip = indel_len + (next_i - (i+1))
        elif char == "+":                                           # parse insertion
            indel_len, next_i = check_indel_len(mapped, (i+1))
            pile_info.skip = indel_len + (next_i - (i+1))
        else:
            # store aligment
            pile_info.record_prev_del_read(low)
            continue

    return

# same as parse_del_only but returns only read alignment information
def parse_read_end_only(mapped):
    skip = 0            # counter
    read_align = 1

    for i,char in enumerate(mapped):    
        # skip_counter if needed
        if skip > 0:
            skip -= 1
            continue
        elif char == "^":       # signify the mapped sequence begins at this position
            skip = 1    # Following ASCII character signify mapping quality. Didn't use this information. 
            continue            # Limitation is not able to use read info in determining error rates
        elif char == "$":
            # last read ends so get rid of last bit in deletion alignment
            read_align ^= 1
            continue
        
        indel_len = 0
        if char == "-":                                             # parse deletion
            indel_len, next_i = check_indel_len(mapped, (i+1))
            skip = indel_len + (next_i - (i+1))
        elif char == "+":                                           # parse insertion
            indel_len, next_i = check_indel_len(mapped, (i+1))
            skip = indel_len + (next_i - (i+1))
        else:
            # store aligment
            read_align <<= 1
            continue

    return read_align

# go through mpileup line and calculate the error_rates
def parse_mpileup_line(pile_info, mpile, low):
    global _low_depth_count
    # pull relevant information from mpileup line 
    # set as variables for readability
    # Note: python is 0-based but samtools mpile is 1-based
    mpile_line = mpile[0].split()

    # for testing depth first
    depth = int(mpile_line[3])

    # initialize an Error object with contig, ref nucleotide and depth info
    pos_erate = Error(mpile_line[0], int(mpile_line[1]) - 1, mpile_line[2], depth)

    if depth < _depth_thres:
        # RECORD
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
        if depth < 1:
            pile_info.set_next_del_to_zero(low)
            return pos_erate
        else:
            pile_info.set_next_del_to_zero(low)
            parse_del_only(pile_info, mpile_line[4], low)
            return pos_erate
    else:
        pass

    # TODO: find better workaround. Ugly way of getting next nt. Temp store next nt for parsing. Will delete after parsing
    try:
        pos_erate.next = mpile[1].split()[2]
    except AttributeError:
        pos_erate.next = "-"

    # get error rates
    parse_pile(pos_erate, mpile_line[4], pile_info, low)

    # reset container for next mpileup line
    pile_info.reset_class()
    pile_info.set_next_del_to_zero(low)        

    # delete store next nt (unnecessary for future calculations)
    del pos_erate.next

    # get next deletion for next position
    parse_del_only(pile_info, mpile_line[4], low)

    return pos_erate

# populate the error rate list with both main and sub erate
def populate_erate_list(pile_info, erate_list, main, sub):
    if main.depth < 1 and sub.depth < 1:
        # both mpile are empty. Same for both
        main.lowErates = sub
        erate_list.append(main)
        if pile_info.lowCounter and pile_info.lowErates:
            # dump if on the right side of a low error region 
            pile_info.add_erate_to_err_rates(erate_list)
            pile_info.reset_low_erate()
    elif main.low and sub.depth > 0:
        # dump when entering low depth region
        if pile_info.lowErates: 
            pile_info.add_erate_to_err_rates(erate_list)
            pile_info.reset_low_erate()

            # Set counter to indicate that recently processed low depth region
            pile_info.set_low_counter()

        # low depth so sub erate becomes main erate
        sub.set_low_depth()
        sub.lowErate = main
        erate_list.append(sub)
    else:
        # add if right end
        if pile_info.lowCounter and len(pile_info.lowErates) > (args.passes * 10 + 1):
            pile_info.add_erate_to_err_rates(erate_list)
            pile_info.reset_low_erate()

        # save the low err info
        pile_info.save_low_erate(args.passes * 10 + 1, sub)
        erate_list.append(main)

# fill hp index and length information. 
# More function calls but wanted it to be separate from error rate parsing
def fill_hp_info(err_rates, last=False):
    index = len(err_rates) - 1
    ref_nt = err_rates[index].ref
    if index == 0:
        err_rates[0].hpStart = 0
        return
    else:
        if ref_nt == err_rates[index-1].ref:
            # is part of homopolymer
            err_rates[index].hpStart = err_rates[index-1].hpStart
        else:
            prev_hp_start = err_rates[index-1].hpStart
            err_rates[prev_hp_start].hpLen = index - prev_hp_start
            err_rates[index].hpStart = index
    
    if last:
        hp_length = 1
        while index != 0:
            if ref_nt != err_rates[index-hp_length].ref:
                start = err_rates[index].hpStart
                err_rates[start].hpLen = hp_length
                break
            else:
                hp_length += 1
    else:
        return
    
    return

# return a list of lists with error rate objects for each position of each contig
def get_mpile_error_rates(main_mpile, sub_mpile=None):
    # error rate storage variable
    err_rates = []

    # object to hold mpileup line info
    container = MpileInfo()

    # iterate over mpiles to get error rates
    # populate the error rate list with Error object with error rates
    if sub_mpile != "":
        # parse at same time
        with open(main_mpile) as main:
            with open(sub_mpile) as sub:
                # get first lines
                main_to_parse, sub_to_parse = next(main), next(sub)
                for mline, sline in zip(main, sub):
                    main_erate = parse_mpileup_line(container, [main_to_parse, mline], False)
                    sub_erate = parse_mpileup_line(container, [sub_to_parse, sline], True)
                    populate_erate_list(container, err_rates, main_erate, sub_erate)
                    fill_hp_info(err_rates)
                    main_to_parse, sub_to_parse = mline, sline

                # parse the last line
                main_erate = parse_mpileup_line(container, [main_to_parse, None], False)
                sub_erate = parse_mpileup_line(container, [sub_to_parse, None], True)
                populate_erate_list(container, err_rates, main_erate, sub_erate)
                fill_hp_info(err_rates, True)
    else:
        with open(main_mpile) as main:
            main_to_parse = next(main)
            for mline in main:
                main_erate = parse_mpileup_line(container, [main_to_parse, mline], False)
                err_rates.append(main_erate)
                fill_hp_info(err_rates)
                main_to_parse = mline

            # parse last line
            main_erate = parse_mpileup_line(container, [main_to_parse, None], False) 
            err_rates.append(main_erate)
            fill_hp_info(err_rates, True)
    
    # finished getting all info
    return err_rates

# retrieve the proper erate to parse. Mostly important for hybrid correction
def retrieve_pos_erate(err_rates, pos, adjust, low):
    if args.low and not adjust and err_rates[pos].low != low:
        try:
            return err_rates[pos].lowErate
        except AttributeError:
            return err_rates[pos]
    else:
        return err_rates[pos]

# adjust read alignment binary depending on which read ends
def adjust_read_alignment(read_align, read_ends):
    b_read_ends = bin(read_ends)
    next_align = list(bin(read_align))
    if len(b_read_ends) != len(next_align): 
        print("BUG: Mismatch in alignment and depth in adjust_read_alignment", file=sys.stderr)
        
    for bit in range(3, len(bin(read_ends))):
        if b_read_ends[bit] == "1":
            # read at this bit ends. Remove from alignment
            next_align[bit] = ""
        else:
            pass
    return int("".join(next_align), 2)

# get read error rates for a particular position with unique position 
# ugly way of doing uniq read alignment only. Not efficient.
def get_uniq_err_rates(pos_erate, err_type, err_len, alignment, uniq_pos_only=False):
    try:
        if err_type == "del":
            if err_len == 1:
                if uniq_pos_only: return pos_erate.del1.reads_w_uniq_err(alignment)
                return pos_erate.del1.uniq_err_erate(alignment), pos_erate.del1.subset_out_uniq_err(alignment)
            elif err_len == 2:
                if uniq_pos_only: return pos_erate.del2.reads_w_uniq_err(alignment)
                return pos_erate.del2.uniq_err_erate(alignment), pos_erate.del2.subset_out_uniq_err(alignment)
            else:
                if uniq_pos_only: return 0
                return pos_erate.fullDel, alignment
        elif err_type == "ins":
            if err_len == 1:
                if uniq_pos_only: return pos_erate.ins1.reads_w_uniq_err(alignment)
                return pos_erate.ins1.uniq_err_erate(alignment), pos_erate.ins1.subset_out_uniq_err(alignment)
            elif err_len == 2:
                if uniq_pos_only: return pos_erate.ins2.reads_w_uniq_err(alignment)
                return pos_erate.ins2.uniq_err_erate(alignment), pos_erate.ins2.subset_out_uniq_err(alignment)
            else:
                if uniq_pos_only: return 0
                return pos_erate.fullIns, alignment
        else:
            print("BUG: Could not determine error type in get_uniq_err_rates. DEBUG.", file=sys.stderr)
    except AttributeError:
        if uniq_pos_only:
            return 0
        else:
            return 0, alignment

# reconstruction of bit alignment using previously ended reads
# faster than bit conversion if less than 1M conversions
def reconstruct_prev_alignment(read_align, read_ends):
    b_read_ends = bin(read_ends)
    next_align = list(bin(read_align))
    align_index = 3
    while align_index < len(b_read_ends):
        if b_read_ends[align_index] == "1":
            # if previous read ended at index, insert a 0 in that position
            next_align[align_index:align_index] = ["0"]
        align_index += 1

    return int("".join(next_align), 2)

# inefficient but these cases should be rare enough that it doesn't matter
def join_single_errs(err_rates, pos_list, low):
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
                    pos1_read_end = retrieve_pos_erate(err_rates, pos1, False, low).readEnd
                    if (pos1_read_end ^ (1 << (pos1_read_end.bit_length() - 1))) > 0:
                        to_compare_v1 = adjust_read_alignment(to_compare_v1, pos1_read_end)

                    for step in range(1, (pos2 - pos1 + 1)):
                        # new variable are for reability
                        # adjust to_compare as we step forward from pos 1
                        pos1_read_end = retrieve_pos_erate(err_rates, pos1+step, False, low).readEnd

                        # pos 1: adjust for new reads
                        bit_diff = pos1_read_end.bit_length() - to_compare_v1.bit_length()
                        if bit_diff > 0: 
                            # only add 0 bits since new reads should not be considered
                            to_compare_v1 <<= bit_diff
                        elif bit_diff < 0:
                            print("BUG: Alignment bit difference in function join_single_errs", file=sys.stderr)
                        else:
                            pass

                        # pos 1: adjust for reads ending
                        # also skip last end correction
                        if (pos1 + step < pos2) and (pos1_read_end ^ (1 << (pos1_read_end.bit_length() - 1))) > 0:
                            to_compare_v1 = adjust_read_alignment(to_compare_v1, pos1_read_end)

                        # adjust to compare as we step backwards from pos 2
                        pos2_read_end = retrieve_pos_erate(err_rates, pos2-step, False, low).readEnd

                        # Pos 2: adjust by read ends in previous position first
                        if (pos2_read_end ^ (1 << (pos2_read_end.bit_length() - 1))) > 0:
                            to_compare_v2 = reconstruct_prev_alignment(to_compare_v2, pos2_read_end)
                        
                        # Pos 2: remove any new reads that were incorporated this round
                        bit_diff = to_compare_v2.bit_length() - pos2_read_end.bit_length()
                        if bit_diff > 0:
                            to_compare_v2 >>= bit_diff
                        elif bit_diff < 0:
                            print("BUG: Alignment bit difference in function join_single_errs", file=sys.stderr)
                        else:
                            pass

                    # make sure that corrected alignment are the same length
                    if (pos_list[first][1]).bit_length() != to_compare_v2.bit_length(): 
                        print("BUG: Alignment bit difference in forward alignment in function join_single_errs", file=sys.stderr)
                    if (pos_list[second][1]).bit_length() != to_compare_v1.bit_length(): 
                        print("BUG: Alignment bit difference in reverse alignment in function join_single_errs", file=sys.stderr)

                    uniq1 = pos_list[first][1] & (to_compare_v2 ^ (1 << (to_compare_v2.bit_length() - 1)))
                    uniq2 = pos_list[second][1] & (to_compare_v1 ^ (1 << (to_compare_v1.bit_length() - 1)))

                    if bin(uniq1).count("1") != bin(uniq2).count("1"): 
                        print("BUG: Mismatch between forward and reverse bit alignments in function join_single_errs", file=sys.stderr)
                    if uniq1 == 0: 
                        # nothing between these two position
                        continue

                    # print("Pos 1: {}, Pos 2: {}, b1: {}, b2: {}". format(pos1, pos2, bin(uniq1), bin(uniq2)))
                    pos_list.append([[pos1, pos2], [uniq1, uniq2]])

                    # correct pos 1 and pos 2 for future correction
                    pos_list[first][1] ^= uniq1
                    pos_list[second][1] ^= uniq2

                    # record to return 
                    pos2_erate = retrieve_pos_erate(err_rates, pos2, False, low)
                    pos_candidate += bin(uniq2).count("1")/pos2_erate.depth
                    uniq_errs ^= uniq2

    # remove anything that doesn't have uniq
    for pos in range(len(pos_list)-1, -1, -1):
        if type(pos_list[pos][0]).__name__ == "list":
            continue
        else:
            align = (1 << (pos_list[pos][1].bit_length() - 1))
            if (pos_list[pos][1] ^ align) == 0:
                del pos_list[pos]

    return pos_candidate, uniq_errs

# add single length errors to error rate
def add_single_length_errors(err_rates, current_pos, current_align, pos_list, low):
    single_err = 0
    for single in pos_list:
        single_pos = single[0]
        to_compare = single[1]
        
        if type(single_pos).__name__ == "list":
            continue
        else:
            # adjust the comparison alignment to current alignment length
            for step in range(single_pos, current_pos+1):
                # new variable are for reability
                # adjust to_compare as we step forward from pos 1
                pos1_read_end = retrieve_pos_erate(err_rates, step, False, low).readEnd

                # pos 1: adjust for new reads
                bit_diff = pos1_read_end.bit_length() - to_compare.bit_length()
                if bit_diff > 0: 
                    # only add 0 bits since new reads should not be considered
                    to_compare <<= bit_diff
                elif bit_diff < 0:
                    print("BUG: Alignment bit difference in function add_single_length_errors", file=sys.stderr)
                else:
                    pass

                # pos 1: adjust for reads ending
                if (step < current_pos) and (pos1_read_end ^ (1 << (pos1_read_end.bit_length() - 1))) > 0:
                    to_compare = adjust_read_alignment(to_compare, pos1_read_end)

            # compare adjusted
            uniq_err = current_align & to_compare
            if uniq_err != 0: uniq_err ^= (1 << (uniq_err.bit_length() - 1))

            # compare and subset out unique single errors
            current_align ^= uniq_err

            # add error rate
            single_pos_erate = retrieve_pos_erate(err_rates, single_pos, False, low)
            single_err += bin(uniq_err).count("1")/err_rates[single_pos].depth

    return single_err, current_align

# In this case, error length represent the error information to adjust
def last_round_error_adjustments(err_rates, err_type, err_len, pos, alignment, low):
    pos_erate = retrieve_pos_erate(err_rates, pos, False, low)
    try:
        if err_len == 1:
            if err_type == "del":
                if ((pos_erate.del1.reads ^ alignment) ^ (1 << (pos_erate.del1.reads.bit_length() - 1))) == 0:
                    del pos_erate.del1
                else:
                    pos_erate.del1.erate -= bin(alignment).count("1")/pos_erate.depth
                    pos_erate.del1.reads ^= alignment
            elif err_type == "ins":
                if ((pos_erate.ins1.reads ^ alignment) ^ (1 << (pos_erate.ins1.reads.bit_length() - 1))) == 0:
                    del pos_erate.ins1
                else:
                    pos_erate.ins1.erate -= bin(alignment).count("1")/pos_erate.depth
                    pos_erate.ins1.reads ^= alignment
        elif err_len == 2:
            if err_type == "del":
                old_reads = pos_erate.del2.reads
                if ((pos_erate.del2.reads ^ alignment) ^ (1 << (pos_erate.del2.reads.bit_length() - 1))) == 0:
                    del pos_erate.del2
                else:
                    pos_erate.del2.erate -= bin(alignment).count("1")/pos_erate.depth
                    pos_erate.del2.reads ^= alignment

                # set aside another err type that can be used in calculation but not in error diff_err selection
                try:
                    if ((pos_erate.del3.reads ^ alignment) ^ (1 << (pos_erate.del3.reads.bit_length() - 1))) == 0:
                        del pos_erate.del3
                    else:
                        pos_erate.del3.erate -= bin(alignment).count("1")/pos_erate.depth
                        pos_erate.del3.reads ^= alignment
                except AttributeError:
                    pos_erate.del3 = ErrType("N", 
                                             bin(alignment).count("1")/pos_erate.depth,
                                             (1 << (old_reads.bit_length() - 1)) | alignment)
            elif err_type == "ins":
                old_reads = pos_erate.ins2.reads
                if ((pos_erate.ins2.reads ^ alignment) ^ (1 << (pos_erate.ins2.reads.bit_length() - 1))) == 0:
                    del pos_erate.ins2
                else:
                    pos_erate.ins2.erate -= bin(alignment).count("1")/pos_erate.depth
                    pos_erate.ins2.reads ^= alignment 
                try:
                    if ((pos_erate.ins3.reads ^ alignment) ^ (1 << (pos_erate.ins3.reads.bit_length() - 1))) == 0:
                        del pos_erate.ins3
                    else:
                        pos_erate.ins3.erate -= bin(alignment).count("1")/pos_erate.depth
                        pos_erate.ins3.reads ^= alignment 
                except AttributeError:
                    pos_erate.ins3 = ErrType("N", 
                                                          bin(alignment).count("1")/pos_erate.depth,
                                                          (1 << (old_reads.bit_length() - 1)) | alignment)
            else:
                print("BUG: Issue with determining error type in last_round_error_adjustments for position {}".format(pos))
        else:
            print("BUG: Issue with determining error length in last_round_error_adjustments for position {}".format(pos))
    except AttributeError:
        print("BUG: No indel attribute found in position. May have been prematurely removed", file=sys.stderr)
    return

def get_singlar_indel_erate_in_range(pos_erate, err_type, err_len, pos, align, last_round, pos_info_storage):
    uniq_align = get_uniq_err_rates(pos_erate, err_type, err_len, align, uniq_pos_only=True)
    pos_candidate, read_align = get_uniq_err_rates(pos_erate, err_type, err_len, align)
    if pos_candidate > 0:
        # store orig positions too
        pos_info_storage.append([pos, uniq_align])

    if last_round:
        if err_len == 1:
            if err_type == "del":
                try:
                    uniq_align_2 = pos_erate.del3.reads_w_uniq_err(read_align)
                    pos_candidate_2 = pos_erate.del3.uniq_err_erate(read_align)
                    read_align = pos_erate.del3.subset_out_uniq_err(read_align)        
                except AttributeError:
                    pos_candidate_2 = 0; uniq_align_2 = 0
            elif err_type == "ins":
                try:
                    uniq_align_2 = pos_erate.ins3.reads_w_uniq_err(read_align)
                    pos_candidate_2 = pos_erate.ins3.uniq_err_erate(read_align)
                    read_align = pos_erate.ins3.subset_out_uniq_err(read_align)
                except AttributeError:
                    pos_candidate_2 = 0; uniq_align_2 = 0
            if pos_candidate_2 > 0:
                pos_info_storage.append([pos * -1, uniq_align_2])
                uniq_align_2 ^= (1 << (uniq_align_2.bit_length()-1))
                pos_candidate += pos_candidate_2
        else:
            uniq_align_2 = 0
    else:
        uniq_align_2 = 0

    return pos_candidate, read_align, (uniq_align ^ uniq_align_2)

# find max indel error rates for range
# last_found option should only be used for finding errors
# having hp = True and find = 0 indicates that it is an adjustment step
def get_best_indel_erate_in_range(err_rates, err_type, err_len, start, end, hp=False, last_round=0, uniq_erates=None, find=0, adjust=False):
    if start < 0: start = 0
    if end > len(err_rates): end = len(err_rates)

    # to keep track of top_non_hp score. 0: position, 1: score
    if adjust: top_non_hp = []

    # Keeps track of position with different errors
    orig_pos_nlen = []
    if last_round: pos_diff_nlen = []

    # keep track of total error rate
    err_candidate = 0
    
    # start off with the original candidate first
    # also determine first read alignment
    # read alignment for determining error rate of unique errors only
    if hp and err_type == "ins" and (start - 1) >= 0:
        orig_start = start - 2
    else:
        orig_start = start - 1

    if find > 0:
        orig_pos = find
        find_low = err_rates[orig_pos].low
        orig_align = (1 << err_rates[orig_pos].readEnd.bit_length()) - 1
        pos_candidate, orig_align, uniq_align = get_singlar_indel_erate_in_range(err_rates[orig_pos], err_type, err_len, orig_pos, orig_align, last_round, orig_pos_nlen)
        err_candidate += pos_candidate

        for reconstruct_pos in range(orig_pos-1, orig_start, -1):
            reconstruct_erate = retrieve_pos_erate(err_rates, reconstruct_pos, adjust, find_low)
            reconstruct_read_end = reconstruct_erate.readEnd

            if (reconstruct_read_end ^ (1 << (reconstruct_read_end.bit_length() - 1))) > 0:
                orig_align = reconstruct_prev_alignment(orig_align, reconstruct_read_end)
            
            # remove any new reads that were incorporated this round
            bit_diff = orig_align.bit_length() - reconstruct_read_end.bit_length()
            if bit_diff > 0:
                orig_align >>= bit_diff
            elif bit_diff < 0:
                print("BUG: Alignment bit difference in function get_best_indel_erate_in_range", file=sys.stderr)
            else:
                pass

        read_align = orig_align     # two copies, one to follow the original alignment (no read addition)
    else:
        find_low = False
        read_align = (1 << err_rates[orig_start + 1].readEnd.bit_length()) - 1

    # get extra error info from previous position only if insertion from previous position is current hp
    if hp and err_type == "ins":
        if (start - 1) >= 0:
            current_pos_erate = retrieve_pos_erate(err_rates, start - 1, adjust, find_low)

            if err_len == 1:
                try:
                    nt = current_pos_erate.ins1.nt
                    current = (nt == err_rates[start].ref)
                except AttributeError:
                    current = False
            else:
                try:
                    nt = current_pos_erate.ins2.nt
                    current = (nt[-1] == err_rates[start].ref) and (nt[0] != current_pos_erate.ref)          # only end needs to match since it will be nt preceding current hp
                except AttributeError:
                    current = False
            
            # parse
            if current:
                # adjust read alignment based on reads that end this position
                pos_candidate, read_align, uniq_align = get_singlar_indel_erate_in_range(current_pos_erate, err_type, err_len, start-1, read_align, 0, orig_pos_nlen)
                err_candidate += pos_candidate
            else:
                pass

            # adjust read alignment based on reads that end this position
            # only adjust if there is a read ending at this position
            if (current_pos_erate.readEnd ^ (1 << current_pos_erate.readEnd.bit_length() - 1)) > 0:
                if find > 0: 
                    adjusted_read_end = current_pos_erate.readEnd >> (current_pos_erate.readEnd.bit_length() - orig_align.bit_length())
                    orig_align = adjust_read_alignment(orig_align, adjusted_read_end) 
                read_align = adjust_read_alignment(read_align, current_pos_erate.readEnd)
            else:
                pass
        else:
            pass
    else:
        pass

    for pos in range(start, end):
        # get proper mpile to assess if using multiple
        current_pos_erate = retrieve_pos_erate(err_rates, pos, adjust, find_low)

        # make sure correct bit size
        bit_diff = current_pos_erate.readEnd.bit_length() - read_align.bit_length()
        if bit_diff > 0:
            # adjust for new read additions. 
            read_align <<= bit_diff
            read_align ^= ((1 << bit_diff) - 1)

            # If for finding errors, then add new bits on found errors up to original position
            if find > 0 and pos <= orig_pos:
                orig_uniq_align = orig_pos_nlen[0][1]
                for reconstruct_pos in range(orig_pos, pos - 1, -1):
                    reconstruct_erate = retrieve_pos_erate(err_rates, reconstruct_pos, adjust, find_low)
                    reconstruct_read_end = reconstruct_erate.readEnd
                    orig_uniq_align = reconstruct_prev_alignment(orig_uniq_align, reconstruct_read_end)
                    
                    # remove any new reads that were incorporated this round
                    bit_diff_2 = orig_uniq_align.bit_length() - reconstruct_read_end.bit_length()
                    if bit_diff_2 > 0:
                        orig_uniq_align >>= bit_diff_2
                    elif bit_diff_2 < 0:
                        print("BUG: Alignment bit difference in correcting alignment during finding errors", file=sys.stderr)
                    else:
                        pass

                # mask out only new bits added this round
                if read_align.bit_length() != orig_uniq_align.bit_length(): 
                    print("BUG: Alignment bit difference when determing remaining intact genomes", file=sys.stderr)
                orig_uniq_align &= ((1 << bit_diff) - 1)
                read_align ^= orig_uniq_align
            else:
                pass
        elif bit_diff < 0:
            if current_pos_erate.depth < _depth_thres:
                pass
            else:
                print("BUG: Unique errors alignment length greater than read depth on contig {} at position {}".format(err_rates[pos].contig, pos), file=sys.stderr)
        else:
            pass

        if current_pos_erate.depth < _depth_thres:
            # just preserve read alignment if it continues beyond low depth region
            if (current_pos_erate.readEnd ^ (1 << current_pos_erate.readEnd.bit_length() - 1)) > 0:
                if find > 0:
                    adjusted_read_end = current_pos_erate.readEnd >> (current_pos_erate.readEnd.bit_length() - orig_align.bit_length())
                    orig_align = adjust_read_alignment(orig_align, adjusted_read_end)            
                read_align = adjust_read_alignment(read_align, current_pos_erate.readEnd)

            # don't parse error rate due to low depth
            # adjust read alignment based on reads that end this position
            # only adjust if there is a read ending at this position
            else:
                pass
            continue
        else:
            pass    # continue as normal

        if err_type == "del":
            if abs(find) == pos:
                uniq_align = 0
            elif not adjust:
                pos_candidate, read_align, uniq_align = get_singlar_indel_erate_in_range(current_pos_erate, err_type, err_len, pos, read_align, last_round, orig_pos_nlen)
                err_candidate += pos_candidate
            else:
            # only consider adjustment error rate if at least one position contains the correct error length 
                if err_len == 1:
                    pos_candidate = 0
                else:
                    pos_candidate, placeholder, uniq_align = get_singlar_indel_erate_in_range(current_pos_erate, err_type, 1, pos, read_align, last_round, orig_pos_nlen)
                err_candidate += current_pos_erate.fullDel - pos_candidate
        elif hp and err_type == "ins":
            # should not enter this area if last round
            # TODO: redo this area since it is a mess
            # don't want to add if the err rate is part of next hp
            uniq_align = get_uniq_err_rates(current_pos_erate, err_type, err_len, read_align, uniq_pos_only=True)
            hp_candidate, sub_align = get_uniq_err_rates(current_pos_erate, err_type, err_len, read_align)
            if abs(find) == pos:
                hp_candidate = 0
            if hp_candidate > 0:
                if err_len == 1:
                    try:
                        err_nt = current_pos_erate.ins1.nt
                    except AttributeError:
                        err_nt = ""
                    if adjust: hp_candidate = current_pos_erate.fullIns
                elif err_len == 2:
                    try:
                        err_nt = current_pos_erate.ins2.nt
                    except AttributeError:
                        err_nt = ""
                    # for adjustments only
                    if adjust: 
                        hp_candidate, placeholder = get_uniq_err_rates(current_pos_erate, err_type, 1, read_align)
                        hp_candidate = (current_pos_erate.fullIns - hp_candidate)
                else:
                    if (err_len - 2) == 1:
                        try:
                            err_nt = current_pos_erate.ins1.nt
                            sub_align = current_pos_erate.ins1.subset_out_uniq_err(sub_align)
                        except AttributeError:
                            err_nt = ""; sub_align = read_align
                    else:
                        try:
                            err_nt = current_pos_erate.ins2.nt
                            sub_align = current_pos_erate.ins2.subset_out_uniq_err(sub_align)
                        except AttributeError:
                            err_nt = ""; sub_align = read_align
                # make sure not next hp for edges
                if (pos + 1) == end:
                    if err_nt == "":
                        err_candidate += 0 
                    elif err_nt[0] == current_pos_erate.ref:
                        err_candidate += hp_candidate
                        read_align = sub_align
                        orig_pos_nlen.append([pos, uniq_align])
                    else:                  
                        next_hp = start + err_rates[start].hpLen
                        try:
                            next_hp_nt = err_rates[next_hp].ref
                        except IndexError:
                            next_hp_nt = ""

                        if err_nt[-1] == next_hp_nt:
                            err_candidate += 0
                        else:
                            if adjust:
                                # matches nothing so treat it as unique
                                if len(top_non_hp) > 0 and top_non_hp[1] >= hp_candidate:
                                    err_candidate += 0
                                elif len(top_non_hp) > 0 and top_non_hp[1] < hp_candidate:
                                    top_non_hp[0] = pos; top_non_hp[1] = hp_candidate
                                else:
                                    top_non_hp.extend([pos, hp_candidate])
                            else:
                                err_candidate += hp_candidate               # matches nothing
                                read_align = sub_align
                                orig_pos_nlen.append([pos, uniq_align])
                else:
                    if err_nt != "":
                        if adjust and err_nt[0] == current_pos_erate.ref and err_nt[0] == err_nt[-1]:
                            err_candidate += hp_candidate
                        elif adjust and (err_nt[0] != current_pos_erate.ref or err_nt[-1] != current_pos_erate.ref):
                            # at least one does not match indicating that an hp is broken and should be treated as separate
                            if len(top_non_hp) > 0 and top_non_hp[1] >= hp_candidate:
                                err_candidate += 0
                            elif len(top_non_hp) > 0 and top_non_hp[1] < hp_candidate:
                                top_non_hp[0] = pos; top_non_hp[1] = hp_candidate
                            else:
                                top_non_hp.extend([pos, hp_candidate])
                        else:
                            err_candidate += hp_candidate                       # matches nothing
                            read_align = sub_align
                            orig_pos_nlen.append([pos, uniq_align])
                    else:
                        err_candidate += 0
            else:
                err_candidate += 0
        elif err_type == "ins":
            if abs(find) == pos:
                uniq_align = 0
            else:
                pos_candidate, read_align, uniq_align = get_singlar_indel_erate_in_range(current_pos_erate, err_type, err_len, pos, read_align, last_round, orig_pos_nlen)
                err_candidate += pos_candidate
        else:
            print("BUG: could not parse error type. Debug", file=sys.stderr)

        # last round calculations
        if last_round: 
            err_len_2 = err_len - (-1)**err_len
            pos_candidate, sub_align = get_uniq_err_rates(current_pos_erate, err_type, err_len_2, read_align)
            
            # keep track of it since something
            if pos_candidate > 0:
                uniq_align = get_uniq_err_rates(current_pos_erate, err_type, err_len_2, read_align, uniq_pos_only=True)
                pos_diff_nlen.append([pos, uniq_align])

                if err_len == 1:
                    err_candidate += pos_candidate
                    read_align = sub_align
                    # to match up with uniq_alignment of err_length 2. Consistency
                    pos_diff_nlen[-1][1] ^= (1 << (uniq_align.bit_length()-1))
                elif err_len == 2:
                    if len(pos_diff_nlen) < 2:
                        pass
                    else:
                        # add potential candidates together
                        pos_candidate, sub_align = join_single_errs(err_rates, pos_diff_nlen, low=find_low)
                        err_candidate += pos_candidate
                        read_align ^= sub_align
                else:
                    print("BUG: cannot parse error length. DEBUG", file=sys.stderr)
            else:
                pass
        else:
            pass

        # adjust read alignment based on reads that end this position
        # only adjust if there is a read ending at this position
        if (current_pos_erate.readEnd ^ (1 << current_pos_erate.readEnd.bit_length() - 1)) > 0:
            if find > 0:
                adjusted_read_end = current_pos_erate.readEnd >> (current_pos_erate.readEnd.bit_length() - orig_align.bit_length())
                orig_align = adjust_read_alignment(orig_align, adjusted_read_end)            
            read_align = adjust_read_alignment(read_align, current_pos_erate.readEnd)
        else:
            pass

    # adjust err_candidate if there any genomes without an indel
    if find > 0 and err_candidate > _indel_thres:
        orig_align = read_align >> (read_align.bit_length() - orig_align.bit_length())
        intact_genomes = (bin(orig_align).count("1") - 1)
        if intact_genomes > 0:
            # print("Position range {} to {} still has {} non indel genome(s)".format(start, end, intact_genomes))
            alt_err_candidate = 1 - (intact_genomes/err_rates[orig_pos].depth)

            # compare the two err_candidates. Take the lowest one
            err_candidate = min(err_candidate, alt_err_candidate)

    # set erate to zero prematurely so they won't be considered when determining error position
    # only if it passes the threshold
    # TODO: clean up. This is a mess
    if last_round and pos_diff_nlen:
        # only re-adjust if pass threshold
        # can't directly modify since slots
        if err_candidate > last_round:
            if err_len == 1:
                for diff_err in pos_diff_nlen:
                    last_round_error_adjustments(err_rates, err_type, 2, diff_err[0], diff_err[1], low=find_low)
                else:
                    pass
            elif err_len == 2:
                for two_pos in pos_diff_nlen:
                    if type(two_pos[0]).__name__ != "list": continue
                    for n in range(2):
                        diff_err = [two_pos[0][n], two_pos[1][n]]
                        last_round_error_adjustments(err_rates, err_type, 1, diff_err[0], diff_err[1], low=find_low)
                    else:
                        pass
                pass
            else:
                print("BUG: cannot parse error length. DEBUG", file=sys.stderr)
        else:
            # couldn't find compensation. Check shorter error length. Only for error length 2
            if err_len == 2:
                # check if the single nt errors can bump error rate above threshold.
                # if so, there is enough evidence for singe nt error
                # too many for loops
                pos_candidate, sub_align = add_single_length_errors(err_rates, pos, read_align, pos_diff_nlen, low=find_low)
                err_candidate += pos_candidate

                # adjust err_candidate if there any genomes without an indel
                if find > 0 and err_candidate > _indel_thres:
                    orig_align = sub_align >> (sub_align.bit_length() - orig_align.bit_length())
                    intact_genomes = (bin(orig_align).count("1") - 1)
                    if intact_genomes > 0:
                        # print("Position range {} to {} still has {} non indel genome(s)".format(start, end, intact_genomes))
                        alt_err_candidate = 1 - (intact_genomes/err_rates[orig_pos].depth)

                        # compare the two err_candidates. Take the lowest one
                        err_candidate = min(err_candidate, alt_err_candidate)

                # recheck if new single nt error is greater than threshold
                if err_candidate > last_round:
                    # go through joined single first
                    for index in range(len(pos_diff_nlen)-1, -1, -1):
                        if type(pos_diff_nlen[index][0]).__name__ != "list": continue
                        for n in range(2):
                            diff_err = [pos_diff_nlen[index][0][n], pos_diff_nlen[index][1][n]]
                            last_round_error_adjustments(err_rates, err_type, 1, diff_err[0], diff_err[1], low=find_low)
                        # delete when done adjusting
                        del pos_diff_nlen[index] 

                    # only readjust error info of length 2 since regular code should adjust errors of length 1
                    for diff_err in orig_pos_nlen:
                        # for consistency with output of join_single_err
                        diff_err[1] ^= (1 << (diff_err[1].bit_length()-1))
                        last_round_error_adjustments(err_rates, err_type, 2, diff_err[0], diff_err[1], low=find_low)

                    # change indel error to negative to signal to upper level function that error 
                    # length was changed from 2 to 1.
                    # keep read_align un changed so we can re adjust erate used when needed
                    err_candidate *= -1

                    # delete orig_pos_nlen since we already  changed error rate
                    del orig_pos_nlen[:]
                else:
                    pass
            else:
                pass
    else:
        pass

    if not adjust and uniq_erates is not None: 
        # can't pass by value and did not want to import deepcopy
        if err_candidate < 0 and err_len == 2:
            # only done when switching for err_len of 2 to 1
            # need to readjust the read_align to get true uniq_err
            for e in pos_diff_nlen:
                if type(e[0]).__name__ == "list": continue
                adjusted_read_align = read_align

                for reconstruct_pos in range(end - 1, e[0] - 1, -1):
                    reconstruct_erate = retrieve_pos_erate(err_rates, reconstruct_pos, adjust, find_low)
                    reconstruct_read_end = reconstruct_erate.readEnd
                    if (reconstruct_read_end ^ (1 << (reconstruct_read_end.bit_length() - 1))) > 0:
                        adjusted_read_align = reconstruct_prev_alignment(adjusted_read_align, reconstruct_read_end)
                    # remove any new reads that were incorporated this round
                    bit_diff = adjusted_read_align.bit_length() - reconstruct_read_end.bit_length()
                    if bit_diff > 0:
                        adjusted_read_align >>= bit_diff
                    elif bit_diff < 0:
                        print("BUG: Alignment bit difference during unique error determination of positions used in error finding", file=sys.stderr)
                    else:
                        pass   
                # get the actual uniq pos 
                current_pos_erate = retrieve_pos_erate(err_rates, e[0], adjust, find_low)
                e[1] = get_uniq_err_rates(current_pos_erate, err_type, err_len_2, adjusted_read_align, uniq_pos_only=True)
                uniq_erates.append(e)
        else:
            for e in orig_pos_nlen:
                uniq_erates.append(e)
    elif adjust and uniq_erates is not None and top_non_hp:
        # save the top non hp to uniq_erates
        # regular list by value copy does not work
        uniq_erates.extend(top_non_hp)
    else:
        pass

    return err_candidate

# determine error position, depth, reference nt and replacement nt
def determine_range_error_info(err_rates, err_type, err_len, start, end, find=0, hp=False, adjust=False, existing_errs=[], check_diff_nlen=False):
    # range adjustment
    if start < 0: start = 0
    if end > len(err_rates): end = len(err_rates)

    if hp and err_type == "ins":
        if (start - 1) >= 0: start -= 1

    if args.low:
        find_low = err_rates[find].low
    else:
        find_low = False

    score_info = {}
    hp_err_pos = {}
    for pos in range(start, end):
        curr_hp = err_rates[pos].hpStart
        curr_hp_nt = err_rates[pos].ref
        curr_hp_len = err_rates[curr_hp].hpLen

        if  (pos + 1) < (curr_hp + curr_hp_len):
            next_hp = curr_hp
            next_hp_nt = curr_hp_nt
            next_hp_len = curr_hp_len
        else:
            # only get next hp information if position is at the junction between two homopolymer sequences
            try:
                next_hp = curr_hp + curr_hp_len
                next_hp_nt = err_rates[next_hp].ref
                next_hp_len = err_rates[next_hp].hpLen
            except IndexError:
                next_hp = 0
                next_hp_nt = ""
                next_hp_len = 0

        hp_length = 0
        hp_match = "none"
        
        # actual position information
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
                        if check_diff_nlen: erate += current_pos_erate.del3.erate
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
                        erate = current_pos_erate.fullDel - current_pos_erate.del1.erate
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
                        if check_diff_nlen: erate += current_pos_erate.ins3.erate
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
                        erate = current_pos_erate.fullIns - current_pos_erate.ins1.erate
                    except AttributeError:
                        erate = current_pos_erate.fullIns
        else:
            print("BUG: Could not access error nucleotide at position {}".format(pos), file=sys.stderr)
        
        # determine how to use erate
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
                # so matching must be the end of the insertion length. Ex. if "A" match "A"
                # if "AC" match "C" since likely next hp can be CCC -> NACCCC
                # want to weigh down cases where +2AC and next hp is AAA -> NACAAA

                if curr_hp_nt == err_nt[0] and next_hp_nt == err_nt[-1]:
                    if curr_hp_len >= next_hp_len:
                        hp_length = curr_hp_len
                        hp_match = "current"
                    else:
                        if adjust: 
                            # need special case for junction insertions when adjusting. Both hp are valid
                            # score would be added to both candidates
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
                    hp_length = 0.5                                                      # set the penalizing score if didn't find anything
            else:                                                                        # either at beginning or end of current hp
                if err_nt[0] == err_nt[-1]:                                              # check if a homopolymer ins/deletion
                    if curr_hp_nt == err_nt[0]:
                        hp_length = curr_hp_len
                        hp_match = "current"
                    else:
                        hp_length = 0.5
                else:                                                                    # Non-hp error implied to break a hp in the middle if insertion. Not possible for deletion
                    if curr_hp_nt == err_nt[0] or curr_hp_nt == err_nt[-1]:              # Penalize score for breaking hp. 
                        hp_length = curr_hp_len/2 
                        hp_match = "current"
                    else:
                        hp_length = 0.25 
        else:
                print("BUG: Can't tell error length to determine error candidate score with length {} at position {}.".format(err_len, pos), file=sys.stderr)

        # calculate score
        if adjust:
            # if getting info for adjustment, you only want info for the current hp
            # next hp will be parsed separately
            # TODO: may change
            if err_type == "ins" and (pos + 1) == end and hp_match == "next":
                # only insertion as deletion always considers current position as well
                continue                                # skip if edge of hp and matches next hp
            else:
                if err_len > 1 and err_nt[0] != err_nt[-1]: 
                    hp_match = "junction"                   # for adjustments, hp matching junctions would be considered as separate
            
            score = erate
        else:
            # default
            score = hp_length * erate                   # proportion of total deletions are in this position
                                                        # the longer the hp at the del position, the higher we want to score this position. 
            if pos in existing_errs: score = 0; #print("Duplicate position found at {}".format(pos))

        # remove position as possibility when it's already being considered as a position of error
        # set the score as 0
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
        else:                                       # unique score
            key = 0 - pos                           # use negative to denote uniqueness in case clash with hp positions
            score_info[key] = score
    
    # get the max score position
    try:
        hp_pos = max(score_info, key=score_info.get)
    except ValueError:
        if adjust: 
            return -1, 0, "-", "-"
        else:
            print("BUG: Could not determine best position of error", file=sys.stderr)



    if hp_pos > 0:
        final_pos = max(hp_err_pos[hp_pos], key=hp_err_pos[hp_pos].get)     # get position of highest score within positions an hp region
    elif hp_pos == 0:
        try: 
            final_pos = max(hp_err_pos[hp_pos], key=hp_err_pos[hp_pos].get)     # try since if unique, there is no -0 so can't tell if a hp pos to consider
        except KeyError:
            final_pos = abs(hp_pos)
    else: 
        final_pos = abs(hp_pos)                                             # if negative, unique pos that is not a homopolymer. Get it directly.

    final_pos_erate = retrieve_pos_erate(err_rates, final_pos, adjust, find_low)

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
                    print("BUG: Could not determine replacement nt at position {}. DEBUG.".format(pos), file=sys.stderr)
                    rep_nt = "-"
        else:
            rep_nt = final_pos_erate.ins2.nt
    else:
        print("BUG: cannot determine error type to pull relevant error info at function determine_range_error_info. Debug", file=sys.stderr)

    if final_pos_erate.pos in existing_errs:
        return -1, final_pos_erate.depth, ref_nt, rep_nt
    else:
        return final_pos_erate.pos, final_pos_erate.depth, ref_nt, rep_nt

# updated version of set_erate_to_zero
# remove error rates that were used to find errors. Only unique ones
def set_erate_to_zero(err_rates, err_type, err_len, uniq_errs, orig_index):
    low = err_rates[orig_index].low
    for erate in uniq_errs:
        index = erate[0]
        erate_to_remove = (bin(erate[1]).count("1") - 1)/err_rates[abs(index)].depth
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
                        print("BUG: Could not find previously split indel of length 2 in function set_erate_to_zero", file=sys.stderr)
            else:
                pos_erate.ins2.remove_uniq_errs(erate[1])
                pos_erate.ins2.erate -= erate_to_remove
        else:
            print("BUG: error type could not be determined in function set_erate_to_zero", file=sys.stderr)

    del uniq_errs[:]
    return

# find errors
def find_errors(err_rates, err_list, nt_range, sub_mpile, last_round=False, existing_errs=[]):    
    erates_used = []    # keep track of erates used to find errors

    # loop through error list to find errors
    for i in range(0, len(err_rates)):
        # avoid calculation if too low depth
        if sub_mpile == "" and err_rates[i].low:
            continue
        elif err_rates[i].low and err_rates[i].depth < _depth_thres:
            continue
        elif err_rates[i].low:
            # adjust threshold for low depth sections substituted by a higher depth mpile
            # consider a different threshold for low depth regions
            thres = args.low_threshold
        else:
            thres = _indel_thres

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

        if nt_range < 0:
            # parse sub since it is the easiest and lower priority
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

        # get the best indel error 
        # bias towards single errors
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

        depth = err_rates[i].depth
        # if, else if: since there should never be both del and ins at the same time
        # else it would be a substitution 
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
                # set the error rate used in this evaluation to zero to ignore for future passes
                # set_erate_to_zero(err_rates, err_type, err_len, i, i+1)
                # keep track to set to zero
                uniq_align = get_uniq_err_rates(err_rates[i], err_type, err_len, read_align, uniq_pos_only=True)
                erates_used.append([i, uniq_align])
                set_erate_to_zero(err_rates, err_type, err_len, erates_used, i)
            else:
                continue
                print("BUG: Did not catch an indel error above threshold on first pass at position {}. DEBUG".format(i), file=sys.stderr)
        elif 0.25 <= indel_err < thres and (indel_err * depth) >= 2:
            # look for candidates
            # should not have any fully pass threshold error rates
            candidate_err = indel_err       # keep a record of the candidate error rate
            if nt_range < 0:
                continue
            elif nt_range == 0:       
                # look at hp only
                # not nt_range(0) should be True
                start = err_rates[i].hpStart
                end = start + err_rates[start].hpLen
            elif nt_range > 0:
                start = i - nt_range
                end = i + nt_range + 1
            else:
                pass

            # only run for every pass after first pass through error rates
            if last_round: 
                last_round_thres = thres
            else:
                last_round_thres = 0
            
            indel_err = get_best_indel_erate_in_range(err_rates, err_type, err_len, start, end, (not nt_range), last_round_thres, uniq_erates=erates_used, find=i)
            
            if indel_err > 1.2: 
                continue
            elif indel_err < 0:
                # found evidence for other error length. Flip length and indel error value.
                err_len -= (-1)**err_len
                indel_err = abs(indel_err)
            elif (not last_round or last_round < 0) and nt_range > 0 and indel_err < thres:
                # below 0 last round very last round: both window shifting and last_round calcs
                # found something nearby but not greater than threshold so lets move the window 
                eused2, eused3 = [], []

                # determine degree of window shifting
                left_start = start - int(nt_range/2); right_start = start + int(nt_range/2)
                left_end = end - int(nt_range/2); right_end = end + int(nt_range/2) 

                # determine left window
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

                # determine right window
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

                # change the starts and end accordingly
                left_end = right + 1; right_start = left
                left_start += 1; right_end -= 1

                # get the different window erates
                left_window_erate = get_best_indel_erate_in_range(err_rates, err_type, err_len, left_start, left_end, (not nt_range), last_round_thres, eused2, find=i)
                if left_window_erate >= 0:
                    # skip calculating right window if left_window already found a potential error
                    right_window_erate = get_best_indel_erate_in_range(err_rates, err_type, err_len, right_start, right_end, (not nt_range), last_round_thres, eused3, find=i)
                else:
                    right_window_erate = 0
                    
                if abs(left_window_erate) > indel_err or abs(right_window_erate) > indel_err:
                    if abs(left_window_erate) > abs(right_window_erate): 
                        indel_err = left_window_erate; start = left_start; end = left_end
                        del erates_used[:]
                        for e in eused2:
                            erates_used.append(e)
                    elif abs(right_window_erate) > abs(left_window_erate):
                        indel_err = right_window_erate; start = right_start; end = right_end
                        del erates_used[:]
                        for e in eused3:
                            erates_used.append(e)
                    else:
                        # if equal take the leftmost one
                        indel_err = left_window_erate; start = left_start; end = left_end
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
                # actually an indel error
                # avoid unnecessary looping unless there is an error to consider. But extra looping.
                    # Second possible conditional: relax threshold a bit for depths higher than expected max genomes mapped
                    # indicates that other parts of the genome are mapping to this area as well so affects positional error rate
                pos, depth, ref_nt, rep_nt = determine_range_error_info(err_rates, err_type, err_len, start, end, i,(not nt_range), adjust=False, existing_errs=existing_errs)
                
                if pos == -1:
                    # since only position is already considered as an error
                    # re-do position determination with different error length consideration
                    pos, depth, ref_nt, rep_nt = determine_range_error_info(err_rates, err_type, err_len, start, end, i, (not nt_range), adjust=False, existing_errs=existing_errs, check_diff_nlen=True)
                    if pos == -1:
                        # still no other candidate. Need to debug.
                        print("BUG: error rate determination is finding an error but no suitable candidate found near position {} within contig {}. DEBUG.".format(err_rates[i].pos, err_rates[i].contig), file=sys.stderr)
                        continue

                set_erate_to_zero(err_rates, err_type, err_len, erates_used, i)
            else:
                del erates_used[:]
                continue
        else:
            continue

        adjusted_i = i - (err_rates[i].pos-pos)         # adjust index if different position
        if pos in existing_errs:
            print("BUG: Collision. Two errors were determined to be on contig {} at position {}".format(err_rates[pos].contig, pos), file=sys.stderr)
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

# get initial Errors from reference and low depth mpileup files
def get_initial_errors(ref_mpile, low_mpile): 
    # get list of lists of positional error rates by contig
    err_rates = get_mpile_error_rates(ref_mpile, low_mpile)

    # determine errors through multiple passes through the error rate list
    # find errors through three passes with each successive pass being less strict
    # Pass 1: pos only (-), pass 2: hp only (0), pass 3: within 5 nt + window shift (1)
    # pass 4: within 10 nt + all indel length(2), pass 5+: window shift + all indel length
    # keep track of errors info and a list of position already determined as error in case of collisions
    found_err = []
    found_err_pos = []  
    for n in range(-1, args.passes):
        nt_range = n * 5
        if n > 2:
            find_errors(err_rates, found_err, nt_range, low_mpile, last_round=-1, existing_errs=found_err_pos)
        elif n > 0:
            find_errors(err_rates, found_err, nt_range, low_mpile, last_round=((n - 1) % 2), existing_errs=found_err_pos)
        else:
            find_errors(err_rates, found_err, nt_range, low_mpile, last_round=0, existing_errs=found_err_pos)

    # last ditch effort to find any residual errors by combining or splitting different n length errors
    # match last round but consider different nlength errors
    # find_errors(err_rates, found_err, (args.passes - 1) * 5, low_mpile, last_round = -1, existing_errs=found_err_pos)    

    # If list, sort
    found_err.sort(key = itemgetter(8,0))

    return found_err

# parse only necessary read info
# assumes that mpile file was in the same order as the other mpiles
# TODO: find a workaround to avoid using index 
def get_region_error_rates(mpile, error, pile_info, nt_range):
    # error information
    err_type = error[0]
    err_contig = error[1]
    err_pos = error[2]
    err_len = error[4]

    err_rates = []

    # set up for deletion error rates
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

        if contig != err_contig or not ((err_pos - end) <= pos <= (err_pos + end)):
            print("""WARNING: Either contig does not match or position {} is not within the range ({} to {}) being assessed.
                    Double check that all mpileup files are the same length and order""".format(pos, start, end), file=sys.stderr)

        pos_erate = Error(contig, pos, mpile_line[2], depth)
        
        # record actual sequence of this position to calculate if region is messy with indels
        if args.messy: pos_erate.coverage = depth

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

        # get error rate
        parse_pile(pos_erate, mpile_line[4], pile_info, get_only=error[0])
        
        # reset
        pile_info.reset_class()
        pile_info.set_next_del_to_zero()        

        # delete storage of next nt
        del pos_erate.next

        # store in growing error list
        err_rates.append(pos_erate)

        # NOTE: hp info on the edge of the region will likely be incorrect, but unlikely to use this edge info.
        # WARNING: No adjustments would be made to account for this edge info
        if end == 1:
            # substitution, no need for hp info. Return prematurely
            return err_rates
        elif index == (end - 1):
            fill_hp_info(err_rates, True)     
        else:
            fill_hp_info(err_rates)

        parse_del_only(pile_info, mpile_line[4])

    # reset variable for next error parsing
    return err_rates            

# get position with error rates greater than threshold
# present as possible adjustment position
def get_adjustment_candidates(error, err_rates, last_pos):
    err_type = error[0]
    err_len = error[4]
    err_candidates = [[],[]]
    start = -1  # initializing

    for index in range(len(err_rates)):
        if start == err_rates[index].hpStart: continue
        
        start = err_rates[index].hpStart
        end = start + err_rates[start].hpLen

        top_non_hp_err1 = []; top_non_hp_err2 = []
        # get best indel erate for both lengths
        err1 = get_best_indel_erate_in_range(err_rates, err_type, 1, start, end, True, uniq_erates=top_non_hp_err1, adjust=True)
        err2 = get_best_indel_erate_in_range(err_rates, err_type, 2, start, end, True, uniq_erates=top_non_hp_err2, adjust=True)

        # get all possible candidates for both length
        # Error length 1
        if top_non_hp_err1 and top_non_hp_err1[1] >= err1 and top_non_hp_err1[1] >= _adjust_thres:
            err_pos, depth, err_ref, err_rep = determine_range_error_info(err_rates, err_type, 1, top_non_hp_err1[0], top_non_hp_err1[0] + 1, hp=True, adjust=True)
            if err_pos == -1:
                # found no position candidate so redo with error of len two and just take the first nt
                err_pos, depth, err_ref, err_rep = determine_range_error_info(err_rates, err_type, 2, top_non_hp_err1[0], top_non_hp_err1[0] + 1, hp=True, adjust=True)

            # save as a candidate if it wasn't previous used as a candidate. 
            # Should minimally slow down script but needed for edge cases
            if err_pos >= 0 not in last_pos:
                err_candidates[0].append([err_type, error[1], err_pos, depth, 1, err_ref[0], err_rep[0], top_non_hp_err1[1]])
        if err1 >= _adjust_thres:
            err_pos, depth, err_ref, err_rep = determine_range_error_info(err_rates, err_type, 1, start, end, hp=True, adjust=True)
            if err_pos == -1:
                err_pos, depth, err_ref, err_rep = determine_range_error_info(err_rates, err_type, 2, start, end, hp=True, adjust=True)
            if err_pos >= 0 and err_pos not in last_pos:
                err_candidates[0].append([err_type, error[1], err_pos, depth, 1, err_ref[0], err_rep[0], err1])

        # Error length 2
        if top_non_hp_err2 and top_non_hp_err2[1] >= err2 and top_non_hp_err2[1] >= (_adjust_thres - 0.1):
            err_pos, depth, err_ref, err_rep = determine_range_error_info(err_rates, err_type, 2, top_non_hp_err2[0], top_non_hp_err2[0] + 1, hp=True, adjust=True)
            if err_pos >= 0 and err_pos not in last_pos:
                err_candidates[1].append([err_type, error[1], err_pos, depth, 2, err_ref, err_rep, top_non_hp_err2[1]])
        if err2 >= (_adjust_thres - 0.1):
            err_pos, depth, err_ref, err_rep = determine_range_error_info(err_rates, err_type, 2, start, end, hp=True, adjust=True)
            if err_pos >= 0 and err_pos not in last_pos:
                err_candidates[1].append([err_type, error[1], err_pos, depth, 2, err_ref, err_rep, err2])

    # return only the useful candidates
    if err_len == 1 and err_candidates[0]:
        return err_candidates[0]
    elif err_len == 2 and err_candidates[1]:
        return err_candidates[1]
    elif err_len == 2 and len(err_candidates[0]) >= 2:
        # for when there is no good evidence of two indel
        # but there is enough evidence for two separate indels within the range
        # will use the top two candidates as correction adjustments
        return err_candidates[0]
    else:
        # no candidates to consider
        return []

# Adjust error based on best candidate
def adjust_error(error, can, len_change=False):
    err_type = error[0]

    # change length if adjusting a two length error to two one length errors
    if len_change:
        error[4] = can[4]

    if error[2] == can[2]:
        # same position
        if error[6] != can[6]:
            # no agreement between reference and read based replacement nt
            # Take read correction error. Really only for insertions
            error[6] = can[6]
        return     
    else:
        # adjust position
        error[2] = can[2]

        # adjust ref or rep nt depending on error type
        if err_type == "del":
            error[5] = can[5]
        elif err_type == "ins":
            error[6] = can[6]
        else:
            print("BUG: Could not determine error type in function adjust_error. DEBUG.", file=sys.stderr)

    # add support number
    error[7] = ";".join([str(error[7]),str(can[7])])
    
    return

# retrieve line offsets so we can pull specific lines without loading entire file
# was lifted from StackOverflow by Adam Rosenfield
# faster than using file.tell()
def _load_read_file_offsets(file):
    line_offset = []
    offset = 0
    try:
        for line in file:
            line_offset.append(offset)
            offset += len(line)
    except UnicodeDecodeError:
        sys.exit("ERROR: Mpileup file is corrupted. Re-try with another mpileup file.")

    return line_offset

# get the chunks of the mpile file based on error position
def _retrieve_mpile_chunk(file, offset, region_size, chunk): 
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

# Read information based adjustments of error found
# Read information should supersede reference information within a certain threshold
def adjust_initial_error(reads_mpile, err_list, adjust_sub=True):
    global _adjustment_count
    global _supported_sub
    # get error rates
    # reads = _load_full_mpile_info(reads_mpile)
    reads = open(reads_mpile)
    offsets = _load_read_file_offsets(reads)
    mpile_chunk = []

    # store mpileup info
    container = MpileInfo()
    
    # create a list for extra errors (when adjusting for length differences)
    # extra list is needed so the already adjusted error does not undergo extra adjustments
    extra = []

    last_adjust_pos = [-1]
    for error in err_list:
        # check for sub support
        if error[0] == "sub":
            if adjust_sub:
                _retrieve_mpile_chunk(reads, offsets[error[8]], 1, mpile_chunk)
                err_rates = get_region_error_rates(mpile_chunk, error, container, 1)
                
                # retrieve error rate of respective substitution
                try:
                    poss_subs = list(err_rates[0].sub.nt.intersection(error[6]))
                except AttributeError:
                    poss_subs = []

                best_sub = "N"; best_sub_erate = 0
                for nt in poss_subs:
                    try:
                        sub_erate = (err_rates[0].subAllErr)[nt]
                    except AttributeError:
                        sub_erate = ({"A": 0, "C": 0, "G": 0,"T": 0})[nt]
                    # store the best sub_erate
                    if sub_erate > best_sub_erate:
                        best_sub_erate = sub_erate; best_nt = nt
                else:
                    # nothing to iterate so set to 0 
                    sub_erate = 0

                if best_sub_erate > 0.1:
                    # keep substitution. And replace the replacement nt with best substitution
                    _supported_sub += 1
                    error[6] = best_nt
                    error[7] = ";".join([str(error[7]),str(best_sub_erate)])
                    del error[-1]
                    # print("Substitution at pos {} has {} support with nucleotde {} from read data".format(error[2], sub_erate, error[6]))
                else:
                    # delete substitution error due to little evidence
                    del error[:]

            # next error because no other adjustments needed
            continue

        # window size
        nt_range = 10
        
        try:
            if error[0] == "del" and (error[8] - 1) > 0:
                _retrieve_mpile_chunk(reads, offsets[(error[8] - nt_range - 1)], (nt_range * 2 + 1), mpile_chunk)
            else:
                _retrieve_mpile_chunk(reads, offsets[(error[8] - nt_range)], (nt_range * 2), mpile_chunk)
        except IndexError:
            print("BUG: Could not retrieve relevant mpileup chunk between indices {} and {}.".format((error[8] - nt_range - 1), (nt_range * 2 + 1)), file=sys.stderr)

        err_rates = get_region_error_rates(mpile_chunk, error, container, nt_range)
        adj_can = get_adjustment_candidates(error, err_rates, last_adjust_pos)

        # adjust errors if possible
        if len(adj_can) == 1:
            _adjustment_count += 1
            print("Adjusting contig {} positon {}".format(error[1], error[2]), file=sys.stderr)
            # print("Adjusting error on contig {} at position {}".format(error[1], error[2]))
            # only one candidate
            adjust_error(error, adj_can[0])
            last_adjust_pos = [adj_can[0][2]]
        elif len(adj_can) > 1:
            print("Adjusting contig {} positon {}".format(error[1], error[2]), file=sys.stderr)
            _adjustment_count += 1
            # print("Adjusting error on contig {} at position {}".format(error[1], error[2]))
            # more than one candidate. Select the best one
            # if one of the candidates is the same as error, use that
            for index in range(len(adj_can)):
                if (adj_can[index][2] == error[2] or adj_can[index][2] == err_rates[err_rates[nt_range].hpStart].pos) and (adj_can[index][6][0] == error[6][0] or adj_can[index][6][-1] == error[6][-1]):
                    best_can = adj_can[index]
                    break
                else:
                    best_can = []
            
            if not best_can:
                best_can = max(adj_can, key=lambda k: k[7])
            
            last_adjust_pos = [best_can[2]]

            if error[4] != adj_can[0][4]:
                # first adjustment
                adjust_error(error, best_can, True)

                # trying to adjust a two length error with one
                # delete the old best candidate and then get best candidate for next adjustment
                to_del = next(index for index in range(len(adj_can)) if adj_can[index] == best_can)
                del adj_can[to_del]
                best_can = max(adj_can, key=lambda k: k[7])
                last_adjust_pos.append(best_can[2])
                
                # add index for further adjustment if not last round of adjustment
                if not adjust_sub: 
                    pos_diff = error[2] - best_can[2]
                    best_can.append((error[8] - pos_diff))

                # add second one length error to the end
                extra.append(best_can)
            else:
                adjust_error(error, best_can)
        else:
            error[7] = ";".join([str(error[7]),"0"])
            pass

        if adjust_sub:
            # delete ind x information as it's not needed anymore
            del error[-1]

    # add extra together 
    err_list += extra
    del extra[:]

    if adjust_sub:
        # remove deleted substitutions. Go backwards to not distrupt list order
        del_count = 0
        for i in range(len(err_list) -1, -1, -1):
            if err_list[i]:
                continue
            else:
                # no error in list. Remove
                del_count += 1
                del err_list[i]
        print("{} substiutions had little read support and were deleted".format(del_count))

    print("{} errors have been adjusted with {}".format(_adjustment_count, reads_mpile))
    _adjustment_count = 0 
    # close file 
    reads.close()

    return

# returns the first position with a differing nt on each side of the current position
def get_flank_nt_pos(gen, pos):
    index = pos                                                 #start with downstream
    pre_pos, suc_pos = pos, pos
    cur_nt = gen[pos]
    direction = 1
    while index < (len(gen) + 1) and index > -1:
        if gen[index] == cur_nt:
            index += direction
        else:
            if index > pos:
                # Downstream position was found. Store position. Reset position to search upstream.
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

#get the sequence and length of homopolymer. Direction determines whether to look upstream or downstream
def get_HP(gen, pos, direction):
    index = pos
    cur_nt = gen[pos]
    hp = ""
    while index < (len(gen) + 1) and index > -1:
        if gen[index] == cur_nt:
            index += direction
        else: 
            index += (direction * -1)   #gets last position where nt was the same if pre, get the position that nt was different if post/suc
            break

    #determine how to subset the genome
    if index < pos:
        hp = gen[index:(pos+1)]
    elif index > pos:
        hp = gen[pos:index+1]
    else:
        hp = gen[pos]

    return hp[0], len(hp)

# Get optional hp data 
def get_diagnostic_hp_data(out, draft, errors):
    # for multiple contigs
    dia_data = open("{}.diagnostic.txt".format(out), "w")
    dia_data.write("Contig\tPosition\tCurrentHp\tCurrentHpLen\tPreHp\tPreHpLen\tPostHp\tPostHpLen\tFor5mer\tRev5mer\tCreateHp\tAdjustHp\tBreakHp\n")
    for contig, err_list in errors.items():
        contig_seq = draft[contig]

        for err in err_list:
            err_type  = err[0]
            err_pos = int(err[2])
            err_len= err[4]

            # only change the central position for insertion
            # only change when inserting nt is next hp 
            # since only error type where nt doesn't involve reference
            if err_type == "ins":
                orig_pos = err_pos
                err_nt = err[6]
                try:
                    if err_nt[0] != contig_seq[err_pos] and err_nt[-1] == contig_seq[(err_pos + 1)]:
                        # next hp
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
            pre_nt, pre_nt_len = get_HP(contig_seq, pre_nt_pos, -1)
            suc_nt, suc_nt_len = get_HP(contig_seq, suc_nt_pos, 1)

            # adjust the error position to last deletion if deletion greater than one
            if err_type == "del" and err_len > 1:
                for_pos = err_pos
                rev_pos = err_pos + err_len
            elif err_type == "ins":
                for_pos = err_pos
                rev_pos = for_pos
            else:
                for_pos = err_pos
                rev_pos = err_pos + 1

            # with error position get the forward and reverse 5-mer 
            try:
                forward_5mer = contig_seq[(err_pos - 5):for_pos]
            except IndexError:
                forward_5mer = contig_seq[0:for_pos]

            try: 
                reverse_5mer = str(Seq(contig_seq[rev_pos:(rev_pos + 5)]).reverse_complement())
            except IndexError:
                reverse_5mer = str(Seq(contig_seq[rev_pos:len(contig_seq)]).reverse_complement())
            
            # break, adjust, create homopolymer information            
            if err_type == "ins":
                pre = contig_seq[orig_pos]
                post = contig_seq[orig_pos + 1]
            else:
                pre = contig_seq[err_pos - 1]
                post = contig_seq[rev_pos]

            create_hp = False
            adjust_hp = False
            break_hp = False        # only considered break if no flanking nt is similar
            # create and adjust
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
                if err_nt == pre_nt and err_nt == suc_nt and err_nt == pre and err_nt == post:
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
                elif (pre == cur_nt or post == cur_nt):
                    if  cur_nt_len == 2:
                        break_hp = True
                    else:
                        adjust_hp = True
                else:
                    pass

            lst_info = [contig, err[2], cur_nt, cur_nt_len, pre_nt[0], pre_nt_len, suc_nt[0], suc_nt_len, forward_5mer, reverse_5mer, create_hp, adjust_hp, break_hp]
            dia_data.write("%s\n" % "\t".join(map(str, lst_info)))

    return
    
# split the error list based on contig
def split_errors_by_contig(errors, keep=False):
    split_err = {}

    if keep:
        err_list = errors[:]            # pass by value
    else:
        err_list = errors               # pass by reference

    # assumes that error list was previously sorted by contig name, then position
    while len(err_list) != 0:
        current_contig = err_list[0][1]
        try:
            next_contig_index = next(index for index in range(len(err_list)) if err_list[index][1] != current_contig)
            split_err[current_contig] = err_list[0:next_contig_index]
            del err_list[0:next_contig_index]
        except StopIteration:
            split_err[current_contig] = err_list[0:len(err_list)]
            del err_list[0:len(err_list)]
            break

    return split_err

# correct the genome by contig
# splitting sequence to a list may be more memory intensive but better for deletion and insertion to not throw off position data
def correct_genome(draft, errors):
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
                    contig_seq[pos] = ""; contig_seq[pos+1] = ""
            elif err_type == "ins":
                contig_seq[pos] += correction[6]
            else:
                print("BUG: Could not determine error type at position {} for sequence correction".format(pos), file=sys.stderr)

            correction_count += 1

        # done for this contig
        draft[contig] = "".join(contig_seq)

    print("Made {} corrections!".format(correction_count))

# final output
def _write_genome_fasta(fasta_name, write_type, gen):
    # genome should be in a dictionary
    with open(fasta_name, write_type) as fasta:
        if type(gen).__name__ == "dict":
            for contig, seq in gen.items():
                fasta.write("%s\n" % ">{}".format(contig))
                if type(seq).__name__ == "list":
                    fasta.write("%s\n" % "".join(seq))
                elif type(seq).__name__ == "str":
                    fasta.write("%s\n"  % seq)
        else:
            print("BUG: Input type is not a dictionary. Could not output genome as fasta", file=sys.stderr)
    return

def main(args):
    # Load errors from an .err file 
    # Compute errors if not loading pre-computed errors
    if args.errors != "":
        errors = _load_prev_error(args.errors)
    else:
        # get reference based initial errors
        errors = get_initial_errors(args.ref, args.low)

        # TODO: set it as extraInfo
        _print_info("{}.02.err".format(args.out), "w", errors)

    if args.adjust < 2:
        # Make error position adjustments based on read information
        for file_num in range(len(args.reads)):
            # resort the errors
            errors.sort(key=itemgetter(1,2))            
            if file_num == (len(args.reads) - 1):
                adjust_initial_error(args.reads[file_num], errors, adjust_sub=True)
            else:
                adjust_initial_error(args.reads[file_num], errors, adjust_sub=False)

        # resort the errors
        errors.sort(key=itemgetter(1,2))

        # TODO: set it as extraInfo
        _print_info("{}.adjusted.02.err".format(args.out), "w", errors)

    # resort the errors, just in case
    errors.sort(key=itemgetter(1,2))

    # have the error list. Correct the genome
    draft = _read_draft()                                            # 1) load
    split_err = split_errors_by_contig(errors)                       # 1.5) split errors by contig 
    
    # write some diagnostic data for hp info of the reads
    get_diagnostic_hp_data(args.out, draft, split_err)

    correct_genome(draft, split_err)                                 # 2) correct
    _write_genome_fasta("{}.fa".format(args.out), "w", draft)        # 3) output


if __name__ == "__main__":
    # retrieve parameters
    args = _parse_Arguments()

    # Global thresholds to maintain
    _indel_thres = args.indel 
    _sub_thres = args.sub 
    _depth_thres = args.depth
    _adjust_thres = args.adjust
    _low_depth_count = 0
    _adjustment_count = 0
    _supported_sub = 0

    # TEMP
    errs_already_found = []

    # run main code
    main(args)

    # temporary 
    print("{} positions were considered low depth and not evaluated".format(_low_depth_count))
    # print("{} adjustments were made".format(_adjustment_count))
    print("{} substitutions have supporting read data".format(_supported_sub))
    print("DONE!")




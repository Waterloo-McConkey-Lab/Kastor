#!/usr/bin/python3
"""
Classes used in Castor.

Note
----
Binary alignments must have a 1 bit at left most position to maintain proper
bit length when bits are converted to 0 bit. I.e. a 64-bit can only store
the information of 63 aligned sequences. Alignment limitations is bypassed
by python's implementation of infinite bit integers

"""

import sys

# TODO: incorporate error objects rather than a raw list
class DetectedError:
    def __init__(self, contig, pos, err_type, depth, length, ref, alt,
                 score):
        self.contig = contig
        self.pos = pos
        self.err_type = err_type
        self.depth = depth
        self.length = length
        self.ref = ref
        self.alt = alt
        self.score = score

    def store_sub_alternatives(self, A=0, C=0, G=0, T=0):
        # input are the score for the respective nucleotide
        self.A = A
        self.C = C
        self.G = G
        self.T = T

        # initialize support information
        self.Asupp = 0
        self.Csupp = 0
        self.Gsupp = 0
        self.Tsupp = 0

    def set_single_sub(self, nt, score):
        if nt == "A":
            self.A = score
        elif nt == "C":
            self.C = score
        elif nt == "G":
            self.G = score
        elif nt == "T":
            self.T = score
        else:
            print("WARNING: Nucleotide does not exist. Nothing was done.")

    # TODO: finish
    def resolve_best_sub_alt(self):
        curr_best_nt = "N"
        curr_best_score = 0
        nt_dict = {"A": self.A, "C": self.C, "G": self.G, "T": self.T}

class PosInfo:
    """General information of an sequence position and any associated errors

    Note
    ----
    This is the parent class of ErrType child class. Specific error
    information will be done via this child class based on error type. To
    avoid initializing more classes than needed, they will be assigned on a
    need basis.

    Attributes
    ----------
    contig : str
        Contig name of the sequence
    pos : int
        Position in contig sequence
    ref : char
        Sequence nucleotide at current position
    depth : int
        Coverage depth of aligned reference genomes
    low : boolean
        Whether coverage depth of primary reference genome exceeds depth
        threshold
    erate : float
        Full error rate including all error types
    readEnd : int
        Binary alignment positions of reference genomes that are ending this
        position
    fullSub : float
        Error rate of substitution errors of all nucleotides
    fullIns : float
        Error rate of deletion errors (requiring insertions) of all lengths
    fullDel : float
        Error rate of insertion errors (requiring deletions) of all lengths

    """
    def __init__(self, contig, pos, ref, depth):
        """Initializes the class with general position information

        Parameters
        ----------
        contig : str
            Contig name of current sequence
        pos : int
            Position in current contig
        ref
            Sequence nucleotide at current position
        depth
            Coverage depth of aligned reference genomes
        """
        self.contig = contig
        self.pos = pos
        self.ref = ref
        self.depth = depth
        self.low = False
        self.erate = 0
        self.readEnd = 1
        self.fullSub = 0
        self.fullIns = 0
        self.fullDel = 0
        self.lineNum = None         # only for initial error detection

    def set_low_depth(self):
        """ Sets low depth indicator

        Returns
        -------
        None

        """
        self.low = True

    def record_read_end(self):
        """Changes last bit to indicate an aligned sequence ends

        Returns
        -------
        None

        """
        self.readEnd ^= 1

    def add_read_alignment(self):
        """Adds an additional sequence (as 0 bit) to the current alignment

        Returns
        -------
        None

        """
        self.readEnd <<= 1

    def get_read_positions(self):
        """Returns alignment information as a binary string

        Returns
        -------
        str
            Sequence alignments as binary

        """
        return bin(self.readEnd)


class ErrorRates:
    """Calculated error information and gen properties

    Note
    ----
    This is the child class of Error parent class. Slots was used to reduce
    memory usage; re-evaluate later to determine if still needed

    Attributes
    ----------
    nt : char
        Suggested nucleotide correction
    erate : float
        Error rate for the specified error and error length
    reads : int
        Binary alignment of sequences that indicate a potential error

    """
    __slots__ = ("nt", "erate", "reads")

    def __init__(self, nt, erate, reads):
        """Initializes the error class with the necessary error information

        Parameters
        ----------
        nt : char
            Suggested nucleotide correction
        erate : float
            Error rate for the specified error
        reads : int
            Binary alignment of sequences that indicate a potential error
        """
        self.nt = nt
        self.erate = erate
        self.reads = reads

    def reads_w_uniq_err(self, full_align):
        """Pull out unique sequence alignments without a previous error

        Sequences pulled out using this function contain an error at this
        specific position but did not contain an error prior to this
        position (within the bounds of the search window).

        Parameters
        ----------
        full_align : int
            Alignment of sequences that previously had an potential error
            detected. 1 bit indicates active sequence; 0 bit indicates that
            the aligned sequence already indicated an error of the same type
            and length.

        Returns
        -------
        int
            Unique alignment of sequences with a error

        """
        if self.reads.bit_length() != full_align.bit_length():
            print(("BUG: Alignment depth do not match between reads at this"
                   "position and full_align during class function "
                   "reads_w_uniq_err. DEBUG"), file=sys.stderr)
        return full_align & self.reads

    def subset_out_uniq_err(self, full_align):
        """Sets sequences with errors found this position as 0 bit

        Parameters
        ----------
        full_align : int
            Alignment of sequences that previously had an potential error
            detected. 1 bit indicates active sequence; 0 bit indicates that
            the aligned sequence already indicated an error of the same type
            and length.

        Returns
        -------
        int
            New alignment of sequences indicating an error

        """
        if self.reads.bit_length() != full_align.bit_length():
            print(("BUG: Alignment depth do not match between reads at this "
                  "position and full_align during class funciton "
                   "subset_out_uniq_err. DEBUG"), file=sys.stderr)

        # get only positions with unique errors
        uniq_err = full_align & self.reads
        uniq_err ^= (1 << (uniq_err.bit_length() - 1))

        return full_align ^ uniq_err

    def uniq_err_erate(self, full_align):
        """Calculate the error rate of sequences with an error this position

        Error rate is calculated by (number of errors) / (depth). Depth is
        considered the bit length - 1.

        Parameters
        ----------
        full_align : int
            Alignment of sequences that previously had an potential error
            detected. 1 bit indicates active sequence; 0 bit indicates that
            the aligned sequence already indicated an error of the same type
            and length.

        Returns
        -------
        float
            calculated error rate
        """
        return ((bin(full_align & self.reads).count("1") - 1)
                / (full_align.bit_length() - 1))

    def remove_uniq_errs(self, sub_align):
        """Sets sequences with errors as 'used' (0-bit)

        This function removes sequences that were 'used' to detect an error.
        Avoids re-detection of these errors in future passes through the
        data.

        Parameters
        ----------
        sub_align : int
            Alignment of sequences that were 'used' in error detection

        Returns
        -------
        None
            reads attributes are directly modified

        """
        if self.reads.bit_length() != sub_align.bit_length():
            print(("BUG: Alignment depth do not match between reads at this "
                  "position and sub_align during class function "
                   "remove_uniq_errs. DEBUG"), file=sys.stderr)
        elif sub_align == 0:
            print(("BUG: remove_uniq_errs class function was called when no"
                   "sequences are indicated to be set to 0. DEBUG"),
                  file=sys.stderr)

        sub_align ^= (1 << (sub_align.bit_length() - 1))
        self.reads ^= sub_align

    def adjust_indel_alignment(self):
        """Set last bit in sequence alignments as 1-bit

        Returns
        -------
        None

        """
        self.reads += 1

    def add_nonindel_alignment(self):
        """Add an additional sequence to the end of the alignment

        Returns
        -------
        None

        """
        self.reads <<= 1


class MpileInfo:
    """Class used to store information between independent positions

    This class is mostly used for storage for easier argument passing
    between functions as error rates are being calculated. It also avoids
    passing multiple lists during function calls as well as groups them up.

    Attributes
    ----------
    match : list
        a list of all the ASCII string that indicates an alignment match at
        the current position
    subStr : list
        a list of all the suggested substitutions at the current position
    insStr : list
        a list of all the suggested insertions at the current position
    prevDel : list
        a list of all suggested deletions starting at the current position
        Information was parsed from the previous position due to how mpileup
        files are structured
    skip : int
        a counter to parse alignment string
    delReads : int
        Alignment of all sequences that suggests a deletion at the parsed
        position
    lowErates : list
        a list that stores, up to n elements, the calculated error rates
        from the supplementary genome alignments until a low depth region
        is reached. These stored error rates allows for the extension of
        low depth boundaries to avoid issues with determining errors at low
        depth boundaries
    subPrevDel : list
        delStr list for the supplementary alignment files
    subDelReads : int
        delReads for the supplementary alignment files
    lowCounter : boolean
        Indicates if a low depth region was reached

    """
    def __init__(self):
        """Initializes the class

        See class docstrings for all attributes

        """
        self.match = []
        self.subStr = []
        self.insStr = []
        self.prevDel = []
        self.skip = 0
        self.delReads = 1   # keep track of which reads have a deletion
        self.lowErates = []
        self.subPrevDel = []
        self.subDelReads = 1
        self.lowCounter = False
        self.lineNum = -1    # keep track of

    def reset_class(self):
        """Resets some class attributes

        This function returns the respective attributes to original state to
        avoid retaining previous position's information

        Returns
        -------
        None

        """
        del self.match[:]
        del self.subStr[:]
        del self.insStr[:]
        self.skip = 0 

    def append_prev_del_info(self, nt, low=False):
        """Retain deletion information for next position

        Parameters
        ----------
        nt : char
            suggested deletion
        low : boolean
            Indicates which attribute to store the deletion info
        Returns
        -------

        """
        if low:
            self.subPrevDel.append(nt)
        else:
            self.prevDel.append(nt)

    def set_next_del_to_zero(self, low=False):
        """Resets deletion information to original state

        This function returns the respective attributes to original state to
        avoid retaining previous position's information. Specific to
        deletions since that is the only information needed to be retained
        between positions

        Parameters
        ----------
        low : boolean
            Indicates which deletion attributes to reset

        Returns
        -------
        None

        """
        if low:
            self.subDelReads = 1
            del self.subPrevDel[:]
        else:
            self.delReads = 1
            del self.prevDel[:]

    def record_prev_del_read(self, low=False):
        """Adds a 0-bit to deletion alignment

        Parameters
        ----------
        low : boolean
            Indicates which attribute to adjust

        Returns
        -------
        None

        """
        if low:
            self.subDelReads <<= 1
        else:
            self.delReads <<= 1
            
    def remove_prev_del_read(self, low=False):
        """Removes the last bit in the deletion alignment

        Parameters
        ----------
        low : boolean
            Indicates which attribute to adjust

        Returns
        -------
        None

        """
        if low:
            self.subDelReads >>= 1
        else:
            self.delReads >>= 1

    def adjust_prev_read_as_del(self, low=False):
        """Flip last bit in deletion alignment

        Parameters
        ----------
        low : boolean
            Indicates which attribute to adjust

        Returns
        -------
        None

        """
        if low:
            self.subDelReads ^= 1
        else:
            self.delReads ^= 1

    def set_low_counter(self):
        """Set low depth counter

        Returns
        -------
        None

        """
        self.lowCounter = True

    def save_low_erate(self, window, pos_erate):
        """Saves supplementary alignment error rate

        Adds the latest calculated error rate from the supplementary
        alignment file to the saved list in case of nearby low depth
        regions. Only the last n elements are retained where n is the window
        size

        Parameters
        ----------
        window : int
            The number of elements to retain in the list
        pos_erate : Error class
            Calculated error rates from the supplementary alignment file

        Returns
        -------
        None

        """

        if len(self.lowErates) > window:
            self.lowErates.append(pos_erate)
            del self.lowErates[0]
        else:
            self.lowErates.append(pos_erate)

    def add_erate_to_err_rates(self, err_list):
        """Dumps all saved supplementary error information to error list

        This class function should only be called during the initial
        encounter of a low depth region. The information is added to
        the Error classes in upstream position under lowErate to avoid
        collision with the primary error information. The function iterates
        backwards (upstream) since there is no set size.

        Parameters
        ----------
        err_list : list
            passed error list containing calculated error rates of previous
            positions

        Returns
        -------
        None
            the error rate list is modified directly

        """
        for n in range(-1, (-1 * len(self.lowErates)) - 1, -1):
            try:
                if hasattr(err_list[n], "lowErate"):
                    # already have an attribute named lowErate. Collision
                    break
                else:
                    err_list[n].lowErate = self.lowErates[n]    
            except IndexError:
                break
            
    def reset_low_erate(self):
        """Resets any low depth information that was previously calculated

        Returns
        -------

        """
        self.lowCounter = False
        del self.lowErates[:]

    def increment_line_num(self):
        self.lineNum += 1
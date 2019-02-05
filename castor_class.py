#!/usr/bin/python3
import pdb
# Data structures for ref_correct.py

# keep track of all errors in a position
class Error:    
    def __init__(self, contig, pos, ref, depth):
        self.contig = contig
        self.pos = pos
        self.ref = ref
        self.depth = depth
        self.low = False
        self.erate = 0
        self.readEnd = 1        # keeps track of which reads end
        self.fullSub = 0
        self.fullIns = 0
        self.fullDel = 0

    def set_low_depth(self):
        self.low = True

    def record_read_end(self):
        self.readEnd ^= 1

    def add_read_alignment(self):
        self.readEnd <<= 1

    def get_read_positions(self):
        return bin(self.readEnd)

# specific error
class ErrType:
    __slots__ = ("nt", "erate", "reads")
    def __init__(self, nt, erate, reads):
        self.nt = nt
        self.erate = erate      # subsetted erate
        self.reads = reads      # where the error containing reads are aligned

    # subset reads with unique errors in a nt region
    def reads_w_uniq_err(self, full_align):
        if self.reads.bit_length() != full_align.bit_length(): pdb.set_trace()
        return (full_align & self.reads)

    def subset_out_uniq_err(self, full_align):
        if self.reads.bit_length() != full_align.bit_length(): pdb.set_trace()

        # get only positions with unique errors
        uniq_err = full_align & self.reads
        uniq_err ^= (1 << (uniq_err.bit_length() - 1))

        return (full_align ^ uniq_err)

    def uniq_err_erate(self, full_align):
        return (bin(full_align & self.reads).count("1") - 1)/(full_align.bit_length() - 1)

    def remove_uniq_errs(self, sub_align):
        if self.reads.bit_length() != sub_align.bit_length(): pdb.set_trace()
        if sub_align == 0: pdb.set_trace()
        sub_align ^= (1 << (sub_align.bit_length() - 1))
        self.reads ^= sub_align

    def adjust_indel_alignment(self):
        self.reads += 1

    def add_nonindel_alignment(self):
        self.reads <<= 1

# class to hold the list involved with parsing mpileup lines
# Avoids passing multiple lists during function calls and groups them up nicely
class MpileInfo:
    def __init__(self):
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

    def reset_class(self):
        del self.match[:]
        del self.subStr[:]
        del self.insStr[:]
        self.skip = 0 

    def append_prev_del_info(self, nt, low=False):
        if low:
            self.subPrevDel.append(nt)
        else:
            self.prevDel.append(nt)

    def set_next_del_to_zero(self, low=False):
        if low:
            self.subDelReads = 1
            del self.subPrevDel[:]
        else:
            self.delReads = 1
            del self.prevDel[:]

    def record_prev_del_read(self, low=False):
        if low:
            self.subDelReads <<= 1
        else:
            self.delReads <<= 1
            
    def remove_prev_del_read(self, low=False):
        if low:
            self.subDelReads >>= 1
        else:
            self.delReads >>= 1

    def adjust_prev_read_as_del(self, low=False):
        if low:
            self.subDelReads ^= 1
        else:
            self.delReads ^= 1

    def set_low_counter(self):
        self.lowCounter = True

    def save_low_erate(self, window, pos_erate):
        # only keep the last n (window size) incase of low erate
        if len(self.lowErates) > (window):
            self.lowErates.append(pos_erate)
            del self.lowErates[0]
        else:
            self.lowErates.append(pos_erate)

    def add_erate_to_err_rates(self, err_list):
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
        self.lowCounter = False
        del self.lowErates[:]

        

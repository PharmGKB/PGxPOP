import numpy as np
import tabix
import os

def get_competitive_haplotype_indices(haplotype_matrix):
    dual_sites = np.where(np.sum(haplotype_matrix, axis=0) > 1)[0]
    competitive_haplotypes = dict()
    for x in haplotype_matrix[:,dual_sites].T:
        comps = set(np.where(x == 1)[0])
        for y in comps:
            if y not in competitive_haplotypes:
                competitive_haplotypes[y] = set()
            competitive_haplotypes[y].update(comps.difference({y}))
    return(competitive_haplotypes)


"""
Normalize chromosome names.  Sometimes chromosomes are listed as just numbers, other times they start with 'chr'. 
This function strips the 'chr' and returns just the chromosome number, if applicable.

Inputs: 
    chr: variable containing chromosome id
Returns:
    chr with 'chr' stripped from beginning.

"""
def clean_chr(chr):
    if isinstance(chr, int):
        return chr
    if chr.startswith('chr'):
        return chr.replace("chr", "")
    return chr

'''
This function parses a line from a VCF into an object with some useful functions.

Inputs:
    line: A string directly from a VCF representing an entire row of call data.
Returns:
    VCFfields object
'''
def parse_vcf_line(line):
    CHROM = 0
    POS = 1
    ID = 2
    REF = 3
    ALT = 4
    QUAL = 5
    FILTER = 6
    INFO = 7
    FORMAT = 8
    calls = 9

    if isinstance(line, list) is True:
        fields = line
    else:
        fields = line.rstrip().split()

    class VCFfields(object):
        def __init__(self):
            self.chrom = None
            self.pos = None
            self.id = None
            self.ref = None
            self.alt = None
            self.qual = None
            self.filter = None
            self.info = {}
            self.format = None
            self.calls = []
            self.gt_index = None

        def print_row(self, chr=True, minimal=False):
            info_print_format = self.format_info()
            calls_print_format = "\t".join(self.calls)
            if chr:
                # If I want to make sure the output has chr on the chromosome, do this.  Kind of messy but works.
                chrom = "chr%s" % clean_chr(self.chrom)
            else:
                chrom = self.chrom
            if minimal is True:
                print("%s\t%s\t%s\t%s\t%s" % (chrom, self.pos, self.id, self.ref, self.alt))
            else:
                print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom, self.pos, self.id, self.ref, self.alt,
                                                                  self.qual, self.filter, info_print_format,
                                                                  self.format, calls_print_format))
        def intake_info(self, info):
            fields = info.split(';')
            for f in fields:
                try:
                    key, value = f.split("=")
                    self.info[key] = value
                except:
                    pass

        def format_info(self):
            output = []
            if len(self.info) == 0:
                return "."
            for key, value in self.info.items():
                output.append("%s=%s" % (key, value))


            return ";".join(output)

        def add_info(self, key, value):
            self.info[key] = value

        def add_allele(self, new):
            if new == self.alt:
                return
            self.alt += ",%s" % new

        def make_multiallelic(self, new_vcf):
            self.add_allele(new_vcf.alt)
            # add the info
            for i in new_vcf.info:
                self.add_info(i, new_vcf.info[i])
            # update the genotypes
            for i in range(len(new_vcf.calls)):
                if "/" in new_vcf.calls[i]:
                    alleles = new_vcf.calls[i].split("/")
                elif "|" in new_vcf.calls[i]:
                    alleles = new_vcf.calls[i].split("|")
                else:
                    print("Unrecognized delimiter!")
                    print(new_vcf.calls[i])
                    exit()

                alt_found = False

                n_alts = str(len(self.alts()))

                if alleles[0] == '1':
                    alleles[0] = n_alts
                    alt_found = True
                if alleles[1] == '1':
                    alleles[1] = n_alts
                    alt_found = True
                if alt_found:
                    self.calls[i] = "/".join(alleles)

        def alts(self):
            return self.alt.split(",")

        def is_multiallelic(self):
            if len(self.alts()) == 1:
                return False
            return True

        def max_insertion(self):
            max_insert = 1
            for a in self.alts():
                if len(a) > max_insert:
                    max_insert = len(a)
            return max_insert

        def is_indel(self):
            if len(self.ref) > 1 or len(self.alt) > 1:
                return True
            return False

    row = VCFfields()
    row.chrom = fields[CHROM]
    row.pos = int(fields[POS])
    row.id = fields[ID]
    row.ref = fields[REF]
    row.alt = fields[ALT]
    row.qual = fields[QUAL]
    row.filter = fields[FILTER]
    row.intake_info(fields[INFO])
    if len(fields) > 8:
        row.format = fields[FORMAT]
        row.calls = fields[calls:]
    return row

'''
Read a VCF and return the subject IDs
'''
def get_vcf_subject_ids(vcf):
    ids = []
    with smart_open(vcf) as f:
        for line in f:
            try:
                line = byte_decoder(line)
            except:
                line = line
            if line.startswith("#CHROM"):
                fields = line.rstrip().split()
                ids = fields[9:]
                break
    return ids


'''
Detect if a file is gzipped before opening.  If it is open with gzip, otherwise open normally.
'''
def smart_open(file):
    if file.endswith('gz'):
        import gzip
        return gzip.open(file)
    return open(file)

'''
Decode byte data into utf-8
'''
def byte_decoder(a):
    return a.decode("utf-8")


'''
Split a genotype string into an array
'''
def split_genotype(gt):
    gt = gt.split(":")[0]
    if "/" in gt:
        alleles = gt.split("/")
        return alleles
    elif "|" in gt:
        alleles = gt.split("|")
        return alleles
    print("Unrecognized delimiter! %s" % gt)
    return gt

'''
Check if a file is phased
'''
def vcf_is_phased(vcf):
    with smart_open(vcf) as f:
        for line in f:
            try:
                line = byte_decoder(line)
            except:
                line = line
            if line.startswith("#"):
                continue
            vcfline = parse_vcf_line(line)

            complete = False
            index = 0
            while not complete:
                gt = vcfline.calls[0].split(':')[index]
                if "|" in gt:
                    f.close()
                    return True
                elif "/" in gt:
                    f.close()
                    return False
                elif index > len(vcfline.calls):
                    complete = True
                else:
                    index += 1

        # If we went through all the genotypes and haven't figure it out, there is something wrong.
        print("Phasing status could not be determined.  VCF may be corrupted.")
        exit(1)

def is_gt_phased(gt):
    if "/" in gt:
        return False
    if "|" in gt:
        return True
    # otherwise I don't know, so we'll say false
    return False

def chr_string(vcf):
    with smart_open(vcf) as f:
        for line in f:
            try:
                line = byte_decoder(line)
            except:
                line = line
            if line.startswith("#"):
                continue
            if line.startswith("chr"):
                f.close()
                return True
            else:
                f.close()
                return False


def fetch_genotypes(vcf, variant, synonym=None):
    tb = tabix.open(vcf)

    # Check whether the chromosomes start with 'chr' or not
    chr_represenation = chr_string(vcf)
    if chr_represenation:
        chr = variant.chromosome
    else:
        chr = variant.clean_chromosome

    # synonym is just an alternate position that the variant may be at.  This is common for INDELs.
    if synonym is not None:
        position = synonym

    else:
        position = variant.position

    records = tb.query("%s" % chr, position-1, position+1)
    #records = tb.query("%s" % chr, position - 1, position + 1)

    n_records = 0
    for r in records:
        n_records += 1

    if n_records > 1:
        strict = True
    else:
        strict = False

    records = tb.query("%s" % chr, position - 1, position + 1)

    for r in records:
        #print("Record!")
        # todo add some additional checks here
        # Otherwise check if there is an rsid, if there is check that it matches the variant file
        # If not just check the position and the alternate alleles
        # If it is an indel, try to match but it might not

        # Check that we fetched the right thing
        valid = True

        # gt_index will store the index of the matched alternate allele
        gt_index = None

        if r[0] != chr:
            #print("Chromosome didn't match")
            valid = False

            if variant.rsid != r[2] and strict is True:
                continue

        if str(r[1]) != str(position):
            valid = False
            #if variant.rsid != r[2] and strict is True:
            #    continue

        if r[3] != variant.ref:
            #print("Ref didnt' match")
            valid = False
            #if variant.rsid != r[2] and strict is True:
            #    continue

        #if variant.position == 42130655:
        #    print("*15!!!!")

        any_alt = False
        #print(r[4])

        found_alts = r[4].split(",")

        for alt in found_alts:
            if alt in variant.alt:
                #gt_index = variant.alt.index(alt)
                gt_index = found_alts.index(alt)
                #print("GT index: %s " % gt_index)
                any_alt = True

                if r[3] == variant.ref:
                    valid = True


        #if any_alt is False:
        #    print("Alt didn't match")



        if valid is False:
            # Check if it's an INDEL
            # todo put this in it's own function
            for a in variant.alt:
                # Check the case where the alt is a deletion
                if a.startswith("del"):
                    # Check if any alt listed satisfies the deletion
                    satisfied, index = del_checker(a, r)
                    if satisfied is True:
                        valid = True
                        gt_index = index

                # Check the case where it is an insertion
                elif len(a) > 1:
                    satisfied, index = ins_checker(a, r)
                    if satisfied is True:
                        valid = True
                        gt_index = index

            # I know this is duplicatd but it works.  The last one wasn't catching cases where a single nucleotide was
            # added
            if variant.type == "INS":
                for a in variant.alt:
                    satisfied, index = ins_checker(a, r)
                    if satisfied is True:
                        valid = True
                        gt_index = index


                # todo add other exceptions where necessary

        if variant.rsid == r[2]:
            # if the IDs match then nothing else matters
            # This will especially catch INDELs with different conventions for ref and alt representation
            #if valid == False:
            valid = True

        if valid is True:
            data = parse_vcf_line(r)
            # If there is only one alternate allele in the definition, set the index
            if len(variant.alt) == 1:
                data.gt_index = gt_index
            return data

    # Recursively try any synonyms - not actually recursive
    if synonym is not None:
        for pos in variant.synonyms:
            syn_result = fetch_genotypes(vcf, variant, pos)
            if syn_result is not None:
                return syn_result


    return None

def ins_checker(original, query):
    # I'm not sure this will work for all insertion representations
    # Basically what I'm doing is adding the added nucleotides from the definition to the ref allele and checking
    # If that matches the listed alternate.
    # I'm sure there will be exceptions to take care of at some point.
    index = 0
    for alt in query[4].split(","):
        if alt == query[3] + original:
            return True, index
        index += 1
    return False, None

def del_checker(original, query):
    # This checks the difference between the ref and the alt and checks if that matches the deletion stated in the
    # definition file.
    deleted_nts = original.strip("del")
    ref = query[3]
    index = 0
    for alt in query[4].split(","):
        if alt in ref:
            diff_nts = ref.replace(alt, '', 1)
            if diff_nts == deleted_nts:
                return True, index
        index += 1

        #r = difflib.ndiff(query[3], alt)
        #for i in r:
        #    print(i)
        #try:
        #    i, s, = difflib.ndiff(query[3], alt)
        #    print(i, s)
        #    if s == "- %s" % deleted_nts:
        #        # Deletion resolved
        #        return True
        #except:
        #    print("Failed")
        #    continue
    return False, None

def fetch_genotype_records(vcf, chromosome, position):
    tb = tabix.open(vcf)
    records = tb.query("%s" % chromosome, position - 1, position + 1)
    for r in records:
        # todo add some additional checks here
        # Otherwise check if there is an rsid, if there is check that it matches the variant file
        # If not just check the position and the alternate alleles
        # If it is an indel, try to match but it might not
        if r[1] == str(position):
            data = parse_vcf_line(r)
            return data
    return None


# Sometimes multiple nucleotides ban be found at a position and have the same function.  In this case
# rather than representing them as multiple letters they can be collapsed into a single letter.  For example
# R is used to represent any purine, so the site could be an A or a G.  This function will map a nucleotide to
# it's IUPAC definition and return a list of nucleotides it could be.
def iupac_nt(nt):
    nts = {
        "A": ["A"],
        "C": ["C"],
        "T": ["T"],
        "G": ["G"],
        "W": ["A", "T"],
        "S": ["C", "G"],
        "M": ["A", "C"],
        "K": ["G", "T"],
        "R": ["A", "G"],
        "Y": ["C", "T"],
        "B": ["C", "G", "T"],
        "D": ["A", "G", "T"],
        "H": ["A", "C", "T"],
        "V": ["A", "C", "G"],
        "N": ["A", "C", "G", "T"],
        "Z": []
    }

    if nt in nts.keys():
        return nts[nt]

    #print("%s not found in nucleotide definition!" % nt)
    return [nt]

pass


def get_definition_file(g):
    definition_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../definition/alleles/")
    filename = "%s_translation.json" % g
    definition_file = os.path.join(definition_dir, filename)
    return definition_file


def welcome_message():
    message = '''
                                                                D\___/\  
                                                                 (0_o)    
                                                                  (V)    
     _________________________________________________________oOo__U__oOo____________
    |                      ___ _ _        ___                                        |
    |        _       _    / __(_) |_ _  _|   \ __ ___ __ ____ _    _       _         | 
    |       (_'-----'_)  | (__| |  _| || | |) / _` \ V  V / _` |  (_'-----'_)        |
    |       (_.'""""._)   \___|_|\__|\_, |___/\__,_|\_/\_/\__, |  (_.'""""._)        |
    |                                |__/                 |___/                      |
    |                                                                                |
    |               v0.0 (pre-pre Alpha) (seriously don't use this)                  |
    |     Written by Adam Lavertu and Greg McInnes with help from PharmGKB.          |
    |  We do not certify that the output from this program are correct in any way.   |
    |                                                                                |
    |________________________________________________________________________________|
    '''

    print(message)
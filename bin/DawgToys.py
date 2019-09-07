import numpy as np
import tabix

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
                    return True
                elif "/" in gt:
                    return False
                elif index > len(vcfline.calls):
                    complete = True
                else:
                    index += 1

        # If we went through all the genotypes and haven't figure it out, there is something wrong.
        print("Phasing status could not be determined.  VCF may be corrupted.")
        exit(1)

def fetch_genotypes(vcf, variant):
    tb = tabix.open(vcf)
    records = tb.query("%s" % variant.chromosome, variant.position - 1, variant.position + 1)
    for r in records:
        # todo add some additional checks here
        # Otherwise check if there is an rsid, if there is check that it matches the variant file
        # If not just check the position and the alternate alleles
        # If it is an indel, try to match but it might not
        if r[1] == str(variant.position):
            data = parse_vcf_line(r)
            return data
    return None

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


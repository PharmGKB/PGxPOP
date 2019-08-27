"""
Greg McInes
Altman Lab
gmcinnes@stanford.edu

A set of commands to work with sequencing data
"""

from __future__ import division
from glib import *
import sys

# Calculate GC content for specified genomic window
def window_gc_content(fasta_file, chr, start, end):
    sequence = fasta_extract(fasta_file, chr, start, end)
    gc = gc_content(sequence)
    return gc

# Calculate the GC content of a sequence of DNA
def gc_content(sequence):
    # Make the entire sequence uppercase
    sequence = sequence.upper()
    # Get the necessary counts
    seq_length = len(sequence)
    if seq_length == 0:
        return 0
    g_count = sequence.count("G")
    c_count = sequence.count("C")

    # Calculate gc content.  Sum of G's and C's divided by sequence length
    gc_content = (g_count + c_count) / seq_length
    return gc_content

# Extract sequence from a fasta file using samtools
def fasta_extract(fasta_file, chr, start, end):
    # Check if samtools is available
    if not which('samtools'):
        print("Samtools not found on path.  Please install to continue.")
        exit(1)
    # Format the region

    #region = "%s:%s-%s" % (clean_chr(chr), start, end)
    region = "%s:%s-%s" % (chr, start, end)
    #print(region)
    #print(fasta_file)
    # Format and execute the command
    command = 'samtools faidx %s %s' % (fasta_file, region)
    output = run_system_cmd(command)
    # Format the response
    # The returned sequence will have a header line (>chr:start-end) followed by newline separated sequence
    seq_list = output.split(b"\n")
    seq = ''
    for s in seq_list[1:]:
        seq += s.decode("utf-8")
    return seq


# Check whether a chromosome is greater than or less than another chromosome
def chr_greater_than_chr(chr1, chr2):

    chr1 = str(chr1)
    chr2 = str(chr2)

    chr1 = strip_chr(chr1)
    chr2 = strip_chr(chr2)

    # Check whether they're the same
    if chr1 == chr2:
        return False

    # Check the non numeric chromosomes.  I know there are more but I'm not accounting for those
    if chr1.lower() == 'x':
        if chr2.lower() == 'y' or chr2.lower() == 'm':
            return False
        return True

    if chr1.lower() == 'y':
        if chr2.lower() == 'x':
            return True

        if chr2.lower() == 'm':
            return False

    if chr2.lower() in ['x', 'y', 'm']:
        return False

    try:
        chr1 = int(chr1)
        chr2 = int(chr2)

        return chr1 > chr2

    except Exception as e:
        print(e)
        exit()

def get_vcf_subject_ids(vcf):
    ids = []
    with g_open(vcf) as f:
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

#todo get subjects and make a dictionary
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

        def print_row(self, chr=True):
            info_print_format = self.format_info()
            calls_print_format = "\t".join(self.calls)
            if chr:
                # If I want to make sure the output has chr on the chromosome, do this.  Kind of messy but works.
                chrom = "chr%s" % strip_chr(self.chrom)
            else:
                chrom = self.chrom
            print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom, self.pos, self.id, self.ref, self.alt,
                                                              self.qual, self.filter, info_print_format, self.format,
                                                              calls_print_format))
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
                #rint(n_alts)

                if alleles[0] == '1':
                    #alleles[0] = str(len(alleles))
                    alleles[0] = n_alts
                    alt_found = True
                if alleles[1] == '1':
                    #alleles[1] = str(len(alleles))
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
                #if len(a) - len(self.ref) > max_insert:
                if len(a) > max_insert:
                    max_insert = len(a) #- len(self.ref)
            return max_insert

        def max_insertion2(self):
            max_insert = 0
            for a in self.alts():
                print(a, file=sys.stderr)
                print(len(a) - len(self.ref), file=sys.stderr)
                if len(a) - len(self.ref) > max_insert:
                    max_insert = len(a) - len(self.ref)
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


#
def bed_extract(bed_file):
    class BedFields(object):
        def __init__(self):
            self.entries = {
                "chrom":[],
                "start":[],
                "end":[]
            }
            self.count = 0
        def add_entry(self, chrom, start, end):
            self.entries["chrom"].append(chrom)
            self.entries["start"].append(int(start))
            self.entries["end"].append(int(end))
            self.count += 1
        def retrieve_entry(self, index):
            # todo, check for valid index
            return self.entries['chrom'][index], self.entries['start'][index], self.entries['end'][index]
    bed_result = BedFields()
    with open(bed_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split()
            bed_result.add_entry(chrom=fields[0], start=fields[1], end=fields[2])
    return bed_result

    def check_range(self, bed, vcf):
        if chr_greater_than_chr(bed[0], vcf[0]):
            return False
        if bed[1] > vcf[1]:
            return False
        return True









def strip_chr(chrom):
    if isinstance(chrom, int):
        return chrom
    if chrom.startswith('chr'):
        return chrom.replace("chr", "")
    return chrom

# Write functions to check if an allele is a pyrimadine or purine
def is_purine(allele):
    purines = ['A', 'G']
    if allele in purines:
        return True
    return False

def is_pyrimadine(allele):
    pyrimadine = ['C', 'T']
    if allele in pyrimadine:
        return True
    return False

# Write functions to check if a variant is a transition or transversion
def is_transition(allele_1, allele_2):
    if (is_purine(allele_1) and is_purine(allele_2)) or (is_pyrimadine(allele_1) and is_pyrimadine(allele_2)):
        return True
    return False

def is_transversion(allele_1, allele_2):
    if (is_purine(allele_1) and is_pyrimadine(allele_2)) or (is_pyrimadine(allele_1) and is_purine(allele_2)):
        return True
    return False

# Check if a variant is an indel
def is_indel(ref, alt):
    if len(ref) != 1 or len(alt) != 1:
        return True
    return False

def indel_size(ref, alt):
    return abs(len(ref) - len(alt))

def is_frameshift(ref, alt):
    if indel_size(ref, alt) % 3 == 0:
        return False
    return True


def high_at(fasta_file, chr, start, end):
    gc = window_gc_content(fasta_file, chr, start, end)
    if gc < 0.25:
        return True
    return False

def high_gc(fasta_file, chr, start, end):
    gc = window_gc_content(fasta_file, chr, start, end)
    if gc > 0.75:
        return True
    return False

# Allele frequency calculations
def is_commom(af):
    if af >= 0.01:
        return True
    return False

def is_rare(af, include_very_rare=True):
    if af == ".":
        return True
    af = float(af)
    # Cutoff for rare from GoT2D paper
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5034897/
    if include_very_rare and af < 0.005:
        return True
    if af >= 0.00001 and af < 0.005:
        return True
    return False

def is_very_rare(af):
    if af == ".":
        return True
    if af < 0.00001:
        return True
    return False

# Conservation
def is_conserved(phylop_score):
    # Super simple threshold here.  Could definitely do better if I did some reading
    if phylop_score > 0.5:
        return True
    return False

def is_dann_pathogenic(dann):
    # Super rough eyeballed cutoff based on this
    # http://www.enlis.com/blog/2015/03/17/the-best-variant-prediction-method-that-no-one-is-using/
    try:
        if float(dann) > 0.9688:
            return True
        else:
            return False
    except:
        return False

def is_fathmm_pathogenic(fathmm):
    # Super rough eyeballed cutoff based on this
    # http://www.enlis.com/blog/2015/03/17/the-best-variant-prediction-method-that-no-one-is-using/
    #if fathmm > 0.8:
    #    return True
    #return False
    if fathmm > 0.8:
        return True
    return False

def is_cadd_pathogenic(cadd):
    # http://cadd.gs.washington.edu/info
    # recommended cutoff
    if cadd > 19.19:
        return True
    return False

def is_gerp_pathogenic(gerp):
    # cutoff from here
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4460046/
    if gerp >= 2:
        return True
    return False

def vcf_info_to_dict(info):
    dict = {}
    info_fields = info.split(';')
    for f in info_fields:
        field = f.split("=")
        if len(field) == 1:
            continue
        dict[field[0]] = field[1]
    return dict

# Functional annotations
def is_exonic(vcf_info):
    if "Func.refGene" in vcf_info:
        function = vcf_info["Func.refGene"]
        if function == 'exonic':
            return True
    return False

def is_intronic(func):
    if func == 'intronic':
        return True
    return False

def is_utr(func):
    if func in ['UTR5', 'UTR3']:
        return True
    return False

def is_intergenic(func):
    if func == 'intergenic':
        return True
    return False

def is_neighboring(func):
    if func in ['upstream', 'downstream']:
        return True
    return False

def is_splicing(func):
    if func == 'splicing':
        return True
    return False

def is_ncrna(func):
    if func in ['ncRNA_exonic', 'ncRNA_intronic']:
        return True
    return False

def is_lof(func):
    if func in ['frameshift_insertion', 'frameshift_deletion', 'stopgain', 'stoploss']:
        return True
    return False

def is_deleterious(dann, func):
    #print(dann, func)
    if is_dann_pathogenic(dann):
        return True
    if is_lof(func):
        return True
    return False

def is_deleterious_2(vcf_info, ref, alt):
    # Check Loftee
    loftee = is_loftee_lof(vcf_info)
    if loftee is True:
        return True
    if loftee is False:
        return False
    # If loftee returns None that means there was no prediction for this variant

    if is_exonic(vcf_info):
        # Check whether it's a frameshift indel.  Loftee should be checking these already.  This is here just in case.
        if is_frameshift(ref, alt):
            return True

        # For SNPs, do majority vote of LRT, MutationAssessor, Provean, VEST, and CADD
        if pgx_coding_deleterious(vcf_info):
            return True

    # If it's non-coding, do majority vote of cadd, dann, and fathmm
    else:
        if pgx_noncoding_deleterious(vcf_info):
            return True

    return False



def pgx_coding_deleterious(vcf_info):
    lrt = is_lrt_deleterious(vcf_info)
    ma = is_mutation_assessor_deleterious(vcf_info)
    cadd = is_cadd_deleterious(vcf_info)
    provean = is_provean_deleterious(vcf_info)
    vest = is_vest_deleterious(vcf_info)
    votes = sum([lrt, ma, cadd, provean, vest])
    if votes >= 3:


        return True
    return False

def pgx_noncoding_deleterious(vcf_info):
    dann = is_dann_deleterious(vcf_info, noncoding=True)
    fathmm = is_fathmm_deleterious(vcf_info, noncoding=True)
    cadd = is_cadd_deleterious(vcf_info, noncoding=True)

    votes = sum([dann, fathmm, cadd])
    if votes >= 2:
        return True
    return False


def is_lrt_deleterious(vcf_info):
    if "LRT_score" in vcf_info:
        LRT_score = vcf_info["LRT_score"]
        if LRT_score == ".":
            return False
        if float(LRT_score) < 0.0025:
            return True
    return False

def is_mutation_assessor_deleterious(vcf_info):
    if "MutationAssessor_score" in vcf_info:
        MutationAssessor_score = vcf_info["MutationAssessor_score"]
        if MutationAssessor_score == ".":
            return False
        if float(MutationAssessor_score) > 2.0566:
            return True
    return False

def is_provean_deleterious(vcf_info):
    if "PROVEAN_score" in vcf_info:
        PROVEAN_score = vcf_info["PROVEAN_score"]
        if PROVEAN_score == ".":
            return False
        if float(PROVEAN_score) < -3.286:
            return True
    return False

def is_vest_deleterious(vcf_info):
    if "VEST3_score" in vcf_info:
        VEST3_score = vcf_info["VEST3_score"]
        if VEST3_score == ".":
            return False
        if float(VEST3_score) > 0.4523:
            return True
    return False

def is_cadd_deleterious(vcf_info, noncoding=False):
    cutoff = 19.19
    if noncoding:
        cutoff = 15

    if "CADD_Phred" in vcf_info:
        CADD_Phred = vcf_info["CADD_Phred"]
        if CADD_Phred == ".":
            return False
        if float(CADD_Phred) > cutoff:
            return True
    return False

def is_dann_deleterious(vcf_info, noncoding=False):
    cutoff = 0.9688
    if noncoding:
        cutoff = 0.9

    if "dann" in vcf_info:
        dann = vcf_info["dann"]
        if dann == ".":
            return False
        if float(dann) > cutoff:
            return True
    return False

def is_fathmm_deleterious(vcf_info, noncoding=False):

    if noncoding:
        cutoff = 0.9
        if "FATHMM_noncoding" in vcf_info:
            CADD_Phred = vcf_info["FATHMM_noncoding"]
            if CADD_Phred == ".":
                return False
            if float(CADD_Phred) > cutoff:
                return True

    else:
        cutoff = 0.3982
        if "FATHMM_coding" in vcf_info:
            CADD_Phred = vcf_info["FATHMM_coding"]
            if CADD_Phred == ".":
                return False
            if float(CADD_Phred) > cutoff:
                return True


    return False

def is_loftee_lof(vcf_info):
    # Get the VEP part
    vep_string = 'Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|REFSEQ_MATCH|SOURCE|GIVEN_REF|USED_REF|BAM_EDIT|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info'
    vep_fields = vep_string.split('|')

    if "CSQ" in vcf_info:
        vep_annotations = [dict(zip(vep_fields, x.split('|'))) for x in vcf_info['CSQ'].split(',')][0]
        #print(vep_annotations)
        if vep_annotations['LoF'] == "HC":
            return True
        if vep_annotations['LoF'] == "LC":
            return False

    return None

def amino_acid_change(vcf_info):
    # Get the VEP part
    vep_string = 'Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|REFSEQ_MATCH|SOURCE|GIVEN_REF|USED_REF|BAM_EDIT|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info'
    vep_fields = vep_string.split('|')

    if "CSQ" in vcf_info:
        vep_annotations = [dict(zip(vep_fields, x.split('|'))) for x in vcf_info['CSQ'].split(',')][0]
        #print(vep_annotations)

        protein_pos = vep_annotations['Protein_position']
        amino_acid = vep_annotations['Amino_acids']

        if "/" in amino_acid:
            aminos = amino_acid.split("/")
            change = "%s%s%s" % (aminos[0], protein_pos, aminos[1])
            return change

        if len(amino_acid) == 1:
            change = "%s%s%s" % (amino_acid, protein_pos, amino_acid)
            return change
        
        return "NA"

def is_methylation_site(vcf_info):
    # Check wgEncodeHaibMethyl450Gm12878SitesRep1 & wgEncodeHaibMethylRrbsGm12878HaibSitesRep1
    # If either are true return True
    if "wgEncodeHaibMethyl450Gm12878SitesRep1" in vcf_info:
        wgEncodeHaibMethyl450Gm12878SitesRep1 = vcf_info["wgEncodeHaibMethyl450Gm12878SitesRep1"]
        if wgEncodeHaibMethyl450Gm12878SitesRep1 != ".":
            return True
    if "wgEncodeHaibMethylRrbsGm12878HaibSitesRep1" in vcf_info:
        wgEncodeHaibMethylRrbsGm12878HaibSitesRep1 = vcf_info["wgEncodeHaibMethylRrbsGm12878HaibSitesRep1"]
        if wgEncodeHaibMethylRrbsGm12878HaibSitesRep1 != ".":
            return True
    return False

def is_tf_binding_site(vcf_info):
    # Check wgEncodeRegTfbsClusteredV3 & tfbsConsSites
    if "tfbsConsSites" in vcf_info:
        tfbsConsSites = vcf_info["tfbsConsSites"]
        if tfbsConsSites != ".":
            return True
    if "wgEncodeRegTfbsClusteredV3" in vcf_info:
        wgEncodeRegTfbsClusteredV3 = vcf_info["wgEncodeRegTfbsClusteredV3"]
        if wgEncodeRegTfbsClusteredV3 != ".":
            return True
    return False

def is_eQTL(vcf_info):
    # Check CYP2D6_eQTLs
    if "eQTL" in vcf_info:
        eQTL = vcf_info["eQTL"]
        if eQTL != ".":
            return True
    return False


def is_dnase_hypersensitivity_site(vcf_info):
    # Check wgEncodeAwgDnaseMasterSites
    if "wgEncodeAwgDnaseMasterSites" in vcf_info:
        wgEncodeAwgDnaseMasterSites = vcf_info["wgEncodeAwgDnaseMasterSites"]
        if wgEncodeAwgDnaseMasterSites != ".":
            return True
    return False


def clean_chr(chr):
    if not chr:
        return chr
    if chr.startswith('chr'):
        return chr[3:]
    return chr

def is_deletion(ref, alt):
    if len(alt) < len(ref):
        return True
    return False
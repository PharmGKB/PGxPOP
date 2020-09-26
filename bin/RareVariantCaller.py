"""
Greg McInes
Altman Lab
gmcinnes@stanford.edu
"""

from DawgToys import parse_vcf_line, get_vcf_subject_ids, get_definition_file, gt_splitter
from Gene import Gene
import argparse
import sys
import os
import tabix

class RareVariantCaller(object):
    def __init__(self, vcf, gene, build='grch38', debug=False):
        self.vcf = vcf
        self.gene = gene
        self.build = build
        self.debug = debug

    def run(self):
        if self.debug:
            print("Extracting rare variants in %s" % self.gene.name)

        # Get the bed regions for the gene
        regions = self.coding_regions()

        if regions is None:
            return None

        # Get the variants from the definition file for this gene
        rare_variants = self.get_rare_variants(regions[self.gene.name])

        return rare_variants

    def get_rare_variants(self, regions):

        # Create dictionary to store results
        sample_variants = {}
        subjects = get_vcf_subject_ids(self.vcf)
        for s in subjects:
             sample_variants[s] = []

        known_vars = self.known_variants()

        # For each region, extract all the variants
        for r in regions:
            variants, keys = self.get_region_variants(r, exclude=known_vars)

            for i in range(len(variants)):
                vcf_row = parse_vcf_line(variants[i])
                for j in range(len(subjects)):
                    s = subjects[j]
                    gt = vcf_row.calls[j].split(':')[0]

                    # If the genotype contains any alternate allele, store it
                    # For best performance SNPs should be split across multiple rows in the VCF
                    # Otherwise it won't be possible to decipher multiallelic sites
                    if '1' in gt or '2' in gt or '3' in gt:

                        sample_variants[s].append(keys[i])

        return sample_variants

    def known_variants(self):
        variants = []

        for v in self.gene.variants:
            chrom = self.gene.variants[v].chromosome
            pos = self.gene.variants[v].position
            ref = self.gene.variants[v].ref

            for a in self.gene.variants[v].alt:
                key = "%s_%s_%s_%s" % (chrom, pos, ref, a)
                variants.append(key)

        return variants

    def coding_regions(self):

        regions = {}

        wd = os.path.dirname(os.path.realpath(__file__))

        if self.build == "grch38":
            bed = os.path.join(wd, "../definition/exome_target_regions.bed")
        elif self.build == 'hg19':
            bed = os.path.join(wd, "../definition/exome_target_regions.hg19.bed")
        else:
            print("Build %s not supported for rare variant discovery.")
            return None

        with open(bed) as f:
            for line in f:
                fields = line.rstrip().split()
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                gene = fields[3]

                if not gene in regions:
                    regions[gene] = []

                r = {"chr": chrom,
                     "start": start,
                     "end": end,
                     "gene": gene}

                regions[gene].append(r)

        return regions

    def get_region_variants(self, region, exclude=None):

        variant_records = []
        keys = []

        chrom = region["chr"]
        start = region["start"]
        end = region["end"]

        tb = tabix.open(self.vcf)

        records = tb.query("%s" % chrom, start, end)

        n_records = 0
        for r in records:
            n_records += 1

            c = r[0]
            p = r[1]
            ref = r[3]
            alt = r[4]

            if alt != ".":

                if exclude is not None:
                    alts = alt.split(",")
                    for a in alts:
                        main_key = "%s_%s_%s_%s" % (c, p, ref, a)

                        # SNP
                        if len(ref) == len(a):

                            key = "%s_%s_%s_%s" % (c, p, ref, a)
                            if key not in exclude:
                                variant_records.append(r)
                                keys.append(main_key)

                        # Deletion
                        elif len(ref) > len(a):

                            deletion = ref[len(a):]


                            key = "%s_%s_%s_del%s" % (c, p, deletion, deletion)

                            if key not in exclude:
                                variant_records.append(r)
                                keys.append(main_key)



                        # Insertion
                        elif len(ref) < len(a):
                            insertion = a[len(ref):]
                            key = "%s_%s_ins_%s" % (c, p, insertion)

                            if key not in exclude:
                                variant_records.append(r)
                                keys.append(main_key)

                else:
                    alts = alt.split(",")
                    for a in alts:
                        main_key = "%s_%s_%s_%s" % (c, p, ref, a)
                        variant_records.append(r)
                        keys.append(main_key)



        return variant_records, keys

"""
Parse the command line
"""
def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'This is a script I wrote')
    parser.add_argument("-f", "--file", help="Input")
    parser.add_argument("-g", "--gene", help="Input")

    parser.add_argument("-d", "--debug", action='store_true', default=False,
                                help="Output debugging messages.  May be very verbose.")
    options = parser.parse_args()
    return options


"""
Main
"""
if __name__ == "__main__":
    options = parse_command_line()
    gene = Gene(get_definition_file(options.gene))
    caller = RareVariantCaller(vcf=options.file, gene=gene, debug=options.debug)
    variants = caller.run()
    print(variants)


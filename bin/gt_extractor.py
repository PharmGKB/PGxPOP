"""
Greg McInes
Altman Lab
gmcinnes@stanford.edu
"""


import sys
import argparse
from DawgToys import fetch_genotypes_2, fetch_genotype_records, get_vcf_subject_ids, parse_vcf_line, split_genotype

class GTExtractor(object):
    def __init__(self, file, chromosome, position, subject=None, carriers_only=False, homs_only=False, debug=False,
                 list=None):
        self.debug = debug
        self.file = file
        self.chr = chromosome
        self.position = position
        self.subject = subject
        self.carriers_only = carriers_only
        self.homs_only = homs_only
        self.list = list

        self.run()

    def run(self):
        # Get the subjects
        subjects = get_vcf_subject_ids(self.file)

        if self.subject is not None:
            if not self.subject in subjects:
                print("Subject not found in VCF")
                exit(1)

        if self.list is not None:
            # iterate over the list, which should look just like a VCF
            with open(self.list) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    vcf = parse_vcf_line(line)
                    chr = vcf.chrom
                    pos = vcf.pos
                    ref = vcf.ref
                    alt = vcf.alt

                    # if you want to pull something out of the INFO to include
                    extras = {
                        "AF": vcf.info['AF'],
                        #"Gene": vcf.info['Gene.refGene'],
                        "Group": "deleterious",
                        "impact": vcf.info["impact"],
                        "high_impact": vcf.info["highimpact"],
                        "Gene": vcf.info["gene"]

                    }

                    self.process_variant(chr, pos, subjects, vcf_row=vcf, extras=extras)
        else:
            self.process_variant(self.chr, self.position, subjects)


    def process_variant(self, chr, position, subjects, ref=None, alt=None, extras=None, vcf_row=None):
        # Extract the variant data

        if vcf_row is not None:
            #vcf_row.print_row(minimal=True)
            variant = fetch_genotypes_2(self.file, vcf_row)
        else:
            variant = fetch_genotype_records(self.file, chr, position)

        if variant is None or variant.alt == ".":
            #print("no variant found")
            return



        #variant.print_row(minimal=True)

        if vcf_row is not None:
            if variant.ref != vcf_row.ref:
                print("REF DIDN'T MATCH")
                print(vcf_row.ref, variant.ref)
                return
            if vcf_row.alt not in variant.alts():
                print("ALTS DIDN'T MATCH")
                print(vcf_row.alt, variant.alts())
                #ariant.print_row()
                #exit()
                return

        if variant is None:
            #print("Position not found in VCF: %s:%s" % (chr, position), file=sys.stderr)
            return
            #exit(1)


        alt_index = 1
        if vcf_row is not None:
            if not variant.ref == vcf_row.ref:
                #print(f"Wrong reference! {ref}", file=sys.stderr)
                return
            if not vcf_row.alt in variant.alts():
                #print("Wrong alt!", file=sys.stderr)
                #print(vcf_row.alt, variant.alts())
                return
            alts = variant.alts()
            #print(alts)
            if len(alts) > 1:
                for i in range(len(alts)):
                    if alts[i] == vcf_row.alt:
                        #print("Found alt!")
                        #print(alts[i], vcf_row.alt)
                        alt_index = i + 1


        #print(alt_index)

        # Get the resultw

        # If a subject is provided, just print that one
        if self.subject is not None:
            index = subjects.index(self.subject)
            print("%s: %s" % (self.subject, variant.calls[index]))
            return

        # If no subject is given, print all the subjects
        for i in range(len(subjects)):
            if self.carriers_only is True:
                gt = split_genotype(variant.calls[i])
                if not str(alt_index) in gt:
                    continue

            elif self.homs_only is True:
                gt = split_genotype(variant.calls[i])
                if not gt[0] in ['1', '2'] and not gt[1] in ['1', '2']:
                    continue
            #else:
            #    print("%s: %s" % (subjects[i], variant.calls[i]))

            gt = variant.calls[i].split(":")[0]
            ac = sum(int(x) for x in split_genotype(variant.calls[i])) / alt_index


            if extras:
                values = []
                for k in extras.keys():
                    values.append(extras[k])

                print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (subjects[i], chr, position, vcf_row.ref, vcf_row.alt, gt, ac, "\t".join(values)))
            else:
                print("%s\t%s\t%s\t%s\t%s" % (subjects[i], chr, position, gt, ac))

"""
Parse the command line
"""
def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'Extract and print genotypes for subjects at a specific position in a VCF')
    parser.add_argument("-f", "--vcf", help="Input")
    parser.add_argument("-c", "--chr", help="Chromosome")
    parser.add_argument("-p", "--position", type=int, help="Position")
    parser.add_argument('-l', "--list", help="List of variants to query")
    parser.add_argument("-s", "--subject", help="Subject ID")
    parser.add_argument("--carriers", action='store_true', help="Return only carriers of an alternate allele")
    parser.add_argument("--homo", action='store_true', help="Return only subjects homozygous for the alternate allele")
    parser.add_argument("-d", "--debug", action='store_true', default=False,
                                help="Output debugging messages.  May be very verbose.")
    options = parser.parse_args()
    return options


"""
Main
"""
if __name__ == "__main__":
    options = parse_command_line()
    GTExtractor(file=options.vcf,
                chromosome=options.chr,
                position=options.position,
                subject=options.subject,
                carriers_only=options.carriers,
                homs_only=options.homo,
                debug=options.debug,
                list=options.list)


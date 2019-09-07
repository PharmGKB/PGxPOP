"""
Greg McInes
Altman Lab
gmcinnes@stanford.edu
"""


import argparse
from DawgToys import fetch_genotype_records, get_vcf_subject_ids, parse_vcf_line, split_genotype

class GTExtractor(object):
    def __init__(self, file, chromosome, position, subject=None, carriers_only=False, homs_only=False, debug=False):
        self.debug = debug
        self.file = file
        self.chr = chromosome
        self.position = int(position)
        self.subject = subject
        self.carriers_only = carriers_only
        self.homs_only = homs_only

        self.run()

    def run(self):
        # Get the subjects
        subjects = get_vcf_subject_ids(self.file)

        if self.subject is not None:
            if not self.subject in subjects:
                print("Subject not found in VCF")
                exit(1)

        # Extract the variant data
        variant = fetch_genotype_records(self.file, self.chr, self.position)

        if variant is None:
            print("Position not found in VCF")
            exit(1)

        # Parse the result
        #vcf_line = parse_vcf_line(variant)

        variant.print_row(minimal=True)

        # Get the result

        # If a subject is provided, just print that one
        if self.subject is not None:
            index = subjects.index(self.subject)
            print("%s: %s" % (self.subject, variant.calls[index]))
            return

        # If no subject is given, print all the subjects
        for i in range(len(subjects)):
            if self.carriers_only is True:
                gt = variant.calls[i]
                if '1' in gt:
                    print("%s: %s" % (subjects[i], variant.calls[i]))
            elif self.homs_only is True:
                gt = variant.calls[i]
                if gt[0] in ['1', '2'] and gt[1] in ['1', '2']:
                    print("%s: %s" % (subjects[i], variant.calls[i]))
            else:
                print("%s: %s" % (subjects[i], variant.calls[i]))

"""
Parse the command line
"""
def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'Extract and print genotypes for subjects at a specific position in a VCF')
    parser.add_argument("-f", "--vcf", help="Input")
    parser.add_argument("-c", "--chr", help="Chromosome")
    parser.add_argument("-p", "--position", help="Position")
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
                debug=options.debug)


"""
Greg McInes
Altman Lab
gmcinnes@stanford.edu
"""

# Check all definitions RSIDs against a VCF file for any positions where
# the RSID is the same but the position is different

import argparse
import sys
import os

from Gene import Gene
from DawgToys import parse_vcf_line

class SynonymChecker(object):
    def __init__(self, file, build='grch38', debug=False):
        self.file = file
        self.build = build
        self.debug = debug

        self.run()

    def run(self):

        # Get all the positions associated with star allele definitions
        definition_rsids = self.get_rsids()

        # Iterate over the VCF and check whether there are any rsids where the position is different in the
        # VCF than it is in the definition
        with open(self.file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                vcf_line = parse_vcf_line(line)
                if vcf_line.id in definition_rsids.keys():
                    print("rsid found: %s" % vcf_line.id)
                    variant = definition_rsids[vcf_line.id]
                    if str(variant["position"]) != str(vcf_line.pos):
                        print("Position mismatch!")
                        vcf_line.print_row()
                        print(definition_rsids[vcf_line.id])
                    # Also check for allele flips
                    if variant["ref"] != vcf_line.ref:
                        print("Ref mismatch!")
                        vcf_line.print_row()
                        print(definition_rsids[vcf_line.id])


    def get_rsids(self):
        variants = {}

        genes = ['CFTR', 'CYP2C9', 'CYP4F2', 'IFNL3', 'TPMT', 'VKORC1', 'NUDT15',
                 'CYP2C19', 'CYP3A5',  'DPYD', 'SLCO1B1', 'UGT1A1', 'CYP2D6', 'CYP2B6', ]

        for g in genes:
            gene_definition = self.get_definition_file(g)
            gene = Gene(gene_definition, build=self.build, debug=self.debug)

            for v in gene.variants:
                if gene.variants[v].rsid != None:
                    variants[gene.variants[v].rsid] = {
                        "position": gene.variants[v].position,
                        "ref": gene.variants[v].ref,
                        "alt": gene.variants[v].alt
                    }

        return variants


    def get_definition_file(self, g):
        definition_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../definition/alleles/")
        filename = "%s_translation.json" % g
        definition_file = os.path.join(definition_dir, filename)
        return definition_file

"""
Parse the command line
"""
def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'Check all definitions RSIDs against a VCF file for any positions '
                      'where the RSID is the same but the position is different')
    parser.add_argument("-f", "--file", help="Input")
    parser.add_argument("-b", "--build", default="grch38", help="Alternate genome build")
    parser.add_argument("-d", "--debug", action='store_true', default=False,
                                help="Output debugging messages.  May be very verbose.")
    options = parser.parse_args()
    return options


"""
Main
"""
if __name__ == "__main__":
    options = parse_command_line()
    SynonymChecker(options.file, options.build, options.debug)


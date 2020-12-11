"""
Greg McInes
Altman Lab
gmcinnes@stanford.edu
"""


import argparse
import os

class SumStats(object):
    def __init__(self, file, prefix, debug=False):
        self.file = file
        self.prefix = prefix
        self.debug = debug

        self.gene_counts = {}
        self.haplotypes = {}
        self.diplotypes = {}
        self.phenotypes = {}
        self.phenotypes_presumptive = {}

        self.run()

    def run(self):
        with open(self.file) as f:
            for line in f:
                fields = line.rstrip().split(',')
                if fields[0] == "sample_id":
                    continue

                gene = fields[1]
                diplotype_raw = fields[2]
                hap_1 = fields[3]
                hap_2 = fields[4]
                phenotype = fields[9]
                phenotype_presumptive = fields[12]

                if "|" in diplotype_raw:
                    diplotype_seperator = "|"
                else:
                    diplotype_seperator = "/"

                diplotype = diplotype_seperator.join(sorted([hap_1, hap_2]))

                if not gene in self.gene_counts:
                    self.gene_counts[gene] = 0
                    self.haplotypes[gene] = {}
                    self.diplotypes[gene] = {}
                    self.phenotypes[gene] = {}
                    self.phenotypes_presumptive[gene] = {}

                self.gene_counts[gene] += 1
                self.update_dict(self.haplotypes[gene], hap_1)
                self.update_dict(self.haplotypes[gene], hap_2)
                self.update_dict(self.diplotypes[gene], diplotype)
                self.update_dict(self.phenotypes[gene], phenotype)
                self.update_dict(self.phenotypes_presumptive[gene], phenotype_presumptive)

        self.print_results()

    def update_dict(self, dict, item):
        if not item in dict:
            dict[item] = 0
        dict[item] += 1

    def print_results(self):
        hap_file = open("%s.haplotypes.tsv" % self.prefix, "w")
        hap_file.write("gene\thaplotype\tcount\n")
        for g in self.haplotypes:
            for h in self.haplotypes[g]:
                hap_file.write("%s\t%s\t%s\n" % (g, h, self.haplotypes[g][h]))
        hap_file.close()

        dip_file = open("%s.diplotypes.tsv" % self.prefix, "w")
        dip_file.write("gene\tdiplotype\tcount\n")
        for g in self.haplotypes:
            for d in self.diplotypes[g]:
                dip_file.write("%s\t%s\t%s\n" % (g, d, self.diplotypes[g][d]))
        dip_file.close()

        phen_file = open("%s.phenotypes.tsv" % self.prefix, "w")
        phen_file.write("gene\tphenotype\tcount\n")
        for g in self.phenotypes_presumptive:
            for p in self.phenotypes_presumptive[g]:
                phen_file.write("%s\t%s\t%s\n" % (g, p, self.phenotypes_presumptive[g][p]))
        phen_file.close()

        #hap_file = open("%s.haplotypes.tsv" % self.prefix, "w")
        #hap_file.write("gene\thaplotype\tcount\n")
        #for g in self.haplotypes:
        #    for h in self.haplotypes[g]:
        #        hap_file.write("%s\t%s\t%s\n" % (g, h, self.haplotypes[g][h]))
        #hap_file.close()

"""
Parse the command line
"""
def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'Calculate summary statistics from PGxPOP output.  Prints to stdout')
    parser.add_argument("-f", "--file", help="Output from PGxPOP run")
    parser.add_argument("-p", "--prefix", help="Output prefix")
    parser.add_argument("-d", "--debug", action='store_true', default=False,
                                help="Output debugging messages.  May be very verbose.")
    options = parser.parse_args()
    return options


"""
Main
"""
if __name__ == "__main__":
    options = parse_command_line()
    SumStats(options.file, options.prefix, options.debug)


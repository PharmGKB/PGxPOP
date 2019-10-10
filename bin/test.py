"""
Greg McInes
Altman Lab
gmcinnes@stanford.edu
"""


import argparse
import os

from CityDawg import CityDawg
from DiplotypeCaller import DiplotypeCaller
from DawgToys import *
from Gene import Gene

class DawgTest(object):
    def __init__(self, file, debug=False):
        self.file = file
        self.debug = debug

        self.run()

    def run(self):

        # Two parts here
        # 1. Check test files defined by PharmCAT
        #    This will test the VCF loader and caller
        # 2. Test all possible combinations of named alleles
        #    This will check for any systematic bias when calling haplotypes from phased and unphased data.


        self.pharmcat_tests()

        self.combo_tests()


    def combo_tests(self):
        print("Running haplotype combination tests")

        genes = ['CFTR', 'CYP2C9', 'CYP2D6', 'CYP4F2', 'IFNL3', 'TPMT', 'VKORC1',
                 'CYP2C19', 'CYP3A5', 'DPYD', 'SLCO1B1', 'UGT1A1']

        for g in genes:
            # Create the gene object
            gene_definition = get_definition_file(g)
            gene = Gene(gene_definition, debug=self.debug)
            for phased in [False, True]:
                dc = DiplotypeCaller(gene, is_phased=phased)
                combo_gene_results = dc.test_gene(gene)
                print("%s error rate (phased = %s): %s" % (g, phased, combo_gene_results['error rate']))
            #print(combo_gene_results.keys())





    def pharmcat_tests(self, gene='all'):
        print("Executing PharmCAT tests")
        tests = self.get_tests()

        failed = []
        missed_alleles = {}
        for t in tests:
            if t['skip'] == '1':
                continue
            result, missed = self.run_pharmcat_test(t)

            if result is not None:
                t['result'] = result
                failed.append(t)

                if t['gene'] not in missed_alleles.keys():
                    missed_alleles[t['gene']] = {}
                for h in missed:
                    if h not in missed_alleles[t['gene']].keys():
                        missed_alleles[t['gene']][h] = 0
                    missed_alleles[t['gene']][h] += 1

        print("--------------------------------------------------------------------")
        print("%s failed tests" % len(failed))
        for f in failed:
            print("%s %s %s" % (f['gene'], f['file'], f['result']))

        print(missed_alleles)


    def run_pharmcat_test(self, test):
        #if self.debug:
        print("--------------------------------------------------------------------")
        print("Gene: %s, File: %s" % (test['gene'], test['file']))

        test_file = self.test_path(test['gene'], test['file'])


        if test['phased'] == '1':
            phased = True
        else:
            phased = False

        failure = None
        missed = []

        try:
            cd = CityDawg(vcf=test_file, gene=test['gene'], phased=phased, debug=self.debug)
            result = cd.process_gene(test['gene'])[0]
        except:
            failure = "parse error"
            print("TEST FAILED")
            print("Parse error")
            return failure, missed


        test_haps = [test['hap1'], test['hap2']]
        result_haps = [result['hap_1'], result['hap_2']]

        test_haps.sort()
        result_haps.sort()

        if test_haps == result_haps:
            print("Test passed")

        else:
            print("TEST FAILED")
            print("Expected: %s" % ",".join(test_haps))
            print("Found: %s" % ",".join(result_haps))
            failure = "Match error"

            missed = []
            for h in test_haps:
                if h not in result_haps:
                    missed.append(h)

        return failure, missed


    def test_path(self, gene, file):
        test_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../test/%s/" % gene)
        test_file = os.path.join(test_dir, file)
        return test_file

    def get_tests(self):
        tests = []
        with open(self.get_pharmcat_test_file()) as f:
            for line in f:
                if line.startswith("gene"):
                    continue
                fields = line.rstrip().split(",")
                new_row = {
                    'gene':fields[0],
                    'file': fields[1],
                    'hap1': fields[2],
                    'hap2': fields[3],
                    'phased': fields[4],
                    'skip': fields[5]
                }
                tests.append(new_row)
        return tests

    def get_pharmcat_test_file(self):
        definition_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../test/")
        filename = "pharmcat_test_cases.csv"
        definition_file = os.path.join(definition_dir, filename)
        return definition_file


"""
Parse the command line
"""
def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'This is a script I wrote')
    parser.add_argument("-f", "--file", help="Input")
    parser.add_argument("-d", "--debug", action='store_true', default=False,
                                help="Output debugging messages.  May be very verbose.")
    options = parser.parse_args()
    return options


"""
Main
"""
if __name__ == "__main__":
    options = parse_command_line()
    DawgTest(options.file, options.debug)


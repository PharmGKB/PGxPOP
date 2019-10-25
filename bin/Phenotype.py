import json
import os
import numpy as np

'''
Define the Gene object

This class will be initiated with a Gene definition json and contain functions related to gene specific tasks.
'''
class Phenotype(object):
    def __init__(self, json_file=None, debug=False):

        # Set global variables
        self.data = None
        self.debug = debug

        self._load_json(json_file)

    def _load_json(self, json_file):
        if json_file is None:
            json_file = self.default_json_file()
        with open(json_file) as f:
            self.data = json.load(f)
        f.close()

    def default_json_file(self):
        definition_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../definition/")
        filename = "gene.phenotypes.json"
        phenotype_file = os.path.join(definition_dir, filename)
        return phenotype_file

    def get_gene(self, gene_name):
        for g in self.data:
            if g["gene"] == gene_name:
                return g
        #print("Gene not found in phenotypes! %s" % gene_name)
        return None

    def get_haplotype_function(self, gene_name, hap):
        # Determine the function of a single haplotype
        gene = self.get_gene(gene_name)
        if gene is None:
            return None
        if hap in gene["haplotypes"].keys():
            return gene["haplotypes"][hap]
        #print("Haplotype %s not found in %s phenotypes!" % (hap, gene_name))
        return None
        #exit(1)

    def get_diplotype_function(self, gene_name, hap1, hap2):

        # Haplotypes with equal scores will now be output as H1A;H1B, so we'll split those
        # and check the function of each separately.  If they have the same function, we
        # can assign the function.  If it is different, then we will report "Conflicting"


        # Determine the function of two haplotypes for a given gene
        hap1_function = self.haplotype_checker(gene_name, hap1)
        hap2_function = self.haplotype_checker(gene_name, hap2)



        if hap1_function is None or hap2_function is None:
            return None

        if hap1_function == "Conflicting" or hap2_function == "Conflicting":
            return "Conflicting"

        hap_functions = sorted([hap1_function, hap2_function])

        # Next we need to iterate over the list of diplotypes and do a list comparison
        gene = self.get_gene(gene_name)
        for d in gene['diplotypes']:
            if sorted(d['diplotype']) == hap_functions:
                return d['phenotype']
        print('Unable to determine function of %s and %s for %s' % (hap1, hap2, gene_name))
        print('%s: %s' % (hap1, hap1_function))
        print('%s: %s' % (hap2, hap2_function))
        exit(1)

    def haplotype_checker(self, gene_name, hap):
        found_functions = set()
        for h in hap.split(";"):
            found_functions.add(self.get_haplotype_function(gene_name, h))
        if len(found_functions) > 1:
            return "Conflicting"
        return list(found_functions)[0]


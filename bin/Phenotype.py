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
        print("Gene not found in phenotypes! %s" % gene_name)
        exit(1)

    def get_haplotype_function(self, gene_name, hap):
        # Determine the function of a single haplotype
        gene = self.get_gene(gene_name)
        if hap in gene["haplotypes"].keys():
            return gene["haplotypes"][hap]
        print("Haplotype %s not found in %s phenotypes!" % (hap, gene_name))
        exit(1)

    def get_diplotype_function(self, gene_name, hap1, hap2):
        # Determine the function of two haplotypes for a given gene
        hap1_function = self.get_haplotype_function(gene_name, hap1)
        hap2_function = self.get_haplotype_function(gene_name, hap2)
        hap_functions = [hap1_function, hap2_function].sort()

        # Next we need to iterate over the list of diplotypes and do a list comparison
        gene = self.get_gene(gene_name)
        for d in gene['diplotypes']:
            if d['diplotype'].sort() == hap_functions:
                return d['phenotype']
        print('Unable to determine function of %s and %s for %s' % (hap1, hap2, gene_name))
        print('%s: %s' % (hap1, hap1_function))
        print('%s: %s' % (hap2, hap2_function))
        exit(1)

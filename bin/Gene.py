import json
from Variant import Variant
from NamedAllele import NamedAllele
import numpy as np

'''
Define the Gene object

This class will be initiated with a Gene definition json and contain functions related to gene specific tasks.
'''
class Gene(object):
    def __init__(self, json_file, build='grch38', debug=False):

        # Set global variables
        self.data = None
        self.build = build
        self.debug = debug

        # Load gene information from JSON file
        self._load_json(json_file)

        self.variants = self.get_variants()
        self.haplotypes = self.get_haplotypes()



        # functions needed
        # haplotype definition matrix
        #   This will return a binary matrix of all the alleles for each haplotype
        #
        # Return a list of variants and positions

    def _load_json(self, json_file):
        with open(json_file) as f:
            self.data = json.load(f)
        f.close()

    def haplotype_matrix(self):
        # contstruct a binary matrix with the alleles along the x axis and haplotypes on the y axis
        # 0's indicate a ref allele and 1's an alternate

        all_hap_alleles = []
        star_index = []
        for h in self.haplotypes:
            star_index.append(self.haplotypes[h].name)
            all_hap_alleles.append(self.haplotypes[h].binary_alleles)

        hap_matrix = np.array(all_hap_alleles)
        return(hap_matrix, star_index)


    def get_variants(self):
        if self.debug:
            print("Formatting variants")
            if self.build != "grch38":
                print("Alternate build detected.  Coordinate conversion may take a few minutes.")
        variants = self.data['variants']
        alleles = self.data['variantAlleles']
        formatted_variants = {}

        index = 0
        for v in variants:
            new_variant = Variant(v, alleles=alleles[index], index=index, build=self.build, debug=self.debug)
            if self.debug:
                new_variant.print_variant()
            formatted_variants[new_variant.index] = new_variant
            index += 1

        return formatted_variants

    def get_haplotypes(self):
        if self.debug:
            print("Formatting haplotypes")
        haplotypes = self.data['namedAlleles']

        formatted_haplotypes = {}

        index = 0
        for h in haplotypes:
            new_haplotype = NamedAllele(h, variants=self.variants, index=index, debug=self.debug)
            if self.debug:
                new_haplotype.print_haplotype()
            formatted_haplotypes[new_haplotype.id] = new_haplotype
            index += 1

        return formatted_haplotypes


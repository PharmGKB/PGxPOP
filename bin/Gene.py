'''
Greg McInes
Altman Lab
gmcinnes@stanford.edu
'''

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

        self.name = self.data['gene']
        self.variants = self.get_variants()
        self.haplotypes = self.get_haplotypes()

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
                print("Alternate build detected.")
        variants = self.data['variants']
        alleles = self.data['variantAlleles']
        formatted_variants = {}

        index = 0
        for v in variants:
            new_variant = Variant(v, alleles=[], index=index, build=self.build, debug=self.debug)
            if self.debug:
                new_variant.print_variant()
            formatted_variants[new_variant.index] = new_variant
            index += 1

        return formatted_variants

    def get_haplotypes(self):
        if self.debug:
            print("Formatting haplotypes")
        haplotypes = self.data['namedAlleles']
        stars = {}

        formatted_haplotypes = {}

        index = 0
        for h in haplotypes:
            new_haplotype = NamedAllele(h, variants=self.variants, index=index, debug=self.debug)
            star = new_haplotype.name

            if star not in stars.keys():
                stars[star] = 0
            else:
                substar = f'{star}%{stars[star]}'
                stars[star] += 1
                new_haplotype.name = substar

            if self.debug:
                new_haplotype.print_haplotype()
            formatted_haplotypes[new_haplotype.id] = new_haplotype
            index += 1

        return formatted_haplotypes

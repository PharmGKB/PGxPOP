import json
import os
import sys

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
            return "Not available"
        if hap in gene["haplotypes"].keys():
            return gene["haplotypes"][hap]
        #print("Haplotype %s not found in %s phenotypes!" % (hap, gene_name))
        return "Not available"
        #exit(1)

    def get_presumptive_haplotype_function(self, gene_name, hap):

        # Determine the function of a combination of haplotypes (+ calls)
        gene = self.get_gene(gene_name)
        if gene is None:
            return None

        #if gene == "DPYD" and hap == "c.1905+1G>A":
        #    return self.get_haplotype_function(gene_name, hap)

        # todo list other known +'s

        split_haps = hap.split("+")

        #print(split_haps)
        #exit()

        if len(split_haps) == 1:
            return self.get_haplotype_function(gene_name, hap)

        if gene_name == "CFTR":
            return self.CFTR_presumptive(split_haps)

        minimum_function = "Not available"

        #print(split_haps)

        min_functions = ["No Function", "Loss of function"]
        if gene_name in ["UGT1A1", "SLCO1B1"]:
            min_functions.append("Decreased Function")
            #min_functions.append("Possible Decreased Function")
            print("adding decreased")

        for h in split_haps:
            f = self.get_haplotype_function(gene_name, h)

            # This could be adjusted to account for a ranking.  Right now only checking for no function
            if f in min_functions:
                minimum_function = f

        return minimum_function

    def CFTR_presumptive(self, split_haps):
        gene = self.get_gene("CFTR")

        all_responsive = True

        for h in split_haps:
            f = self.get_haplotype_function("CFTR", h)
            #print(h, f)
            if f != "ivacaftor responsive":
                all_responsive = False

        if all_responsive is True:
            return "ivacaftor responsive"
        else:
            return "Not available"




    def function_rank(self):
        functions = {
            "No function": 0,
            "Loss of function": 0,
            "Decreased function": 1,
            "Probable Decreased Function": 2,
            "Possible Decreased function": 2,
            "Normal function": 3,
            "Increased function": 4

        }

        return functions

    def get_diplotype_function(self, gene_name, hap1, hap2, presumptive=False):

        # Haplotypes with equal scores will now be output as H1A;H1B, so we'll split those
        # and check the function of each separately.  If they have the same function, we
        # can assign the function.  If it is different, then we will report "Conflicting"


        # Determine the function of two haplotypes for a given gene
        hap1_function = self.haplotype_checker(gene_name, hap1, presumptive)
        hap2_function = self.haplotype_checker(gene_name, hap2, presumptive)


        if hap1_function is None or hap2_function is None:
            return None

        if hap1_function == "Not available" or hap2_function == "Not available":
            return "Not available"


        if hap1_function == "Conflicting" or hap2_function == "Conflicting":
            return "Conflicting"

        hap_functions = sorted([hap1_function, hap2_function])

        # Next we need to iterate over the list of diplotypes and do a list comparison
        gene = self.get_gene(gene_name)
        for d in gene['diplotypes']:
            if sorted(d['diplotype']) == hap_functions:
                return d['phenotype']
        print('Unable to determine function of %s and %s for %s' % (hap1, hap2, gene_name), file=sys.stderr)
        #print('%s: %s' % (hap1, hap1_function))
        #print('%s: %s' % (hap2, hap2_function))
        return "Not available"
        #exit(1)


    def haplotype_checker(self, gene_name, hap, presumptive=False):
        found_functions = set()
        for h in hap.split(";"):
            if presumptive is True:
                found_functions.add(self.get_presumptive_haplotype_function(gene_name, h))
            else:
                found_functions.add(self.get_haplotype_function(gene_name, h))
        if len(found_functions) > 1:
            return "Conflicting"
        return list(found_functions)[0]

    def get_activity_score(self, gene_name, hap1, hap2, presumptive=False):

        hap1_as = self.get_hap_activity_score(gene_name, hap1, presumptive)
        hap2_as = self.get_hap_activity_score(gene_name, hap2, presumptive)

        if None in [hap1_as, hap2_as]:
            return None

        return hap1_as + hap2_as

    def get_hap_activity_score(self, gene_name, hap, presumptive=False):
        gene = self.get_gene(gene_name)

        if gene is None:
            return None

        if "haplotype_as" not in list(gene.keys()):
            return None

        if not presumptive:
            if hap in list(gene['haplotype_as'].keys()):
                return gene['haplotype_as'][hap]
            else:
                return None

        else:
            if hap in list(gene['haplotype_as'].keys()):
                return gene['haplotype_as'][hap]
            if self.get_presumptive_haplotype_function(gene_name, hap) == "No Function":
                return 0







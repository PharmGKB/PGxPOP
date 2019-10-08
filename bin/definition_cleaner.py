"""
Greg McInes
Altman Lab
gmcinnes@stanford.edu
"""

# This script will make any necessary changes to the allele definition files.

# Convert all INDELs to null in the reference sequence
# Split apart the TA repeat in UGT1A1
# Add known synonyms

import json
import argparse
from Gene import Gene

# todo put this in a file
synonyms = {
    "CYP2C9": {94949282: [94949281]}
}

class DefinitionCleaner(object):
    def __init__(self, file, debug=False):
        self.file = file
        self.debug = debug

        # Load the data into a global object
        self.data = self.get_data()
        self.gene = Gene(self.file, debug=self.debug)

        self.run()

    def run(self):

        # Handle special cases
        self.special_cases()

        # Add synonynms
        self.add_synonyms()

        # Fix all the INDELs
        self.update_indels()

        # print the output
        print(json.dumps(self.data, indent=2))

    def special_cases(self):

        # todo
        # This should be an expansion of the INDEL in UGT1A1
        # Leaving blank for now

        pass

    def add_synonyms(self):
        current_gene = self.data["gene"]
        if not current_gene in synonyms.keys():
            return
        gene_syns = synonyms[current_gene]
        for pos in gene_syns.keys():
            for i in self.gene.variants:
                if self.gene.variants[i].position == pos:
                    self.data["variants"][i]["synonyms"] = gene_syns[pos]

    def update_indels(self):
        # Get the indices of all INDELs in the file
        for i in self.gene.variants:
            if self.gene.variants[i].type in ["INS", "DEL"]:
                #self.gene.variants[i].print_variant()
                # If an INDEL is identified update the index for the reference named allele
                j = 0
                for n in self.data["namedAlleles"]:
                    if n['name'] in ["*1", "Reference"]:
                        # Set the index of the INDEL to None
                        self.data["namedAlleles"][j]['alleles'][i] = None
                    j += 1

    def get_data(self):
        with open(self.file) as json_file:
            data = json.load(json_file)
        return data





"""
Parse the command line
"""
def parse_command_line():
    parser = argparse.ArgumentParser(
        description = 'This script will make any necessary changes to the allele definition files.')
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
    DefinitionCleaner(options.file, options.debug)


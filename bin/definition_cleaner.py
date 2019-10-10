"""
Greg McInes
Altman Lab
gmcinnes@stanford.edu
"""

# This script will make any necessary changes to the allele definition files.

# Convert all INDELs to null in the reference sequence
# Split apart the TA repeat in UGT1A1
# Add known synonyms
# Add hg19 position

'''
CYP2D6 almost needs it's own script
Due to the whole gene deletion (*5), every single variant is of type deletion in the file, which doesn't make sense
and breaks things on our end.  These need to be changed to the actual type of variant they represent.

Additionally, for *4 IUPAC nucleotide representations are used, perhaps overused. For instance an R is used for a 
purine.  This is not handled at the moment, so I manually changed all them to the nucleotide in the PharmVar definition.
'''

import json
import argparse
from pyliftover import LiftOver
from Gene import Gene

# todo put this in a file
synonyms = {
    "CYP2C9": {94949282: [94949281]},
    "CYP2D6": {42129084: [42129083],
               42128249: [42128248],
               42128242: [42128241],
               42128218: [42128211],
               42128201: [42128200],
               42128174: [42128173],
               42127963: [42127962],
               42127846: [42127845]}
}

# CFTR hg19 synonym
# 117199645: 117199644

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

        # Add hg19 position
        self.add_alternate_build()

        # print the output
        print(json.dumps(self.data, indent=2))

    def special_cases(self):

        # todo
        # Expand the INDEL in UGT1A1 into 3 different variants

        # For CYP3A5, *3 needs to have the R and Y changed to G and T

        #if self.data["gene"] == "CYP2D6":
        #    self.CYP2D6()

        pass

    def CYP2D6(self):
        if self.debug:
            print("Fixing CYP2D6")
        # Fix all the types
        self.fix_type()
        exit()
        pass


    def fix_type(self):
        for i in self.gene.variants:
            print(self.data["variants"][i]["type"])
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

    def add_alternate_build(self):
        # Don't forget to update the synonyms too
        for i in self.gene.variants:
            # lift the chromosome and position
            chr = self.gene.variants[i].chromosome
            pos = self.gene.variants[i].position
            new_chr, new_pos = self.liftover(chr, pos)

            new_synonyms = []
            if "synonyms" in self.data["variants"][i].keys():
                for s in self.data["variants"][i]["synonyms"]:
                    new_chr, s_pos = self.liftover(chr, pos)
                    new_synonyms.append(s_pos)

            alt_build_info = {
                'build':'hg19',
                'chromosome': new_chr,
                'position': new_pos,
                'synonyms': new_synonyms
            }

            self.data["variants"][i]["hg19"] = alt_build_info

    def liftover(self, chromosome, position, build='hg19'):

        # todo
        # Not sure what the failure mode of this tool is.  Will probably need to write a try catch eventually
        # Changing the chromosome and position messes up the key as well.  Could probably fix that.  But i don't have
        # the ref and alt alleles on hand and I don't want to parse them out of chromosomeHgvsName.

        lo = LiftOver('hg38', build)
        lifted = lo.convert_coordinate(chromosome, position)

        new_chromosome = lifted[0][0]
        new_position = lifted[0][1]

        if self.debug:
            print("%s %s -> %s %s" % (chromosome, position, new_chromosome, new_position))

        return new_chromosome, new_position




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


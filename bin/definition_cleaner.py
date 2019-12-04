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
***FIXED***
CYP2D6 almost needs it's own script
Due to the whole gene deletion (*5), every single variant is of type deletion in the file, which doesn't make sense
and breaks things on our end.  These need to be changed to the actual type of variant they represent.
*6 has an extra variant in it S486T
*3 extra: N166D
***********

***FIXED***
Additionally, for *4 IUPAC nucleotide representations are used, perhaps overused. For instance an R is used for a 
purine.  This is not handled at the moment, so I manually changed all them to the nucleotide in the PharmVar definition.
***********

***FIXED***
CYP2B6, I changed the multiallelic sites to be single alts.  This isn't a problem in the UKBB, but the other 
alts are common enough that this should be fixed.

g.40991390C>A/T to g.40991390C>T
g.41004406G>A/C/T to g.41004406G>T
g.41010006G>T/A/C to g.41010006G>C
***********
 
Add allele flips for 
CYP2C19: rs3758581
CYP3A5: rs776746

# CFTR hg19 synonym - note that this synonym is not in hg38
# 117199645: 117199644
'''

import json
import argparse
from pyliftover import LiftOver
from Gene import Gene
from DawgToys import iupac_nt
import copy
import myvariant

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

        # Get ref and alts from MyVariant
        #self.update_variants()

        # Add synonynms
        self.add_synonyms()

        # Expand out the wobble nucleotides
        self.dewobbler()

        # Fix all the INDELs
        self.update_indels()

        # Add hg19 position
        self.add_alternate_build()

        # Handle special cases
        self.special_cases()

        # print the output
        print(json.dumps(self.data, indent=2))

    def special_cases(self):

        # todo
        # Expand the INDEL in UGT1A1 into 3 different variants

        # For CYP3A5, *3 needs to have the R and Y changed to G and T

        if self.data["gene"] == "CYP2D6":
            self.CYP2D6()

        pass

    def dewobbler(self):
        # go through each named allele
        new_alleles = []

        for allele in self.data["namedAlleles"]:
            cyp2d6_exon_9_conv = False
            cyp2d6_exon_9_conv_copy = None
            if self.data['gene'] == "CYP2D6" and allele["name"] == "*4":
                cyp2d6_exon_9_conv = True

            branched_alleles = [[]]
            #print(allele['name'])
            for v in range(len(allele['alleles'])):


                #print("Current variant: %s" % v)
                expanded = iupac_nt(allele['alleles'][v])



                #print("Expanded: %s" % expanded)
                if len(expanded) > 1:

                    # create a copy of the original for safe keeping
                    if cyp2d6_exon_9_conv is True and cyp2d6_exon_9_conv_copy is None:
                        # check which variant we're at
                        position = self.data['variants'][v]['position']
                        if position == 42126663:
                            # We're in the conversion.  Stop splitting.
                            #print("We're in the conversion!")
                            cyp2d6_exon_9_conv_copy = copy.deepcopy(branched_alleles)

                    if cyp2d6_exon_9_conv_copy is not None:
                        for i in range(len(branched_alleles)):
                            branched_alleles[i].append(expanded[0])
                        for i in range(len(cyp2d6_exon_9_conv_copy)):
                            cyp2d6_exon_9_conv_copy[i].append(expanded[1])

                    else:
                        new_branches = []
                        for nt in expanded:

                            # create a new copy of all branches and append the new variant
                            new_branch = copy.deepcopy(branched_alleles)

                            #print("original: %s" % branched_alleles)
                            #print(new_branch)
                            for b in range(len(new_branch)):
                                new_branch[b].append(nt)
                                #print(new_branch)

                            new_branches = new_branches + new_branch

                        branched_alleles = new_branches

                else:
                    for i in range(len(branched_alleles)):
                        branched_alleles[i].append(expanded[0])

                    if cyp2d6_exon_9_conv_copy is not None:
                        for i in range(len(cyp2d6_exon_9_conv_copy)):
                            cyp2d6_exon_9_conv_copy[i].append(expanded[0])

            if cyp2d6_exon_9_conv_copy is not None:
                branched_alleles = branched_alleles + cyp2d6_exon_9_conv_copy

            # add the new branches to namedAlleles
            if len(branched_alleles) > 1:
                #print("Created %s branches" % len(branched_alleles))
                for i in range(len(branched_alleles)):
                    new_allele = copy.deepcopy(allele)
                    new_allele['alleles'] = branched_alleles[i]
                    new_allele['id'] = new_allele['id'] + '.%s' % i
                    #print(new_allele)
                    new_alleles.append(new_allele)

        self.data["namedAlleles"] = self.data["namedAlleles"] + new_alleles
        #print(self.data["namedAlleles"])





    def CYP2D6(self):
        if self.debug:
            print("Fixing CYP2D6")

        self.variant_splitter('rs72549356')

        # Fix all the types - no longer need to do this after ryan removed *5
        #self.fix_type()

        #exit()
        #pass

    def variant_splitter(self, rsid):
        # Find the variant in the variants

        index = None

        for i in range(len(self.data['variants'])):
            variant = self.data['variants'][i]
            if variant['rsid'] == rsid:
                index = i
                break

        if index is None:
            print("Variant not found in definition: %s" % rsid)
            exit(1)

        # Update the variant definition
        variant = self.gene.variants[index]
        alts = variant.alt

        if variant.type == "INS":
            for i in range(len(alts)):
                alts[i] = "ins" + alts[i]

        new_variants = []
        for alt in variant.alt:
            new_variant = copy.deepcopy(self.data['variants'][index])

            if variant.type == "INS":

                new_alt = new_variant['chromosomeHgvsName'].split('ins')[0] + alt
                new_variant['chromosomeHgvsName'] = new_alt

            new_variants.append(new_variant)

        del self.data['variants'][index]

        for i in range(len(new_variants)):

            self.data['variants'].insert(index+i, new_variants[i])

        #print(self.data['variants'][index])
        #print(self.data['variants'][index+1])

        # fix variantAlleles
        #print(self.data['variantAlleles'][index])
        del self.data['variantAlleles'][index]
        for i in range(len(alts)):
            new_alts = [alts[i], variant.ref]

            self.data['variantAlleles'].insert(index + i, new_alts)

        #print(self.data['variantAlleles'][index])
        #print(self.data['variantAlleles'][index+1])

        # ok now go through all the named allele definitinos
        for i in range(len(self.data['namedAlleles'])):
            #print(self.data['namedAlleles'][i]['name'])
            if self.data['namedAlleles'][i]['name'] == "*1":
                # skip *1, just because
                self.data['namedAlleles'][i]['alleles'].insert(index, None)
                continue

            if self.data['namedAlleles'][i]['alleles'][index] is None:
                #just add another one
                self.data['namedAlleles'][i]['alleles'].insert(index, None)

            else:
                found_alt = self.data['namedAlleles'][i]['alleles'][index]

                del self.data['namedAlleles'][i]['alleles'][index]
                for j in range(len(alts)):
                    if found_alt == alts[j]:
                        #print('inserting %s' % found_alt)
                        self.data['namedAlleles'][i]['alleles'].insert(index+j, found_alt)
                    else:
                        self.data['namedAlleles'][i]['alleles'].insert(index+j, None)

                #print(self.data['namedAlleles'][i]['alleles'][index])
                #print(self.data['namedAlleles'][i]['alleles'][index+1])

                #exit()

        #exit()


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


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
from DawgToys import iupac_nt


nonmulti_multis = {
    42130692: "A",
    42129819: "T",
    42129809: "C",
    42127941: "A",
    42126663: "C",
    42126660: "C",
    42126647: "T",
    42126636: "A",
    42126635: "G",
    42126633: "G",
    42126627: "C",
    42126624: "T",
    42126623: "C",
    42129075: "T",
    42130715: "T",
    42129042: "C",
    42126611: "G"
}

# 42129071 needs to be split into two.  with A and G as alts.
# all the others that need to be split should have a / in the definition

class CYP2D6Cleaner(object):
    def __init__(self, file, debug=False):
        self.file = file
        self.debug = debug

        # Load the data into a global object
        self.data = self.get_data()

        self.run()

    def run(self):

        # change all the types to the actual type instead of "del"
        self.fix_type()

        # Remove IUPAC codes from definitions (only where unnecessary)
        self.remove_iupac()

        # Split all the multiallelic variants apart
        split_loci = self.fix_variant_defs()

        # delete *5
        self.remove_s5()

        # now update the allele variant lists with the split loci
        self.update_allele_variants(split_loci)

        # print the output
        print(json.dumps(self.data, indent=2))

    def remove_s5(self):
        for i in range(len(self.data['namedAlleles'])):
            if self.data['namedAlleles'][i]["name"] == "*5":
                del self.data['namedAlleles'][i]
                return

    def update_allele_variants(self, split_loci):

        for i in range(len(self.data['namedAlleles'])):
            offset = 0
            copy = self.data['namedAlleles'][i]["alleles"]
            # for each loci that was split, split the according loci into two.
            # in cases where the allele for that star allele matches one of the split loci
            # make sure it ends up in the right index
            for s in split_loci.keys():
                # this is the position in the original allele list
                #print(self.data['namedAlleles'][i]["alleles"][s])
                #print(split_loci[s])

                loci_alleles = []
                for l in split_loci[s]:
                    loci_alleles.append(self.data["variants"][l]["alt"])

                current_allele = self.data['namedAlleles'][i]["alleles"][s + offset]

                if current_allele in loci_alleles:
                    #print(self.data['namedAlleles'][i]['name'])
                    #print(s)
                    # make sure it goes in the right place
                    #print(self.data['namedAlleles'][i]["alleles"][s])

                    #print(split_loci[s])
                    #print(loci_alleles)

                    # Get index where we need to put the matching allele
                    index = loci_alleles.index(current_allele)
                    new = [None] * len(loci_alleles)
                    new = self.update_variant_list(new, index, [current_allele])
                    #print(new)

                    #print(copy)
                    #print(s + offset)
                    copy = self.update_variant_list(copy, s + offset, new)
                    #print(copy)

                else:
                    #print("no changes")
                    # otherwise just fill the new space with Nones

                    copy = self.update_variant_list(copy, s + offset, [None] * len(split_loci[s]))

                offset += len(split_loci[s]) - 1

            #print(len(copy))




    def fix_type(self):
        for i in range(len(self.data['variants'])):
            if "ins" in self.data["variants"][i]['chromosomeHgvsName']:
                self.data["variants"][i]["type"] = "INS"
            elif "del" in self.data["variants"][i]['chromosomeHgvsName']:
                self.data["variants"][i]["type"] = "DEL"
            elif ">" in self.data["variants"][i]['chromosomeHgvsName']:
                self.data["variants"][i]["type"] = "SNP"

    def remove_iupac(self):
        # will need to first get the index of the variant in the variants
        for v in nonmulti_multis.keys():
            index = self.get_variant_index_by_pos(v)

            for i in range(len(self.data['namedAlleles'])):
                if self.data['namedAlleles'][i]["alleles"][index] in ["W", "S", "M", "K", "R", "Y"]:
                    #print("fixing %s" % self.data['namedAlleles'][i]["alleles"][index])
                    self.data['namedAlleles'][i]["alleles"][index] = nonmulti_multis[v]

    def fix_variant_defs(self):
        # for any snps with multiple alts (either by iupac or by /) split them apart
        # also any indels
        split_loci = {}

        variants_copy = self.data['variants'].copy()

        for i in range(len(self.data['variants'])):
            copy_index = i + len(split_loci)
            if "/" in self.data["variants"][i]['chromosomeHgvsName']:

                #print("-----------------------------------------------------------------------")
                #print(self.data["variants"][i])

                # This will not work if there are more than two alts.  Thankfully that is very rare.
                if self.data["variants"][i]['type'] == "SNP":
                    head, alt2 = self.data["variants"][i]['chromosomeHgvsName'].split("/")
                    head2 = head[0:-1] + alt2

                    copy_1 = self.data["variants"][i].copy()
                    copy_1['chromosomeHgvsName'] = head
                    copy_1['alt'] = head[-1]

                    copy_2 = self.data["variants"][i].copy()

                    copy_2['chromosomeHgvsName'] = head2
                    copy_2['alt'] = head2[-1]

                    new = [copy_1, copy_2]

                    variants_copy = self.update_variant_list(variants_copy, copy_index, new)

                    split_loci[i] = [copy_index, copy_index+1]

                elif self.data["variants"][i]['type'] == "INS":
                    head, alt2 = self.data["variants"][i]['chromosomeHgvsName'].split("/")
                    head2 = head.split("ins")[0] + alt2

                    copy_1 = self.data["variants"][i].copy()
                    copy_1['chromosomeHgvsName'] = head
                    copy_1['alt'] = 'ins' + head.split("ins")[-1]

                    copy_2 = self.data["variants"][i].copy()
                    copy_2['chromosomeHgvsName'] = head2
                    copy_2['alt'] = 'ins' + head2.split("ins")[-1]

                    new = [copy_1, copy_2]

                    variants_copy = self.update_variant_list(variants_copy, copy_index, new)

                    split_loci[i] = [copy_index, copy_index + 1]

            # Split IUPAC notation apart
            if self.data["variants"][i]['chromosomeHgvsName'][-1] not in ['A', 'C', 'T', 'G']:
                iupac = self.data["variants"][i]['chromosomeHgvsName'][-1]
                alternate_nts = iupac_nt(iupac)
                new = []
                original = self.data["variants"][i]['chromosomeHgvsName']
                split_loci[i] = []
                index = 0
                for alt in alternate_nts:
                    copy = self.data["variants"][i].copy()
                    copy['chromosomeHgvsName'] = original[0:-1] + alt
                    copy['alt'] = copy['chromosomeHgvsName'][-1]
                    new.append(copy)

                    split_loci[i].append(copy_index + index)
                    index += 1

                variants_copy = self.update_variant_list(variants_copy, copy_index, new)

        #print(split_loci)
        #print(len(self.data["variants"]))
        #print(len(variants_copy))

        self.data["variants"] = variants_copy
        return split_loci

    def update_variant_list(self, variants, index, new_variants):
        # Remove a variant and replace it with new variants in the same position in the list
        del variants[index]
        for v in reversed(new_variants):
            #print(v)
            variants.insert(index, v)
        return variants

    def get_variant_index_by_pos(self, pos):
        for i in range(len(self.data['variants'])):
            if self.data["variants"][i]['position'] == pos:
                return i
        return None

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
    CYP2D6Cleaner(options.file, options.debug)


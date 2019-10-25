

class NamedAllele(object):
    def __init__(self, data, variants, index=0, debug=False):

        # Set global variables
        self.debug = debug
        self.name = None
        self.id = None
        self.function = None
        self.alleles = None
        self.population_frequency = None
        self.index = index
        self.variants = variants
        self.binary_alleles = None
        self.key_list = None

        self._setup(data)

    def _setup(self, data):
        self.name = data["name"]
        self.id = data["id"]
        #self.function = data["function"]
        self.alleles = data["alleles"]
        self.populationFrequency = data["populationFrequency"]

        self._alleles_to_binary()

    def _alleles_to_binary(self):
        self.binary_alleles = []
        self.key_list = []


        # We'll create a binary matrix by iterating over the variants
        # Variants with multiple alternates will get multiple columns
        for v in range(len(self.variants)):
            # Get the variant ref and alts
            variant = self.variants[v]
            ref = variant.ref
            alts = variant.alt
            flipped = variant.flipped

            if flipped is True:
                temp = ref
                ref = alts[0]
                alts = [temp]

            #variant.print_variant()


            # If there are multiple alts iterate over them and check if the listed allele matches the alt
            for a in alts:
                # If it matches the ref, just add a 0
                if self.alleles is None:
                    # No allele listed
                    self.binary_alleles.append(0)
                elif self.alleles[v] == ref:
                    # Matches the reference
                    self.binary_alleles.append(0)
                elif self.alleles[v] == a:
                    # Matches the listed alternate
                    self.binary_alleles.append(1)

                # If it is an INDEL we have to do this differently
                # I am going to set all reference INDELs to "null" in the definition files
                # That way, if it is an indel and not null, then it must be a 1
                elif variant.type in ["INS", "DEL"]:
                    if not self.alleles[v]:
                        # it's an indel loci but a null value was listed
                        self.binary_alleles.append(0)
                    else:
                        # Something was there, assuming alt allele
                        self.binary_alleles.append(1)

                # Null values indicate a ref call
                elif not self.alleles[v]:
                    self.binary_alleles.append(0)

                else:
                    # No match, just make it a zero
                    if self.debug:
                        print("NO MATCH FOUND!")
                        variant.print_variant()
                        print(variant.alleles)
                        print(alts)
                        print(self.alleles[v])
                        print("Current alt: %s" % a)
                    #exit()
                    self.binary_alleles.append(0)

                # Add a variant key of the alt for every column added
                self.key_list.append("%s_%s" % (v, a))


    def print_haplotype(self):
        print("%s %s %s" % (self.name, self.id, self.function))

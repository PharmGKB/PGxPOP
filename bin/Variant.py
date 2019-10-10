import json
from DawgToys import clean_chr, iupac_nt

'''
Define the Variant object
'''
class Variant(object):
    def __init__(self, data, alleles, index=0, build='grch38', debug=False):

        # Initialize global variables
        self.clean_chromosome = None
        self.chromosome = None
        self.position = None
        self.rsid = None
        self.chromosomeHgvsName = None
        self.geneHgvsName = None
        self.proteinNote = None
        self.resourceNote = None
        self.type = None
        self.referenceRepeat = None
        self.key = None
        self.alleles = alleles
        self.ref = None
        self.alt = []
        self.synonyms = []
        self.index = index
        self.build = build.lower()
        self.debug = debug

        self._setup(data, alleles)

    def _setup(self, data, alleles):

        #todo
        # The definition files don't explicitly say what the ref and alt alleles are.  It's in the chromosomeHgvsName
        # but that is annoying to parse.

        self.chromosome = data['chromosome']
        self.clean_chromosome = clean_chr(self.chromosome)
        self.position = data['position']
        self.rsid = data['rsid']
        self.chromosomeHgvsName = data['chromosomeHgvsName']
        self.geneHgvsName = data['geneHgvsName']
        self.proteinNote = data['proteinNote']
        self.resourceNote = data['resourceNote']
        self.type = data['type']
        self.referenceRepeat = data['referenceRepeat']
        self.key = "%s.%s" % (self.chromosome, self.chromosomeHgvsName)
        self.keys = self.chromosomeHgvsName.split("; ")

        if "synonyms" in data.keys():
            self.synonyms = data["synonyms"]

        self._parse_alleles()

        if self.build != 'grch38':
            if self.debug:
                print("Alternate genome build found.  Converting positions.")
            if self.build not in ['hg19']:
                print("Build not stored.  Coordinate conversion may take a while: %s" % self.build)
                self.liftover()
            else:
                self.chromosome = data[self.build]['chromosome']
                self.position = data[self.build]['position']
                self.synonyms = data[self.build]['synonyms']

    def _parse_alleles(self):
        if self.type == "SNP":
            hgvs = self.chromosomeHgvsName.split(";")
            for i in hgvs:
                i = i.strip()
                ref, alt = i[-3:].split('>')
                if self.ref != ref:
                    self.ref = ref
                self.alt = self.alt + iupac_nt(alt)
            return

        if self.type =="DEL":
            hgvs = self.chromosomeHgvsName.split(";")
            for i in hgvs:
                i = i.strip()
                ref = i.split('del')[1]
                if self.ref != ref:
                    self.ref = ref
                self.alt.append("del" + ref)
            return

        if self.type =="INS":
            hgvs = self.chromosomeHgvsName.split(";")
            for i in hgvs:
                i = i.strip()
                alt = i.split('ins')[1]
                self.alt.append(alt)
            self.ref = "ins"
            return

        print("Allele parsing of %s not yet supported" % self.type)
        print(self.chromosomeHgvsName)
        exit()



    def liftover(self):

        # todo
        # Not sure what the failure mode of this tool is.  Will probably need to write a try catch eventually
        # Changing the chromosome and position messes up the key as well.  Could probably fix that.  But i don't have
        # the ref and alt alleles on hand and I don't want to parse them out of chromosomeHgvsName.

        from pyliftover import LiftOver
        lo = LiftOver('hg38', self.build)
        lifted = lo.convert_coordinate(self.chromosome, self.position)

        self.chromosome = lifted[0][0]
        self.position = lifted[0][1]

    def print_variant(self):
        print("%s %s %s %s %s %s" % (self.key, self.chromosome, self.position, self.rsid, self.ref, self.alt))

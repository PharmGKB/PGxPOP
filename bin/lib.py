import tabix
import numpy as np
from gseq import parse_vcf_line, get_vcf_subject_ids, split_genotype

class GenotypeParser(object):
    def __init__(self, vcf, gene, debug=False):
        # Given a VCF and a Gene definition extract all the necessary genotypes and construct haplotype matrices
        self.vcf = vcf
        self.vcf_header = get_vcf_subject_ids(vcf)
        self.gene = gene
        self.debug = debug
        
        # Placeholder for needed property
        self.is_phased = None
        
    def get_sample_index(self, sample_id):
        sample_index = [x for x,i in enumerate(self.vcf_header) if i == sample_id][0]
        return(sample_index)
        
    def haplotype_matrices(self):

        subjects = get_vcf_subject_ids(self.vcf)
        subject_hap_lists = {}





        null_row = []
        for s in subjects:
            subject_hap_lists[s] = [[], []]
            null_row.append(0)

        # For each variant in the gene
        variants = self.gene.variants

        # Variants not called
        uncalled = []

        # This will be a list of lists.  Where the inner lists represent the allele status for every variant in the
        # definition table.  Column represent subjects.  Same order as in the VCF.  Each added row should be of the
        # same length as the number of samples.
        # The left and right lists represent the left and right sides of the genotype
        all_alleles = ([], [])


        # Extract the variant with tabix
        for v in range(len(variants)):
            variant = variants[v]
            genotypes = self.fetch_genotypes(variant)

            # if nothing found, make a row of all zeroes
            # also save all the variants that were not callable
            if genotypes is None:
                uncalled.append(variant)
                # add the same number of null rows as there are genotypes
                for i in range(len(variant.alt)):
                    all_alleles[0].append(null_row)
                    all_alleles[1].append(null_row)
                continue

            # dictionary of alts to store genotypes in
            alt_alleles = {}
            for a in variant.alt:
                alt_alleles[a] = [[], []]

            # Dictionary of allele indices in the VCF
            alt_indices = {}
            for a in variant.alt:
                # Get the index of the alt in the listed genotypes.
                if a in genotypes.alts():
                    index = genotypes.alts().index(a)
                    alt_indices[a] = index + 1

                else:
                    alt_indices[a] = None

            # Loop over every genotype in the row, 1 for each subject
            for s in range(len(subjects)):
                # Get the genotype and split it into an array
                gt = split_genotype(genotypes.calls[s])

                # For each allele in the genotype (2) figure out if it is a ref call, or if it's an alt, which one.
                for g in range(len(gt)):
                    allele = gt[g]
                    if allele == '0':
                        for a in alt_alleles.keys():
                            alt_alleles[a][g].append(0)
                    else:
                        # Determine which alterate allele
                        found = False
                        for alt in alt_indices:
                            if str(alt_indices[alt]) == allele:
                                if self.debug:
                                    print("Found alt call")
                                    print(gt)

                                alt_alleles[alt][g].append(1)

                                # Update the other alleles
                                for other_alt in alt_indices:
                                    if other_alt != alt:
                                        alt_alleles[other_alt][g].append(0)

                                found = True
                        if found is not True:
                            if self.debug:
                                print("Didn't find alt call")
                                print(gt)
                                print(genotypes.ref, genotypes.alts())
                            # Allele does not exist in definitions.  Mark as 0.
                            for a in alt_alleles.keys():
                                alt_alleles[a][g].append(0)

            for a in variant.alt:

                if len(alt_alleles[a][0]) != len(subjects) or len(alt_alleles[a][1]) != len(subjects):
                    print("Incorrect number of alleles!")
                    #print(alt_alleles)
                    print(len(alt_alleles[a][0]))
                    print(len(alt_alleles[a][1]))
                    print("Should be %s" % len(subjects))
                    variant.print_variant()
                    exit()

                all_alleles[0].append(alt_alleles[a][0])
                all_alleles[1].append(alt_alleles[a][1])

        left_hap = np.array(all_alleles[0])
        right_hap = np.array(all_alleles[1])

        return((left_hap, right_hap))


        # construct a list of allele matches for each alternate allele.  Do this in a single pass
        # Make sure that you match the order of the alts in the gene definition


        pass

    def fetch_genotypes(self, variant):
        tb = tabix.open(self.vcf)
        records = tb.query("%s" % variant.clean_chromosome, variant.position-1, variant.position+1)

        for r in records:
            # todo add some additional checks here
            # Otherwise check if there is an rsid, if there is check that it matches the variant file
            # If not just check the position and the alternate alleles
            # If it is an indel, try to match but it might not
            if r[1] == str(variant.position):
                if self.debug:
                    print("Matching variant found in VCF")
                data = parse_vcf_line(r)
                return data

        return None
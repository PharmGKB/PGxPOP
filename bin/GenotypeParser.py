import tabix
import numpy as np
from DawgToys import parse_vcf_line, get_vcf_subject_ids, split_genotype, vcf_is_phased, fetch_genotypes

class GenotypeParser(object):
    def __init__(self, vcf, gene, debug=False):
        # Given a VCF and a Gene definition extract all the necessary genotypes and construct haplotype matrices
        self.vcf = vcf
        self.vcf_header = get_vcf_subject_ids(vcf)
        self.gene = gene
        self.debug = debug
        self.is_phased = vcf_is_phased(vcf)
        
    def get_sample_index(self, sample_id):
        sample_index = [x for x,i in enumerate(self.vcf_header) if i == sample_id][0]
        return(sample_index)
        
    def haplotype_matrices(self, batch_mode=False):

        subjects = get_vcf_subject_ids(self.vcf)
        # IS THIS USED ANYWHERE ELSE? DON'T UNDERSTAND WHY THIS IS CREATED?
        # subject_hap_lists = {}
        #
        #
        null_row = []
        for s in subjects:
        #     subject_hap_lists[s] = [[], []]
            null_row.append(0)

        if batch_mode:
            samps_per_batch = int(np.floor(np.sqrt(len(subjects))))
            batches = [subjects[x:x+samps_per_batch] for x in range(0, len(subjects)-samps_per_batch)]
            batches.append(subjects[len(subjects)-samps_per_batch:])
        else:
            batches = [subjects]

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
        for batch in batches:
            for v in range(len(variants)):
                variant = variants[v]
                if self.debug:

                    variant.print_variant()
                genotypes = fetch_genotypes(self.vcf, variant)

                # if nothing found, make a row of all zeroes
                # also save all the variants that were not callable
                if genotypes is None:
                    uncalled.append(variant)
                    if self.debug:
                        print("Variant not found")
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
                for i,s in enumerate(batch):
                    # Get the genotype and split it into an array
                    gt = split_genotype(genotypes.calls[i])

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
                                    #if self.debug:
                                    #    print("Found alt call")
                                    #    print(gt)

                                    alt_alleles[alt][g].append(1)

                                    # Update the other alleles
                                    for other_alt in alt_indices:
                                        if other_alt != alt:
                                            alt_alleles[other_alt][g].append(0)

                                    found = True
                            if found is not True:
                                # if there is only one alternate allele, we can assume it is right because of the extra
                                # checks when pulling the genotype.  Otherwise, we'll need to add more exceptions to catch it.
                                # This is only true for INDELs.

                                if len(alt_alleles.keys()) == 1:
                                    alt_alleles[list(alt_alleles.keys())[0]][g].append(1)

                                else:
                                    if self.debug:
                                        print("Didn't find alt call")
                                        print(alt_alleles.keys())
                                        print(gt)
                                        print(genotypes.ref, genotypes.alts())
                                        exit()
                                    # Allele does not exist in definitions.  Mark as 0.
                                    for a in alt_alleles.keys():
                                        alt_alleles[a][g].append(0)

                for a in variant.alt:

                    if len(alt_alleles[a][0]) != len(batch) or len(alt_alleles[a][1]) != len(batch):
                        print("Incorrect number of alleles!")
                        #print(alt_alleles)
                        print(len(alt_alleles[a][0]))
                        print(len(alt_alleles[a][1]))
                        print("Should be %s" % len(batch))
                        variant.print_variant()
                        exit()


                    if self.debug:
                        print("%s calls found in A: %s" % (a, np.sum(alt_alleles[a][0])))
                        print("%s calls found in B: %s" % (a, np.sum(alt_alleles[a][1])))


                    all_alleles[0].append(alt_alleles[a][0])
                    all_alleles[1].append(alt_alleles[a][1])

            left_hap = np.array(all_alleles[0])
            right_hap = np.array(all_alleles[1])

            yield((left_hap, right_hap))


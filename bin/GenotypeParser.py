import numpy as np
from DawgToys import parse_vcf_line, get_vcf_subject_ids, split_genotype, vcf_is_phased, fetch_genotypes, is_gt_phased

class GenotypeParser(object):
    def __init__(self, vcf, gene, debug=False):
        # Given a VCF and a Gene definition extract all the necessary genotypes and construct haplotype matrices
        self.vcf = vcf
        self.vcf_header = get_vcf_subject_ids(vcf)
        self.gene = gene
        self.phase_matrix = None
        self.debug = debug
        #self.is_phased = vcf_is_phased(vcf) # no need for this with the command line flag.
        self.is_phased = False
        
    def get_sample_index(self, sample_id):
        sample_index = [x for x,i in enumerate(self.vcf_header) if i == sample_id][0]
        return(sample_index)
        
    def haplotype_matrices(self, batch_mode=False):

        subjects = get_vcf_subject_ids(self.vcf)

        null_row = []
        null_phase_row = []
        for s in subjects:
            null_row.append(0)
            null_phase_row.append(1)

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

        # This will also be a list of lists, which will later be converted to a matrix
        phase_index = []

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

                        phase_index.append(null_phase_row)

                    continue

                # dictionary of alts to store genotypes in
                alt_alleles = {}
                for a in variant.alt:
                    alt_alleles[a] = [[], []]

                # Dictionary of allele indices in the VCF.  This will be the genotype (e.g. 0/1) that we look for
                # For example most alternate alleles will get a value of 1 if they are the only alt.
                alt_indices = {}
                for a in variant.alt:
                #for a in genotypes.alt:
                    #print("Variant alt: %s" % a)
                    # Get the index of the alt in the listed genotypes.

                    #print(genotypes.gt_index)

                    if genotypes.gt_index is not None:
                        alt_indices[a] = genotypes.gt_index + 1

                    elif a in genotypes.alts():
                        #print("ALT: %s:" % a)
                        index = genotypes.alts().index(a)
                        alt_indices[a] = index + 1

                    else:
                        alt_indices[a] = None
                    #print(alt_indices[a])

                # Loop over every genotype in the row, 1 for each subject
                row_phasing_data = []
                for i, s in enumerate(batch):
                    # Get the genotype and split it into an array
                    #print(genotypes.calls[i])
                    gt = split_genotype(genotypes.calls[i])

                    # Check the phasing status first.  Add a zero if it is phased, one if it's not
                    if is_gt_phased(gt):
                        row_phasing_data.append(0)
                    else:
                        row_phasing_data.append(1)

                    # For each allele in the genotype (2) figure out if it is a ref call, or if it's an alt, which one.
                    #for g in range(len(gt)): # for some reason there are 0/0/0 genotypes in the pharmcat test files
                    for g in range(2):
                        allele = gt[g]
                        if allele == '0':
                            for a in alt_alleles.keys():
                                alt_alleles[a][g].append(0)
                        else:
                            # Determine which alterate allele
                            found = False
                            for alt in alt_indices:
                                #print(alt)
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
                                # if there is only one alternate allele, we can assume it is right because of the extra
                                # checks when pulling the genotype.  Otherwise, we'll need to add more exceptions to catch it.
                                # This is only true for INDELs.

                                # This is resolved by including gt_index in the parser
                                #if len(alt_alleles.keys()) == 1:
                                #    pass
                                    #alt_alleles[list(alt_alleles.keys())[0]][g].append(1)

                                #else:
                                #if self.debug:
                                #    print("Didn't find alt call")
                                #    print(alt_alleles.keys())
                                #    print(gt)
                                #    print(genotypes.ref, genotypes.alts())
                                #    exit()
                                # Allele does not exist in definitions.  Mark as 0.
                                for a in alt_alleles.keys():
                                    alt_alleles[a][g].append(0)

                phase_index.append(row_phasing_data)

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

            phase_matrix = np.array(phase_index)

            yield((left_hap, right_hap), phase_matrix)

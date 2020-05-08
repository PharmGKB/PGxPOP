'''
Greg McInes
Altman Lab
gmcinnes@stanford.edu
'''

'''
This script creates matrices of variant calls from a VCF for a gene region.  The output is binary matrices where each
 column is a variant that is included in the haplotype definition and each row is a sample.  Two matrices are output,
 one for each strand.
'''

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
        self.is_phased = False
        self.sample_variants = {}
        self.called = []
        self.variant_list = []

    '''
    Find the index of a sample in the VCF
    '''
    def get_sample_index(self, sample_id):
        sample_index = [x for x, i in enumerate(self.vcf_header) if i == sample_id][0]
        return sample_index

    '''
    Create the matrices
    '''
    def haplotype_matrices(self, batch_mode=False):

        # Get the subjects from the file that we will be operating on
        subjects = get_vcf_subject_ids(self.vcf)

        # Define some defaults to output in the case where no variant is identified
        null_row = []
        null_phase_row = []
        for s in subjects:
            null_row.append(0)
            null_phase_row.append(1)
            self.sample_variants[s] = [[], []]

        # Batch mode setup
        if batch_mode:
            samps_per_batch = int(np.floor(np.sqrt(len(subjects))))
            batches = [subjects[x:x+samps_per_batch] for x in range(0, len(subjects)-samps_per_batch)]
            batches.append(subjects[len(subjects)-samps_per_batch:])
        else:
            batches = [subjects]

        # Get all the variants that are defined for the gene we are working with
        variants = self.gene.variants

        # All alleles identified
        # This will be a list of lists.  Where the inner lists represent the allele status for every variant in the
        # definition table.  Column represent subjects.  Same order as in the VCF.  Each added row should be of the
        # same length as the number of samples.
        # The left and right lists represent the left and right sides of the genotype
        all_alleles = ([], [])

        # Initiate variables to store indices for phasing and no calls
        phase_index = []
        no_call_index = []

        # Batch processing
        for batch in batches:

            # Iterate over all variants and process one by one
            for v in range(len(variants)):
                variant = variants[v]

                if self.debug:
                    variant.print_variant()

                # Extract the variant with tabix
                genotypes = fetch_genotypes(self.vcf, variant)

                # if nothing found, make a row of all zeroes
                # also save all the variants that were not callable
                if genotypes is None:
                    if self.debug:
                        print("Variant not found")
                    # add the same number of null rows as there are genotypes
                    for i in range(len(variant.alt)):
                        all_alleles[0].append(null_row)
                        all_alleles[1].append(null_row)
                        phase_index.append(null_phase_row)
                        no_call_index.append(null_phase_row)
                        new_key = "%s.g.%s%s>%s" % (variant.chromosome, variant.position, variant.ref, variant.alt[i])
                        self.variant_list.append(new_key)
                    continue

                # We found something, so add the variant to the called list
                self.called.append((variant,v))

                # Create a dictionary of alts to store genotypes in
                alt_alleles = {}
                n_alts = 0

                # Since there can be multiallelic sites, we need to iterate over all possible alternates
                # Create matrices for each alt found
                for a in variant.alt:
                    alt_alleles[a] = [[], []]
                    new_key = "%s.g.%s%s>%s" % (variant.chromosome, variant.position, variant.ref, a)
                    self.variant_list.append(new_key)
                    n_alts += 1

                # Dictionary of allele indices in the VCF.  This will be the genotype (e.g. 0/1) that we look for
                # For example most alternate alleles will get a value of 1 if they are the only alt.
                alt_indices = {}
                for a in variant.alt:

                    # Get the index of the alt in the listed genotypes.
                    if genotypes.gt_index is not None:
                        alt_indices[a] = genotypes.gt_index + 1

                    elif a in genotypes.alts():
                        index = genotypes.alts().index(a)
                        alt_indices[a] = index + 1

                    else:
                        alt_indices[a] = None

                # Loop over every genotype in the row, 1 for each subject
                row_phasing_data = []
                row_no_call_data = []

                for i in range(n_alts):
                    row_phasing_data.append([])
                    row_no_call_data.append([])

                for i, s in enumerate(batch):
                    # Get the genotype and split it into an array
                    gt = split_genotype(genotypes.calls[i])

                    # Update the no-call matrix.
                    if "." in gt:
                        if self.debug:
                            print('No call found: %s:%s' % (variant.key, genotypes.calls[i]))
                        for a in range(n_alts):
                            row_no_call_data[a].append(1)
                    else:
                        for a in range(n_alts):
                            row_no_call_data[a].append(0)

                    # Check the phasing status first.  Add a zero if it is phased, one if it's not
                    if is_gt_phased(genotypes.calls[i]):
                        for a in range(n_alts):
                            row_phasing_data[a].append(0)
                    else:
                        for a in range(n_alts):
                            row_phasing_data[a].append(1)

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
                                if str(alt_indices[alt]) == allele:
                                    if self.debug:
                                        print("Found alt call")

                                    alt_alleles[alt][g].append(1)
                                    self.sample_variants[s][g].append(variant.key)

                                    # Update the other alleles
                                    for other_alt in alt_indices:
                                        if other_alt != alt:
                                            alt_alleles[other_alt][g].append(0)

                                    found = True
                            if found is not True:
                                # if there is only one alternate allele, we can assume it is right because of the extra
                                # checks when pulling the genotype.  Otherwise, we'll need to add more exceptions to catch it.
                                # This is only true for INDELs.

                                # Allele does not exist in definitions.  Mark as 0.
                                for a in alt_alleles.keys():
                                    alt_alleles[a][g].append(0)

                for i in range(n_alts):
                    phase_index.append(row_phasing_data[i])
                    no_call_index.append(row_no_call_data[i])

                for a in variant.alt:

                    # Check that we ended up with the proper number of alleles in the matrix
                    if len(alt_alleles[a][0]) != len(batch) or len(alt_alleles[a][1]) != len(batch):
                        print("Incorrect number of alleles!")
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

            # Extract the final matrices
            left_hap = np.array(all_alleles[0])
            right_hap = np.array(all_alleles[1])

            # Prepare the phase and no call matrices
            phase_matrix = np.array(phase_index)
            no_call_matrix = np.array(no_call_index)

            # We yield all these outputs like this to enable parallelization
            yield((left_hap, right_hap), phase_matrix, self.sample_variants, self.variant_list, no_call_matrix)

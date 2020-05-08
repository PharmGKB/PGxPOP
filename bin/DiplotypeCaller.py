"""
DiplotypeCaller.py
Written by Adam Lavertu
05-07-2020

DiplotypeCaller takes in a gene object, as created by Gene.py. Then given a matrix of an individual's haplotypes, the caller
identifies haplotypes that match the pgx haplotypes in the gene object, scores those matches, and returns the top scoring match.
"""
import re

import numpy as np
from itertools import combinations,product

np.seterr(divide='ignore', invalid='ignore')

class DiplotypeCaller(object):

    # Instantiate DiplotypeCaller object
    # Inputs:
    #   gene - Gene.py object
    #   is_phased - Boolean indicating whether the genotype data is phased, if false shuffled haplotype scoring is used
    #   get_partials - Boolean indicating whether partial matches to PGx alleles should be returned
    def __init__(self, gene, is_phased=False, get_partials=True):
        # Get haplotype matrix and associated star allele names from the gene object
        self.hap_matrix, self.stars = gene.haplotype_matrix()
        self.gene = gene

        # Store total number of defined positions for each gene haplotype
        self.hap_alleles_nums = np.sum(self.hap_matrix, axis =1)

        # Find and store the reference allele
        self.ref_allele = self.stars[np.where(self.hap_alleles_nums == 0)[0][0]]

        # Store phased flag
        self.is_phased = is_phased

        #Store partials flag
        self.get_partial_matches = get_partials


        # Create dictionary object mapping star alleles to their corresponding haplotype vectors
        self.star2dip = dict()
        self.variant_list = [var.key for var in gene.variants.values()]
        for j,x in enumerate(self.stars):
            self.star2dip[x] = self.hap_matrix[j,:]

    # check_for_partial_haps
    # Inputs:
    #   comb - list of haplotype combinations
    def check_for_partial_haps(self, comb):
        outPartials = []
        # Subtract assigned haplotypes from the sample haplotype,
        # find variants that are not part of the haplotype allele definitions
        for x in np.where(comb[1] - np.sum([self.star2dip.get(s) for s in comb[0]], axis=0) == 1)[0]:
            outPartials.append(self.variant_list[x])
        # Return list of additional PGx variants that are not part of the haplotype
        return(outPartials)

    # call_diplotype
    # Inputs:
    #   diplotype - list of two vectors, each representing one haplotype of the sample
    #   uncalled - list of uncalled positions to exclude from analysis
    #   phase_vector (optional) - vector indicating phased status of each position (1 = phased, 0 = unphased)
    # Outputs:
    #   string of diplotype calls for the given sample, list of uncallable alleles
    def call_diplotype(self, diplotype, uncalled, phase_vector = None):
        # If data is phased, simply pass through the given diplotype
        if self.is_phased:
            possib_diplotypes = [diplotype]
        # If data is unphased, get list of all possible combinations of the observed variants
        else:
            possib_diplotypes = self.get_possible_diplotypes(*diplotype, phase_vector)

        # Score all possible diplotypes
        dips, dip_scores = self.score_diplotypes(possib_diplotypes)
        
        # Create list of uncallable star alleles, alleles where at least one relevant position is unobserved
        uncallable = [self.stars[x] for x in np.where(np.sum(np.multiply(self.hap_matrix, uncalled.transpose()), axis=1) > 0)[0]]
        uncallable = list(set(uncallable))

        # Find top scoring diplotype
        top_dip_score = np.max(dip_scores)

        out_dips = set()
        # Find all top scoring diplotypes and format the corresponding string output
        for combs in [dips[x] for x in np.where(dip_scores == top_dip_score)[0]]:
            # If the star alleles on a single strand have overlapping variants, adjust the output strings
            strand1 = "or".join(self.process_for_overlaps(combs[0]))
            strand2 = "or".join(self.process_for_overlaps(combs[1]))

            # If data is phased, output with pipe
            if self.is_phased:
                star_dip = "|".join([strand1,strand2])
            else:
                # If data is not phased, sort strings and output with forward slash
                try:
                    star_dip = "/".join(sorted([strand1,strand2], key=lambda a: float(a[1:])))
                except:
                    star_dip = "/".join(sorted([strand1,strand2]))
            out_dips.add(star_dip)
                        
        return([";".join(out_dips), uncallable])

    # score_diplotypes
    # For each possible diplotypes in the nested list of hap_sets, score the star alleles and return
    # Inputs:
    #   hap_sets - nested list of diplotypes, with one pair of haplotypes per nested list
    # Outputs:
    #   list of diplotypes and their associated scores
    def score_diplotypes(self, hap_sets):
        dips = []
        dip_scores = []
        for haps in hap_sets:
            hap1_results = self.get_max_star_alleles(haps[0])
            hap2_results = self.get_max_star_alleles(haps[1])
            tot_dip_score = hap1_results[0] + hap2_results[0]
            dip_scores.append(tot_dip_score)
            dips.append([hap1_results[1:], hap2_results[1:]])
        return([dips, dip_scores])

    # get_possible_diplotypes
    # Input:
    #   haplo1 - vector representation of one haplotype
    #   haplo2 - vector representation of the other haplotype
    #   phase_vector - vector representation of phased state, only unphased positions will be included in novel
    #                   haplotype combinations
    # Output:
    # hap_sets - list of possible haplotype combinations given phasing status
    def get_possible_diplotypes(self, haplo1, haplo2, phase_vector=None):
        # If data is_phased, find positions that are not phased and available for rearrangement
        if self.is_phased:
            temp_hap1 = np.multiply(haplo1,phase_vector)
            temp_hap2 = np.multiply(haplo2,phase_vector)
            alt_sites = np.where(temp_hap1 != temp_hap2)[0]
        else:
            alt_sites = np.where(haplo1 != haplo2)[0]
        hap_sets = []

        # For each alternate allele site create new haplotypes by shuffling the haplotypes at that site
        if len(alt_sites) > 1:
            combos = []
            for i in range((len(alt_sites)**2)//2):
                b = bin(i)[2:].zfill(len(alt_sites))
                combos.append(np.array([int(x) for x in list(b)]))
            for c in combos:
                hap1 = np.zeros(haplo1.shape)
                hap2 = np.zeros(haplo1.shape)
                alt_sites_new = set([alt_sites[x] for x in np.where(c == 1)[0]])
                for x in alt_sites:
                    if x in alt_sites_new:
                        hap1[x] = 1
                    else:
                        hap2[x] = 1
                hap_sets.append([hap1, hap2])
        else:
            hap_sets.append([haplo1, haplo2])
        return(hap_sets)

    # get_max_star_alleles
    # Input:
    #   hap - a single haplotype vector
    # Output:
    # top_score - highest score achieved by one or more pgx haplotypes (total number of matching positions)
    # alleles - list of alleles that achieve that top score
    # hap - input haplotype for downstream functions
    def get_max_star_alleles(self,hap):
        # Calculate the number of hits per haplotype
        hap_hit_nums = np.sum(np.multiply(self.hap_matrix,hap), axis=1)
        
        # Calculate the proportion of alleles matched per haplotype
        hap_prop = (hap_hit_nums/self.hap_alleles_nums)
        
        # If no haplotype matches 100%, return reference
        if np.nanmax(hap_prop) != 1.0:
            top_score = 0.0
            alleles = [self.ref_allele]
        
        else:
            
            top_score = 1.0
            
            # Find all alleles with a 100% match
            matching_alleles = np.argwhere(hap_prop == 1.0)[:,0]
            alleles = {self.stars[x] for x in np.argwhere(hap_prop == 1.0)[:,0]}
            
            # Check for overlapping haplotypes, i.e. those with positions that sum to more than 1
            pot_dips = self.hap_matrix[matching_alleles,:]
            pot_sum = np.sum(pot_dips, axis=0)
            toRemove = set()
            if np.nanmax(pot_sum) > 1:
                
                # Find which positions have overlapping alleles
                cols_oi = np.argwhere(pot_sum > 1)
                
                for pos in cols_oi:
                    # Get haplotypes corresponding to the overlapping alleles
                    overlap_coords = np.squeeze(np.asarray(np.argwhere(self.hap_matrix[:,pos] == 1)[:,0]))
                    
                    # Find the position with the maximum matching and remove the rest
                    overlap_scores = sorted([(self.hap_alleles_nums[x],self.stars[x]) for x in overlap_coords if self.stars[x] in alleles], reverse=True)
                    
                    # Make sure we're not removing duplicate star alleles
                    for x in overlap_scores[1:]:
                        # Check that the allele is not the same as the lead allele and is not tied
                        # for overlap score with the lead allele
                        if x[1] != overlap_scores[0][1] and x[0] != overlap_scores[0][0]:
                            toRemove.add(x[1]) 
                    
            # Filter alleles for the overlaps then return the remaining matching alleles
            alleles = list(alleles.difference(toRemove))
            if "*1" in alleles and len(alleles) > 1:
                _ = alleles.remove("*1")

        return([top_score, alleles, hap])

    # process_for_overlaps
    # Input:
    #   hap_result - a single haplotype vector
    # Output:
    # outSet - set of haplotypes that overlap within the selected haplotypes
    def process_for_overlaps(self, hap_result):
        # Create matrix of the identified haplotypes
        hit_gen_mat = np.matrix([self.star2dip.get(x) for x in hap_result[0]])

        # Find positions not overlapping with the identified haplotypes
        overlap = np.where(np.squeeze(np.array(np.any(np.dot(hit_gen_mat, hit_gen_mat.transpose()) > 1, axis=1))))[0]

        # Find locations where there's an alternate allele and it's already covered by the identified haplotypes
        non_overlap = np.where(np.squeeze(np.array(np.all(np.dot(hit_gen_mat, hit_gen_mat.transpose()) <= 1, axis=1))))[0]
        add_haps = [hap_result[0][j] for j in non_overlap]

        outSet = set()
        partial_vars_hap = []
        # Find overlapping haplotypes and format for output
        if len(overlap) > 0:
            for sep_allele in overlap:
                subset = [hap_result[0][sep_allele]] + add_haps
                if self.get_partial_matches:
                    partial_vars_hap = self.check_for_partial_haps([subset,hap_result[1]])
                subset = [a.split(f"%")[0] for a in subset]
                subset = list(set(subset))
                x = "+".join(sorted(subset, key=lambda a: a[1:]) + partial_vars_hap)
                outSet.add(x)
        else:
            if self.get_partial_matches:
                partial_vars_hap = self.check_for_partial_haps([add_haps,hap_result[1]])
            subset = [a.split(f"%")[0] for a in add_haps]
            subset = list(set(subset))
            x = "+".join(sorted(subset, key=lambda a: a[1:]) + partial_vars_hap)
            outSet.add(x)
        return(outSet)

    # Test function
    def test_gene(self, gene):
        haps, stars = gene.haplotype_matrix()
        indices = [x for x in range(len(stars))]
        true_dips = []
        true_haps = []
        for hap1, hap2 in combinations(indices, 2):
            if self.is_phased:
                true_dips.append("|".join([stars[hap1], stars[hap2]]))
            else:
                true_dips.append("/".join([stars[hap1], stars[hap2]]))
            true_haps.append([haps[hap1], haps[hap2]])
#         len(fake_dips)

        pred_dips = []
        for hap_set in true_haps:
            pred_dips.append(self.call_diplotype(hap_set))
        
        misses = []
        for x in range(len(pred_dips)):
            if self.is_phased and true_dips[x] != pred_dips[x]:
                misses.append([true_dips[x], pred_dips[x]])
            else:
                t = set(re.split(r"[/\|]+", true_dips[x]))
                p = set(re.split(r"[/\|]+", pred_dips[x]))
                if len(p.difference(t)) > 0:
                    misses.append([true_dips[x], pred_dips[x]])
                
        return({"error rate":len(misses)/len(true_dips), "misses":misses, "true_dips":true_dips, "pred_dips":pred_dips})

        

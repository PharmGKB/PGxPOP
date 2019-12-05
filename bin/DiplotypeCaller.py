import re

import numpy as np
from itertools import combinations,product

np.seterr(divide='ignore', invalid='ignore')

class DiplotypeCaller(object):
    def __init__(self, gene, is_phased=False, get_partials=True):
        self.hap_matrix, self.stars = gene.haplotype_matrix()
        self.gene = gene
        self.hap_alleles_nums = np.sum(self.hap_matrix, axis =1)
        self.ref_allele = self.stars[np.where(self.hap_alleles_nums == 0)[0][0]]
        self.is_phased = is_phased
        self.star2dip = dict()
        self.variant_list = None
        for j,x in enumerate(self.stars):
            self.star2dip[x] = self.hap_matrix[j,:]
        self.get_partial_matches = get_partials

    def check_for_partial_haps(self, comb):
        outPartials = []
        for x in np.where(comb[1] - np.sum([self.star2dip.get(s) for s in comb[0]], axis=0) == 1)[0]:
            outPartials.append(self.variant_list[x])
        return(outPartials)

    def call_diplotype(self, diplotype, phase_vector = None):
        if self.is_phased:
            possib_diplotypes = [diplotype]
        else:
            possib_diplotypes = self.get_possible_diplotypes(*diplotype, phase_vector)
            
        dips, dip_scores = self.score_diplotypes(possib_diplotypes)
        top_dip_score = np.max(dip_scores)
        out_dips = set()
#         print(dips)
        for combs in [dips[x] for x in np.where(dip_scores == top_dip_score)[0]]:
            if self.get_partial_matches:
                partial_vars_hap1 = self.check_for_partial_haps(combs[0])
                partial_vars_hap2 = self.check_for_partial_haps(combs[1])
            else:
                partial_vars_hap1 = []
                partial_vars_hap2 = []
            # Add partial matches to list of output alleles
            try:
                x = "+".join(sorted(combs[0][0]))
                x = "+".join(["*" + str(z) for z in sorted([int(z) for z in re.split("\*|\+",x) if len(z) > 0])] + partial_vars_hap1)

                # Access star allele call
                y = "+".join(sorted(combs[1][0]))
                y = "+".join(["*" + str(z) for z in sorted([int(z) for z in re.split("\*|\+",y) if len(z) > 0])] + partial_vars_hap2)
            except:
                # Access star allele calls
                print(partial_vars_hap1)
                print(partial_vars_hap2)
                x = "+".join(sorted(combs[0][0]) + partial_vars_hap1)
                y = "+".join(sorted(combs[1][0]) + partial_vars_hap2)
                                                            
            if self.is_phased:
                star_dip = "|".join([x,y])
            else:
                try:
                    star_dip = "/".join(sorted([x,y], key=lambda a: float(a[1:])))
                except:
                    star_dip = "/".join(sorted([x,y]))
            out_dips.add(star_dip) 
                        
        return(";".join(out_dips))
    
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
    
    
    def get_possible_diplotypes(self, haplo1, haplo2, phase_vector=None):
        if self.is_phased:
            temp_hap1 = np.multiply(haplo1,phase_vector)
            temp_hap2 = np.multiply(haplo2,phase_vector)
            alt_sites = np.where(temp_hap1 != temp_hap2)[0]
        else:
            alt_sites = np.where(haplo1 != haplo2)[0]
        hap_sets = []

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
        
        

        

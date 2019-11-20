import re

import numpy as np
from itertools import combinations,product

np.seterr(divide='ignore', invalid='ignore')

class DiplotypeCaller(object):
    def __init__(self, gene, is_phased=False):
        self.hap_matrix, self.stars = gene.haplotype_matrix()
        self.hap_alleles_nums = np.sum(self.hap_matrix, axis =1)
        self.ref_allele = self.stars[np.where(self.hap_alleles_nums == 0)[0][0]]
        self.is_phased = is_phased

    
        
#     def get_overlapping_haplotypes(self):
#         hap_map_dict = dict()
#         temp_hap_mat = np.matrix(self.hap_matrix)
#         pos_oi = np.where(temp_hap_mat==1)

#         var_sets = dict()
#         for i in range(len(pos_oi[0])):
#             if pos_oi[1][i] not in var_sets:
#                 var_sets[pos_oi[1][i]] = set()
#             var_sets[pos_oi[1][i]].add(pos_oi[0][i])

#         var_oi = sorted([k for k,v in var_sets.items() if len(v) > 1])
#         for j in range(temp_hap_mat[:,var_oi].shape[0]):
#             pos_num = np.sum(temp_hap_mat[j,var_oi])
#             if pos_num > 1:
#                 pos_vars = dict()
#                 for pos_var in np.where(temp_hap_mat[j,:] > 0)[1]:
#                     pos_vars[pos_var] = list(np.where(temp_hap_mat[:,pos_var] == 1)[0])

#                 comb_components = list(pos_vars.values())

#                 for combs in product(*comb_components):
#                     combs = sorted(list(set(combs)))
#                     if len(combs) != 1:
#                         if np.all(np.sum(temp_hap_mat[combs,:], axis=0) == temp_hap_mat[j,:]):
#                             hap_map_dict[";".join([self.stars[k] for k in combs])] = self.stars[j]


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
#             print(combs)
            try:
                x = "+".join(sorted(combs[0], key=lambda a: float(a[1:])))
                y = "+".join(sorted(combs[1], key=lambda a: float(a[1:])))
            except:
                x = "+".join(sorted(combs[0]))
                y = "+".join(sorted(combs[1]))
                
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
            dips.append([hap1_results[1], hap2_results[1]])
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
            
            # Check for overlapping haplotypes, i.e. those with positions that sum to more than 1
            pot_dips = self.hap_matrix[matching_alleles,:]
            pot_sum = np.sum(pot_dips, axis=0)
            toRemove = set()
            if np.nanmax(pot_sum) > 1:
                # Find which positions have overlapping alleles
                cols_oi = np.argwhere(pot_sum > 1)[0]
                
                for pos in cols_oi:
                    # Get haplotypes corresponding to the overlapping alleles
                    overlap_coords = np.squeeze(np.asarray(np.argwhere(self.hap_matrix[:,pos] == 1)[:,0]))
                    
                    # Find the position with the maximum matching and remove the rest
                    overlap_scores = sorted([(self.hap_alleles_nums[x],self.stars[x]) for x in overlap_coords], reverse=True)
                    # MAke sure we're not removing duplicate star alleles
                    toRemove = {x[1] for x in overlap_scores[1:] if x[1] != overlap_scores[0][1]}
                    
            # Filter alleles for the overlaps then return the remaining matching alleles
            alleles = {self.stars[x] for x in np.argwhere(hap_prop == 1.0)[:,0]}
            alleles = list(alleles.difference(toRemove))
            
        return([top_score, alleles])
    
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
        
        

        

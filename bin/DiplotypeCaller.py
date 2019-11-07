import re

import numpy as np
from itertools import combinations

np.seterr(divide='ignore', invalid='ignore')

class DiplotypeCaller(object):
    def __init__(self, gene, is_phased=False):
        self.hap_matrix, self.stars = gene.haplotype_matrix()
        self.hap_alleles_nums = np.sum(self.hap_matrix, axis =1)
        self.ref_allele = self.stars[np.where(self.hap_alleles_nums == 0)[0][0]]
        self.is_phased = is_phased
        # For mixed phasing, give binary vector of position specific phasing info
        # 1 = not-phased, 0 = phased. Multiply by the sumed hap vector
        
    def call_diplotype(self, diplotype, phase_vector = None):
        if self.is_phased:
            possib_diplotypes = [diplotype]
        else:
            possib_diplotypes = self.get_possible_diplotypes(*diplotype, phase_vector)
            
        dips, dip_scores = self.score_diplotypes(possib_diplotypes)
        top_dip_score = np.max(dip_scores)
        out_dips = set()

        for combs in [dips[x] for x in np.where(dip_scores == top_dip_score)[0]]:
            for x in combs[0]:
                for y in combs[1]:
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
        hap_hit_nums = np.sum(np.multiply(self.hap_matrix,hap), axis=1)
        hap_prop = (hap_hit_nums/self.hap_alleles_nums)
        if np.nanmax(hap_prop) != 1.0:
            top_score = 0.0
            alleles = [self.ref_allele]
        else:
            hit_locs = np.zeros(hap_prop.shape)
            hit_locs[np.argwhere(hap_prop == 1.0)] = 1
            hap_scores = np.multiply(hit_locs,self.hap_alleles_nums)
            top_score = np.max(hap_scores)
            alleles = [self.stars[x] for x in np.argwhere(hap_scores == top_score)[:,0]]
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
        
        
        
        
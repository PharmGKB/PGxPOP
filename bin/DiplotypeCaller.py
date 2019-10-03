import numpy as np

class DiplotypeCaller(object):
    def __init__(self, gene, is_phased = False):
        self.hap_matrix, self.stars = gene.haplotype_matrix()
        self.hap_alleles_nums = np.sum(self.hap_matrix, axis =1)
        self.ref_allele = self.stars[np.where(self.hap_alleles_nums == 0)[0][0]]
        self.is_phased = is_phased
        
    def call_diplotype(self, diplotype):
        if self.is_phased:
            possib_haplotypes = diplotype
        else:
            possib_haplotypes = self.get_possible_haplotypes(*diplotype)
        
        dips, dip_scores = self.score_diplotypes(possib_haplotypes)
        top_dip_score = np.max(dip_scores)
        out_dips = list()
        for combs in [dips[x] for x in np.where(dip_scores == top_dip_score)[0]]:
            for x in combs[0]:
                for y in combs[1]:
                    star_dip = "/".join(sorted([x,y]))
                    out_dips.append(star_dip)
        
        return(out_dips[0])
    
    def score_diplotypes(self,hap_sets):
        dips = []
        dip_scores = []
        for haps in hap_sets:
            hap1_results = self.get_max_star_alleles(haps[0])
            hap2_results = self.get_max_star_alleles(haps[1])
            tot_dip_score = hap1_results[0] + hap2_results[0]
            dip_scores.append(tot_dip_score)
            dips.append([hap1_results[1], hap2_results[1]])
        return([dips, dip_scores])
    
    
    def get_possible_haplotypes(self, haplo1, haplo2):
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
        hap_scores = np.multiply((hap_hit_nums/self.hap_alleles_nums),hap_hit_nums)
        top_score = np.nanmax(hap_scores)
        if top_score == 0.0:
            alleles = [self.ref_allele]
        else:
            alleles = [self.stars[x] for x in np.where(hap_scores == top_score)[0]]
        return([top_score, alleles])

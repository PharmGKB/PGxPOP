import numpy as np
    
def get_competitive_haplotype_indices(haplotype_matrix):
    dual_sites = np.where(np.sum(haplotype_matrix, axis=0) > 1)[0]
    competitive_haplotypes = dict()
    for x in haplotype_matrix[:,dual_sites].T:
        comps = set(np.where(x == 1)[0])
        for y in comps:
            if y not in competitive_haplotypes:
                competitive_haplotypes[y] = set()
            competitive_haplotypes[y].update(comps.difference({y}))
    return(competitive_haplotypes)


import numpy as np
import json

class ExceptionCaller(object):
    def __init__(self, gene_exception_file, gene, is_phased = False):
        self.data = json.load(open(gene_exception_file))
        self.is_phased = is_phased
        self.override = self.data['overide_all']
        self.affected_diplotypes = self.data['affected_diplotypes']
        self.build_exception_rules(gene)
        
    def build_exception_rules(self, gene):
        hapCond_to_index_mapping = {}
        for haplotype in self.data['rules']:
            if 'rsids' in self.data['rules'][haplotype]:
                for rsid in self.data['rules'][haplotype]['rsids']:
                    hapCond_to_index_mapping[rsid] = None
                for ind, var in gene.variants.items():
                    if var.rsid in hapCond_to_index_mapping:
                        hapCond_to_index_mapping[var.rsid] = ind
                        
        self.hap_defs = {}
        for hap, rules in self.data['rules'].items():
            temp = {}
            for v in rules.values():
                for var, info in v.items():
                    temp[hapCond_to_index_mapping[var]] = int(info)
            self.hap_defs[hap] = temp
            
            
    def call_samples(self, sample_order, gt_matrices):
        sample_dips = {}
        for gt_mat in gt_matrices:
            for hap, crits in self.hap_defs.items():
                for crit, val in crits.items():
                    
                    for samp in np.where(gt_mat[0][crit] == val)[0]:
                        sid = sample_order[samp]
                        if sid not in sample_dips:
                            sample_dips[sid] = [hap, None]
                        else:
                            sample_dips[sid][0] = hap
                            
                    for samp in np.where(gt_mat[1][crit] == val)[0]:
                        sid = sample_order[samp]
                        if sid not in sample_dips:
                            sample_dips[sid] = [None, hap]
                        else:
                            sample_dips[sid][1] = hap
        
        out_calls = dict()
        print(sample_dips)
        if self.is_phased:
            for sample in sample_dips:
                out_calls[sample] = "|".join(sample_dips[sample])
        else:
            for sample in sample_dips:
                out_calls[sample] = "/".join(sorted(sample_dips[sample]))
                
        return(out_calls)

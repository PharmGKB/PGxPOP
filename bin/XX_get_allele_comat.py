import sys
import pathlib
from collections import Counter
import numpy as np

import Gene
from GenotypeParser import GenotypeParser
from DawgToys import *

gene_name = sys.argv[1]
vcf = sys.argv[2]
oFilePath = sys.argv[3]

gene = Gene.Gene("../definition/alleles/"+gene_name+  "_translation.json")

haps, stars = gene.haplotype_matrix()

gp = GenotypeParser(vcf, gene)

variant_list = list()
for variant in gene.variants.values():
    for a in variant.alt:
        new_key = "%s.g.%s%s>%s" % (variant.chromosome, variant.position, variant.ref, a)
        #print("ALT: %s" % a)
        variant_list.append(new_key)
        
def variant_gen(gp):
    gt_matrices = gp.haplotype_matrices()
    for gt_mat, phase_matrix, sample_vars, variant_list in gt_matrices:
        for samp in range(gt_mat[0].shape[1]):
            yield(gt_mat[0][:,samp])
            yield(gt_mat[1][:,samp])
            
coOcc_data = np.zeros((len(variant_list), len(variant_list)))
for hap in variant_gen(gp):
    var_data = np.where(hap == 1)[0]
    if len(var_data) != 0:
        for i in var_data:
            coOcc_data[i,:] += hap
            
with open(oFilePath, "w+") as outFile:
    _ = outFile.write(",".join(variant_list) + "\n")
    for samp in range(coOcc_data.shape[0]):
        _ = outFile.write(variant_list[samp] + "," + ",".join([str(int(x)) for x in coOcc_data[samp,:]]) + "\n")
            
            
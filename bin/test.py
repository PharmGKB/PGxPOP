
import Gene
from timeit import default_timer as timer
from datetime import timedelta

# your code ...



start = timer()

gene_def = "/Users/gmcinnes/sherlock/mount/CityDawg/definition/alleles/CYP2C9_translation.json"
gene_def = "/oak/stanford/groups/rbaltman/gmcinnes/projects/CityDawg/definition/alleles/CYP2C9_translation.json"

gene = Gene.Gene(gene_def, build='hg19', debug=True)


hap_matrix, stars = gene.haplotype_matrix()
end = timer()

print("Haplotype matrix")
print(hap_matrix.shape)
print(timedelta(seconds=end-start))



from lib import GenotypeParser

start = timer()
vcf = '/Users/gmcinnes/sherlock/mount/CityDawg/data/ALL.chr10_GRCh38.genotypes.20170504.vcf.gz'
vcf = '/oak/stanford/groups/rbaltman/ukbiobank/pgx/imputed/pgx_chr10.vcf.gz'

gp = GenotypeParser(vcf, gene, debug=False)
gt_matrics = gp.haplotype_matrices()
end = timer()


print("GT matrices")
print(gt_matrics[0].shape)
print(gt_matrics[1].shape)
print(timedelta(seconds=end-start))

import numpy as np
print(np.sum(gt_matrics[0]))
print(np.sum(gt_matrics[1]))
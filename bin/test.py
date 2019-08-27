
import Gene

gene = Gene.Gene("/Users/gmcinnes/sherlock/mount/CityDawg/definition/alleles/CYP2C9_translation.json", debug=True)

hap_matrix, stars = gene.haplotype_matrix()

print(hap_matrix)

from lib import GenotypeParser

vcf = '/Users/gmcinnes/sherlock/mount/CityDawg/data/ALL.chr10_GRCh38.genotypes.20170504.vcf.gz'
gp = GenotypeParser(vcf, gene, debug=True)
gp.haplotype_matrices()
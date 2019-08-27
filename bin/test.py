
import Gene

gene = Gene.Gene("/Users/gmcinnes/sherlock/mount/CityDawg/definition/alleles/CYP2C9_translation.json", debug=True)

hap_matrix, stars = gene.haplotype_matrix()

print(hap_matrix)


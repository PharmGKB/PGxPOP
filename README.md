![PGX_POP_logo](logo/PGxPop_logo.png)
## A population-scale pharmacogenetic allele and phenotype caller  
### Greg McInnes and Adam Lavertu  

If you find our code useful, please cite [our paper](https://doi.org/10.1002/cpt.2122).

Pharmacogenetics at scale: An analysis of the UK Biobank<br />
Greg McInnes, Adam Lavertu, Katrin Sangkuhl, Teri E. Klein, Michelle Whirl-Carrillo, Russ B. Altman<br />
Clinical Pharmacology and Therapeutics 2020/11/25; doi: https://doi.org/10.1002/cpt.2122<br />

PGxPOP is a population-scale PGx allele caller designed to handle 100,000s of samples. Input is a phased VCF file, that has been indexed with [tabix](http://www.htslib.org/doc/tabix.html). 
Uses PGx allele definitions as developed by the original [PharmCAT](https://github.com/PharmGKB/PharmCAT) effort.

## Installation:

PGxPOP requires Python 3.6 or greater

First, clone this repository to your local directory.

> git clone https://github.com/PharmGKB/PGxPOP.git  
> cd PGxPOP  
> pip install -r requirements.txt

Once all the necessary requirement files have been installed, run using the following command:

> python bin/PGxPOP.py --vcf <path/to/tabixed_vcf> -g <gene_name> --phased -o <path/to/output_dir>

You can also specify a genome build with the build flag:

> python bin/PGxPOP.py --vcf <path/to/tabixed_vcf> -g <gene_name> --phased --build hg19 -o <path/to/output_dir>



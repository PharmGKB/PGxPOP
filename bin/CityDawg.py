"""
Greg McInnes and Adam Lavertu
Altman Lab
gmcinnes@stanford.edu
"""

import os
import argparse

# import multiprocessing as mp
import numpy as np

from timeit import default_timer as timer
from datetime import timedelta

import Gene
from DawgToys import *
from GenotypeParser import GenotypeParser
from DiplotypeCaller import *
from ExceptionCaller import ExceptionCaller
from Phenotype import Phenotype

class CityDawg(object):
    def __init__(self, vcf, gene, phased=False, build='grch38', output=None,
                 debug=False, batch_mode=False):
        self.vcf = vcf
        self.gene = gene
        self.phased = phased
        self.build = build
        self.debug = debug
        self.batch_mode = batch_mode
        self.output = output

    def run(self):

        # Get the genes we want to run on
        genes = self.get_genes()

        all_results = []
        # For each gene
        for g in genes:
            results = self.process_gene(g)
            all_results = all_results + results

        # Print results to file
        self.print_results(all_results)

    def process_gene(self, g):
        gene = self.get_gene(g)
        gt_matrices = self.get_gt_matrices(gene)
        diplotypes, sample_variants = self.get_calls(gene, gt_matrices)
        phenotypes = self.get_phenotypes(gene, diplotypes, sample_variants)
        return phenotypes
    
    def get_exception_file(self, g):
        definition_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../definition/alleles/")
        filename = "%s_exceptions.json" % g
        if os.path.exists(os.path.join(definition_dir, filename)):
            exception_file = os.path.join(definition_dir, filename)
            return(exception_file)
        else:
            return(False)

    def get_genes(self):
        genes = ['CFTR', 'CYP2C9', 'CYP2D6', 'CYP4F2', 'IFNL3', 'TPMT', 'VKORC1',
                 'CYP2C19', 'CYP3A5',  'DPYD', 'SLCO1B1', 'UGT1A1', 'CYP2B6', 'NUDT15']


        if self.gene == 'all':
            return genes

        if not self.gene in genes:
            print("Selected gene not available.  Please choose from list:")
            print(",".join(genes))
            exit(1)

        return [self.gene]

    def get_gene(self, g):
        # Get the definition file
        gene_definition = get_definition_file(g)
        gene = Gene.Gene(gene_definition, build=self.build, debug=self.debug)
        return gene

    def get_gt_matrices(self, gene):
        if self.debug:
            print("Extracting genotype matrices")

        extraction_start_time = timer()
        gp = GenotypeParser(self.vcf, gene, debug=self.debug)
        gt_matrices = gp.haplotype_matrices(batch_mode=self.batch_mode)
        extraction_end_time = timer()

        if self.debug:
            print("Genotype extraction finished")
            #print("Genotype matrix shapes: %s, %s" % (gt_matrices[0].shape, gt_matrices[1].shape))
            print("Execution time: %s" % timedelta(seconds=extraction_end_time - extraction_start_time))

        return gt_matrices

    def get_calls(self, gene, gt_matrices):

        g = gene.name

        if self.debug:
            print("Calling diplotypes")
            print("Checking for exception file")
        gene_exceptions = self.get_exception_file(g)

        if gene_exceptions:
            ec = ExceptionCaller(gene_exceptions, gene, is_phased=self.phased)
            ec_override = ec.override
            if self.debug:
                print("Exception file found")
        else:
            ec_override = False
            if self.debug:
                print("No exception file found")

        diplotype_caller_start_time = timer()

        sample_variants = {}
        if ec_override:
            if self.debug:
                print("Exception override, calling with ExceptionCaller")
            sample_ids = get_vcf_subject_ids(self.vcf)
            sample_calls, sample_variants = ec.call_samples(sample_ids, gt_matrices)

        else:
            dipCal = DiplotypeCaller(gene, is_phased=self.phased)
            sample_ids = get_vcf_subject_ids(self.vcf)
            sample_calls = {}
            for gt_mat, phase_matrix, sample_vars, variant_list in gt_matrices:
                dipCal.variant_list = variant_list
                for samp in range(gt_mat[0].shape[1]):
                    cd_call = dipCal.call_diplotype([gt_mat[0][:, samp], gt_mat[1][:, samp]], phase_matrix[:,samp])

                    cd_call = self.clean_up_call(cd_call)

                    sample_calls[sample_ids[samp]] = cd_call
                    sample_variants.update(sample_vars)

        diplotype_caller_end_time = timer()

        if self.debug:
            print("Diplotype calling finished")
            print("Execution time: %s" % timedelta(seconds=diplotype_caller_end_time - diplotype_caller_start_time))

            for k, v in sample_calls.items():
                print("%s: %s" % (k, v))

        return sample_calls, sample_variants

    # This is hacky, but wobbles weren't working so I named them differently, like *2%1, *2%2, etc are all *2.
    # As a result of this a *2.1 will get called a *2+*2.1.  So I merge all identical calls here.  We can remove this
    # but I just want to get calls for Katrin that she won't complain about.
    def clean_up_call(self, call):

        if "|" in call:
            calls = call.split("|")
            delim = "|"
        elif "/" in call:
            calls = call.split("/")
            delim = "/"
        else:
            return call
        final_calls = []
        for c in calls:
            subs = c.split("+")
            core_star = []
            for s in subs:
                core = s.split('%')[0]
                if not core in core_star:
                    core_star.append(core)
            final_calls.append("+".join(core_star))
        #final_calls.sort()
        return f'{delim.join(final_calls)}'

    def get_phenotypes(self, gene, diplotypes, sample_variants):
        g = gene.name

        # Load phenotype file into object
        phenotypes = Phenotype()

        results = []

        for sample, diplotype in diplotypes.items():
            haps = split_genotype(diplotype)
            haplotype_data = {}
            for i in range(len(haps)):
                h = haps[i]
                h_function = phenotypes.get_haplotype_function(g, h)
                haplotype_data["hap_%s" % i] = h
                haplotype_data["hap_%s_function" % i] = h_function

            phenotype = phenotypes.get_diplotype_function(g, haps[0], haps[1])

            results.append({
                "sample": sample,
                "gene": g,
                "diplotype": diplotype,
                "hap_1": haplotype_data["hap_0"],
                "hap_2": haplotype_data["hap_1"],
                "hap_1_function": haplotype_data["hap_0_function"],
                "hap_2_function": haplotype_data["hap_1_function"],
                "hap_1_variants": sample_variants[sample][0],
                "hap_2_variants": sample_variants[sample][1],
                "phenotype": phenotype
            })

        return results


    def print_results(self, results):
        f = open(self.output, "w")
        f.write("sample_id,gene,diplotype,hap_1,hap_2,hap_1_function,hap_2_function,hap_1_variants,hap_2_variants,phenotype\n")
        for r in results:
            if self.debug:
                print("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" % (r["sample"], r["gene"], r["diplotype"], r["hap_1"], r["hap_2"],
                                                r["hap_1_function"], r["hap_2_function"], ";".join(r["hap_1_variants"]),
                                                         ";".join(r["hap_2_variants"]), r["phenotype"]))
            f.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (r["sample"], r["gene"], r["diplotype"], r["hap_1"],
                                                         r["hap_2"], r["hap_1_function"], r["hap_2_function"],
                                                         ";".join(r["hap_1_variants"]), ";".join(r["hap_2_variants"]),
                                                         r["phenotype"]))


"""
Parse the command line
"""
def parse_command_line():
    welcome_message()
    parser = argparse.ArgumentParser(
        description = 'CityDawg determines star allele haplotypes for samples in a VCF file and outputs predicted '
                      'pharmacogenetic phenotypes.')
    parser.add_argument("-f", "--vcf", help="Input VCF")
    parser.add_argument("-g", "--gene", default='all', help="Gene to run.  Select from [].  Run all by default.")
    parser.add_argument("--phased", action='store_true', default=False, help="Data is phased.  Will try to determine phasing status "
                                                              "from VCF by default.")
    parser.add_argument("--build", default='grch38', help="Select build genome reference.  By default CityDawg assumes "
                                                            "GRCh38.")
    parser.add_argument("-d", "--debug", action='store_true', default=False,
                                help="Output debugging messages.  May be very verbose.")
    parser.add_argument("-b", "--batch", action='store_true', default=False,
                        help="Fragment into batched sample runs. Suggested for runs with more than 10k samples.")
    parser.add_argument("-o", "--output", default="citydawg_results.txt", help="Output file")
    options = parser.parse_args()
    return options


"""
Main
"""
if __name__ == "__main__":
    options = parse_command_line()
    cd = CityDawg(vcf=options.vcf, gene=options.gene, phased=options.phased, build=options.build, debug=options.debug,
             batch_mode=options.batch, output=options.output)

    cd.run()







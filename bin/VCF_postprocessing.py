#!/usr/bin/env python
import argparse
import sys
try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3

# local imports
from PipelineLib.VariantAnalysis import postprocessing

class DummyLogger:
    def info(self,log_string):
        sys.stderr.write("INFO: {}\n".format(log_string))
    def warning(self,log_string):
        sys.stderr.write("WARNING: {}\n".format(log_string))
    def error(self,log_string):
        sys.stderr.write("ERROR: {}\n".format(log_string))
    def critical(self,log_string):
        sys.stderr.write("CRITICAL: {}\n".format(log_string))
    def debug(self,log_string):
        sys.stderr.write("DEBUG: {}\n".format(log_string))

def main():
    parser=argparse.ArgumentParser(description='VCF postprocessor')
    parser.add_argument('-i','--input',help='Input VCF file name; if unset, stdin is assumed', metavar='VALUE',dest='input_filename',default=None)
    parser.add_argument('-o','--output',help='Output VCF file name; if unset, stdout is assumed', metavar='VALUE',dest='output_filename',default=None)
    parser.add_argument('--ploidy',help='genome ploidy',required=False,default=2,dest='ploidy',type=int)
    parser.add_argument('--variant-calling-region-bed',help='Regions to cover in variant calling',default=None,required=False,dest='variant_calling_bed')
    parser.add_argument('--reference-fasta',help='reference sequence in FASTA format.',dest='reference_fasta',required=False,default=None)
    parser.add_argument('--call-all-positions',help='Call all positions (in target bed file if provided) - so far mpileup and Freebayes only',action='store_true',default=False,required=False,dest='call_all_positions')
    parser.add_argument('--genotype-min-count-het-call-threshold',help='Post-filter het calls genotypes for >= total fraction [default: deactivated with -1] (only applies when --genotype-min-count-filter is used for ploidy=2; non-effective for polyploids)', metavar='VALUE',dest='genotype_min_count_het_call_threshold',default=-1,type=float)
    parser.add_argument('--genotype-min-coverage-threshold',help='Minimum read coverage to keep a genotype value', metavar='VALUE',dest='min_read_threshold',type=int,default=0)
    parser.add_argument('--genotype-min-quality-threshold',help='Minimum genotype quality to keep a genotype value', metavar='VALUE',dest='min_genotype_quality',type=float,default=0.0)
    parser.add_argument('--remove-monomorphic-variants',help='Remove monomorphic loci during post-filtering', dest='remove_monomorphic',default=False,action='store_true')
    args=parser.parse_args()
    
    # First input comes as an argument / stdin
    current_input=args.input_filename
    
    # Append uncalled positions? (need a BED, FASTA and argument)
    if (not args.variant_calling_bed is None) and (not args.reference_fasta is None) and (args.call_all_positions):
        current_output=StringIO()
        postprocessing.append_uncalled_positions_to_VCF(input_vcf=current_input,output_vcf=current_output,target_bed=args.variant_calling_bed,reference_fasta=args.reference_fasta)
        current_output.seek(0)
        current_input=current_output

    # Append IDs? (need a BED)
    if (not args.variant_calling_bed is None):
        current_output=StringIO()
        postprocessing.append_BED_IDs_to_VCF(input_vcf=current_input,output_vcf=current_output,target_bed=args.variant_calling_bed)
        current_output.seek(0)
        current_input=current_output
    
    # Finally, apply read count filter
    if args.ploidy > 2:
        postprocessing.polyploid_Freebayes_VCF_allele_count_hard_filter(input_vcf=current_input, output_vcf=args.output_filename,min_total_count=args.min_read_threshold,logger=DummyLogger(),sample_string="all",remove_monomorphic=args.remove_monomorphic,min_genotype_quality=args.min_genotype_quality)
    else:
        postprocessing.Freebayes_VCF_allele_count_hard_filter(input_vcf=current_input, output_vcf=args.output_filename,min_total_count=args.min_read_threshold,min_het_ab=args.genotype_min_count_het_call_threshold,logger=DummyLogger(),sample_string="all",remove_monomorphic=args.remove_monomorphic,min_genotype_quality=args.min_genotype_quality,ploidy=args.ploidy)
    
    # Done
    return()
    
if __name__ == "__main__":
    main()


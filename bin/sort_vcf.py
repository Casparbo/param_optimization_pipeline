#!/usr/bin/env python3
import subprocess
import argparse

import natsort


def _natsort_vcf_samples(input_vcf_filename,output_vcf_filename):
		# Get sample names from VCF:
		command = "grep -m 1 ^#CHROM "+input_vcf_filename+" | cut -f 10- | sed 's/\t/,/g'"
		vcf_sample_header = subprocess.check_output(args=command,shell=True,universal_newlines=True)
		samples = vcf_sample_header.strip().split(",")

		# Natsort the samples and retrieve their respective index positions
		sorted_samples = natsort.natsorted(samples)
		sorted_sample_indices = []
		for sample in sorted_samples:
				sorted_sample_indices.append(samples.index(sample))

		# Parse file and sort
		infile = open(input_vcf_filename,mode="r")
		outfile = open(output_vcf_filename,mode="w")
		current_line = infile.readline().strip()

		# Parse VCF
		while current_line != "":
				if current_line[0:2] == "##":
						# Header, keep as is
						output_line = current_line
				else:
						# Re-sort
						line_parts = current_line.split("\t")
						output_line = ("\t").join(line_parts[0:9])
						for new_sample_index in sorted_sample_indices:
								output_line += "\t"+ line_parts[new_sample_index+9]

				# Write output
				outfile.write(output_line+"\n")

				# Next line
				current_line = infile.readline().rstrip()

		# Done
		infile.close()
		outfile.close()
		del infile
		del outfile
		return


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("input_file")
	parser.add_argument("output_file")
	args = parser.parse_args()

	_natsort_vcf_samples(args.input_file, args.output_file)


if __name__ == '__main__':
	main()

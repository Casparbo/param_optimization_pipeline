#!/usr/bin/env python3
import argparse
import os

import pandas as pd
import numpy as np


def assemble_sort_df(chromosomes, positions, refs, alts, data_columns):
	"""build a pandas DataFrame out of the given data and sort values by locus"""
	# fill useless columns with zero values
	empty_column = np.zeros_like(positions)
	columns = {"#CHROM": chromosomes, "POS": positions, "ID": empty_column, "REF": refs, "ALT": alts, "QUAL": empty_column,
				"FILTER": empty_column, "INFO": empty_column, "FORMAT": empty_column}
	info_df = pd.DataFrame(data=columns)

	sorted_data_columns = data_columns.sort_index(axis=1)

	df = pd.concat([info_df, sorted_data_columns], axis=1)
	df.sort_values(by=["#CHROM", "POS"], inplace=True)
	df.set_index(["#CHROM", "POS"], inplace=True)

	return df


def switch_call(ref, alt, call):
	"""convert call format to vcf (0/1, 0/0, etc. instead of A/A, A/G, etc)"""
	if call == ref:
		return 0
	elif call == alt:
		return 1
	elif call in ["A", "T", "G", "C"]:
		return 2
	else:
		return call


def txt_to_vcf(df):
	"""switch call nomenclature to fit with vcf format"""
	for i, s in enumerate(df.iloc):
		for j, elem in enumerate(s[3:]):
			if type(elem) is not str or elem=="?/?":
				df.iloc[i, j+3] = "."
				continue
			a, b = elem.split("/")
			df.iloc[i, j+3] = f"{switch_call(s['REF'], s['ALT'], a)}/{switch_call(s['REF'], s['ALT'], b)}"

	return df


def rename_columns(name):
	"""remove suffix from column names so that they are equal to the column names of the sample"""
	name = "".join(name.split("_")[:2])
	return name[:6] + "-" + name[6:]


def determine_alt_base(df):
	"""determine the alt value to allow conversion to vcf-like format"""
	alts_list = []

	for s in df.iloc:
		ref = s["ID"]
		alt = "."

		for elem in s[5:]:
			a, b = elem.split("/")
			if a != ref:
				alt = a
				break
			if b != ref:
				alt = b
				break

		alts_list.append(alt)

	return pd.Series(alts_list)


def parse_txt(input_file):
	"""parse reference file"""
	# parse file with pandas
	txt_df = pd.read_table(input_file, delim_whitespace=True)

	# get the relevant information
	chromosomes = txt_df["Chromsome"]
	positions = txt_df["Stop"]
	refs = txt_df["ID"]
	alts = determine_alt_base(txt_df)

	# get the data
	data_columns = txt_df.iloc[:, 6:]

	# rename data columns to fit with names of vcf file by removing the '_'-seperated suffix
	data_columns.rename(rename_columns, axis=1, inplace=True)

	df = assemble_sort_df(chromosomes, positions, refs, alts, data_columns)
	df = txt_to_vcf(df)

	return df


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("input_file")
	args = parser.parse_args()

	vcf_df = parse_txt(args.input_file)

	with open("converted.vcf", "w") as f:
		vcf_df.to_csv(f, sep="\t")


if __name__ == '__main__':
	main()
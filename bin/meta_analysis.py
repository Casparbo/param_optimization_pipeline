#!/usr/bin/env python3

import argparse

import pandas as pd


def read_data_frames(filelist):
	"""read metadata files and create DataFrames to store them in, return them as a list"""
	df_list = []

	for f in filelist:
		df = pd.read_csv(f)
		# the last two columns are percentage and absolute, so they dont belong in the index, everything else does
		df.set_index(df.columns.to_list()[:-2], inplace=True)
		df_list.append(df)

	return df_list


def combineFrames(df_list):
	"""combine list of metadata DataFrames into a sample DataFrame, a both DataFrame, and a ref DataFrame"""
	combined_df = pd.DataFrame()

	# combine everything into one dataframe
	for df in df_list:
		combined_df = pd.concat([combined_df, df])

	# split off sample
	sample_df = combined_df.xs("sample", level="origin")

	# split off "both"
	both_df = combined_df.xs("both", level="origin")

	# split off ref
	ref_df = combined_df.xs("ref", level="origin")
	# remove param columns since they dont influence the ref
	for n in ref_df.index.names:
		if n != "category":
			ref_df = ref_df.droplevel(level=n)

	ref_df.drop_duplicates(inplace=True)

	return sample_df, both_df, ref_df


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("filelist", nargs="*")
	args = parser.parse_args()

	df_list = read_data_frames(args.filelist)
	sample_df, both_df, ref_df = combineFrames(df_list)


if __name__ == '__main__':
	main()
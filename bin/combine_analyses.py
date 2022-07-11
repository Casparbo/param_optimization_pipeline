#!/usr/bin/env python3

import argparse

import pandas as pd


def read_data_frames(filelist):
	"""read metadata files and create DataFrames to store them in, return them as a list"""
	df_list = []

	for f in filelist:
		df = pd.read_csv(f)
		# the last two columns are percentage and absolute, so they dont belong in the index, everything else does
		df.set_index(df.columns.to_list()[0], inplace=True)
		df_list.append(df)

	return df_list


def combineFrames(df_list):
	"""combine list of metadata DataFrames"""
	combined_df = pd.DataFrame()

	if len(df_list) == 1:
		return df_list[0].xs("average").to_frame().swapaxes("index", "columns")

	# combine everything into one dataframe
	for df in df_list:
		combined_df = pd.concat([combined_df, df])

	return combined_df.xs("average")


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("filelist", nargs="*")
	args = parser.parse_args()

	df_list = read_data_frames(args.filelist)
	average_df = combineFrames(df_list)

	with open("combined.vcf", "w") as f:
		average_df.to_csv(f)


if __name__ == '__main__':
	main()

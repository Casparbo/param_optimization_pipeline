#!/usr/bin/env python3

import argparse
import os

import pandas as pd


def read_data_frames(filelist):
	"""read metadata files and create DataFrames to store them in, return them as a list"""
	df_list = []

	for f in filelist:
		df = pd.read_csv(f, index_col=0, header=[0, 1])
		df_list.append(df)

	return df_list


def combineFrames(df_list):
	"""combine list of metadata DataFrames"""
	if len(df_list) == 1:
		return df_list[0].xs("average").to_frame().swapaxes("index", "columns")

	df_buckets = dict()
	# sort dfs into buckets of the same sample
	for df in df_list:
		sample_name = df.columns[0][0]
		if sample_name in df_buckets:
			df_buckets[sample_name].append(df)
		else:
			df_buckets[sample_name] = [df]

	# combine all dfs from the same sample
	for sample_name in df_buckets:
		combined_df = pd.DataFrame()

		for df in df_buckets[sample_name]:
			combined_df = pd.concat([combined_df, df])

		df_buckets[sample_name] = combined_df.xs("average")

	return df_buckets


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("dirlist", nargs="*")
	args = parser.parse_args()

	filelist = [os.path.join(d, f) for d in args.dirlist for f in os.listdir(d)]

	df_list = read_data_frames(filelist)
	df_buckets = combineFrames(df_list)

	for sample_name, df in df_buckets.items():
		with open(f"{sample_name}_combined.csv", "w") as f:
			df.to_csv(f)


if __name__ == '__main__':
	main()

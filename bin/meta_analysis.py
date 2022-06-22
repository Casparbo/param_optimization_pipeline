#!/usr/bin/env python3

import argparse

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


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


def plot_f1_score(df):
	"""create a scatterplot/swarmplot of the f1 score over params, depending on if they are numbers or not"""
	fig, axs = plt.subplots(2, 3, sharey=True, figsize=(20, 20))
	for i, column in enumerate(df.columns.to_list()[:6]):
		ax = axs.flatten()[i]
		plot_df = pd.concat([df[column], df["f1-score"]], axis=1)
		#print(df[column].astype("category").cat.codes)
		#print(df[column].astype("category").cat.codes.corr(df["f1-score"]))

		if plot_df[column].dtype in ["int64", "float6"]:
			sns.scatterplot(x=column, y="f1-score", data=plot_df, ax=ax)
		else:
			sns.swarmplot(x=column, y="f1-score", data=plot_df, ax=ax)

		ax.set_ylim((0, 1))

	return fig


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("filelist", nargs="*")
	args = parser.parse_args()

	df_list = read_data_frames(args.filelist)
	average_df = combineFrames(df_list)
	figs = [plot_f1_score(average_df)]

	for i, f in enumerate(figs):
		f.savefig(f"metadata_{i}.png")


if __name__ == '__main__':
	main()
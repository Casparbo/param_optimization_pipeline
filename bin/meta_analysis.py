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


def create_3d_param_plot(df, params):
	fig = plt.figure(figsize=(20, 20))
	ax = fig.add_subplot(projection="3d")
	numeric_params = []

	for p in params:
		if df[p].dtype in ["int64", "float64"]:
			numeric_params.append(df[p])
		else:
			numeric_params.append(df[p].astype("category").cat.codes)

	ax.scatter(numeric_params[0], numeric_params[1], df["f1-score"])
	
	ax.set_xlabel(params[0])
	ax.set_ylabel(params[1])
	ax.set_zlabel("f1-score")
	ax.set_zlim((0, 1))

	# set axes labels in case one of them was categorical
	ax.axes.set_xticklabels(df[params[0]])
	ax.axes.set_xticks(numeric_params[0])
	ax.axes.set_yticklabels(df[params[1]])
	ax.axes.set_yticks(numeric_params[1])

	return fig


def plot_f1_score(df):
	"""create a scatterplot/swarmplot of the f1 score over params, depending on if they are numbers or not
	also plot the two params with the biggest correlation in a 3d plot
	"""
	fig, axs = plt.subplots(2, 3, sharey=True, figsize=(20, 20))
	param_corrs = {}
	for i, column in enumerate(df.columns.to_list()[:6]):
		ax = axs.flatten()[i]
		plot_df = pd.concat([df[column], df["f1-score"]], axis=1)

		if plot_df[column].dtype in ["int64", "float64"]:
			sns.scatterplot(x=column, y="f1-score", data=plot_df, ax=ax)
			key = abs((df[column].corr(df["f1-score"])))
		else:
			sns.stripplot(x=column, y="f1-score", data=plot_df, ax=ax)
			key = abs((df[column].astype("category").cat.codes.corr(df["f1-score"])))

		param_corrs[key] = column

		ax.set_ylim((0, 1))

	largest = dict(sorted(param_corrs.items()))

	fig_3d = create_3d_param_plot(df, list(largest.values())[-2:])

	return fig, fig_3d


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("filelist", nargs="*")
	args = parser.parse_args()

	df_list = read_data_frames(args.filelist)
	average_df = combineFrames(df_list)
	figs, fig_3d = plot_f1_score(average_df)

	figs.savefig(f"metadata.png")
	fig_3d.savefig("metadata_3d.png")


if __name__ == '__main__':
	main()
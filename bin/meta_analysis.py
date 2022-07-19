#!/usr/bin/env python3

import argparse
import functools

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf
import seaborn as sns
import numpy as np


def per_sample_stats(df_list, n):
	df_list = [df.copy() for df in df_list]
	top_list = []
	for df in df_list:
		df.columns = df.columns.droplevel(0)
		df.sort_values("f1-score", ascending=False, inplace=True)

		# get only params of top n f1-scores
		top = df.head(n).iloc[:, :-3].reset_index(drop=True)
		top_list.append(top)

	merged = functools.reduce(lambda left, right: left.merge(right, how="inner"), top_list)
	
	return merged


def concat_dfs(df_list):
	for df in df_list:
		df.columns = df.columns.droplevel(0)

	combined = functools.reduce(lambda left, right: pd.concat([left, right]), df_list)
	
	return combined


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
	ax.axes.set_xticks(numeric_params[0])
	ax.axes.set_xticklabels(df[params[0]])
	ax.axes.set_yticks(numeric_params[1])
	ax.axes.set_yticklabels(df[params[1]])

	return fig


def plot_f1_score(df, title):
	"""create a scatterplot/swarmplot of the f1 score over params, depending on if they are numbers or not
	also plot the two params with the biggest correlation in a 3d plot
	"""
	# last three columns are sensitivity, specificity, and f1-score, all others are params
	fig, axs = plt.subplots(len(df.columns[:-3])//3 + 1, 3, sharey=True, figsize=(20, 20))
	fig.suptitle(title, fontsize=20)
	param_corrs = {}
	for i, column in enumerate(df.columns.to_list()[:-3]):
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
	fig_3d.suptitle(title, fontsize=20)

	return fig, fig_3d


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--confusion_vars", nargs="+")
	parser.add_argument("--confusion_vars_missing", nargs="+")
	parser.add_argument("--top_params_thresh", type=int)
	args = parser.parse_args()

	confusion_vars_list = [pd.read_csv(file, index_col=0, header=[0, 1]) for file in args.confusion_vars]
	confusion_vars_missing_list = [pd.read_csv(file, index_col=0, header=[0, 1]) for file in args.confusion_vars]

	top_params = per_sample_stats(confusion_vars_list, args.top_params_thresh)
	top_params_missing = per_sample_stats(confusion_vars_missing_list, args.top_params_thresh)

	confusion_vars = concat_dfs(confusion_vars_list)
	confusion_vars_missing = concat_dfs(confusion_vars_missing_list)

	figs, fig_3d = plot_f1_score(confusion_vars, "F1-score over params, excluding missing data")
	figs_missing, fig_3d_missing = plot_f1_score(confusion_vars_missing, "F1-score over params, including missing data")

	with open("top_params_no_missing.csv", "w") as f:
		top_params.to_csv(f)

	with open("top_params_missing.csv", "w") as f:
		top_params_missing.to_csv(f)

	figs.savefig("metadata_no_missing.png")
	fig_3d.savefig("metadata_3d_no_missing.png")

	figs_missing.savefig("metadata_missing.png")
	fig_3d_missing.savefig("metadata_3d_missing.png")


if __name__ == '__main__':
	main()

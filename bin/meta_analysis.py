#!/usr/bin/env python3

import argparse

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


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
	# last three columns are sensitivity, specificity, and f1-score, all others are params
	fig, axs = plt.subplots(len(df.columns[:-3])//3 + 1, 3, sharey=True, figsize=(20, 20))
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

	return fig, fig_3d


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("confusion_vars")
	parser.add_argument("confusion_vars_missing")
	args = parser.parse_args()

	confusion_vars = pd.read_csv(args.confusion_vars)
	confusion_vars_missing = pd.read_csv(args.confusion_vars_missing)

	figs, fig_3d = plot_f1_score(confusion_vars)
	figs_missing, fig_3d_missing = plot_f1_score(confusion_vars_missing)

	figs.savefig(f"metadata.png")
	fig_3d.savefig("metadata_3d.png")

	figs_missing.savefig("metadata_missing.png")
	fig_3d_missing.savefig("metadata_3d_missing.png")


if __name__ == '__main__':
	main()

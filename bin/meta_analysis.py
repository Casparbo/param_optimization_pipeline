#!/usr/bin/env python3

import argparse
import functools

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns


def table_from_df(df):
	"""create a 2d-list out of a pandas dataframe, including index, excluding column names"""
	cell_text = []
	for row in range(len(df)):
		cell_text.append([df.index[row]] + df.iloc[row].to_list())

	return cell_text


def param_occurence_in_percentage(df_list, percentage):
	df_list = [df.copy() for df in df_list]
	top_list = []

	# prep dfs
	for df in df_list:
		df.columns = df.columns.droplevel(0)
		# sort by param values first in case f1-scores occur multiple times
		df.sort_values(df.columns.to_list()[:-3], ascending=False, inplace=True)
		df.sort_values("f1-score", ascending=False, inplace=True)
	
	thresh_f1_score = max([df["f1-score"].max() for df in df_list]) * percentage / 100

	# find param combinations above threshold of max
	for df in df_list:
		# get only params of top percentage f1-scores
		#thresh_f1_score = df["f1-score"].max() * percentage / 100
		top = df[df["f1-score"] > thresh_f1_score].reset_index(drop=True).drop(["sensitivity", "specificity"], axis=1)
		top_list.append(top)

	combined = functools.reduce(lambda left, right: pd.concat([left, right]), top_list)
	count_df = combined.groupby(combined.columns.to_list()[:-1], as_index=False).size()
	count_df = count_df.rename({"size": percentage}, axis="columns").sort_values(percentage, ascending=False).reset_index(drop=True)

	return count_df


def per_sample_stats(df_list, threshold, title):
	occurence_list = []
	for percentage in range(101):
		occurence_list.append(param_occurence_in_percentage(df_list, percentage))

	occurence_df = functools.reduce(lambda left, right: pd.concat([left, right[right.columns[-1]]], axis=1), occurence_list)

	# too many results lead to overcrowded graphs
	top_params = occurence_df[occurence_df[threshold] > 0].head(20)
	cutoffs = pd.DataFrame()
	param_length = 0

	for i, row in top_params.iterrows():
		params = row[:-101]
		data = row[-101:]
		param_length = len(params)
		changes = pd.Series([perc for perc in data.index[1:] if data[perc] != data[perc-1]])
		new_row = pd.concat([params, changes])
		cutoffs = cutoffs.append(new_row, ignore_index=True)

	data_length = len(cutoffs.columns[param_length:])
	cutoffs.columns = cutoffs.columns[:param_length].to_list() + [f"occurs {data_length - x} times" for x in range(data_length)]

	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(50, 20))
	fig.suptitle(title, fontsize=20)
	fig.tight_layout()
	fig.subplots_adjust(top=0.8)

	legend_patches = []

	# ignore actual param values in plotting, just occurences
	for column, color in zip(cutoffs.iloc[:, param_length:].columns[::-1], sns.color_palette()):
		sns.barplot(x=cutoffs.index, y=cutoffs[column], color=color, ax=ax1)
		legend_patches.append(mpatches.Patch(color=color, label=column))

	ax1.legend(handles=legend_patches)
	ax1.set_ylabel("percentage")
	ax1.set_ylim((90, 100))

	ax2.axis("off")
	table = ax2.table(cellText=table_from_df(cutoffs), colLabels=["index"]+cutoffs.columns.to_list(), loc="center")
	table.auto_set_font_size(False)
	table.set_fontsize(10)
	table.scale(1, 1.5)

	return fig


def concat_dfs(df_list):
	"""concatenate dfs into one, removing the sample name"""
	for df in df_list:
		df.columns = df.columns.droplevel(0)

	combined = functools.reduce(lambda left, right: pd.concat([left, right]), df_list)
	
	return combined


def create_3d_param_plot(df, params):
	"""create a 3d plot of the f1-score over the two params"""
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

	top_params = per_sample_stats(confusion_vars_list, args.top_params_thresh,
		"Top percentage where param combinations occur, excluding missing data")
	top_params_missing = per_sample_stats(confusion_vars_missing_list, args.top_params_thresh,
		"Top percentage where param combinations occur, including missing data")

	confusion_vars = concat_dfs(confusion_vars_list)
	confusion_vars_missing = concat_dfs(confusion_vars_missing_list)

	figs, fig_3d = plot_f1_score(confusion_vars, "F1-score over params, excluding missing data")
	figs_missing, fig_3d_missing = plot_f1_score(confusion_vars_missing, "F1-score over params, including missing data")

	top_params.savefig("top_params_no_missing.png", bbox_inches="tight")
	top_params_missing.savefig("top_params_missing.png", bbox_inches="tight")

	figs.savefig("metadata_no_missing.png")
	fig_3d.savefig("metadata_3d_no_missing.png")

	figs_missing.savefig("metadata_missing.png")
	fig_3d_missing.savefig("metadata_3d_missing.png")


if __name__ == '__main__':
	main()

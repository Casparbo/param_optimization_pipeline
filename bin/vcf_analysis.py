#!/usr/bin/env python3
import argparse

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def assemble_sort_df(chromosomes, positions, refs, alts, data_columns):
	"""build a pandas DataFrame out of the given data and sort values by locus"""
	columns = {"chromosome": chromosomes, "position": positions, "ref": refs, "alt": alts}
	info_df = pd.DataFrame(data=columns)

	sorted_data_columns = data_columns.sort_index(axis=1)

	df = pd.concat([info_df, sorted_data_columns], axis=1)
	df.sort_values(by=["chromosome", "position"], inplace=True)
	df.set_index(["chromosome", "position"], inplace=True)

	return df


def vcf_get_data(data):
	"""get the genome type out of the vcf-data-column"""
	elems = data.split(":")

	return elems[0]


def parse_vcf(input_file):
	"""parse a vcf file and create a pandas DataFrame out of it"""
	# count header lines to be ignored by pandas
	header_count = 0
	with open(input_file) as f:
		while f.readline()[:2] == "##":
			header_count += 1

	# parse file with pandas, ignoring header
	vcf_df = pd.read_table(input_file, skiprows=header_count, delim_whitespace=True)

	# get the relevant information
	chromosomes = vcf_df["#CHROM"]
	positions = vcf_df["POS"]
	refs = vcf_df["REF"]
	alts = vcf_df["ALT"]

	# get the data and extract the relevant part
	data_columns = vcf_df.iloc[:, 9:]
	new_data_columns = data_columns.applymap(vcf_get_data)

	return assemble_sort_df(chromosomes, positions, refs, alts, new_data_columns)


def remove_differing_alts(sample_df, ref_df):
	different_alts = 0

	# remove rows where alts are different between ref and sample
	for chrom, pos in sample_df.index:
		if (chrom, pos) not in ref_df.index:
			sample_df.drop(index=(chrom, pos), inplace=True)
		elif sample_df.loc[chrom].loc[pos]["alt"] != ref_df.loc[chrom].loc[pos]["alt"]:
			different_alts += 1
			ref_df.drop(index=(chrom, pos), inplace=True)
			sample_df.drop(index=(chrom, pos), inplace=True)

	return different_alts


def count_states(sample_df):
	# count states row- and position-wise
	sample_state_count = pd.DataFrame()
	position_state_count = pd.DataFrame()
	
	for column in sample_df.columns[2:]:
		sample_state_count = pd.concat([sample_state_count, sample_df[column].value_counts(dropna=False)], axis=1)

	for i, row in sample_df.iterrows():
		position_state_count = pd.concat([position_state_count, row[2:].value_counts(dropna=False)], axis=1)

	normalize = lambda s: s.apply(lambda e: e/s.sum())

	sample_state_norm = sample_state_count.apply(normalize, axis=1)
	position_state_norm = position_state_count.apply(normalize, axis=1)

	combine_count_norm = lambda x, y: x.combine(y, lambda a, b: f"{a}|{b:%}")

	sample_state_df = sample_state_count.combine(sample_state_norm, combine_count_norm)
	position_state_df = position_state_count.combine(position_state_norm, combine_count_norm)

	return sample_state_df, position_state_df


def count_differences(sample_df, ref_df):
	# build difference DataFrame and count variations
	comb_df = pd.DataFrame()
	sample_count_df = pd.DataFrame()

	# build combination dataFrame
	for column in sample_df.columns[2:]:
		new_column = sample_df[column].combine(ref_df[column], lambda x, y: f"{y}>{x}", fill_value=".")
		comb_df[column] = new_column

	# count combinations per sample
	for column in comb_df.columns:
		sample_count_df = pd.concat([sample_count_df, comb_df[column].value_counts(dropna=False)], axis=1)

	total_counts = sample_count_df.sum(axis=1).sort_values(ascending=False).astype("int32")
	norm_total_counts = total_counts.apply(lambda val: (val*100)/total_counts.sum())
	total_count_df = pd.concat([total_counts, norm_total_counts], axis=1, keys=("absolute", "percentage"))

	return total_count_df


def quantify_differences(sample_df, ref_df):
	"""collect data about the differences in calls between ref and sample"""
	
	# check if everything is the same
	if sample_df.equals(ref_df):
		print("Dataframes are equal!")
		return

	different_alts = remove_differing_alts(sample_df, ref_df)

	sample_state_dfs = count_states(sample_df)
	ref_state_dfs = count_states(ref_df)

	total_count_df = count_differences(sample_df, ref_df)
	
	return total_count_df, sample_state_dfs, ref_state_dfs, different_alts


def build_matrix(count_df):
	"""build a difference matrix out of the result from quantify_differences"""
	ref_idx = set()
	sample_idx = set()

	for idx in count_df.index:
		a, b = idx.split(">")
		ref_idx.add(a)
		sample_idx.add(b)

	perc_df = pd.DataFrame(columns=sorted(sample_idx), index=sorted(ref_idx))
	abs_df = pd.DataFrame(columns=sorted(sample_idx), index=sorted(ref_idx))

	# percentage matrix
	for idx, elem in zip(count_df.index, count_df["percentage"]):
		a, b = idx.split(">")
		perc_df[b][a] = elem

	# absolute value matrix
	for idx, elem in zip(count_df.index, count_df["absolute"]):
		a, b = idx.split(">")
		abs_df[b][a] = elem

	return perc_df, abs_df


def build_heatmap(matrix_df):
	"""build a heatmap out of a difference matrix and save as heatmap.png"""
	matrix_df = matrix_df.astype("float32")
	
	fig = plt.figure(figsize=(10, 10))

	heatmap = sns.heatmap(matrix_df, cmap="Reds", annot=True)
	heatmap.set(xlabel="Sample", ylabel="Reference")
	
	fig.add_axes(heatmap)

	return fig


def parse_param_string(param_string):
	"""parse the param string into a list of params and a list of names"""
	words = param_string.split("_")

	params = [w for w in words if not w.isupper()]
	param_names = [f"param{i}" for i in range(len(params))]

	return params, param_names
	

def calc_metadata(percentage_df, absolute_df, different_alts, params, param_names):
	"""calculate total number of dots, hetero calls and homo calls in sample and ref"""
	# sums of dots
	sample_dots_perc = percentage_df["."].sum()
	sample_dots_abs = absolute_df["."].sum()
	ref_dots_perc = percentage_df.loc["."].sum()
	ref_dots_abs = absolute_df.loc["."].sum()

	# sums of hets
	sample_hets_perc = percentage_df["0/1"].sum()
	sample_hets_abs = absolute_df["0/1"].sum()
	ref_hets_perc = percentage_df.loc["0/1"].sum()
	ref_hets_abs = absolute_df.loc["0/1"].sum()

	# sums of homs
	sample_homs_perc = percentage_df["1/1"].sum() + percentage_df["0/0"].sum()
	sample_homs_abs = absolute_df["1/1"].sum() + absolute_df["0/0"].sum()
	ref_homs_perc = percentage_df.loc["1/1"].sum() + percentage_df.loc["0/0"].sum()
	ref_homs_abs = absolute_df.loc["1/1"].sum() + absolute_df.loc["0/0"].sum()

	idx = pd.MultiIndex.from_tuples(tuples = [params + ["sample", "dots"],
												params + ["sample", "hets"],
												params + ["sample", "homs"],
												params + ["ref", "dots"],
												params + ["ref", "hets"],
												params + ["ref", "homs"],
												params + ["both", "different alts"]],
									names = param_names + ["origin", "category"])
	
	# leave a zero at the end because percentages are calculated without different_alts value
	perc_data = pd.Series(data=[sample_dots_perc, sample_hets_perc, sample_homs_perc,
					ref_dots_perc, ref_hets_perc, ref_homs_perc, 0], index=idx, name="percentage")
	abs_data = pd.Series([sample_dots_abs, sample_hets_abs, sample_homs_abs,
					ref_dots_abs, ref_hets_abs, ref_homs_abs, different_alts], index=idx, name="absolute")
	metadata_df = pd.concat((perc_data, abs_data), axis=1)

	return metadata_df


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("input_file")
	parser.add_argument("comparison_file")
	parser.add_argument("param_string")
	args = parser.parse_args()

	sample_df = parse_vcf(args.input_file)
	ref_df = parse_vcf(args.comparison_file)	

	analysis_result, sample_state_dfs, ref_state_dfs, different_alts = quantify_differences(sample_df, ref_df)
			
	percentage_matrix, absolute_matrix = build_matrix(analysis_result)

	heatmap = build_heatmap(percentage_matrix)

	params, param_names = parse_param_string(args.param_string)
	metadata = calc_metadata(percentage_matrix, absolute_matrix, different_alts, params, param_names)

	with open("input_sample_counts.csv", "w") as f:
		sample_state_dfs[0].to_csv(f)

	with open("sample_position_counts.csv", "w") as f:
		sample_state_dfs[1].to_csv(f)

	with open("ref_sample_counts.csv", "w") as f:
		ref_state_dfs[0].to_csv(f)

	with open("ref_position_counts.csv", "w") as f:
		ref_state_dfs[1].to_csv(f)

	with open("metadata.csv", "w") as f:
		metadata.to_csv(f)

	heatmap.savefig("heatmap.png")

	with open("matrix.csv", "w") as f:
		percentage_matrix.to_csv(f)


if __name__ == '__main__':
	main()

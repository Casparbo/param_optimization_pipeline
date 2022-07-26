#!/usr/bin/env python3
import argparse

import pandas as pd


def assemble_sort_df(chromosomes, positions, refs, alts, data_columns):
	"""build a pandas DataFrame out of the given data and sort values by locus"""
	columns = {"chromosome": chromosomes, "position": positions, "ref": refs, "alt": alts}
	info_df = pd.DataFrame(data=columns)

	sorted_data_columns = data_columns.sort_index(axis=1)

	df = pd.concat([info_df, sorted_data_columns], axis=1)
	df.sort_values(by=["chromosome", "position"], inplace=True)
	df.set_index(["chromosome", "position"], inplace=True)

	#remove duplicates
	df.drop_duplicates(inplace=True)

	return df


def vcf_get_data(data):
	"""get the genome type out of the vcf-data-column"""
	if type(data) is not str:
		return data
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
	"""remove rows in both dfs where alts are different between ref and sample"""
	different_alts = 0
	sample_to_drop = set()
	ref_to_drop = set()

	for chrom, pos in sample_df.index:
		if (chrom, pos) not in ref_df.index:
			sample_to_drop.add((chrom, pos))
		# some variant callers make different lines for multiple alt, resulting in a series being returned here
		elif type(sample_df.loc[chrom].loc[pos]["alt"]) is not str:
			if sample_df.loc[chrom].loc[pos]["alt"].all() != ref_df.loc[chrom].loc[pos]["alt"]:
				sample_to_drop.add((chrom, pos))
				ref_to_drop.add((chrom, pos))
		elif sample_df.loc[chrom].loc[pos]["alt"] != ref_df.loc[chrom].loc[pos]["alt"]:
			sample_to_drop.add((chrom, pos))
			ref_to_drop.add((chrom, pos))

	for chrom, pos in sample_to_drop:
		sample_df.drop(index=(chrom, pos), inplace=True)

	for chrom, pos in ref_to_drop:
		ref_df.drop(index=(chrom, pos), inplace=True)

	return len(ref_df)


def count_states(sample_df):
	"""count states row- and position-wise"""
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


def build_state_dfs(sample_df, ref_df):
	"""collect data about the differences in calls between ref and sample"""
	
	# check if everything is the same
	if sample_df.equals(ref_df):
		print("Dataframes are equal!")
		return

	sample_state_dfs = count_states(sample_df)
	ref_state_dfs = count_states(ref_df)
	
	return sample_state_dfs, ref_state_dfs


def build_transition_df(sample_df, ref_df):
	"""build difference DataFrame and count variations"""
	comb_df = pd.DataFrame()
	sample_count_df = pd.DataFrame()

	# build combination dataFrame
	for column in sample_df.columns[2:]:
		new_column = sample_df[column].combine(ref_df[column], lambda x, y: f"{y}>{x}", fill_value=".")
		comb_df[column] = new_column

	# count combinations per sample
	for column in comb_df.columns:
		sample_count_df = pd.concat([sample_count_df, comb_df[column].value_counts(dropna=False)], axis=1)

	total_counts = sample_count_df#.sum(axis=1).sort_values(ascending=False).astype("int32")
	
	# fill missing combinations, except for points
	calls = set()
	for idx in total_counts.index:
		ref, samp = idx.split(">")
		calls.add(ref)
		calls.add(samp)
	
	for t in [f"{r}>{s}" for r in calls for s in calls]:
		if t not in total_counts.index:
			fill_df = pd.DataFrame(index=[t], data={col: 0 for col in total_counts.columns})
			total_counts = pd.concat([total_counts, fill_df])

	total_counts.fillna(0, inplace=True)

	percentage_total_counts = total_counts.apply(lambda col:  col.apply(lambda val: (val*100)/col.sum()))

	return total_counts, percentage_total_counts


def build_transition_matrices(count_df):
	"""build a transition matrix out of the transition_df"""
	ref_idx = set()
	sample_idx = set()

	for idx in count_df.index:
		a, b = idx.split(">")
		ref_idx.add(a)
		sample_idx.add(b)

	# each column (=sample) gets its own confusion matrix
	matrices = [pd.DataFrame(columns=sorted(sample_idx), index=sorted(ref_idx)) for col in count_df.columns]

	for i, col in enumerate(count_df.columns):
		for idx in count_df.index:
			a, b = idx.split(">")
			matrices[i][b][a] = count_df[col][idx]

	return matrices


def calc_confusion_variables(matrix_df, params, param_names, include_missing=False):
	"""calculate sensitivity, specificity and f1-score of the transition matrix"""
	# omit "." calls?
	if include_missing:
		no_points_matrix = matrix_df
	else:
		no_points_matrix = matrix_df.drop(".", axis=0).drop(".", axis=1)

	confusion_df = pd.DataFrame(index=no_points_matrix.index.to_list() + ["average"],
								columns=param_names+["sensitivity", "specificity", "f1-score"])
	confusion_df.index.name = "call"
	
	# calc variables for each call
	for idx in no_points_matrix.index:
		tp = no_points_matrix.loc[idx].loc[idx]
		p = no_points_matrix.loc[idx].sum()
		sensitivity = tp/p if p != 0 else 0

		tn = no_points_matrix.drop(idx, axis=0).drop(idx, axis=1).to_numpy().sum()
		n = no_points_matrix.drop(idx, axis=0).to_numpy().sum()
		specificity = tn/n if n!= 0 else 0

		fp = no_points_matrix[idx].drop(idx).sum()
		precision = tp/(tp+fp) if tp != 0 or fp != 0 else 0
		
		fn = no_points_matrix.loc[idx].drop(idx)
		f1 = (2*precision*sensitivity)/(sensitivity+precision) if sensitivity != 0 and precision != 0 else 0

		confusion_df.loc[idx] = params + [sensitivity, specificity, f1]

	# average is sum / len -1 because "average" is already part of the index and is thus included in len
	sensitivity_avg = confusion_df["sensitivity"].sum()/(len(confusion_df["sensitivity"]) - 1)
	specificity_avg = confusion_df["specificity"].sum()/(len(confusion_df["specificity"]) - 1)
	f1_avg = confusion_df["f1-score"].sum()/(len(confusion_df["f1-score"]) - 1)
	confusion_df.loc["average"] = params + [sensitivity_avg, specificity_avg, f1_avg]

	return confusion_df


def parse_param_string(param_string, param_names):
	"""parse the param string into a list of params and a list of names"""
	words = param_string.split("_")

	params = [w for w in words if not w.isupper()]
	param_names = param_names + [f"param{i}" for i in range(len(params) - len(param_names))]

	return params, param_names


def sum_by_call(df, call, loc=False):
	if loc and (call in df.index):
		return df.loc[call].sum()
	elif (not loc) and (call in df):
		return df[call].sum()
	else:
		return 0
	

def calc_metadata(percentage_df, absolute_df, different_alts, params, param_names):
	"""calculate total number of dots, hetero calls and homo calls in sample and ref, as well as the difference between them (delta)"""
	# sums of dots
	sample_dots_perc = sum_by_call(percentage_df, ".")
	sample_dots_abs = sum_by_call(absolute_df, ".")
	ref_dots_perc = sum_by_call(percentage_df, ".", True)
	ref_dots_abs = sum_by_call(absolute_df, ".", True)

	delta_dots_perc = abs(ref_dots_perc - sample_dots_perc)
	delta_dots_abs = abs(ref_dots_abs - sample_dots_abs)

	sample_hets_perc = sum_by_call(percentage_df, "0/1")
	sample_hets_abs = sum_by_call(absolute_df, "0/1")
	ref_hets_perc = sum_by_call(percentage_df, "0/1", True)
	ref_hets_abs = sum_by_call(absolute_df, "0/1", True)

	delta_hets_perc = abs(ref_hets_perc - sample_hets_perc)
	delta_hets_abs = abs(ref_hets_abs - sample_hets_abs)

	# sums of homs
	sample_homs_perc = sum_by_call(percentage_df, "1/1")+ sum_by_call(percentage_df, "0/0")
	sample_homs_abs = sum_by_call(absolute_df, "1/1") + sum_by_call(absolute_df, "0/0")

	ref_homs_perc = sum_by_call(percentage_df, "1/1", True) + sum_by_call(percentage_df, "0/0", True)
	ref_homs_abs = sum_by_call(absolute_df, "1/1", True) + sum_by_call(absolute_df, "0/0", True)

	delta_homs_perc = abs(ref_homs_perc - sample_homs_perc)
	delta_homs_abs = abs(ref_homs_abs - sample_homs_abs)

	idx = pd.MultiIndex.from_tuples(tuples = [params + ["sample", "dots"],
												params + ["sample", "hets"],
												params + ["sample", "homs"],
												params + ["ref", "dots"],
												params + ["ref", "hets"],
												params + ["ref", "homs"],
												params + ["delta", "dots"],
												params + ["delta", "hets"],
												params + ["delta", "homs"],
												params + ["both", "different alts"]],
									names = param_names + ["origin", "category"])
	
	# leave a zero at the end because percentages are calculated without different_alts value
	perc_data = pd.Series(data=[sample_dots_perc, sample_hets_perc, sample_homs_perc,
					ref_dots_perc, ref_hets_perc, ref_homs_perc,
					delta_dots_perc, delta_hets_perc, delta_homs_perc, 0], index=idx, name="percentage")
	abs_data = pd.Series([sample_dots_abs, sample_hets_abs, sample_homs_abs,
					ref_dots_abs, ref_hets_abs, ref_homs_abs,
					delta_dots_abs, delta_hets_abs, delta_homs_abs, different_alts], index=idx, name="absolute")
	metadata_df = pd.concat((perc_data, abs_data), axis=1)

	return metadata_df


def add_sample_name_to_columns(sample_names, dfs):
	"""make columns into a multiindex to preserve the sample name"""
	for sn, df in zip(sample_names, dfs):
		mux = pd.MultiIndex.from_tuples([(sn, col) for col in df.columns])
		df.columns = mux


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("input_file")
	parser.add_argument("comparison_file")
	parser.add_argument("param_string")
	parser.add_argument("param_names", nargs="*")
	args = parser.parse_args()

	params, param_names = parse_param_string(args.param_string, args.param_names)
	sample_df = parse_vcf(args.input_file)
	ref_df = parse_vcf(args.comparison_file)

	different_alts = remove_differing_alts(sample_df, ref_df)
	sample_state_dfs, ref_state_dfs = build_state_dfs(sample_df, ref_df)
	abs_transition_count_df, percentage_transition_count_df = build_transition_df(sample_df, ref_df)

	abs_confusion_matrices = build_transition_matrices(abs_transition_count_df)
	percentage_confusion_matrices = build_transition_matrices(percentage_transition_count_df)

	metadata = [calc_metadata(pm, am, different_alts, params, param_names) for (pm, am) in zip(percentage_confusion_matrices, abs_confusion_matrices)]

	confusion_vars = [calc_confusion_variables(m, params, param_names) for m in abs_confusion_matrices]
	confusion_vars_missing = [calc_confusion_variables(m, params, param_names, include_missing=True) for m in abs_confusion_matrices]

	# add sample names to confusion vars
	add_sample_name_to_columns(abs_transition_count_df.columns, confusion_vars)
	add_sample_name_to_columns(abs_transition_count_df.columns, confusion_vars_missing)

	with open("input_sample_counts.csv", "w") as f:
		sample_state_dfs[0].to_csv(f)

	with open("sample_position_counts.csv", "w") as f:
		sample_state_dfs[1].to_csv(f)

	with open("ref_sample_counts.csv", "w") as f:
		ref_state_dfs[0].to_csv(f)

	with open("ref_position_counts.csv", "w") as f:
		ref_state_dfs[1].to_csv(f)

	# column names of transition_dfs are sample names
	for i, sample_name in enumerate(abs_transition_count_df.columns):
		with open(f"{sample_name}_metadata.csv", "w") as f:
			metadata[i].to_csv(f)

		with open(f"{sample_name}_matrix.csv", "w") as f:
			percentage_confusion_matrices[i].to_csv(f)

		with open(f"{sample_name}_confusion_vars.csv", "w") as f:
			confusion_vars[i].to_csv(f)

		with open(f"{sample_name}_confusion_vars_missing.csv", "w") as f:
			confusion_vars_missing[i].to_csv(f)


if __name__ == '__main__':
	main()

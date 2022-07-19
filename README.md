# Pipeline for parameteroptimization in high-throughput-data
This pipeline executes parameter-sweeps in variant callers.

## Requirements
- Nextflow
- Apptainer or Singularity
- Conda

## Set-up
Quickly run the pipeline through ```nextflow run <repo-name>```.
If you want to have the code, you can do the following:
1. Install Requirements
2. Clone Repo
3. Run ```nextflow run main.nf```

## Options
Use the *-stub* option to start a minimal run. This is useful for testing and debugging.

## Output
An output directory will be created according to the variant caller config file. In it, there will be a number of 2d-scatterplots, with the f1-score plotted over each parameter. There will also be another file containing a 3d-scatterplot, with the f1-score plotted over the two parameters that have the biggest correlation with it.
In addition, there will be a file containing the top performing parameter combinations.
These results are generated twice each, once disregarding missing-data-calls and once including them.
A directory _report_ will be created. It contains data regarding the ressource usage of the run.

## Adding a custom dataset
Profiles may be used. A dataset defines the following params:
- fasta - fasta-file
- bedFile - bed-file
- bamFile - bam-file or file containing a list of bam-files
- comparison - reference-vcf

## Adding a custom variant caller
Per default, there is an implementation for variant calling using freebayes. It is possible to add additional custom variant callers.
The existing implementation for freebayes may be used as a reference.
The file with the implementation of the variant caller must include a workflow _callVariants_, which will be imported by the main workflow.
This workflow must return a vcf file and the *param-string*.
The param string contains the values of all sweeped parameters, seperated by \_. Their order should correspond to the order of names defined in the config.
The variant caller config defines the following params:
- outdir - the output directory for final results
- workflow - name of the file that defines the callVariants workflow
- the values and ranges for all parameters that need to be tested
- paramNames, which includes the names of all tested parameters, seperated by space
- container - path to the container which contains the variant caller

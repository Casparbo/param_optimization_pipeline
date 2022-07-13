# Pipeline for parameteroptimization in high-throughput-data
With this pipeline you can execute parameter-sweeps in variant callers. The results are presented in scatterplots.

## Requirements
- Nextflow
- Apptainer or Singularity
- Conda

## Quick set-up
You can quickly run the pipeline through ```nextflow run <repo-name>```.
If you want to have the code, you can do the following:
1. Install Requirements
2. Clone Repo
3. Run ```nextflow run main.nf -stub```, all necessary containers and environments will be downloaded. There should be no error messages.
4. Run ```nextflow run main.nf```

## Output
A directory _output_ will be created. In it, you will find a number of 2d-scatterplots, with the f1-score plotted over each parameter. You will also find another file containing a 3d-scatterplot, with the f1-score plotted over the two parameters that have the biggest correlation with it.
These results be there twice each, once disregarding missing-data-calls and once including them.
A directory _report_ will be created. It contains data regarding the ressource usage of the run.

## Adding custom dataset
Profiles may be used. A dataset defines the following params:
- fasta - fasta-file
- bedFile - bed-file
- bamFile - bam-file or file containing a list of bam-files
- comparison - reference-vcf

## Adding a custom variant caller
Per default, there is an implementation for variant calling using freebayes and bcftools. You can add your own implementation of a variant caller of your choice.
The file with the implementation of the variant caller must include a workflow _callVariants_, which will be imported by the main workflow.
The variant caller config defines the following params:
- outdir - the output directory for final results
- workflow - name of the file that defines the callVariants workflow
- the values and ranges for all parameters that need to be tested
- paramNames, which includes the names of all tested parameters, seperated by space
- container - path to the container which contains the variant caller

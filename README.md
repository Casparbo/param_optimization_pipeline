# Pipeline for parameteroptimization in high-throughput-data
## Requirements
- Nextflow
- Apptainer or Singularity
- Conda

## Quickstart
1. Install Requirements
2. Clone Repo
3. Run ```nextflow run main.nf -stub```

## Output
A directory _output_ will be created. It contains the results.
A directory _report_ will be created. It contains data regarding the ressource usage of the run.

## Adding custom datasets
Profiles may be used. A dataset defines the following params:
- fasta - fasta-file
- bedFile - bed-file
- bamFile - bam-file or file containing a list of bam-files
- comparison - reference-vcf

## Adding a custom variant caller
The variant caller must include a workflow _callVariants_, which will be imported by the main workflow.
The variant caller config defines the following params:
- outdir - the output directory for final results
- workflow - name of the file that defines the callVariants workflow
- the values and ranges for all parameters that need to be tested
- paramNames, which includes the names of all tested parameters, seperated by space
- container - path to the container which contains the variant caller

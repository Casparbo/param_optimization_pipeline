# Pipeline for parameteroptimization in high-throughput-data
## Config
Two config files need to be imported: One for the datasets and one for the variant caller.

### Datasets
Profiles may be used. A dataset defines the following params:
- fasta - fasta-file
- bedFile - bed-file
- bamFile - bam-file or file containing a list of bam-files
- comparison - reference-vcf

### Variant caller
The variant caller must include a workflow _callVariants_, which will be imported by the main workflow.
The variant caller config defines the following params:
- outdir - the output directory for final results
- workflow - name of the file that defines the callVariants workflow
- the values and ranges for all parameters that need to be tested
- paramNames, which includes the names of all tested parameters, seperated by space
- container - path to the container which contains the variant caller

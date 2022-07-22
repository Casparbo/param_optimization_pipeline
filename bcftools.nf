#!/usr/bin/env nextflow

process splitBedFile {
  input:
  path bedFile

  output:
  path "split*.bed"

  script:
  """
  lines=\$(wc -l < $bedFile)
  for ((i=1;i<=lines;i+=50))
  do
    sed -n \$i,\$((\$i + 49))p $bedFile > split\$i.bed
  done
  """

  stub:
  """
  sed -n 1,50p $bedFile > split0.bed
  """
}

process bcftools {
  container = params.container

  input:
  path fasta
  path fastaIndex

  each path(bedFile)

  path bamfile
  path bamindex

  each minMQ
  each minBQ

  output:
  tuple val(minMQ), val(minBQ), path("out.vcf")

  script:
  """
  bcftools mpileup \
  -Ou \
  --min-MQ $minMQ \
  --min-BQ $minBQ \
  -f $fasta \
  --regions-file $bedFile \
  $bamfile \
  | bcftools call -c -Ov -o out.vcf -f GQ
  """
}

process sortVCFsamples {
  conda "natsort"

  input:
  tuple val(minMQ), val(minBQ), path(vcf)

  output:
  tuple val(minMQ), val(minBQ), path("sorted.vcf")

  script:
  """
  sort_vcf.py $vcf sorted.vcf
  """
}

process combineVCF {
  container "https://depot.galaxyproject.org/singularity/vcftools:0.1.16--pl5321hd03093a_7"
  input:
  tuple val(minMQ), val(minBQ), path(vcf, stageAs: "*.vcf")

  output:
  publishDir "${params.outdir}/${minMQ}_${minBQ}"
  tuple val("${minMQ}_${minBQ}"), path("combined.vcf")

  script:
  """
  vcf-concat $vcf > combined.vcf
  """
}

process lgcPostProcessing {
  conda "natsort biopython"
  
  input:
  tuple val(minMQ), val(minBQ), path(vcf)

  each hetCorrectFilter
  each minCoverageThresh

  output:
  publishDir "${params.outdir}/${minMQ}_${minBQ}_${hetCorrectFilter}_${minCoverageThresh}", mode:"copy", enabled: workflow.stubRun
  tuple val("${minMQ}_${minBQ}_${hetCorrectFilter}_${minCoverageThresh}"), path("filtered.vcf")

  script:
  """
  VCF_postprocessing.py \
    --input $vcf \
    --output filtered.vcf \
    --genotype-min-count-het-call-threshold $hetCorrectFilter \
    --genotype-min-coverage-threshold $minCoverageThresh
  """
}

workflow callVariants {
	take:
    fasta
    fastaIndex
    bedFile
    bamlist
    bamindex

  main:
    minMQ = Channel.from(params.minMQ)
    minBQ = Channel.from(params.minBQ)

    hetCorrectFilter = Channel.from(params.hetCorrectFilter)
    minCoverageThresh = Channel.from(params.minCoverageThresh)

    if(workflow.stubRun) {
      minMQ = minMQ.first().concat(minMQ.last())
      minBQ = minBQ.first().concat(minBQ.last())

      hetCorrectFilter = hetCorrectFilter.first().concat(hetCorrectFilter.last())
      minCoverageThresh = minCoverageThresh.first().concat(minCoverageThresh.last())
    }

    splitBedFile(bedFile)
    
    bcftools(fasta, fastaIndex, splitBedFile.out.flatMap(), bamlist.collect(), bamindex.collect(),
              minMQ, minBQ)
    sortVCFsamples(bcftools.out)
    sortVCFsamples.out.collect()
    //group by the first 2 elements of tuple (which are the parameters)
    combineVCF(sortVCFsamples.out.groupTuple(by: 0..1))
    //lgcPostProcessing(combineVCF.out, hetCorrectFilter, minCoverageThresh)

  emit:
    combineVCF.out
}

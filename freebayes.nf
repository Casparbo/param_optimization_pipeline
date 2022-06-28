#!/usr/bin/env nextflow

def determine_version() {
  name = params.container.split("/").last()
  name_no_appendix = name.split("--").first()
  version_nr = name_no_appendix.split(":").last()

  return version_nr
}

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
    sed -n \$i,\$((\$i + 50))p $bedFile > split\$i.bed
  done
  """

  stub:
  """
  sed -n 1,50p $bedFile > split0.bed
  """
}

process freebayes {
  container = params.container

  input:
  path fasta
  path fastaIndex

  each path(bedFile)

  path bamfile
  path bamindex

  each minQsum
  each readMismatchLimit
  each minAlternateFraction
  each noMnps
  each noComplex

  output:
  tuple val(minQsum), val(readMismatchLimit), val(minAlternateFraction), val(noMnps), val(noComplex), path("out.vcf")

  script:
  """
  freebayes \
  --report-genotype-likelihood-max \
  --exclude-unobserved-genotypes \
  --genotype-qualities \
  --strict-vcf \
  --ploidy 2 \
  --min-base-quality 20 \
  --min-coverage 4 \
  --min-alternate-count 2 \
  --mismatch-base-quality-threshold 10 \
  --report-monomorphic \
  --min-supporting-allele-qsum $minQsum \
  --read-mismatch-limit $readMismatchLimit \
  --min-alternate-fraction $minAlternateFraction \
  ${if(noMnps) {
      return'--no-mnps'
    } else {
      return ''
    }
  } \
  ${if(noComplex) {
      return'--no-complex'
    } else {
      return ''
    }
  } \
  --fasta-reference $fasta \
  --target $bedFile \
  $bamfile > out.vcf
  """
}

process sortVCFsamples {
  conda "natsort pandas"

  input:
  tuple val(minQsum), val(readMismatchLimit), val(minAlternateFraction), val(noMnps), val(noComplex), path(vcf)

  output:
  tuple val(minQsum), val(readMismatchLimit), val(minAlternateFraction), val(noMnps), val(noComplex), path("sorted.vcf")

  script:
  """
  sort_vcf.py $vcf sorted.vcf
  """
}

process combineVCF {
  container "https://depot.galaxyproject.org/singularity/vcftools:0.1.16--pl5321hd03093a_7"
  input:
  tuple val(minQsum), val(readMismatchLimit), val(minAlternateFraction), val(noMnps), val(noComplex), path(vcf, stageAs: "*.vcf")

  output:
  //publishDir "${params.outdir}/FREEBAYES_${determine_version()}_${minQsum}_${readMismatchLimit}_${minAlternateFraction}_${noMnps}_${noComplex}", mode:"copy"
  tuple val(minQsum), val(readMismatchLimit), val(minAlternateFraction), val(noMnps), val(noComplex), path("combined.vcf")

  script:
  """
  vcf-concat $vcf > combined.vcf
  """
}

process lgcPostProcessing {
  conda "natsort biopython"
  
  input:
  tuple val(minQsum), val(readMismatchLimit), val(minAlternateFraction), val(noMnps), val(noComplex), path(vcf)

  each hetCorrectFilter
  each minCoverageThresh

  output:
  publishDir "${params.outdir}/FREEBAYES_${determine_version()}_${minQsum}_${readMismatchLimit}_${minAlternateFraction}_${noMnps}_${noComplex}_${hetCorrectFilter}_${minCoverageThresh}", mode:"copy"
  tuple val("FREEBAYES_${determine_version()}_${minQsum}_${readMismatchLimit}_${minAlternateFraction}_${noMnps}_${noComplex}_${hetCorrectFilter}_${minCoverageThresh}"), path("filtered.vcf")

  script:
  """
  /illumina/software/Development/Pipelines-SA-dev_Docker-py3/VCF_postprocessing.py \
    --input $vcf \
    --output filtered.vcf \
    --genotype-min-coverage-threshold $minCoverageThresh \
    --genotype-min-count-het-call-threshold $hetCorrectFilter
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
    minQsum = Channel.from(params.minQsum)
    readMismatchLimit = Channel.from(params.readMismatchLimit)
    minAlternateFraction = Channel.from(params.minAlternateFraction)
    noMnps = Channel.from(params.noMnps)
    noComplex = Channel.from(params.noComplex)

    hetCorrectFilter = Channel.from(params.hetCorrectFilter)
    minCoverageThresh = Channel.from(params.minCoverageThresh)

    if(workflow.stubRun) {
      minQsum = minQsum.first().mix(minQsum.last())
      readMismatchLimit = readMismatchLimit.first().mix(readMismatchLimit.last())
      minAlternateFraction = minAlternateFraction.first().mix(minAlternateFraction.last())
      noMnps = noMnps.first().mix(noMnps.last())
      noComplex = noComplex.first().mix(noComplex.last())

      hetCorrectFilter = hetCorrectFilter.first().mix(hetCorrectFilter.last())
      moinCoverageThresh = minCoverageThresh.first().mix(minCoverageThresh.last())
    }

    splitBedFile(bedFile)
    
    freebayes(fasta, fastaIndex, splitBedFile.out.flatMap(), bamlist.collect(), bamindex.collect(),
    minQsum, readMismatchLimit, minAlternateFraction, noMnps, noComplex)
    sortVCFsamples(freebayes.out)
    sortVCFsamples.out.collect()
    //group by the first 5 elements of tuple (which are the parameters)
    combineVCF(sortVCFsamples.out.groupTuple(by: 0..4))
    lgcPostProcessing(combineVCF.out, hetCorrectFilter, minCoverageThresh)

  emit:
    lgcPostProcessing.out
}
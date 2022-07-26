#!/usr/bin/env nextflow

def createBamChannels(bamfile) {
  //check if input ist bamfile or bamlist, then add files to Channel
  bamlist = Channel.empty()
  bamindex = Channel.empty()

  fileExtension = bamfile.name.substring(bamfile.name.length()-3)
  
  if(fileExtension.equals("bam")) {
    bamlist << file(bamfile)
    bamindex << file(bamfile + ".{bai, csi}")
  } else {
    bamfile.eachLine {line ->
      if (line) {
        bamlist << file(line)
        bamindex << file(line + ".{bai, csi}")
      }
    }
  }

  return [bamlist, bamindex]
}

process txtToVcf {
  conda "pandas"
  
  input:
  path txt

  output:
  path "converted.vcf"

  script:
  """
  txt_to_vcf.py $txt
  """
}

process vcfPandas {
  //debug true
  conda "pandas seaborn"
  
  input:
  tuple val(paramString), path(vcf), path(comparison), val(paramNames)

  output:
  publishDir "${params.outdir}/$paramString", mode:"copy", enabled: workflow.stubRun
  
  path "*_confusion_vars.csv", emit: confusionVars
  path "*_confusion_vars_missing.csv", emit: confusionVarsMissing
  // also publish all other generated csv
  path "*.csv"

  script:
  """
  vcf_analysis.py $vcf $comparison $paramString $paramNames
  """
}

// to avoid a too long input for next process, collect files into directory and then output that
// preserve idx for next process
process collectOutputToDirectory {
  input:
  tuple path(input_files, stageAs: "*.csv"), val(idx)

  output:
  tuple path("filedir"), val(idx)

  script:
  """
  mkdir filedir
  ln $input_files filedir
  """
}

process combineAnalyses {
  debug true
  conda "pandas"
  beforeScript "ulimit -Ss unlimited"
  input:
  tuple path(fileList, stageAs: "*filedir"), val(idx)

  output:  
  tuple path("*_combined.csv"), val(idx)

  script:
  """
  combine_analyses.py $fileList
  """
}

process metaAnalysis {
  debug true
  conda "pandas matplotlib seaborn"
  input:
  tuple path(confusionVarsNoMissing, stageAs: "*noMissing.csv"), path(confusionVarsMissing, stageAs: "*missing.csv")

  output:
  publishDir "${params.outdir}", mode: "copy"
  path "metadata_*.png"
  path "top_params_*.png"

  script:
  """
  meta_analysis.py --confusion_vars $confusionVarsNoMissing --confusion_vars_missing $confusionVarsMissing --top_params_thresh 20
  """
}

include {callVariants} from "./${params.workflow}"

workflow {
  fasta = Channel.fromPath(params.fasta)
  fastaIndex = Channel.fromPath(params.fasta + ".fai")
  bedFile = Channel.fromPath(params.bedFile)

  (bamlist, bamindex) = createBamChannels(file(params.bamfile))

  comparison = Channel.fromPath(params.comparison)
  // if comparison is a txt file, convert it to vcf and move on
  if(params.comparison[-3..-1] == "txt") {
    comparison = txtToVcf(comparison)
  }

  paramNames = Channel.from(params.paramNames)

  callVariants(fasta, fastaIndex, bedFile, bamlist, bamindex)

  vcfPandas(callVariants.out.combine(comparison).combine(paramNames))

  // add number so that we can later tell which one is with and without missing data (0: without, 1: with)
  confusionBuckets = vcfPandas.out.confusionVars.flatten().collate(50).map{[it, 0]}
  confusionBucketsMissing = vcfPandas.out.confusionVarsMissing.flatten().collate(50).map{[it, 1]}
  confusionBucketsBoth = confusionBuckets.concat(confusionBucketsMissing)

  collectOutputToDirectory(confusionBucketsBoth)

  combineAnalyses(collectOutputToDirectory.out.groupTuple(by: 1))

  // group by added number (0/1), then fork and recombine the two paths so that the one without missing data comes first
  bothCombined = combineAnalyses.out.groupTuple(by: 1).branch{no_missing: it[1] == 0; missing: it[1] == 1}
  bothPaths = bothCombined.no_missing.concat(bothCombined.missing).collect{it[0]}

  metaAnalysis(bothPaths)
}

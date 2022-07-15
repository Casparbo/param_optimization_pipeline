#!/usr/bin/env nextflow

def createBamChannels(bamfile) {
  //check if input ist bamfile or bamlist, then add files to Channel
  bamlist = Channel.empty()
  bamindex = Channel.empty()

  fileExtension = bamfile.name.substring(bamfile.name.length()-3)
  
  if(fileExtension.equals("bam")) {
    bamlist << file(bamfile)
    bamindex << file(bamfile + ".csi")
  } else {
    bamfile.eachLine {line ->
      if (line) {
        bamlist << file(line)
        bamindex << file(line + ".csi")
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
  
  path "metadata.csv", emit: metadata
  path "confusion_vars.csv", emit: confusionVars
  path "confusion_vars_missing.csv", emit: confusionVarsMissing
  path "matrix.csv"
  path "heatmap.png"
  path "input_sample_counts.csv"
  path "sample_position_counts.csv"
  path "ref_sample_counts.csv"
  path "ref_position_counts.csv"

  script:
  """
  vcf_analysis.py $vcf $comparison $paramString $paramNames
  """
}

process combineAnalyses {
  conda "pandas"
  beforeScript "ulimit -Ss unlimited"
  input:
  tuple path(fileList, stageAs: "*.csv"), val(idx)

  output:
  tuple path("combined.vcf"), val(idx)

  script:
  """
  combine_analyses.py $fileList
  """
}

process metaAnalysis {
  debug true
  conda "pandas matplotlib seaborn"
  input:
  path(combinedConfusionVars, stageAs: "*confusionVars.csv")

  output:
  publishDir "${params.outdir}", mode: "copy"
  path "metadata.pdf"
  path "metadata_3d.pdf"

  script:
  """
  meta_analysis.py $combinedConfusionVars
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
  confusionVarsFiles = vcfPandas.out.confusionVars.collect().concat(Channel.from(0)).toList()
  confusionVarsMissingFiles = vcfPandas.out.confusionVarsMissing.collect().concat(Channel.from(1)).toList()
  confusionVarsFilesBoth = confusionVarsFiles.concat(confusionVarsMissingFiles)

  combineAnalyses(confusionVarsFilesBoth)

  // fork and recombine the two paths so that the one without missing data comes first
  bothCombined = combineAnalyses.out.branch{no_missing: it[1] == 0; missing: it[1] == 1}
  bothPaths = bothCombined.no_missing.concat(bothCombined.missing).collect{it[0]}

  metaAnalysis(bothPaths)
}

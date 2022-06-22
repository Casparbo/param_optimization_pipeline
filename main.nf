#!/usr/bin/env nextflow

def bedFileToRegionStrings(bedFile) {
  regionList = Channel.empty()

  if(workflow.stubRun) {
    splt = bedFile.readLines()[0].split("\t")
    regionString = "${splt[0]}:${splt[1]}..${splt[2]}"
    regionList << regionString

    return regionList
  }

  bedFile.eachLine {str ->
    splt = str.split("\t")
    regionString = "${splt[0]}:${splt[1]}..${splt[2]}"
    regionList << regionString
  }

  return regionList
}

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
  debug true
  conda "pandas seaborn"
  
  input:
  tuple val(paramString), path(vcf), path(comparison)

  output:
  publishDir "${params.outdir}/$paramString", mode:"copy"
  path "metadata.csv", emit: metadata
  path "confusion_vars.csv", emit: confusionVars
  path "matrix.csv"
  path "heatmap.png"
  path "input_sample_counts.csv"
  path "sample_position_counts.csv"
  path "ref_sample_counts.csv"
  path "ref_position_counts.csv"

  script:
  """
  vcf_analysis.py $vcf $comparison $paramString
  """
}

process metaAnalysis {
  debug true
  conda "pandas matplotlib seaborn"
  input:
  path(metadataList, stageAs: "*metadata.csv")

  output:
  publishDir "${params.outdir}", mode: "copy"
  path "metadata_*.png"

  script:
  """
  meta_analysis.py $metadataList
  """
}

include {callVariants} from "./${params.workflow}"

workflow {
  fasta = Channel.fromPath(params.fasta)
  fastaIndex = Channel.fromPath(params.fasta + ".fai")

  bedFile = file(params.bedFile)

  targets = bedFileToRegionStrings(bedFile)

  (bamlist, bamindex) = createBamChannels(file(params.bamfile))

  comparison = Channel.fromPath(params.comparison)
  // if comparison is a txt file, convert it to vcf and move on
  if(params.comparison[-3..-1] == "txt") {
    comparison = txtToVcf(comparison)
  }

  callVariants(fasta, fastaIndex, targets, bamlist, bamindex)

  vcfPandas(callVariants.out.combine(comparison))

  metaAnalysis(vcfPandas.out.confusionVars.collect())
}

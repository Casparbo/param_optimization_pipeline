#!/usr/bin/env python
import sys 
import os
import subprocess
import tempfile
import natsort
import copy
import shutil
from Bio import SeqIO
# Py2/3 compatibility
try:
    basestring
except NameError:
    basestring = str

def _check_file_object(filename,mode="r"):
    if filename is None:
        if mode=="r":
            return sys.stdin, False
        elif mode=="w":
            return sys.stdout, False
    elif isinstance(filename, basestring):
        return open(filename,mode), True
    else:
        # assume file object
        return filename, False

def _calculate_AO(format_data,sample_fields):
    AO = {}
    for sample_field in sample_fields:
        sample_data = sample_field.split(":")
        if (len(sample_data) == len(format_data)): # otherwise it is truncated
            # count AO
            if sample_data[format_data.index("AO")] != ".":
                all_counts = sample_data[format_data.index("AO")].split(",")
                for i in range(0,len(all_counts)):
                    if not i in AO:
                        AO[i] = 0
                    AO[i] += int(all_counts[i])
    return ",".join([str(AO[allele_index]) for allele_index in sorted(list(AO.keys()))])

def remove_uncalled_alt_alleles(input_vcf,output_vcf,logger,sample_string):
    
    # removes ALTs and associated fields from the VCF
    infile,close_input_on_exit = _check_file_object(filename=input_vcf,mode="r")
    outfile,close_output_on_exit = _check_file_object(filename=output_vcf,mode="w")
    current_line = infile.readline().strip()
    
    # parse input
    sample_fields_to_update_total_indexing = []
    sample_fields_to_update_alt_indexing = []
    info_fields_to_update_total_indexing = []
    info_fields_to_update_alt_indexing = []
    line_count = 0
    multiallelic_line_count = 0
    updated_line_count = 0
    while current_line != "":
        if current_line[0] == "#":
            line_parts = current_line.replace("#","").split("=",1)
            # parse INFO and FORMAT annotation fields
            if line_parts[0] in ["FORMAT","INFO"]:
                try:
                    line_data = {}
                    # parse into dict
                    for field in line_parts[1][1:-2].split(","):
                        field_data = field.split("=")
                        line_data[field_data[0]] = field_data[1]
                        if field_data[0] == "Number":
                            break
                except:
                    pass # ignore parsing errors for now
                
            # Check which fields contain allele data, and mark them for updates
            if line_parts[0] == "INFO":
                if line_data["Number"] == "A":
                    info_fields_to_update_alt_indexing.append(line_data["ID"])
                elif line_data["Number"] == "R":
                    info_fields_to_update_total_indexing.append(line_data["ID"])
            elif line_parts[0] == "FORMAT":
                if line_data["Number"] == "A":
                    sample_fields_to_update_alt_indexing.append(line_data["ID"])
                elif line_data["Number"] == "R":
                    sample_fields_to_update_total_indexing.append(line_data["ID"])
            
        else:
            # parse the line
            line_parts = current_line.split("\t")
            alt_alleles = line_parts[4].split(",")
            format_data = line_parts[8].split(":")
            allele_sequences = [line_parts[3]] + alt_alleles
            try:
                alt_alleles.remove(".")
            except:
                pass
            if len(alt_alleles) > 1:
                multiallelic_line_count += 1
            # Check which alt alleles are called
            sample_fields = line_parts[9:]
            called_alleles = set(['0']) # always include the ref allele
            for sample_field in sample_fields:
                called_alleles.update(set(sample_field.split(":")[0].split("/")))
            try:
                called_alleles.remove(".")
            except:
                pass
            
            # Are all alleles called?
            if len(called_alleles) != len(alt_alleles) + 1:
                # no, let's remove some
                all_alleles_to_keep = sorted([int(x) for x in called_alleles])
                alt_alleles_to_keep = sorted([int(x)-1 for x in called_alleles if x != '0'])
                allele_sequences_to_keep = [allele_sequences[index] for index in all_alleles_to_keep]
                
                # update ALT field
                new_alt_alleles = [alt_alleles[index] for index in alt_alleles_to_keep]
                if len(new_alt_alleles) == 0:
                    new_alt_alleles = ["."]
                line_parts[4] = ",".join(new_alt_alleles)
                
                # Update sample fields
                for sample_index in range(0,len(sample_fields)):
                    sample_data = sample_fields[sample_index].split(":")
                    if (len(sample_data) == len(format_data)): # otherwise it is truncated
                        # recode GT indices
                        if sample_data[format_data.index("GT")] != ".":
                            sample_data[format_data.index("GT")] = "/".join([str(allele_sequences_to_keep.index(allele_sequences[int(index)])) for index in sample_data[format_data.index("GT")].split("/")])
                        
                        # update fields with alt indexing
                        for sample_field in sample_fields_to_update_alt_indexing:
                            if len(alt_alleles_to_keep) > 0:
                                sample_data[format_data.index(sample_field)] = ",".join([sample_data[format_data.index(sample_field)].split(",")[index] for index in alt_alleles_to_keep])
                            else:
                                sample_data[format_data.index(sample_field)] = "."
                        # update fields with all allele indexing
                        for sample_field in sample_fields_to_update_total_indexing:
                            if len(all_alleles_to_keep) > 0:
                                sample_data[format_data.index(sample_field)] = ",".join([sample_data[format_data.index(sample_field)].split(",")[index] for index in all_alleles_to_keep])
                            else:
                                sample_data[format_data.index(sample_field)] = "."
                    
                    # update sample data
                    sample_fields[sample_index] = ":".join(sample_data)
                
                # update line
                line_parts = line_parts[0:9] + sample_fields
                
                # update INFO fields
                info_data = line_parts[7].split(";")
                for info_index in range(0,len(info_data)):
                    info_field_data = info_data[info_index].split("=")
                    if (len(info_field_data) > 0): # otherwise it is truncated
                        if (info_field_data[1] != "."):
                            updated = False
                            if info_field_data[0] == "AO":
                                info_field_data[1] = _calculate_AO(format_data,sample_fields)
                            elif info_field_data[0] == "MinAF":
                                info_field_data[1] = info_field_data[1].split(",")[0]
                            elif info_field_data[0] in info_fields_to_update_alt_indexing:
                                updated = True
                                if len(alt_alleles_to_keep) > 0:
                                    info_field_data[1] = ",".join([info_field_data[1].split(",")[index] for index in alt_alleles_to_keep])
                                else:
                                    info_field_data[1] = "."
                            elif info_field_data[0] in info_fields_to_update_total_indexing:
                                updated = True
                                if len(alt_alleles_to_keep) > 0:
                                    info_field_data[1] = ",".join([info_field_data[1].split(",")[index] for index in all_alleles_to_keep])
                                else:
                                    info_field_data[1] = "."
                            
                            # update field
                            info_data[info_index] = "=".join(info_field_data)
                # update into line
                line_parts[7] = ";".join(info_data)
                
                # update counter
                updated_line_count += 1
            
            # Re-construct output
            current_line = "\t".join(line_parts)
            line_count += 1
        
        # Write output
        outfile.write(current_line+"\n")
        
        # Next line
        current_line = infile.readline().rstrip()
    
    logger.debug(str(line_count) + " variants (of which " + str(multiallelic_line_count) + " multiallelic); " + str(updated_line_count) + " updated")
    
    # Clean up
    if close_output_on_exit:
        outfile.close()
        del outfile
    if close_input_on_exit:
        infile.close()
        del infile
    
    # Done
    return

def passthrough_filter_VCF(input_vcf,output_vcf,passthrough_filters):
    # filters VCF for anything but passthrough filters
    infile,close_input_on_exit = _check_file_object(filename=input_vcf,mode="r")
    outfile,close_output_on_exit = _check_file_object(filename=output_vcf,mode="w")
    current_line = infile.readline().strip()
    passthrough_filters.append("PASS")
    
    # parse input
    while current_line != "":
        keep = True
        if current_line[0] == "#":
            # Propagate info
            output_line = current_line
        else:
            # Parse the VCF data line and check filter
            line_parts = current_line.split("\t")
            filter_field = line_parts[6]
            
            filter_hits = 0
            for filter_name in filter_field.split(";"):
                if filter_name == ".":
                    line_parts[6] = "PASS"
                    current_line = "\t".join(line_parts)
                elif not filter_name in passthrough_filters:
                    filter_hits += 1
            
            # Construct filter field
            if filter_hits > 0:
                keep = False
        
        # Write output
        if keep:
            outfile.write(current_line+"\n")
        
        # Next line
        current_line = infile.readline().rstrip()
    
    # Clean up
    if close_output_on_exit:
        outfile.close()
        del outfile
    if close_input_on_exit:
        infile.close()
        del infile
    
    # Done
    return


def apply_padding_correction_to_VCF(input_vcf,output_vcf,locus_padding_length):
    # filters VCF for anything but passthrough filters
    infile,close_input_on_exit = _check_file_object(filename=input_vcf,mode="r")
    outfile,close_output_on_exit = _check_file_object(filename=output_vcf,mode="w")
    current_line = infile.readline().rstrip()
    
    # parse input
    while current_line != "":
        if current_line[0] == "#":
            # Propagate info
            output_line = current_line
        else:
            # Parse the VCF data and apply position correction
            line_parts = current_line.split("\t")
            output_fields = [line_parts[0]]
            output_fields.append(str(int(line_parts[1])-locus_padding_length))
            output_fields.extend(line_parts[2:])
            output_line = "\t".join(output_fields)
            
        outfile.write(output_line+"\n")
        
        # Next line
        current_line = infile.readline().rstrip()
    
    # Clean up
    if close_output_on_exit:
        outfile.close()
        del outfile
    if close_input_on_exit:
        infile.close()
        del infile
    
    # Done
    return

def Freebayes_VCF_background_variant_filter(input_vcf,output_vcf,background_samples,logger,sample_string):
    # inspects background sample genotypes, filters variants if any allele thereof if observed in any of the background samples
    infile,close_input_on_exit = _check_file_object(filename=input_vcf,mode="r")
    outfile,close_output_on_exit = _check_file_object(filename=output_vcf,mode="w")
    current_line = infile.readline().strip()
    
    # statistic counters
    parsed_variant_count = 0
    passed_variant_count = 0
    missing_background_count = 0
    inconsistent_background_count = 0
    background_variant_count = 0
    background_variant_only_count = 0
    
    # parse input
    while current_line != "":
        if current_line[0] == "#":
            # Propagate info
            output_line = current_line
            if current_line.startswith("#CHROM"):
                # Get samples
                line_parts = current_line.split("\t")
                samples = line_parts[9:]
                foreground_samples = copy.deepcopy(samples)
                # check if all background samples are there
                for background_sample in background_samples:
                    foreground_samples.remove(background_sample)
                    if background_sample not in samples:
                        logger.warning("Freebayes_VCF_background_variant_filter: Background {} cannot be found in VCF file".format(background_sample))
                
                # Add filter descriptions
                original_output_line=output_line
                output_line = "##FILTER=<ID=MissingBackground,Description=\"All replicates for background samples have missing data\">\n"
                output_line += "##FILTER=<ID=InconsistentBackground,Description=\"Mismatching genotypes found between replicates of background samples\">\n"
                output_line += "##FILTER=<ID=BackgroundVariant,Description=\"ALT allele present in at least one background sample\">\n"
                output_line += "##FILTER=<ID=BackgroundVariantOnly,Description=\"ALT allele only present in background samples\">\n"
                output_line+=original_output_line
                # => PASS all if the allele isn't observed in any of the background samples, and at least one background sample hold data
        else:
            # Parse the VCF data line and check filter
            line_parts = current_line.split("\t")
            variant_id = line_parts[0]+"_"+line_parts[1]
            format_parts = line_parts[8].split(":")
            sample_data = line_parts[9:]
            genotype_field_index = format_parts.index("GT")
            filter_field = line_parts[6]
            info_field = line_parts[7]
            alt_field = line_parts[5]
            
            # Checking consistency of parent GT calls
            background_gts=[sample_data[samples.index(background_sample)].split(":")[genotype_field_index] for background_sample in background_samples]
            foreground_gts=[sample_data[samples.index(foreground_sample)].split(":")[genotype_field_index] for foreground_sample in foreground_samples]
            present_background_gts = list(set(background_gts))
            present_foreground_gts = list(set(foreground_gts))
            if ("." in present_background_gts):
                present_background_gts.remove(".")
                background_replicates_complete=False # there are missing genotypes among background replicates
            else:
                background_replicates_complete=True
            if ("." in present_foreground_gts):
                present_foreground_gts.remove(".")

            # Extract alleles
            alt_alleles = [str(x) for x in range(1,len(alt_field.split(","))+1)]
            present_background_alleles = list(set([item for sublist in [genotype.split("/") for genotype in present_background_gts] for item in sublist]))
            present_foreground_alleles = list(set([item for sublist in [genotype.split("/") for genotype in present_foreground_gts] for item in sublist]))
            
            # Check for available genotypes for this parent
            if len(present_background_gts) == 0:
                background_missing=True
            else:
                background_missing=False
            
            # Check for consistency
            if len(present_background_gts) > 1:
                background_consistent=False
            else:
                background_consistent=True
            
            # Check if alt allele(s) exist(s) in background samples            
            background_variant=False
            for alt_allele in alt_alleles:
                if alt_allele in present_background_alleles:
                    background_variant=True
                    break
            
            # Check if alt allele exists in background samples only
            background_variant_only=True
            for alt_allele in alt_alleles:
                if alt_allele in present_foreground_alleles:
                    background_variant_only=False
                    break
            
            # Apply criteria
            filter_hits = []
            if filter_field != "PASS":
                filter_hits=filter_field.split(";")
            
            if background_missing:
                filter_hits.append("MissingBackground")
                missing_background_count += 1
            if not background_consistent:
                filter_hits.append("InconsistentBackground")
                inconsistent_background_count += 1
            if background_variant:
                filter_hits.append("BackgroundVariant")
                background_variant_count += 1
            if background_variant_only:
                filter_hits.append("BackgroundVariantOnly")
                background_variant_only_count += 1
            
            # Construct filter field
            if len(filter_hits) == 0:
                # Passed!
                passed_variant_count+=1
                filter_field = "PASS"
            else:
                filter_field = ";".join(filter_hits)
            
            # Construct output line
            output_line = ("\t".join(line_parts[0:6]))+"\t"+filter_field+"\t"+info_field+"\t"+("\t".join(line_parts[8:]))
        
        # Stats
        parsed_variant_count += 1
        if (parsed_variant_count % 1000 == 0):
            logger.info(sample_string+" - Freebayes_VCF_background_variant_filter() RUNNING - scanned %d variants; PASS: %d, completely missing background GTs: %d, inconsistent background GTs: %d, variants present in background samples: %d, variants only present in background samples: %d" % (parsed_variant_count,passed_variant_count,missing_background_count,inconsistent_background_count,background_variant_count,background_variant_only_count))
        
        # Write output
        outfile.write(output_line+"\n")
        
        # Next line
        current_line = infile.readline().strip()
    
    # Done
    logger.info(sample_string+" - Freebayes_VCF_background_variant_filter() DONE - scanned %d variants; PASS: %d, completely missing background GTs: %d, inconsistent background GTs: %d, variants present in background samples: %d, variants only present in background samples: %d" % (parsed_variant_count,passed_variant_count,missing_background_count,inconsistent_background_count,background_variant_count,background_variant_only_count))
    
    # Clean up
    if close_output_on_exit:
        outfile.close()
        del outfile
    if close_input_on_exit:
        infile.close()
        del infile
    
    # Done
    return

def Freebayes_VCF_parental_allele_filter(info_vcf_filename,output_vcf_filename,parents,logger,sample_string,generation=1):
    # inspects parental genotypes, filters missing or inconsistent data
    infile,close_input_on_exit = _check_file_object(filename=info_vcf_filename,mode="r")
    outfile,close_output_on_exit = _check_file_object(filename=output_vcf_filename,mode="w")
    current_line = infile.readline().strip()
    
    parsed_variant_count = 0
    passed_variant_count = 0
    missing_parents_count = 0
    missing_parental_replicate_count = 0
    inconsistent_parents_count = 0
    non_parental_alt_allele_count = 0
    parents_monomorphic_count = 0
    parents_identical_count = 0
    
    while current_line != "":
        if current_line[0] == "#":
            # Propagate info
            output_line = current_line
            if current_line.startswith("#CHROM"):
                # Get samples
                line_parts = current_line.split("\t")
                samples = line_parts[9:]
                # check if all parents are there
                for parent in parents[0]+parents[1]:
                    if parent not in samples:
                        logger.warning("Freebayes_VCF_parental_allele_filter: Parent {} cannot be found in VCF file".format(parent))
                
                # Add filter descriptions
                original_output_line=output_line
                output_line = "##FILTER=<ID=MissingParentalReplicate,Description=\"Parental replicate missing\">\n"
                output_line += "##FILTER=<ID=MissingParent,Description=\"All replicates for a parent missing\">\n"
                output_line += "##FILTER=<ID=InconsistentParents,Description=\"Mismatching genotypes found between replicates of parental samples\">\n"
                output_line += "##FILTER=<ID=NonParentalAltAllele,Description=\"ALT allele called in any sample without being present in parental alleles\">\n"
                output_line += "##FILTER=<ID=MonomorphicParents,Description=\"All parents carry only one allele\">\n"
                output_line += "##FILTER=<ID=IdenticalParents,Description=\"All parents have the same genotype\">\n"
                output_line += "##FILTER=<ID=HomozygousParents,Description=\"All parents are homozygous, so the alleles will not segregate in an F1 population\">\n"
                
                output_line+=original_output_line
        else:
            # Parse the VCF data line and check filter
            line_parts = current_line.split("\t")
            variant_id = line_parts[0]+"_"+line_parts[1]
            format_parts = line_parts[8].split(":")
            sample_data = line_parts[9:]
            genotype_field_index = format_parts.index("GT")
            filter_field = line_parts[6]
            info_field = line_parts[7]
            alt_field = line_parts[5]
            
            # Checking consistency of parent GT calls
            parent_replicates_complete=True
            parents_complete=True
            parent_replicates_consistent=True
            parents_monomorphic=False
            parents_identical=False
            parent_het_count=0
            all_parent_gts=[]
            for parent_replicates in parents:
                parent_gts=[sample_data[samples.index(parent_rep)].split(":")[genotype_field_index] for parent_rep in parent_replicates]
                
                for gt in parent_gts:
                    if len(list(set(gt.split("/")))) > 1:
                        parent_het_count+=1
                
                all_parent_gts+=parent_gts
                present_parent_gts = list(set(parent_gts))
                if ("." in parent_gts):
                    if (len(parent_gts) > 1):
                        parent_replicates_complete=False
                    present_parent_gts.remove(".")
                # Check for available genotypes for this parent
                if len(present_parent_gts) == 0:
                    parents_complete=False
                # Check for consistency
                if len(present_parent_gts) > 1:
                    parent_replicates_consistent=False
            
            # check if all parental calls are identical
            present_all_parent_gts = list(set(all_parent_gts))
            if "." in present_all_parent_gts:
                present_all_parent_gts.remove(".")
            all_parent_alleles=list(set([allele for alleles in present_all_parent_gts for allele in alleles.split("/")]))
            if len(set(all_parent_alleles))==1:
                parents_monomorphic=True
            if len(set(present_all_parent_gts))==1:
                parents_identical=True
            
            # Apply criteria
            if not parents_complete:
                if filter_field == "PASS":
                    filter_field = ""
                else:
                    filter_field += ";"
                filter_field += "MissingParent"
                missing_parents_count+=1
            
            if not parent_replicates_complete:
                if filter_field == "PASS":
                    filter_field = ""
                else:
                    filter_field += ";"
                filter_field += "MissingParentalReplicate"
                missing_parental_replicate_count+=1
            
            if not parent_replicates_consistent:
                if filter_field == "PASS":
                    filter_field = ""
                else:
                    filter_field += ";"
                filter_field += "InconsistentParents"
                inconsistent_parents_count+=1
            
            # check if ALT alleles match parental alleles
            alt_alleles = [str(x) for x in range(1,len(alt_field.split(","))+1)]
            for alt_allele in alt_alleles:
                if alt_allele not in all_parent_alleles:
                    # Apply criteria
                    if filter_field == "PASS":
                        filter_field = ""
                    else:
                        filter_field += ";"
                    filter_field += "NonParentalAltAllele"
                    non_parental_alt_allele_count+=1
                break
            
            if parents_complete:
                if (generation == 1) and ((not parents_monomorphic) and parents_identical):
                    # het parents will segregate in F1
                    pass
                else:
                    # Apply monomorphic parent criterium
                    if parents_monomorphic:
                        if filter_field == "PASS":
                            filter_field = ""
                        else:
                            filter_field += ";"
                        filter_field += "MonomorphicParents"
                        parents_monomorphic_count+=1
                    
                    # Apply identical parent criterium
                    if parents_identical:
                        if filter_field == "PASS":
                            filter_field = ""
                        else:
                            filter_field += ";"
                        filter_field += "IdenticalParents"
                        parents_identical_count+=1
                    else:
                        if (generation == 1) and (parent_het_count == 0):
                            if filter_field == "PASS":
                                filter_field = ""
                            else:
                                filter_field += ";"
                            filter_field += "HomozygousParents"
            
            # Passed?
            if filter_field == "PASS":
                passed_variant_count+=1
            
            # Construct output line
            output_line = ("\t".join(line_parts[0:6]))+"\t"+filter_field+"\t"+info_field+"\t"+("\t".join(line_parts[8:]))
        
        # Stats
        parsed_variant_count += 1
        if (parsed_variant_count % 1000 == 0):
            logger.info(sample_string+" - Freebayes_VCF_parental_allele_filter() RUNNING - scanned %d variants; PASS: %d, parents with missing GTs: %d, parents with missing replicate GT data: %d, parents with inconsistent GTs: %d, ALT alleles not matching parent GTs: %d, monomorphic parents: %d, identical parents: %d" % (parsed_variant_count,passed_variant_count,missing_parents_count,missing_parental_replicate_count,inconsistent_parents_count,non_parental_alt_allele_count,parents_monomorphic_count,parents_identical_count))
        
        # Write output
        outfile.write(output_line+"\n")
        
        # Next line
        current_line = infile.readline().strip()
    
    # Done
    logger.info(sample_string+" - Freebayes_VCF_parental_allele_filter() DONE - scanned %d variants; PASS: %d, parents with missing GTs: %d, parents with missing replicate GT data: %d, parents with inconsistent GTs: %d, ALT alleles not matching parent GTs: %d, monomorphic parents: %d, identical parents: %d" % (parsed_variant_count,passed_variant_count,missing_parents_count,missing_parental_replicate_count,inconsistent_parents_count,non_parental_alt_allele_count,parents_monomorphic_count,parents_identical_count))
    
    # Clean up
    if close_output_on_exit:
        outfile.close()
        del outfile
    if close_input_on_exit:
        infile.close()
        del infile
    # Done
    return

def append_BED_IDs_to_VCF(input_vcf,output_vcf,target_bed):
    # IDs?
    variant_calling_bed_id_filename = target_bed.rsplit(".")[0]+"_ids.txt"
    if not os.path.exists(variant_calling_bed_id_filename):
        variant_calling_bed_id_file = open(variant_calling_bed_id_filename,"w")
        for line_fields in [x.rstrip().split("\t") for x in open(target_bed).readlines()]:
            if len(line_fields) > 2:
                # use existing IDs
                variant_calling_bed_id_file.write("%s_%d\t%s\n" % (line_fields[0],int(line_fields[1])+1,line_fields[3]))
        variant_calling_bed_id_file.close()
    
    # Check if variant calling BED includes IDs
    bed_identifiers = dict([line_fields for line_fields in [current_line.rstrip().split("\t") for current_line in open(variant_calling_bed_id_filename).readlines()]])
    
    # VCF
    outfile,close_output_on_exit = _check_file_object(filename=output_vcf,mode="w")
    infile,close_input_on_exit  = _check_file_object(filename=input_vcf,mode="r")
    
    current_line = infile.readline()
    while current_line != "":
        if current_line[0] == "#":
            output_line = current_line.rstrip()
        else:
            line_fields = current_line.rstrip().split("\t")
            if ((line_fields[0]+"_"+line_fields[1]) in bed_identifiers.keys()):
                line_fields[2] = bed_identifiers[line_fields[0]+"_"+line_fields[1]]
            output_line = "\t".join(line_fields)
            
        # write output
        outfile.write(output_line+"\n")
        
        # next
        current_line = infile.readline()
    
    # done
    if close_output_on_exit:
        outfile.close()
        del outfile
    if close_input_on_exit:
        infile.close()
        del infile
    return

def append_uncalled_positions_to_VCF(input_vcf,output_vcf,target_bed,reference_fasta="",reference_base_fasta=""):
    # VCF
    outfile,close_output_on_exit = _check_file_object(filename=output_vcf,mode="w")
    infile,close_input_on_exit  = _check_file_object(filename=input_vcf,mode="r")
    
    # Read input:
    vcf_input_lines = infile.readlines()
    vcf_input_line_fields = [input_line.rstrip("\n").split("\t") for input_line in vcf_input_lines if input_line[0] != "#"]
    
    # Write VCF header to output and extract samples to process
    for input_line in vcf_input_lines:
        if input_line[0] == "#":
            outfile.write(input_line)
            if (input_line.find("#CHROM") == 0):
                samples = input_line.split("\t")[9:]
        else:
            # stop when first content line has been reached
            break
    
    # Check if all positions were examined (positions with no alignments will be missing); add missing data to such positions
    tmp_bed = tempfile.NamedTemporaryFile(mode='w', suffix='.bed',delete=False)
    tmp_bed.close()
    child = subprocess.Popen(args="bedtools intersect -v -a "+target_bed+" -b stdin > "+tmp_bed.name,stdin=subprocess.PIPE,shell=True,universal_newlines=True)
    child.communicate(input="\n".join(vcf_input_lines))
    # Extract ref bases:
    tmp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta',delete=False)
    tmp_fasta.close()
    if reference_fasta != "":
        subprocess.call(args="bedtools getfasta -fi "+reference_fasta+" -fo "+tmp_fasta.name+" -bed "+ tmp_bed.name,shell=True)
    elif reference_base_fasta != "":
        shutil.copyfile(reference_base_fasta,tmp_fasta.name)
    # TODO: test
    # Make VCF for missing positions
    ploidy=2
    uncalled_bed_file = open(tmp_bed.name,"r")
    uncalled_ref_bases = SeqIO.parse(tmp_fasta.name,"fasta")
    vcf_uncalled_lines = []
    for bed_line in uncalled_bed_file.readlines():
        ref_base = next(uncalled_ref_bases)
        bed_fields = bed_line.strip().split("\t")
        vcf_line = []
        vcf_line.append(bed_fields[0])
        vcf_line.append(str(int(bed_fields[1])+1))
        vcf_line.append(bed_fields[0]+"_"+str(int(bed_fields[1])+1))
        vcf_line.append(str(ref_base.seq).upper())
        vcf_line.append(".")
        vcf_line.append("0")
        vcf_line.append(".")
        vcf_line.append("DP=0;DPB=0;EPPR=0;GTI=0;MQMR=0;NS=0;NUMALT=0;ODDS=0;PAIREDR=0;PQR=0;PRO=0;QR=0;RO=0;RPPR=0;MOTIF="+str(ref_base.seq).upper())
        vcf_line.append("GT:GQ:DP:DPR:RO:QR:AO:QA")
        for sample in samples:
            vcf_line.append(("./"*(ploidy-1))+".:.:.:.:.:.:.:.")
        vcf_uncalled_lines.append(vcf_line)
    uncalled_bed_file.close()
    
    # Merge and sort VCF content
    from operator import itemgetter
    vcf_input_line_fields.extend(vcf_uncalled_lines)
    sorted_vcf_input_line_fields = natsort.natsorted(vcf_input_line_fields, key=itemgetter(0,1))
    for current_line in sorted_vcf_input_line_fields:
        outfile.write(("\t".join(current_line))+"\n")
    
    # done
    os.remove(tmp_fasta.name)
    os.remove(tmp_bed.name)
    if close_output_on_exit:
        outfile.close()
        del outfile
    if close_input_on_exit:
        infile.close()
        del infile
    return
    

def Freebayes_VCF_allele_count_hard_filter(input_vcf,output_vcf,min_total_count,min_het_ab,logger,sample_string,remove_monomorphic=True,min_genotype_quality=0.0,ploidy=2):
    ## Applies hard filters on Freebayes calls:
    # min_total_count: samples with total read counts below this threshold gets their genotype set to .
    # min_het_ab: samples with a het call and AB below this threshold gets a hom genotype instead
    
    logger.info(sample_string+" - Applying hard allele count and het AB filters")
    
    # Scan VCF:
    outfile,close_output_on_exit = _check_file_object(filename=output_vcf,mode="w")
    infile,close_input_on_exit  = _check_file_object(filename=input_vcf,mode="r")
    
    current_line = infile.readline().strip()
    samples = []
    truncated_genotypes = 0
    forced_ref_genotypes = 0
    forced_alt_genotypes = 0
    parsed_variants = 0
    removed_variants = 0
    while current_line != "":
        include_line=True
        if current_line[0:1] == "#":
            # Propagate info
            output_line = current_line
            if current_line[0:6] == "#CHROM":
                # Get samples
                line_parts = current_line.split("\t")
                samples = line_parts[9:]
                # Add MinAF info field description:
                output_line ="##INFO=<ID=MinAF,Number=A,Type=Float,Description=\"Estimated minor allele frequency in the range (0,0.5]\">\n"+output_line
        else:
            # Parse the VCF data line and check filters
            line_parts = current_line.split("\t")
            variant_id = line_parts[0]+"_"+line_parts[1]
            info_fields = line_parts[7].split(";")
            info_field_data = {}
            for info_field in info_fields:
                info_field_parts = info_field.split("=")
                info_field_data[info_field_parts[0]] = info_field_parts[1]
            format_parts = line_parts[8].split(":")
            sample_data = line_parts[9:]
            genotype_field_index = format_parts.index("GT")
            GQ_field_index = format_parts.index("GQ")
            RO_field_index = format_parts.index("RO")
            AO_field_index = format_parts.index("AO")
            
            # re-calculate these
            line_DP = 0
            line_AO = 0
            line_RO = 0
            
            # Check all samples:
            updated_sample_data = []
            sample_genotypes = []
            for sample_name in samples:
                current_sample_data = sample_data[samples.index(sample_name)].split(":")
                sample_genotype_data = current_sample_data[genotype_field_index].replace(".","").split("/")
                # Only work on non-missing data
                if sample_genotype_data[0] != "":
                    # Collect counts
                    sample_ref_count = int(sample_data[samples.index(sample_name)].split(":")[RO_field_index])
                    sample_alt_counts = sample_data[samples.index(sample_name)].split(":")[AO_field_index].replace(".","").split(",")
                    total_count = sample_ref_count
                    for sample_alt_count in sample_alt_counts:
                        if sample_alt_count != "":
                            total_count += int(sample_alt_count)
                    
                    # Check filters
                    if (total_count < min_total_count) or (float(current_sample_data[GQ_field_index]) < min_genotype_quality):
                        # Low sample total coverage; set data as missing:
                        updated_sample_data.append(".:.:"+(":".join(current_sample_data[2:])))
                        truncated_genotypes += 1
                        sample_genotypes.append(None) # missing
                    else:
                        # add to line counts
                        line_DP += total_count
                        line_RO += sample_ref_count
                        line_AO += sum([int(sample_alt_count) for sample_alt_count in sample_alt_counts if sample_alt_count != ""])
                        
                        # check and update genotypes
                        if ploidy == 1:
                            if (min_het_ab > 0.0) and (sample_genotype_data[0] != "0"):
                                # Check call 
                                if float(sample_alt_counts[0])/float(total_count) < min_het_ab:
                                    # low AB, set ref homozygote
                                    updated_sample_data.append("0:"+(":".join(current_sample_data[1:])))
                                    forced_ref_genotypes += 1
                                    sample_genotypes.append(0) # hom ref
                                else:
                                    # High AB, keep alt homozygote
                                    updated_sample_data.append(sample_genotype_data[0]+":"+(":".join(current_sample_data[1:])))
                                    sample_genotypes.append(1) # keep alt
                            else:
                                # Do not modify current data
                                updated_sample_data.append(":".join(current_sample_data))
                                sample_genotypes.append(int(sample_genotype_data[0])) 
                        else:
                            if len(set(sample_genotype_data)) == 2:
                                if min_het_ab > 0.0:
                                    # Biallelic heterozygote, check het AB:
                                    if float(sample_alt_counts[0])/float(total_count) < min_het_ab:
                                        # low AB, set ref homozygote
                                        updated_sample_data.append("0/0:"+(":".join(current_sample_data[1:])))
                                        forced_ref_genotypes += 1
                                        sample_genotypes.append(0) # hom ref
                                    elif float(sample_alt_counts[0])/float(total_count) > (1-min_het_ab):
                                        # high AB, set alt homozygote
                                        updated_sample_data.append("1/1:"+(":".join(current_sample_data[1:])))
                                        forced_alt_genotypes += 1
                                        sample_genotypes.append(1) # hom alt
                                    else:
                                        # AB OK:
                                        updated_sample_data.append(":".join(current_sample_data))
                                        sample_genotypes.append(2) # het
                                else:
                                    # Do not modify current data
                                    updated_sample_data.append(":".join(current_sample_data))
                                    sample_genotypes.append(2) # het
                            elif len(set(sample_genotype_data)) == 1:
                                # Homozygote
                                sample_genotypes.append(int(list(set(sample_genotype_data))[0])) # hom
                                updated_sample_data.append(":".join(current_sample_data))
                            else:
                                # Multiallelic heterozygote, assume all OK
                                updated_sample_data.append(":".join(current_sample_data))
                                sample_genotypes.append(2) # het
                else:
                    # Propagate missing data
                    if len(current_sample_data[1:]) > 0:
                        # there are further fields with data to include
                        updated_sample_data.append(":".join(["."]+current_sample_data[1:]))
                    else:
                        updated_sample_data.append(".")
                    sample_genotypes.append(None) # missing
            
            # Calculate new NS
            info_field_data["NS"] = len(sample_genotypes)-sample_genotypes.count(None)
            # Calculate new MinAF
            if (info_field_data["NS"] > 0):
                if len(samples) == 1:
                    info_field_data["MinAF"] = float(sample_genotypes.count(1)+(sample_genotypes.count(2)/2))/ploidy
                else:
                    info_field_data["MinAF"] = float((2*min(sample_genotypes.count(0),sample_genotypes.count(1)))+sample_genotypes.count(2))/float(2*info_field_data["NS"])
            else:
                info_field_data["MinAF"] = 0
            
            # New DP, RO and AO
            info_field_data["DP"] = line_DP
            info_field_data["RO"] = line_RO
            info_field_data["AO"] = line_AO
            
            logger.debug("Freebayes_VCF_allele_count_hard_filter() info_field_data[\"NS\"]: "+str(info_field_data["NS"]))
            logger.debug("Freebayes_VCF_allele_count_hard_filter() info_field_data[\"MinAF\"]: "+str(info_field_data["MinAF"]))
            logger.debug("Freebayes_VCF_allele_count_hard_filter() sample_genotypes: "+str(sample_genotypes))
            # Remove non-polymorphic loci (MinAF = 0 or hom-count = 0)
            if remove_monomorphic:
                if ((sample_genotypes.count(0)+sample_genotypes.count(1)) == 0) or (info_field_data["MinAF"] == 0):
                    include_line=False
                    removed_variants+=1
            
            # Construct output line
            line_parts[7] = ";".join([key+"="+str(info_field_data[key]) for key in info_field_data]) # replace INFO
            output_line = ("\t".join(line_parts[0:9]))+"\t"+("\t".join(updated_sample_data))
            
            # Stats
            parsed_variants += 1
            if (parsed_variants % 1000 == 0):
                logger.info(sample_string+" - Freebayes_VCF_allele_count_hard_filter() RUNNING - scanned %d variants; truncated genotypes: %d, forced reference homozygotes: %d, forced alternate homozygotes: %d, removed variants: %d" % (parsed_variants,truncated_genotypes,forced_ref_genotypes,forced_alt_genotypes,removed_variants))
            
        # Write output
        if include_line:
            outfile.write(output_line+"\n")
        
        # Next line
        current_line = infile.readline().strip()
    
    # Done
    logger.info(sample_string+" - Freebayes_VCF_allele_count_hard_filter() DONE - scanned %d variants; truncated genotypes: %d, forced reference homozygotes: %d, forced alternate homozygotes: %d, removed variants: %d" % (parsed_variants,truncated_genotypes,forced_ref_genotypes,forced_alt_genotypes,removed_variants))
    
    # Clean up
    if close_output_on_exit:
        outfile.close()
        del outfile
    if close_input_on_exit:
        infile.close()
        del infile
    
    # Done
    return


def polyploid_Freebayes_VCF_allele_count_hard_filter(input_vcf,output_vcf,min_total_count,logger,sample_string,remove_monomorphic=True,min_genotype_quality=0.0):
    ## Applies hard filters on Freebayes calls:
    # min_total_count: samples with total read counts below this threshold gets their genotype set to .
    
    logger.info(sample_string+" - Applying hard allele count filter on '"+str(input_vcf)+"'.")
    # Scan VCF:
    outfile,close_output_on_exit = _check_file_object(filename=output_vcf,mode="w")
    infile,close_input_on_exit  = _check_file_object(filename=input_vcf,mode="r")
    
    current_line = infile.readline().strip()
    samples = []
    truncated_genotypes = 0
    forced_ref_genotypes = 0
    forced_alt_genotypes = 0
    parsed_variants = 0
    removed_variants = 0
    while current_line != "":
        include_line=True
        if current_line[0:1] == "#":
            # Propagate info
            output_line = current_line
            if current_line[0:6] == "#CHROM":
                # Get samples
                line_parts = current_line.split("\t")
                samples = line_parts[9:]
                # Add MinAF info field description:
                output_line ="##INFO=<ID=MinAF,Number=A,Type=Float,Description=\"Estimated minor allele frequency in the range (0,0.5]\">\n"+output_line
        else:
            # Parse the VCF data line and check filters
            line_parts = current_line.split("\t")
            variant_id = line_parts[0]+"_"+line_parts[1]
            info_fields = line_parts[7].split(";")
            info_field_data = {}
            for info_field in info_fields:
                info_field_parts = info_field.split("=")
                info_field_data[info_field_parts[0]] = info_field_parts[1]
            format_parts = line_parts[8].split(":")
            sample_data = line_parts[9:]
            genotype_field_index = format_parts.index("GT")
            GQ_field_index = format_parts.index("GQ")
            RO_field_index = format_parts.index("RO")
            AO_field_index = format_parts.index("AO")
            
            # re-calculate these
            line_DP = 0
            line_AO = 0
            line_RO = 0
            
            # Check all samples:
            updated_sample_data = []
            sample_genotypes = []
            current_variant_updated = False
            for sample_name in samples:
                current_sample_data = sample_data[samples.index(sample_name)].split(":")
                sample_genotype_data = current_sample_data[genotype_field_index].replace(".","").split("/")
                # Only work on non-missing data
                if sample_genotype_data[0] != "":
                    # Collect counts
                    sample_ref_count = int(current_sample_data[RO_field_index])
                    sample_alt_counts = current_sample_data[AO_field_index].replace(".","").split(",")
                    total_count = sample_ref_count
                    for sample_alt_count in sample_alt_counts:
                        if sample_alt_count != "":
                            total_count += int(sample_alt_count)
                    
                    # Check filters
                    if (total_count < min_total_count) or (float(current_sample_data[GQ_field_index]) < min_genotype_quality):
                        # Low sample total coverage; set data as missing:
                        updated_sample_data.append(".:.:"+(":".join(current_sample_data[2:])))
                        truncated_genotypes += 1
                        sample_genotypes.append(None) # missing
                        current_variant_updated = True
                    else:
                        # all OK
                        updated_sample_data.append(":".join(current_sample_data))
                        # hom/het
                        sample_genotypes.append([int(x) for x in sample_genotype_data]) # propagate
                        # add to line counts
                        line_DP += total_count
                        line_RO += sample_ref_count
                        line_AO += sum([int(sample_alt_count) for sample_alt_count in sample_alt_counts if sample_alt_count != ""])
                else:
                    # Propagate missing data
                    if len(current_sample_data[1:]) > 0:
                        # there are further fields with data to include
                        updated_sample_data.append(":".join(["."]+current_sample_data[1:]))
                    else:
                        updated_sample_data.append(".")
                    sample_genotypes.append(None) # missing
            
            # Calculate new NS
            info_field_data["NS"] = len(sample_genotypes)-sample_genotypes.count(None)
            
            # Calculate new AF and MinAF
            current_allele_counts = {}
            for sample_genotype in sample_genotypes:
                if not sample_genotype is None:
                    for allele in sample_genotype:
                        if not allele in current_allele_counts:
                            current_allele_counts[allele]=0
                        current_allele_counts[allele]+=1 
            
            total_alleles = sum([x for x in current_allele_counts.values()])
            if (line_parts[4] != ".") and total_alleles > 0:
                current_af = []
                for current_allele_index in range(0,1+len(line_parts[4].replace(";",",").split(","))):
                    if current_allele_index in current_allele_counts:
                        current_af.append(current_allele_counts[current_allele_index]/float(total_alleles))
                    else:
                        current_af.append(0.00)
                
                info_field_data["AF"] = ",".join([str(x) for x in current_af[1:]])
                info_field_data["MinAF"] = min(current_af)
                
            else:
                info_field_data["MinAF"] = 0.00
                info_field_data["AF"] = "0.00"
            
            # Update DP, RO and AO
            info_field_data["DP"] = line_DP
            info_field_data["RO"] = line_RO
            info_field_data["AO"] = line_AO
            
            # Truncate certain INFO fields; we cannot easily recalculate them
            if current_variant_updated:
                info_field_data["AC"] = "."
                info_field_data["AN"] = "."
                info_field_data["AB"] = "."
                info_field_data["ABP"] = "."
                info_field_data["DPRA"] = "."
                info_field_data["NUMALT"] = "."
                info_field_data["MEANALT"] = "."
                info_field_data["MQM"] = "."
                info_field_data["PAIRED"] = "."
                info_field_data["PAIREDR"] = "."
            
            # Remove non-polymorphic loci (MinAF = 0)
            if remove_monomorphic:
                if (info_field_data["MinAF"] == 0):
                    include_line=False
                    removed_variants+=1
            
            # Construct output line
            line_parts[7] = ";".join([key+"="+str(info_field_data[key]) for key in info_field_data]) # replace INFO
            output_line = ("\t".join(line_parts[0:9]))+"\t"+("\t".join(updated_sample_data))
            
            # Stats
            parsed_variants += 1
            if (parsed_variants % 1000 == 0):
                logger.info(sample_string+" - polyploid_Freebayes_VCF_allele_count_hard_filter() RUNNING - scanned %d variants; truncated genotypes: %d, removed variants: %d" % (parsed_variants,truncated_genotypes,removed_variants))
            
        # Write output
        if include_line:
            outfile.write(output_line+"\n")
        
        # Next line
        current_line = infile.readline().strip()
    
    # Done
    logger.info(sample_string+" - polyploid_Freebayes_VCF_allele_count_hard_filter() DONE - scanned %d variants; truncated genotypes: %d, removed variants: %d" % (parsed_variants,truncated_genotypes,removed_variants))
    
    # Clean up
    if close_output_on_exit:
        outfile.close()
        del outfile
    if close_input_on_exit:
        infile.close()
        del infile
    
    # Done
    return
    

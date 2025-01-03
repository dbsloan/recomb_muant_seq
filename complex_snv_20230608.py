import sys
import csv
csv.field_size_limit(sys.maxsize)
import re
import os
import copy
import itertools
from itertools import chain
import subprocess


from Bio.Blast.Applications import NcbiblastnCommandline


def complex_snv():

	input_sam = sys.argv[1]
	input_fasta = sys.argv[2]
	output_complex_muts = sys.argv[3]
	output_snv_muts = sys.argv[4]
	output_depth = sys.argv[5]
	rep_ID =sys.argv[6]
	excludefromend = int(sys.argv[7])
	recomb_flank=int(sys.argv[8])
	input_plastid_coord=sys.argv[9]
	input_mito_coord=sys.argv[10]
	output_total_cov_and_snvs=sys.argv[11]
	output_total_cov_and_dinucs=sys.argv[12]
	output_region=sys.argv[13]
	output_strand=sys.argv[14]
	output_filtered_stats=sys.argv[15]
	output_trinuc=sys.argv[16]
	output_trinuc_strand_aware=sys.argv[17]
	

	
	'''
	This script is designed to parse a sam_depth file and report dinucleotide and other interesting mutations. It includes a blastnn based check for numt derived muts.
	It is called with:
	
	python3 complex_snv_20230608.py	
	followed by these arrguments, in order, seperated by spaces
	
	1. input sam file
	2. fasta with identical chrom names to those used in sam generation
	3. output list of dinucleotide mutations
	4. output list of single nucleotide mutations
	5. custom depth file (like samtools depth, but calculated using only 'non-discard' reads)
	6. replicate ID
	7. exclude from end (integer)
	8. recomb_flank (use 100)
	9. plastid coord file (download here: https://github.com/dbsloan/duplexseq/tree/master/example_files)
	10. mito coord file, download from above link
	11. output_total_cov_and_snvs
	12. output_total_cov_and_dinucs
	13. output_region
	14. output_strand
	
	
	note blast db must be created for the references fasta (used for sam generation):
	makeblastdb -in refs_cp28673mod.fas -parse_seqids -dbtype nucl -out refs_cp28673mod_db
	
	and for the contamination fasta (containing all numts and nuclear genome):
	makeblastdb -in tair10_nuclear_and_numt.fas -parse_seqids -dbtype nucl -out tair10_nuclear_and_numt.db

	example call(s):
	
	python3 complex_snv_20230608.py mt_MSH1CS3372_mut_2.DCS.filt.test.sam refs_cp28673mod.fas mt_MSH1CS3372_mut_2.comp_mut.txt  mt_MSH1CS3372_mut_2.snv.txt mt_MSH1CS3372_mut_2.cus_depth.txt    mt_MSH1CS3372_mut_2 9 100 NC_000932.mod28673.coordmap.txt NC_037304.coordmap.txt mt_MSH1CS3372_mut_2.snvs_and_cov.txt mt_MSH1CS3372_mut_2.dinucs_and_cov.txt mt_MSH1CS3372_mut_2.snvs_region.txt mt_MSH1CS3372_mut_2.snvs_strand.txt mt_MSH1CS3372_mut_2.filtered_stats.txt  mt_MSH1CS3372_mut_2.trinuc_freq.txt mt_MSH1CS3372_mut_2.trinuc_freq_strand_aware.txt


	'''
	
	try: 
		fasta_handle=open(input_fasta)
	except:
		return 'unable to open input fasta'
		
	fasta_reader = csv.reader(fasta_handle, delimiter = '\t')
	
	try: 
		sam_handle=open(input_sam)
	except:
		return 'unable to open input sam'
		
	sam_reader = csv.reader(sam_handle, delimiter = '\t')
	
	try: 
		plastid_coord_handle=open(input_plastid_coord)
	except:
		return 'unable to open input plastid coord'
		
	plastid_coord_reader = csv.reader(plastid_coord_handle, delimiter = '\t')
	
	try: 
		mito_coord_handle=open(input_mito_coord)
	except:
		return 'unable to open input mito coord'
		
	mito_coord_reader = csv.reader(mito_coord_handle, delimiter = '\t')
		
	try: 
		out_handle_comp=open(output_complex_muts, 'w')
	except:
		return 'unable to create output comp mut'
	
	try: 
		out_handle_snv=open(output_snv_muts, 'w')
	except:
		return 'unable to create output snv'
	 
	try: 
		out_handle_depth=open(output_depth, 'w')
	except:
		return 'unable to create output depth'

	try: 
		out_handle_total_cov_and_snvs=open(output_total_cov_and_snvs, 'w')
	except:
		return 'unable to create output total_cov_and_snvs'

	try: 
		out_handle_total_cov_and_dinucs=open(output_total_cov_and_dinucs, 'w')
	except:
		return 'unable to create output total_cov_and_dinucs'
	
	try:
		out_handle_region=open(output_region, 'w')
	except:
		return 'unable to create output output_region'
	
	try:
		out_handle_strand=open(output_strand, 'w')
	except:
		return 'unable to create output output_strand'

	try:
		out_handle_filtered=open(output_filtered_stats, 'w')
	except:
		return 'unable to create output filtered stats'

	try:
		out_handle_trinuc=open(output_trinuc, 'w')
	except:
		return 'unable to create output trinuc'
		
	
	try:
		out_handle_trinuc_strnd_aware=open(output_trinuc_strand_aware, 'w')
	except:
		return 'unable to create output trinuc'
					
	
# Read in the Fasta file and turn it into a dictionary where all the keys are the fasta headers (chrom names) and 
# All the values are matrixes with 8 columns with the following attributes:
# The matrix is as many rows an the length of the sequence
# The matrix contains 10 columns: 1:position, 2:reference base, 3:coverage of reference(#), 4-8: info from the coordmap file - like genomic region, strand and gene name. 9 and 10 are as of yet empty.
	seq_dict={}
	#tempseq=''
	tempchrom='ignore'
	tempseq_mat=[]
	pos_counter=0
	fasta_dict={}
	temp_seq=''
	for line in fasta_reader:
		if line[0][0]=='>':
			#print(line)
			seq_dict[tempchrom]=seq_dict.get(tempchrom,tempseq_mat)
			fasta_dict[tempchrom]=fasta_dict.get(tempchrom,temp_seq)
			tempseq_mat=[]
			temp_seq=''
			tempchrom=line[0][1:]
			#print(tempchrom)
			#print(1)
			pos_counter=0
		else:
			for char in line[0]:
				if char != '\n':
					tempseq_mat.append([0]*10)
					tempseq_mat[pos_counter][0]=pos_counter+1
					tempseq_mat[pos_counter][1]=char
					pos_counter +=1
					temp_seq+=char
	else:
		seq_dict[tempchrom]=seq_dict.get(tempchrom,tempseq_mat)
		fasta_dict[tempchrom]=fasta_dict.get(tempchrom,temp_seq)
		
		
	


	dict_keys= ['mito','plastid']
	
	pointer={'mito':mito_coord_reader,'plastid':plastid_coord_reader}
	
	for key in dict_keys:
		#if key == 'mito':
		index=-2
		for line in pointer[key]:
			index+=1
			if index>-1:
				#
				seq_dict[key][index][3]=line[2]
				seq_dict[key][index][4]=line[3]
				seq_dict[key][index][5]=line[4]
				seq_dict[key][index][6]=line[5]
				seq_dict[key][index][7]=line[6]
		
	
	trinuc_count={}
		
	for chrom in seq_dict:
		
		if chrom != 'ignore':
		

			trinuc_count[chrom]=trinuc_count.get(chrom,{})
			trinuc=''
			oldtrinuc=''
			first='off'
			for line in seq_dict[chrom]:

				
				if len(trinuc) <3:
					trinuc+=line[1]
				elif len(trinuc) ==3:
					oldtrinuc=trinuc
					
					trinuc=trinuc[1:]+line[1]
				
				pos=int(line[0])
					
				if len(trinuc)==3:
					if(int(line[0])) > 2:
				
						trinuc_count[chrom][(pos-3,pos-2,pos-1)]=trinuc_count[chrom].get((pos-3,pos-2,pos-1),[oldtrinuc,line[3],'',0,0,0,'NA'])

						if line[3] != 'intergenic':
							trinuc_count[chrom][(pos-3,pos-2,pos-1)][6]=line[5]


	
	snvcount=0
	dinuc_count=0
	twosnvsnotadj=0
	tri_nuc_count=0
	dinuc_plussnv=0
	threesnvs_notadj=0
	other=0
	temp_dict={}



	recom_hit_dict={} ### this will be a complicated dictionary of dictionaries for storing info about mutations explained by recombinations
	SNV_hit_dict={}
	di_nuc_hit_dict={}
	
	snv_recomb_count=0
	contamcount=0
	contam_count_snv=0
	contam_count_di=0
	contam_count_two_snvs=0
	contam_count_tri=0
									
	
	
	dinuc_recomb=0
	for line in sam_reader:
		contamhit=''
		
		recombhit=''
		
		if line[0][0]!='@':                                            #header line

			if line[1] != 4: #column 1 is 4 if the read didn't map
				if re.fullmatch("[0-9]+M{1}", line[5]): #only do all this stuff beneath if the cigar sting is #M
					if line[-2][5:].isnumeric(): #no SNPs in MDZ field
					
	
						#add one to the depth counter (column 4)
						for nt in range(int(line[3])-1,(  int(line[3])+len(line[9]) -excludefromend)):
							if nt >= int(line[3]) + excludefromend:
								
								seq_dict[line[2]][nt-1][2]+=1
					else: 
					
						mdz= ["".join(x) for _, x in itertools.groupby(line[-2][5:], key=str.isdigit)]
						if int(mdz[0])>excludefromend and int(mdz[-1]) > excludefromend:
						
							# test= subprocess.call(,shell=True)
# 							print(test)
							query = line[9]
							blastn_temp=NcbiblastnCommandline(db='tair10_nuclear_and_numt.db',outfmt=7)
							out,err=blastn_temp(stdin=query)
							hits = [i for i in (out.split('\n')) if i.startswith('Query_1')]
							for hit in hits:
								if hit.split('\t')[2]=='100.000' and int(hit.split('\t')[3])>max(int(mdz[0])+2,int(mdz[-1])+2):
									#print(out)
									#print(line)
									contamhit='yes'
									contamcount+=1
									if len(mdz)==3:
										contam_count_snv+=1
									elif len(mdz)==5:
										if mdz[2]=='0':
											contam_count_di+=1
										else:
											contam_count_two_snvs+=1
										
									elif len(mdz)==7:
										contam_count_tri+=1
									
									
									
							if contamhit=='':
								if len(mdz)==3:
									### these are single nucleotide variants
									### check for recombination mediated artifacts by simulating a read - with x bps of flanking seq (x = user specified)on either side that lacks the mutation. If the mut allele has a WT copy elsewhere in the reference, it's more likely the allele is a result fo recombination - not mutation 
									leftflank=max(int(line[3]) + int(mdz[0]) - recomb_flank-1,0)
									rightflank=min(int(line[3]) + int(mdz[0]) + recomb_flank,len(fasta_dict[line[2]]))

									query=fasta_dict[line[2]][leftflank:rightflank]
									blastn_temp=NcbiblastnCommandline(db='refs_cp28673mod_db',outfmt=5,word_size=11,evalue=1e-10 )
									out,err=blastn_temp(stdin=query)
									
									key=''
									
									#print(out)
									
									### before conducting the blast- check if this allele has already been identified as recombination derrived:
									alt_allele= line[9][int(mdz[0])]
									#print(alt_allele)
									
									recom_hit_dict[line[2]]=recom_hit_dict.get(line[2],{})
									if int(line[3]) + int(mdz[0]) in recom_hit_dict[line[2]]:
										if alt_allele in recom_hit_dict[line[2]][int(line[3]) + int(mdz[0])]:
											recombhit='yes'
											#print('already flagged')
										else:
											print('multiple snv types at positions explained by recombination')
									else:
										
										### first i do a little gymnastics to ignore the hits for the other reference genome (like when a plastid mapping snp has a hit in the mito genome)
										for blasthit in out.split('<Hit_id>'):
											
											if blasthit.split('<')[0]==line[2]:
												#print(blasthit)
												### then I do a check of each blast hit - to see if there is a different location in the genome where the 'mut' allele is part of the WT reference sequence.
												for hsp_num in blasthit.split('<Hsp_num>'):
													for entry in hsp_num.split('\n'):
														entry=entry.replace(' ','')
														entry_list=re.split('>|<',entry)
														entry_list=[i for i in entry_list if i != '']
														if len(entry_list)>0:	
													
															if entry_list[0].isnumeric():
																### each unique blast hit gets a numeric ID. I use these IDs as dictionary keys, to efficiently store all the info for a given blast hit
																key=entry_list[0]
																temp_dict={}
																temp_dict[key]=temp_dict.get(key,{})

															else:
																if key!='':
																	if len(entry_list)>1:
																		#here is where I populate the blast hit dictionary. 
																
																		temp_dict[key][entry_list[0]]=temp_dict[key].get(entry_list[0],entry_list[1])
																		
													
													if key !='': ### This is a hacked way of deetermining we have a occupied temp_dict (containing a blast hit) Below I assess if the putative SNP could be explained by recombination (putative mut exists as WT allele elsewhere in genome)
														ref_coord=int(temp_dict[key]['Hsp_query-from'])
														target_position=min(recomb_flank+1,(int(line[3])+int(mdz[0])))
														alt_allele= line[9][int(mdz[0])]
														for i in range(0,len(temp_dict[key]['Hsp_qseq'])):
															substring=temp_dict[key]['Hsp_qseq'][i:i+1]
															#print(substring)
															if ref_coord==target_position:
																if temp_dict[key]['Hsp_hseq'][i:i+1]==alt_allele:
# 																	print(substring)
# 																	print(alt_allele)
# 																	print(temp_dict[key]['Hsp_hseq'][i:i+1])
# 																	print(hsp_num)
# 																	print(query)
																	recombhit='yes'
																	recom_hit_dict[line[2]]=recom_hit_dict.get(line[2],{})
																	recom_hit_dict[line[2]][int(line[3]) + int(mdz[0])]=recom_hit_dict[line[2]].get(int(line[3]) + int(mdz[0]),{})
																	recom_hit_dict[line[2]][int(line[3]) + int(mdz[0])][alt_allele]=recom_hit_dict[line[2]][int(line[3]) + int(mdz[0])].get(alt_allele,[key,temp_dict[key]['Hsp_hit-from'],temp_dict[key]['Hsp_align-len']])
															ref_coord+=1
															if substring== '-':
																target_position+=1
														
									if recombhit=='':
										
										SNV_hit_dict[line[2]]=SNV_hit_dict.get(line[2],{})
										SNV_hit_dict[line[2]][int(line[3]) + int(mdz[0])] = SNV_hit_dict[line[2]].get(int(line[3]) + int(mdz[0]),{})
										SNV_hit_dict[line[2]][int(line[3]) + int(mdz[0])][alt_allele]=SNV_hit_dict[line[2]][int(line[3]) + int(mdz[0])].get(alt_allele,[mdz[1],0])
										
										SNV_hit_dict[line[2]][int(line[3]) + int(mdz[0])][alt_allele][1] +=1
										snvcount+=1		
										
										
										
										
										for nt in range(int(line[3])-1,(  int(line[3])+len(line[9]) -excludefromend)):
											if nt >= int(line[3]) + excludefromend:
												if nt != int(line[3]) + int(mdz[0]):
													seq_dict[line[2]][nt-1][2]+=1	
									else:
										snv_recomb_count+=1		

								elif len(mdz)==5:
									if mdz[2]=='0':
										###these are dinucleotide mutations
										
										### check for recombination mediated artifacts by simulating a read - with x bps of flanking seq (x = user specified)on either side that lacks the mutation. If the mut allele has a WT copy elsewhere in the reference, it's more likely the allele is a result fo recombination - not mutation 
										leftflank=max(int(line[3]) + int(mdz[0]) - recomb_flank-1,0)
										rightflank=min(int(line[3]) + int(mdz[0]) + recomb_flank+1,len(fasta_dict[line[2]]))

										query=fasta_dict[line[2]][leftflank:rightflank]
										blastn_temp=NcbiblastnCommandline(db='refs_cp28673mod_db',outfmt=5,word_size=11,evalue=1e-10)
										out,err=blastn_temp(stdin=query)

										key=''
									
										### before conducting the blast- check if this allele has already been identified as recombination derrived:
										alt_allele= line[9][int(mdz[0]):int(mdz[0])+2]
										
										recom_hit_dict[line[2]]=recom_hit_dict.get(line[2],{})
										if int(line[3]) + int(mdz[0]) in recom_hit_dict[line[2]]:
											if alt_allele in recom_hit_dict[line[2]][int(line[3]) + int(mdz[0])]:
												recombhit='yes'
												#print('already flagged')
											else:
												print('multiple snv types at positions explained by recombination')
										else:
											### first i do a little gymnastics to ignore the hits for the other reference genome (like when a plastid mapping snp has a hit in the mito genome)
											for blasthit in out.split('<Hit_id>'):
												if blasthit.split('<')[0]==line[2]:
													### then I do a check of each blast hit - to see if there is a different location in the genome where the 'mut' allele is part of the WT reference sequence.
													for hsp_num in blasthit.split('<Hsp_num>'):
														for entry in hsp_num.split('\n'):
															entry=entry.replace(' ','')
															entry_list=re.split('>|<',entry)
															entry_list=[i for i in entry_list if i != '']
															if len(entry_list)>0:	
													
																if entry_list[0].isnumeric():
																	### each unique blast hit gets a numeric ID. I use these IDs as dictionary keys, to efficiently store all the info for a given blast hit
																	key=entry_list[0]
																	temp_dict={}
																	temp_dict[key]=temp_dict.get(key,{})

																else:
																	if key!='':
																		if len(entry_list)>1:
																			#here is where I populate the blast hit dictionary. 
																			temp_dict[key][entry_list[0]]=temp_dict[key].get(entry_list[0],entry_list[1])
													
													
														if key !='': ### This is a hacked way of deetermining we have a occupied temp_dict (containing a blast hit) 
														###Below I assess if the putative SNP could be explained by recombination (putative mut exists as WT allele elsewhere in genome)
															ref_coord=int(temp_dict[key]['Hsp_query-from'])
															target_position=min(recomb_flank+1,(int(line[3])+int(mdz[0])))
															alt_allele= line[9][int(mdz[0]):int(mdz[0])+2]
															for i in range(0,len(temp_dict[key]['Hsp_qseq'])):
																substring=temp_dict[key]['Hsp_qseq'][i:i+2]
																if ref_coord==target_position:
																	if temp_dict[key]['Hsp_hseq'][i:i+2]==alt_allele:

																		recombhit='yes'
																		recom_hit_dict[line[2]]=recom_hit_dict.get(line[2],{})
																		recom_hit_dict[line[2]][int(line[3]) + int(mdz[0])]=recom_hit_dict[line[2]].get(int(line[3]) + int(mdz[0]),{})
																		recom_hit_dict[line[2]][int(line[3]) + int(mdz[0])][alt_allele]=recom_hit_dict[line[2]][int(line[3]) + int(mdz[0])].get(alt_allele,[key,temp_dict[key]['Hsp_hit-from'],temp_dict[key]['Hsp_align-len']])
																			

																			
																ref_coord+=1
																if substring[0]=='-':
																	target_position+=1
														
										if recombhit=='':
											di_nuc_hit_dict[line[2]]=di_nuc_hit_dict.get(line[2],{})
											di_nuc_hit_dict[line[2]][int(line[3]) + int(mdz[0])] = di_nuc_hit_dict[line[2]].get(int(line[3]) + int(mdz[0]),{})
											di_nuc_hit_dict[line[2]][int(line[3]) + int(mdz[0])][alt_allele]=di_nuc_hit_dict[line[2]][int(line[3]) + int(mdz[0])].get(alt_allele,[mdz[1]+mdz[3],0])
											di_nuc_hit_dict[line[2]][int(line[3]) + int(mdz[0])][alt_allele][1]+=1
											dinuc_count+=1						
											
											
											
								
										
										
										
											for nt in range(int(line[3])-1,(  int(line[3])+len(line[9]) -excludefromend)):
												if nt >= int(line[3]) + excludefromend:
													if nt != int(line[3]) + int(mdz[0]) and nt != int(line[3]) + int(mdz[0]) +1 :
														seq_dict[line[2]][nt-1][2]+=1				
										else:
											dinuc_recomb+=1
										

# 											
										
									else:
									
										
										###two non adj snvs
										### check for recombination mediated artifacts by simulating a read - with x bps of flanking seq (x = user specified)on either side that lacks the mutation. If the mut allele has a WT copy elsewhere in the reference, it's more likely the allele is a result fo recombination - not mutation 
										#print(line)
										leftflank=max(int(line[3]) + int(mdz[0]) - recomb_flank-1,0)
										
										rightflank=min(int(line[3]) + int(mdz[0]) +int(mdz[2]) + recomb_flank+1,len(fasta_dict[line[2]]))

										query=fasta_dict[line[2]][leftflank:rightflank]
										blastn_temp=NcbiblastnCommandline(db='refs_cp28673mod_db',outfmt=5)
										out,err=blastn_temp(stdin=query)

										key=''
									
										### before conducting the blast- check if this allele has already been identified as recombination derrived:
										alt_allele_list= [line[9][int(mdz[0])],line[9][int(mdz[0])+1+int(mdz[2])]]
										#print(alt_allele_list)
										two_SNVs_hit_list=['','']
									# 	recom_hit_dict[line[2]]=recom_hit_dict.get(line[2],{})
# 										if int(line[3]) + int(mdz[0]) in recom_hit_dict[line[2]]:
# 											if alt_allele in recom_hit_dict[line[2]][int(line[3]) + int(mdz[0])]:
# 												recombhit='yes'
# 												#print('already flagged')
# 											else:
# 												print('multiple snv types at positions explained by recombination')
# 										else:
										### first i do a little gymnastics to ignore the hits for the other reference genome (like when a plastid mapping snp has a hit in the mito genome)
										for blasthit in out.split('<Hit_id>'):
											if blasthit.split('<')[0]==line[2]:
												### then I do a check of each blast hit - to see if there is a different location in the genome where the 'mut' allele is part of the WT reference sequence.
												for hsp_num in blasthit.split('<Hsp_num>'):
													for entry in hsp_num.split('\n'):
														entry=entry.replace(' ','')
														entry_list=re.split('>|<',entry)
														entry_list=[i for i in entry_list if i != '']
														if len(entry_list)>0:	
												
															if entry_list[0].isnumeric():
																### each unique blast hit gets a numeric ID. I use these IDs as dictionary keys, to efficiently store all the info for a given blast hit
																key=entry_list[0]
																temp_dict={}
																temp_dict[entry[0]]=temp_dict.get(entry[0],{})
															else:
																if key!='':
																	if len(entry_list)>1:
																		#here is where I populate the blast hit dictionary. 
																		temp_dict[key][entry_list[0]]=temp_dict[key].get(entry_list[0],entry_list[1])
												
												
													if key !='': ### This is a hacked way of deetermining we have a occupied temp_dict (containing a blast hit) 
													###Below I assess if the putative SNP could be explained by recombination (putative mut exists as WT allele elsewhere in genome)
														ref_coord=int(temp_dict[key]['Hsp_query-from'])
														target_positions_list=[min(recomb_flank+1,(int(line[3])+int(mdz[0]))), recomb_flank+2+int(mdz[2])]
														
														

														for i in range(0,len(temp_dict[key]['Hsp_qseq'])):
															substring=temp_dict[key]['Hsp_qseq'][i:i+2]
															index=-1
															
															for target_position in target_positions_list:
																index+=1

																
																if ref_coord==target_position:
																	if temp_dict[key]['Hsp_hseq'][i:i+1]==alt_allele_list[index]:

																		two_SNVs_hit_list[index]='yes'
																		recom_hit_dict[line[2]]=recom_hit_dict.get(line[2],{})
																		recom_hit_dict[line[2]][int(line[3]) + int(mdz[0])]=recom_hit_dict[line[2]].get(int(line[3]) + int(mdz[0]),{})
																		recom_hit_dict[line[2]][int(line[3]) + int(mdz[0])][alt_allele_list[index]]=recom_hit_dict[line[2]][int(line[3]) + int(mdz[0])].get(alt_allele_list[index],[key,temp_dict[key]['Hsp_hit-from'],temp_dict[key]['Hsp_align-len']])
															ref_coord+=1
															if substring=='-':
																target_position+=1
													
										if two_SNVs_hit_list[0]=='' and two_SNVs_hit_list[1]=='':
									
											twosnvsnotadj+=1
										elif two_SNVs_hit_list[0]=='' or two_SNVs_hit_list[1]=='' :
											print('one SNV, one recomb')
	
										
								elif len(mdz)==7:
									if mdz[2]=='0' and mdz[4]=='0':
										tri_nuc_count+=1
									elif mdz[2]=='0' or mdz[4]=='0':
										dinuc_plussnv+=1
									else:
										threesnvs_notadj+=1
														

								else:
									other+=1

							else:
								contamcount+=1
	
	type_count={}
	cov_count={}
	region_count={}
	strand_count={}
	total_cov={}
	all_pos_snv=['A>C','A>G','A>T','C>A','C>G','C>T','G>A','G>C','G>T','T>A','T>C','T>G']


	rep_ID_list=rep_ID.split('_')
	out_handle_filtered.write(rep_ID +'\t' +rep_ID_list[0] +'\t' +rep_ID_list[1] +'\t' +rep_ID_list[2] +'\t' +rep_ID_list[3] +'\t' + 'snv_contam_count' + '\t' +str(contam_count_snv) +'\n' )
	out_handle_filtered.write(rep_ID +'\t' +rep_ID_list[0] +'\t' +rep_ID_list[1] +'\t' +rep_ID_list[2] +'\t' +rep_ID_list[3] +'\t' + 'dinuc_contam_count' + '\t' +str(contam_count_di) +'\n' )
	out_handle_filtered.write(rep_ID +'\t' +rep_ID_list[0] +'\t' +rep_ID_list[1] +'\t' +rep_ID_list[2] +'\t' +rep_ID_list[3] +'\t' + 'two_nonadj_snv_contam_count' + '\t' +str(contam_count_two_snvs) +'\n' )
	out_handle_filtered.write(rep_ID +'\t' +rep_ID_list[0] +'\t' +rep_ID_list[1] +'\t' +rep_ID_list[2] +'\t' +rep_ID_list[3] +'\t' + 'three_snv_contam_count' + '\t' +str(contam_count_tri) +'\n' )
	out_handle_filtered.write(rep_ID +'\t' +rep_ID_list[0] +'\t' +rep_ID_list[1] +'\t' +rep_ID_list[2] +'\t' +rep_ID_list[3] +'\t' + '4plus_snv_contam_count' + '\t' +str(contamcount-contam_count_tri-contam_count_two_snvs-contam_count_di-contam_count_snv) +'\n' )
	out_handle_filtered.write(rep_ID +'\t' +rep_ID_list[0] +'\t' +rep_ID_list[1] +'\t' +rep_ID_list[2] +'\t' +rep_ID_list[3] +'\t' + 'snv_recomb_count' + '\t' +str(snv_recomb_count) +'\n' )
	out_handle_filtered.write(rep_ID +'\t' +rep_ID_list[0] +'\t' +rep_ID_list[1] +'\t' +rep_ID_list[2] +'\t' +rep_ID_list[3] +'\t' + 'dinuc_recomb_count' + '\t' +str(dinuc_recomb) +'\n' )
	out_handle_filtered.write(rep_ID +'\t' +rep_ID_list[0] +'\t' +rep_ID_list[1] +'\t' +rep_ID_list[2] +'\t' +rep_ID_list[3] +'\t' + 'snv_count' + '\t' +str(snvcount) +'\n' )
	out_handle_filtered.write(rep_ID +'\t' +rep_ID_list[0] +'\t' +rep_ID_list[1] +'\t' +rep_ID_list[2] +'\t' +rep_ID_list[3] +'\t' + 'dinuc_count' + '\t' +str(dinuc_count) +'\n' )
	
	
	#print(seq_dict)
	for chrom in seq_dict:
		
		if chrom != 'ignore':
		
			type_count[chrom]=type_count.get(chrom,{})
			cov_count[chrom]=cov_count.get(chrom,{})
			region_count[chrom]=region_count.get(chrom,{})
			strand_count[chrom]=strand_count.get(chrom,{})
			total_cov[chrom]=total_cov.get(chrom,0)


			for line in seq_dict[chrom]:



	
				pos=int(line[0])
				
				
				setter=(pos-1,pos,pos+1)
				if setter in trinuc_count[chrom]:
					trinuc_count[chrom][setter][3]+=line[2]

				
			
				cov_count[chrom][line[1]]=cov_count[chrom].get(line[1],[0,{}])
				cov_count[chrom][line[1]][0]+= line[2]
			
				region_count[chrom][line[3]]=region_count[chrom].get(line[3],{})
				region_count[chrom][line[3]][line[1]]=region_count[chrom][line[3]].get(line[1],[0,{}])
				region_count[chrom][line[3]][line[1]][0]+=line[2]
			
			
				strand_count[chrom][line[3]]=strand_count[chrom].get(line[3],{})
				strand_count[chrom][line[3]][line[1]]=strand_count[chrom][line[3]].get(line[1],{})
				strand_count[chrom][line[3]][line[1]][line[5]]=strand_count[chrom][line[3]][line[1]].get(line[5],[0,{}])
				strand_count[chrom][line[3]][line[1]][line[5]][0]+=line[2]
			
				total_cov[chrom]+=line[2]
			
				for pos_snv in all_pos_snv:
					if pos_snv[0] == line[1]:
						type_count[chrom][pos_snv]=type_count[chrom].get(pos_snv,[0,0])
						cov_count[chrom][line[1]][1][pos_snv]=cov_count[chrom][line[1]][1].get(pos_snv,[0,0])
						region_count[chrom][line[3]][line[1]][1][pos_snv]=region_count[chrom][line[3]][line[1]][1].get(pos_snv,[0,0])
						strand_count[chrom][line[3]][line[1]][line[5]][1][pos_snv]=strand_count[chrom][line[3]][line[1]][line[5]][1].get(pos_snv,[0,0])
			
			
				if chrom in SNV_hit_dict:
				
					if line[0] in SNV_hit_dict[chrom]:
						if len(SNV_hit_dict[chrom][line[0]]) ==1:
							for alt_allele in SNV_hit_dict[chrom][line[0]]:
								newkey=line[1] + '>' + alt_allele
						
								type_count[chrom][newkey][0]+=1
								type_count[chrom][newkey][1]+=SNV_hit_dict[chrom][line[0]][alt_allele][1]
						
								cov_count[chrom][line[1]][1][newkey][0]+=1
								cov_count[chrom][line[1]][1][newkey][1]+=SNV_hit_dict[chrom][line[0]][alt_allele][1]
						
								region_count[chrom][line[3]][line[1]][1][newkey][0]+=1
								region_count[chrom][line[3]][line[1]][1][newkey][1]+=SNV_hit_dict[chrom][line[0]][alt_allele][1]
					
								strand_count[chrom][line[3]][line[1]][line[5]][1][newkey][0]+=1
								strand_count[chrom][line[3]][line[1]][line[5]][1][newkey][1]+=SNV_hit_dict[chrom][line[0]][alt_allele][1]
						
	
						else:
							print(line)
							print(SNV_hit_dict[chrom][line[0]])
							print('above site has multiple mutations')
					


	
	for chrom in SNV_hit_dict:
		for pos in SNV_hit_dict[chrom]:
			post=int(pos)
			setter=(post-1,post,post+1)
			if setter in trinuc_count[chrom]:
				#print('hit')
				if len(SNV_hit_dict[chrom][pos])==1:
					
					for alt_allele in SNV_hit_dict[chrom][pos]:
						#print(SNV_hit_dict[chrom][pos])
						#print(pos)
						trinuc_count[chrom][setter][2]=alt_allele
						trinuc_count[chrom][setter][4]+=1
						trinuc_count[chrom][setter][5]+=SNV_hit_dict[chrom][pos][alt_allele][1]
						#print(trinuc_count[chrom][setter])
				else:
					
					print('above site has multiple mutations')
			
	
	

	tri_dict= {'AGA':'TCT','AGC':'GCT','AGG':'CCT','AGT':'ACT','CGA':'TCG','CGC':'GCG','CGG':'CCG','CGT':'ACG','GGA':'TCC','GGC':'GCC','GGG':'CCC','GGT':'ACC','TGA':'TCA','TGC':'GCA','TGG':'CCA','TGT':'ACA','ATA':'TAT','ATC':'GAT','ATG':'CAT','ATT':'AAT','CTA':'TAG','CTC':'GAG','CTG':'CAG','CTT':'AAG','GTA':'TAC','GTC':'GAC','GTG':'CAC','GTT':'AAC','TTA':'TAA' ,'TTC':'GAA' ,'TTG':'CAA' ,'TTT' :'AAA'}
	tri_dict2= {'mito': {'AAA':[0,{}], 'AAC':[0,{}],'AAG':[0,{}], 'AAT':[0,{}],'CAA':[0,{}],'CAC':[0,{}],'CAG':[0,{}],'CAT':[0,{}],'GAA':[0,{}],'GAC':[0,{}],'GAG':[0,{}],'GAT':[0,{}],'TAA':[0,{}],'TAC':[0,{}],'TAG':[0,{}],'TAT':[0,{}],'ACA':[0,{}],'ACC':[0,{}],'ACG':[0,{}],'ACT':[0,{}],'CCA':[0,{}],'CCC':[0,{}],'CCG':[0,{}],'CCT':[0,{}],'GCA':[0,{}],'GCC':[0,{}],'GCG':[0,{}],'GCT':[0,{}],'TCA':[0,{}],'TCC':[0,{}],'TCG':[0,{}],'TCT':[0,{}]}, 'plastid':{'AAA':[0,{}], 'AAC':[0,{}],'AAG':[0,{}], 'AAT':[0,{}],'CAA':[0,{}],'CAC':[0,{}],'CAG':[0,{}],'CAT':[0,{}],'GAA':[0,{}],'GAC':[0,{}],'GAG':[0,{}],'GAT':[0,{}],'TAA':[0,{}],'TAC':[0,{}],'TAG':[0,{}],'TAT':[0,{}],'ACA':[0,{}],'ACC':[0,{}],'ACG':[0,{}],'ACT':[0,{}],'CCA':[0,{}],'CCC':[0,{}],'CCG':[0,{}],'CCT':[0,{}],'GCA':[0,{}],'GCC':[0,{}],'GCG':[0,{}],'GCT':[0,{}],'TCA':[0,{}],'TCC':[0,{}],'TCG':[0,{}],'TCT':[0,{}]}  }
	cg_tri_dict_trans={'ACA':'TGT','ACC':'GGT','ACG':'CGT','ACT':'AGT','CCA':'TGG','CCC':'GGG','CCG':'CGG','CCT':'AGG','GCA':'TGC','GCC':'GGC','GCG':'CGC','GCT':'AGC','TCA':'TGA','TCC':'GGA','TCG':'CGA','TCT':'AGA','AGA':'TCT','AGC':'GCT','AGG':'CCT','AGT':'ACT','CGA':'TCG','CGC':'GCG','CGG':'CCG','CGT':'ACG','GGA':'TCC','GGC':'GCC','GGG':'CCC','GGT':'ACC','TGA':'TCA','TGC':'GCA','TGG':'CCA','TGT':'ACA'}
	at_sub=['AT>CG','AT>GC','AT>TA']
	cg_sub=['CG>AT','CG>GC','CG>TA']
	

	
	for chrom in tri_dict2:
		for tri in tri_dict2[chrom]:
			if tri[1] == 'A':
				for sub in at_sub:
					tri_dict2[chrom][tri][1][sub]=tri_dict2[chrom][tri][1].get(sub,[0,0])
			elif tri[1] == 'C':
				for sub in cg_sub:
					tri_dict2[chrom][tri][1][sub]=tri_dict2[chrom][tri][1].get(sub,[0,0])

	
	cd2={'A>C':'AT>CG','A>G':'AT>GC','A>T':'AT>TA','C>A':'CG>AT','C>G':'CG>GC','C>T':'CG>TA','T>G':'AT>CG','T>C':'AT>GC','T>A':'AT>TA','G>T':'CG>AT','G>C':'CG>GC','G>A':'CG>TA'}

	# trinuc_sum={}
	
	trinuc_sa={'mito': {'CDS': {'ACA': [0, 0], 'ACC': [0, 0], 'ACG': [0, 0], 'ACT': [0, 0], 'CCA': [0, 0], 'CCC': [0, 0], 'CCG': [0, 0], 'CCT': [0, 0], 'GCA': [0, 0], 'GCC': [0, 0], 'GCG': [0, 0], 'GCT': [0, 0], 'TCA': [0, 0], 'TCC': [0, 0], 'TCG': [0, 0], 'TCT': [0, 0], 'AGA': [0, 0], 'AGC': [0, 0], 'AGG': [0, 0], 'AGT': [0, 0], 'CGA': [0, 0], 'CGC': [0, 0], 'CGG': [0, 0], 'CGT': [0, 0], 'GGA': [0, 0], 'GGC': [0, 0], 'GGG': [0, 0], 'GGT': [0, 0], 'TGA': [0, 0], 'TGC': [0, 0], 'TGG': [0, 0], 'TGT': [0, 0]}, 'intron': {'ACA': [0, 0], 'ACC': [0, 0], 'ACG': [0, 0], 'ACT': [0, 0], 'CCA': [0, 0], 'CCC': [0, 0], 'CCG': [0, 0], 'CCT': [0, 0], 'GCA': [0, 0], 'GCC': [0, 0], 'GCG': [0, 0], 'GCT': [0, 0], 'TCA': [0, 0], 'TCC': [0, 0], 'TCG': [0, 0], 'TCT': [0, 0], 'AGA': [0, 0], 'AGC': [0, 0], 'AGG': [0, 0], 'AGT': [0, 0], 'CGA': [0, 0], 'CGC': [0, 0], 'CGG': [0, 0], 'CGT': [0, 0], 'GGA': [0, 0], 'GGC': [0, 0], 'GGG': [0, 0], 'GGT': [0, 0], 'TGA': [0, 0], 'TGC': [0, 0], 'TGG': [0, 0], 'TGT': [0, 0]}, 'tRNA': {'ACA': [0, 0], 'ACC': [0, 0], 'ACG': [0, 0], 'ACT': [0, 0], 'CCA': [0, 0], 'CCC': [0, 0], 'CCG': [0, 0], 'CCT': [0, 0], 'GCA': [0, 0], 'GCC': [0, 0], 'GCG': [0, 0], 'GCT': [0, 0], 'TCA': [0, 0], 'TCC': [0, 0], 'TCG': [0, 0], 'TCT': [0, 0], 'AGA': [0, 0], 'AGC': [0, 0], 'AGG': [0, 0], 'AGT': [0, 0], 'CGA': [0, 0], 'CGC': [0, 0], 'CGG': [0, 0], 'CGT': [0, 0], 'GGA': [0, 0], 'GGC': [0, 0], 'GGG': [0, 0], 'GGT': [0, 0], 'TGA': [0, 0], 'TGC': [0, 0], 'TGG': [0, 0], 'TGT': [0, 0]}, 'rRNA': {'ACA': [0, 0], 'ACC': [0, 0], 'ACG': [0, 0], 'ACT': [0, 0], 'CCA': [0, 0], 'CCC': [0, 0], 'CCG': [0, 0], 'CCT': [0, 0], 'GCA': [0, 0], 'GCC': [0, 0], 'GCG': [0, 0], 'GCT': [0, 0], 'TCA': [0, 0], 'TCC': [0, 0], 'TCG': [0, 0], 'TCT': [0, 0], 'AGA': [0, 0], 'AGC': [0, 0], 'AGG': [0, 0], 'AGT': [0, 0], 'CGA': [0, 0], 'CGC': [0, 0], 'CGG': [0, 0], 'CGT': [0, 0], 'GGA': [0, 0], 'GGC': [0, 0], 'GGG': [0, 0], 'GGT': [0, 0], 'TGA': [0, 0], 'TGC': [0, 0], 'TGG': [0, 0], 'TGT': [0, 0]}}, 'plastid': {'tRNA': {'ACA': [0, 0], 'ACC': [0, 0], 'ACG': [0, 0], 'ACT': [0, 0], 'CCA': [0, 0], 'CCC': [0, 0], 'CCG': [0, 0], 'CCT': [0, 0], 'GCA': [0, 0], 'GCC': [0, 0], 'GCG': [0, 0], 'GCT': [0, 0], 'TCA': [0, 0], 'TCC': [0, 0], 'TCG': [0, 0], 'TCT': [0, 0], 'AGA': [0, 0], 'AGC': [0, 0], 'AGG': [0, 0], 'AGT': [0, 0], 'CGA': [0, 0], 'CGC': [0, 0], 'CGG': [0, 0], 'CGT': [0, 0], 'GGA': [0, 0], 'GGC': [0, 0], 'GGG': [0, 0], 'GGT': [0, 0], 'TGA': [0, 0], 'TGC': [0, 0], 'TGG': [0, 0], 'TGT': [0, 0]}, 'CDS': {'ACA': [0, 0], 'ACC': [0, 0], 'ACG': [0, 0], 'ACT': [0, 0], 'CCA': [0, 0], 'CCC': [0, 0], 'CCG': [0, 0], 'CCT': [0, 0], 'GCA': [0, 0], 'GCC': [0, 0], 'GCG': [0, 0], 'GCT': [0, 0], 'TCA': [0, 0], 'TCC': [0, 0], 'TCG': [0, 0], 'TCT': [0, 0], 'AGA': [0, 0], 'AGC': [0, 0], 'AGG': [0, 0], 'AGT': [0, 0], 'CGA': [0, 0], 'CGC': [0, 0], 'CGG': [0, 0], 'CGT': [0, 0], 'GGA': [0, 0], 'GGC': [0, 0], 'GGG': [0, 0], 'GGT': [0, 0], 'TGA': [0, 0], 'TGC': [0, 0], 'TGG': [0, 0], 'TGT': [0, 0]}, 'intron': {'ACA': [0, 0], 'ACC': [0, 0], 'ACG': [0, 0], 'ACT': [0, 0], 'CCA': [0, 0], 'CCC': [0, 0], 'CCG': [0, 0], 'CCT': [0, 0], 'GCA': [0, 0], 'GCC': [0, 0], 'GCG': [0, 0], 'GCT': [0, 0], 'TCA': [0, 0], 'TCC': [0, 0], 'TCG': [0, 0], 'TCT': [0, 0], 'AGA': [0, 0], 'AGC': [0, 0], 'AGG': [0, 0], 'AGT': [0, 0], 'CGA': [0, 0], 'CGC': [0, 0], 'CGG': [0, 0], 'CGT': [0, 0], 'GGA': [0, 0], 'GGC': [0, 0], 'GGG': [0, 0], 'GGT': [0, 0], 'TGA': [0, 0], 'TGC': [0, 0], 'TGG': [0, 0], 'TGT': [0, 0]}, 'rRNA': {'ACA': [0, 0], 'ACC': [0, 0], 'ACG': [0, 0], 'ACT': [0, 0], 'CCA': [0, 0], 'CCC': [0, 0], 'CCG': [0, 0], 'CCT': [0, 0], 'GCA': [0, 0], 'GCC': [0, 0], 'GCG': [0, 0], 'GCT': [0, 0], 'TCA': [0, 0], 'TCC': [0, 0], 'TCG': [0, 0], 'TCT': [0, 0], 'AGA': [0, 0], 'AGC': [0, 0], 'AGG': [0, 0], 'AGT': [0, 0], 'CGA': [0, 0], 'CGC': [0, 0], 'CGG': [0, 0], 'CGT': [0, 0], 'GGA': [0, 0], 'GGC': [0, 0], 'GGG': [0, 0], 'GGT': [0, 0], 'TGA': [0, 0], 'TGC': [0, 0], 'TGG': [0, 0], 'TGT': [0, 0]}}}

	for chrom in trinuc_count:

	
		for setter in trinuc_count[chrom]:
			if len(trinuc_count[chrom][setter][0])==3:
				trinuc=trinuc_count[chrom][setter][0]
				if trinuc[1]=='C' or trinuc[1]=='G':
					if trinuc_count[chrom][setter][1] != 'intergenic':
						if trinuc_count[chrom][setter][6]== 1:
							trinuc_sa[chrom][trinuc_count[chrom][setter][1]][trinuc][0]+=trinuc_count[chrom][setter][3]
						else:
							trinuc_sa[chrom][trinuc_count[chrom][setter][1]][cg_tri_dict_trans[trinuc]][0]+=trinuc_count[chrom][setter][3]
						if (trinuc[1]=='C'and trinuc_count[chrom][setter][2]=='T') or (trinuc[1]=='G'and trinuc_count[chrom][setter][2]=='A'):

							if trinuc_count[chrom][setter][6]== '1':
								trinuc_sa[chrom][trinuc_count[chrom][setter][1]][trinuc][1]+=trinuc_count[chrom][setter][4]

							elif trinuc_count[chrom][setter][6]=='-1':
								trinuc_sa[chrom][trinuc_count[chrom][setter][1]][cg_tri_dict_trans[trinuc]][1]+=trinuc_count[chrom][setter][4]

						
				if trinuc[1]=='A' or trinuc[1]=='C': 
					tri_dict2[chrom][trinuc][0]+=trinuc_count[chrom][setter][3]
					if trinuc_count[chrom][setter][2]!='':
						
						tri_dict2[chrom][trinuc][1][cd2[trinuc[1]+'>'+trinuc_count[chrom][setter][2]]][0]+=trinuc_count[chrom][setter][4]
						tri_dict2[chrom][trinuc][1][cd2[trinuc[1]+'>'+trinuc_count[chrom][setter][2]]][1]+=trinuc_count[chrom][setter][5]

							
							
							
				else: 
					tri_dict2[chrom][tri_dict[trinuc]][0]+=trinuc_count[chrom][setter][3]
					if trinuc_count[chrom][setter][2]!='':
						tri_dict2[chrom][tri_dict[trinuc]][1][cd2[trinuc[1]+'>'+trinuc_count[chrom][setter][2]]][0]+=trinuc_count[chrom][setter][4]
						tri_dict2[chrom][tri_dict[trinuc]][1][cd2[trinuc[1]+'>'+trinuc_count[chrom][setter][2]]][1]+=trinuc_count[chrom][setter][5]

	
	for chrom in trinuc_sa:
		for region in trinuc_sa[chrom]:
			for tri in trinuc_sa[chrom][region]:
				if tri[1]=='C':
					out_handle_trinuc_strnd_aware.write(rep_ID +'\t' +rep_ID_list[0] +'\t' +rep_ID_list[1] +'\t'+rep_ID_list[2]+'\t' +rep_ID_list[3]+ '\t'+ '\t' + chrom + '\t' + region + '\t' + tri + '\t' + 'C' +'\t' + tri  + '\t' + str(trinuc_sa[chrom][region][tri][0])+ '\t' + str(trinuc_sa[chrom][region][tri][1]) +'\n' )
				else:
					out_handle_trinuc_strnd_aware.write(rep_ID +'\t' +rep_ID_list[0] +'\t' +rep_ID_list[1] +'\t'+rep_ID_list[2]+'\t' +rep_ID_list[3]+ '\t'+ '\t' + chrom + '\t' + region + '\t' + tri_dict[tri] + '\t' + 'G'+ '\t' + tri  + '\t' + str(trinuc_sa[chrom][region][tri][0])+ '\t' + str(trinuc_sa[chrom][region][tri][1])+'\n' )
	
	for chrom in tri_dict2:
		for trinuc in tri_dict2[chrom]:
			for sub in tri_dict2[chrom][trinuc][1]:
				out_handle_trinuc.write(rep_ID +'\t' +rep_ID_list[0] +'\t' +rep_ID_list[1] +'\t'+rep_ID_list[2]+'\t' +rep_ID_list[3]+ '\t' + chrom + '\t' + trinuc+ '\t' + str(tri_dict2[chrom][trinuc][0]) +'\t'     + sub + '\t' + str(tri_dict2[chrom][trinuc][1][sub][0]) + '\t' + str(tri_dict2[chrom][trinuc][1][sub][1]) +'\n'   )
						

	for chrom in SNV_hit_dict:
		for pos in SNV_hit_dict[chrom]:
			for alt_allele in SNV_hit_dict[chrom][pos]:
				out_handle_snv.write(rep_ID +'\t' +rep_ID_list[0] +'\t' +rep_ID_list[1] +'\t' +rep_ID_list[2] +'\t' +rep_ID_list[3] +'\t' + chrom + '\t' + str(pos) + '\t' + str(SNV_hit_dict[chrom][pos][alt_allele][0]) +'\t' + alt_allele +'\t' + str(SNV_hit_dict[chrom][pos][alt_allele][1]) +'\t' + str(seq_dict[chrom][pos][2]) +'\n')
				
				
				
	for chrom in di_nuc_hit_dict:
		for pos in di_nuc_hit_dict[chrom]:
			for alt_allele in di_nuc_hit_dict[chrom][pos]:
				out_handle_comp.write(rep_ID +'\t' +rep_ID_list[0] +'\t' +rep_ID_list[1] +'\t' +rep_ID_list[2] +'\t' +rep_ID_list[3] +'\t' + chrom + '\t' + str(pos) + '\t'  + str(di_nuc_hit_dict[chrom][pos][alt_allele][0]) +'\t' + alt_allele +'\t' + str(di_nuc_hit_dict[chrom][pos][alt_allele][1]) +'\t' + str(seq_dict[chrom][pos][2]) +'\n' )
				
				
	
				
	for chrom in cov_count:
		for nt in cov_count[chrom]:
			for snv in cov_count[chrom][nt][1]:
				out_handle_total_cov_and_snvs.write(rep_ID +'\t' +rep_ID_list[0] +'\t' +rep_ID_list[1] +'\t' +rep_ID_list[2] +'\t' +rep_ID_list[3] +'\t' + chrom + '\t' + nt  + '\t' + snv + '\t' + str(cov_count[chrom][nt][0]) + '\t' + str(cov_count[chrom][nt][1][snv][0]) + '\t' + str(cov_count[chrom][nt][1][snv][1]) + '\t' + str(cov_count[chrom][nt][1][snv][0]/cov_count[chrom][nt][0])+ '\t' + str(cov_count[chrom][nt][1][snv][1]/cov_count[chrom][nt][0]) + '\n'   )
				
	for chrom in total_cov:
		if chrom in di_nuc_hit_dict:
			out_handle_total_cov_and_dinucs.write(rep_ID +'\t' +rep_ID_list[0] +'\t' +rep_ID_list[1] +'\t' +rep_ID_list[2] +'\t' +rep_ID_list[3] +'\t' + chrom + '\t' + str(len(di_nuc_hit_dict[chrom])) + '\t'  + str(total_cov[chrom]) + '\t'  + str(len(di_nuc_hit_dict[chrom])/total_cov[chrom]) +'\n'  )
		else:
			out_handle_total_cov_and_dinucs.write(rep_ID +'\t' +rep_ID_list[0] +'\t' +rep_ID_list[1] +'\t' +rep_ID_list[2] +'\t' +rep_ID_list[3] +'\t' + chrom + '\t' + str(0) + '\t'  + str(total_cov[chrom]) + '\t'  + str(0) +'\n'  )
		
				
	for chrom in region_count:
		for region in region_count[chrom]:
			for nt in region_count[chrom][region]:
				for snv in region_count[chrom][region][nt][1]:
					out_handle_region.write(rep_ID +'\t' +rep_ID_list[0] +'\t' +rep_ID_list[1] +'\t' +rep_ID_list[2] +'\t' +rep_ID_list[3] +'\t' + chrom + '\t' + region  + '\t' + nt +'\t'   + snv + '\t' + str(region_count[chrom][region][nt][0]) +'\t' + str(region_count[chrom][region][nt][1][snv][0]) +'\t' + str(region_count[chrom][region][nt][1][snv][1])  +'\t' + str(region_count[chrom][region][nt][1][snv][0]/region_count[chrom][region][nt][0])  +'\t' + str(region_count[chrom][region][nt][1][snv][1]/region_count[chrom][region][nt][0]) +      '\n'  )
					
	strand_corrected={}
	comp_dict={'A':'T','C':'G','G':'C','T':'A'}
	
	SNV_comp_dict={'A>C':'T>G','A>G':'T>C','A>T':'T>A','C>A':'G>T','C>G':'G>C','C>T':'G>A','G>A':'C>T','G>C':'C>G','G>T':'C>A','T>A':'A>T','T>C':'A>G','T>G':'A>C'}
	for chrom in strand_count:
		for region in strand_count[chrom]:
			if region != 'intergenic':
				for nt in strand_count[chrom][region]:
					for strand in strand_count[chrom][region][nt]:
						snv_count=0
						for snv in strand_count[chrom][region][nt][strand][1]:
							snv_count+=1
							if strand == '1':
								strand_corrected[chrom]=strand_corrected.get(chrom,{})
								strand_corrected[chrom][region]=strand_corrected[chrom].get(region,{})
								strand_corrected[chrom][region][nt]=strand_corrected[chrom][region].get(nt,[0,{}])
								if snv_count<2:
									strand_corrected[chrom][region][nt][0]+=strand_count[chrom][region][nt][strand][0]
								strand_corrected[chrom][region][nt][1][snv]=strand_corrected[chrom][region][nt][1].get(snv,[0,0])
								strand_corrected[chrom][region][nt][1][snv][0]+=strand_count[chrom][region][nt][strand][1][snv][0]
								strand_corrected[chrom][region][nt][1][snv][1]+=strand_count[chrom][region][nt][strand][1][snv][1]
							else:
								strand_corrected[chrom]=strand_corrected.get(chrom,{})
								strand_corrected[chrom][region]=strand_corrected[chrom].get(region,{})
								strand_corrected[chrom][region][comp_dict[nt]]= strand_corrected[chrom][region].get(comp_dict[nt],[0,{}])
								if snv_count<2:
									strand_corrected[chrom][region][comp_dict[nt]][0]+=strand_count[chrom][region][nt][strand][0]
								strand_corrected[chrom][region][comp_dict[nt]][1][SNV_comp_dict[snv]]=strand_corrected[chrom][region][comp_dict[nt]][1].get(SNV_comp_dict[snv],[0,0])	
								strand_corrected[chrom][region][comp_dict[nt]][1][SNV_comp_dict[snv]][0]+=strand_count[chrom][region][nt][strand][1][snv][0]
								strand_corrected[chrom][region][comp_dict[nt]][1][SNV_comp_dict[snv]][1]+=strand_count[chrom][region][nt][strand][1][snv][1]
					
		
				
	for chrom in strand_corrected:
		for region in strand_corrected[chrom]:
			for nt in strand_corrected[chrom][region]:
				for snv in strand_corrected[chrom][region][nt][1]:
					out_handle_strand.write(rep_ID +'\t' +rep_ID_list[0] +'\t' +rep_ID_list[1] +'\t' +rep_ID_list[2] +'\t' +rep_ID_list[3] +'\t' + chrom +'\t' + region +'\t' +nt + '\t' + snv +'\t' + str(strand_corrected[chrom][region][nt][0]) + '\t' + str(strand_corrected[chrom][region][nt][1][snv][0]) + '\t' + str(strand_corrected[chrom][region][nt][1][snv][1])  + '\t' + str(strand_corrected[chrom][region][nt][1][snv][0]/strand_corrected[chrom][region][nt][0])+ '\t' + str(strand_corrected[chrom][region][nt][1][snv][1]/strand_corrected[chrom][region][nt][0]) +'\n'  )
								
	for chrom in seq_dict:
		if chrom != 'ignore':
			for line in seq_dict[chrom]:
				liner=''
				liner+=chrom
				liner+='\t'
				count=0
				for col in line:
					count+=1
					if count < 8:
						liner += str(col)
						liner +='\t'
					elif count == 8:
						liner +=str(col)
						liner +='\n'
				
						
				out_handle_depth.write(liner)
								
	

	
if __name__=='__main__':
	print(complex_snv())	
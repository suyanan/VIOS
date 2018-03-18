#-*- coding:utf-8 -*-

#from .config_paras import *

import os
import django

from .config_paras import *
from .show_results_blastn import *
from matplotlib.pyplot import *
from collections import Counter

if __name__ == '__main__':
    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    os.environ['DJANGO_SETTINGS_MODULE'] = 'VIOS.settings'
    django.setup()


def files_reads_conter(querySample):
	"""for fig1-3 data sources"""
	#raw fastq file
	reads_raw_datas_file = os.path.join(trimmomatic_path,querySample+'_trimlog.txt')
	
	#after AC, four fastq file
	reads_after_qc_P_file = os.path.join(trimmomatic_path,querySample+'_output_forward_paired.fastq')
	reads_after_qc_FU_file = os.path.join(trimmomatic_path,querySample+'_output_forward_unpaired.fastq')
	reads_after_qc_RU_file = os.path.join(trimmomatic_path,querySample+'_output_reverse_unpaired.fastq')
	reads_after_qc_se_output_file = os.path.join(trimmomatic_path,querySample+'_output.fastq')

	#after map GH_HG fastq file
	reads_after_map_hg38_P_file = os.path.join(bowtie2_path,querySample+'_pair_host.1.fastq')
	reads_after_map_hg38_U_file = os.path.join(bowtie2_path,querySample+map_gb_hg_singleAndunpair_fastq_file_suffix)

	#after map NT_VIRUSE fastq file(first map aginst Virus,can use the viruses_SAM file as the fourth figure file)
	reads_after_map_viruses_P_file = os.path.join(bowtie2_path,querySample+'_pair_VIRUS.1.fastq')
	reads_after_map_viruses_U_file = os.path.join(bowtie2_path,querySample+map_nt_viruses_singleAndunpair_fastq_file_suffix)
	reads_after_map_viruses_SAM_file = os.path.join(bowtie2_path,querySample+map_nt_viruses_singleAndunpair_sam_file_suffix)
	
	reads_counter_file = os.path.join(result_total_sta_path,'pipeline_reads_sta.txt')
	if os.path.exists(reads_after_qc_P_file):
		os.system('wc -l %s %s %s %s %s %s %s %s %s > %s'
			%(reads_raw_datas_file,
				reads_after_qc_P_file,reads_after_qc_FU_file,reads_after_qc_RU_file,
				reads_after_map_hg38_P_file,reads_after_map_hg38_U_file,
				reads_after_map_viruses_P_file,reads_after_map_viruses_U_file,reads_after_map_viruses_SAM_file,
				reads_counter_file))
		with open(reads_counter_file,'r')  as reads_counter_file_obj:
			reads_counter_lines_list = reads_counter_file_obj.readlines()
			reads_raw_datas_numbers = int(reads_counter_lines_list[0].split(' ')[0:-1][-1])
			reads_after_qc_numbers = (int(reads_counter_lines_list[1].split(' ')[0:-1][-1])*2+int(reads_counter_lines_list[2].split(' ')[0:-1][-1])+int(reads_counter_lines_list[3].split(' ')[0:-1][-1]))/4
			reads_after_map_hg38_numbers = (int(reads_counter_lines_list[4].split(' ')[0:-1][-1])*2+int(reads_counter_lines_list[5].split(' ')[0:-1][-1]))/4
			reads_after_map_viruses_numbers = (int(reads_counter_lines_list[6].split(' ')[0:-1][-1])*2+int(reads_counter_lines_list[7].split(' ')[0:-1][-1]))/4
			reads_after_map_viruses_SAM_numbers = int(reads_counter_lines_list[8].split(' ')[0:-1][-1])##not equal to reads_after_map_viruses_numbers
			return reads_raw_datas_numbers, reads_after_qc_numbers, reads_after_map_hg38_numbers, reads_after_map_viruses_numbers,reads_after_map_viruses_SAM_numbers
	else:
		os.system('wc -l %s %s %s %s %s > %s'
			%(reads_raw_datas_file,
				reads_after_qc_se_output_file,
				reads_after_map_hg38_U_file,
				reads_after_map_viruses_U_file,reads_after_map_viruses_SAM_file,
				reads_counter_file))
		with open(reads_counter_file,'r') as reads_counter_file_obj:
			reads_counter_lines_list = reads_counter_file_obj.readlines()
			reads_raw_datas_numbers = int(reads_counter_lines_list[0].split(' ')[0:-1][-1])
			reads_after_qc_numbers = int(reads_counter_lines_list[1].split(' ')[0:-1][-1])/4
			reads_after_map_hg38_numbers = int(reads_counter_lines_list[2].split(' ')[0:-1][-1])/4
			reads_after_map_viruses_numbers = int(reads_counter_lines_list[3].split(' ')[0:-1][-1])/4
			reads_after_map_viruses_SAM_numbers = int(reads_counter_lines_list[4].split(' ')[0:-1][-1])##not equal to reads_after_map_viruses_numbers
			return reads_raw_datas_numbers, reads_after_qc_numbers, reads_after_map_hg38_numbers, reads_after_map_viruses_numbers,reads_after_map_viruses_SAM_numbers

		

def nt_filter_target_accessionLength(fromBigFile,queryTheAccession):
    '''
    too slow ,need to modify
    '''
    with open(fromBigFile, 'r') as fromBigFile_obj:
        lines_list = fromBigFile_obj.readlines()
        lines_list_len = len(lines_list)
        keys_accession_list = []
        keys_accession_index_list = []
        for i in range(len(lines_list)):
            if lines_list[i][0] == '>':
                keys_accession_list.append(lines_list[i].strip('\n'))
                keys_accession_index_list.append(i)
        keys_accession_index_list.append(lines_list_len)

        values_seq_list = []
        for i in range(1, len(keys_accession_index_list)):
            values_seq_one_list = []
            for j in range(keys_accession_index_list[i - 1] + 1, keys_accession_index_list[i]):
                values_seq_one = lines_list[j].strip('\n')
                values_seq_one_list.append(values_seq_one)                
            values_seq_list.append(values_seq_one_list)  # before values_seq_one_list = [], append
 

        for i in range(len(keys_accession_list)):
        	accession = keys_accession_list[i].split(' ')[0].split('>')[-1]
        	if accession == queryTheAccession:
        		accession_length_sum = []
        		for j in range(len(values_seq_list[i])):
        			one_line_length = len(values_seq_list[i][j])
        			accession_length_sum.append(one_line_length)
        		return sum(accession_length_sum)


def pipeline_reads_sta_to_grafic(sequenceMethod,sampleListResults,querySample,ntVirusesFamilyName):

	fig = figure(figsize=(15,5),dpi=80)
	p2 = subplot(1,1,1)

	for i in range(len(sampleListResults)):
		if sampleListResults[i] == querySample and sequenceMethod == 'pe':
			reads_after_map_viruses_P_file = os.path.join(bowtie2_path,querySample+'_pair_VIRUS.1.fastq')
			reads_after_map_viruses_U_file = os.path.join(bowtie2_path,querySample+map_nt_viruses_singleAndunpair_fastq_file_suffix)

			forward_file = samplesRawDataForwardFiles[i]
			with open(forward_file,'r') as forward_file_obj:
				reads_raw_datas_numbers = 2*(len(forward_file_obj.readlines()))/4
			with open(reads_after_map_viruses_P_file,'r') as reads_after_qc_P_file_obj:
				reads_after_map_viruses_pair_numbers = 2*(len(reads_after_qc_P_file_obj.readlines()))/4
			with open(reads_after_map_viruses_U_file,'r') as reads_after_map_viruses_U_file_obj:
				reads_after_map_viruses_unpair_numbers = (len(reads_after_map_viruses_U_file_obj.readlines()))/4
				reads_after_map_viruses_numbers = reads_after_map_viruses_pair_numbers+reads_after_map_viruses_unpair_numbers
				print reads_raw_datas_numbers
				print reads_after_map_viruses_numbers

				reads_map_viruses_list = [reads_after_map_viruses_numbers,reads_raw_datas_numbers-reads_after_map_viruses_numbers]
				print reads_map_viruses_list
				label_reads_map_viruses_list = ['viruses','others']
				explode = [0.05,0]
				p2.pie(reads_map_viruses_list,explode=explode,shadow=True,labels=label_reads_map_viruses_list,startangle=0,labeldistance=1.1,pctdistance=0.6)
				p2.axis('equal')		
		
				#subplots_adjust(left=-0.1,right=1.1,wspace=-0.6)
				title_str = "the reads change of %s smaple after mapping"%(querySample)
				text(-1.5,1.2,title_str,fontsize=14)  #1.larger,righter; 2.larger,upper

				fig_file = os.path.join(result_total_sta_path,'sta_reads_'+querySample+'.png')
				fig.savefig(fig_file)
				show()

		if sampleListResults[i] == querySample and sequenceMethod == 'se':
			reads_after_map_viruses_U_file = os.path.join(bowtie2_path,querySample+map_nt_viruses_singleAndunpair_fastq_file_suffix)

			single_file = samplesRawDataSingleFiles[i]

			with open(single_file,'r') as single_file_obj:
				reads_raw_datas_numbers = len(single_file_obj.readlines())/4
			with open(reads_after_map_viruses_U_file,'r') as reads_after_map_viruses_U_file_obj:
				reads_after_map_viruses_unpair_numbers = (len(reads_after_map_viruses_U_file_obj.readlines()))/4
				reads_after_map_viruses_numbers = reads_after_map_viruses_unpair_numbers
				print reads_raw_datas_numbers
				print reads_after_map_viruses_numbers

				reads_map_viruses_list = [reads_after_map_viruses_numbers,reads_raw_datas_numbers-reads_after_map_viruses_numbers]
				label_reads_map_viruses_list = ['viruses','others']
				explode = [0.05,0]
				p2.pie(reads_map_viruses_list,explode=explode,shadow=True,labels=label_reads_map_viruses_list,startangle=0,labeldistance=1.1,pctdistance=0.6)
				p2.axis('equal')		

				#subplots_adjust(left=-0.1,right=1.1,wspace=-0.6)
				title_str = "the reads change of %s smaple after mapping"%(querySample)
				text(-1.5,1.2,title_str,fontsize=14)  #1.larger,righter; 2.larger,upper

				fig_file = os.path.join(result_total_sta_path,'sta_reads_'+querySample+'.png')
				fig.savefig(fig_file)
				show()

def pipeline_identifer_fig(querySample,ntVirusesFamilyName):
	fig = figure(figsize=(15,5),dpi=80)

	##fig4: reads statistics about top 10 virus in SAM file.
	reads_after_map_viruses_SAM_file = os.path.join(bowtie2_path,querySample+map_nt_viruses_singleAndunpair_sam_file_suffix)
	with open(reads_after_map_viruses_SAM_file,'r') as reads_after_map_viruses_SAM_file_obj:
		reads_lines = reads_after_map_viruses_SAM_file_obj.readlines()
		reads_after_map_viruses_SAM_numbers = len(reads_lines)
		reads_fig4_data = []
		for i in range(len(reads_lines)):
			each_accesion = reads_lines[i].split('\t')[2]
			reads_fig4_data.append(each_accesion)
		#print reads_fig4_data  ##[..., 'AB256208.1', 'AB256497.1', 'AB256497.1', '*']

		counter = Counter(reads_fig4_data)
		#print counter  ##Counter({'KF429752.1': 13, ..., 'FJ943620.1': 1})	#<class 'collections.Counter'>
		counter_most = counter.most_common(setShowTopLines-5)  ##!!!get most common(improper!!!), not map aginst bacteria, results may include bacteria.
		#print counter_most  ##[('KF429752.1', 13), ('KF268195.1', 12), ('KF429748.1', 9), ('JN860680.1', 8), ('KF268315.1', 8), ('AY601636.1', 8), ('JX423386.1', 8), ('HC084950.1', 7), ('KF268120.1', 7), ('FJ025928.1', 7)]	<type 'list'>
		counter_most_dict = dict(counter_most)
		#print counter_most_dict  ##{'KF268120.1': 7, 'KF429748.1': 9, 'KF429752.1': 13, 'FJ025928.1': 7, 'KF268315.1': 8, 'JN860680.1': 8, 'HC084950.1': 7, 'KF268195.1': 12, 'AY601636.1': 8, 'JX423386.1': 8}
		
		counter_most_dict_sorted_list = sorted(counter_most_dict.items(),key=lambda item:item[1],reverse = True) #[('KM209255.1', 5868), ('KM209277.1', 1973), 

		accession_data = []
		counter_data = []
		for i in range(len(counter_most_dict_sorted_list)):
			accession_data.append(counter_most_dict_sorted_list[i][0])
			counter_data.append(counter_most_dict_sorted_list[i][1])
			#print accession_data  #just counter the 3rd colomn, not 7th colomn(mayby include the accession string)

		##change accesion to annotation
		annotation_data = []
		for i in range(len(accession_data)):
			annotation = search_nt_title(ntVirusesFamilyName,accession_data[i])
			annotation_data.append(annotation)

		annotation_data.append('other viruses')
		counter_data.append(reads_after_map_viruses_SAM_numbers-sum(counter_data))	
		#print annotation_data  #['CP002121.1', 'KP279748.1', 'CP016318.1', 'KX494979.1', 'AM420293.1', 'FN997652.1', 'CP008698.1', 'CP000233.1', 'AB256208.1', 'KY065497.1', 'others viruses']
		#print counter_data  #[10335, 3922, 195, 3715, 432, 291, 264, 346, 317, 1619, 1599]

		#p4 = subplot(2,1,2)  #the before 3 figs was (2,1,1)
		p4 = subplot(1,1,1)  #the before 3 figs was (2,1,1)
		explode = [0,0,0,0,0,0.1]
		#p4.pie(counter_data,shadow=True,labels=annotation_data,autopct='%2.1f%%',startangle=0,labeldistance=1.1,pctdistance=0.6)
		p4.pie(counter_data,explode=explode,shadow=True,autopct='%2.1f%%',startangle=0,labeldistance=1.1,pctdistance=0.6)
		p4.axis('equal')
		subplots_adjust(left=-0.1,right=1.1,wspace=-0.6)
		#legend(annotation_data)
		legend(annotation_data,loc='bottom,left',bbox_to_anchor=(0.6,0.40))	#1value larger righter,2value largerer upper
		#legend(annotation_data,loc='bottom,left',bbox_to_anchor=(0.53,0.07))
		#subplots_adjust(bottom=0.2)
		#tight_layout(h_pad=2)
		title_str = "%s reads changes while processing"%(querySample)
		#text(1.6,-1.4,title_str,fontsize=14)
		text(-1.0,1.1,title_str,fontsize=14)  #0,0 center, 1value larger righter, 2value larger upper

		fig_file = os.path.join(result_total_sta_path,'sta_identifer_'+querySample+'.png')
		fig.savefig(fig_file)
		show()


def pipeline_identifer_list(querySample,ntVirusesFamilyName):
	##statistics:reads_after_map_viruses_SAM_numbers
	#reads_raw_datas_numbers, reads_after_qc_numbers, reads_after_map_hg38_numbers, reads_after_map_viruses_numbers, reads_after_map_viruses_SAM_numbers = files_reads_conter(querySample)

	reads_after_map_viruses_SAM_file = os.path.join(bowtie2_path,querySample+map_nt_viruses_singleAndunpair_sam_file_suffix)

	with open(reads_after_map_viruses_SAM_file,'r') as reads_after_map_viruses_SAM_file_obj:
		reads_lines = reads_after_map_viruses_SAM_file_obj.readlines()
		reads_after_map_viruses_SAM_numbers = len(reads_lines)
		reads_fig4_data = []
		for i in range(len(reads_lines)):
			each_accesion = reads_lines[i].split('\t')[2]
			reads_fig4_data.append(each_accesion)
		#print reads_fig4_data  ##[..., 'AB256208.1', 'AB256497.1', 'AB256497.1', '*']
	

		counter = Counter(reads_fig4_data)
		#print counter  ##Counter({'KF429752.1': 13, ..., 'FJ943620.1': 1})	#<class 'collections.Counter'>
		
		#counter_most = counter.most_common(setShowTopLines+20)  ##!!!get most common(improper!!!), not map aginst bacteria, results may include bacteria.
		counter_most = counter.most_common(5) 
		#counter_most = counter.most_common(100) 
		#print counter_most  ##[('KF429752.1', 13), ('KF268195.1', 12), ('KF429748.1', 9), ('JN860680.1', 8), ('KF268315.1', 8), ('AY601636.1', 8), ('JX423386.1', 8), ('HC084950.1', 7), ('KF268120.1', 7), ('FJ025928.1', 7)]	<type 'list'>
		counter_most_dict = dict(counter_most)
		#print counter_most_dict  ##{'KF268120.1': 7, 'KF429748.1': 9, 'KF429752.1': 13, 'FJ025928.1': 7, 'KF268315.1': 8, 'JN860680.1': 8, 'HC084950.1': 7, 'KF268195.1': 12, 'AY601636.1': 8, 'JX423386.1': 8}
		
		counter_most_dict_sorted_list = sorted(counter_most_dict.items(),key=lambda item:item[1],reverse = True) #[('KM209255.1', 5868), ('KM209277.1', 1973), 

		##TABLE: 4 COLOMN (accession annotation number ratio)
		accession_data = []
		counter_data = []
		for i in range(len(counter_most_dict_sorted_list)):
			accession_data.append(counter_most_dict_sorted_list[i][0])
			counter_data.append(counter_most_dict_sorted_list[i][1])
			#print accession_data  #just counter the 3rd colomn, not 7th colomn(mayby include the accession string)
		
		##change accesion to annotation
		annotation_data = []
		for i in range(len(accession_data)):
			annotation = search_nt_title(ntVirusesFamilyName,accession_data[i])
			annotation_data.append(annotation)

		'''
		##resutlt.html show_results_coverage method can't match the "the left" like accession
		accession_data.append('the left')
		annotation_data.append('left viruses')
		counter_data.append(reads_after_map_viruses_SAM_numbers-sum(counter_data))	
		#print annotation_data  #['CP002121.1', 'KP279748.1', 'CP016318.1', 'KX494979.1', 'AM420293.1', 'FN997652.1', 'CP008698.1', 'CP000233.1', 'AB256208.1', 'KY065497.1', 'others viruses']
		#print counter_data  #[10335, 3922, 195, 3715, 432, 291, 264, 346, 317, 1619, 1599]
		'''
		counterRatio_data = []
		for i in range(len(counter_data)):
			eachRatio = format(float(counter_data[i])/float(reads_after_map_viruses_SAM_numbers),'.4f')
			counterRatio_data.append(eachRatio)

		identifer_data2D = []
		for i in range(len(counter_data)):
			each_line = [accession_data[i],annotation_data[i],counter_data[i],counterRatio_data[i]]
			identifer_data2D.append(each_line)
		return identifer_data2D


def eachAccessionCoverageFig_fromSAMFile(querySample,ntVirusesFamilyName,queryTheAccession):
	"""SAM File to target the accession from results.html page. EXAMPLE:AB289986.1 """
	##time_clock:1.01194   time_time:2.02953600883
	reads_after_map_viruses_file = os.path.join(bowtie2_path,querySample+map_nt_viruses_singleAndunpair_sam_file_suffix)
	
	family_path = 'nt_%s'%(ntVirusesFamilyName)
	family_file_name = 'nt_%s.fasta'%(ntVirusesFamilyName)

	nt_viruses_family_path = blastDB_path+family_path
	nt_viruses_family_file = os.path.join(nt_viruses_family_path,family_file_name)

	nt_viruses_family_seqLengthList_file_name = 'nt_%s_seqLengthList'%(ntVirusesFamilyName)
	nt_viruses_family_seqLengthList_file = os.path.join(nt_viruses_family_path,nt_viruses_family_seqLengthList_file_name)

	accession_length = 0
	with open(nt_viruses_family_seqLengthList_file,'r') as nt_viruses_family_seqLengthList_file_obj:
		linesList = nt_viruses_family_seqLengthList_file_obj.readlines()
		for i in range(len(linesList)):
			eachLine = linesList[i]
			if eachLine.strip('\n').split('\t')[0] == queryTheAccession:
				accession_length = int(eachLine.strip('\n').split('\t')[1])
	print accession_length
	#accession_length = nt_filter_target_accessionLength(nt_viruses_family_file,queryTheAccession)
	baseNum_atPos_list = [0 for i in range(accession_length)]  #[0,0,...,0]  3141 bases
	print len(baseNum_atPos_list)
	
	with open(reads_after_map_viruses_file,'r') as reads_after_map_viruses_file_obj:
		reads_lines = reads_after_map_viruses_file_obj.readlines()
		reads_posStart_list = [] #the 4th column		
		reads_posEnd_list = []  #each_posStart+each_cigar+1   ##reads_cigar_list = []  #the 6th column
		for i in range(len(reads_lines)):
			each_accession = reads_lines[i].split('\t')[2]
			if each_accession == queryTheAccession:
				each_posStart = int(reads_lines[i].split('\t')[3])
				reads_posStart_list.append(each_posStart)
				#each_cigar = int(reads_lines[i].split('\t')[5])  #simplify to just number 65, not 65M!!!at beginning
				#each_posEnd = each_posStart+each_cigar  #[60, 52, 56, 33, 56, 64, 89]
				##EXAMAPLE: '8M8I135M'  MIDNSHP=X nine cigar character with number before.
				each_cigar = reads_lines[i].split('\t')[5]
				cigarPattern = re.compile(r'\d+|M|I|D|N|S|H|P|=|X')
				cigar_list = cigarPattern.findall(each_cigar)  #['8','M','8','I','135','M']
				cigar_ref_list = []
				for j in range(len(cigar_list)):
					if cigar_list[j] == 'M' or cigar_list[j] == 'D' or cigar_list[j] == 'N' or cigar_list[j] == '=' or cigar_list[j] == 'X':
						cigar_ref_list.append(int(cigar_list[j-1]))
				cigar_ref_counter = sum(cigar_ref_list)
				each_posEnd = each_posStart+cigar_ref_counter
				reads_posEnd_list.append(each_posEnd)
		#print reads_posStart_list  #[2800, 2834, 2730, 2823, 2757, 2848, 2777]
		#print reads_posEnd_list    #[2860, 2886, 2786, 2856, 2813, 2912, 2866]  ##2800/2801/.../2869 sixty bases were seqed!2860 not!
		for i in range(len(reads_posStart_list)):
			pos_start = reads_posStart_list[i]
			pos_end = reads_posEnd_list[i]
			#for j in range(pos_start-1,pos_end):  #[pos_start,pos_end]  ##IndexError: list index out of range
			for j in range(pos_start-1,min(pos_end,accession_length)):
				baseNum_atPos_list[j] += 1  #elemnts number is accession_length, but maybe beyond the reference(0,accession_length), sequenced base exists.
	print len(baseNum_atPos_list)

	fig,ax = subplots(figsize=(10,5))
	x = range(1,accession_length+1)
	y = baseNum_atPos_list

	plot(x,y)
	fill_between(x,0,y)

	title_str = "%s base coverage with mappting aginst %s"%(querySample,queryTheAccession)
	title(title_str,fontsize=14)
	xlabel("base position")
	ylabel("coverage")
	xlim(1,accession_length)
	ylim(ymin=0)
	grid(color='lightskyblue')

	output_path = bowtie2_path_basecoverage + querySample
	if not os.path.exists(output_path):
	    os.system('mkdir %s'%output_path)

	fig_file = os.path.join(output_path,queryTheAccession+'.png')
	fig.savefig(fig_file)
	show()

def eachAccessionCoverageFig_fromSAMFile1(querySample,ntVirusesFamilyName,queryTheAccession):
	"""SAM File to target the accession from results.html page. EXAMPLE:AB289986.1 """
	##time_clock:17.894935  time_time:22.2830429077
	reads_after_map_viruses_file = os.path.join(bowtie2_path,querySample+map_nt_viruses_singleAndunpair_sam_file_suffix)
	
	family_path = 'nt_%s'%(ntVirusesFamilyName)
	family_file_name = 'nt_%s.fasta'%(ntVirusesFamilyName)

	nt_viruses_family_path = blastDB_path+family_path
	nt_viruses_family_file = os.path.join(nt_viruses_family_path,family_file_name)

	nt_viruses_family_seqLengthList_file_name = 'nt_%s_seqLengthList'%(ntVirusesFamilyName)
	nt_viruses_family_seqLengthList_file = os.path.join(nt_viruses_family_path,nt_viruses_family_seqLengthList_file_name)

	accession_length = 0
	with open(nt_viruses_family_seqLengthList_file,'r') as nt_viruses_family_seqLengthList_file_obj:
		linesList = nt_viruses_family_seqLengthList_file_obj.readlines()
		for i in range(len(linesList)):
			eachLine = linesList[i]
			if eachLine.strip('\n').split('\t')[0] == queryTheAccession:
				accession_length = int(eachLine.strip('\n').split('\t')[1])
	print accession_length
	#accession_length = nt_filter_target_accessionLength(nt_viruses_family_file,queryTheAccession)
	
	with open(reads_after_map_viruses_file,'r') as reads_after_map_viruses_file_obj:
		reads_lines = reads_after_map_viruses_file_obj.readlines()
		reads_posStart_list = [] #the 4th column		
		reads_posEnd_list = []  #each_posStart+each_cigar+1   ##reads_cigar_list = []  #the 6th column
		for i in range(len(reads_lines)):
			each_accession = reads_lines[i].split('\t')[2]
			if each_accession == queryTheAccession:
				each_posStart = int(reads_lines[i].split('\t')[3])
				reads_posStart_list.append(each_posStart)
				#each_cigar = int(reads_lines[i].split('\t')[5])  #simplify to just number 65, not 65M!!!at beginning
				#each_posEnd = each_posStart+each_cigar  #[60, 52, 56, 33, 56, 64, 89]
				##EXAMAPLE: '8M8I135M'  MIDNSHP=X nine cigar character with number before.
				each_cigar = reads_lines[i].split('\t')[5]
				cigarPattern = re.compile(r'\d+|M|I|D|N|S|H|P|=|X')
				cigar_list = cigarPattern.findall(each_cigar)  #['8','M','8','I','135','M']
				cigar_ref_list = []
				for j in range(len(cigar_list)):
					if cigar_list[j] == 'M' or cigar_list[j] == 'D' or cigar_list[j] == 'N' or cigar_list[j] == '=' or cigar_list[j] == 'X':
						cigar_ref_list.append(int(cigar_list[j-1]))
				cigar_ref_counter = sum(cigar_ref_list)
				each_posEnd = each_posStart+cigar_ref_counter
				reads_posEnd_list.append(each_posEnd)
		#print reads_posStart_list  #[2800, 2834, 2730, 2823, 2757, 2848, 2777]
		#print reads_posEnd_list    #[2860, 2886, 2786, 2856, 2813, 2912, 2866]  ##2800/2801/.../2859 sixty bases were seqed!2860 not!
		
		posAndNum_counter_list = [] #1D  //not 2D #[[2800,2801,...,2859],[2834,2835,...,2885],...]
		for i in range(len(reads_posStart_list)):
			pos_start = reads_posStart_list[i]
			pos_end = reads_posEnd_list[i]
			each_startAndENd_pair_list = range(pos_start,pos_end)  #[2800,2801,...,2859]
			#posAndNum_counter_list.append(each_startAndENd_pair_list)  
			posAndNum_counter_list = posAndNum_counter_list+each_startAndENd_pair_list
		#print posAndNum_counter_list  
	
		posAndNum_counter = Counter(posAndNum_counter_list)
		posAndNum_counter_most = posAndNum_counter.most_common() 
		posAndNum_counter_most_dict = dict(posAndNum_counter_most)
		posAndNum_counter_most_dict_sorted_list = sorted(posAndNum_counter_most_dict.items(),key=lambda item:item[0],reverse = False) #shengxu pailie [('KM209255.1', 5868), ('KM209277.1', 1973), 
		#print posAndNum_counter_most_dict_sorted_list

		pos_data = []
		num_data = []
		for i in range(len(posAndNum_counter_most_dict_sorted_list)):
			pos_data.append(posAndNum_counter_most_dict_sorted_list[i][0])
			num_data.append(posAndNum_counter_most_dict_sorted_list[i][1])
		#print pos_data
		#print num_data


		fig,ax = subplots(figsize=(10,5))
		x = pos_data
		y = num_data

		plot(x,y)
		fill_between(x,0,y)

		title_str = "%s base coverage with mappting aginst %s"%(querySample,queryTheAccession)
		title(title_str,fontsize=14)
		xlabel("base position")
		ylabel("coverage")
		xlim(1,accession_length)
		ylim(ymin=0)
		grid(color='lightskyblue')

		output_path = bowtie2_path_basecoverage + querySample
		if not os.path.exists(output_path):
		    os.system('mkdir %s'%output_path)

		fig_file = os.path.join(output_path,queryTheAccession+'.png')
		fig.savefig(fig_file)
		show()
	



####CONTIGS
def pipeline_contigs_sta_to_grafic(ntVirusesFamilyName,querySample,blastnColumn,queryTheAccession):
	'''
	just show one target accession's contigs in blastn output.
	##paras:
	#blastPathSorted:blast_path_sorted=result_path+'blast_result_sorted_by_%s/'%(colNumSortedStr)
	#querySample: query_sample = 'jiayan'
	return queryTheAccession's all lines.
	'''
	##accession lenth from method nt_filter_target_accessionLength:time_clock:8.693668  time_time:8.47259616852
	##accession lenth from file nt_viruses_final_seqLengthList:time_clock:0.578304  time_time:2.52949690819

	    
	family_path = 'nt_%s'%(ntVirusesFamilyName)
	family_file_name = 'nt_%s.fasta'%(ntVirusesFamilyName)

	nt_viruses_family_path = blastDB_path+family_path
	nt_viruses_family_file = os.path.join(nt_viruses_family_path,family_file_name)

	nt_viruses_family_seqLengthList_file_name = 'nt_%s_seqLengthList'%(ntVirusesFamilyName)
	nt_viruses_family_seqLengthList_file = os.path.join(nt_viruses_family_path,nt_viruses_family_seqLengthList_file_name)

	accession_length = 0
	with open(nt_viruses_family_seqLengthList_file,'r') as nt_viruses_family_seqLengthList_file_obj:
		linesList = nt_viruses_family_seqLengthList_file_obj.readlines()
		for i in range(len(linesList)):
			eachLine = linesList[i]
			if eachLine.strip('\n').split('\t')[0] == queryTheAccession:
				accession_length = int(eachLine.strip('\n').split('\t')[1])
	print accession_length

	#accession_length = nt_filter_target_accessionLength(nt_viruses_family_file,queryTheAccession)
	#print accession_length

	query_blastn_file = target_blastn_query_file(querySample,blastnColumn)
	#print query_blastn_file
	with open(query_blastn_file,'r') as query_blastn_file_obj:
		queryLinesList = query_blastn_file_obj.readlines()
		#print queryLinesList  #['NODE_1_length_33877_cov_26.581575\tKJ883521.1\t99.917\t33907\t23\t5\t1\t33907\t108\t34009\t0.0\t62489\n',...
		
		target_identity_list = []
		target_length_list = []
		target_refStartPos_list = []
		target_refEndPos_list = []
		for i in range(len(queryLinesList)):
			accession = queryLinesList[i].split('\t')[1]
			if queryTheAccession == accession:
				identity = queryLinesList[i].split('\t')[2]
				length = queryLinesList[i].split('\t')[3]
				refStart = queryLinesList[i].split('\t')[8]
				refEnd = queryLinesList[i].split('\t')[9]
				refStartPos = min([refStart,refEnd])
				refEndPos = max([refStart,refEnd])
				target_identity_list.append(identity)
				target_length_list.append(length)
				target_refStartPos_list.append(refStartPos)
				target_refEndPos_list.append(refEndPos)
		'''
		print target_identity_list
		print target_length_list
		print target_refStartPos_list
		print target_refEndPos_list
		'''

		#fig = figure(figsize=(10,4),dpi=80)
		fig,ax = subplots(figsize=(10,2))
		hlines(0,0,accession_length,color='skyblue',linewidth=10.0)
		for i in range(len(target_length_list)):
			hlines(0.2,int(target_refStartPos_list[i]),int(target_refEndPos_list[i]),color='r',linewidth=4.0)
		title_str = "%s contigs with mapping aginst %s"%(querySample,queryTheAccession)
		title(title_str,fontsize=14)
		subplots_adjust(bottom=0.2,top=0.8)
		ax.spines['top'].set_visible(False)
		ax.spines['left'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.set_yticks([])  #no visible the y ticks

		output_path = blast_contigs_againt_ref_path + querySample		
		if not os.path.exists(output_path):
			os.system('mkdir %s'%output_path)
		fig_file = os.path.join(output_path,queryTheAccession+'_contigs.png')
		#print fig_file
		fig.savefig(fig_file)
		show()

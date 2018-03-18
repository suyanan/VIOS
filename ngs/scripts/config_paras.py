#-*- coding:utf-8 -*-
import os
import time
import sys
import re
import string
import django
import xlwt
import xlrd
import numpy
import math


from VIOS import settings

"""0.before : some paras need to change while another environment settings"""
#db_path : set the real directory of database
#db_path = "/media/analysis/2f81e2c0-4bf1-4dca-b4a6-a946b51a5145/ngs_website/database"

db_path = os.environ.get("DATABASE_HOME")  #get the environment variable:DATABASE_HOME
#print '============%s'%(db_path)

#num_threads : based on the number(CPU*cores) of user's platform , used for database/blastn/bowtie2
num_threads = 8 #24CPU*2cores
#num_threads = 500 #max_num_threads_batch1=24*2*12=576
#num_threads = 180 #max_num_threads_bigmem1=8*12*2=192
#num_threads = 8 #4CPU*2cores



"""1.basic directory setting"""
BASE_DIR = settings.BASE_DIR  #/root/ngs_website/VIOS
MEDIA_ROOT = settings.MEDIA_ROOT  #/root/ngs_website/VIOS/media


genomesDB_path = db_path+"/genomes/"

blastDB_path = db_path+'/blastDB/'
nt_path=blastDB_path+'nt/'
nt_file=os.path.join(nt_path,'nt')  #raw nt file: nt.gz --> nt  (fastq_file)
nt_index_header_file=os.path.join(nt_path,'index_nt_title')


#root_path = MEDIA_ROOT+"/"+time.strftime("%Y")+"/"+time.strftime("%m")+"/"+time.strftime("%d")+"/"  ##/root/ngs_website/VIOS/media/2017/05/10/
root_path = MEDIA_ROOT+"/"+time.strftime("%Y")+"/"+time.strftime("%m")+"/"  ##/root/ngs_website/VIOS/media/2017/05/

#raw data directory
rawData_path = root_path+"raw_datas/"  #/root/ngs_website/VIOS/media/2017/05/raw_datas/


##function:detemine blastn's parameter :  -max_target_seqs
##os.system('grep "^>" -c nt')
def paraMaxTargetSeqs(refFormatFileAP):
    refFormatFile = open(refFormatFileAP,'r')
    refFormatLinesList = refFormatFile.readlines()
    counter = 0
    for i in range(len(refFormatLinesList)):
        if refFormatLinesList[i][0] == '>':
            counter += 1
    return counter
#max_target_seqs = paraMaxTargetSeqs(nt_path)  ##big memory : server
#max_target_seqs = 40426590  ##Also, use can in the index_nt_title


"""2.database handler"""
#nt_viruses_family_name_list = ['viruses_final',]  #0 means viruses full database
#nt_viruses_family_name_readable_list = ['viruses_final',]  #0 means viruses full database

#nt_viruses_family_name_list = ['viruses_final','viruses_full']  #0 means viruses full database
#nt_viruses_family_name_readable_list = ['viruses_final','viruses_full']  #0 means viruses full database

nt_viruses_family_name_list = ['viruses_final','viruses_beforeRemovalLongSeqs','viruses_beforeRemovalOfHg38','viruses_full']  #0 means viruses full database
nt_viruses_family_name_readable_list = ['smart_virus_database','viruses_beforeRemovalLongSeqs','viruses_beforeRemovalOfHg38','viruses_full']  #0 means viruses full database

#nt_viruses_family_name_list = ['viruses_final','viruses_full','bacterias_full','viruses_Adenoviridae','viruses_Coronaviridae','viruses_Filoviridae','viruses_Flaviviridae','viruses_Orthomyxoviridae']  #0 means viruses full database
#nt_viruses_family_name_readable_list = ['viruses','viruses_full','bacterias_full','Adenovirus','SARS','Ebola_virus','Zika_virus','Influenza']  #0 means viruses full database

def get_nt_viruses_family_results(ntVirusesFamilyNameList):
    nt_viruses_family_file_list = []
    nt_viruses_family_path_list = []
    bt2_index_base_name_list = []
    nt_viruses_family_index_header_file_list = []
    nt_viruses_family_seqLengthList_file_list = []
    for i in range(len(ntVirusesFamilyNameList)):
        family_path = 'nt_%s'%(ntVirusesFamilyNameList[i])
        family_file_name = 'nt_%s.fasta'%(ntVirusesFamilyNameList[i])

        nt_viruses_family_path = blastDB_path+family_path
        nt_viruses_family_file = os.path.join(nt_viruses_family_path,family_file_name)
        nt_viruses_family_path_list.append(nt_viruses_family_path)
        nt_viruses_family_file_list.append(nt_viruses_family_file)
            
        bt2_index_base_name = 'bt2_index_%s'%(ntVirusesFamilyNameList[i])
        bt2_index_base_name_list.append(bt2_index_base_name)

        nt_viruses_family_index_header_file_name = 'nt_%s_index_title'%(ntVirusesFamilyNameList[i])  #0_index_title / Adenoviridae_index_title
        nt_viruses_family_index_header_file = os.path.join(nt_viruses_family_path,nt_viruses_family_index_header_file_name)
        nt_viruses_family_index_header_file_list.append(nt_viruses_family_index_header_file)

        nt_viruses_family_seqLengthList_file_name = 'nt_%s_seqLengthList'%(ntVirusesFamilyNameList[i])
        nt_viruses_family_seqLengthList_file = os.path.join(nt_viruses_family_path,nt_viruses_family_seqLengthList_file_name)
        nt_viruses_family_seqLengthList_file_list.append(nt_viruses_family_seqLengthList_file)
    return nt_viruses_family_file_list,nt_viruses_family_path_list,bt2_index_base_name_list,nt_viruses_family_index_header_file_list,nt_viruses_family_seqLengthList_file_list
    '''
    print nt_viruses_family_file_list  #['/media/analysis/../ngs_website/database/blastDB/nt_viruses_Adenoviridae/nt_viruses_Adenoviridae_len1000.fasta',
    print nt_viruses_family_path_list  #['/media/analysis/../ngs_website/database/blastDB/nt_viruses_Adenoviridae',
    print bt2_index_base_name_list  #['bt2_index_Adenoviridae', 
    print nt_viruses_family_index_header_file_list  #['/media/analysis/../ngs_website/database/blastDB/nt_viruses_Adenoviridae/Adenoviridae_index_title',
    '''


genome_hosts_name_list = []
genome_hosts_name_readable_list = []
#genome_hosts_name_list = ['hg38','mm10','rn6','galGal5','papAnu2','panTro5','gorGor5','rheMac8','felCat8','susScr3','oryCun2','oviAri3','myoLuc2']  #
#genome_hosts_name_readable_list = ['human','house_mouse','Norway_rat','chicken','baboon','chimpanzee','gorilla','rhesus_monkey','cat','pig','rabbit','sheep','microbat']
##add two to test:
#'bosTau8','equCab2',
#'cow','horse',

def get_genome_hosts_results(genomeHostsNameList):
    genome_hosts_file_list = []
    genome_hosts_path_list = []
    bt2_index_base_name_list = []
    for i in range(len(genomeHostsNameList)):
        path_name = genomeHostsNameList[i]
        file_name = '%s.2bit'%(genomeHostsNameList[i])
        
        genome_hosts_path = genomesDB_path+path_name
        genome_hosts_file = os.path.join(genome_hosts_path,file_name)
        genome_hosts_path_list.append(genome_hosts_path)
        genome_hosts_file_list.append(genome_hosts_file)
            
        bt2_index_base_name = 'bt2_index_%s'%(genomeHostsNameList[i])
        bt2_index_base_name_list.append(bt2_index_base_name)     
    return genome_hosts_file_list,genome_hosts_path_list,bt2_index_base_name_list

genome_hosts_name_adder_list = []
genome_hosts_name_readable_adder_list = []
genome_host_adder_file = os.path.join(db_path,'genome_host_adder.log')
if os.path.exists(genome_host_adder_file):
    with open(genome_host_adder_file,'r') as genome_host_adder_file_obj:
        genome_hosts_adder_list = genome_host_adder_file_obj.readlines()
        for i in range(len(genome_hosts_adder_list)):
            each_genome_host = genome_hosts_adder_list[i].strip('\n')
            each_genome_host_name = each_genome_host.split('\t')[0]
            each_genome_host_name_readable = each_genome_host.split('\t')[1]
            genome_hosts_name_adder_list.append(each_genome_host_name)
            genome_hosts_name_readable_adder_list.append(each_genome_host_name_readable)

genome_hosts_updated_list = []
genome_hosts_name_updated_list = genome_hosts_name_list+genome_hosts_name_adder_list
genome_hosts_name_readable_updated_list = genome_hosts_name_readable_list+genome_hosts_name_readable_adder_list



"""3.sample handleing"""
if os.path.exists(rawData_path):
    raw_data_files = os.listdir(rawData_path)
    raw_data_files.sort()

    ####pe sample datas
    raw_data_forward_name_list = []
    raw_data_reverse_name_list = []
    samplesRawDataForwardFiles = []
    samplesRawDataReverseFiles = []
    for i in range(len(raw_data_files)):
        raw_data_files_tuple = os.path.splitext(raw_data_files[i])
        raw_data_files_tuple_tuple = os.path.splitext(raw_data_files_tuple[0])
        if raw_data_files_tuple[1] == '.fastq' and len(raw_data_files)%2 == 0:  #.fastq file is ok , and raw_data_files_tuple_tuple[1] == '' avoid 'sample.v2.fastq', and len(raw_data_files)%2 == 0 show PE,not SE
            if i%2 == 0 :
                raw_data_forward_name = raw_data_files_tuple_tuple[0]
                raw_data_forward_name_list.append(raw_data_forward_name)
                raw_data_forward_file = os.path.join(rawData_path,raw_data_files[i])           
                samplesRawDataForwardFiles.append(raw_data_forward_file)
                
            else:
                raw_data_reverse_name = raw_data_files_tuple_tuple[0]
                raw_data_reverse_name_list.append(raw_data_reverse_name)
                raw_data_reverse_file = os.path.join(rawData_path,raw_data_files[i])
                samplesRawDataReverseFiles.append(raw_data_reverse_file)
        if raw_data_files_tuple[1] == '.gz' and raw_data_files_tuple_tuple[1] == '.fastq' and len(raw_data_files)%2 == 0:  #.fastq.gz file is ok
            if i%2 == 0 :
                raw_data_forward_name = raw_data_files_tuple_tuple[0]
                raw_data_forward_name_list.append(raw_data_forward_name)
                raw_data_forward_file = os.path.join(rawData_path,raw_data_files[i])           
                samplesRawDataForwardFiles.append(raw_data_forward_file)
                
            else:
                raw_data_reverse_name = raw_data_files_tuple_tuple[0]
                raw_data_reverse_name_list.append(raw_data_reverse_name)
                raw_data_reverse_file = os.path.join(rawData_path,raw_data_files[i])
                samplesRawDataReverseFiles.append(raw_data_reverse_file)

    sample_list = []

    for i in range(len(raw_data_forward_name_list)):
        for j in range(len(raw_data_forward_name_list[i])):
            if raw_data_forward_name_list[i][j] != raw_data_reverse_name_list[i][j]:
                if raw_data_forward_name_list[i][j-1] == '_' or raw_data_forward_name_list[i][j-1] == '.':  #sample_1.fastq
                    sample_name = raw_data_forward_name_list[i][0:j-1]
                else:  #sample1.fastq
                    sample_name = raw_data_forward_name_list[i][0:j]
                sample_list.append(sample_name)
    print "pair end sequencing"
    print sample_list

    ####se sample datas
    sample_se_list = []
    samplesRawDataSingleFiles = []
    for i in range(len(raw_data_files)):
        raw_data_files_tuple = os.path.splitext(raw_data_files[i])
        raw_data_files_tuple_tuple = os.path.splitext(raw_data_files_tuple[0])
        if raw_data_files_tuple[1] == '.fastq':  #.fastq file is ok , and raw_data_files_tuple_tuple[1] == '' avoid 'sample.v2.fastq'
            sample_se_name = raw_data_files_tuple_tuple[0]
            sample_se_list.append(sample_se_name)
            raw_data_single_file = os.path.join(rawData_path,raw_data_files[i])
            samplesRawDataSingleFiles.append(raw_data_single_file)
        if raw_data_files_tuple[1] == '.gz' and raw_data_files_tuple_tuple[1] == '.fastq':  #.fastq.gz file is ok
            sample_se_name = raw_data_files_tuple_tuple[0]
            sample_se_list.append(sample_se_name)
            raw_data_single_file = os.path.join(rawData_path,raw_data_files[i])
            samplesRawDataSingleFiles.append(raw_data_single_file)

    print "single end sequencing"
    print sample_se_list
    

"""4.result directory AND pipelines'paras"""
#result directory
result_path = root_path+'results/'
result_total_sta_path = result_path+'0.TOTAL_STA'

##trimmomatic directory and paras
#trimmomatic_path=result_path+'trimmomatic_result/'


##bowtie2 directory and paras
bowtie2_path=result_path+'1.bowtie2_result/'
bowtie2IndexHG = 'hg38bt2Index'
bowtie2IndexVirus = 'virusbt2Index'

#some common files ,need to handle next,like reads stat.
map_gb_hg_singleAndunpair_fastq_file_suffix = '_unpair_host.fastq'
map_nt_viruses_singleAndunpair_fastq_file_suffix = '_unpair_VIRUS.fastq'  #for SE PE , bt2 results
map_nt_viruses_singleAndunpair_sam_file_suffix = '_VIRUS.sam'

##velvet directory and paras
velvet_path=result_path+'2.velvet_result/'

kmerMin = 31
#kmerMax = 35  #small value for test!!!
kmerMax = 127
kmerStep = 2


##blastn directory and paras
blast_path=result_path+'3.blast_result/'



"""4.show results and user query settings on website"""
##4.1.show blastn result

#4.1.0.set the user's query settings of blastn for test!!!
#query_sample = 'jiayan'
#query_sample_kmer = '31'

#4.1.1.set the sorting several methods:score/identisties/evalue, default score
#change to default top 10 lines!!!???

#for method of sort_blast_result_by_blastnColumn(blastnColumn)
colNumSorted_list = [2,3,11]  #identities = 2,alingnmentLength = 3,score = 11
blastnLengthThreshold = 300

#for name of directory of blast_path_sorted
colNumSortedStr_list = ['identity','length','score']

#4.1.2.blastn_sorted directory
blast_path_sorted_list = []
for i in range(len(colNumSortedStr_list)):
    blast_path_sorted=result_path+'3.blast_result_sorted_by_%s/'%(colNumSortedStr_list[i])
    blast_path_sorted_list.append(blast_path_sorted)

blast_contigs_againt_ref_path = result_path+'3.blast_result_contigs/'

#4.1.3.top lines
#setShowTopLines = 2  #small value for test!!!
setShowTopLines = 10


##4.2.show bowtie2 result
nt_viruses_the_accession_path = db_path +'/nt_viruses_the_accession/'
bowtie2_path_basecoverage = result_path+'1.bowtie2_result_coverage/'

##4.3.other paras, set by admin
sizehint = 209715200  #200M = 200*1024*1024, the B as the danwei



#-*- coding:utf-8 -*-

from .config_paras import *
from .db_nt_process import *


if __name__ == '__main__':
    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    os.environ['DJANGO_SETTINGS_MODULE'] = 'VIOS.settings'
    django.setup()

def search_nt_title(ntVirusesFamilyName,queryRefName):
    ##paras:
    #ntVirusesFamilyName to target the IndexPathFile:nt_index_path=blastDB_path+'nt/'+'index_nt_title'
    #queryRefName: queryAssessionStr = ref (show_top_lines_blastn method)

    ####2D list : blockLinesList2D , element is EACH BLOCK from readlines(sizehint)
    family_path = 'nt_%s'%(ntVirusesFamilyName)
    family_file_name = 'nt_%s.fasta'%(ntVirusesFamilyName)
   
    nt_viruses_family_path = blastDB_path+family_path
    nt_viruses_family_file = os.path.join(nt_viruses_family_path,family_file_name)
        
    bt2_index_base_name = 'bt2_index_%s'%(ntVirusesFamilyName)

    nt_viruses_family_index_header_file_name = 'nt_%s_index_title'%(ntVirusesFamilyName)  #0_index_title / Adenoviridae_index_title
    nt_viruses_family_index_header_file = os.path.join(nt_viruses_family_path,nt_viruses_family_index_header_file_name)

    
    with open(nt_viruses_family_index_header_file,'r') as indexNTFile:
        linesList = indexNTFile.readlines()
        for i in range(len(linesList)):
            eachline = linesList[i].strip('\n')
            accession = eachline.split('\t')[0]
            annotation = eachline.split('\t')[1]
            if accession == queryRefName:
                return annotation
            if accession != queryRefName and i == len(linesList):##no found and read the file final line
                return 'NO FOUND in assigned reference'


def target_blastn_query_file(querySample,blastnColumn):
    ##paras:
    #blastnColumn, to get the blast_path_sorted
    #querySample: query_sample = 'jiayan'
    ##return the target blast_path_sorted/sample.out file        
    
    if blastnColumn == 2:
        blast_path_sorted=blast_path_sorted_list[0]  ##judge the directory name by blastnColumn
    if blastnColumn == 3:
        blast_path_sorted=blast_path_sorted_list[1]
    if blastnColumn == 11:
        blast_path_sorted=blast_path_sorted_list[2]
    #print querySample
    #if os.path.exists(blast_path_sorted):
    blastFiles = os.listdir(blast_path_sorted)
    for i in range(len(blastFiles)):
        blastnSample = blastFiles[i].split('.')[0]
        if querySample == blastnSample:
            blastFileAP = os.path.join(blast_path_sorted,blastFiles[i])
            return blastFileAP

def show_blastn_top_lines(ntVirusesFamilyName,querySample,blastnColumn):
    '''
    method to get top_lines of blastn resutls.
    Next,first line's ref,targetAnnotation can used for bowtie2 first.
    And,from blastn targetAnnotation URL can used for contigs sta.
    paras:
    blastnColumn to target blastPathSorted:blast_path_sorted=result_path+'blast_result_sorted_by_%s/'%(colNumSortedStr)
    querySample: query_sample = 'jiayan'
    ntVirusesFamilyName: transfer to method search_nt_title
    return all setShowTopLines lines list.
    '''
    
    query_blastn_file = target_blastn_query_file(querySample,blastnColumn)
    #print query_blastn_file
    
    with open(query_blastn_file,'r') as query_blastn_file_obj:       
        queryLinesList = query_blastn_file_obj.readlines()
        queryEachLineFilterList2D = []
        ##check len(queryLinesList) < setShowTopLines
        #setShowTopLines_min = min(len(queryLinesList),setShowTopLines)
        setShowTopLines_min = min(len(queryLinesList),100)
        #setShowTopLines_min = len(queryLinesList)
        for k in range(setShowTopLines_min):##len(queryLinesList) refer to all lines
            queryEachLineFilterList = []

            ref = queryLinesList[k].split('\t')[1]
            identity = float(queryLinesList[k].split('\t')[2])
            alignment_length = int(queryLinesList[k].split('\t')[3])
            pos_query_start = int(queryLinesList[k].split('\t')[6])
            pos_query_end = int(queryLinesList[k].split('\t')[7])
            pos_ref_start = int(queryLinesList[k].split('\t')[8])
            pos_ref_end = int(queryLinesList[k].split('\t')[9])
            e_value = float(queryLinesList[k].split('\t')[10])
            bit_score = float(queryLinesList[k].split('\t')[11][:-2])  #\n

            targetAnnotation = search_nt_title(ntVirusesFamilyName,ref)
            
            queryEachLineFilterList = [ref,targetAnnotation,identity,alignment_length,pos_query_start,pos_query_end,pos_ref_start,pos_ref_end,e_value,bit_score]
            queryEachLineFilterList2D.append(queryEachLineFilterList)
        #print queryEachLineFilterList2D  ##all lines, to show up in blastn_resutls.
        return queryEachLineFilterList2D

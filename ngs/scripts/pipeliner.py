#-*- coding:utf-8 -*-

from .config_paras import *
from .db_nt_process import *



if __name__ == '__main__':
    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    os.environ['DJANGO_SETTINGS_MODULE'] = 'VIOS.settings'
    django.setup()


##method : (use in pipeliner.py)sort the blastn result accoding to the para 'colNumSorted'
def sort_blast_result_by_blastnColumn(blastnColumn):
    ##filter by blastnColumn(int): colomn 2.identity 3.length(filter below 300,before all filter this) 11.score
    ##sorted the all blast result, use the different colomn    

    #if os.path.exists(blast_path):
    blastOutFiles = os.listdir(blast_path)
    #print blastOutFiles
    blastOutFilesAP = []
    for i in range(len(blastOutFiles)):
        blastOutFileAP = '%s%s'%(blast_path,blastOutFiles[i])
        if os.path.getsize(blastOutFileAP) != 0:
            blastOutFilesAP.append(blastOutFileAP)


    for i in range(len(blastOutFilesAP)):
        eachBlastOutFile = open(blastOutFilesAP[i],'r')
        try:
            fileLinesList = eachBlastOutFile.read().strip('\n').split('\n')
            fileLinesListFormat = []  #big  list including all lines of one file
            for j in range(len(fileLinesList)):
                fileEachLine = fileLinesList[j].split('\t')
                fileEachLine[2] = float(fileEachLine[2])  #third column : percent identity
                fileEachLine[3] = int(fileEachLine[3])  ##4th column : alignment length
                fileEachLine[4] = int(fileEachLine[4])
                fileEachLine[5] = int(fileEachLine[5])
                fileEachLine[6] = int(fileEachLine[6])
                fileEachLine[7] = int(fileEachLine[7])
                fileEachLine[8] = int(fileEachLine[8])
                fileEachLine[9] = int(fileEachLine[9])
                fileEachLine[10] = float(fileEachLine[10])  
                fileEachLine[11] = float(fileEachLine[11])  #12nd clolun : bit socre
                if fileEachLine[3] >= blastnLengthThreshold:  ##first get above alignment length value,then sort by three column
                    fileLinesListFormat.append(fileEachLine)
            #-------sorted to the formatted list---------
            #print sorted(fileLinesListFormat,cmp = lambda x,y:cmp(x[colNumSorted],y[colNumSorted],reverse = True))
            fileLinesList_New = sorted(fileLinesListFormat,key = lambda fileLinesListFormat:fileLinesListFormat[blastnColumn],reverse = True) #reverse = True bi
            
            #eachBlastOutFile_s = '%ssorted_%s'%(blast_path_sorted,blastOutFiles[i])
            if blastnColumn == 2:
                blast_path_sorted=blast_path_sorted_list[0]  ##judge the directory name by blastnColumn
            if blastnColumn == 3:
                blast_path_sorted=blast_path_sorted_list[1]
            if blastnColumn == 11:
                blast_path_sorted=blast_path_sorted_list[2]
            eachBlastOutFile_s = '%s%s'%(blast_path_sorted,blastOutFiles[i])
            #if os.path.exists(eachBlastOutFile_s):

            eachBlastOutFile_New = open(eachBlastOutFile_s,'w+')
            #-------write the formatted_list to the new file----------
            try:
                for ii in range(len(fileLinesList_New)):
                #for ii in range(len(fileLinesList_New)):
                    for jj in range(len(fileLinesList_New[ii])):
                        if jj != len(fileLinesList_New[ii])-1:
                            fileLinesList_New[ii][jj] = str(fileLinesList_New[ii][jj])+'\t'
                        else:
                            fileLinesList_New[ii][jj] = str(fileLinesList_New[ii][jj])+'\n'
                    eachBlastOutFile_New.writelines(fileLinesList_New[ii])  #zhu hang xie ru
            finally:
                eachBlastOutFile_New.close()
        finally:
            eachBlastOutFile.close()

##method : justMapAginstVirusWithoutQC , NO HOST GENOME
def process_pipeliner_pe_NoHostGenome(ntVirusesFamilyName):
    ##rebuild : use all genome host
    
    ##rebuild : use ntVirusesFamilyName to replace nt_viruses_family_name_list[i], used for mapping with bowtie2 on nt_viruses_ and blastn
    family_path = 'nt_%s'%(ntVirusesFamilyName)
    family_file_name = 'nt_%s.fasta'%(ntVirusesFamilyName)
    
    nt_viruses_family_path = blastDB_path+family_path
    nt_viruses_family_file = os.path.join(nt_viruses_family_path,family_file_name)
        
    nt_virus_bt2_index_base_name = 'bt2_index_%s'%(ntVirusesFamilyName)

    nt_viruses_family_index_header_file_name = 'nt_%s_index_title'%(ntVirusesFamilyName)  #0_index_title / Adenoviridae_index_title
    nt_viruses_family_index_header_file = os.path.join(nt_viruses_family_path,nt_viruses_family_index_header_file_name)

    #print sample_list

    time_total_start = time.clock()
    time_total_start1 = time.time()

    pipelinerSta = xlwt.Workbook()
    sheetPipelinerStaClock = pipelinerSta.add_sheet('pipeline_clock_statistics',cell_overwrite_ok=True)
    sheetPipelinerStaTime = pipelinerSta.add_sheet('pipeline_time_statistics',cell_overwrite_ok=True)
    #rowTitle = ['qc_trimmomatic','map_hp38','map_genebank_virus','map_nt_virus','assembly_velvet','map_blastn','total_time']
    #rowTitle = ['qc_trimmomatic','map_hp38','map_nt_virus','assembly_velvet','map_blastn','total_time']
    rowTitle = ['map_%s'%ntVirusesFamilyName,'assembly_velvet','map_blastn','total_time']
    for i in range(0,len(rowTitle)):
        sheetPipelinerStaClock.write(0,i,rowTitle[i]) #first row
        sheetPipelinerStaTime.write(0,i,rowTitle[i])



    print '==========map and align : AGINST NT_viruses reference,keep mapped reads'

    #bowtie2IndexNTVirusFamily = nt_viruses_family_name_list[0]
    #bowtie2IndexNTVirusFamily = bt2_index_base_name_list[0]  #'bt2_index_Adenoviridae'

    timer_bowtie2_nt_virus_start = time.clock()
    timer_bowtie2_nt_virus_start1 = time.time()

    outputFastqPair1FileNTVirus_list = []
    outputFastqPair2FileNTVirus_list = []
    outputFastqUnpairFileNTVirus_list = []
    for i in range(len(sample_list)):
        outputSAMFileNTVirus = bowtie2_path+sample_list[i]+map_nt_viruses_singleAndunpair_sam_file_suffix
        outputFastqPairFileNTVirus = bowtie2_path+sample_list[i]+'_pair_VIRUS.fastq'
        outputFastqPair1FileNTVirus = bowtie2_path+sample_list[i]+'_pair_VIRUS.1.fastq'
        outputFastqPair2FileNTVirus = bowtie2_path+sample_list[i]+'_pair_VIRUS.2.fastq'
        
        outputFastqUnpairFileNTVirus = bowtie2_path+sample_list[i]+map_nt_viruses_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -1 %s -2 %s --no-head --no-unal -S %s --al-conc %s --al %s -p %d'
            %(nt_viruses_family_path,nt_virus_bt2_index_base_name,samplesRawDataForwardFiles[i],samplesRawDataReverseFiles[i],outputSAMFileNTVirus,outputFastqPairFileNTVirus,outputFastqUnpairFileNTVirus,num_threads))

        outputFastqPair1FileNTVirus_list.append(outputFastqPair1FileNTVirus)
        outputFastqPair2FileNTVirus_list.append(outputFastqPair2FileNTVirus)
        outputFastqUnpairFileNTVirus_list.append(outputFastqUnpairFileNTVirus)
        
    timer_bowtie2_nt_virus_end = time.clock()
    timer_bowtie2_nt_virus_end1 = time.time()
    bowtie2_nt_virus_time = timer_bowtie2_nt_virus_end - timer_bowtie2_nt_virus_start
    bowtie2_nt_virus_time1 = timer_bowtie2_nt_virus_end1 - timer_bowtie2_nt_virus_start1
    '''
    print '@@@@@----------timer:bowtie2_nt_virus_time_clock:'
    print bowtie2_nt_virus_time    #18.209088
    print '@@@@@----------timer:bowtie2_nt_virus_time_time:'
    print bowtie2_nt_virus_time1    #1678.66673303
    '''
    sheetPipelinerStaClock.write(1,0,bowtie2_nt_virus_time)
    sheetPipelinerStaTime.write(1,0,bowtie2_nt_virus_time1)

    print '==========assembly : velvet with alignment data above'

    timer_velvet_start = time.clock()
    timer_velvet_start1 = time.time()

    for i in range(len(sample_list)):
        #inputPE = velvet_result+sample_list[i]+'_pe.fastq'
        #os.system('%s %s %s %s' %(shuffleSequences_fastq_path,outputFastq1AfterHGandVirus[i],outputFastq2AfterHGandVirus[i],inputPE))
        for kmer in range(kmerMin,kmerMax,kmerStep):
            #print sample_list[i]+'************************************'
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_list[i],kmer)
            os.system('velveth %s %d -shortPaired -fastq -separate %s %s -short -fastq %s'%(eachKmerDir,kmer,outputFastqPair1FileNTVirus_list[i],outputFastqPair2FileNTVirus_list[i],outputFastqUnpairFileNTVirus_list[i]))
            os.system('velvetg %s -exp_cov auto -cov_cutoff auto ' %(eachKmerDir))


            eachKmerFiles = os.listdir(eachKmerDir)
            eachKmerFilesAP = []
            for j in range(len(eachKmerFiles)):
                eachKmerFilesAP.append(eachKmerDir+'/'+eachKmerFiles[j])
            for j in range(len(eachKmerFilesAP)):
                #pass
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Graph2' : os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/LastGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/PreGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Roadmaps': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Sequences': os.remove(eachKmerFilesAP[j])

    timer_velvet_end = time.clock()
    timer_velvet_end1 = time.time()
    velvet_time = timer_velvet_end - timer_velvet_start
    velvet_time1 = timer_velvet_end1 - timer_velvet_start1
    '''
    print '@@@@@----------timer:velvet_time_clock:'
    print velvet_time    #0.034675
    print '@@@@@----------timer:velvet_time_time:'
    print velvet_time1    #1.79084706306
    '''
    sheetPipelinerStaClock.write(1,1,velvet_time)
    sheetPipelinerStaTime.write(1,1,velvet_time1)

    ##================statistics velvet : pick up the proper kmer , accoding to N50 max========
    velvetSta = xlwt.Workbook()
    sheetVelvetSta = velvetSta.add_sheet('velvet_statistics', cell_overwrite_ok=True)
    kmerNum = (kmerMax-kmerMin)/2  #(35-31)/2=2
    row0 = ['sampleNums','kmer','medianCovDepth','nodesNum','contigN50','contigMax','genomeLen','usedReads','totalReads']
    for i in range(0,len(row0)):
        sheetVelvetSta.write(0,i,row0[i])  #first row

    #first cow is sample_list
    #second cow is kmer (kmerMin to kmerMax , kmerStep)

    for i in range(0,len(sample_list)):
        sheetVelvetSta.write(kmerNum*i+1,0,sample_list[i])   #first column is the sampleNums

    colKmer = range(kmerMin,kmerMax,kmerStep)
    for i in range(0,len(sample_list)*kmerNum,kmerNum): #i : every column
        for k in range(0,kmerNum):
            sheetVelvetSta.write(i+k+1,1,colKmer[k])   #second column is the k-mer
            k += 1

    logLast2Str = []
    MedianCovDepth = []

    logLast1Str = []
    nodesNum = []
    contigN50 = []
    contigMax = []
    genomeLen = []
    usedReads = []
    totalReads = []
    for i in range(len(sample_list)):
        for kmer in range(kmerMin,kmerMax,kmerStep):
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_list[i],kmer)
            Log_file = os.path.join(eachKmerDir,'Log')
            logs = open(Log_file,'r').readlines()
            #logLast2line = logs[len(logs)-2]
            #if logLast2line != '':
            #Median coverage depth = 96.434783  #not the logs[26] , to avoid the task again with adding the log_information to Log File
            #logLast2Str = logLast2line.split(' ')
            #print logLast2Str
            #print logLast2Str[4][:-2]
            #MedianCovDepth.append(string.atof(logLast2Str[4][:-2]))
            logLast1line = logs[len(logs)-1]  #Final graph has 199 nodes and n50 of 47325, max 166302, total 1570595, using 1470371/1487666 reads
            logLast1Str = logLast1line.split(' ')
            #['Final', 'graph', 'has', '172', 'nodes', 'and', 'n50', 'of', '46917,', 'max', '237233,', 'total', '1570012,', 'using', '1478778/1487666', 'reads\n']
            #print logLast1Str
            nodesNum.append(string.atoi(logLast1Str[3], 10))
            contigN50.append(string.atoi(logLast1Str[8][:-1],10))
            contigMax.append(string.atoi(logLast1Str[10][:-1],10))
            genomeLen.append(string.atoi(logLast1Str[12][:-1],10))
            usedReads.append(string.atoi(logLast1Str[14].split('/')[0], 10))
            totalReads.append(string.atoi(logLast1Str[14].split('/')[1], 10))
    #print contigN50
    for i in range(0,len(sample_list)*kmerNum):
        #sheetVelvetSta.write(i+1,2,MedianCovDepth[i])
        sheetVelvetSta.write(i+1,3,nodesNum[i])
        sheetVelvetSta.write(i+1,4,contigN50[i])
        sheetVelvetSta.write(i+1,5,contigMax[i])
        sheetVelvetSta.write(i+1,6,genomeLen[i])
        sheetVelvetSta.write(i+1,7,usedReads[i])
        sheetVelvetSta.write(i+1,8,totalReads[i])

    velvetSta.save(os.path.join(result_total_sta_path,'velvet_log_sta.xls'))

    contigN50_max_value_list = []  #all_samples
    contigN50_max_index_list = []
    for i in range(len(sample_list)):
        each_sample_contigN50_list = contigN50[kmerNum*i:kmerNum*(i+1)]
        #print each_sample_contigN50_list
        contigN50_max_value_list.append(max(each_sample_contigN50_list))
        contigN50_max_index_list.append(kmerMin+(each_sample_contigN50_list.index(max(each_sample_contigN50_list)))*kmerStep)

    print contigN50_max_value_list
    print contigN50_max_index_list

    print '=============blast:use the contigs.fa to map AGINST NT_viruses'
    #1.first:format the local nt fasta_database file
    #os.system('makeblastdb -dbtype nucl -input_type fasta -in %s' %(nt_path))  ###big memory : server

    #after makeblastdb,generate some other indexed files besides nt file,
    #use the indexed nt file(with some other indexed files) as the following input ref db '-db'.

    #2.second:blastn

    timer_blastn_start = time.clock()
    timer_blastn_start1 = time.time()

    #-outfmt 6
    for i in range(len(sample_list)):
        kmer = contigN50_max_index_list[i]
        print kmer
        eachKmerDir = '%s%s_%d'%(velvet_path,sample_list[i],kmer) #/public1/home/yefq/project_adenovirus/results/velvet_result/jiayan_31
        #print eachKmerDir
        eachContigFile = os.path.join(eachKmerDir,'contigs.fa')
        #print eachContigFile
        #estimate the contigs.fa is empty or not.
        if os.path.getsize(eachContigFile) != 0:
            eachBlastResult = os.path.join(blast_path,sample_list[i]+'.out')
            os.system('blastn -db %s -query %s -outfmt 6 -dust no -num_threads %d -perc_identity 80 -word_size 20 -max_target_seqs %d -evalue 0.0000001 -out %s'
                %(nt_viruses_family_file,eachContigFile,num_threads,max_target_seqs_nt_viruses,eachBlastResult))

    timer_blastn_end = time.clock()
    timer_blastn_end1 = time.time()
    blastn_time = timer_blastn_end - timer_blastn_start
    blastn_time1 = timer_blastn_end1 - timer_blastn_start1
    '''
    print '@@@@@----------timer:blastn_time_clock:'
    print blastn_time  #0.287237
    print '@@@@@----------timer:blastn_time_time:'
    print blastn_time1  #30.0745470524
    '''
    sheetPipelinerStaClock.write(1,2,blastn_time)
    sheetPipelinerStaTime.write(1,2,blastn_time1)

    ##3.third:sorted the blast result, use the different colomn
    for i in range(len(colNumSorted_list)):
        sort_blast_result_by_blastnColumn(colNumSorted_list[i])


    time_total_end = time.clock()
    time_total_end1 = time.time()
    

    total_time = time_total_end - time_total_start
    total_time1 = time_total_end1 - time_total_start1
    '''
    print '@@@@@-----------timer:total_time_clock:'
    print total_time  #23.519992
    print '@@@@@-----------timer:total_time_time:'
    print total_time1  #2177.6757071
    '''
    sheetPipelinerStaClock.write(1,3,total_time)
    sheetPipelinerStaTime.write(1,3,total_time1)

    pipelinerSta.save(os.path.join(result_total_sta_path,'pipeline_time_sta.xls'))

##method : justMapAginstVirusWithoutQC , NO HOST GENOME
def process_pipeliner_se(ntVirusesFamilyName):
    ##rebuild : use all genome host
    
    ##rebuild : use ntVirusesFamilyName to replace nt_viruses_family_name_list[i], used for mapping with bowtie2 on nt_viruses_ and blastn
    family_path = 'nt_%s'%(ntVirusesFamilyName)
    family_file_name = 'nt_%s.fasta'%(ntVirusesFamilyName)
    
    nt_viruses_family_path = blastDB_path+family_path
    nt_viruses_family_file = os.path.join(nt_viruses_family_path,family_file_name)
        
    nt_virus_bt2_index_base_name = 'bt2_index_%s'%(ntVirusesFamilyName)

    nt_viruses_family_index_header_file_name = 'nt_%s_index_title'%(ntVirusesFamilyName)  #0_index_title / Adenoviridae_index_title
    nt_viruses_family_index_header_file = os.path.join(nt_viruses_family_path,nt_viruses_family_index_header_file_name)

    
    #print sample_se_list

    time_total_start = time.clock()
    time_total_start1 = time.time()

    pipelinerSta = xlwt.Workbook()
    sheetPipelinerStaClock = pipelinerSta.add_sheet('pipeline_clock_statistics',cell_overwrite_ok=True)
    sheetPipelinerStaTime = pipelinerSta.add_sheet('pipeline_time_statistics',cell_overwrite_ok=True)
    
    rowTitle = ['map_%s'%ntVirusesFamilyName,'assembly_velvet','map_blastn','total_time']
    for i in range(0,len(rowTitle)):
        sheetPipelinerStaClock.write(0,i,rowTitle[i]) #first row
        sheetPipelinerStaTime.write(0,i,rowTitle[i])



    print '==========map and align : AGINST NT_viruses reference,keep mapped reads'

    #bowtie2IndexNTVirusFamily = nt_viruses_family_name_list[0]
    #bowtie2IndexNTVirusFamily = bt2_index_base_name_list[0]  #'bt2_index_Adenoviridae'

    timer_bowtie2_nt_virus_start = time.clock()
    timer_bowtie2_nt_virus_start1 = time.time()

    afterNTVIRUSfastqFiles =[]
    for i in range(len(sample_se_list)):
        outputSAMFileNTVirus = bowtie2_path+sample_se_list[i]+map_nt_viruses_singleAndunpair_sam_file_suffix
        outputFASTQFileNTVirus = bowtie2_path+sample_se_list[i]+map_nt_viruses_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -U %s --no-head --no-unal -S %s --al %s -p %d'
            %(nt_viruses_family_path,nt_virus_bt2_index_base_name,samplesRawDataSingleFiles[i],outputSAMFileNTVirus,outputFASTQFileNTVirus,num_threads))

        afterNTVIRUSfastqFiles.append(outputFASTQFileNTVirus)


    timer_bowtie2_nt_virus_end = time.clock()
    timer_bowtie2_nt_virus_end1 = time.time()
    bowtie2_nt_virus_time = timer_bowtie2_nt_virus_end - timer_bowtie2_nt_virus_start
    bowtie2_nt_virus_time1 = timer_bowtie2_nt_virus_end1 - timer_bowtie2_nt_virus_start1
    '''
    print '@@@@@----------timer:bowtie2_nt_virus_time_clock:'
    print bowtie2_nt_virus_time    #18.209088
    print '@@@@@----------timer:bowtie2_nt_virus_time_time:'
    print bowtie2_nt_virus_time1    #1678.66673303
    '''
    sheetPipelinerStaClock.write(1,0,bowtie2_nt_virus_time)
    sheetPipelinerStaTime.write(1,0,bowtie2_nt_virus_time1)


    print '==========assembly : velvet with alignment data above'

    timer_velvet_start = time.clock()
    timer_velvet_start1 = time.time()

    minContigLength = 100
    for i in range(len(sample_se_list)):
        for kmer in range(kmerMin,kmerMax,kmerStep):
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_se_list[i],kmer)
            os.system('velveth %s %d -short -fastq %s' %(eachKmerDir,kmer,afterNTVIRUSfastqFiles[i]))
            os.system('velvetg %s -exp_cov auto -cov_cutoff auto -min_contig_lgth %d' %(eachKmerDir,minContigLength))

            eachKmerFiles = os.listdir(eachKmerDir)
            eachKmerFilesAP = []
            for j in range(len(eachKmerFiles)):
                eachKmerFilesAP.append(eachKmerDir+'/'+eachKmerFiles[j])
            for j in range(len(eachKmerFilesAP)):
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Graph2' : os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/LastGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/PreGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Roadmaps': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Sequences': os.remove(eachKmerFilesAP[j])

    timer_velvet_end = time.clock()
    timer_velvet_end1 = time.time()
    velvet_time = timer_velvet_end - timer_velvet_start
    velvet_time1 = timer_velvet_end1 - timer_velvet_start1
    '''
    print '@@@@@----------timer:velvet_time_clock:'
    print velvet_time    #0.034675
    print '@@@@@----------timer:velvet_time_time:'
    print velvet_time1    #1.79084706306
    '''
    sheetPipelinerStaClock.write(1,1,velvet_time)
    sheetPipelinerStaTime.write(1,1,velvet_time1)

    ##================statistics velvet : pick up the proper kmer , accoding to N50 max========
    velvetSta = xlwt.Workbook()
    sheetVelvetSta = velvetSta.add_sheet('velvet_statistics', cell_overwrite_ok=True)
    kmerNum = (kmerMax-kmerMin)/2  #(35-31)/2=2
    row0 = ['sampleNums','kmer','medianCovDepth','nodesNum','contigN50','contigMax','genomeLen','usedReads','totalReads']
    for i in range(0,len(row0)):
        sheetVelvetSta.write(0,i,row0[i])  #first row

    #first cow is sample_se_list
    #second cow is kmer (kmerMin to kmerMax , kmerStep)

    for i in range(0,len(sample_se_list)):
        sheetVelvetSta.write(kmerNum*i+1,0,sample_se_list[i])   #first column is the sampleNums

    colKmer = range(kmerMin,kmerMax,kmerStep)
    for i in range(0,len(sample_se_list)*kmerNum,kmerNum): #i : every column
        for k in range(0,kmerNum):
            sheetVelvetSta.write(i+k+1,1,colKmer[k])   #second column is the k-mer
            k += 1

    logLast2Str = []
    MedianCovDepth = []

    logLast1Str = []
    nodesNum = []
    contigN50 = []
    contigMax = []
    genomeLen = []
    usedReads = []
    totalReads = []
    for i in range(len(sample_se_list)):
        for kmer in range(kmerMin,kmerMax,kmerStep):
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_se_list[i],kmer)
            Log_file = os.path.join(eachKmerDir,'Log')
            logs = open(Log_file,'r').readlines()
            #logLast2line = logs[len(logs)-2]
            #if logLast2line != '':
            #Median coverage depth = 96.434783  #not the logs[26] , to avoid the task again with adding the log_information to Log File
            #logLast2Str = logLast2line.split(' ')
            #print logLast2Str
            #print logLast2Str[4][:-2]
            #MedianCovDepth.append(string.atof(logLast2Str[4][:-2]))
            logLast1line = logs[len(logs)-1]  #Final graph has 199 nodes and n50 of 47325, max 166302, total 1570595, using 1470371/1487666 reads
            logLast1Str = logLast1line.split(' ')
            #['Final', 'graph', 'has', '172', 'nodes', 'and', 'n50', 'of', '46917,', 'max', '237233,', 'total', '1570012,', 'using', '1478778/1487666', 'reads\n']
            #print logLast1Str
            nodesNum.append(string.atoi(logLast1Str[3], 10))
            contigN50.append(string.atoi(logLast1Str[8][:-1],10))
            contigMax.append(string.atoi(logLast1Str[10][:-1],10))
            genomeLen.append(string.atoi(logLast1Str[12][:-1],10))
            usedReads.append(string.atoi(logLast1Str[14].split('/')[0], 10))
            totalReads.append(string.atoi(logLast1Str[14].split('/')[1], 10))
    #print contigN50
    for i in range(0,len(sample_se_list)*kmerNum):
        #sheetVelvetSta.write(i+1,2,MedianCovDepth[i])
        sheetVelvetSta.write(i+1,3,nodesNum[i])
        sheetVelvetSta.write(i+1,4,contigN50[i])
        sheetVelvetSta.write(i+1,5,contigMax[i])
        sheetVelvetSta.write(i+1,6,genomeLen[i])
        sheetVelvetSta.write(i+1,7,usedReads[i])
        sheetVelvetSta.write(i+1,8,totalReads[i])

    velvetSta.save(os.path.join(result_total_sta_path,'velvet_log_sta.xls'))

    ##check the max contigN50 to pick up ,then blastn
    contigN50_max_value_list = []  #all_samples
    contigN50_max_index_list = []
    for i in range(len(sample_se_list)):
        each_sample_contigN50_list = contigN50[kmerNum*i:kmerNum*(i+1)]
        #print each_sample_contigN50_list
        contigN50_max_value_list.append(max(each_sample_contigN50_list))
        contigN50_max_index_list.append(kmerMin+(each_sample_contigN50_list.index(max(each_sample_contigN50_list)))*kmerStep)

    #print contigN50_max_value_list
    #print contigN50_max_index_list

    print '=============blast:use the contigs.fa to map AGINST NT_viruses'
    #1.first:format the local nt fasta_database file
    #os.system('makeblastdb -dbtype nucl -input_type fasta -in %s' %(nt_path))  ###big memory : server

    #after makeblastdb,generate some other indexed files besides nt file,
    #use the indexed nt file(with some other indexed files) as the following input ref db '-db'.

    #2.second:blastn

    timer_blastn_start = time.clock()
    timer_blastn_start1 = time.time()

                 
    for i in range(len(sample_se_list)):
        kmer = contigN50_max_index_list[i]
        eachKmerDir = '%s%s_%d'%(velvet_path,sample_se_list[i],kmer) #/public1/home/yefq/project_adenovirus/results/velvet_result/jiayan_31
        #print eachKmerDir
        eachContigFile = os.path.join(eachKmerDir,'contigs.fa')
        #print eachContigFile
        if os.path.getsize(eachContigFile) != 0:
            eachBlastResult = os.path.join(blast_path,sample_se_list[i]+'.out')
            os.system('blastn -db %s -query %s -outfmt 6 -dust no -num_threads %d -perc_identity 80 -word_size 20 -max_target_seqs %d -evalue 0.0000001 -out %s'
                %(nt_viruses_family_file,eachContigFile,num_threads,max_target_seqs_nt_viruses,eachBlastResult))

    timer_blastn_end = time.clock()
    timer_blastn_end1 = time.time()
    blastn_time = timer_blastn_end - timer_blastn_start
    blastn_time1 = timer_blastn_end1 - timer_blastn_start1
    '''
    print '@@@@@----------timer:blastn_time_clock:'
    print blastn_time  #0.287237
    print '@@@@@----------timer:blastn_time_time:'
    print blastn_time1  #30.0745470524
    '''
    sheetPipelinerStaClock.write(1,2,blastn_time)
    sheetPipelinerStaTime.write(1,2,blastn_time1)

    ##3.third:sorted the blast result, use the different colomn
    for i in range(len(colNumSorted_list)):
        sort_blast_result_by_blastnColumn(colNumSorted_list[i])


    time_total_end = time.clock()
    time_total_end1 = time.time()
    

    total_time = time_total_end - time_total_start
    total_time1 = time_total_end1 - time_total_start1
    '''
    print '@@@@@-----------timer:total_time_clock:'
    print total_time  #23.519992
    print '@@@@@-----------timer:total_time_time:'
    print total_time1  #2177.6757071
    '''
    sheetPipelinerStaClock.write(1,3,total_time)
    sheetPipelinerStaTime.write(1,3,total_time1)

    pipelinerSta.save(os.path.join(result_total_sta_path,'pipeline_time_sta.xls'))



##method : withGenomeHost -----**---------
def process_pipeliner_pe(genomeHostName,ntVirusesFamilyName):
    ##rebuild : use all genome host
    genome_host_path = genomesDB_path+genomeHostName
    genome_host_file = os.path.join(genome_host_path,genomeHostName+'.fa')

    genome_host_bt2_index_base_name = 'bt2_index_%s'%(genomeHostName)
    
    ##rebuild : use ntVirusesFamilyName to replace nt_viruses_family_name_list[i], used for mapping with bowtie2 on nt_viruses_ and blastn
    family_path = 'nt_%s'%(ntVirusesFamilyName)
    family_file_name = 'nt_%s.fasta'%(ntVirusesFamilyName)
    
    nt_viruses_family_path = blastDB_path+family_path
    nt_viruses_family_file = os.path.join(nt_viruses_family_path,family_file_name)
        
    nt_virus_bt2_index_base_name = 'bt2_index_%s'%(ntVirusesFamilyName)

    nt_viruses_family_index_header_file_name = 'nt_%s_index_title'%(ntVirusesFamilyName)  #0_index_title / Adenoviridae_index_title
    nt_viruses_family_index_header_file = os.path.join(nt_viruses_family_path,nt_viruses_family_index_header_file_name)

    #print sample_list

    time_total_start = time.clock()
    time_total_start1 = time.time()

    pipelinerSta = xlwt.Workbook()
    sheetPipelinerStaClock = pipelinerSta.add_sheet('pipeline_clock_statistics',cell_overwrite_ok=True)
    sheetPipelinerStaTime = pipelinerSta.add_sheet('pipeline_time_statistics',cell_overwrite_ok=True)
    #rowTitle = ['qc_trimmomatic','map_hp38','map_genebank_virus','map_nt_virus','assembly_velvet','map_blastn','total_time']
    #rowTitle = ['qc_trimmomatic','map_hp38','map_nt_virus','assembly_velvet','map_blastn','total_time']
    rowTitle = ['map_%s'%ntVirusesFamilyName,'assembly_velvet','map_blastn','total_time']
    for i in range(0,len(rowTitle)):
        sheetPipelinerStaClock.write(0,i,rowTitle[i]) #first row
        sheetPipelinerStaTime.write(0,i,rowTitle[i])



    print '==========map and align : AGINST NT_viruses reference,keep mapped reads'

    #bowtie2IndexNTVirusFamily = nt_viruses_family_name_list[0]
    #bowtie2IndexNTVirusFamily = bt2_index_base_name_list[0]  #'bt2_index_Adenoviridae'

    timer_bowtie2_nt_virus_start = time.clock()
    timer_bowtie2_nt_virus_start1 = time.time()

    outputFastqPair1FileNTVirus_list = []
    outputFastqPair2FileNTVirus_list = []
    outputFastqUnpairFileNTVirus_list = []
    for i in range(len(sample_list)):
        outputSAMFileNTVirus = bowtie2_path+sample_list[i]+map_nt_viruses_singleAndunpair_sam_file_suffix
        outputFastqPairFileNTVirus = bowtie2_path+sample_list[i]+'_pair_VIRUS.fastq'
        outputFastqPair1FileNTVirus = bowtie2_path+sample_list[i]+'_pair_VIRUS.1.fastq'
        outputFastqPair2FileNTVirus = bowtie2_path+sample_list[i]+'_pair_VIRUS.2.fastq'
        
        outputFastqUnpairFileNTVirus = bowtie2_path+sample_list[i]+map_nt_viruses_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -1 %s -2 %s --no-head --no-unal -S %s --al-conc %s --al %s -p %d'
            %(nt_viruses_family_path,nt_virus_bt2_index_base_name,samplesRawDataForwardFiles[i],samplesRawDataReverseFiles[i],outputSAMFileNTVirus,outputFastqPairFileNTVirus,outputFastqUnpairFileNTVirus,num_threads))

        outputFastqPair1FileNTVirus_list.append(outputFastqPair1FileNTVirus)
        outputFastqPair2FileNTVirus_list.append(outputFastqPair2FileNTVirus)
        outputFastqUnpairFileNTVirus_list.append(outputFastqUnpairFileNTVirus)
        
    timer_bowtie2_nt_virus_end = time.clock()
    timer_bowtie2_nt_virus_end1 = time.time()
    bowtie2_nt_virus_time = timer_bowtie2_nt_virus_end - timer_bowtie2_nt_virus_start
    bowtie2_nt_virus_time1 = timer_bowtie2_nt_virus_end1 - timer_bowtie2_nt_virus_start1
    '''
    print '@@@@@----------timer:bowtie2_nt_virus_time_clock:'
    print bowtie2_nt_virus_time    #18.209088
    print '@@@@@----------timer:bowtie2_nt_virus_time_time:'
    print bowtie2_nt_virus_time1    #1678.66673303
    '''
    sheetPipelinerStaClock.write(1,0,bowtie2_nt_virus_time)
    sheetPipelinerStaTime.write(1,0,bowtie2_nt_virus_time1)

    print '==========assembly : velvet with alignment data above'

    timer_velvet_start = time.clock()
    timer_velvet_start1 = time.time()

    for i in range(len(sample_list)):
        #inputPE = velvet_result+sample_list[i]+'_pe.fastq'
        #os.system('%s %s %s %s' %(shuffleSequences_fastq_path,outputFastq1AfterHGandVirus[i],outputFastq2AfterHGandVirus[i],inputPE))
        for kmer in range(kmerMin,kmerMax,kmerStep):
            #print sample_list[i]+'************************************'
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_list[i],kmer)
            os.system('velveth %s %d -shortPaired -fastq -separate %s %s -short -fastq %s'%(eachKmerDir,kmer,outputFastqPair1FileNTVirus_list[i],outputFastqPair2FileNTVirus_list[i],outputFastqUnpairFileNTVirus_list[i]))
            os.system('velvetg %s -exp_cov auto -cov_cutoff auto ' %(eachKmerDir))


            eachKmerFiles = os.listdir(eachKmerDir)
            eachKmerFilesAP = []
            for j in range(len(eachKmerFiles)):
                eachKmerFilesAP.append(eachKmerDir+'/'+eachKmerFiles[j])
            for j in range(len(eachKmerFilesAP)):
                #pass
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Graph2' : os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/LastGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/PreGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Roadmaps': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Sequences': os.remove(eachKmerFilesAP[j])

    timer_velvet_end = time.clock()
    timer_velvet_end1 = time.time()
    velvet_time = timer_velvet_end - timer_velvet_start
    velvet_time1 = timer_velvet_end1 - timer_velvet_start1
    '''
    print '@@@@@----------timer:velvet_time_clock:'
    print velvet_time    #0.034675
    print '@@@@@----------timer:velvet_time_time:'
    print velvet_time1    #1.79084706306
    '''
    sheetPipelinerStaClock.write(1,1,velvet_time)
    sheetPipelinerStaTime.write(1,1,velvet_time1)

    ##================statistics velvet : pick up the proper kmer , accoding to N50 max========
    velvetSta = xlwt.Workbook()
    sheetVelvetSta = velvetSta.add_sheet('velvet_statistics', cell_overwrite_ok=True)
    kmerNum = (kmerMax-kmerMin)/2  #(35-31)/2=2
    row0 = ['sampleNums','kmer','medianCovDepth','nodesNum','contigN50','contigMax','genomeLen','usedReads','totalReads']
    for i in range(0,len(row0)):
        sheetVelvetSta.write(0,i,row0[i])  #first row

    #first cow is sample_list
    #second cow is kmer (kmerMin to kmerMax , kmerStep)

    for i in range(0,len(sample_list)):
        sheetVelvetSta.write(kmerNum*i+1,0,sample_list[i])   #first column is the sampleNums

    colKmer = range(kmerMin,kmerMax,kmerStep)
    for i in range(0,len(sample_list)*kmerNum,kmerNum): #i : every column
        for k in range(0,kmerNum):
            sheetVelvetSta.write(i+k+1,1,colKmer[k])   #second column is the k-mer
            k += 1

    logLast2Str = []
    MedianCovDepth = []

    logLast1Str = []
    nodesNum = []
    contigN50 = []
    contigMax = []
    genomeLen = []
    usedReads = []
    totalReads = []
    for i in range(len(sample_list)):
        for kmer in range(kmerMin,kmerMax,kmerStep):
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_list[i],kmer)
            Log_file = os.path.join(eachKmerDir,'Log')
            logs = open(Log_file,'r').readlines()
            #logLast2line = logs[len(logs)-2]
            #if logLast2line != '':
            #Median coverage depth = 96.434783  #not the logs[26] , to avoid the task again with adding the log_information to Log File
            #logLast2Str = logLast2line.split(' ')
            #print logLast2Str
            #print logLast2Str[4][:-2]
            #MedianCovDepth.append(string.atof(logLast2Str[4][:-2]))
            logLast1line = logs[len(logs)-1]  #Final graph has 199 nodes and n50 of 47325, max 166302, total 1570595, using 1470371/1487666 reads
            logLast1Str = logLast1line.split(' ')
            #['Final', 'graph', 'has', '172', 'nodes', 'and', 'n50', 'of', '46917,', 'max', '237233,', 'total', '1570012,', 'using', '1478778/1487666', 'reads\n']
            #print logLast1Str
            nodesNum.append(string.atoi(logLast1Str[3], 10))
            contigN50.append(string.atoi(logLast1Str[8][:-1],10))
            contigMax.append(string.atoi(logLast1Str[10][:-1],10))
            genomeLen.append(string.atoi(logLast1Str[12][:-1],10))
            usedReads.append(string.atoi(logLast1Str[14].split('/')[0], 10))
            totalReads.append(string.atoi(logLast1Str[14].split('/')[1], 10))
    #print contigN50
    for i in range(0,len(sample_list)*kmerNum):
        #sheetVelvetSta.write(i+1,2,MedianCovDepth[i])
        sheetVelvetSta.write(i+1,3,nodesNum[i])
        sheetVelvetSta.write(i+1,4,contigN50[i])
        sheetVelvetSta.write(i+1,5,contigMax[i])
        sheetVelvetSta.write(i+1,6,genomeLen[i])
        sheetVelvetSta.write(i+1,7,usedReads[i])
        sheetVelvetSta.write(i+1,8,totalReads[i])

    velvetSta.save(os.path.join(result_total_sta_path,'velvet_log_sta.xls'))

    contigN50_max_value_list = []  #all_samples
    contigN50_max_index_list = []
    for i in range(len(sample_list)):
        each_sample_contigN50_list = contigN50[kmerNum*i:kmerNum*(i+1)]
        #print each_sample_contigN50_list
        contigN50_max_value_list.append(max(each_sample_contigN50_list))
        contigN50_max_index_list.append(kmerMin+(each_sample_contigN50_list.index(max(each_sample_contigN50_list)))*kmerStep)

    print contigN50_max_value_list
    print contigN50_max_index_list

    print '=============blast:use the contigs.fa to map AGINST NT_viruses'
    #1.first:format the local nt fasta_database file
    #os.system('makeblastdb -dbtype nucl -input_type fasta -in %s' %(nt_path))  ###big memory : server

    #after makeblastdb,generate some other indexed files besides nt file,
    #use the indexed nt file(with some other indexed files) as the following input ref db '-db'.

    #2.second:blastn

    timer_blastn_start = time.clock()
    timer_blastn_start1 = time.time()

    #-outfmt 6
    for i in range(len(sample_list)):
        kmer = contigN50_max_index_list[i]
        print kmer
        eachKmerDir = '%s%s_%d'%(velvet_path,sample_list[i],kmer) #/public1/home/yefq/project_adenovirus/results/velvet_result/jiayan_31
        #print eachKmerDir
        eachContigFile = os.path.join(eachKmerDir,'contigs.fa')
        #print eachContigFile
        #estimate the contigs.fa is empty or not.
        if os.path.getsize(eachContigFile) != 0:
            eachBlastResult = os.path.join(blast_path,sample_list[i]+'.out')
            os.system('blastn -db %s -query %s -outfmt 6 -dust no -num_threads %d -perc_identity 80 -word_size 20 -max_target_seqs %d -evalue 0.0000001 -out %s'
                %(nt_viruses_family_file,eachContigFile,num_threads,max_target_seqs_nt_viruses,eachBlastResult))

    timer_blastn_end = time.clock()
    timer_blastn_end1 = time.time()
    blastn_time = timer_blastn_end - timer_blastn_start
    blastn_time1 = timer_blastn_end1 - timer_blastn_start1
    '''
    print '@@@@@----------timer:blastn_time_clock:'
    print blastn_time  #0.287237
    print '@@@@@----------timer:blastn_time_time:'
    print blastn_time1  #30.0745470524
    '''
    sheetPipelinerStaClock.write(1,2,blastn_time)
    sheetPipelinerStaTime.write(1,2,blastn_time1)

    ##3.third:sorted the blast result, use the different colomn
    for i in range(len(colNumSorted_list)):
        sort_blast_result_by_blastnColumn(colNumSorted_list[i])


    time_total_end = time.clock()
    time_total_end1 = time.time()
    

    total_time = time_total_end - time_total_start
    total_time1 = time_total_end1 - time_total_start1
    '''
    print '@@@@@-----------timer:total_time_clock:'
    print total_time  #23.519992
    print '@@@@@-----------timer:total_time_time:'
    print total_time1  #2177.6757071
    '''
    sheetPipelinerStaClock.write(1,3,total_time)
    sheetPipelinerStaTime.write(1,3,total_time1)

    pipelinerSta.save(os.path.join(result_total_sta_path,'pipeline_time_sta.xls'))

##method : withGenomeHost
def process_pipeliner_se_withGenomeHost(genomeHostName,ntVirusesFamilyName):
    ##rebuild : use all genome host
    genome_host_path = genomesDB_path+genomeHostName
    genome_host_file = os.path.join(genome_host_path,genomeHostName+'.fa')

    genome_host_bt2_index_base_name = 'bt2_index_%s'%(genomeHostName)
    
    ##rebuild : use ntVirusesFamilyName to replace nt_viruses_family_name_list[i], used for mapping with bowtie2 on nt_viruses_ and blastn
    family_path = 'nt_%s'%(ntVirusesFamilyName)
    family_file_name = 'nt_%s.fasta'%(ntVirusesFamilyName)
    
    nt_viruses_family_path = blastDB_path+family_path
    nt_viruses_family_file = os.path.join(nt_viruses_family_path,family_file_name)
        
    nt_virus_bt2_index_base_name = 'bt2_index_%s'%(ntVirusesFamilyName)

    nt_viruses_family_index_header_file_name = 'nt_%s_index_title'%(ntVirusesFamilyName)  #0_index_title / Adenoviridae_index_title
    nt_viruses_family_index_header_file = os.path.join(nt_viruses_family_path,nt_viruses_family_index_header_file_name)

    
    #print sample_se_list

    time_total_start = time.clock()
    time_total_start1 = time.time()

    pipelinerSta = xlwt.Workbook()
    sheetPipelinerStaClock = pipelinerSta.add_sheet('pipeline_clock_statistics',cell_overwrite_ok=True)
    sheetPipelinerStaTime = pipelinerSta.add_sheet('pipeline_time_statistics',cell_overwrite_ok=True)
    
    rowTitle = ['map_%s'%ntVirusesFamilyName,'assembly_velvet','map_blastn','total_time']
    for i in range(0,len(rowTitle)):
        sheetPipelinerStaClock.write(0,i,rowTitle[i]) #first row
        sheetPipelinerStaTime.write(0,i,rowTitle[i])



    print '==========map and align : AGINST NT_viruses reference,keep mapped reads'

    #bowtie2IndexNTVirusFamily = nt_viruses_family_name_list[0]
    #bowtie2IndexNTVirusFamily = bt2_index_base_name_list[0]  #'bt2_index_Adenoviridae'

    timer_bowtie2_nt_virus_start = time.clock()
    timer_bowtie2_nt_virus_start1 = time.time()

    afterNTVIRUSfastqFiles =[]
    for i in range(len(sample_se_list)):
        outputSAMFileNTVirus = bowtie2_path+sample_se_list[i]+map_nt_viruses_singleAndunpair_sam_file_suffix
        outputFASTQFileNTVirus = bowtie2_path+sample_se_list[i]+map_nt_viruses_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -U %s --no-head --no-unal -S %s --al %s -p %d'
            %(nt_viruses_family_path,nt_virus_bt2_index_base_name,samplesRawDataSingleFiles[i],outputSAMFileNTVirus,outputFASTQFileNTVirus,num_threads))

        afterNTVIRUSfastqFiles.append(outputFASTQFileNTVirus)


    timer_bowtie2_nt_virus_end = time.clock()
    timer_bowtie2_nt_virus_end1 = time.time()
    bowtie2_nt_virus_time = timer_bowtie2_nt_virus_end - timer_bowtie2_nt_virus_start
    bowtie2_nt_virus_time1 = timer_bowtie2_nt_virus_end1 - timer_bowtie2_nt_virus_start1
    '''
    print '@@@@@----------timer:bowtie2_nt_virus_time_clock:'
    print bowtie2_nt_virus_time    #18.209088
    print '@@@@@----------timer:bowtie2_nt_virus_time_time:'
    print bowtie2_nt_virus_time1    #1678.66673303
    '''
    sheetPipelinerStaClock.write(1,0,bowtie2_nt_virus_time)
    sheetPipelinerStaTime.write(1,0,bowtie2_nt_virus_time1)


    print '==========assembly : velvet with alignment data above'

    timer_velvet_start = time.clock()
    timer_velvet_start1 = time.time()

    minContigLength = 100
    for i in range(len(sample_se_list)):
        for kmer in range(kmerMin,kmerMax,kmerStep):
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_se_list[i],kmer)
            os.system('velveth %s %d -short -fastq %s' %(eachKmerDir,kmer,afterNTVIRUSfastqFiles[i]))
            os.system('velvetg %s -exp_cov auto -cov_cutoff auto -min_contig_lgth %d' %(eachKmerDir,minContigLength))

            eachKmerFiles = os.listdir(eachKmerDir)
            eachKmerFilesAP = []
            for j in range(len(eachKmerFiles)):
                eachKmerFilesAP.append(eachKmerDir+'/'+eachKmerFiles[j])
            for j in range(len(eachKmerFilesAP)):
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Graph2' : os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/LastGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/PreGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Roadmaps': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Sequences': os.remove(eachKmerFilesAP[j])

    timer_velvet_end = time.clock()
    timer_velvet_end1 = time.time()
    velvet_time = timer_velvet_end - timer_velvet_start
    velvet_time1 = timer_velvet_end1 - timer_velvet_start1
    '''
    print '@@@@@----------timer:velvet_time_clock:'
    print velvet_time    #0.034675
    print '@@@@@----------timer:velvet_time_time:'
    print velvet_time1    #1.79084706306
    '''
    sheetPipelinerStaClock.write(1,1,velvet_time)
    sheetPipelinerStaTime.write(1,1,velvet_time1)

    ##================statistics velvet : pick up the proper kmer , accoding to N50 max========
    velvetSta = xlwt.Workbook()
    sheetVelvetSta = velvetSta.add_sheet('velvet_statistics', cell_overwrite_ok=True)
    kmerNum = (kmerMax-kmerMin)/2  #(35-31)/2=2
    row0 = ['sampleNums','kmer','medianCovDepth','nodesNum','contigN50','contigMax','genomeLen','usedReads','totalReads']
    for i in range(0,len(row0)):
        sheetVelvetSta.write(0,i,row0[i])  #first row

    #first cow is sample_se_list
    #second cow is kmer (kmerMin to kmerMax , kmerStep)

    for i in range(0,len(sample_se_list)):
        sheetVelvetSta.write(kmerNum*i+1,0,sample_se_list[i])   #first column is the sampleNums

    colKmer = range(kmerMin,kmerMax,kmerStep)
    for i in range(0,len(sample_se_list)*kmerNum,kmerNum): #i : every column
        for k in range(0,kmerNum):
            sheetVelvetSta.write(i+k+1,1,colKmer[k])   #second column is the k-mer
            k += 1

    logLast2Str = []
    MedianCovDepth = []

    logLast1Str = []
    nodesNum = []
    contigN50 = []
    contigMax = []
    genomeLen = []
    usedReads = []
    totalReads = []
    for i in range(len(sample_se_list)):
        for kmer in range(kmerMin,kmerMax,kmerStep):
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_se_list[i],kmer)
            Log_file = os.path.join(eachKmerDir,'Log')
            logs = open(Log_file,'r').readlines()
            #logLast2line = logs[len(logs)-2]
            #if logLast2line != '':
            #Median coverage depth = 96.434783  #not the logs[26] , to avoid the task again with adding the log_information to Log File
            #logLast2Str = logLast2line.split(' ')
            #print logLast2Str
            #print logLast2Str[4][:-2]
            #MedianCovDepth.append(string.atof(logLast2Str[4][:-2]))
            logLast1line = logs[len(logs)-1]  #Final graph has 199 nodes and n50 of 47325, max 166302, total 1570595, using 1470371/1487666 reads
            logLast1Str = logLast1line.split(' ')
            #['Final', 'graph', 'has', '172', 'nodes', 'and', 'n50', 'of', '46917,', 'max', '237233,', 'total', '1570012,', 'using', '1478778/1487666', 'reads\n']
            #print logLast1Str
            nodesNum.append(string.atoi(logLast1Str[3], 10))
            contigN50.append(string.atoi(logLast1Str[8][:-1],10))
            contigMax.append(string.atoi(logLast1Str[10][:-1],10))
            genomeLen.append(string.atoi(logLast1Str[12][:-1],10))
            usedReads.append(string.atoi(logLast1Str[14].split('/')[0], 10))
            totalReads.append(string.atoi(logLast1Str[14].split('/')[1], 10))
    #print contigN50
    for i in range(0,len(sample_se_list)*kmerNum):
        #sheetVelvetSta.write(i+1,2,MedianCovDepth[i])
        sheetVelvetSta.write(i+1,3,nodesNum[i])
        sheetVelvetSta.write(i+1,4,contigN50[i])
        sheetVelvetSta.write(i+1,5,contigMax[i])
        sheetVelvetSta.write(i+1,6,genomeLen[i])
        sheetVelvetSta.write(i+1,7,usedReads[i])
        sheetVelvetSta.write(i+1,8,totalReads[i])

    velvetSta.save(os.path.join(result_total_sta_path,'velvet_log_sta.xls'))

    ##check the max contigN50 to pick up ,then blastn
    contigN50_max_value_list = []  #all_samples
    contigN50_max_index_list = []
    for i in range(len(sample_se_list)):
        each_sample_contigN50_list = contigN50[kmerNum*i:kmerNum*(i+1)]
        #print each_sample_contigN50_list
        contigN50_max_value_list.append(max(each_sample_contigN50_list))
        contigN50_max_index_list.append(kmerMin+(each_sample_contigN50_list.index(max(each_sample_contigN50_list)))*kmerStep)

    #print contigN50_max_value_list
    #print contigN50_max_index_list

    print '=============blast:use the contigs.fa to map AGINST NT_viruses'
    #1.first:format the local nt fasta_database file
    #os.system('makeblastdb -dbtype nucl -input_type fasta -in %s' %(nt_path))  ###big memory : server

    #after makeblastdb,generate some other indexed files besides nt file,
    #use the indexed nt file(with some other indexed files) as the following input ref db '-db'.

    #2.second:blastn

    timer_blastn_start = time.clock()
    timer_blastn_start1 = time.time()

                 
    for i in range(len(sample_se_list)):
        kmer = contigN50_max_index_list[i]
        eachKmerDir = '%s%s_%d'%(velvet_path,sample_se_list[i],kmer) #/public1/home/yefq/project_adenovirus/results/velvet_result/jiayan_31
        #print eachKmerDir
        eachContigFile = os.path.join(eachKmerDir,'contigs.fa')
        #print eachContigFile
        if os.path.getsize(eachContigFile) != 0:
            eachBlastResult = os.path.join(blast_path,sample_se_list[i]+'.out')
            os.system('blastn -db %s -query %s -outfmt 6 -dust no -num_threads %d -perc_identity 80 -word_size 20 -max_target_seqs %d -evalue 0.0000001 -out %s'
                %(nt_viruses_family_file,eachContigFile,num_threads,max_target_seqs_nt_viruses,eachBlastResult))

    timer_blastn_end = time.clock()
    timer_blastn_end1 = time.time()
    blastn_time = timer_blastn_end - timer_blastn_start
    blastn_time1 = timer_blastn_end1 - timer_blastn_start1
    '''
    print '@@@@@----------timer:blastn_time_clock:'
    print blastn_time  #0.287237
    print '@@@@@----------timer:blastn_time_time:'
    print blastn_time1  #30.0745470524
    '''
    sheetPipelinerStaClock.write(1,2,blastn_time)
    sheetPipelinerStaTime.write(1,2,blastn_time1)

    ##3.third:sorted the blast result, use the different colomn
    for i in range(len(colNumSorted_list)):
        sort_blast_result_by_blastnColumn(colNumSorted_list[i])


    time_total_end = time.clock()
    time_total_end1 = time.time()
    

    total_time = time_total_end - time_total_start
    total_time1 = time_total_end1 - time_total_start1
    '''
    print '@@@@@-----------timer:total_time_clock:'
    print total_time  #23.519992
    print '@@@@@-----------timer:total_time_time:'
    print total_time1  #2177.6757071
    '''
    sheetPipelinerStaClock.write(1,3,total_time)
    sheetPipelinerStaTime.write(1,3,total_time1)

    pipelinerSta.save(os.path.join(result_total_sta_path,'pipeline_time_sta.xls'))



##method : justMapAginstVirusWithoutQC
def process_pipeliner_pe_justMapAginstVirus(genomeHostName,ntVirusesFamilyName):
    ##rebuild : use all genome host
    genome_host_path = genomesDB_path+genomeHostName
    genome_host_file = os.path.join(genome_host_path,genomeHostName+'.fa')

    genome_host_bt2_index_base_name = 'bt2_index_%s'%(genomeHostName)
    
    ##rebuild : use ntVirusesFamilyName to replace nt_viruses_family_name_list[i], used for mapping with bowtie2 on nt_viruses_ and blastn
    family_path = 'nt_%s'%(ntVirusesFamilyName)
    family_file_name = 'nt_%s.fasta'%(ntVirusesFamilyName)
    
    nt_viruses_family_path = blastDB_path+family_path
    nt_viruses_family_file = os.path.join(nt_viruses_family_path,family_file_name)
        
    nt_virus_bt2_index_base_name = 'bt2_index_%s'%(ntVirusesFamilyName)

    nt_viruses_family_index_header_file_name = 'nt_%s_index_title'%(ntVirusesFamilyName)  #0_index_title / Adenoviridae_index_title
    nt_viruses_family_index_header_file = os.path.join(nt_viruses_family_path,nt_viruses_family_index_header_file_name)

    #print sample_list

    time_total_start = time.clock()
    time_total_start1 = time.time()

    pipelinerSta = xlwt.Workbook()
    sheetPipelinerStaClock = pipelinerSta.add_sheet('pipeline_clock_statistics',cell_overwrite_ok=True)
    sheetPipelinerStaTime = pipelinerSta.add_sheet('pipeline_time_statistics',cell_overwrite_ok=True)
    #rowTitle = ['qc_trimmomatic','map_hp38','map_genebank_virus','map_nt_virus','assembly_velvet','map_blastn','total_time']
    #rowTitle = ['qc_trimmomatic','map_hp38','map_nt_virus','assembly_velvet','map_blastn','total_time']
    rowTitle = ['qc_trimmomatic','map_%s'%ntVirusesFamilyName,'assembly_velvet','map_blastn','total_time']
    for i in range(0,len(rowTitle)):
        sheetPipelinerStaClock.write(0,i,rowTitle[i]) #first row
        sheetPipelinerStaTime.write(0,i,rowTitle[i])


    
    print '==========qc:trimmomatic'
    timer_trimmomatic_start = time.clock()
    timer_trimmomatic_start1 = time.time()

    trimmomatic_results_prefix = []
    for i in range(len(sample_list)):
        trimmomatic_result_prefix = trimmomatic_path+sample_list[i]
        trimmomatic_results_prefix.append(trimmomatic_result_prefix)
    #print trimmomatic_results_prefix ##/public1/home/yefq/project_adenovirus/results/trimmomatic_result/jiayan

    outputFPFiles = []
    outputFUFiles = []
    outputRPFiles = []
    outputRUFiles = []
    for i in range(len(sample_list)):
        #print sample_list[i]
        trimlogFile=trimmomatic_results_prefix[i]+'_trimlog.txt'
        outputFPFile=trimmomatic_results_prefix[i]+'_output_forward_paired.fastq'
        outputFUFile=trimmomatic_results_prefix[i]+'_output_forward_unpaired.fastq'
        outputRPFile=trimmomatic_results_prefix[i]+'_output_reverse_paired.fastq'
        outputRUFile=trimmomatic_results_prefix[i]+'_output_reverse_unpaired.fastq'

        os.system('java -jar /usr/share/java/trimmomatic.jar PE -threads %d -trimlog %s %s %s %s %s %s %s ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'
        %(num_threads,trimlogFile,samplesRawDataForwardFiles[i],samplesRawDataReverseFiles[i],outputFPFile,outputFUFile,outputRPFile,outputRUFile))

        outputFPFiles.append(outputFPFile)
        outputRPFiles.append(outputRPFile)
        outputFUFiles.append(outputFUFile)
        outputRUFiles.append(outputRUFile)

    timer_trimmomatic_end = time.clock()
    timer_trimmomatic_end1 = time.time()
    trimmomatic_time = timer_trimmomatic_end - timer_trimmomatic_start
    trimmomatic_time1 = timer_trimmomatic_end1 - timer_trimmomatic_start1
    '''
    print '@@@@@----------timer:trimmomatic_time_clock:'
    print trimmomatic_time  #4.822289
    print '@@@@@----------timer:trimmomatic_time_time:'
    print trimmomatic_time1  #466.976895094
    '''
    sheetPipelinerStaClock.write(1,0,trimmomatic_time)
    sheetPipelinerStaTime.write(1,0,trimmomatic_time1)


    print '==========map and align : AGINST NT_viruses reference,keep mapped reads'

    #bowtie2IndexNTVirusFamily = nt_viruses_family_name_list[0]
    #bowtie2IndexNTVirusFamily = bt2_index_base_name_list[0]  #'bt2_index_Adenoviridae'

    timer_bowtie2_nt_virus_start = time.clock()
    timer_bowtie2_nt_virus_start1 = time.time()

    outputFastqPair1FileNTVirus_list = []
    outputFastqPair2FileNTVirus_list = []
    outputFastqUnpairFileNTVirus_list = []
    for i in range(len(sample_list)):
        outputSAMFileNTVirus = bowtie2_path+sample_list[i]+map_nt_viruses_singleAndunpair_sam_file_suffix
        outputFastqPairFileNTVirus = bowtie2_path+sample_list[i]+'_pair_VIRUS.fastq'
        outputFastqPair1FileNTVirus = bowtie2_path+sample_list[i]+'_pair_VIRUS.1.fastq'
        outputFastqPair2FileNTVirus = bowtie2_path+sample_list[i]+'_pair_VIRUS.2.fastq'
        
        outputFastqUnpairFileNTVirus = bowtie2_path+sample_list[i]+map_nt_viruses_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -1 %s -2 %s -U %s,%s --no-head --no-unal -S %s --al-conc %s --al %s -p %d'
            %(nt_viruses_family_path,nt_virus_bt2_index_base_name,outputFPFiles[i],outputRPFiles[i],outputFUFiles[i],outputRUFiles[i],outputSAMFileNTVirus,outputFastqPairFileNTVirus,outputFastqUnpairFileNTVirus,num_threads))

        outputFastqPair1FileNTVirus_list.append(outputFastqPair1FileNTVirus)
        outputFastqPair2FileNTVirus_list.append(outputFastqPair2FileNTVirus)
        outputFastqUnpairFileNTVirus_list.append(outputFastqUnpairFileNTVirus)
        
    timer_bowtie2_nt_virus_end = time.clock()
    timer_bowtie2_nt_virus_end1 = time.time()
    bowtie2_nt_virus_time = timer_bowtie2_nt_virus_end - timer_bowtie2_nt_virus_start
    bowtie2_nt_virus_time1 = timer_bowtie2_nt_virus_end1 - timer_bowtie2_nt_virus_start1
    '''
    print '@@@@@----------timer:bowtie2_nt_virus_time_clock:'
    print bowtie2_nt_virus_time    #18.209088
    print '@@@@@----------timer:bowtie2_nt_virus_time_time:'
    print bowtie2_nt_virus_time1    #1678.66673303
    '''
    sheetPipelinerStaClock.write(1,1,bowtie2_nt_virus_time)
    sheetPipelinerStaTime.write(1,1,bowtie2_nt_virus_time1)

    print '==========assembly : velvet with alignment data above'

    timer_velvet_start = time.clock()
    timer_velvet_start1 = time.time()

    for i in range(len(sample_list)):
        #inputPE = velvet_result+sample_list[i]+'_pe.fastq'
        #os.system('%s %s %s %s' %(shuffleSequences_fastq_path,outputFastq1AfterHGandVirus[i],outputFastq2AfterHGandVirus[i],inputPE))
        for kmer in range(kmerMin,kmerMax,kmerStep):
            #print sample_list[i]+'************************************'
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_list[i],kmer)
            os.system('velveth %s %d -shortPaired -fastq -separate %s %s -short -fastq %s'%(eachKmerDir,kmer,outputFastqPair1FileNTVirus_list[i],outputFastqPair2FileNTVirus_list[i],outputFastqUnpairFileNTVirus_list[i]))
            os.system('velvetg %s -exp_cov auto -cov_cutoff auto ' %(eachKmerDir))


            eachKmerFiles = os.listdir(eachKmerDir)
            eachKmerFilesAP = []
            for j in range(len(eachKmerFiles)):
                eachKmerFilesAP.append(eachKmerDir+'/'+eachKmerFiles[j])
            for j in range(len(eachKmerFilesAP)):
                #pass
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Graph2' : os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/LastGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/PreGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Roadmaps': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Sequences': os.remove(eachKmerFilesAP[j])

    timer_velvet_end = time.clock()
    timer_velvet_end1 = time.time()
    velvet_time = timer_velvet_end - timer_velvet_start
    velvet_time1 = timer_velvet_end1 - timer_velvet_start1
    '''
    print '@@@@@----------timer:velvet_time_clock:'
    print velvet_time    #0.034675
    print '@@@@@----------timer:velvet_time_time:'
    print velvet_time1    #1.79084706306
    '''
    sheetPipelinerStaClock.write(1,2,velvet_time)
    sheetPipelinerStaTime.write(1,2,velvet_time1)

    ##================statistics velvet : pick up the proper kmer , accoding to N50 max========
    velvetSta = xlwt.Workbook()
    sheetVelvetSta = velvetSta.add_sheet('velvet_statistics', cell_overwrite_ok=True)
    kmerNum = (kmerMax-kmerMin)/2  #(35-31)/2=2
    row0 = ['sampleNums','kmer','medianCovDepth','nodesNum','contigN50','contigMax','genomeLen','usedReads','totalReads']
    for i in range(0,len(row0)):
        sheetVelvetSta.write(0,i,row0[i])  #first row

    #first cow is sample_list
    #second cow is kmer (kmerMin to kmerMax , kmerStep)

    for i in range(0,len(sample_list)):
        sheetVelvetSta.write(kmerNum*i+1,0,sample_list[i])   #first column is the sampleNums

    colKmer = range(kmerMin,kmerMax,kmerStep)
    for i in range(0,len(sample_list)*kmerNum,kmerNum): #i : every column
        for k in range(0,kmerNum):
            sheetVelvetSta.write(i+k+1,1,colKmer[k])   #second column is the k-mer
            k += 1

    logLast2Str = []
    MedianCovDepth = []

    logLast1Str = []
    nodesNum = []
    contigN50 = []
    contigMax = []
    genomeLen = []
    usedReads = []
    totalReads = []
    for i in range(len(sample_list)):
        for kmer in range(kmerMin,kmerMax,kmerStep):
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_list[i],kmer)
            Log_file = os.path.join(eachKmerDir,'Log')
            logs = open(Log_file,'r').readlines()
            #logLast2line = logs[len(logs)-2]
            #if logLast2line != '':
            #Median coverage depth = 96.434783  #not the logs[26] , to avoid the task again with adding the log_information to Log File
            #logLast2Str = logLast2line.split(' ')
            #print logLast2Str
            #print logLast2Str[4][:-2]
            #MedianCovDepth.append(string.atof(logLast2Str[4][:-2]))
            logLast1line = logs[len(logs)-1]  #Final graph has 199 nodes and n50 of 47325, max 166302, total 1570595, using 1470371/1487666 reads
            logLast1Str = logLast1line.split(' ')
            #['Final', 'graph', 'has', '172', 'nodes', 'and', 'n50', 'of', '46917,', 'max', '237233,', 'total', '1570012,', 'using', '1478778/1487666', 'reads\n']
            #print logLast1Str
            nodesNum.append(string.atoi(logLast1Str[3], 10))
            contigN50.append(string.atoi(logLast1Str[8][:-1],10))
            contigMax.append(string.atoi(logLast1Str[10][:-1],10))
            genomeLen.append(string.atoi(logLast1Str[12][:-1],10))
            usedReads.append(string.atoi(logLast1Str[14].split('/')[0], 10))
            totalReads.append(string.atoi(logLast1Str[14].split('/')[1], 10))
    #print contigN50
    for i in range(0,len(sample_list)*kmerNum):
        #sheetVelvetSta.write(i+1,2,MedianCovDepth[i])
        sheetVelvetSta.write(i+1,3,nodesNum[i])
        sheetVelvetSta.write(i+1,4,contigN50[i])
        sheetVelvetSta.write(i+1,5,contigMax[i])
        sheetVelvetSta.write(i+1,6,genomeLen[i])
        sheetVelvetSta.write(i+1,7,usedReads[i])
        sheetVelvetSta.write(i+1,8,totalReads[i])

    velvetSta.save(os.path.join(result_total_sta_path,'velvet_log_sta.xls'))

    contigN50_max_value_list = []  #all_samples
    contigN50_max_index_list = []
    for i in range(len(sample_list)):
        each_sample_contigN50_list = contigN50[kmerNum*i:kmerNum*(i+1)]
        #print each_sample_contigN50_list
        contigN50_max_value_list.append(max(each_sample_contigN50_list))
        contigN50_max_index_list.append(kmerMin+(each_sample_contigN50_list.index(max(each_sample_contigN50_list)))*kmerStep)

    #print contigN50_max_value_list
    #print contigN50_max_index_list

    print '=============blast:use the contigs.fa to map AGINST NT_viruses'
    #1.first:format the local nt fasta_database file
    #os.system('makeblastdb -dbtype nucl -input_type fasta -in %s' %(nt_path))  ###big memory : server

    #after makeblastdb,generate some other indexed files besides nt file,
    #use the indexed nt file(with some other indexed files) as the following input ref db '-db'.

    #2.second:blastn

    timer_blastn_start = time.clock()
    timer_blastn_start1 = time.time()

    #-outfmt 6
    for i in range(len(sample_list)):
        kmer = contigN50_max_index_list[i]
        eachKmerDir = '%s%s_%d'%(velvet_path,sample_list[i],kmer) #/public1/home/yefq/project_adenovirus/results/velvet_result/jiayan_31
        #print eachKmerDir
        eachContigFile = os.path.join(eachKmerDir,'contigs.fa')
        #print eachContigFile
        #estimate the contigs.fa is empty or not.
        if os.path.getsize(eachContigFile) != 0:
            eachBlastResult = os.path.join(blast_path,sample_list[i]+'.out')
            os.system('blastn -db %s -query %s -outfmt 6 -dust no -num_threads %d -perc_identity 80 -word_size 20 -max_target_seqs %d -evalue 0.0000001 -out %s'
                %(nt_viruses_family_file,eachContigFile,num_threads,max_target_seqs_nt_viruses,eachBlastResult))

    timer_blastn_end = time.clock()
    timer_blastn_end1 = time.time()
    blastn_time = timer_blastn_end - timer_blastn_start
    blastn_time1 = timer_blastn_end1 - timer_blastn_start1
    '''
    print '@@@@@----------timer:blastn_time_clock:'
    print blastn_time  #0.287237
    print '@@@@@----------timer:blastn_time_time:'
    print blastn_time1  #30.0745470524
    '''
    sheetPipelinerStaClock.write(1,3,blastn_time)
    sheetPipelinerStaTime.write(1,3,blastn_time1)

    ##3.third:sorted the blast result, use the different colomn
    for i in range(len(colNumSorted_list)):
        sort_blast_result_by_blastnColumn(colNumSorted_list[i])


    time_total_end = time.clock()
    time_total_end1 = time.time()
    

    total_time = time_total_end - time_total_start
    total_time1 = time_total_end1 - time_total_start1
    '''
    print '@@@@@-----------timer:total_time_clock:'
    print total_time  #23.519992
    print '@@@@@-----------timer:total_time_time:'
    print total_time1  #2177.6757071
    '''
    sheetPipelinerStaClock.write(1,4,total_time)
    sheetPipelinerStaTime.write(1,4,total_time1)

    pipelinerSta.save(os.path.join(result_total_sta_path,'pipeline_time_sta.xls'))

##method : justMapAginstVirusWithoutQC
def process_pipeliner_se_justMapAginstVirus(genomeHostName,ntVirusesFamilyName):
    ##rebuild : use all genome host
    genome_host_path = genomesDB_path+genomeHostName
    genome_host_file = os.path.join(genome_host_path,genomeHostName+'.fa')

    genome_host_bt2_index_base_name = 'bt2_index_%s'%(genomeHostName)
    
    ##rebuild : use ntVirusesFamilyName to replace nt_viruses_family_name_list[i], used for mapping with bowtie2 on nt_viruses_ and blastn
    family_path = 'nt_%s'%(ntVirusesFamilyName)
    family_file_name = 'nt_%s.fasta'%(ntVirusesFamilyName)
    
    nt_viruses_family_path = blastDB_path+family_path
    nt_viruses_family_file = os.path.join(nt_viruses_family_path,family_file_name)
        
    nt_virus_bt2_index_base_name = 'bt2_index_%s'%(ntVirusesFamilyName)

    nt_viruses_family_index_header_file_name = 'nt_%s_index_title'%(ntVirusesFamilyName)  #0_index_title / Adenoviridae_index_title
    nt_viruses_family_index_header_file = os.path.join(nt_viruses_family_path,nt_viruses_family_index_header_file_name)

    
    #print sample_se_list

    time_total_start = time.clock()
    time_total_start1 = time.time()

    pipelinerSta = xlwt.Workbook()
    sheetPipelinerStaClock = pipelinerSta.add_sheet('pipeline_clock_statistics',cell_overwrite_ok=True)
    sheetPipelinerStaTime = pipelinerSta.add_sheet('pipeline_time_statistics',cell_overwrite_ok=True)
    
    rowTitle = ['qc_trimmomatic','map_%s'%ntVirusesFamilyName,'map_%s'%genomeHostName,'assembly_velvet','map_blastn','total_time']
    for i in range(0,len(rowTitle)):
        sheetPipelinerStaClock.write(0,i,rowTitle[i]) #first row
        sheetPipelinerStaTime.write(0,i,rowTitle[i])

    
    print '==========qc:trimmomatic'
    timer_trimmomatic_start = time.clock()
    timer_trimmomatic_start1 = time.time()
    
    MINLEN = 10
    LEADING = 0
    TRAILING = 0    
    outputSingleFiles = []
    for i in range(len(sample_se_list)):
        trimlogFile=trimmomatic_path+sample_se_list[i]+'_trimlog.txt'
        outputSingleFile=trimmomatic_path+sample_se_list[i]+'_output.fastq'

        os.system('java -jar /usr/share/java/trimmomatic.jar SE -threads %d -trimlog %s %s %s ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-SE.fa:2:30:10 LEADING:%d TRAILING:%d SLIDINGWINDOW:4:15 MINLEN:%d'
        %(num_threads,trimlogFile,samplesRawDataSingleFiles[i],outputSingleFile,LEADING,TRAILING,MINLEN))

        outputSingleFiles.append(outputSingleFile)

    timer_trimmomatic_end = time.clock()
    timer_trimmomatic_end1 = time.time()
    trimmomatic_time = timer_trimmomatic_end - timer_trimmomatic_start
    trimmomatic_time1 = timer_trimmomatic_end1 - timer_trimmomatic_start1
    '''
    print '@@@@@----------timer:trimmomatic_time_clock:'
    print trimmomatic_time  #4.822289
    print '@@@@@----------timer:trimmomatic_time_time:'
    print trimmomatic_time1  #466.976895094
    '''
    sheetPipelinerStaClock.write(1,0,trimmomatic_time)
    sheetPipelinerStaTime.write(1,0,trimmomatic_time1)



    print '==========map and align : AGINST NT_viruses reference,keep mapped reads'

    #bowtie2IndexNTVirusFamily = nt_viruses_family_name_list[0]
    #bowtie2IndexNTVirusFamily = bt2_index_base_name_list[0]  #'bt2_index_Adenoviridae'

    timer_bowtie2_nt_virus_start = time.clock()
    timer_bowtie2_nt_virus_start1 = time.time()

    afterNTVIRUSfastqFiles =[]
    for i in range(len(sample_se_list)):
        outputSAMFileNTVirus = bowtie2_path+sample_se_list[i]+map_nt_viruses_singleAndunpair_sam_file_suffix
        outputFASTQFileNTVirus = bowtie2_path+sample_se_list[i]+map_nt_viruses_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -U %s --no-head --no-unal -S %s --al %s -p %d'
            %(nt_viruses_family_path,nt_virus_bt2_index_base_name,outputSingleFiles[i],outputSAMFileNTVirus,outputFASTQFileNTVirus,num_threads))

        afterNTVIRUSfastqFiles.append(outputFASTQFileNTVirus)


    timer_bowtie2_nt_virus_end = time.clock()
    timer_bowtie2_nt_virus_end1 = time.time()
    bowtie2_nt_virus_time = timer_bowtie2_nt_virus_end - timer_bowtie2_nt_virus_start
    bowtie2_nt_virus_time1 = timer_bowtie2_nt_virus_end1 - timer_bowtie2_nt_virus_start1
    '''
    print '@@@@@----------timer:bowtie2_nt_virus_time_clock:'
    print bowtie2_nt_virus_time    #18.209088
    print '@@@@@----------timer:bowtie2_nt_virus_time_time:'
    print bowtie2_nt_virus_time1    #1678.66673303
    '''
    sheetPipelinerStaClock.write(1,1,bowtie2_nt_virus_time)
    sheetPipelinerStaTime.write(1,1,bowtie2_nt_virus_time1)


    print '==========assembly : velvet with alignment data above'

    timer_velvet_start = time.clock()
    timer_velvet_start1 = time.time()

    minContigLength = 100
    for i in range(len(sample_se_list)):
        for kmer in range(kmerMin,kmerMax,kmerStep):
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_se_list[i],kmer)
            os.system('velveth %s %d -short -fastq %s' %(eachKmerDir,kmer,afterNTVIRUSfastqFiles[i]))
            os.system('velvetg %s -exp_cov auto -cov_cutoff auto -min_contig_lgth %d' %(eachKmerDir,minContigLength))

            eachKmerFiles = os.listdir(eachKmerDir)
            eachKmerFilesAP = []
            for j in range(len(eachKmerFiles)):
                eachKmerFilesAP.append(eachKmerDir+'/'+eachKmerFiles[j])
            for j in range(len(eachKmerFilesAP)):
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Graph2' : os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/LastGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/PreGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Roadmaps': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Sequences': os.remove(eachKmerFilesAP[j])

    timer_velvet_end = time.clock()
    timer_velvet_end1 = time.time()
    velvet_time = timer_velvet_end - timer_velvet_start
    velvet_time1 = timer_velvet_end1 - timer_velvet_start1
    '''
    print '@@@@@----------timer:velvet_time_clock:'
    print velvet_time    #0.034675
    print '@@@@@----------timer:velvet_time_time:'
    print velvet_time1    #1.79084706306
    '''
    sheetPipelinerStaClock.write(1,2,velvet_time)
    sheetPipelinerStaTime.write(1,2,velvet_time1)

    ##================statistics velvet : pick up the proper kmer , accoding to N50 max========
    velvetSta = xlwt.Workbook()
    sheetVelvetSta = velvetSta.add_sheet('velvet_statistics', cell_overwrite_ok=True)
    kmerNum = (kmerMax-kmerMin)/2  #(35-31)/2=2
    row0 = ['sampleNums','kmer','medianCovDepth','nodesNum','contigN50','contigMax','genomeLen','usedReads','totalReads']
    for i in range(0,len(row0)):
        sheetVelvetSta.write(0,i,row0[i])  #first row

    #first cow is sample_se_list
    #second cow is kmer (kmerMin to kmerMax , kmerStep)

    for i in range(0,len(sample_se_list)):
        sheetVelvetSta.write(kmerNum*i+1,0,sample_se_list[i])   #first column is the sampleNums

    colKmer = range(kmerMin,kmerMax,kmerStep)
    for i in range(0,len(sample_se_list)*kmerNum,kmerNum): #i : every column
        for k in range(0,kmerNum):
            sheetVelvetSta.write(i+k+1,1,colKmer[k])   #second column is the k-mer
            k += 1

    logLast2Str = []
    MedianCovDepth = []

    logLast1Str = []
    nodesNum = []
    contigN50 = []
    contigMax = []
    genomeLen = []
    usedReads = []
    totalReads = []
    for i in range(len(sample_se_list)):
        for kmer in range(kmerMin,kmerMax,kmerStep):
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_se_list[i],kmer)
            Log_file = os.path.join(eachKmerDir,'Log')
            logs = open(Log_file,'r').readlines()
            #logLast2line = logs[len(logs)-2]
            #if logLast2line != '':
            #Median coverage depth = 96.434783  #not the logs[26] , to avoid the task again with adding the log_information to Log File
            #logLast2Str = logLast2line.split(' ')
            #print logLast2Str
            #print logLast2Str[4][:-2]
            #MedianCovDepth.append(string.atof(logLast2Str[4][:-2]))
            logLast1line = logs[len(logs)-1]  #Final graph has 199 nodes and n50 of 47325, max 166302, total 1570595, using 1470371/1487666 reads
            logLast1Str = logLast1line.split(' ')
            #['Final', 'graph', 'has', '172', 'nodes', 'and', 'n50', 'of', '46917,', 'max', '237233,', 'total', '1570012,', 'using', '1478778/1487666', 'reads\n']
            #print logLast1Str
            nodesNum.append(string.atoi(logLast1Str[3], 10))
            contigN50.append(string.atoi(logLast1Str[8][:-1],10))
            contigMax.append(string.atoi(logLast1Str[10][:-1],10))
            genomeLen.append(string.atoi(logLast1Str[12][:-1],10))
            usedReads.append(string.atoi(logLast1Str[14].split('/')[0], 10))
            totalReads.append(string.atoi(logLast1Str[14].split('/')[1], 10))
    #print contigN50
    for i in range(0,len(sample_se_list)*kmerNum):
        #sheetVelvetSta.write(i+1,2,MedianCovDepth[i])
        sheetVelvetSta.write(i+1,3,nodesNum[i])
        sheetVelvetSta.write(i+1,4,contigN50[i])
        sheetVelvetSta.write(i+1,5,contigMax[i])
        sheetVelvetSta.write(i+1,6,genomeLen[i])
        sheetVelvetSta.write(i+1,7,usedReads[i])
        sheetVelvetSta.write(i+1,8,totalReads[i])

    velvetSta.save(os.path.join(result_total_sta_path,'velvet_log_sta.xls'))

    ##check the max contigN50 to pick up ,then blastn
    contigN50_max_value_list = []  #all_samples
    contigN50_max_index_list = []
    for i in range(len(sample_se_list)):
        each_sample_contigN50_list = contigN50[kmerNum*i:kmerNum*(i+1)]
        #print each_sample_contigN50_list
        contigN50_max_value_list.append(max(each_sample_contigN50_list))
        contigN50_max_index_list.append(kmerMin+(each_sample_contigN50_list.index(max(each_sample_contigN50_list)))*kmerStep)

    #print contigN50_max_value_list
    #print contigN50_max_index_list

    print '=============blast:use the contigs.fa to map AGINST NT_viruses'
    #1.first:format the local nt fasta_database file
    #os.system('makeblastdb -dbtype nucl -input_type fasta -in %s' %(nt_path))  ###big memory : server

    #after makeblastdb,generate some other indexed files besides nt file,
    #use the indexed nt file(with some other indexed files) as the following input ref db '-db'.

    #2.second:blastn

    timer_blastn_start = time.clock()
    timer_blastn_start1 = time.time()

                 
    for i in range(len(sample_se_list)):
        kmer = contigN50_max_index_list[i]
        eachKmerDir = '%s%s_%d'%(velvet_path,sample_se_list[i],kmer) #/public1/home/yefq/project_adenovirus/results/velvet_result/jiayan_31
        #print eachKmerDir
        eachContigFile = os.path.join(eachKmerDir,'contigs.fa')
        #print eachContigFile
        if os.path.getsize(eachContigFile) != 0:
            eachBlastResult = os.path.join(blast_path,sample_se_list[i]+'.out')
            os.system('blastn -db %s -query %s -outfmt 6 -dust no -num_threads %d -perc_identity 80 -word_size 20 -max_target_seqs %d -evalue 0.0000001 -out %s'
                %(nt_viruses_family_file,eachContigFile,num_threads,max_target_seqs_nt_viruses,eachBlastResult))

    timer_blastn_end = time.clock()
    timer_blastn_end1 = time.time()
    blastn_time = timer_blastn_end - timer_blastn_start
    blastn_time1 = timer_blastn_end1 - timer_blastn_start1
    '''
    print '@@@@@----------timer:blastn_time_clock:'
    print blastn_time  #0.287237
    print '@@@@@----------timer:blastn_time_time:'
    print blastn_time1  #30.0745470524
    '''
    sheetPipelinerStaClock.write(1,3,blastn_time)
    sheetPipelinerStaTime.write(1,3,blastn_time1)

    ##3.third:sorted the blast result, use the different colomn
    for i in range(len(colNumSorted_list)):
        sort_blast_result_by_blastnColumn(colNumSorted_list[i])


    time_total_end = time.clock()
    time_total_end1 = time.time()
    

    total_time = time_total_end - time_total_start
    total_time1 = time_total_end1 - time_total_start1
    '''
    print '@@@@@-----------timer:total_time_clock:'
    print total_time  #23.519992
    print '@@@@@-----------timer:total_time_time:'
    print total_time1  #2177.6757071
    '''
    sheetPipelinerStaClock.write(1,4,total_time)
    sheetPipelinerStaTime.write(1,4,total_time1)

    pipelinerSta.save(os.path.join(result_total_sta_path,'pipeline_time_sta.xls'))




##method : mapAginstVirusFirstwithoutQC **
def process_pipeliner_pe_mapAginstVirusFirstwithoutQC(genomeHostName,ntVirusesFamilyName):
    ##rebuild : use all genome host
    genome_host_path = genomesDB_path+genomeHostName
    genome_host_file = os.path.join(genome_host_path,genomeHostName+'.fa')

    genome_host_bt2_index_base_name = 'bt2_index_%s'%(genomeHostName)
    
    ##rebuild : use ntVirusesFamilyName to replace nt_viruses_family_name_list[i], used for mapping with bowtie2 on nt_viruses_ and blastn
    family_path = 'nt_%s'%(ntVirusesFamilyName)
    family_file_name = 'nt_%s.fasta'%(ntVirusesFamilyName)
    
    nt_viruses_family_path = blastDB_path+family_path
    nt_viruses_family_file = os.path.join(nt_viruses_family_path,family_file_name)
        
    nt_virus_bt2_index_base_name = 'bt2_index_%s'%(ntVirusesFamilyName)

    nt_viruses_family_index_header_file_name = 'nt_%s_index_title'%(ntVirusesFamilyName)  #0_index_title / Adenoviridae_index_title
    nt_viruses_family_index_header_file = os.path.join(nt_viruses_family_path,nt_viruses_family_index_header_file_name)

    #print sample_list

    time_total_start = time.clock()
    time_total_start1 = time.time()

    pipelinerSta = xlwt.Workbook()
    sheetPipelinerStaClock = pipelinerSta.add_sheet('pipeline_clock_statistics',cell_overwrite_ok=True)
    sheetPipelinerStaTime = pipelinerSta.add_sheet('pipeline_time_statistics',cell_overwrite_ok=True)
    #rowTitle = ['qc_trimmomatic','map_hp38','map_genebank_virus','map_nt_virus','assembly_velvet','map_blastn','total_time']
    #rowTitle = ['qc_trimmomatic','map_hp38','map_nt_virus','assembly_velvet','map_blastn','total_time']
    rowTitle = ['map_%s'%ntVirusesFamilyName,'map_%s'%genomeHostName,'assembly_velvet','map_blastn','total_time']
    for i in range(0,len(rowTitle)):
        sheetPipelinerStaClock.write(0,i,rowTitle[i]) #first row
        sheetPipelinerStaTime.write(0,i,rowTitle[i])


    print '==========map and align : AGINST NT_viruses reference,keep mapped reads'

    #bowtie2IndexNTVirusFamily = nt_viruses_family_name_list[0]
    #bowtie2IndexNTVirusFamily = bt2_index_base_name_list[0]  #'bt2_index_Adenoviridae'

    timer_bowtie2_nt_virus_start = time.clock()
    timer_bowtie2_nt_virus_start1 = time.time()

    outputFastqPair1FileNTVirus_list = []
    outputFastqPair2FileNTVirus_list = []
    outputFastqUnpairFileNTVirus_list = []
    for i in range(len(sample_list)):
        outputSAMFileNTVirus = bowtie2_path+sample_list[i]+map_nt_viruses_singleAndunpair_sam_file_suffix
        outputFastqPairFileNTVirus = bowtie2_path+sample_list[i]+'_pair_VIRUS.fastq'
        outputFastqPair1FileNTVirus = bowtie2_path+sample_list[i]+'_pair_VIRUS.1.fastq'
        outputFastqPair2FileNTVirus = bowtie2_path+sample_list[i]+'_pair_VIRUS.2.fastq'
        
        outputFastqUnpairFileNTVirus = bowtie2_path+sample_list[i]+map_nt_viruses_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -1 %s -2 %s --no-head --no-unal -S %s --al-conc %s --al %s -p %d'
            %(nt_viruses_family_path,nt_virus_bt2_index_base_name,samplesRawDataForwardFiles[i],samplesRawDataReverseFiles[i],outputSAMFileNTVirus,outputFastqPairFileNTVirus,outputFastqUnpairFileNTVirus,num_threads))

        outputFastqPair1FileNTVirus_list.append(outputFastqPair1FileNTVirus)
        outputFastqPair2FileNTVirus_list.append(outputFastqPair2FileNTVirus)
        outputFastqUnpairFileNTVirus_list.append(outputFastqUnpairFileNTVirus)
        
    timer_bowtie2_nt_virus_end = time.clock()
    timer_bowtie2_nt_virus_end1 = time.time()
    bowtie2_nt_virus_time = timer_bowtie2_nt_virus_end - timer_bowtie2_nt_virus_start
    bowtie2_nt_virus_time1 = timer_bowtie2_nt_virus_end1 - timer_bowtie2_nt_virus_start1
    '''
    print '@@@@@----------timer:bowtie2_nt_virus_time_clock:'
    print bowtie2_nt_virus_time    #18.209088
    print '@@@@@----------timer:bowtie2_nt_virus_time_time:'
    print bowtie2_nt_virus_time1    #1678.66673303
    '''
    sheetPipelinerStaClock.write(1,0,bowtie2_nt_virus_time)
    sheetPipelinerStaTime.write(1,0,bowtie2_nt_virus_time1)


    print '==========map and align : AGINST hg38 reference,keep unmapped reads'

    timer_bowtie2_hg38_start = time.clock()
    timer_bowtie2_hg38_start1 = time.time()

    outputFastqPair1FileGBHG_list = []
    outputFastqPair2FileGBHG_list = []
    outputFastqUnpairFileGBHG_list = []

    for i in range(len(sample_list)):
        outputSAMFileNTGBHG = bowtie2_path+sample_list[i]+'_host.sam'
        outputFastqPairFileGBHG = bowtie2_path+sample_list[i]+'_pair_host.fastq'
        outputFastqPair1FileGBHG = bowtie2_path+sample_list[i]+'_pair_host.1.fastq'
        outputFastqPair2FileGBHG = bowtie2_path+sample_list[i]+'_pair_host.2.fastq'
        outputFastqUnpairFileGBHG = bowtie2_path+sample_list[i]+map_gb_hg_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -1 %s -2 %s -U %s --no-head --no-unal -S %s --un-conc %s --un %s -p %d'
        %(genome_host_path,genome_host_bt2_index_base_name,outputFastqPair1FileNTVirus_list[i],outputFastqPair2FileNTVirus_list[i],outputFastqUnpairFileNTVirus_list[i],outputSAMFileNTGBHG,outputFastqPairFileGBHG,outputFastqUnpairFileGBHG,num_threads))

        outputFastqPair1FileGBHG_list.append(outputFastqPair1FileGBHG)
        outputFastqPair2FileGBHG_list.append(outputFastqPair2FileGBHG)
        outputFastqUnpairFileGBHG_list.append(outputFastqUnpairFileGBHG)

    timer_bowtie2_hg38_end = time.clock()
    timer_bowtie2_hg38_end1 = time.time()
    bowtie2_hg38_time = timer_bowtie2_hg38_end - timer_bowtie2_hg38_start
    bowtie2_hg38_time1 = timer_bowtie2_hg38_end1 - timer_bowtie2_hg38_start1
    '''
    print '@@@@@----------timer:bowtie2_hg38_time_clock:'
    print bowtie2_hg38_time    #17.103814
    print '@@@@@----------timer:bowtie2_hg38_time_time:'
    print bowtie2_hg38_time1    #1568.48754096
    '''
    sheetPipelinerStaClock.write(1,1,bowtie2_hg38_time)
    sheetPipelinerStaTime.write(1,1,bowtie2_hg38_time1)


    print '==========assembly : velvet with alignment data above'

    timer_velvet_start = time.clock()
    timer_velvet_start1 = time.time()

    for i in range(len(sample_list)):
        #inputPE = velvet_result+sample_list[i]+'_pe.fastq'
        #os.system('%s %s %s %s' %(shuffleSequences_fastq_path,outputFastq1AfterHGandVirus[i],outputFastq2AfterHGandVirus[i],inputPE))
        for kmer in range(kmerMin,kmerMax,kmerStep):
            #print sample_list[i]+'************************************'
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_list[i],kmer)
            os.system('velveth %s %d -shortPaired -fastq -separate %s %s -short -fastq %s'%(eachKmerDir,kmer,outputFastqPair1FileGBHG_list[i],outputFastqPair2FileGBHG_list[i],outputFastqUnpairFileGBHG_list[i]))
            os.system('velvetg %s -exp_cov auto -cov_cutoff auto ' %(eachKmerDir))

            eachKmerFiles = os.listdir(eachKmerDir)
            eachKmerFilesAP = []
            for j in range(len(eachKmerFiles)):
                eachKmerFilesAP.append(eachKmerDir+'/'+eachKmerFiles[j])
            for j in range(len(eachKmerFilesAP)):
                #pass
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Graph2' : os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/LastGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/PreGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Roadmaps': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Sequences': os.remove(eachKmerFilesAP[j])

    timer_velvet_end = time.clock()
    timer_velvet_end1 = time.time()
    velvet_time = timer_velvet_end - timer_velvet_start
    velvet_time1 = timer_velvet_end1 - timer_velvet_start1
    '''
    print '@@@@@----------timer:velvet_time_clock:'
    print velvet_time    #0.034675
    print '@@@@@----------timer:velvet_time_time:'
    print velvet_time1    #1.79084706306
    '''
    sheetPipelinerStaClock.write(1,2,velvet_time)
    sheetPipelinerStaTime.write(1,2,velvet_time1)

    ##================statistics velvet : pick up the proper kmer , accoding to N50 max========
    velvetSta = xlwt.Workbook()
    sheetVelvetSta = velvetSta.add_sheet('velvet_statistics', cell_overwrite_ok=True)
    kmerNum = (kmerMax-kmerMin)/2  #(35-31)/2=2
    row0 = ['sampleNums','kmer','medianCovDepth','nodesNum','contigN50','contigMax','genomeLen','usedReads','totalReads']
    for i in range(0,len(row0)):
        sheetVelvetSta.write(0,i,row0[i])  #first row

    #first cow is sample_list
    #second cow is kmer (kmerMin to kmerMax , kmerStep)

    for i in range(0,len(sample_list)):
        sheetVelvetSta.write(kmerNum*i+1,0,sample_list[i])   #first column is the sampleNums

    colKmer = range(kmerMin,kmerMax,kmerStep)
    for i in range(0,len(sample_list)*kmerNum,kmerNum): #i : every column
        for k in range(0,kmerNum):
            sheetVelvetSta.write(i+k+1,1,colKmer[k])   #second column is the k-mer
            k += 1

    logLast2Str = []
    MedianCovDepth = []

    logLast1Str = []
    nodesNum = []
    contigN50 = []
    contigMax = []
    genomeLen = []
    usedReads = []
    totalReads = []
    for i in range(len(sample_list)):
        for kmer in range(kmerMin,kmerMax,kmerStep):
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_list[i],kmer)
            Log_file = os.path.join(eachKmerDir,'Log')
            logs = open(Log_file,'r').readlines()
            #logLast2line = logs[len(logs)-2]
            #if logLast2line != '':
            #Median coverage depth = 96.434783  #not the logs[26] , to avoid the task again with adding the log_information to Log File
            #logLast2Str = logLast2line.split(' ')
            #print logLast2Str
            #print logLast2Str[4][:-2]
            #MedianCovDepth.append(string.atof(logLast2Str[4][:-2]))
            logLast1line = logs[len(logs)-1]  #Final graph has 199 nodes and n50 of 47325, max 166302, total 1570595, using 1470371/1487666 reads
            logLast1Str = logLast1line.split(' ')
            #['Final', 'graph', 'has', '172', 'nodes', 'and', 'n50', 'of', '46917,', 'max', '237233,', 'total', '1570012,', 'using', '1478778/1487666', 'reads\n']
            #print logLast1Str
            nodesNum.append(string.atoi(logLast1Str[3], 10))
            contigN50.append(string.atoi(logLast1Str[8][:-1],10))
            contigMax.append(string.atoi(logLast1Str[10][:-1],10))
            genomeLen.append(string.atoi(logLast1Str[12][:-1],10))
            usedReads.append(string.atoi(logLast1Str[14].split('/')[0], 10))
            totalReads.append(string.atoi(logLast1Str[14].split('/')[1], 10))
    #print contigN50
    for i in range(0,len(sample_list)*kmerNum):
        #sheetVelvetSta.write(i+1,2,MedianCovDepth[i])
        sheetVelvetSta.write(i+1,3,nodesNum[i])
        sheetVelvetSta.write(i+1,4,contigN50[i])
        sheetVelvetSta.write(i+1,5,contigMax[i])
        sheetVelvetSta.write(i+1,6,genomeLen[i])
        sheetVelvetSta.write(i+1,7,usedReads[i])
        sheetVelvetSta.write(i+1,8,totalReads[i])

    velvetSta.save(os.path.join(result_total_sta_path,'velvet_log_sta.xls'))

    contigN50_max_value_list = []  #all_samples
    contigN50_max_index_list = []
    for i in range(len(sample_list)):
        each_sample_contigN50_list = contigN50[kmerNum*i:kmerNum*(i+1)]
        #print each_sample_contigN50_list
        contigN50_max_value_list.append(max(each_sample_contigN50_list))
        contigN50_max_index_list.append(kmerMin+(each_sample_contigN50_list.index(max(each_sample_contigN50_list)))*kmerStep)

    #print contigN50_max_value_list
    #print contigN50_max_index_list

    print '=============blast:use the contigs.fa to map AGINST NT_viruses'
    #1.first:format the local nt fasta_database file
    #os.system('makeblastdb -dbtype nucl -input_type fasta -in %s' %(nt_path))  ###big memory : server

    #after makeblastdb,generate some other indexed files besides nt file,
    #use the indexed nt file(with some other indexed files) as the following input ref db '-db'.

    #2.second:blastn

    timer_blastn_start = time.clock()
    timer_blastn_start1 = time.time()

    #-outfmt 6
    for i in range(len(sample_list)):
        kmer = contigN50_max_index_list[i]
        eachKmerDir = '%s%s_%d'%(velvet_path,sample_list[i],kmer) #/public1/home/yefq/project_adenovirus/results/velvet_result/jiayan_31
        #print eachKmerDir
        eachContigFile = os.path.join(eachKmerDir,'contigs.fa')
        #print eachContigFile
        #estimate the contigs.fa is empty or not.
        if os.path.getsize(eachContigFile) != 0:
            eachBlastResult = os.path.join(blast_path,sample_list[i]+'.out')
            os.system('blastn -db %s -query %s -outfmt 6 -dust no -num_threads %d -perc_identity 80 -word_size 20 -max_target_seqs %d -evalue 0.0000001 -out %s'
                %(nt_viruses_family_file,eachContigFile,num_threads,max_target_seqs_nt_viruses,eachBlastResult))

    timer_blastn_end = time.clock()
    timer_blastn_end1 = time.time()
    blastn_time = timer_blastn_end - timer_blastn_start
    blastn_time1 = timer_blastn_end1 - timer_blastn_start1
    '''
    print '@@@@@----------timer:blastn_time_clock:'
    print blastn_time  #0.287237
    print '@@@@@----------timer:blastn_time_time:'
    print blastn_time1  #30.0745470524
    '''
    sheetPipelinerStaClock.write(1,3,blastn_time)
    sheetPipelinerStaTime.write(1,3,blastn_time1)

    ##3.third:sorted the blast result, use the different colomn
    for i in range(len(colNumSorted_list)):
        sort_blast_result_by_blastnColumn(colNumSorted_list[i])


    time_total_end = time.clock()
    time_total_end1 = time.time()
    

    total_time = time_total_end - time_total_start
    total_time1 = time_total_end1 - time_total_start1
    '''
    print '@@@@@-----------timer:total_time_clock:'
    print total_time  #23.519992
    print '@@@@@-----------timer:total_time_time:'
    print total_time1  #2177.6757071
    '''
    sheetPipelinerStaClock.write(1,4,total_time)
    sheetPipelinerStaTime.write(1,4,total_time1)

    pipelinerSta.save(os.path.join(result_total_sta_path,'pipeline_time_sta.xls'))

##method : mapAginstVirusFirstwithoutQC
def process_pipeliner_se_mapAginstVirusFirstwithoutQC(genomeHostName,ntVirusesFamilyName):
    ##rebuild : use all genome host
    genome_host_path = genomesDB_path+genomeHostName
    genome_host_file = os.path.join(genome_host_path,genomeHostName+'.fa')

    genome_host_bt2_index_base_name = 'bt2_index_%s'%(genomeHostName)
    
    ##rebuild : use ntVirusesFamilyName to replace nt_viruses_family_name_list[i], used for mapping with bowtie2 on nt_viruses_ and blastn
    family_path = 'nt_%s'%(ntVirusesFamilyName)
    family_file_name = 'nt_%s.fasta'%(ntVirusesFamilyName)
    
    nt_viruses_family_path = blastDB_path+family_path
    nt_viruses_family_file = os.path.join(nt_viruses_family_path,family_file_name)
        
    nt_virus_bt2_index_base_name = 'bt2_index_%s'%(ntVirusesFamilyName)

    nt_viruses_family_index_header_file_name = 'nt_%s_index_title'%(ntVirusesFamilyName)  #0_index_title / Adenoviridae_index_title
    nt_viruses_family_index_header_file = os.path.join(nt_viruses_family_path,nt_viruses_family_index_header_file_name)

    
    #print sample_se_list

    time_total_start = time.clock()
    time_total_start1 = time.time()

    pipelinerSta = xlwt.Workbook()
    sheetPipelinerStaClock = pipelinerSta.add_sheet('pipeline_clock_statistics',cell_overwrite_ok=True)
    sheetPipelinerStaTime = pipelinerSta.add_sheet('pipeline_time_statistics',cell_overwrite_ok=True)
    
    rowTitle = ['map_%s'%ntVirusesFamilyName,'map_%s'%genomeHostName,'assembly_velvet','map_blastn','total_time']
    for i in range(0,len(rowTitle)):
        sheetPipelinerStaClock.write(0,i,rowTitle[i]) #first row
        sheetPipelinerStaTime.write(0,i,rowTitle[i])



    print '==========map and align : AGINST hg38 reference,keep unmapped reads'

    timer_bowtie2_hg38_start = time.clock()
    timer_bowtie2_hg38_start1 = time.time()

    afterHGfastqFiles =[]
    for i in range(len(sample_se_list)):
        outputSAMFileGBHG = bowtie2_path+sample_se_list[i]+'_host.sam'
        outputFASTQFileGBHG = bowtie2_path+sample_se_list[i]+map_gb_hg_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -U %s --no-head --no-unal -S %s --un %s -p %d'
            %(genome_host_path,genome_host_bt2_index_base_name,samplesRawDataSingleFiles[i],outputSAMFileGBHG,outputFASTQFileGBHG,num_threads))

        afterHGfastqFiles.append(outputFASTQFileGBHG)



    timer_bowtie2_hg38_end = time.clock()
    timer_bowtie2_hg38_end1 = time.time()
    bowtie2_hg38_time = timer_bowtie2_hg38_end - timer_bowtie2_hg38_start
    bowtie2_hg38_time1 = timer_bowtie2_hg38_end1 - timer_bowtie2_hg38_start1
    '''
    print '@@@@@----------timer:bowtie2_hg38_time_clock:'
    print bowtie2_hg38_time    #17.103814
    print '@@@@@----------timer:bowtie2_hg38_time_time:'
    print bowtie2_hg38_time1    #1568.48754096
    '''
    sheetPipelinerStaClock.write(1,0,bowtie2_hg38_time)
    sheetPipelinerStaTime.write(1,0,bowtie2_hg38_time1)


    print '==========map and align : AGINST NT_viruses reference,keep mapped reads'

    #bowtie2IndexNTVirusFamily = nt_viruses_family_name_list[0]
    #bowtie2IndexNTVirusFamily = bt2_index_base_name_list[0]  #'bt2_index_Adenoviridae'

    timer_bowtie2_nt_virus_start = time.clock()
    timer_bowtie2_nt_virus_start1 = time.time()

    afterNTVIRUSfastqFiles =[]
    for i in range(len(sample_se_list)):
        outputSAMFileNTVirus = bowtie2_path+sample_se_list[i]+map_nt_viruses_singleAndunpair_sam_file_suffix
        outputFASTQFileNTVirus = bowtie2_path+sample_se_list[i]+map_nt_viruses_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -U %s --no-head --no-unal -S %s --al %s -p %d'
            %(nt_viruses_family_path,nt_virus_bt2_index_base_name,afterHGfastqFiles[i],outputSAMFileNTVirus,outputFASTQFileNTVirus,num_threads))

        afterNTVIRUSfastqFiles.append(outputFASTQFileNTVirus)


    timer_bowtie2_nt_virus_end = time.clock()
    timer_bowtie2_nt_virus_end1 = time.time()
    bowtie2_nt_virus_time = timer_bowtie2_nt_virus_end - timer_bowtie2_nt_virus_start
    bowtie2_nt_virus_time1 = timer_bowtie2_nt_virus_end1 - timer_bowtie2_nt_virus_start1
    '''
    print '@@@@@----------timer:bowtie2_nt_virus_time_clock:'
    print bowtie2_nt_virus_time    #18.209088
    print '@@@@@----------timer:bowtie2_nt_virus_time_time:'
    print bowtie2_nt_virus_time1    #1678.66673303
    '''
    sheetPipelinerStaClock.write(1,1,bowtie2_nt_virus_time)
    sheetPipelinerStaTime.write(1,1,bowtie2_nt_virus_time1)


    print '==========assembly : velvet with alignment data above'

    timer_velvet_start = time.clock()
    timer_velvet_start1 = time.time()

    minContigLength = 100
    for i in range(len(sample_se_list)):
        for kmer in range(kmerMin,kmerMax,kmerStep):
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_se_list[i],kmer)
            os.system('velveth %s %d -short -fastq %s' %(eachKmerDir,kmer,afterNTVIRUSfastqFiles[i]))
            os.system('velvetg %s -exp_cov auto -cov_cutoff auto -min_contig_lgth %d' %(eachKmerDir,minContigLength))

            eachKmerFiles = os.listdir(eachKmerDir)
            eachKmerFilesAP = []
            for j in range(len(eachKmerFiles)):
                eachKmerFilesAP.append(eachKmerDir+'/'+eachKmerFiles[j])
            for j in range(len(eachKmerFilesAP)):
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Graph2' : os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/LastGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/PreGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Roadmaps': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Sequences': os.remove(eachKmerFilesAP[j])

    timer_velvet_end = time.clock()
    timer_velvet_end1 = time.time()
    velvet_time = timer_velvet_end - timer_velvet_start
    velvet_time1 = timer_velvet_end1 - timer_velvet_start1
    '''
    print '@@@@@----------timer:velvet_time_clock:'
    print velvet_time    #0.034675
    print '@@@@@----------timer:velvet_time_time:'
    print velvet_time1    #1.79084706306
    '''
    sheetPipelinerStaClock.write(1,2,velvet_time)
    sheetPipelinerStaTime.write(1,2,velvet_time1)

    ##================statistics velvet : pick up the proper kmer , accoding to N50 max========
    velvetSta = xlwt.Workbook()
    sheetVelvetSta = velvetSta.add_sheet('velvet_statistics', cell_overwrite_ok=True)
    kmerNum = (kmerMax-kmerMin)/2  #(35-31)/2=2
    row0 = ['sampleNums','kmer','medianCovDepth','nodesNum','contigN50','contigMax','genomeLen','usedReads','totalReads']
    for i in range(0,len(row0)):
        sheetVelvetSta.write(0,i,row0[i])  #first row

    #first cow is sample_se_list
    #second cow is kmer (kmerMin to kmerMax , kmerStep)

    for i in range(0,len(sample_se_list)):
        sheetVelvetSta.write(kmerNum*i+1,0,sample_se_list[i])   #first column is the sampleNums

    colKmer = range(kmerMin,kmerMax,kmerStep)
    for i in range(0,len(sample_se_list)*kmerNum,kmerNum): #i : every column
        for k in range(0,kmerNum):
            sheetVelvetSta.write(i+k+1,1,colKmer[k])   #second column is the k-mer
            k += 1

    logLast2Str = []
    MedianCovDepth = []

    logLast1Str = []
    nodesNum = []
    contigN50 = []
    contigMax = []
    genomeLen = []
    usedReads = []
    totalReads = []
    for i in range(len(sample_se_list)):
        for kmer in range(kmerMin,kmerMax,kmerStep):
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_se_list[i],kmer)
            Log_file = os.path.join(eachKmerDir,'Log')
            logs = open(Log_file,'r').readlines()
            #logLast2line = logs[len(logs)-2]
            #if logLast2line != '':
            #Median coverage depth = 96.434783  #not the logs[26] , to avoid the task again with adding the log_information to Log File
            #logLast2Str = logLast2line.split(' ')
            #print logLast2Str
            #print logLast2Str[4][:-2]
            #MedianCovDepth.append(string.atof(logLast2Str[4][:-2]))
            logLast1line = logs[len(logs)-1]  #Final graph has 199 nodes and n50 of 47325, max 166302, total 1570595, using 1470371/1487666 reads
            logLast1Str = logLast1line.split(' ')
            #['Final', 'graph', 'has', '172', 'nodes', 'and', 'n50', 'of', '46917,', 'max', '237233,', 'total', '1570012,', 'using', '1478778/1487666', 'reads\n']
            #print logLast1Str
            nodesNum.append(string.atoi(logLast1Str[3], 10))
            contigN50.append(string.atoi(logLast1Str[8][:-1],10))
            contigMax.append(string.atoi(logLast1Str[10][:-1],10))
            genomeLen.append(string.atoi(logLast1Str[12][:-1],10))
            usedReads.append(string.atoi(logLast1Str[14].split('/')[0], 10))
            totalReads.append(string.atoi(logLast1Str[14].split('/')[1], 10))
    #print contigN50
    for i in range(0,len(sample_se_list)*kmerNum):
        #sheetVelvetSta.write(i+1,2,MedianCovDepth[i])
        sheetVelvetSta.write(i+1,3,nodesNum[i])
        sheetVelvetSta.write(i+1,4,contigN50[i])
        sheetVelvetSta.write(i+1,5,contigMax[i])
        sheetVelvetSta.write(i+1,6,genomeLen[i])
        sheetVelvetSta.write(i+1,7,usedReads[i])
        sheetVelvetSta.write(i+1,8,totalReads[i])

    velvetSta.save(os.path.join(result_total_sta_path,'velvet_log_sta.xls'))

    ##check the max contigN50 to pick up ,then blastn
    contigN50_max_value_list = []  #all_samples
    contigN50_max_index_list = []
    for i in range(len(sample_se_list)):
        each_sample_contigN50_list = contigN50[kmerNum*i:kmerNum*(i+1)]
        #print each_sample_contigN50_list
        contigN50_max_value_list.append(max(each_sample_contigN50_list))
        contigN50_max_index_list.append(kmerMin+(each_sample_contigN50_list.index(max(each_sample_contigN50_list)))*kmerStep)

    #print contigN50_max_value_list
    #print contigN50_max_index_list

    print '=============blast:use the contigs.fa to map AGINST NT_viruses'
    #1.first:format the local nt fasta_database file
    #os.system('makeblastdb -dbtype nucl -input_type fasta -in %s' %(nt_path))  ###big memory : server

    #after makeblastdb,generate some other indexed files besides nt file,
    #use the indexed nt file(with some other indexed files) as the following input ref db '-db'.

    #2.second:blastn

    timer_blastn_start = time.clock()
    timer_blastn_start1 = time.time()

                 
    for i in range(len(sample_se_list)):
        kmer = contigN50_max_index_list[i]
        eachKmerDir = '%s%s_%d'%(velvet_path,sample_se_list[i],kmer) #/public1/home/yefq/project_adenovirus/results/velvet_result/jiayan_31
        #print eachKmerDir
        eachContigFile = os.path.join(eachKmerDir,'contigs.fa')
        #print eachContigFile
        if os.path.getsize(eachContigFile) != 0:
            eachBlastResult = os.path.join(blast_path,sample_se_list[i]+'.out')
            os.system('blastn -db %s -query %s -outfmt 6 -dust no -num_threads %d -perc_identity 80 -word_size 20 -max_target_seqs %d -evalue 0.0000001 -out %s'
                %(nt_viruses_family_file,eachContigFile,num_threads,max_target_seqs_nt_viruses,eachBlastResult))

    timer_blastn_end = time.clock()
    timer_blastn_end1 = time.time()
    blastn_time = timer_blastn_end - timer_blastn_start
    blastn_time1 = timer_blastn_end1 - timer_blastn_start1
    '''
    print '@@@@@----------timer:blastn_time_clock:'
    print blastn_time  #0.287237
    print '@@@@@----------timer:blastn_time_time:'
    print blastn_time1  #30.0745470524
    '''
    sheetPipelinerStaClock.write(1,3,blastn_time)
    sheetPipelinerStaTime.write(1,3,blastn_time1)

    ##3.third:sorted the blast result, use the different colomn
    for i in range(len(colNumSorted_list)):
        sort_blast_result_by_blastnColumn(colNumSorted_list[i])


    time_total_end = time.clock()
    time_total_end1 = time.time()
    

    total_time = time_total_end - time_total_start
    total_time1 = time_total_end1 - time_total_start1
    '''
    print '@@@@@-----------timer:total_time_clock:'
    print total_time  #23.519992
    print '@@@@@-----------timer:total_time_time:'
    print total_time1  #2177.6757071
    '''
    sheetPipelinerStaClock.write(1,4,total_time)
    sheetPipelinerStaTime.write(1,4,total_time1)

    pipelinerSta.save(os.path.join(result_total_sta_path,'pipeline_time_sta.xls'))


##method : mapAginstVirusFirst
def process_pipeliner_pe_mapAginstVirusFirst(genomeHostName,ntVirusesFamilyName):
    ##rebuild : use all genome host
    genome_host_path = genomesDB_path+genomeHostName
    genome_host_file = os.path.join(genome_host_path,genomeHostName+'.fa')

    genome_host_bt2_index_base_name = 'bt2_index_%s'%(genomeHostName)
    
    ##rebuild : use ntVirusesFamilyName to replace nt_viruses_family_name_list[i], used for mapping with bowtie2 on nt_viruses_ and blastn
    family_path = 'nt_%s'%(ntVirusesFamilyName)
    family_file_name = 'nt_%s.fasta'%(ntVirusesFamilyName)
    
    nt_viruses_family_path = blastDB_path+family_path
    nt_viruses_family_file = os.path.join(nt_viruses_family_path,family_file_name)
        
    nt_virus_bt2_index_base_name = 'bt2_index_%s'%(ntVirusesFamilyName)

    nt_viruses_family_index_header_file_name = 'nt_%s_index_title'%(ntVirusesFamilyName)  #0_index_title / Adenoviridae_index_title
    nt_viruses_family_index_header_file = os.path.join(nt_viruses_family_path,nt_viruses_family_index_header_file_name)

    #print sample_list

    time_total_start = time.clock()
    time_total_start1 = time.time()

    pipelinerSta = xlwt.Workbook()
    sheetPipelinerStaClock = pipelinerSta.add_sheet('pipeline_clock_statistics',cell_overwrite_ok=True)
    sheetPipelinerStaTime = pipelinerSta.add_sheet('pipeline_time_statistics',cell_overwrite_ok=True)
    #rowTitle = ['qc_trimmomatic','map_hp38','map_genebank_virus','map_nt_virus','assembly_velvet','map_blastn','total_time']
    #rowTitle = ['qc_trimmomatic','map_hp38','map_nt_virus','assembly_velvet','map_blastn','total_time']
    rowTitle = ['qc_trimmomatic','map_%s'%ntVirusesFamilyName,'map_%s'%genomeHostName,'assembly_velvet','map_blastn','total_time']
    for i in range(0,len(rowTitle)):
        sheetPipelinerStaClock.write(0,i,rowTitle[i]) #first row
        sheetPipelinerStaTime.write(0,i,rowTitle[i])


    
    print '==========qc:trimmomatic'
    timer_trimmomatic_start = time.clock()
    timer_trimmomatic_start1 = time.time()

    trimmomatic_results_prefix = []
    for i in range(len(sample_list)):
        trimmomatic_result_prefix = trimmomatic_path+sample_list[i]
        trimmomatic_results_prefix.append(trimmomatic_result_prefix)
    #print trimmomatic_results_prefix ##/public1/home/yefq/project_adenovirus/results/trimmomatic_result/jiayan

    outputFPFiles = []
    outputFUFiles = []
    outputRPFiles = []
    outputRUFiles = []
    for i in range(len(sample_list)):
        #print sample_list[i]
        trimlogFile=trimmomatic_results_prefix[i]+'_trimlog.txt'
        outputFPFile=trimmomatic_results_prefix[i]+'_output_forward_paired.fastq'
        outputFUFile=trimmomatic_results_prefix[i]+'_output_forward_unpaired.fastq'
        outputRPFile=trimmomatic_results_prefix[i]+'_output_reverse_paired.fastq'
        outputRUFile=trimmomatic_results_prefix[i]+'_output_reverse_unpaired.fastq'

        os.system('java -jar /usr/share/java/trimmomatic.jar PE -threads %d -trimlog %s %s %s %s %s %s %s ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'
        %(num_threads,trimlogFile,samplesRawDataForwardFiles[i],samplesRawDataReverseFiles[i],outputFPFile,outputFUFile,outputRPFile,outputRUFile))

        outputFPFiles.append(outputFPFile)
        outputRPFiles.append(outputRPFile)
        outputFUFiles.append(outputFUFile)
        outputRUFiles.append(outputRUFile)

    timer_trimmomatic_end = time.clock()
    timer_trimmomatic_end1 = time.time()
    trimmomatic_time = timer_trimmomatic_end - timer_trimmomatic_start
    trimmomatic_time1 = timer_trimmomatic_end1 - timer_trimmomatic_start1
    '''
    print '@@@@@----------timer:trimmomatic_time_clock:'
    print trimmomatic_time  #4.822289
    print '@@@@@----------timer:trimmomatic_time_time:'
    print trimmomatic_time1  #466.976895094
    '''
    sheetPipelinerStaClock.write(1,0,trimmomatic_time)
    sheetPipelinerStaTime.write(1,0,trimmomatic_time1)


    print '==========map and align : AGINST NT_viruses reference,keep mapped reads'

    #bowtie2IndexNTVirusFamily = nt_viruses_family_name_list[0]
    #bowtie2IndexNTVirusFamily = bt2_index_base_name_list[0]  #'bt2_index_Adenoviridae'

    timer_bowtie2_nt_virus_start = time.clock()
    timer_bowtie2_nt_virus_start1 = time.time()

    outputFastqPair1FileNTVirus_list = []
    outputFastqPair2FileNTVirus_list = []
    outputFastqUnpairFileNTVirus_list = []
    for i in range(len(sample_list)):
        outputSAMFileNTVirus = bowtie2_path+sample_list[i]+map_nt_viruses_singleAndunpair_sam_file_suffix
        outputFastqPairFileNTVirus = bowtie2_path+sample_list[i]+'_pair_VIRUS.fastq'
        outputFastqPair1FileNTVirus = bowtie2_path+sample_list[i]+'_pair_VIRUS.1.fastq'
        outputFastqPair2FileNTVirus = bowtie2_path+sample_list[i]+'_pair_VIRUS.2.fastq'
        
        outputFastqUnpairFileNTVirus = bowtie2_path+sample_list[i]+map_nt_viruses_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -1 %s -2 %s -U %s,%s --no-head --no-unal -S %s --al-conc %s --al %s -p %d'
            %(nt_viruses_family_path,nt_virus_bt2_index_base_name,outputFPFiles[i],outputRPFiles[i],outputFUFiles[i],outputRUFiles[i],outputSAMFileNTVirus,outputFastqPairFileNTVirus,outputFastqUnpairFileNTVirus,num_threads))

        outputFastqPair1FileNTVirus_list.append(outputFastqPair1FileNTVirus)
        outputFastqPair2FileNTVirus_list.append(outputFastqPair2FileNTVirus)
        outputFastqUnpairFileNTVirus_list.append(outputFastqUnpairFileNTVirus)
        
    timer_bowtie2_nt_virus_end = time.clock()
    timer_bowtie2_nt_virus_end1 = time.time()
    bowtie2_nt_virus_time = timer_bowtie2_nt_virus_end - timer_bowtie2_nt_virus_start
    bowtie2_nt_virus_time1 = timer_bowtie2_nt_virus_end1 - timer_bowtie2_nt_virus_start1
    '''
    print '@@@@@----------timer:bowtie2_nt_virus_time_clock:'
    print bowtie2_nt_virus_time    #18.209088
    print '@@@@@----------timer:bowtie2_nt_virus_time_time:'
    print bowtie2_nt_virus_time1    #1678.66673303
    '''
    sheetPipelinerStaClock.write(1,1,bowtie2_nt_virus_time)
    sheetPipelinerStaTime.write(1,1,bowtie2_nt_virus_time1)


    print '==========map and align : AGINST hg38 reference,keep unmapped reads'

    timer_bowtie2_hg38_start = time.clock()
    timer_bowtie2_hg38_start1 = time.time()

    outputFastqPair1FileGBHG_list = []
    outputFastqPair2FileGBHG_list = []
    outputFastqUnpairFileGBHG_list = []

    for i in range(len(sample_list)):
        outputSAMFileNTGBHG = bowtie2_path+sample_list[i]+'_host.sam'
        outputFastqPairFileGBHG = bowtie2_path+sample_list[i]+'_pair_host.fastq'
        outputFastqPair1FileGBHG = bowtie2_path+sample_list[i]+'_pair_host.1.fastq'
        outputFastqPair2FileGBHG = bowtie2_path+sample_list[i]+'_pair_host.2.fastq'
        outputFastqUnpairFileGBHG = bowtie2_path+sample_list[i]+map_gb_hg_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -1 %s -2 %s -U %s --no-head --no-unal -S %s --un-conc %s --un %s -p %d'
        %(genome_host_path,genome_host_bt2_index_base_name,outputFastqPair1FileNTVirus_list[i],outputFastqPair2FileNTVirus_list[i],outputFastqUnpairFileNTVirus_list[i],outputSAMFileNTGBHG,outputFastqPairFileGBHG,outputFastqUnpairFileGBHG,num_threads))

        outputFastqPair1FileGBHG_list.append(outputFastqPair1FileGBHG)
        outputFastqPair2FileGBHG_list.append(outputFastqPair2FileGBHG)
        outputFastqUnpairFileGBHG_list.append(outputFastqUnpairFileGBHG)

    timer_bowtie2_hg38_end = time.clock()
    timer_bowtie2_hg38_end1 = time.time()
    bowtie2_hg38_time = timer_bowtie2_hg38_end - timer_bowtie2_hg38_start
    bowtie2_hg38_time1 = timer_bowtie2_hg38_end1 - timer_bowtie2_hg38_start1
    '''
    print '@@@@@----------timer:bowtie2_hg38_time_clock:'
    print bowtie2_hg38_time    #17.103814
    print '@@@@@----------timer:bowtie2_hg38_time_time:'
    print bowtie2_hg38_time1    #1568.48754096
    '''
    sheetPipelinerStaClock.write(1,2,bowtie2_hg38_time)
    sheetPipelinerStaTime.write(1,2,bowtie2_hg38_time1)


    print '==========assembly : velvet with alignment data above'

    timer_velvet_start = time.clock()
    timer_velvet_start1 = time.time()

    for i in range(len(sample_list)):
        #inputPE = velvet_result+sample_list[i]+'_pe.fastq'
        #os.system('%s %s %s %s' %(shuffleSequences_fastq_path,outputFastq1AfterHGandVirus[i],outputFastq2AfterHGandVirus[i],inputPE))
        for kmer in range(kmerMin,kmerMax,kmerStep):
            #print sample_list[i]+'************************************'
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_list[i],kmer)
            os.system('velveth %s %d -shortPaired -fastq -separate %s %s -short -fastq %s'%(eachKmerDir,kmer,outputFastqPair1FileGBHG_list[i],outputFastqPair2FileGBHG_list[i],outputFastqUnpairFileGBHG_list[i]))
            os.system('velvetg %s -exp_cov auto -cov_cutoff auto ' %(eachKmerDir))

            eachKmerFiles = os.listdir(eachKmerDir)
            eachKmerFilesAP = []
            for j in range(len(eachKmerFiles)):
                eachKmerFilesAP.append(eachKmerDir+'/'+eachKmerFiles[j])
            for j in range(len(eachKmerFilesAP)):
                #pass
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Graph2' : os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/LastGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/PreGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Roadmaps': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Sequences': os.remove(eachKmerFilesAP[j])

    timer_velvet_end = time.clock()
    timer_velvet_end1 = time.time()
    velvet_time = timer_velvet_end - timer_velvet_start
    velvet_time1 = timer_velvet_end1 - timer_velvet_start1
    '''
    print '@@@@@----------timer:velvet_time_clock:'
    print velvet_time    #0.034675
    print '@@@@@----------timer:velvet_time_time:'
    print velvet_time1    #1.79084706306
    '''
    sheetPipelinerStaClock.write(1,3,velvet_time)
    sheetPipelinerStaTime.write(1,3,velvet_time1)

    ##================statistics velvet : pick up the proper kmer , accoding to N50 max========
    velvetSta = xlwt.Workbook()
    sheetVelvetSta = velvetSta.add_sheet('velvet_statistics', cell_overwrite_ok=True)
    kmerNum = (kmerMax-kmerMin)/2  #(35-31)/2=2
    row0 = ['sampleNums','kmer','medianCovDepth','nodesNum','contigN50','contigMax','genomeLen','usedReads','totalReads']
    for i in range(0,len(row0)):
        sheetVelvetSta.write(0,i,row0[i])  #first row

    #first cow is sample_list
    #second cow is kmer (kmerMin to kmerMax , kmerStep)

    for i in range(0,len(sample_list)):
        sheetVelvetSta.write(kmerNum*i+1,0,sample_list[i])   #first column is the sampleNums

    colKmer = range(kmerMin,kmerMax,kmerStep)
    for i in range(0,len(sample_list)*kmerNum,kmerNum): #i : every column
        for k in range(0,kmerNum):
            sheetVelvetSta.write(i+k+1,1,colKmer[k])   #second column is the k-mer
            k += 1

    logLast2Str = []
    MedianCovDepth = []

    logLast1Str = []
    nodesNum = []
    contigN50 = []
    contigMax = []
    genomeLen = []
    usedReads = []
    totalReads = []
    for i in range(len(sample_list)):
        for kmer in range(kmerMin,kmerMax,kmerStep):
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_list[i],kmer)
            Log_file = os.path.join(eachKmerDir,'Log')
            logs = open(Log_file,'r').readlines()
            #logLast2line = logs[len(logs)-2]
            #if logLast2line != '':
            #Median coverage depth = 96.434783  #not the logs[26] , to avoid the task again with adding the log_information to Log File
            #logLast2Str = logLast2line.split(' ')
            #print logLast2Str
            #print logLast2Str[4][:-2]
            #MedianCovDepth.append(string.atof(logLast2Str[4][:-2]))
            logLast1line = logs[len(logs)-1]  #Final graph has 199 nodes and n50 of 47325, max 166302, total 1570595, using 1470371/1487666 reads
            logLast1Str = logLast1line.split(' ')
            #['Final', 'graph', 'has', '172', 'nodes', 'and', 'n50', 'of', '46917,', 'max', '237233,', 'total', '1570012,', 'using', '1478778/1487666', 'reads\n']
            #print logLast1Str
            nodesNum.append(string.atoi(logLast1Str[3], 10))
            contigN50.append(string.atoi(logLast1Str[8][:-1],10))
            contigMax.append(string.atoi(logLast1Str[10][:-1],10))
            genomeLen.append(string.atoi(logLast1Str[12][:-1],10))
            usedReads.append(string.atoi(logLast1Str[14].split('/')[0], 10))
            totalReads.append(string.atoi(logLast1Str[14].split('/')[1], 10))
    #print contigN50
    for i in range(0,len(sample_list)*kmerNum):
        #sheetVelvetSta.write(i+1,2,MedianCovDepth[i])
        sheetVelvetSta.write(i+1,3,nodesNum[i])
        sheetVelvetSta.write(i+1,4,contigN50[i])
        sheetVelvetSta.write(i+1,5,contigMax[i])
        sheetVelvetSta.write(i+1,6,genomeLen[i])
        sheetVelvetSta.write(i+1,7,usedReads[i])
        sheetVelvetSta.write(i+1,8,totalReads[i])

    velvetSta.save(os.path.join(result_total_sta_path,'velvet_log_sta.xls'))

    contigN50_max_value_list = []  #all_samples
    contigN50_max_index_list = []
    for i in range(len(sample_list)):
        each_sample_contigN50_list = contigN50[kmerNum*i:kmerNum*(i+1)]
        #print each_sample_contigN50_list
        contigN50_max_value_list.append(max(each_sample_contigN50_list))
        contigN50_max_index_list.append(kmerMin+(each_sample_contigN50_list.index(max(each_sample_contigN50_list)))*kmerStep)

    #print contigN50_max_value_list
    #print contigN50_max_index_list

    print '=============blast:use the contigs.fa to map AGINST NT_viruses'
    #1.first:format the local nt fasta_database file
    #os.system('makeblastdb -dbtype nucl -input_type fasta -in %s' %(nt_path))  ###big memory : server

    #after makeblastdb,generate some other indexed files besides nt file,
    #use the indexed nt file(with some other indexed files) as the following input ref db '-db'.

    #2.second:blastn

    timer_blastn_start = time.clock()
    timer_blastn_start1 = time.time()

    #-outfmt 6
    for i in range(len(sample_list)):
        kmer = contigN50_max_index_list[i]
        eachKmerDir = '%s%s_%d'%(velvet_path,sample_list[i],kmer) #/public1/home/yefq/project_adenovirus/results/velvet_result/jiayan_31
        #print eachKmerDir
        eachContigFile = os.path.join(eachKmerDir,'contigs.fa')
        #print eachContigFile
        #estimate the contigs.fa is empty or not.
        if os.path.getsize(eachContigFile) != 0:
            eachBlastResult = os.path.join(blast_path,sample_list[i]+'.out')
            os.system('blastn -db %s -query %s -outfmt 6 -dust no -num_threads %d -perc_identity 80 -word_size 20 -max_target_seqs %d -evalue 0.0000001 -out %s'
                %(nt_viruses_family_file,eachContigFile,num_threads,max_target_seqs_nt_viruses,eachBlastResult))

    timer_blastn_end = time.clock()
    timer_blastn_end1 = time.time()
    blastn_time = timer_blastn_end - timer_blastn_start
    blastn_time1 = timer_blastn_end1 - timer_blastn_start1
    '''
    print '@@@@@----------timer:blastn_time_clock:'
    print blastn_time  #0.287237
    print '@@@@@----------timer:blastn_time_time:'
    print blastn_time1  #30.0745470524
    '''
    sheetPipelinerStaClock.write(1,4,blastn_time)
    sheetPipelinerStaTime.write(1,4,blastn_time1)

    ##3.third:sorted the blast result, use the different colomn
    for i in range(len(colNumSorted_list)):
        sort_blast_result_by_blastnColumn(colNumSorted_list[i])


    time_total_end = time.clock()
    time_total_end1 = time.time()
    

    total_time = time_total_end - time_total_start
    total_time1 = time_total_end1 - time_total_start1
    '''
    print '@@@@@-----------timer:total_time_clock:'
    print total_time  #23.519992
    print '@@@@@-----------timer:total_time_time:'
    print total_time1  #2177.6757071
    '''
    sheetPipelinerStaClock.write(1,5,total_time)
    sheetPipelinerStaTime.write(1,5,total_time1)

    pipelinerSta.save(os.path.join(result_total_sta_path,'pipeline_time_sta.xls'))

##method : mapAginstVirusFirst
def process_pipeliner_se_mapAginstVirusFirst(genomeHostName,ntVirusesFamilyName):
    ##rebuild : use all genome host
    genome_host_path = genomesDB_path+genomeHostName
    genome_host_file = os.path.join(genome_host_path,genomeHostName+'.fa')

    genome_host_bt2_index_base_name = 'bt2_index_%s'%(genomeHostName)
    
    ##rebuild : use ntVirusesFamilyName to replace nt_viruses_family_name_list[i], used for mapping with bowtie2 on nt_viruses_ and blastn
    family_path = 'nt_%s'%(ntVirusesFamilyName)
    family_file_name = 'nt_%s.fasta'%(ntVirusesFamilyName)
    
    nt_viruses_family_path = blastDB_path+family_path
    nt_viruses_family_file = os.path.join(nt_viruses_family_path,family_file_name)
        
    nt_virus_bt2_index_base_name = 'bt2_index_%s'%(ntVirusesFamilyName)

    nt_viruses_family_index_header_file_name = 'nt_%s_index_title'%(ntVirusesFamilyName)  #0_index_title / Adenoviridae_index_title
    nt_viruses_family_index_header_file = os.path.join(nt_viruses_family_path,nt_viruses_family_index_header_file_name)

    
    #print sample_se_list

    time_total_start = time.clock()
    time_total_start1 = time.time()

    pipelinerSta = xlwt.Workbook()
    sheetPipelinerStaClock = pipelinerSta.add_sheet('pipeline_clock_statistics',cell_overwrite_ok=True)
    sheetPipelinerStaTime = pipelinerSta.add_sheet('pipeline_time_statistics',cell_overwrite_ok=True)
    
    rowTitle = ['qc_trimmomatic','map_%s'%ntVirusesFamilyName,'map_%s'%genomeHostName,'assembly_velvet','map_blastn','total_time']
    for i in range(0,len(rowTitle)):
        sheetPipelinerStaClock.write(0,i,rowTitle[i]) #first row
        sheetPipelinerStaTime.write(0,i,rowTitle[i])

    
    print '==========qc:trimmomatic'
    timer_trimmomatic_start = time.clock()
    timer_trimmomatic_start1 = time.time()
    
    MINLEN = 10
    LEADING = 0
    TRAILING = 0    
    outputSingleFiles = []
    for i in range(len(sample_se_list)):
        trimlogFile=trimmomatic_path+sample_se_list[i]+'_trimlog.txt'
        outputSingleFile=trimmomatic_path+sample_se_list[i]+'_output.fastq'

        os.system('java -jar /usr/share/java/trimmomatic.jar SE -threads %d -trimlog %s %s %s ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-SE.fa:2:30:10 LEADING:%d TRAILING:%d SLIDINGWINDOW:4:15 MINLEN:%d'
        %(num_threads,trimlogFile,samplesRawDataSingleFiles[i],outputSingleFile,LEADING,TRAILING,MINLEN))

        outputSingleFiles.append(outputSingleFile)

    timer_trimmomatic_end = time.clock()
    timer_trimmomatic_end1 = time.time()
    trimmomatic_time = timer_trimmomatic_end - timer_trimmomatic_start
    trimmomatic_time1 = timer_trimmomatic_end1 - timer_trimmomatic_start1
    '''
    print '@@@@@----------timer:trimmomatic_time_clock:'
    print trimmomatic_time  #4.822289
    print '@@@@@----------timer:trimmomatic_time_time:'
    print trimmomatic_time1  #466.976895094
    '''
    sheetPipelinerStaClock.write(1,0,trimmomatic_time)
    sheetPipelinerStaTime.write(1,0,trimmomatic_time1)


    print '==========map and align : AGINST hg38 reference,keep unmapped reads'

    timer_bowtie2_hg38_start = time.clock()
    timer_bowtie2_hg38_start1 = time.time()

    afterHGfastqFiles =[]
    for i in range(len(sample_se_list)):
        outputSAMFileGBHG = bowtie2_path+sample_se_list[i]+'_host.sam'
        outputFASTQFileGBHG = bowtie2_path+sample_se_list[i]+map_gb_hg_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -U %s --no-head --no-unal -S %s --un %s -p %d'
            %(genome_host_path,genome_host_bt2_index_base_name,outputSingleFiles[i],outputSAMFileGBHG,outputFASTQFileGBHG,num_threads))

        afterHGfastqFiles.append(outputFASTQFileGBHG)



    timer_bowtie2_hg38_end = time.clock()
    timer_bowtie2_hg38_end1 = time.time()
    bowtie2_hg38_time = timer_bowtie2_hg38_end - timer_bowtie2_hg38_start
    bowtie2_hg38_time1 = timer_bowtie2_hg38_end1 - timer_bowtie2_hg38_start1
    '''
    print '@@@@@----------timer:bowtie2_hg38_time_clock:'
    print bowtie2_hg38_time    #17.103814
    print '@@@@@----------timer:bowtie2_hg38_time_time:'
    print bowtie2_hg38_time1    #1568.48754096
    '''
    sheetPipelinerStaClock.write(1,1,bowtie2_hg38_time)
    sheetPipelinerStaTime.write(1,1,bowtie2_hg38_time1)


    print '==========map and align : AGINST NT_viruses reference,keep mapped reads'

    #bowtie2IndexNTVirusFamily = nt_viruses_family_name_list[0]
    #bowtie2IndexNTVirusFamily = bt2_index_base_name_list[0]  #'bt2_index_Adenoviridae'

    timer_bowtie2_nt_virus_start = time.clock()
    timer_bowtie2_nt_virus_start1 = time.time()

    afterNTVIRUSfastqFiles =[]
    for i in range(len(sample_se_list)):
        outputSAMFileNTVirus = bowtie2_path+sample_se_list[i]+map_nt_viruses_singleAndunpair_sam_file_suffix
        outputFASTQFileNTVirus = bowtie2_path+sample_se_list[i]+map_nt_viruses_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -U %s --no-head --no-unal -S %s --al %s -p %d'
            %(nt_viruses_family_path,nt_virus_bt2_index_base_name,afterHGfastqFiles[i],outputSAMFileNTVirus,outputFASTQFileNTVirus,num_threads))

        afterNTVIRUSfastqFiles.append(outputFASTQFileNTVirus)


    timer_bowtie2_nt_virus_end = time.clock()
    timer_bowtie2_nt_virus_end1 = time.time()
    bowtie2_nt_virus_time = timer_bowtie2_nt_virus_end - timer_bowtie2_nt_virus_start
    bowtie2_nt_virus_time1 = timer_bowtie2_nt_virus_end1 - timer_bowtie2_nt_virus_start1
    '''
    print '@@@@@----------timer:bowtie2_nt_virus_time_clock:'
    print bowtie2_nt_virus_time    #18.209088
    print '@@@@@----------timer:bowtie2_nt_virus_time_time:'
    print bowtie2_nt_virus_time1    #1678.66673303
    '''
    sheetPipelinerStaClock.write(1,2,bowtie2_nt_virus_time)
    sheetPipelinerStaTime.write(1,2,bowtie2_nt_virus_time1)


    print '==========assembly : velvet with alignment data above'

    timer_velvet_start = time.clock()
    timer_velvet_start1 = time.time()

    minContigLength = 100
    for i in range(len(sample_se_list)):
        for kmer in range(kmerMin,kmerMax,kmerStep):
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_se_list[i],kmer)
            os.system('velveth %s %d -short -fastq %s' %(eachKmerDir,kmer,afterNTVIRUSfastqFiles[i]))
            os.system('velvetg %s -exp_cov auto -cov_cutoff auto -min_contig_lgth %d' %(eachKmerDir,minContigLength))

            eachKmerFiles = os.listdir(eachKmerDir)
            eachKmerFilesAP = []
            for j in range(len(eachKmerFiles)):
                eachKmerFilesAP.append(eachKmerDir+'/'+eachKmerFiles[j])
            for j in range(len(eachKmerFilesAP)):
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Graph2' : os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/LastGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/PreGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Roadmaps': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Sequences': os.remove(eachKmerFilesAP[j])

    timer_velvet_end = time.clock()
    timer_velvet_end1 = time.time()
    velvet_time = timer_velvet_end - timer_velvet_start
    velvet_time1 = timer_velvet_end1 - timer_velvet_start1
    '''
    print '@@@@@----------timer:velvet_time_clock:'
    print velvet_time    #0.034675
    print '@@@@@----------timer:velvet_time_time:'
    print velvet_time1    #1.79084706306
    '''
    sheetPipelinerStaClock.write(1,3,velvet_time)
    sheetPipelinerStaTime.write(1,3,velvet_time1)

    ##================statistics velvet : pick up the proper kmer , accoding to N50 max========
    velvetSta = xlwt.Workbook()
    sheetVelvetSta = velvetSta.add_sheet('velvet_statistics', cell_overwrite_ok=True)
    kmerNum = (kmerMax-kmerMin)/2  #(35-31)/2=2
    row0 = ['sampleNums','kmer','medianCovDepth','nodesNum','contigN50','contigMax','genomeLen','usedReads','totalReads']
    for i in range(0,len(row0)):
        sheetVelvetSta.write(0,i,row0[i])  #first row

    #first cow is sample_se_list
    #second cow is kmer (kmerMin to kmerMax , kmerStep)

    for i in range(0,len(sample_se_list)):
        sheetVelvetSta.write(kmerNum*i+1,0,sample_se_list[i])   #first column is the sampleNums

    colKmer = range(kmerMin,kmerMax,kmerStep)
    for i in range(0,len(sample_se_list)*kmerNum,kmerNum): #i : every column
        for k in range(0,kmerNum):
            sheetVelvetSta.write(i+k+1,1,colKmer[k])   #second column is the k-mer
            k += 1

    logLast2Str = []
    MedianCovDepth = []

    logLast1Str = []
    nodesNum = []
    contigN50 = []
    contigMax = []
    genomeLen = []
    usedReads = []
    totalReads = []
    for i in range(len(sample_se_list)):
        for kmer in range(kmerMin,kmerMax,kmerStep):
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_se_list[i],kmer)
            Log_file = os.path.join(eachKmerDir,'Log')
            logs = open(Log_file,'r').readlines()
            #logLast2line = logs[len(logs)-2]
            #if logLast2line != '':
            #Median coverage depth = 96.434783  #not the logs[26] , to avoid the task again with adding the log_information to Log File
            #logLast2Str = logLast2line.split(' ')
            #print logLast2Str
            #print logLast2Str[4][:-2]
            #MedianCovDepth.append(string.atof(logLast2Str[4][:-2]))
            logLast1line = logs[len(logs)-1]  #Final graph has 199 nodes and n50 of 47325, max 166302, total 1570595, using 1470371/1487666 reads
            logLast1Str = logLast1line.split(' ')
            #['Final', 'graph', 'has', '172', 'nodes', 'and', 'n50', 'of', '46917,', 'max', '237233,', 'total', '1570012,', 'using', '1478778/1487666', 'reads\n']
            #print logLast1Str
            nodesNum.append(string.atoi(logLast1Str[3], 10))
            contigN50.append(string.atoi(logLast1Str[8][:-1],10))
            contigMax.append(string.atoi(logLast1Str[10][:-1],10))
            genomeLen.append(string.atoi(logLast1Str[12][:-1],10))
            usedReads.append(string.atoi(logLast1Str[14].split('/')[0], 10))
            totalReads.append(string.atoi(logLast1Str[14].split('/')[1], 10))
    #print contigN50
    for i in range(0,len(sample_se_list)*kmerNum):
        #sheetVelvetSta.write(i+1,2,MedianCovDepth[i])
        sheetVelvetSta.write(i+1,3,nodesNum[i])
        sheetVelvetSta.write(i+1,4,contigN50[i])
        sheetVelvetSta.write(i+1,5,contigMax[i])
        sheetVelvetSta.write(i+1,6,genomeLen[i])
        sheetVelvetSta.write(i+1,7,usedReads[i])
        sheetVelvetSta.write(i+1,8,totalReads[i])

    velvetSta.save(os.path.join(result_total_sta_path,'velvet_log_sta.xls'))

    ##check the max contigN50 to pick up ,then blastn
    contigN50_max_value_list = []  #all_samples
    contigN50_max_index_list = []
    for i in range(len(sample_se_list)):
        each_sample_contigN50_list = contigN50[kmerNum*i:kmerNum*(i+1)]
        #print each_sample_contigN50_list
        contigN50_max_value_list.append(max(each_sample_contigN50_list))
        contigN50_max_index_list.append(kmerMin+(each_sample_contigN50_list.index(max(each_sample_contigN50_list)))*kmerStep)

    #print contigN50_max_value_list
    #print contigN50_max_index_list

    print '=============blast:use the contigs.fa to map AGINST NT_viruses'
    #1.first:format the local nt fasta_database file
    #os.system('makeblastdb -dbtype nucl -input_type fasta -in %s' %(nt_path))  ###big memory : server

    #after makeblastdb,generate some other indexed files besides nt file,
    #use the indexed nt file(with some other indexed files) as the following input ref db '-db'.

    #2.second:blastn

    timer_blastn_start = time.clock()
    timer_blastn_start1 = time.time()

                 
    for i in range(len(sample_se_list)):
        kmer = contigN50_max_index_list[i]
        eachKmerDir = '%s%s_%d'%(velvet_path,sample_se_list[i],kmer) #/public1/home/yefq/project_adenovirus/results/velvet_result/jiayan_31
        #print eachKmerDir
        eachContigFile = os.path.join(eachKmerDir,'contigs.fa')
        #print eachContigFile
        if os.path.getsize(eachContigFile) != 0:
            eachBlastResult = os.path.join(blast_path,sample_se_list[i]+'.out')
            os.system('blastn -db %s -query %s -outfmt 6 -dust no -num_threads %d -perc_identity 80 -word_size 20 -max_target_seqs %d -evalue 0.0000001 -out %s'
                %(nt_viruses_family_file,eachContigFile,num_threads,max_target_seqs_nt_viruses,eachBlastResult))

    timer_blastn_end = time.clock()
    timer_blastn_end1 = time.time()
    blastn_time = timer_blastn_end - timer_blastn_start
    blastn_time1 = timer_blastn_end1 - timer_blastn_start1
    '''
    print '@@@@@----------timer:blastn_time_clock:'
    print blastn_time  #0.287237
    print '@@@@@----------timer:blastn_time_time:'
    print blastn_time1  #30.0745470524
    '''
    sheetPipelinerStaClock.write(1,4,blastn_time)
    sheetPipelinerStaTime.write(1,4,blastn_time1)

    ##3.third:sorted the blast result, use the different colomn
    for i in range(len(colNumSorted_list)):
        sort_blast_result_by_blastnColumn(colNumSorted_list[i])


    time_total_end = time.clock()
    time_total_end1 = time.time()
    

    total_time = time_total_end - time_total_start
    total_time1 = time_total_end1 - time_total_start1
    '''
    print '@@@@@-----------timer:total_time_clock:'
    print total_time  #23.519992
    print '@@@@@-----------timer:total_time_time:'
    print total_time1  #2177.6757071
    '''
    sheetPipelinerStaClock.write(1,5,total_time)
    sheetPipelinerStaTime.write(1,5,total_time1)

    pipelinerSta.save(os.path.join(result_total_sta_path,'pipeline_time_sta.xls'))



##method : mapAginstHumanFirstwithoutQC  **
def process_pipeliner_pe_mapAginstHumanFirstwithoutQC(genomeHostName,ntVirusesFamilyName):
    ##rebuild : use all genome host
    genome_host_path = genomesDB_path+genomeHostName
    genome_host_file = os.path.join(genome_host_path,genomeHostName+'.fa')

    genome_host_bt2_index_base_name = 'bt2_index_%s'%(genomeHostName)
    
    ##rebuild : use ntVirusesFamilyName to replace nt_viruses_family_name_list[i], used for mapping with bowtie2 on nt_viruses_ and blastn
    family_path = 'nt_%s'%(ntVirusesFamilyName)
    family_file_name = 'nt_%s.fasta'%(ntVirusesFamilyName)
    
    nt_viruses_family_path = blastDB_path+family_path
    nt_viruses_family_file = os.path.join(nt_viruses_family_path,family_file_name)
        
    nt_virus_bt2_index_base_name = 'bt2_index_%s'%(ntVirusesFamilyName)

    nt_viruses_family_index_header_file_name = 'nt_%s_index_title'%(ntVirusesFamilyName)  #0_index_title / Adenoviridae_index_title
    nt_viruses_family_index_header_file = os.path.join(nt_viruses_family_path,nt_viruses_family_index_header_file_name)

    #print sample_list

    time_total_start = time.clock()
    time_total_start1 = time.time()

    pipelinerSta = xlwt.Workbook()
    sheetPipelinerStaClock = pipelinerSta.add_sheet('pipeline_clock_statistics',cell_overwrite_ok=True)
    sheetPipelinerStaTime = pipelinerSta.add_sheet('pipeline_time_statistics',cell_overwrite_ok=True)
    #rowTitle = ['qc_trimmomatic','map_hp38','map_genebank_virus','map_nt_virus','assembly_velvet','map_blastn','total_time']
    #rowTitle = ['qc_trimmomatic','map_hp38','map_nt_virus','assembly_velvet','map_blastn','total_time']
    rowTitle = ['map_%s'%genomeHostName,'map_%s'%ntVirusesFamilyName,'assembly_velvet','map_blastn','total_time']
    for i in range(0,len(rowTitle)):
        sheetPipelinerStaClock.write(0,i,rowTitle[i]) #first row
        sheetPipelinerStaTime.write(0,i,rowTitle[i])



    print '==========map and align : AGINST hg38 reference,keep unmapped reads'

    timer_bowtie2_hg38_start = time.clock()
    timer_bowtie2_hg38_start1 = time.time()

    outputFastqPair1FileGBHG_list = []
    outputFastqPair2FileGBHG_list = []
    outputFastqUnpairFileGBHG_list = []

    for i in range(len(sample_list)):
        outputSAMFileNTGBHG = bowtie2_path+sample_list[i]+'_host.sam'
        outputFastqPairFileGBHG = bowtie2_path+sample_list[i]+'_pair_host.fastq'
        outputFastqPair1FileGBHG = bowtie2_path+sample_list[i]+'_pair_host.1.fastq'
        outputFastqPair2FileGBHG = bowtie2_path+sample_list[i]+'_pair_host.2.fastq'
        outputFastqUnpairFileGBHG = bowtie2_path+sample_list[i]+map_gb_hg_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -1 %s -2 %s --no-head --no-unal -S %s --un-conc %s --un %s -p %d'
        %(genome_host_path,genome_host_bt2_index_base_name,samplesRawDataForwardFiles[i],samplesRawDataReverseFiles[i],outputSAMFileNTGBHG,outputFastqPairFileGBHG,outputFastqUnpairFileGBHG,num_threads))

        outputFastqPair1FileGBHG_list.append(outputFastqPair1FileGBHG)
        outputFastqPair2FileGBHG_list.append(outputFastqPair2FileGBHG)
        outputFastqUnpairFileGBHG_list.append(outputFastqUnpairFileGBHG)

    timer_bowtie2_hg38_end = time.clock()
    timer_bowtie2_hg38_end1 = time.time()
    bowtie2_hg38_time = timer_bowtie2_hg38_end - timer_bowtie2_hg38_start
    bowtie2_hg38_time1 = timer_bowtie2_hg38_end1 - timer_bowtie2_hg38_start1
    '''
    print '@@@@@----------timer:bowtie2_hg38_time_clock:'
    print bowtie2_hg38_time    #17.103814
    print '@@@@@----------timer:bowtie2_hg38_time_time:'
    print bowtie2_hg38_time1    #1568.48754096
    '''
    sheetPipelinerStaClock.write(1,0,bowtie2_hg38_time)
    sheetPipelinerStaTime.write(1,0,bowtie2_hg38_time1)


    print '==========map and align : AGINST NT_viruses reference,keep mapped reads'

    #bowtie2IndexNTVirusFamily = nt_viruses_family_name_list[0]
    #bowtie2IndexNTVirusFamily = bt2_index_base_name_list[0]  #'bt2_index_Adenoviridae'

    timer_bowtie2_nt_virus_start = time.clock()
    timer_bowtie2_nt_virus_start1 = time.time()

    outputFastqPair1FileNTVirus_list = []
    outputFastqPair2FileNTVirus_list = []
    outputFastqUnpairFileNTVirus_list = []
    for i in range(len(sample_list)):
        outputSAMFileNTVirus = bowtie2_path+sample_list[i]+map_nt_viruses_singleAndunpair_sam_file_suffix
        outputFastqPairFileNTVirus = bowtie2_path+sample_list[i]+'_pair_VIRUS.fastq'
        outputFastqPair1FileNTVirus = bowtie2_path+sample_list[i]+'_pair_VIRUS.1.fastq'
        outputFastqPair2FileNTVirus = bowtie2_path+sample_list[i]+'_pair_VIRUS.2.fastq'
        
        outputFastqUnpairFileNTVirus = bowtie2_path+sample_list[i]+map_nt_viruses_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -1 %s -2 %s -U %s --no-head --no-unal -S %s --al-conc %s --al %s -p %d'
            %(nt_viruses_family_path,nt_virus_bt2_index_base_name,outputFastqPair1FileGBHG_list[i],outputFastqPair2FileGBHG_list[i],outputFastqUnpairFileGBHG_list[i],outputSAMFileNTVirus,outputFastqPairFileNTVirus,outputFastqUnpairFileNTVirus,num_threads))

        outputFastqPair1FileNTVirus_list.append(outputFastqPair1FileNTVirus)
        outputFastqPair2FileNTVirus_list.append(outputFastqPair2FileNTVirus)
        outputFastqUnpairFileNTVirus_list.append(outputFastqUnpairFileNTVirus)
        
    timer_bowtie2_nt_virus_end = time.clock()
    timer_bowtie2_nt_virus_end1 = time.time()
    bowtie2_nt_virus_time = timer_bowtie2_nt_virus_end - timer_bowtie2_nt_virus_start
    bowtie2_nt_virus_time1 = timer_bowtie2_nt_virus_end1 - timer_bowtie2_nt_virus_start1
    '''
    print '@@@@@----------timer:bowtie2_nt_virus_time_clock:'
    print bowtie2_nt_virus_time    #18.209088
    print '@@@@@----------timer:bowtie2_nt_virus_time_time:'
    print bowtie2_nt_virus_time1    #1678.66673303
    '''
    sheetPipelinerStaClock.write(1,1,bowtie2_nt_virus_time)
    sheetPipelinerStaTime.write(1,1,bowtie2_nt_virus_time1)


    print '==========assembly : velvet with alignment data above'

    timer_velvet_start = time.clock()
    timer_velvet_start1 = time.time()

    for i in range(len(sample_list)):
        #inputPE = velvet_result+sample_list[i]+'_pe.fastq'
        #os.system('%s %s %s %s' %(shuffleSequences_fastq_path,outputFastq1AfterHGandVirus[i],outputFastq2AfterHGandVirus[i],inputPE))
        for kmer in range(kmerMin,kmerMax,kmerStep):
            #print sample_list[i]+'************************************'
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_list[i],kmer)
            os.system('velveth %s %d -shortPaired -fastq -separate %s %s -short -fastq %s'%(eachKmerDir,kmer,outputFastqPair1FileNTVirus_list[i],outputFastqPair2FileNTVirus_list[i],outputFastqUnpairFileNTVirus_list[i]))
            os.system('velvetg %s -exp_cov auto -cov_cutoff auto ' %(eachKmerDir))

            eachKmerFiles = os.listdir(eachKmerDir)
            eachKmerFilesAP = []
            for j in range(len(eachKmerFiles)):
                eachKmerFilesAP.append(eachKmerDir+'/'+eachKmerFiles[j])
            for j in range(len(eachKmerFilesAP)):
                #pass
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Graph2' : os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/LastGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/PreGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Roadmaps': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Sequences': os.remove(eachKmerFilesAP[j])

    timer_velvet_end = time.clock()
    timer_velvet_end1 = time.time()
    velvet_time = timer_velvet_end - timer_velvet_start
    velvet_time1 = timer_velvet_end1 - timer_velvet_start1
    '''
    print '@@@@@----------timer:velvet_time_clock:'
    print velvet_time    #0.034675
    print '@@@@@----------timer:velvet_time_time:'
    print velvet_time1    #1.79084706306
    '''
    sheetPipelinerStaClock.write(1,2,velvet_time)
    sheetPipelinerStaTime.write(1,2,velvet_time1)

    ##================statistics velvet : pick up the proper kmer , accoding to N50 max========
    velvetSta = xlwt.Workbook()
    sheetVelvetSta = velvetSta.add_sheet('velvet_statistics', cell_overwrite_ok=True)
    kmerNum = (kmerMax-kmerMin)/2  #(35-31)/2=2
    row0 = ['sampleNums','kmer','medianCovDepth','nodesNum','contigN50','contigMax','genomeLen','usedReads','totalReads']
    for i in range(0,len(row0)):
        sheetVelvetSta.write(0,i,row0[i])  #first row

    #first cow is sample_list
    #second cow is kmer (kmerMin to kmerMax , kmerStep)

    for i in range(0,len(sample_list)):
        sheetVelvetSta.write(kmerNum*i+1,0,sample_list[i])   #first column is the sampleNums

    colKmer = range(kmerMin,kmerMax,kmerStep)
    for i in range(0,len(sample_list)*kmerNum,kmerNum): #i : every column
        for k in range(0,kmerNum):
            sheetVelvetSta.write(i+k+1,1,colKmer[k])   #second column is the k-mer
            k += 1

    logLast2Str = []
    MedianCovDepth = []

    logLast1Str = []
    nodesNum = []
    contigN50 = []
    contigMax = []
    genomeLen = []
    usedReads = []
    totalReads = []
    for i in range(len(sample_list)):
        for kmer in range(kmerMin,kmerMax,kmerStep):
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_list[i],kmer)
            Log_file = os.path.join(eachKmerDir,'Log')
            logs = open(Log_file,'r').readlines()
            #logLast2line = logs[len(logs)-2]
            #if logLast2line != '':
            #Median coverage depth = 96.434783  #not the logs[26] , to avoid the task again with adding the log_information to Log File
            #logLast2Str = logLast2line.split(' ')
            #print logLast2Str
            #print logLast2Str[4][:-2]
            #MedianCovDepth.append(string.atof(logLast2Str[4][:-2]))
            logLast1line = logs[len(logs)-1]  #Final graph has 199 nodes and n50 of 47325, max 166302, total 1570595, using 1470371/1487666 reads
            logLast1Str = logLast1line.split(' ')
            #['Final', 'graph', 'has', '172', 'nodes', 'and', 'n50', 'of', '46917,', 'max', '237233,', 'total', '1570012,', 'using', '1478778/1487666', 'reads\n']
            #print logLast1Str
            nodesNum.append(string.atoi(logLast1Str[3], 10))
            contigN50.append(string.atoi(logLast1Str[8][:-1],10))
            contigMax.append(string.atoi(logLast1Str[10][:-1],10))
            genomeLen.append(string.atoi(logLast1Str[12][:-1],10))
            usedReads.append(string.atoi(logLast1Str[14].split('/')[0], 10))
            totalReads.append(string.atoi(logLast1Str[14].split('/')[1], 10))
    #print contigN50
    for i in range(0,len(sample_list)*kmerNum):
        #sheetVelvetSta.write(i+1,2,MedianCovDepth[i])
        sheetVelvetSta.write(i+1,3,nodesNum[i])
        sheetVelvetSta.write(i+1,4,contigN50[i])
        sheetVelvetSta.write(i+1,5,contigMax[i])
        sheetVelvetSta.write(i+1,6,genomeLen[i])
        sheetVelvetSta.write(i+1,7,usedReads[i])
        sheetVelvetSta.write(i+1,8,totalReads[i])

    velvetSta.save(os.path.join(result_total_sta_path,'velvet_log_sta.xls'))

    contigN50_max_value_list = []  #all_samples
    contigN50_max_index_list = []
    for i in range(len(sample_list)):
        each_sample_contigN50_list = contigN50[kmerNum*i:kmerNum*(i+1)]
        #print each_sample_contigN50_list
        contigN50_max_value_list.append(max(each_sample_contigN50_list))
        contigN50_max_index_list.append(kmerMin+(each_sample_contigN50_list.index(max(each_sample_contigN50_list)))*kmerStep)

    #print contigN50_max_value_list
    #print contigN50_max_index_list

    print '=============blast:use the contigs.fa to map AGINST NT_viruses'
    #1.first:format the local nt fasta_database file
    #os.system('makeblastdb -dbtype nucl -input_type fasta -in %s' %(nt_path))  ###big memory : server

    #after makeblastdb,generate some other indexed files besides nt file,
    #use the indexed nt file(with some other indexed files) as the following input ref db '-db'.

    #2.second:blastn

    timer_blastn_start = time.clock()
    timer_blastn_start1 = time.time()

    #-outfmt 6
    for i in range(len(sample_list)):
        kmer = contigN50_max_index_list[i]
        eachKmerDir = '%s%s_%d'%(velvet_path,sample_list[i],kmer) #/public1/home/yefq/project_adenovirus/results/velvet_result/jiayan_31
        #print eachKmerDir
        eachContigFile = os.path.join(eachKmerDir,'contigs.fa')
        #print eachContigFile
        #estimate the contigs.fa is empty or not.
        if os.path.getsize(eachContigFile) != 0:
            eachBlastResult = os.path.join(blast_path,sample_list[i]+'.out')
            os.system('blastn -db %s -query %s -outfmt 6 -dust no -num_threads %d -perc_identity 80 -word_size 20 -max_target_seqs %d -evalue 0.0000001 -out %s'
                %(nt_viruses_family_file,eachContigFile,num_threads,max_target_seqs_nt_viruses,eachBlastResult))

    timer_blastn_end = time.clock()
    timer_blastn_end1 = time.time()
    blastn_time = timer_blastn_end - timer_blastn_start
    blastn_time1 = timer_blastn_end1 - timer_blastn_start1
    '''
    print '@@@@@----------timer:blastn_time_clock:'
    print blastn_time  #0.287237
    print '@@@@@----------timer:blastn_time_time:'
    print blastn_time1  #30.0745470524
    '''
    sheetPipelinerStaClock.write(1,3,blastn_time)
    sheetPipelinerStaTime.write(1,3,blastn_time1)

    ##3.third:sorted the blast result, use the different colomn
    for i in range(len(colNumSorted_list)):
        sort_blast_result_by_blastnColumn(colNumSorted_list[i])


    time_total_end = time.clock()
    time_total_end1 = time.time()
    

    total_time = time_total_end - time_total_start
    total_time1 = time_total_end1 - time_total_start1
    '''
    print '@@@@@-----------timer:total_time_clock:'
    print total_time  #23.519992
    print '@@@@@-----------timer:total_time_time:'
    print total_time1  #2177.6757071
    '''
    sheetPipelinerStaClock.write(1,4,total_time)
    sheetPipelinerStaTime.write(1,4,total_time1)

    pipelinerSta.save(os.path.join(result_total_sta_path,'pipeline_time_sta.xls'))

##method : mapAginstHumanFirstwithoutQC
def process_pipeliner_se_mapAginstHumanFirstwithoutQC(genomeHostName,ntVirusesFamilyName):
    ##rebuild : use all genome host
    genome_host_path = genomesDB_path+genomeHostName
    genome_host_file = os.path.join(genome_host_path,genomeHostName+'.fa')

    genome_host_bt2_index_base_name = 'bt2_index_%s'%(genomeHostName)
    
    ##rebuild : use ntVirusesFamilyName to replace nt_viruses_family_name_list[i], used for mapping with bowtie2 on nt_viruses_ and blastn
    family_path = 'nt_%s'%(ntVirusesFamilyName)
    family_file_name = 'nt_%s.fasta'%(ntVirusesFamilyName)
    
    nt_viruses_family_path = blastDB_path+family_path
    nt_viruses_family_file = os.path.join(nt_viruses_family_path,family_file_name)
        
    nt_virus_bt2_index_base_name = 'bt2_index_%s'%(ntVirusesFamilyName)

    nt_viruses_family_index_header_file_name = 'nt_%s_index_title'%(ntVirusesFamilyName)  #0_index_title / Adenoviridae_index_title
    nt_viruses_family_index_header_file = os.path.join(nt_viruses_family_path,nt_viruses_family_index_header_file_name)

    
    #print sample_se_list

    time_total_start = time.clock()
    time_total_start1 = time.time()

    pipelinerSta = xlwt.Workbook()
    sheetPipelinerStaClock = pipelinerSta.add_sheet('pipeline_clock_statistics',cell_overwrite_ok=True)
    sheetPipelinerStaTime = pipelinerSta.add_sheet('pipeline_time_statistics',cell_overwrite_ok=True)
    
    rowTitle = ['map_%s'%genomeHostName,'map_%s'%ntVirusesFamilyName,'assembly_velvet','map_blastn','total_time']
    for i in range(0,len(rowTitle)):
        sheetPipelinerStaClock.write(0,i,rowTitle[i]) #first row
        sheetPipelinerStaTime.write(0,i,rowTitle[i])

    

    print '==========map and align : AGINST NT_viruses reference,keep mapped reads'

    #bowtie2IndexNTVirusFamily = nt_viruses_family_name_list[0]
    #bowtie2IndexNTVirusFamily = bt2_index_base_name_list[0]  #'bt2_index_Adenoviridae'

    timer_bowtie2_nt_virus_start = time.clock()
    timer_bowtie2_nt_virus_start1 = time.time()

    afterNTVIRUSfastqFiles =[]
    for i in range(len(sample_se_list)):
        outputSAMFileNTVirus = bowtie2_path+sample_se_list[i]+map_nt_viruses_singleAndunpair_sam_file_suffix
        outputFASTQFileNTVirus = bowtie2_path+sample_se_list[i]+map_nt_viruses_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -U %s --no-head --no-unal -S %s --al %s -p %d'
            %(nt_viruses_family_path,nt_virus_bt2_index_base_name,samplesRawDataSingleFiles[i],outputSAMFileNTVirus,outputFASTQFileNTVirus,num_threads))

        afterNTVIRUSfastqFiles.append(outputFASTQFileNTVirus)


    timer_bowtie2_nt_virus_end = time.clock()
    timer_bowtie2_nt_virus_end1 = time.time()
    bowtie2_nt_virus_time = timer_bowtie2_nt_virus_end - timer_bowtie2_nt_virus_start
    bowtie2_nt_virus_time1 = timer_bowtie2_nt_virus_end1 - timer_bowtie2_nt_virus_start1
    '''
    print '@@@@@----------timer:bowtie2_nt_virus_time_clock:'
    print bowtie2_nt_virus_time    #18.209088
    print '@@@@@----------timer:bowtie2_nt_virus_time_time:'
    print bowtie2_nt_virus_time1    #1678.66673303
    '''
    sheetPipelinerStaClock.write(1,0,bowtie2_nt_virus_time)
    sheetPipelinerStaTime.write(1,0,bowtie2_nt_virus_time1)

    print '==========map and align : AGINST hg38 reference,keep unmapped reads'

    timer_bowtie2_hg38_start = time.clock()
    timer_bowtie2_hg38_start1 = time.time()

    afterHGfastqFiles =[]
    for i in range(len(sample_se_list)):
        outputSAMFileGBHG = bowtie2_path+sample_se_list[i]+'_host.sam'
        outputFASTQFileGBHG = bowtie2_path+sample_se_list[i]+map_gb_hg_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -U %s --no-head --no-unal -S %s --un %s -p %d'
            %(genome_host_path,genome_host_bt2_index_base_name,afterNTVIRUSfastqFiles[i],outputSAMFileGBHG,outputFASTQFileGBHG,num_threads))
        
        afterHGfastqFiles.append(outputFASTQFileGBHG)



    timer_bowtie2_hg38_end = time.clock()
    timer_bowtie2_hg38_end1 = time.time()
    bowtie2_hg38_time = timer_bowtie2_hg38_end - timer_bowtie2_hg38_start
    bowtie2_hg38_time1 = timer_bowtie2_hg38_end1 - timer_bowtie2_hg38_start1
    '''
    print '@@@@@----------timer:bowtie2_hg38_time_clock:'
    print bowtie2_hg38_time    #17.103814
    print '@@@@@----------timer:bowtie2_hg38_time_time:'
    print bowtie2_hg38_time1    #1568.48754096
    '''
    sheetPipelinerStaClock.write(1,1,bowtie2_hg38_time)
    sheetPipelinerStaTime.write(1,1,bowtie2_hg38_time1)

    print '==========assembly : velvet with alignment data above'

    timer_velvet_start = time.clock()
    timer_velvet_start1 = time.time()

    minContigLength = 100
    for i in range(len(sample_se_list)):
        for kmer in range(kmerMin,kmerMax,kmerStep):
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_se_list[i],kmer)
            os.system('velveth %s %d -short -fastq %s' %(eachKmerDir,kmer,afterHGfastqFiles[i]))
            os.system('velvetg %s -exp_cov auto -cov_cutoff auto -min_contig_lgth %d' %(eachKmerDir,minContigLength))

            eachKmerFiles = os.listdir(eachKmerDir)
            eachKmerFilesAP = []
            for j in range(len(eachKmerFiles)):
                eachKmerFilesAP.append(eachKmerDir+'/'+eachKmerFiles[j])
            for j in range(len(eachKmerFilesAP)):
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Graph2' : os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/LastGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/PreGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Roadmaps': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Sequences': os.remove(eachKmerFilesAP[j])

    timer_velvet_end = time.clock()
    timer_velvet_end1 = time.time()
    velvet_time = timer_velvet_end - timer_velvet_start
    velvet_time1 = timer_velvet_end1 - timer_velvet_start1
    '''
    print '@@@@@----------timer:velvet_time_clock:'
    print velvet_time    #0.034675
    print '@@@@@----------timer:velvet_time_time:'
    print velvet_time1    #1.79084706306
    '''
    sheetPipelinerStaClock.write(1,2,velvet_time)
    sheetPipelinerStaTime.write(1,2,velvet_time1)

    ##================statistics velvet : pick up the proper kmer , accoding to N50 max========
    velvetSta = xlwt.Workbook()
    sheetVelvetSta = velvetSta.add_sheet('velvet_statistics', cell_overwrite_ok=True)
    kmerNum = (kmerMax-kmerMin)/2  #(35-31)/2=2
    row0 = ['sampleNums','kmer','medianCovDepth','nodesNum','contigN50','contigMax','genomeLen','usedReads','totalReads']
    for i in range(0,len(row0)):
        sheetVelvetSta.write(0,i,row0[i])  #first row

    #first cow is sample_se_list
    #second cow is kmer (kmerMin to kmerMax , kmerStep)

    for i in range(0,len(sample_se_list)):
        sheetVelvetSta.write(kmerNum*i+1,0,sample_se_list[i])   #first column is the sampleNums

    colKmer = range(kmerMin,kmerMax,kmerStep)
    for i in range(0,len(sample_se_list)*kmerNum,kmerNum): #i : every column
        for k in range(0,kmerNum):
            sheetVelvetSta.write(i+k+1,1,colKmer[k])   #second column is the k-mer
            k += 1

    logLast2Str = []
    MedianCovDepth = []

    logLast1Str = []
    nodesNum = []
    contigN50 = []
    contigMax = []
    genomeLen = []
    usedReads = []
    totalReads = []
    for i in range(len(sample_se_list)):
        for kmer in range(kmerMin,kmerMax,kmerStep):
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_se_list[i],kmer)
            Log_file = os.path.join(eachKmerDir,'Log')
            logs = open(Log_file,'r').readlines()
            #logLast2line = logs[len(logs)-2]
            #if logLast2line != '':
            #Median coverage depth = 96.434783  #not the logs[26] , to avoid the task again with adding the log_information to Log File
            #logLast2Str = logLast2line.split(' ')
            #print logLast2Str
            #print logLast2Str[4][:-2]
            #MedianCovDepth.append(string.atof(logLast2Str[4][:-2]))
            logLast1line = logs[len(logs)-1]  #Final graph has 199 nodes and n50 of 47325, max 166302, total 1570595, using 1470371/1487666 reads
            logLast1Str = logLast1line.split(' ')
            #['Final', 'graph', 'has', '172', 'nodes', 'and', 'n50', 'of', '46917,', 'max', '237233,', 'total', '1570012,', 'using', '1478778/1487666', 'reads\n']
            #print logLast1Str
            nodesNum.append(string.atoi(logLast1Str[3], 10))
            contigN50.append(string.atoi(logLast1Str[8][:-1],10))
            contigMax.append(string.atoi(logLast1Str[10][:-1],10))
            genomeLen.append(string.atoi(logLast1Str[12][:-1],10))
            usedReads.append(string.atoi(logLast1Str[14].split('/')[0], 10))
            totalReads.append(string.atoi(logLast1Str[14].split('/')[1], 10))
    #print contigN50
    for i in range(0,len(sample_se_list)*kmerNum):
        #sheetVelvetSta.write(i+1,2,MedianCovDepth[i])
        sheetVelvetSta.write(i+1,3,nodesNum[i])
        sheetVelvetSta.write(i+1,4,contigN50[i])
        sheetVelvetSta.write(i+1,5,contigMax[i])
        sheetVelvetSta.write(i+1,6,genomeLen[i])
        sheetVelvetSta.write(i+1,7,usedReads[i])
        sheetVelvetSta.write(i+1,8,totalReads[i])

    velvetSta.save(os.path.join(result_total_sta_path,'velvet_log_sta.xls'))

    ##check the max contigN50 to pick up ,then blastn
    contigN50_max_value_list = []  #all_samples
    contigN50_max_index_list = []
    for i in range(len(sample_se_list)):
        each_sample_contigN50_list = contigN50[kmerNum*i:kmerNum*(i+1)]
        #print each_sample_contigN50_list
        contigN50_max_value_list.append(max(each_sample_contigN50_list))
        contigN50_max_index_list.append(kmerMin+(each_sample_contigN50_list.index(max(each_sample_contigN50_list)))*kmerStep)

    #print contigN50_max_value_list
    #print contigN50_max_index_list

    print '=============blast:use the contigs.fa to map AGINST NT_viruses'
    #1.first:format the local nt fasta_database file
    #os.system('makeblastdb -dbtype nucl -input_type fasta -in %s' %(nt_path))  ###big memory : server

    #after makeblastdb,generate some other indexed files besides nt file,
    #use the indexed nt file(with some other indexed files) as the following input ref db '-db'.

    #2.second:blastn

    timer_blastn_start = time.clock()
    timer_blastn_start1 = time.time()

                 
    for i in range(len(sample_se_list)):
        kmer = contigN50_max_index_list[i]
        eachKmerDir = '%s%s_%d'%(velvet_path,sample_se_list[i],kmer) #/public1/home/yefq/project_adenovirus/results/velvet_result/jiayan_31
        #print eachKmerDir
        eachContigFile = os.path.join(eachKmerDir,'contigs.fa')
        #print eachContigFile
        if os.path.getsize(eachContigFile) != 0:
            eachBlastResult = os.path.join(blast_path,sample_se_list[i]+'.out')
            os.system('blastn -db %s -query %s -outfmt 6 -dust no -num_threads %d -perc_identity 80 -word_size 20 -max_target_seqs %d -evalue 0.0000001 -out %s'
                %(nt_viruses_family_file,eachContigFile,num_threads,max_target_seqs_nt_viruses,eachBlastResult))

    timer_blastn_end = time.clock()
    timer_blastn_end1 = time.time()
    blastn_time = timer_blastn_end - timer_blastn_start
    blastn_time1 = timer_blastn_end1 - timer_blastn_start1
    '''
    print '@@@@@----------timer:blastn_time_clock:'
    print blastn_time  #0.287237
    print '@@@@@----------timer:blastn_time_time:'
    print blastn_time1  #30.0745470524
    '''
    sheetPipelinerStaClock.write(1,3,blastn_time)
    sheetPipelinerStaTime.write(1,3,blastn_time1)

    ##3.third:sorted the blast result, use the different colomn
    for i in range(len(colNumSorted_list)):
        sort_blast_result_by_blastnColumn(colNumSorted_list[i])


    time_total_end = time.clock()
    time_total_end1 = time.time()
    

    total_time = time_total_end - time_total_start
    total_time1 = time_total_end1 - time_total_start1
    '''
    print '@@@@@-----------timer:total_time_clock:'
    print total_time  #23.519992
    print '@@@@@-----------timer:total_time_time:'
    print total_time1  #2177.6757071
    '''
    sheetPipelinerStaClock.write(1,4,total_time)
    sheetPipelinerStaTime.write(1,4,total_time1)

    pipelinerSta.save(os.path.join(result_total_sta_path,'pipeline_time_sta.xls'))

##method : mapAginstHumanFirst
def process_pipeliner_pe_mapAginstHumanFirst(genomeHostName,ntVirusesFamilyName):
    ##rebuild : use all genome host
    genome_host_path = genomesDB_path+genomeHostName
    genome_host_file = os.path.join(genome_host_path,genomeHostName+'.fa')

    genome_host_bt2_index_base_name = 'bt2_index_%s'%(genomeHostName)
    
    ##rebuild : use ntVirusesFamilyName to replace nt_viruses_family_name_list[i], used for mapping with bowtie2 on nt_viruses_ and blastn
    family_path = 'nt_%s'%(ntVirusesFamilyName)
    family_file_name = 'nt_%s.fasta'%(ntVirusesFamilyName)
    
    nt_viruses_family_path = blastDB_path+family_path
    nt_viruses_family_file = os.path.join(nt_viruses_family_path,family_file_name)
        
    nt_virus_bt2_index_base_name = 'bt2_index_%s'%(ntVirusesFamilyName)

    nt_viruses_family_index_header_file_name = 'nt_%s_index_title'%(ntVirusesFamilyName)  #0_index_title / Adenoviridae_index_title
    nt_viruses_family_index_header_file = os.path.join(nt_viruses_family_path,nt_viruses_family_index_header_file_name)

    #print sample_list

    time_total_start = time.clock()
    time_total_start1 = time.time()

    pipelinerSta = xlwt.Workbook()
    sheetPipelinerStaClock = pipelinerSta.add_sheet('pipeline_clock_statistics',cell_overwrite_ok=True)
    sheetPipelinerStaTime = pipelinerSta.add_sheet('pipeline_time_statistics',cell_overwrite_ok=True)
    #rowTitle = ['qc_trimmomatic','map_hp38','map_genebank_virus','map_nt_virus','assembly_velvet','map_blastn','total_time']
    #rowTitle = ['qc_trimmomatic','map_hp38','map_nt_virus','assembly_velvet','map_blastn','total_time']
    rowTitle = ['qc_trimmomatic','map_%s'%genomeHostName,'map_%s'%ntVirusesFamilyName,'assembly_velvet','map_blastn','total_time']
    for i in range(0,len(rowTitle)):
        sheetPipelinerStaClock.write(0,i,rowTitle[i]) #first row
        sheetPipelinerStaTime.write(0,i,rowTitle[i])


    
    print '==========qc:trimmomatic'
    timer_trimmomatic_start = time.clock()
    timer_trimmomatic_start1 = time.time()

    trimmomatic_results_prefix = []
    for i in range(len(sample_list)):
        trimmomatic_result_prefix = trimmomatic_path+sample_list[i]
        trimmomatic_results_prefix.append(trimmomatic_result_prefix)
    #print trimmomatic_results_prefix ##/public1/home/yefq/project_adenovirus/results/trimmomatic_result/jiayan

    outputFPFiles = []
    outputFUFiles = []
    outputRPFiles = []
    outputRUFiles = []
    for i in range(len(sample_list)):
        #print sample_list[i]
        trimlogFile=trimmomatic_results_prefix[i]+'_trimlog.txt'
        outputFPFile=trimmomatic_results_prefix[i]+'_output_forward_paired.fastq'
        outputFUFile=trimmomatic_results_prefix[i]+'_output_forward_unpaired.fastq'
        outputRPFile=trimmomatic_results_prefix[i]+'_output_reverse_paired.fastq'
        outputRUFile=trimmomatic_results_prefix[i]+'_output_reverse_unpaired.fastq'

        os.system('java -jar /usr/share/java/trimmomatic.jar PE -threads %d -trimlog %s %s %s %s %s %s %s ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'
        %(num_threads,trimlogFile,samplesRawDataForwardFiles[i],samplesRawDataReverseFiles[i],outputFPFile,outputFUFile,outputRPFile,outputRUFile))

        outputFPFiles.append(outputFPFile)
        outputRPFiles.append(outputRPFile)
        outputFUFiles.append(outputFUFile)
        outputRUFiles.append(outputRUFile)

    timer_trimmomatic_end = time.clock()
    timer_trimmomatic_end1 = time.time()
    trimmomatic_time = timer_trimmomatic_end - timer_trimmomatic_start
    trimmomatic_time1 = timer_trimmomatic_end1 - timer_trimmomatic_start1
    '''
    print '@@@@@----------timer:trimmomatic_time_clock:'
    print trimmomatic_time  #4.822289
    print '@@@@@----------timer:trimmomatic_time_time:'
    print trimmomatic_time1  #466.976895094
    '''
    sheetPipelinerStaClock.write(1,0,trimmomatic_time)
    sheetPipelinerStaTime.write(1,0,trimmomatic_time1)


    print '==========map and align : AGINST hg38 reference,keep unmapped reads'

    timer_bowtie2_hg38_start = time.clock()
    timer_bowtie2_hg38_start1 = time.time()

    outputFastqPair1FileGBHG_list = []
    outputFastqPair2FileGBHG_list = []
    outputFastqUnpairFileGBHG_list = []

    for i in range(len(sample_list)):
        outputSAMFileNTGBHG = bowtie2_path+sample_list[i]+'_host.sam'
        outputFastqPairFileGBHG = bowtie2_path+sample_list[i]+'_pair_host.fastq'
        outputFastqPair1FileGBHG = bowtie2_path+sample_list[i]+'_pair_host.1.fastq'
        outputFastqPair2FileGBHG = bowtie2_path+sample_list[i]+'_pair_host.2.fastq'
        outputFastqUnpairFileGBHG = bowtie2_path+sample_list[i]+map_gb_hg_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -1 %s -2 %s -U %s,%s --no-head --no-unal -S %s --un-conc %s --un %s -p %d'
        %(genome_host_path,genome_host_bt2_index_base_name,outputFPFiles[i],outputRPFiles[i],outputFUFiles[i],outputRUFiles[i],outputSAMFileNTGBHG,outputFastqPairFileGBHG,outputFastqUnpairFileGBHG,num_threads))

        outputFastqPair1FileGBHG_list.append(outputFastqPair1FileGBHG)
        outputFastqPair2FileGBHG_list.append(outputFastqPair2FileGBHG)
        outputFastqUnpairFileGBHG_list.append(outputFastqUnpairFileGBHG)

    timer_bowtie2_hg38_end = time.clock()
    timer_bowtie2_hg38_end1 = time.time()
    bowtie2_hg38_time = timer_bowtie2_hg38_end - timer_bowtie2_hg38_start
    bowtie2_hg38_time1 = timer_bowtie2_hg38_end1 - timer_bowtie2_hg38_start1
    '''
    print '@@@@@----------timer:bowtie2_hg38_time_clock:'
    print bowtie2_hg38_time    #17.103814
    print '@@@@@----------timer:bowtie2_hg38_time_time:'
    print bowtie2_hg38_time1    #1568.48754096
    '''
    sheetPipelinerStaClock.write(1,1,bowtie2_hg38_time)
    sheetPipelinerStaTime.write(1,1,bowtie2_hg38_time1)


    print '==========map and align : AGINST NT_viruses reference,keep mapped reads'

    #bowtie2IndexNTVirusFamily = nt_viruses_family_name_list[0]
    #bowtie2IndexNTVirusFamily = bt2_index_base_name_list[0]  #'bt2_index_Adenoviridae'

    timer_bowtie2_nt_virus_start = time.clock()
    timer_bowtie2_nt_virus_start1 = time.time()

    outputFastqPair1FileNTVirus_list = []
    outputFastqPair2FileNTVirus_list = []
    outputFastqUnpairFileNTVirus_list = []
    for i in range(len(sample_list)):
        outputSAMFileNTVirus = bowtie2_path+sample_list[i]+map_nt_viruses_singleAndunpair_sam_file_suffix
        outputFastqPairFileNTVirus = bowtie2_path+sample_list[i]+'_pair_VIRUS.fastq'
        outputFastqPair1FileNTVirus = bowtie2_path+sample_list[i]+'_pair_VIRUS.1.fastq'
        outputFastqPair2FileNTVirus = bowtie2_path+sample_list[i]+'_pair_VIRUS.2.fastq'
        
        outputFastqUnpairFileNTVirus = bowtie2_path+sample_list[i]+map_nt_viruses_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -1 %s -2 %s -U %s --no-head --no-unal -S %s --al-conc %s --al %s -p %d'
            %(nt_viruses_family_path,nt_virus_bt2_index_base_name,outputFastqPair1FileGBHG_list[i],outputFastqPair2FileGBHG_list[i],outputFastqUnpairFileGBHG_list[i],outputSAMFileNTVirus,outputFastqPairFileNTVirus,outputFastqUnpairFileNTVirus,num_threads))

        outputFastqPair1FileNTVirus_list.append(outputFastqPair1FileNTVirus)
        outputFastqPair2FileNTVirus_list.append(outputFastqPair2FileNTVirus)
        outputFastqUnpairFileNTVirus_list.append(outputFastqUnpairFileNTVirus)
        
    timer_bowtie2_nt_virus_end = time.clock()
    timer_bowtie2_nt_virus_end1 = time.time()
    bowtie2_nt_virus_time = timer_bowtie2_nt_virus_end - timer_bowtie2_nt_virus_start
    bowtie2_nt_virus_time1 = timer_bowtie2_nt_virus_end1 - timer_bowtie2_nt_virus_start1
    '''
    print '@@@@@----------timer:bowtie2_nt_virus_time_clock:'
    print bowtie2_nt_virus_time    #18.209088
    print '@@@@@----------timer:bowtie2_nt_virus_time_time:'
    print bowtie2_nt_virus_time1    #1678.66673303
    '''
    sheetPipelinerStaClock.write(1,2,bowtie2_nt_virus_time)
    sheetPipelinerStaTime.write(1,2,bowtie2_nt_virus_time1)


    print '==========assembly : velvet with alignment data above'

    timer_velvet_start = time.clock()
    timer_velvet_start1 = time.time()

    for i in range(len(sample_list)):
        #inputPE = velvet_result+sample_list[i]+'_pe.fastq'
        #os.system('%s %s %s %s' %(shuffleSequences_fastq_path,outputFastq1AfterHGandVirus[i],outputFastq2AfterHGandVirus[i],inputPE))
        for kmer in range(kmerMin,kmerMax,kmerStep):
            #print sample_list[i]+'************************************'
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_list[i],kmer)
            os.system('velveth %s %d -shortPaired -fastq -separate %s %s -short -fastq %s'%(eachKmerDir,kmer,outputFastqPair1FileNTVirus_list[i],outputFastqPair2FileNTVirus_list[i],outputFastqUnpairFileNTVirus_list[i]))
            os.system('velvetg %s -exp_cov auto -cov_cutoff auto ' %(eachKmerDir))

            eachKmerFiles = os.listdir(eachKmerDir)
            eachKmerFilesAP = []
            for j in range(len(eachKmerFiles)):
                eachKmerFilesAP.append(eachKmerDir+'/'+eachKmerFiles[j])
            for j in range(len(eachKmerFilesAP)):
                #pass
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Graph2' : os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/LastGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/PreGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Roadmaps': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Sequences': os.remove(eachKmerFilesAP[j])

    timer_velvet_end = time.clock()
    timer_velvet_end1 = time.time()
    velvet_time = timer_velvet_end - timer_velvet_start
    velvet_time1 = timer_velvet_end1 - timer_velvet_start1
    '''
    print '@@@@@----------timer:velvet_time_clock:'
    print velvet_time    #0.034675
    print '@@@@@----------timer:velvet_time_time:'
    print velvet_time1    #1.79084706306
    '''
    sheetPipelinerStaClock.write(1,3,velvet_time)
    sheetPipelinerStaTime.write(1,3,velvet_time1)

    ##================statistics velvet : pick up the proper kmer , accoding to N50 max========
    velvetSta = xlwt.Workbook()
    sheetVelvetSta = velvetSta.add_sheet('velvet_statistics', cell_overwrite_ok=True)
    kmerNum = (kmerMax-kmerMin)/2  #(35-31)/2=2
    row0 = ['sampleNums','kmer','medianCovDepth','nodesNum','contigN50','contigMax','genomeLen','usedReads','totalReads']
    for i in range(0,len(row0)):
        sheetVelvetSta.write(0,i,row0[i])  #first row

    #first cow is sample_list
    #second cow is kmer (kmerMin to kmerMax , kmerStep)

    for i in range(0,len(sample_list)):
        sheetVelvetSta.write(kmerNum*i+1,0,sample_list[i])   #first column is the sampleNums

    colKmer = range(kmerMin,kmerMax,kmerStep)
    for i in range(0,len(sample_list)*kmerNum,kmerNum): #i : every column
        for k in range(0,kmerNum):
            sheetVelvetSta.write(i+k+1,1,colKmer[k])   #second column is the k-mer
            k += 1

    logLast2Str = []
    MedianCovDepth = []

    logLast1Str = []
    nodesNum = []
    contigN50 = []
    contigMax = []
    genomeLen = []
    usedReads = []
    totalReads = []
    for i in range(len(sample_list)):
        for kmer in range(kmerMin,kmerMax,kmerStep):
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_list[i],kmer)
            Log_file = os.path.join(eachKmerDir,'Log')
            logs = open(Log_file,'r').readlines()
            #logLast2line = logs[len(logs)-2]
            #if logLast2line != '':
            #Median coverage depth = 96.434783  #not the logs[26] , to avoid the task again with adding the log_information to Log File
            #logLast2Str = logLast2line.split(' ')
            #print logLast2Str
            #print logLast2Str[4][:-2]
            #MedianCovDepth.append(string.atof(logLast2Str[4][:-2]))
            logLast1line = logs[len(logs)-1]  #Final graph has 199 nodes and n50 of 47325, max 166302, total 1570595, using 1470371/1487666 reads
            logLast1Str = logLast1line.split(' ')
            #['Final', 'graph', 'has', '172', 'nodes', 'and', 'n50', 'of', '46917,', 'max', '237233,', 'total', '1570012,', 'using', '1478778/1487666', 'reads\n']
            #print logLast1Str
            nodesNum.append(string.atoi(logLast1Str[3], 10))
            contigN50.append(string.atoi(logLast1Str[8][:-1],10))
            contigMax.append(string.atoi(logLast1Str[10][:-1],10))
            genomeLen.append(string.atoi(logLast1Str[12][:-1],10))
            usedReads.append(string.atoi(logLast1Str[14].split('/')[0], 10))
            totalReads.append(string.atoi(logLast1Str[14].split('/')[1], 10))
    #print contigN50
    for i in range(0,len(sample_list)*kmerNum):
        #sheetVelvetSta.write(i+1,2,MedianCovDepth[i])
        sheetVelvetSta.write(i+1,3,nodesNum[i])
        sheetVelvetSta.write(i+1,4,contigN50[i])
        sheetVelvetSta.write(i+1,5,contigMax[i])
        sheetVelvetSta.write(i+1,6,genomeLen[i])
        sheetVelvetSta.write(i+1,7,usedReads[i])
        sheetVelvetSta.write(i+1,8,totalReads[i])

    velvetSta.save(os.path.join(result_total_sta_path,'velvet_log_sta.xls'))

    contigN50_max_value_list = []  #all_samples
    contigN50_max_index_list = []
    for i in range(len(sample_list)):
        each_sample_contigN50_list = contigN50[kmerNum*i:kmerNum*(i+1)]
        #print each_sample_contigN50_list
        contigN50_max_value_list.append(max(each_sample_contigN50_list))
        contigN50_max_index_list.append(kmerMin+(each_sample_contigN50_list.index(max(each_sample_contigN50_list)))*kmerStep)

    #print contigN50_max_value_list
    #print contigN50_max_index_list

    print '=============blast:use the contigs.fa to map AGINST NT_viruses'
    #1.first:format the local nt fasta_database file
    #os.system('makeblastdb -dbtype nucl -input_type fasta -in %s' %(nt_path))  ###big memory : server

    #after makeblastdb,generate some other indexed files besides nt file,
    #use the indexed nt file(with some other indexed files) as the following input ref db '-db'.

    #2.second:blastn

    timer_blastn_start = time.clock()
    timer_blastn_start1 = time.time()

    #-outfmt 6
    for i in range(len(sample_list)):
        kmer = contigN50_max_index_list[i]
        eachKmerDir = '%s%s_%d'%(velvet_path,sample_list[i],kmer) #/public1/home/yefq/project_adenovirus/results/velvet_result/jiayan_31
        #print eachKmerDir
        eachContigFile = os.path.join(eachKmerDir,'contigs.fa')
        #print eachContigFile
        #estimate the contigs.fa is empty or not.
        if os.path.getsize(eachContigFile) != 0:
            eachBlastResult = os.path.join(blast_path,sample_list[i]+'.out')
            os.system('blastn -db %s -query %s -outfmt 6 -dust no -num_threads %d -perc_identity 80 -word_size 20 -max_target_seqs %d -evalue 0.0000001 -out %s'
                %(nt_viruses_family_file,eachContigFile,num_threads,max_target_seqs_nt_viruses,eachBlastResult))

    timer_blastn_end = time.clock()
    timer_blastn_end1 = time.time()
    blastn_time = timer_blastn_end - timer_blastn_start
    blastn_time1 = timer_blastn_end1 - timer_blastn_start1
    '''
    print '@@@@@----------timer:blastn_time_clock:'
    print blastn_time  #0.287237
    print '@@@@@----------timer:blastn_time_time:'
    print blastn_time1  #30.0745470524
    '''
    sheetPipelinerStaClock.write(1,4,blastn_time)
    sheetPipelinerStaTime.write(1,4,blastn_time1)

    ##3.third:sorted the blast result, use the different colomn
    for i in range(len(colNumSorted_list)):
        sort_blast_result_by_blastnColumn(colNumSorted_list[i])


    time_total_end = time.clock()
    time_total_end1 = time.time()
    

    total_time = time_total_end - time_total_start
    total_time1 = time_total_end1 - time_total_start1
    '''
    print '@@@@@-----------timer:total_time_clock:'
    print total_time  #23.519992
    print '@@@@@-----------timer:total_time_time:'
    print total_time1  #2177.6757071
    '''
    sheetPipelinerStaClock.write(1,5,total_time)
    sheetPipelinerStaTime.write(1,5,total_time1)

    pipelinerSta.save(os.path.join(result_total_sta_path,'pipeline_time_sta.xls'))

##method : mapAginstHumanFirst
def process_pipeliner_se_mapAginstHumanFirst(genomeHostName,ntVirusesFamilyName):
    ##rebuild : use all genome host
    genome_host_path = genomesDB_path+genomeHostName
    genome_host_file = os.path.join(genome_host_path,genomeHostName+'.fa')

    genome_host_bt2_index_base_name = 'bt2_index_%s'%(genomeHostName)
    
    ##rebuild : use ntVirusesFamilyName to replace nt_viruses_family_name_list[i], used for mapping with bowtie2 on nt_viruses_ and blastn
    family_path = 'nt_%s'%(ntVirusesFamilyName)
    family_file_name = 'nt_%s.fasta'%(ntVirusesFamilyName)
    
    nt_viruses_family_path = blastDB_path+family_path
    nt_viruses_family_file = os.path.join(nt_viruses_family_path,family_file_name)
        
    nt_virus_bt2_index_base_name = 'bt2_index_%s'%(ntVirusesFamilyName)

    nt_viruses_family_index_header_file_name = 'nt_%s_index_title'%(ntVirusesFamilyName)  #0_index_title / Adenoviridae_index_title
    nt_viruses_family_index_header_file = os.path.join(nt_viruses_family_path,nt_viruses_family_index_header_file_name)

    
    #print sample_se_list

    time_total_start = time.clock()
    time_total_start1 = time.time()

    pipelinerSta = xlwt.Workbook()
    sheetPipelinerStaClock = pipelinerSta.add_sheet('pipeline_clock_statistics',cell_overwrite_ok=True)
    sheetPipelinerStaTime = pipelinerSta.add_sheet('pipeline_time_statistics',cell_overwrite_ok=True)
    
    rowTitle = ['qc_trimmomatic','map_%s'%genomeHostName,'map_%s'%ntVirusesFamilyName,'assembly_velvet','map_blastn','total_time']
    for i in range(0,len(rowTitle)):
        sheetPipelinerStaClock.write(0,i,rowTitle[i]) #first row
        sheetPipelinerStaTime.write(0,i,rowTitle[i])

    
    print '==========qc:trimmomatic'
    timer_trimmomatic_start = time.clock()
    timer_trimmomatic_start1 = time.time()
    
    MINLEN = 10
    LEADING = 0
    TRAILING = 0    
    outputSingleFiles = []
    for i in range(len(sample_se_list)):
        trimlogFile=trimmomatic_path+sample_se_list[i]+'_trimlog.txt'
        outputSingleFile=trimmomatic_path+sample_se_list[i]+'_output.fastq'

        os.system('java -jar /usr/share/java/trimmomatic.jar SE -threads %d -trimlog %s %s %s ILLUMINACLIP:/usr/share/trimmomatic/TruSeq3-SE.fa:2:30:10 LEADING:%d TRAILING:%d SLIDINGWINDOW:4:15 MINLEN:%d'
        %(num_threads,trimlogFile,samplesRawDataSingleFiles[i],outputSingleFile,LEADING,TRAILING,MINLEN))

        outputSingleFiles.append(outputSingleFile)

    timer_trimmomatic_end = time.clock()
    timer_trimmomatic_end1 = time.time()
    trimmomatic_time = timer_trimmomatic_end - timer_trimmomatic_start
    trimmomatic_time1 = timer_trimmomatic_end1 - timer_trimmomatic_start1
    '''
    print '@@@@@----------timer:trimmomatic_time_clock:'
    print trimmomatic_time  #4.822289
    print '@@@@@----------timer:trimmomatic_time_time:'
    print trimmomatic_time1  #466.976895094
    '''
    sheetPipelinerStaClock.write(1,0,trimmomatic_time)
    sheetPipelinerStaTime.write(1,0,trimmomatic_time1)


    print '==========map and align : AGINST NT_viruses reference,keep mapped reads'

    #bowtie2IndexNTVirusFamily = nt_viruses_family_name_list[0]
    #bowtie2IndexNTVirusFamily = bt2_index_base_name_list[0]  #'bt2_index_Adenoviridae'

    timer_bowtie2_nt_virus_start = time.clock()
    timer_bowtie2_nt_virus_start1 = time.time()

    afterNTVIRUSfastqFiles =[]
    for i in range(len(sample_se_list)):
        outputSAMFileNTVirus = bowtie2_path+sample_se_list[i]+map_nt_viruses_singleAndunpair_sam_file_suffix
        outputFASTQFileNTVirus = bowtie2_path+sample_se_list[i]+map_nt_viruses_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -U %s --no-head --no-unal -S %s --al %s -p %d'
            %(nt_viruses_family_path,nt_virus_bt2_index_base_name,outputSingleFiles[i],outputSAMFileNTVirus,outputFASTQFileNTVirus,num_threads))

        afterNTVIRUSfastqFiles.append(outputFASTQFileNTVirus)


    timer_bowtie2_nt_virus_end = time.clock()
    timer_bowtie2_nt_virus_end1 = time.time()
    bowtie2_nt_virus_time = timer_bowtie2_nt_virus_end - timer_bowtie2_nt_virus_start
    bowtie2_nt_virus_time1 = timer_bowtie2_nt_virus_end1 - timer_bowtie2_nt_virus_start1
    '''
    print '@@@@@----------timer:bowtie2_nt_virus_time_clock:'
    print bowtie2_nt_virus_time    #18.209088
    print '@@@@@----------timer:bowtie2_nt_virus_time_time:'
    print bowtie2_nt_virus_time1    #1678.66673303
    '''
    sheetPipelinerStaClock.write(1,1,bowtie2_nt_virus_time)
    sheetPipelinerStaTime.write(1,1,bowtie2_nt_virus_time1)

    print '==========map and align : AGINST hg38 reference,keep unmapped reads'

    timer_bowtie2_hg38_start = time.clock()
    timer_bowtie2_hg38_start1 = time.time()

    afterHGfastqFiles =[]
    for i in range(len(sample_se_list)):
        outputSAMFileGBHG = bowtie2_path+sample_se_list[i]+'_host.sam'
        outputFASTQFileGBHG = bowtie2_path+sample_se_list[i]+map_gb_hg_singleAndunpair_fastq_file_suffix
        
        os.system('bowtie2 -x %s/%s -q -U %s --no-head --no-unal -S %s --un %s -p %d'
            %(genome_host_path,genome_host_bt2_index_base_name,afterNTVIRUSfastqFiles[i],outputSAMFileGBHG,outputFASTQFileGBHG,num_threads))
        
        afterHGfastqFiles.append(outputFASTQFileGBHG)



    timer_bowtie2_hg38_end = time.clock()
    timer_bowtie2_hg38_end1 = time.time()
    bowtie2_hg38_time = timer_bowtie2_hg38_end - timer_bowtie2_hg38_start
    bowtie2_hg38_time1 = timer_bowtie2_hg38_end1 - timer_bowtie2_hg38_start1
    '''
    print '@@@@@----------timer:bowtie2_hg38_time_clock:'
    print bowtie2_hg38_time    #17.103814
    print '@@@@@----------timer:bowtie2_hg38_time_time:'
    print bowtie2_hg38_time1    #1568.48754096
    '''
    sheetPipelinerStaClock.write(1,2,bowtie2_hg38_time)
    sheetPipelinerStaTime.write(1,2,bowtie2_hg38_time1)

    print '==========assembly : velvet with alignment data above'

    timer_velvet_start = time.clock()
    timer_velvet_start1 = time.time()

    minContigLength = 100
    for i in range(len(sample_se_list)):
        for kmer in range(kmerMin,kmerMax,kmerStep):
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_se_list[i],kmer)
            os.system('velveth %s %d -short -fastq %s' %(eachKmerDir,kmer,afterHGfastqFiles[i]))
            os.system('velvetg %s -exp_cov auto -cov_cutoff auto -min_contig_lgth %d' %(eachKmerDir,minContigLength))

            eachKmerFiles = os.listdir(eachKmerDir)
            eachKmerFilesAP = []
            for j in range(len(eachKmerFiles)):
                eachKmerFilesAP.append(eachKmerDir+'/'+eachKmerFiles[j])
            for j in range(len(eachKmerFilesAP)):
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Graph2' : os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/LastGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/PreGraph': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Roadmaps': os.remove(eachKmerFilesAP[j])
                if os.path.exists(eachKmerFilesAP[j]) and eachKmerFilesAP[j] == eachKmerDir+'/Sequences': os.remove(eachKmerFilesAP[j])

    timer_velvet_end = time.clock()
    timer_velvet_end1 = time.time()
    velvet_time = timer_velvet_end - timer_velvet_start
    velvet_time1 = timer_velvet_end1 - timer_velvet_start1
    '''
    print '@@@@@----------timer:velvet_time_clock:'
    print velvet_time    #0.034675
    print '@@@@@----------timer:velvet_time_time:'
    print velvet_time1    #1.79084706306
    '''
    sheetPipelinerStaClock.write(1,3,velvet_time)
    sheetPipelinerStaTime.write(1,3,velvet_time1)

    ##================statistics velvet : pick up the proper kmer , accoding to N50 max========
    velvetSta = xlwt.Workbook()
    sheetVelvetSta = velvetSta.add_sheet('velvet_statistics', cell_overwrite_ok=True)
    kmerNum = (kmerMax-kmerMin)/2  #(35-31)/2=2
    row0 = ['sampleNums','kmer','medianCovDepth','nodesNum','contigN50','contigMax','genomeLen','usedReads','totalReads']
    for i in range(0,len(row0)):
        sheetVelvetSta.write(0,i,row0[i])  #first row

    #first cow is sample_se_list
    #second cow is kmer (kmerMin to kmerMax , kmerStep)

    for i in range(0,len(sample_se_list)):
        sheetVelvetSta.write(kmerNum*i+1,0,sample_se_list[i])   #first column is the sampleNums

    colKmer = range(kmerMin,kmerMax,kmerStep)
    for i in range(0,len(sample_se_list)*kmerNum,kmerNum): #i : every column
        for k in range(0,kmerNum):
            sheetVelvetSta.write(i+k+1,1,colKmer[k])   #second column is the k-mer
            k += 1

    logLast2Str = []
    MedianCovDepth = []

    logLast1Str = []
    nodesNum = []
    contigN50 = []
    contigMax = []
    genomeLen = []
    usedReads = []
    totalReads = []
    for i in range(len(sample_se_list)):
        for kmer in range(kmerMin,kmerMax,kmerStep):
            eachKmerDir = '%s%s_%d'%(velvet_path,sample_se_list[i],kmer)
            Log_file = os.path.join(eachKmerDir,'Log')
            logs = open(Log_file,'r').readlines()
            #logLast2line = logs[len(logs)-2]
            #if logLast2line != '':
            #Median coverage depth = 96.434783  #not the logs[26] , to avoid the task again with adding the log_information to Log File
            #logLast2Str = logLast2line.split(' ')
            #print logLast2Str
            #print logLast2Str[4][:-2]
            #MedianCovDepth.append(string.atof(logLast2Str[4][:-2]))
            logLast1line = logs[len(logs)-1]  #Final graph has 199 nodes and n50 of 47325, max 166302, total 1570595, using 1470371/1487666 reads
            logLast1Str = logLast1line.split(' ')
            #['Final', 'graph', 'has', '172', 'nodes', 'and', 'n50', 'of', '46917,', 'max', '237233,', 'total', '1570012,', 'using', '1478778/1487666', 'reads\n']
            #print logLast1Str
            nodesNum.append(string.atoi(logLast1Str[3], 10))
            contigN50.append(string.atoi(logLast1Str[8][:-1],10))
            contigMax.append(string.atoi(logLast1Str[10][:-1],10))
            genomeLen.append(string.atoi(logLast1Str[12][:-1],10))
            usedReads.append(string.atoi(logLast1Str[14].split('/')[0], 10))
            totalReads.append(string.atoi(logLast1Str[14].split('/')[1], 10))
    #print contigN50
    for i in range(0,len(sample_se_list)*kmerNum):
        #sheetVelvetSta.write(i+1,2,MedianCovDepth[i])
        sheetVelvetSta.write(i+1,3,nodesNum[i])
        sheetVelvetSta.write(i+1,4,contigN50[i])
        sheetVelvetSta.write(i+1,5,contigMax[i])
        sheetVelvetSta.write(i+1,6,genomeLen[i])
        sheetVelvetSta.write(i+1,7,usedReads[i])
        sheetVelvetSta.write(i+1,8,totalReads[i])

    velvetSta.save(os.path.join(result_total_sta_path,'velvet_log_sta.xls'))

    ##check the max contigN50 to pick up ,then blastn
    contigN50_max_value_list = []  #all_samples
    contigN50_max_index_list = []
    for i in range(len(sample_se_list)):
        each_sample_contigN50_list = contigN50[kmerNum*i:kmerNum*(i+1)]
        #print each_sample_contigN50_list
        contigN50_max_value_list.append(max(each_sample_contigN50_list))
        contigN50_max_index_list.append(kmerMin+(each_sample_contigN50_list.index(max(each_sample_contigN50_list)))*kmerStep)

    #print contigN50_max_value_list
    #print contigN50_max_index_list

    print '=============blast:use the contigs.fa to map AGINST NT_viruses'
    #1.first:format the local nt fasta_database file
    #os.system('makeblastdb -dbtype nucl -input_type fasta -in %s' %(nt_path))  ###big memory : server

    #after makeblastdb,generate some other indexed files besides nt file,
    #use the indexed nt file(with some other indexed files) as the following input ref db '-db'.

    #2.second:blastn

    timer_blastn_start = time.clock()
    timer_blastn_start1 = time.time()

                 
    for i in range(len(sample_se_list)):
        kmer = contigN50_max_index_list[i]
        eachKmerDir = '%s%s_%d'%(velvet_path,sample_se_list[i],kmer) #/public1/home/yefq/project_adenovirus/results/velvet_result/jiayan_31
        #print eachKmerDir
        eachContigFile = os.path.join(eachKmerDir,'contigs.fa')
        #print eachContigFile
        if os.path.getsize(eachContigFile) != 0:
            eachBlastResult = os.path.join(blast_path,sample_se_list[i]+'.out')
            os.system('blastn -db %s -query %s -outfmt 6 -dust no -num_threads %d -perc_identity 80 -word_size 20 -max_target_seqs %d -evalue 0.0000001 -out %s'
                %(nt_viruses_family_file,eachContigFile,num_threads,max_target_seqs_nt_viruses,eachBlastResult))

    timer_blastn_end = time.clock()
    timer_blastn_end1 = time.time()
    blastn_time = timer_blastn_end - timer_blastn_start
    blastn_time1 = timer_blastn_end1 - timer_blastn_start1
    '''
    print '@@@@@----------timer:blastn_time_clock:'
    print blastn_time  #0.287237
    print '@@@@@----------timer:blastn_time_time:'
    print blastn_time1  #30.0745470524
    '''
    sheetPipelinerStaClock.write(1,4,blastn_time)
    sheetPipelinerStaTime.write(1,4,blastn_time1)

    ##3.third:sorted the blast result, use the different colomn
    for i in range(len(colNumSorted_list)):
        sort_blast_result_by_blastnColumn(colNumSorted_list[i])


    time_total_end = time.clock()
    time_total_end1 = time.time()
    

    total_time = time_total_end - time_total_start
    total_time1 = time_total_end1 - time_total_start1
    '''
    print '@@@@@-----------timer:total_time_clock:'
    print total_time  #23.519992
    print '@@@@@-----------timer:total_time_time:'
    print total_time1  #2177.6757071
    '''
    sheetPipelinerStaClock.write(1,5,total_time)
    sheetPipelinerStaTime.write(1,5,total_time1)

    pipelinerSta.save(os.path.join(result_total_sta_path,'pipeline_time_sta.xls'))

#-*- coding:utf-8 -*-

from .config_paras import *

if __name__ == '__main__':
    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    os.environ['DJANGO_SETTINGS_MODULE'] = 'VIOS.settings'
    django.setup()



##index the nt file : 20170614(download)
def nt_database_pre_process(ntFile):
	os.system('makeblastdb -dbtype nucl -input_type fasta -in %s' %(ntFile))
	#ntVirusesLen1000File_dir = os.path.dirname(ntVirusesLen1000File)
	#os.system('bowtie2-build -f -q --threads %d %s %s/%s' %(numThreads,ntVirusesLen1000File,ntVirusesLen1000File_dir,bt2_index_base))  #%s %s%s : ref bt2_index_base
	os.system('samtools faidx %s'%(ntFile))


#nt_database_pre_process(nt_file)


##index the nt_index_title
#def nt_index_title() #the nt file is too large, so in big_server to get the log/ouput by using "print .."


def nt_viruses_onceforall(numThreads,ntVirusesFamilyNameList):
	"""index the nt_viruses _family and generate index_header file"""
	nt_viruses_family_file_list,nt_viruses_family_path_list,bt2_index_base_name_list,nt_viruses_family_index_header_file_list,nt_viruses_family_seqLengthList_file_list = get_nt_viruses_family_results(ntVirusesFamilyNameList)
	#print nt_viruses_family_file_list
	for i in range(len(ntVirusesFamilyNameList)):
		print '===========%s'%(ntVirusesFamilyNameList[i])
		#format the local nt_ file to blastn , maxFileSize = '2GB'  #BLAST options error: max_file_sz must be < 2 GiB
		#os.system('makeblastdb -dbtype nucl -input_type fasta -max_file_sz 2GB -in %s' %(nt_viruses_family_file_list[i]))
		'''
		os.system('makeblastdb -dbtype nucl -input_type fasta -in %s' %(nt_viruses_family_file_list[i]))
		os.system('bowtie2-build -f -q --threads %d %s %s/%s' %(numThreads,nt_viruses_family_file_list[i],nt_viruses_family_path_list[i],bt2_index_base_name_list[i]))
		os.system('samtools faidx %s'%(nt_viruses_family_file_list[i]))
		'''

		##generate nt_index_header file
		
		nt_viruses_family_file_obj = open(nt_viruses_family_file_list[i], 'r')
		nt_viruses_family_index_header_file_obj = open(nt_viruses_family_index_header_file_list[i],'w+')
		try:
			while True:
			    eachLine = nt_viruses_family_file_obj.readline()
			    if eachLine == '':  #read the end of file : EOF
			        break;
			    if eachLine[0] == '>':
			        #header = eachLine.split('>')[-1]
			        # >V01351.1 Sea urchin fragment, 3' to the actin gene in <SPAC01>
			        header = eachLine[1:-1]
			        assession = header.split(' ')[0]
			        annotation = header[header.index(' ') + 1:]
			        for i in range(len(annotation) - 1):
			            if not (re.compile(r'\w|,|\'|<|>|\.|-|_|\(|\)|\+|:')).match(annotation[i]):
			                annotation = annotation.replace(annotation[i], ' ')
			        nt_viruses_family_index_header_file_obj.write(assession+'\t'+annotation+'\n')
		finally:
			nt_viruses_family_file_obj.close()
			nt_viruses_family_index_header_file_obj.close()
		
		##generage nt_viruses_family_seqLengthList_file
		with open(nt_viruses_family_file_list[i], 'r') as nt_viruses_family_file_obj:
			lines_list = nt_viruses_family_file_obj.readlines()
			lines_list_len = len(lines_list)

			keys_header_list = []
			keys_header_index_list = []  #'>' at line number
			keys_accession_list = []
			for j in range(len(lines_list)):
			    if lines_list[j][0] == '>':
			    	each_header = lines_list[j].strip('\n')
			        keys_header_list.append(each_header)
			        keys_header_index_list.append(j)
			        accession = each_header.split(' ')[0].split('>')[-1]
			        keys_accession_list.append(accession)
			keys_header_index_list.append(lines_list_len)
			#print keys_header_list
			#print keys_header_index_list
			#print keys_accession_list


			values_seq_list = []  #2D:each elemnt is values_seq_one_list
			for j in range(1, len(keys_header_index_list)):
			    values_seq_one_list = []
			    for k in range(keys_header_index_list[j - 1] + 1, keys_header_index_list[j]):
			        values_seq_one = lines_list[k].strip('\n')
			        values_seq_one_list.append(values_seq_one)                
			    values_seq_list.append(values_seq_one_list)  # before values_seq_one_list = [], append
			#print values_seq_list

			values_seqNum_list = []  #2D
			for j in range(len(values_seq_list)):
				oneSeq_list = values_seq_list[j]
				oneSeq_eachLineLength_list = []
				for k in range(len(oneSeq_list)):
					oneSeq_eachLineLength = len(oneSeq_list[k])
					oneSeq_eachLineLength_list.append(oneSeq_eachLineLength)
				#print oneSeq_eachLineLength_list
				oneSeqNum = sum(oneSeq_eachLineLength_list)
				values_seqNum_list.append(oneSeqNum)
			#print values_seqNum_list

		with open(nt_viruses_family_seqLengthList_file_list[i],'w+') as nt_viruses_family_seqLengthList_file_obj:
			for j in range(len(keys_accession_list)):
				accession = keys_accession_list[j]
				seqNum = values_seqNum_list[j]
				nt_viruses_family_seqLengthList_file_obj.writelines(accession+'\t'+str(seqNum)+'\n')


			'''
			for i in range(len(keys_accession_list)):
				accession = keys_accession_list[i].split(' ')[0].split('>')[-1]
				if accession == queryTheAccession:
					accession_length_sum = []
					for j in range(len(values_seq_list[i])):
						one_line_length = len(values_seq_list[i][j])
						accession_length_sum.append(one_line_length)
					return sum(accession_length_sum)
			'''


#nt_viruses_family_name_list = ['viruses_final','viruses_full','bacterias_full','viruses_Adenoviridae','viruses_Coronaviridae','viruses_Filoviridae','viruses_Flaviviridae','viruses_Orthomyxoviridae']  #0 means viruses full database
max_target_seqs_nt_viruses_list = [120211, 3817, 6224600, 9894, 2462, 65183, 420380]  #while employ nt_viruses_onceforall , you can change it.
max_target_seqs_nt_viruses = max(max_target_seqs_nt_viruses_list)  #420380, use the max as the paras of 'blastn -max_target_seqs'


def genome_hosts_onceforall(numThreads,genomeHostsNameList):
	"""index the genome_hosts"""
	genome_hosts_file_list,genome_hosts_path_list,bt2_index_base_name_list = get_genome_hosts_results(genomeHostsNameList)
	for i in range(len(genomeHostsNameList)):
		print '===========%s'%(genomeHostsNameList[i])
		genome_hosts_fasta_file = os.path.join(genome_hosts_path_list[i],genomeHostsNameList[i]+'.fa')
		#print genome_hosts_fasta_file
		if os.path.getsize(genome_hosts_file_list[i]):  #some genome not in USCS with 2bit format, but can find in Genome browser and ftp download fasta format
			os.system('twoBitToFa %s %s' %(genome_hosts_file_list[i],genome_hosts_fasta_file))
		os.system('bowtie2-build -f -q --threads %d %s %s/%s'%(numThreads,genome_hosts_fasta_file,genome_hosts_path_list[i],bt2_index_base_name_list[i]))
		

def genome_host_adder(numThreads,genomeHostName,genomeHostNameReadable):
	path_name = genomeHostName
	file_name = '%s.2bit'%(genomeHostName)
	genome_host_path = genomesDB_path+path_name
	if not os.path.exists(genome_host_path):
		os.system('mkdir %s'%(genome_host_path))
	genome_host_file = os.path.join(genome_host_path,file_name)
	bt2_index_base_name = 'bt2_index_%s'%(genomeHostName)

	print '===========%s'%(genomeHostName)
	genome_host_adder_file = os.path.join(db_path,'genome_host_adder.log')
	genome_hosts_name_adder_list = []
	if os.path.exists(genome_host_adder_file):#if file exists,check. if no exists just pass
		with open(genome_host_adder_file,'r') as genome_host_adder_file_obj:
			genome_hosts_adder_list = genome_host_adder_file_obj.readlines()
			for i in range(len(genome_hosts_adder_list)):
				each_genome_host = genome_hosts_adder_list[i].strip('\n')
				each_genome_host_name = each_genome_host.split('\t')[0]
				genome_hosts_name_adder_list.append(each_genome_host_name)

	###no repeate submit the genome name by users
	genome_hosts_name_updated_list = genome_hosts_name_list+genome_hosts_name_adder_list
	if not genomeHostName in genome_hosts_name_updated_list:  
		os.system('rsync -aP rsync://hgdownload.soe.ucsc.edu/goldenPath/%s/bigZips/%s.2bit %s' %(genomeHostName,genomeHostName,genome_host_path))
		genome_host_fasta_file = os.path.join(genome_host_path,genomeHostName+'.fa')
		
		##rsync download fail, there is no user's host name in USCS
		if not os.path.exists(genome_host_file):
			os.rmdir(genome_host_path)
		if os.path.exists(genome_host_path) and os.path.getsize(genome_host_file):#file exists and not empty
			os.system('twoBitToFa %s %s' %(genome_host_file,genome_host_fasta_file))
			os.system('bowtie2-build -f -q --threads %d %s %s/%s'%(numThreads,genome_host_fasta_file,genome_host_path,bt2_index_base_name))
			#if all is successful,remember the genome_host_name_readable
			genome_host_adder_file = os.path.join(db_path,'genome_host_adder.log')
			with open(genome_host_adder_file,'a+') as genome_host_adder_file_obj:
				genome_hosts_adder_list = [genomeHostName+'\t',genomeHostNameReadable+'\n']
				genome_host_adder_file_obj.writelines(genome_hosts_adder_list)


		

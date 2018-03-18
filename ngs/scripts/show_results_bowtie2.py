#-*- coding:utf-8 -*-

from .config_paras import *
from .db_nt_process import *

from matplotlib.pyplot import *

if __name__ == '__main__':
    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    os.environ['DJANGO_SETTINGS_MODULE'] = 'VIOS.settings'
    django.setup()


def nt_filter_from_accession(fromBigFile,queryTheAccession):
    '''
    get the fasta file, from reference(fromBigFile) to filter the target accession's fasta.
    then index the refence file, but if there is accession path ,this method can pass(maybe not achived fully, avoid this)
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
 
        file_dict = dict(zip(keys_accession_list, values_seq_list))


    toFilterFilePath = os.path.join(nt_viruses_the_accession_path,queryTheAccession)
    if not os.path.exists(toFilterFilePath):
        os.mkdir(toFilterFilePath)
    toFilterFile = os.path.join(toFilterFilePath,queryTheAccession+'.fasta')
    with open(toFilterFile, 'w+') as toFilterFile_obj:
        for j in range(len(keys_accession_list)):
            accession = keys_accession_list[j].split(' ')[0].split('>')[-1]
            if accession == queryTheAccession:
                target_seq_one_list = file_dict.get(keys_accession_list[j])  # if KeyError ,return None
                toFilterFile_obj.write(keys_accession_list[j] + '\n')
                for k in range(len(target_seq_one_list)):
                    target_seq_one_line = target_seq_one_list[k] + '\n'
                    toFilterFile_obj.write(target_seq_one_line)

def mpileup_to_coverage_grafic(mpileupFile,targetAccession,querySample):
    with open(mpileupFile, 'r') as mpileup_file_obj:
        lines_list = mpileup_file_obj.readlines()
        pos_list = []
        cov_list = []
        for i in range(len(lines_list)):
            each_line = lines_list[i]
            pos = int(each_line.split('\t')[1])
            cov = int(each_line.split('\t')[3])
            if cov == 0:
                cov_log = 0
            else:
                cov_log = math.log(cov,10)#log(cov,10)  #base is 10
            pos_list.append(pos)
            #cov_list.append(cov)
            cov_list.append(cov_log)
        min_pos = pos_list[0]
        max_pos = pos_list[len(pos_list)-1]

    #fig = figure(figsize=(10,5),dpi=80)
    fig,ax = subplots(figsize=(10,5))
    x = pos_list
    y = cov_list

    plot(x,y)
    fill_between(x,0,y)

    title_str = "%s base coverage with mappting aginst %s"%(querySample,targetAccession)
    title(title_str,fontsize=14)
    xlabel("base position")
    ylabel("coverage")
    xlim(min_pos,max_pos)
    ylim(ymin=0)
    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    grid(color='lightskyblue')

    fig_file = os.path.join(os.path.dirname(mpileupFile),targetAccession+'.png')
    fig.savefig(fig_file)
    show()


def show_bowtie2_base_coverage(ntVirusesFamilyName,targetAccession,querySample):
    '''0.nt_filter_from_accession method :
    from reference to filter the target accession's fasta,and index the fasta to database directory 
    pick up the targetAccession's fasta file as reference , use all filter aginst GH NT database output files to map lonely.
    '''
    family_path = 'nt_%s'%(ntVirusesFamilyName)
    family_file_name = 'nt_%s.fasta'%(ntVirusesFamilyName)
    
    nt_viruses_family_path = blastDB_path+family_path
    nt_viruses_family_file = os.path.join(nt_viruses_family_path,family_file_name)
    
    to_filter_path = os.path.join(nt_viruses_the_accession_path,targetAccession)
    #print to_filter_path
    to_filter_file = os.path.join(to_filter_path,targetAccession+'.fasta')
    bt2_index_base_name = 'bt2_index_%s'%(targetAccession)
    nt_filter_from_accession(nt_viruses_family_file,targetAccession)  #...blastDB/nt_viruses_the_accession/KJ883521.1.fasta

    os.system('bowtie2-build -f -q --threads %d %s %s/%s' %(num_threads,to_filter_file,to_filter_path,bt2_index_base_name))  #%s %s%s : ref bt2_index_base
    os.system('samtools faidx %s'%(to_filter_file))


    ##use all filtered(after qc/virus steps output files) to map the targetAccession as the one reference
    outputFastqPair1FileNTVirus = bowtie2_path+querySample+'_pair_VIRUS.1.fastq'        
    outputFastqPair2FileNTVirus = bowtie2_path+querySample+'_pair_VIRUS.2.fastq'
    outputFastqUnpairFileNTVirus = bowtie2_path+querySample+map_nt_viruses_singleAndunpair_fastq_file_suffix
  
    output_path = bowtie2_path_basecoverage + querySample
    if not os.path.exists(output_path):
        os.system('mkdir %s'%output_path)
    '''else:
        for each_out_put_file in os.listdir(output_path):
            os.remove(os.path.join(output_path,each_out_put_file))'''
    
    output_file_base_name = querySample+'_'+targetAccession
    SAMFileAccession = os.path.join(output_path,output_file_base_name+'.sam')	    
    if os.path.exists(outputFastqPair1FileNTVirus) and os.path.exists(outputFastqPair2FileNTVirus): #if sequenceMethod == 'pe':
        os.system('bowtie2 -x %s/%s -q -1 %s -2 %s -U %s -S %s -p %d'%(to_filter_path,bt2_index_base_name,outputFastqPair1FileNTVirus,outputFastqPair2FileNTVirus,outputFastqUnpairFileNTVirus,SAMFileAccession,num_threads))
    else:
        os.system('bowtie2 -x %s/%s -q -U %s -S %s -p %d'%(to_filter_path,bt2_index_base_name,outputFastqUnpairFileNTVirus,SAMFileAccession,num_threads))
        
    BAMFileAccession = os.path.join(output_path,output_file_base_name+'.bam')
    BAMSortedFileAccession = os.path.join(output_path,output_file_base_name+'.sorted.bam')
    MPILEUPFileAccession = os.path.join(output_path,output_file_base_name+'.mpileup')

    os.system('samtools view -bS %s > %s' % (SAMFileAccession,BAMFileAccession))
    '''
    os.system('samtools sort -m 2G %s -o %s' % (BAMFileAccession,BAMSortedFileAccession))
    os.system('samtools index %s' % (BAMSortedFileAccession))
    os.system('samtools mpileup -d 10000000000 -f %s %s > %s' % (to_filter_file,BAMSortedFileAccession,MPILEUPFileAccession))
    '''
    ##the coverage of each base:
    mpileup_to_coverage_grafic(MPILEUPFileAccession,targetAccession,querySample)


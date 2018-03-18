#from __future__ import absolute_import
from django.shortcuts import render,redirect,render_to_response
from django.http import JsonResponse
from django.views import View
from django.contrib.auth.decorators import login_required
from django.utils.decorators import method_decorator

from .forms import SequencingFilesForm,LoginUsernameForm
from .models import SequencingFiles,LoginUsername

from .scripts.db_nt_process import *

from .scripts.pipeliner import *
from .scripts.show_results_bowtie2 import *
from .scripts.show_results_blastn import *


from .scripts.show_results_total_time import *
from .scripts.show_results_total_reads import *

def mkdir_path():
    if os.path.exists(rawData_path):
        rawDataFiles = os.listdir(rawData_path)
    if not os.path.exists(MEDIA_ROOT):
        os.system('mkdir %s'%MEDIA_ROOT)
        os.system('mkdir %s'%(MEDIA_ROOT+"/"+time.strftime("%Y")))
        os.system('mkdir %s'%(MEDIA_ROOT+"/"+time.strftime("%Y")+"/"+time.strftime("%m")))
    if not os.path.exists(nt_viruses_the_accession_path):
        os.system('mkdir %s'%nt_viruses_the_accession_path)
    if not os.path.exists(result_path):
        os.system('mkdir %s'%result_path)
    if not os.path.exists(result_total_sta_path):
        os.system('mkdir %s'%result_total_sta_path)
    if not os.path.exists(bowtie2_path):
        os.system('mkdir %s'%bowtie2_path)
    if not os.path.exists(bowtie2_path_basecoverage):
        os.system('mkdir %s'%bowtie2_path_basecoverage)
    if not os.path.exists(velvet_path):
        os.system('mkdir %s'%velvet_path)
    if not os.path.exists(blast_path):
        os.system('mkdir %s'%blast_path)
    for i in range(len(blast_path_sorted_list)):
        if not os.path.exists(blast_path_sorted_list[i]):
            os.system('mkdir %s'%blast_path_sorted_list[i])
    if not os.path.exists(blast_contigs_againt_ref_path):
        os.system('mkdir %s'%blast_contigs_againt_ref_path)

##start each query, execute this method , make sure path must create
mkdir_path()


##@method_decorator(login_required,name='dispatch')
class ProgressBarUploadView(View):
    def get(self, request):
        SequencingFiles_list = SequencingFiles.objects.all()
        return render(self.request, 'ngs/progress_bar_upload.html', {'SequencingFiles': SequencingFiles_list})

    def post(self, request):
        #time.sleep(3)  # You don't need this line. This is just to delay the process so you can see the progress bar testing locally.
        form = SequencingFilesForm(self.request.POST, self.request.FILES)
        if form.is_valid():
            SequencingFiles = form.save()
            #print SequencingFiles
            #print SequencingFiles.username

            '''#write to database
            db = SequencingFiles()
            db.username = SequencingFiles.username
            db.SeqFiles = SequencingFiles.SeqFiles
            db.uploaded_at = SequencingFiles.uploaded_at
            db.save()'''

            data = {'is_valid': True, 'user_name':SequencingFiles.username,'name': SequencingFiles.SeqFiles.name, 'url': SequencingFiles.SeqFiles.url}
            #print data  #{'url': '/media/2017/05/04/hello', 'is_valid': True, 'name': u'2017/05/04/hello'}
            #print JsonResponse(data)
        else:
            data = {'is_valid': False}
        return JsonResponse(data)  #An HTTP response class that consumes data to be serialized to JSON.


def addLoginUsername(request):
	if request.method == "POST":
		loginusername_form = LoginUsernameForm(request.POST,request.FILES)
		if loginusername_form.is_valid():
			login_username = loginusername_form.cleaned_data['username']

			#write to database
			loginusername_db = LoginUsername()
			loginusername_db.username = login_username
			loginusername_db.save()

			#if all this is successful(get form data and write data to the database)
			return render_to_response('ngs/progress_bar_upload.html',{'login_username':login_username})
	else:
	    loginusername_form = LoginUsernameForm()
	return render(request,'ngs/addloginusername.html',{'loginusername_form':loginusername_form})


def clear_database(request):
    for sequencing_files in SequencingFiles.objects.all():
        sequencing_files.SeqFiles.delete()
        sequencing_files.delete()
    return redirect(request.POST.get('next'))



####database
class DatabaseView(View):
    def get(self,request):
        return render(self.request,'ngs/database.html')

def database_nt_viruses_handler(request):
    if request.method == 'POST':
        print '===========database nt viruses handler=============='
        nt_viruses_onceforall(num_threads,nt_viruses_family_name_list)               
    return redirect(request.POST.get('nt_viruses'))

def database_genome_hosts_handler(request):
    if request.method == 'POST':
        print '===========database genome hosts handler============='
        genome_hosts_onceforall(num_threads,genome_hosts_name_list)        
    return redirect(request.POST.get('genome_hosts'))
def database_genome_host_adder(request):
    if request.method == 'POST':
        print '===========database genome host adder==============='
        genome_host_name = request.POST.get('genome_host_name')
        genome_host_name_readable = request.POST.get('genome_host_name_readable')
        
        genome_host_adder(num_threads,genome_host_name,genome_host_name_readable)        
    return redirect(request.POST.get('genome_host_adder'))

####process
class ProcessView(View):
    def get(self,request):
        return render(self.request,'ngs/process.html',
            {'genome_hosts_name_readable_updated_list':genome_hosts_name_readable_updated_list,
            'nt_viruses_family_name_readable_list':nt_viruses_family_name_readable_list})

#database:NO HOST
def process_pipeline_NoHostGeome(request):
    if request.method == 'POST':        
        print '============process_pipeline============'
        genome_host_name_readable = request.POST.get('genome_hosts')
        nt_virus_family_name_readable = request.POST.get('nt_virus_family')  ##Adenoviridae   value="Adenoviridae"(user choose)
        sequence_method = request.POST.get('pipeliner')
        #print nt_virus_family_name_readable

        reference_list = []
        '''
        for i in range(len(genome_hosts_name_readable_updated_list)):
            if genome_host_name_readable == genome_hosts_name_readable_updated_list[i]:
                genome_host_name = genome_hosts_name_updated_list[i]
                reference_list.append(genome_host_name)
        '''
        for i in range(len(nt_viruses_family_name_readable_list)):
            if nt_virus_family_name_readable == nt_viruses_family_name_readable_list[i]:
                nt_virus_family_name = nt_viruses_family_name_list[i]
                reference_list.append(nt_virus_family_name)
        #print reference_list

        if sequence_method == 'pe':
            process_pipeliner_pe(reference_list[0])
        if sequence_method == 'se':
            process_pipeliner_se(reference_list[0])

        results_paras_transfer_log_file = os.path.join(result_total_sta_path,'paras_transfer_log.txt')
        with open(results_paras_transfer_log_file,'w+') as results_paras_transfer_log_file_obj:
            paras_list = [nt_virus_family_name+'\n',sequence_method+'\n']
            results_paras_transfer_log_file_obj.writelines(paras_list)

    return redirect(request.POST.get('process_pipeliner'))

#database:with Host Genome
def process_pipeline(request):
    if request.method == 'POST':        
        print '============process_pipeline============'
        genome_host_name_readable = request.POST.get('genome_hosts')
        nt_virus_family_name_readable = request.POST.get('nt_virus_family')  ##Adenoviridae   value="Adenoviridae"(user choose)
        sequence_method = request.POST.get('pipeliner')
        #print nt_virus_family_name_readable

        reference_list = []
        
        for i in range(len(genome_hosts_name_readable_updated_list)):
            if genome_host_name_readable == genome_hosts_name_readable_updated_list[i]:
                genome_host_name = genome_hosts_name_updated_list[i]
                reference_list.append(genome_host_name)
        
        for i in range(len(nt_viruses_family_name_readable_list)):
            if nt_virus_family_name_readable == nt_viruses_family_name_readable_list[i]:
                nt_virus_family_name = nt_viruses_family_name_list[i]
                reference_list.append(nt_virus_family_name)
        #print reference_list

        if sequence_method == 'pe':
            process_pipeliner_pe(reference_list[0],reference_list[1])
        if sequence_method == 'se':
            process_pipeliner_se(reference_list[0],reference_list[1])

        results_paras_transfer_log_file = os.path.join(result_total_sta_path,'paras_transfer_log.txt')
        with open(results_paras_transfer_log_file,'w+') as results_paras_transfer_log_file_obj:
            paras_list = [nt_virus_family_name+'\n',sequence_method+'\n']
            results_paras_transfer_log_file_obj.writelines(paras_list)

    return redirect(request.POST.get('process_pipeliner'))


####results
class ResultsView(View):
    def get(self,request):
        results_paras_transfer_log_file = os.path.join(result_total_sta_path,'paras_transfer_log.txt')
        if os.path.exists(results_paras_transfer_log_file):            
            with open(results_paras_transfer_log_file,'r') as results_paras_transfer_log_file_obj:
                paras_list = results_paras_transfer_log_file_obj.readlines()
                #print paras_list
                sequence_method = paras_list[1].strip('\n')
                sample_list_results = []
                if sequence_method == 'pe':
                    sample_list_results = sample_list
                if sequence_method == 'se':
                    sample_list_results = sample_se_list
                #print sample_list_results
                #print paras_list[0].strip('\n')
        else:
            sample_list_results = []
        return render(self.request,'ngs/results.html',{'sample_list_results':sample_list_results,})


def show_results_total_time(request):
    if request.method == 'POST':
        print '=============total time statics==============='
        results_paras_transfer_log_file = os.path.join(result_total_sta_path,'paras_transfer_log.txt')
        with open(results_paras_transfer_log_file,'r') as results_paras_transfer_log_file_obj:
            paras_list = results_paras_transfer_log_file_obj.readlines()
            sequence_method = paras_list[1].strip('\n')
            sample_list_results = []
            if sequence_method == 'pe':
                sample_list_results = sample_list
            if sequence_method == 'se':
                sample_list_results = sample_se_list
            #print sample_list_results

        result_total_time_file = os.path.join(result_total_sta_path,'pipeline_time_sta.xls')
        
        pipeline_time_sta_to_grafic(sample_list_results,result_total_time_file)
    return render(request,'ngs/results.html',{'sample_list_results':sample_list_results,})


def show_results_blastn(request):
    if request.method == 'POST':
        search_time_clock_start = time.clock()
        search_time_time_start = time.time()

        print '=============blastn_results================='
        results_paras_transfer_log_file = os.path.join(result_total_sta_path,'paras_transfer_log.txt')
        with open(results_paras_transfer_log_file,'r') as results_paras_transfer_log_file_obj:
            paras_list = results_paras_transfer_log_file_obj.readlines()
            sequence_method = paras_list[1].strip('\n')
            sample_list_results = []
            if sequence_method == 'pe':
                sample_list_results = sample_list
            if sequence_method == 'se':
                sample_list_results = sample_se_list

            nt_virus_family = paras_list[0].strip('\n')

        query_sample = request.POST.get('query_sample')

        blastn_column = int(request.POST.get('blast_result_order'))

        target_top_lines = show_blastn_top_lines(nt_virus_family,query_sample,blastn_column)
        #queryEachLineFilterList = [ref,targetAnnotation,identity,alignment_length,pos_query_start,pos_query_end,pos_ref_start,pos_ref_end,e_value,bit_score]
        table_header = ['accession','annotation','identity','alignment_length','pos_query_start','pos_query_end','pos_ref_start','pos_ref_end','e_value','bit_score']
        #print target_top_lines
                   
        ##write list locally.        
        sample_contig_path = blast_contigs_againt_ref_path + query_sample
         
        if not os.path.exists(sample_contig_path):
            os.system('mkdir %s'%sample_contig_path)

        #['identity','length','score'] #[2,3,11]
        if blastn_column == 11:
            sortedBy = 'score'
        if blastn_column == 3:
            sortedBy ='length'
        if blastn_column == 2:
            sortedBy = 'identity'
        sample_contig_file =  os.path.join(sample_contig_path,'sta_contig_%s_%s.txt'%(sortedBy,query_sample))
        with open(sample_contig_file,'w+') as sample_contig_file_obj:
            for i in range(len(table_header)):
                if i != len(table_header)-1:
                    sample_contig_file_obj.writelines(table_header[i]+'\t')
                else:
                    sample_contig_file_obj.writelines(table_header[i]+'\n')

            for i in range(len(target_top_lines)):
                for j in range(len(target_top_lines[i])):
                    if j != len(target_top_lines[i])-1:
                        sample_contig_file_obj.writelines(str(target_top_lines[i][j])+'\t')
                    else:
                        sample_contig_file_obj.writelines(str(target_top_lines[i][j])+'\n')


        search_time_clock_end = time.clock()
        search_time_time_end = time.time()
        time_search_clock = search_time_clock_end - search_time_clock_start
        time_search_time = search_time_time_end - search_time_time_start
        '''
        print '-------------time_search_clock:'
        print time_search_clock
        print '-------------time_search_time:'
        print time_search_time
        '''
        
    #return redirect(request.POST.get('top_lines'))
    return render(request,'ngs/results.html',
        {'sample_chosed':query_sample,'family_chosed':nt_virus_family,'column_chosed':blastn_column,
        'sample_list_results':sample_list_results,'table_header':table_header,'targetTopLines':target_top_lines})

def show_results_contigs_sta(request,sample,nt_virus_family,blastn_column,accession):
    print '=============contigs_usage_sta================='
    results_paras_transfer_log_file = os.path.join(result_total_sta_path,'paras_transfer_log.txt')
    with open(results_paras_transfer_log_file,'r') as results_paras_transfer_log_file_obj:
        paras_list = results_paras_transfer_log_file_obj.readlines()
        sequence_method = paras_list[1].strip('\n')
        sample_list_results = []
        if sequence_method == 'pe':
            sample_list_results = sample_list
        if sequence_method == 'se':
            sample_list_results = sample_se_list        

        #nt_virus_family = paras_list[0].strip('\n')    

    pipeline_contigs_sta_to_grafic(nt_virus_family,sample,int(blastn_column),accession)

    #return redirect(request.POST.get('top_lines'))   
    return render(request,'ngs/results.html',{'sample_list_results':sample_list_results,})




def show_results_total_reads(request):
    if request.method == 'POST':
        print '=============reads_usage_sta================'
        results_paras_transfer_log_file = os.path.join(result_total_sta_path,'paras_transfer_log.txt')
        with open(results_paras_transfer_log_file,'r') as results_paras_transfer_log_file_obj:
            paras_list = results_paras_transfer_log_file_obj.readlines()
            sequence_method = paras_list[1].strip('\n')
            sample_list_results = []
            if sequence_method == 'pe':
                sample_list_results = sample_list
            if sequence_method == 'se':
                sample_list_results = sample_se_list

            nt_virus_family = paras_list[0].strip('\n')

        query_sample = request.POST.get('query_sample')   

        
        '''results_fig_temp_path = result_path + 'results_fig_temp'
        if not os.path.exists(results_fig_temp_path):
            os.mkdir(results_fig_temp_path)'''

        pipeline_reads_sta_to_grafic(sequence_method,sample_list_results,query_sample,nt_virus_family)
        
    return render(request,'ngs/results.html',{'sample_list_results':sample_list_results,})

def show_results_identifer_fig(request):
    if request.method == 'POST':
        print '=============reads_usage_sta================'
        results_paras_transfer_log_file = os.path.join(result_total_sta_path,'paras_transfer_log.txt')
        with open(results_paras_transfer_log_file,'r') as results_paras_transfer_log_file_obj:
            paras_list = results_paras_transfer_log_file_obj.readlines()
            sequence_method = paras_list[1].strip('\n')
            sample_list_results = []
            if sequence_method == 'pe':
                sample_list_results = sample_list
            if sequence_method == 'se':
                sample_list_results = sample_se_list

            nt_virus_family = paras_list[0].strip('\n')

        query_sample = request.POST.get('query_sample')   

        pipeline_identifer_fig(query_sample,nt_virus_family)
        
    return render(request,'ngs/results.html',{'sample_list_results':sample_list_results,})

def show_results_identifer_list(request):
    if request.method == 'POST':
        print '=============identifer================'
        results_paras_transfer_log_file = os.path.join(result_total_sta_path,'paras_transfer_log.txt')
        with open(results_paras_transfer_log_file,'r') as results_paras_transfer_log_file_obj:
            paras_list = results_paras_transfer_log_file_obj.readlines()
            sequence_method = paras_list[1].strip('\n')
            sample_list_results = []
            if sequence_method == 'pe':
                sample_list_results = sample_list
            if sequence_method == 'se':
                sample_list_results = sample_se_list

            nt_virus_family = paras_list[0].strip('\n')

        query_sample = request.POST.get('query_sample')  

        identifer_table_header = ['accession','annotation','readsNumber','readsRatio']

        target_identifer_list = pipeline_identifer_list(query_sample,nt_virus_family)
        
        ##write list locally.        
        sample_identifer_path = bowtie2_path_basecoverage + query_sample
         
        if not os.path.exists(sample_identifer_path):
            os.system('mkdir %s'%sample_identifer_path)
        
        sample_identifer_file =  os.path.join(sample_identifer_path,'sta_identifer_%s.txt'%(query_sample))
        with open(sample_identifer_file,'w+') as sample_identifer_file_obj:
            sample_identifer_file_obj.writelines(identifer_table_header[0]+'\t'+identifer_table_header[1]+'\t'+identifer_table_header[2]+'\t'+identifer_table_header[3]+'\n')
            for i in range(len(target_identifer_list)):
                each_line = target_identifer_list[i][0]+'\t'+target_identifer_list[i][1]+'\t'+str(target_identifer_list[i][2])+'\t'+str(target_identifer_list[i][3])+'\n'
                sample_identifer_file_obj.writelines(each_line)
        
    return render(request,'ngs/results.html',
        {'sample_chosed':query_sample,'family_chosed':nt_virus_family,
        'sample_list_results':sample_list_results,
        'identifer_table_header':identifer_table_header,
        'targetIdentiferList':target_identifer_list,
        })


def show_results_coverage_SAMsta(request,sample,nt_virus_family,accession):
    print '=============mapAgainstVirus_SAM_STATISTICS====================='
    results_paras_transfer_log_file = os.path.join(result_total_sta_path,'paras_transfer_log.txt')
    with open(results_paras_transfer_log_file,'r') as results_paras_transfer_log_file_obj:
        paras_list = results_paras_transfer_log_file_obj.readlines()
        sequence_method = paras_list[1].strip('\n')
        sample_list_results = []
        if sequence_method == 'pe':
            sample_list_results = sample_list
        if sequence_method == 'se':
            sample_list_results = sample_se_list

        #nt_virus_family = paras_list[0].strip('\n')
           
    eachAccessionCoverageFig_fromSAMFile(sample,nt_virus_family,accession)  #choose the method of two circle not Counter . due to time consumption

    return render(request,'ngs/results.html',{'sample_list_results':sample_list_results,})


def show_results_coverage_sta(request,sample,nt_virus_family,accession):
    bt2_time_clock_start = time.clock()
    bt2_time_time_start = time.time()
    print '=============bowtie2_results====================='
    results_paras_transfer_log_file = os.path.join(result_total_sta_path,'paras_transfer_log.txt')
    with open(results_paras_transfer_log_file,'r') as results_paras_transfer_log_file_obj:
        paras_list = results_paras_transfer_log_file_obj.readlines()
        sequence_method = paras_list[1].strip('\n')
        sample_list_results = []
        if sequence_method == 'pe':
            sample_list_results = sample_list
        if sequence_method == 'se':
            sample_list_results = sample_se_list


        #nt_virus_family = paras_list[0].strip('\n')
           
    show_bowtie2_base_coverage(nt_virus_family,accession,sample)
    bt2_time_clock_end = time.clock()
    bt2_time_time_end = time.time()      
    bt2_time_clock = bt2_time_clock_end - bt2_time_clock_start
    bt2_time_time = bt2_time_time_end - bt2_time_time_start
    print 'time_clock:%s' %(bt2_time_clock)
    print 'time_time:%s' %(bt2_time_time)
    return render(request,'ngs/results.html',{'sample_list_results':sample_list_results,})



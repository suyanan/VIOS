from django.conf.urls import url

from .scripts.config_paras import *
from . import views
from django.views.generic import TemplateView  #A view that renders a template.

app_name = 'ngs'

urlpatterns = [
    url(r'^$', TemplateView.as_view(template_name='ngs/pipelines.html'), name='ngs_pipelines'),
    url(r'^loginusername/$', views.addLoginUsername, name='loginusername'),
    url(r'^clear/$', views.clear_database, name='clear_database'),
    url(r'^progress-bar-upload/$', views.ProgressBarUploadView.as_view(), name='progress_bar_upload'),
    
    
    url(r'^database/$', views.DatabaseView.as_view(), name='database'),
    url(r'^database_nt_viruses/$',views.database_nt_viruses_handler, name='database_nt_viruses'),
    url(r'^database_genome_hosts/$',views.database_genome_hosts_handler, name='database_genome_hosts'),
    url(r'^database_genome_host_adder/$',views.database_genome_host_adder, name='database_genome_host_adder'),

    
    url(r'^process/$', views.ProcessView.as_view(), name='process'),
    url(r'^process_pipeline/$', views.process_pipeline, name='process_pipeline'),
    
    
    url(r'^results/$', views.ResultsView.as_view(), name='results'),
        
    url(r'^results/total_time/$', views.show_results_total_time, name='show_results_total_time'),

    url(r'^results/blastn/$', views.show_results_blastn, name='show_results_blastn'),     
    url(r'^results/contigs_sta/(?P<sample>\w+)/(?P<nt_virus_family>\w+)/(?P<blastn_column>\d+)/(?P<accession>\w+\.\d)/$', views.show_results_contigs_sta, name='show_results_contigs'),

    url(r'^results/total_reads/$', views.show_results_total_reads, name='show_results_total_reads'),  
    url(r'^results/which_identifer_fig/$', views.show_results_identifer_fig, name='show_results_identifer_fig'),
    url(r'^results/which_identifer_list/$', views.show_results_identifer_list, name='show_results_identifer_list'),
    url(r'^results/coverage_sta/(?P<sample>\w+)/(?P<nt_virus_family>\w+)/(?P<accession>\w+\.\d)/$', views.show_results_coverage_sta, name='show_results_coverage'),
    url(r'^results/coverage_sam_sta/(?P<sample>\w+)/(?P<nt_virus_family>\w+)/(?P<accession>\w+\.\d)/$', views.show_results_coverage_SAMsta, name='show_results_coverageSAM'),
]

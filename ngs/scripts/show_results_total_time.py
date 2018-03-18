#-*- coding:utf-8 -*-

#from .config_paras import *

from .config_paras import *
from matplotlib.pyplot import *

if __name__ == '__main__':
    sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    os.environ['DJANGO_SETTINGS_MODULE'] = 'VIOS.settings'
    django.setup()


def pipeline_time_sta_to_grafic(sampleListResults,resultTotalTimeFile):
	resultTotalTimeFile_bk = xlrd.open_workbook(resultTotalTimeFile)
	time_sta_sheet = resultTotalTimeFile_bk.sheet_by_index(1)  #time.time()
	nrows = time_sta_sheet.nrows
	ncols = time_sta_sheet.ncols
	x_data = ['map_host','map_viruses','assembly','blastn','total']
	#x_data = ['map_viruses','map_host','assembly','blastn','total']
	#x_data = ['map_viruses','assembly','blastn','total']
	y_data = time_sta_sheet.row_values(1)  #second row
	x_data.pop()
	y_data.pop()
	#print x_data
	#print y_data


	fig = figure(figsize=(10,4),dpi=80)
	

	ind = numpy.arange(len(y_data))
	#ind = numpy.arange(6)
	width = 0.5

	bar(ind,y_data,width,color='skyblue')

	
	title('Histogram of process time of %d samples'%(len(sampleListResults)),fontsize=14)
	#xlabel("process")
	ylabel("time/s")
	
	xticks(ind,x_data)
	#yticks(0,max(y_data),100)
	#ylim(ymin=0,ymax=max(y_data))
	

	fig_file = os.path.join(result_total_sta_path,'sta_time.png')
	#print fig_file
	fig.savefig(fig_file)
	show()

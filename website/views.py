from django.shortcuts import render

# Create your views here.

def manual(request):
	manual_msg = 'READ ME'
	return render(request,'website/manual.html',{'manual_message':manual_msg})

def reference(request):
	reference_msg = 'References'
	return render(request,'website/reference.html',{'reference_message':reference_msg})
'''
def faq(request):
	faq_msg = 'FAQ'
	return render(request,'website/faq.html',{'faq_message':faq_msg})
'''

def contact(request):
	contact_msg = 'Contact us'
	return render(request,'website/contact.html',{'contact_message':contact_msg})

def vios_intro(request):
	return render(request,'website/includes/vios_intro.html')

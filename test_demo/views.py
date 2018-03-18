from django.shortcuts import render,get_object_or_404
from django.urls import reverse_lazy
from django.http import JsonResponse
from django.views.generic import ListView,DetailView
from django.views.generic.edit import FormView,CreateView, UpdateView, DeleteView

from .models import Publisher,Book,Author
from .forms import ContactForm

# Create your views here.


class PublisherList(ListView):
    model = Publisher
    context_object_name = 'my_favorite_publishers'


class PublisherDetail(DetailView):
    context_object_name = 'publisher'
    queryset = Publisher.objects.all()

    '''model = Publisher
    def get_context_data(self,**kwargs):
        context = super(PublisherDetail,self).get_context_data(**kwargs)
        context['book_list'] = Book.objects.all()
        return context'''

class BookList(ListView):
    queryset = Book.objects.order_by('publication_date')
    context_object_name = 'book_list'

class YaziBookList(ListView):
    context_object_name = 'bookOnePublisher_list'
    queryset = Book.objects.filter(publisher__name='yazi')
    template_name = 'test_demo/yazi_list.html'

class PublisherBookList(ListView):
    context_object_name = 'publisherbook_list'
    template_name = 'test_demo/books_by_publisher.html'
    def get_queryset(self):
        self.publisher = get_object_or_404(Publisher, name=self.args[0])
        return Book.objects.filter(publisher=self.publisher)


class AuthorCreate(CreateView):
    model = Author
    context_object_name = 'author_list'
    fields = ['name']
class AuthorUpdate(UpdateView):
    model = Author
    fields = ['name']
class AuthorDelete(DeleteView):
    model = Author
    #fields = ['name']
    success_url = reverse_lazy('author-list')


#FormView:A view for displaying a form, and rendering a template response.
class ContactView(FormView):
    template_name = 'contact.html'
    form_class = ContactForm
    success_url = '/thanks'

    def form_valid(self,form):  #This method is called when valid form data has been POSTed.
        form.send_email()
        return super(ContactView,self).form_valid(form)

from django.conf.urls import url

from . import views
from django.views.generic import TemplateView

app_name = 'test_demo'

urlpatterns = [
    url(r'^publishers/$', views.PublisherList.as_view(), name='publishers'),
    #url(r'^publishers/$', views.PublisherDetail.as_view(), name='publishers'),
    url(r'^books/$', views.BookList.as_view(), name='books'),
    url(r'^yazibooks/$', views.YaziBookList.as_view(), name='yazibooks'),
    url(r'^eachpublisherbooks/([\w-]+)/$', views.PublisherBookList.as_view()),
    url(r'author/add/$', views.AuthorCreate.as_view(), name='author-add'),
    url(r'author/(?P<pk>[0-9]+)/$', views.AuthorUpdate.as_view(), name='author-update'),
    url(r'author/(?P<pk>[0-9]+)/delete/$', views.AuthorDelete.as_view(), name='author-delete'),

]

from django.conf.urls import url

from . import views

app_name = 'website'

urlpatterns = [
    url(r'^manual/$',views.manual,name='manual'),
    url(r'^reference/$',views.reference,name='reference'),
    #url(r'^faq/$',views.faq,name='faq'),
    url(r'^contact$',views.contact,name='contact'),
    url(r'^vios_intro$',views.vios_intro,name='vios_intro'),
]

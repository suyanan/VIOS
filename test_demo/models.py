from __future__ import unicode_literals
from django.urls import reverse
from django.db import models

# Create your models here.

class Publisher(models.Model):
    name = models.CharField(max_length=30)
    address = models.CharField(max_length=50)
    city = models.CharField(max_length=60)
    state_province = models.CharField(max_length=30)
    country = models.CharField(max_length=50)
    website = models.URLField()
    class Meta:
        ordering = ["-name"]
    def __unicode__(self):
        return self.name

class Author(models.Model):
    salutation = models.CharField(max_length=10)
    name = models.CharField(max_length=200)
    email = models.EmailField()
    headshot = models.ImageField(upload_to='author_headshots')
    def __unicode__(self): # __unicode__ on Python 2
        return self.name
    def get_absolute_url(self):
        return reverse('author-detail', kwargs={'pk': self.pk})

class Book(models.Model):
    title = models.CharField(max_length=100)
    authors = models.ManyToManyField('Author')
    publisher = models.ForeignKey(Publisher, on_delete=models.CASCADE)
    publication_date = models.DateField()
    def __unicode__(self): # __unicode__ on Python 2
        return self.title

class Animal(models.Model):
    name = models.CharField(max_length=100)
    sound = models.CharField(max_length=100)
    def speak(self):
        uni_string = 'The '+self.name +' says "'+ self.sound+'"'
        return str(uni_string)

from __future__ import unicode_literals

from django.db import models
from django.core.files.storage import FileSystemStorage
from VIOS import settings
import time
import os
# Create your models heres.

class LoginUsername(models.Model):
    username = models.CharField(max_length=255, blank=True)

class OverwriteStorage(FileSystemStorage):
    def get_available_name(self, name,max_length=None):
        if self.exists(name):
            os.remove(os.path.join(settings.MEDIA_ROOT, name))
        return name


def user_directory_path(instance,filename):
	#dirpath = time.strftime("%Y")+"/"+time.strftime("%m")+"/"+time.strftime("%d")+"/"+"raw_datas"+instance.username
    dirpath = time.strftime("%Y")+"/"+time.strftime("%m")+"/"+"raw_datas"+instance.username
    print '88888888'
    print dirpath
    return '{0}/{1}/'.format(dirpath,filename)

class SequencingFiles(models.Model):
    username = models.CharField(max_length=255, blank=True)
    SeqFiles = models.FileField(upload_to=user_directory_path,storage=OverwriteStorage())
    uploaded_at = models.DateTimeField(auto_now_add=True)

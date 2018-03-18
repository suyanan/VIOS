from django import forms

from .models import SequencingFiles,LoginUsername

class SequencingFilesForm(forms.ModelForm):
    class Meta:
	    model = SequencingFiles
	    fields = ('SeqFiles', )

class LoginUsernameForm(forms.ModelForm):
    class Meta:
	    model = LoginUsername
	    fields = ('username', )

# -*- coding: utf-8 -*-
# Generated by Django 1.10.3 on 2017-04-28 02:23
from __future__ import unicode_literals

from django.db import migrations, models
import ngs.models


class Migration(migrations.Migration):

    dependencies = [
        ('ngs', '0002_auto_20170410_0700'),
    ]

    operations = [
        migrations.CreateModel(
            name='LoginUsername',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('username', models.CharField(blank=True, max_length=255)),
            ],
        ),
        migrations.AlterField(
            model_name='sequencingfiles',
            name='SeqFiles',
            field=models.FileField(storage=ngs.models.OverwriteStorage(), upload_to=ngs.models.user_directory_path),
        ),
    ]

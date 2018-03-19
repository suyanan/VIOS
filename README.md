# VIOS
Virus Identification On Site

[TOC]

![VIOS HOME PAGE](https://github.com/suyanan/VIOS/raw/master/INTRO/HOME-main.png)

## Pre-requirement
- Java 1.8
- Python 2.7
- Perl 5
- Python 3rd libraries

  > xlwt xlrd biopython numpy ReportLab matplotlib django1.11.3 django_extensions


There are mainly three steps, including short reads' mapping with Bowtie2, assembling with Velvet, contigs' mapping with blast+.
Deloy the bioinformatics softwares , and add the commond to PATH.
- [CD-HIT](http://weizhongli-lab.org/cd-hit/)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [velvet](https://www.ebi.ac.uk/~zerbino/velvet/)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

## DEMO DATA
adenovirus two samples:
-jiayan
-liublood

You need to upload the raw datas, after uploading the files were in the directory of "$PATH/VIOS/media/2018/02/raw_datas/".
You can get the them from [BaiduYun_demoData](https://pan.baidu.com/s/1hqsjfqi8-1ANb_ZmS8Sqgg).

## Database preparation
BEFORE, set the **DATABASE_HOME** in your environment variables.
```
export DATABASE_HOME=$YOUR_PATH/database
##genome directory
$DATABASE_HOME/genomes
##virus gene/genome directory
$DATABASE_HOME/blastDB
```
1. genome

   Go to the modul 用户自有参考数据库索引 on the web.
   Then , You can search species genome from [UCSC](http://hgdownload.soe.ucsc.edu/downloads.html) to get host "genome name" you wanna map aginst to remove the host genetics informations from the raw data reads. And, name the "genome name" with the readable name in the "host name" on the web.

   For example, you can add "hg38" and "human" in the input box attached a label with "genome name" and "host name", respectively.

   There are also some common species
   - model animal in the lab
   - poultry and livestock
   - primates (closer to human relations)

       | genome name |   host name   |
       |:-----------:|:-------------:|
       |    hg38     |     human     |
       |    mm10     |  house_mouse  |
       |     rn6     |  Norway_rat   |
       |   galGal5   |    chicken    |
       |   felCat8   |      cat      |
       |   susScr3   |      pig      |
       |   oryCun2   |    rabbit     |
       |   oviAri3   |     sheep     |
       |   bosTau8   |      cow      |
       |   equCab2   |     horse     |
       |   myoLuc2   |   microbat    |
       |   papAnu2   |    baboon     |
       |   panTro5   |  chimpanzee   |
       |   gorGor5   |    gorilla    |
       |   rheMac8   | rhesus_monkey |
2. virus gene/genome

   - Step1. pre-handle the virus database referring to the paper;
   - Step2. name it with "nt_viruses_final.fasta";
   - Step3. put the file into the **$DATABASE_HOME/blastDB/nt_viruses/final** directory;
   - Step4. Go to the the modelue 用户自有参考数据库索引, press the "Update Viruses" button.

   Besides, you can also custom the virus 's some **family** like Filoviridae which Ebola belongs to. The process method is like the virus handling method above, in addition, modify the variables(nt_viruses_family_name_list,nt_viruses_family_name_readable_list) in the **/VIOS/ngs/scripts/config_paras.py** script.

   Well, you can email at Su Yanan <suyanan1991@163.com> to get the virus database, I'll share with you through cloud storage.
   You can download nt_viruses_final database from [BaiduYun_database](https://pan.baidu.com/s/1lSpOnhK-uga0AEWrFpeWEw).

## Deploy the Django Project **VIOS**
  - **Step1. Config MYSQL database**
    ```
    $apt install mysql-server mysql-client libmysqlclient-dev  #REMEMBER ROOT PASSWORD
    $myslq -u root -p
    mysql>CREATE DATABSAE vios_db;
    mysql>USE vios_db;
    ```
    Besides, you can change MySQL to other database like Postgres.

  - **Step2. Clone the VIOS project code locally**
    ```
    #clone the VIOS codes locally:
    git clone git://github.com/suyanan/VIOS.git

    #Modify the /VIOS/VIOS/settings.py file
    ##replace those variables:USER/PASSWORD
    DATABASES = {
        'default': {
            'ENGINE': 'django.db.backends.mysql',
            'NAME': 'vios_db',
            'HOST':'127.0.0.1',
            'PORT':'3306',
            'USER':'root',
            'PASSWORD':'1234aaaa',
        }
    }
    ##ADD your IP
    ALLOWED_HOSTS = [
        u'127.0.0.1',
        u'localhost',
        u'your PC IP',
    ]
    ```

  - **Step3. start the VIOS project**
    ```
    $python /VIOS/manage.py runserver 0.0.0.0:8000

    #apply the migrations for VIOS
    $python manage.py migrate
    ```
    Open the [VIOS_dress](127.0.0.1:8000) with link **127.0.0.1:8000**

  - **Step4. ANALYSIS: database**
  - **Step5. ANALYSIS: pipeline**
  - **Step6. ANALYSIS: visualization**


## TOOL'S FUNCTIONDEMO ON BROWSER SITE

  ![ANALYSIS-uploadFile](https://github.com/suyanan/VIOS/raw/master/INTRO/ANALYSIS-1-uploadFiles.png)

  ![ANALYSIS-process](https://github.com/suyanan/VIOS/raw/master/INTRO/ANALYSIS-2-pipeline.png)

  ![ANALYSIS-results](https://github.com/suyanan/VIOS/raw/master/INTRO/ANALYSIS-3-results-a.png)

  ![ANALYSIS-updateDatabase](https://github.com/suyanan/VIOS/raw/master/INTRO/ANALYSIS-4-updateDatabase.png)


<br>
**COPY RIGHT!**

2018.02.09


---
The End

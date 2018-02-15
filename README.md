# Fettuccine_MSc_Group_Project
- Ensure you have the following modules install on your pc (by any means necessary e.g using pip install)
- Flask
- Pandas (Provides pytz, numpy)
- Pytz
- MySQLdb
- Numpy
- Biopython
- Pylab
- Matplotlib
- Pygraphviz (and graphviz, not via pip)

- You also need to download MAMP or AAMP for linux (which does the MYSQL hosting for you) otherwise there will be no database
- https://www.mamp.info/en/

# Complete
- Gantt Chart
- Architecture
- UI
- All Pages complete


link:
https://academic.oup.com/gbe/article/8/12/3485/2680041


# Useful Commands:
- git clone (link)
- git add (file)
- git commit -m "comment goes here"
- git push
- git pull

# README FOR INSTALL

------------------------------
USE python 32bit or 64bit. (different library files will be provided for each version)

Use virtualenv at your own discretion

Open up CMD terminal and run:
pip install -r /path/to/requirements.txt
------------------------------

Download the DB files from:

- Dropbox: https://www.dropbox.com/s/0ajx2p0t9wbe3ib/genome_data_export.sql?dl=0

Install Your Own MYSQL server of Choice (we Used MAMP for windows to make it simple to turn on / off)

Download https://www.graphviz.org/download/ V2.38 and install

The software was ran and tested using MySQK 5.7.20 in Windows 32-bit, 64-bit and Linux 64-bit.

Copy the MySQLdb (MySQL python connector) and pygraphviz folders + pylab files to your C:\Python27\Lib\site-packages\ or equivalent folders for your OS or if using virtualenv. (These files have been provided in the compressed zip format)

PyGraphviz May not work on 64bit python, python 64bit is reccomended as the 32bit will cause memory errors if you try to upload files > 300MB (however pygraphviz works on 32bit)

If Pygraphviz does not load, you are still able to run the website as it is not needed for the main sections of the site
(The SECOND custom tree will not load on the webpage after submission of ph file, however you can visit the url custom_tree2 > "localhost:5000/custom_tree2" to see how it would look if it had rendered properly with the functional plugin)

# MYSQL
-----------

Create a database called genome_data:

> CREATE DATABASE genome_data;

Import the SQL file into MySQL using:

> Path/to/Mysql.exe -u root -p genome_data < "path/to/genome_data.sql"

User can also create an identical database from scratch by
following the instructions  inside SQL_db_instructions folder.
This file also contains random data from our database that can
populate the new db.

---------------
# RUNNING FLASK
-----------

To run the software go to Website_Design and execute the commands based on OS:

# run flask on WINDOWS using:
OPEN COMMAND PROMPT, GO TO CORRECT DIRECTORY of file (cd ..\..)
- set FLASK_APP=project.py
- set FLASK_DEBUG=1
- flask run

# run flask on LINUX using:
- export FLASK_APP=project.py
- export FLASK_DEBUG=1
- flask run

----------------------

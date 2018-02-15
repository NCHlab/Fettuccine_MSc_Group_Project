from mod_tables.serverside.serverside_table import ServerSideTable
from mod_tables.serverside import herv_table_schemas
from mod_tables.serverside import l1_table_schemas
from flask import Flask, send_from_directory, render_template, request, url_for, flash
from datetime import datetime
from pytz import timezone
import os
import MySQLdb
import pandas as pd
import numpy as py
from Bio import SeqIO
import re
import csv
import hashlib
import json
#Connect to the datbase
try:
	connection = MySQLdb.connect(host="localhost",
						port=3306,
						user="root",
						passwd="root",
						db="genome_data")#genome_data
except:
	print "Failed to connect to db!!"
	print "Please ensure you have started your MySQL server and the db name is correct"
	print "and that the port is set to 3306 in your server"
	pass

#Get the information from the database
cur = connection.cursor()
cur.execute("SELECT peptide_id, family, superfamily, predicted_protein, sequence FROM `herv_prot_seqs`")
row_headers_erv=[x[0] for x in cur.description]
ervs = cur.fetchall()
json_data_erv=[]
cur.execute("SELECT peptide_id, family, protein, sequence FROM `l1_prot_seqs`")
row_headers_l1=[x[0] for x in cur.description]
l1s = cur.fetchall()
json_data_l1=[]
for result in l1s:
    json_data_l1.append(dict(zip(row_headers_l1,result)))
for result in ervs:
    json_data_erv.append(dict(zip(row_headers_erv,result)))


DATA_SAMPLE_ERV = json_data_erv
DATA_SAMPLE_L1 = json_data_l1

#return the data for serverside processing
class TableBuilder(object):

	def collect_data_serverside_erv(self, request):
		columns = herv_table_schemas.SERVERSIDE_TABLE_COLUMNS
		return ServerSideTable(request, DATA_SAMPLE_ERV, columns).output_result()

	def collect_data_serverside_l1(self, request):
		columns = l1_table_schemas.SERVERSIDE_TABLE_COLUMNS
		return ServerSideTable(request, DATA_SAMPLE_L1, columns).output_result()

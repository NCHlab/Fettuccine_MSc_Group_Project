from flask import Flask, send_from_directory, render_template, request, url_for
from datetime import datetime
from pytz import timezone
import os
import MySQLdb
import pandas as pd
import numpy as py
#import fasta
from Bio import SeqIO
import re
import xml.etree.ElementTree as ET
#from pyfastaq import sequences

app = Flask(__name__)
APP_ROOT = os.path.dirname(os.path.abspath(__file__))

# Connect to the database
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


cur = connection.cursor()


@app.route("/")
def indexpage():
	now = datetime.now(timezone('Europe/London'))
	return render_template("index.html", time = now)
	#return "The time is {}".format(now)

@app.route("/family_table")
def family_table():
	global HF
	global LF
	global DATA2
	global DATA3
	cur.execute("SELECT * FROM RT_data")
	table_family = cur.fetchall()
	d2=[]
	for row in table_family:
		d2.append({'Family':row[1], 'name': row[2], 'details': row[3], 'sequence': row[4]})
	df2=pd.DataFrame(d2)
	DATA2=df2.to_html()
	df2_fam= df2[df2.Family == 'HERV']
	df3_fam= df2[df2.Family == 'LINE1']
	HF=df2_fam.to_html()
	LF=df3_fam.to_html()

	cur.execute("SELECT * FROM prelim_data")
	table_family2 = cur.fetchall()
	d4=[]
	for row in table_family2:
		d4.append({'genoName':row[1], 'genoStart': row[2], 'genoEnd': row[3], 'genoLeft': row[4],
		 'strand': row[5], 'repName': row[6], 'repClass': row[7], 'repFamily': row[8],
		  'repStart': row[9], 'repEnd': row[10], 'repLeft': row[11]})
	df4=pd.DataFrame(d4)
	DATA3=df4.to_html()
	return render_template("family_table.html", data3=DATA3, data2=DATA2)

@app.route('/family_table_LINE1')
def family_table_LINE1():
    Sum="variable passed on"
    return render_template('family_table_LINE1.html',data2=LF)

@app.route('/family_table_HERV')
def family_table_HERV():
    Sum="variable passed on"
    return render_template('family_table_HERV.html',data2=HF)

@app.route('/family_table_O')
def family_table_O():
    Sum="variable passed on"
    return render_template('family_table_O.html',data3=DATA3)

@app.route("/distribution")
def distribution():
    return render_template("distribution.html")

@app.route("/AA_seq_list")
def AA_seq_list():
	return render_template("AA_seq_list.html")

@app.route("/relationship_AA")
def relationship_AA():
	return render_template("relationship_AA.html")

@app.route("/peptide_seq_ident", methods=["GET","POST"])
def peptide_seq_ident():
	global empty_error
	global rows_count
	global result_seq
	global no_match
	result_seq=""
	rows_count = ""
	error_empty2 = "This is an empty file! Please upload a populated FASTA file"
	no_match= "No Match was found"

	#fastaseq = request.form["fasta_content"]
	# If data has been submitted to the page i.e uploaded, then the POST method engages
	if request.method == "POST":
		if request.form["fasta_content"] != "":
			fastaseq = request.form["fasta_content"]
			rows_count = cur.execute("SELECT Family, Sequence FROM prelim2 WHERE Sequence = %s", fastaseq)
			cur.execute("SELECT Family, Sequence FROM prelim2 WHERE Sequence = %s", fastaseq)
			result_seq =  cur.fetchall()
			if not cur.rowcount:
			  return render_template("peptide_seq_ident.html", result_family=no_match)
			else:
			  return render_template("peptide_seq_ident.html", result_family=result_seq[0][0], result_seq1=result_seq[0][1])

		elif 'file2' in request.files:
			target = os.path.join(APP_ROOT, "sequence_ident/")

			if not os.path.isdir(target):
				os.mkdir(target)

			  # Saves the file to the target folder explicitly mentioned earlier
			for file in request.files.getlist("file2"):
				filename = file.filename
				destination = "/".join([target, filename])
				file.save(destination)
			  # Goes through fasta file
			global recordID


			seqfile = SeqIO.parse(filename, "fasta")
			if os.stat(filename).st_size == 0:
				return render_template("peptide_seq_ident.html", empty = error_empty2)
			else:
				for record in SeqIO.parse(filename, "fasta"):
					recordID = record.seq
					rows_count = cur.execute("SELECT Family, Sequence FROM prelim2 WHERE Sequence = %s", recordID)
					cur.execute("SELECT Family, Sequence FROM prelim2 WHERE Sequence = %s", recordID)
					result_seq =  cur.fetchall()
					if not cur.rowcount:
					  return render_template("peptide_seq_ident.html", result_family=no_match)
					else:
					  return render_template("peptide_seq_ident.html", result_family=result_seq[0][0], result_seq1=result_seq[0][1])
		elif request.form["fasta_content"] == "":
			return render_template("peptide_seq_ident.html", empty = error_empty2)
		  #return render_template("peptide_seq_ident.html", result_family=result_seq, result_seq1=result_seq, emptyfile=error_empty)
	  #return render_template("peptide_seq_ident.html", result_family=result_seq[0][0], result_seq1=result_seq[0][1], emptyfile=error_empty)
	else:
	  	return render_template("peptide_seq_ident.html")


@app.route("/upload_peptide", methods=["GET","POST"])
def upload_peptide():
	global filename2
	global pepseq

	if request.method == "POST":
		target = os.path.join(APP_ROOT, "uploaded/")

		if not os.path.isdir(target):
			os.mkdir(target)

		for file in request.files.getlist("file"):
			filename2 = file.filename
			destination = "/".join([target, filename2])
			file.save(destination)
		pepseq=[]
		tree = ET.parse("Galaxy18.xml")
		root = tree.getroot()
		#root = ET.fromstring(Galaxy18_as_string)
		for Peptide in root.findall("Peptide"):
			pepseq.append(Peptide.find("PeptideSequence").text)
		#for i in range(0,4):
		#	pepseq.append(i)

		return render_template("uploaded.html", matches1 = pepseq)
		#return uploaded() #render_template("uploaded.html")
	else:
		return render_template("upload_peptide.html")

@app.route("/uploaded")
def uploaded():

	return render_template("uploaded.html", matches1 = pepseq)

@app.route("/expression_atlas")
def atlas():
	return render_template("expression_atlas.html")

@app.route("/documentation")
def documentation():
	return render_template("documentation.html")

@app.route("/about_us")
def about_us():
	return render_template("about_us.html")

###############################################

@app.route("/profile/<name>")
def profile(name):
    return render_template("profile2.html", name=name)


@app.route('/user/<username>')
def show_user_profile(username):
	return render_template("style.html", username=username)
    # show the user profile for that user
    #return 'User %s' % username

@app.route('/login', methods=['GET', 'POST'])
def login():
    if request.method == 'POST':
        pass
    else:
        show_the_login_form()

@app.route("/visualise")
def visualise():
    return "lolnope"



#@app.route("/static")
#def static1():
#	url_for('static', filename='style.css')

@app.route('/<path:path>')
def static_file(path):
    return app.send_static_file(path)

#if __name__ == "__main__":
#	app.run(debug=True)

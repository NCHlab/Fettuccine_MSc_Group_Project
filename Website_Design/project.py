from flask import Flask, send_from_directory, render_template, request, url_for, flash
from datetime import datetime
from pytz import timezone
import os
import MySQLdb
import pandas as pd
import numpy as py
from Bio import SeqIO
import re
#import xml.etree.ElementTree as ET
#from pyfastaq import sequences
#import fasta

app = Flask(__name__)
APP_ROOT = os.path.dirname(os.path.abspath(__file__))
#app.secret_key = 'SUPER_SECRET_KEY_BIO_PROJECT'

ALLOWED_EXTENSIONS = set(["xml", "mzid", "mzTab", "mztab" ,"fasta"])
ALLOWED_EX_XML = set(["xml", "mzid"])
ALLOWED_EX_MZTAB = set(["mzTab", "mztab"])
ALLOWED_EX_FASTA = set(["fasta"])

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


cur = connection.cursor() # May need to open it in each function instead of globally


@app.route("/")
def indexpage():
	now = datetime.now(timezone('Europe/London'))
	return render_template("index.html", time = now)
	#return "The time is {}".format(now)
@app.route("/family_table")

def family_table():
	cur.execute("SELECT Family, Repeat_Name, Counts FROM `HERV_Groupby_Families`")
	HERV_GbF_rows = cur.fetchall()
	return render_template("family_table.html", data=HERV_GbF_rows)

@app.route('/family_table_LINE1')
def family_table_LINE1():
	cur.execute("SELECT Repeat_Name, Counts FROM `L1_Groupby_Repeats`")
	HERV_LbR_rows = cur.fetchall()
	return render_template("family_table_LINE1.html", data=HERV_LbR_rows)

@app.route("/family_table_HERV")
def family_table_HERV():
	cur.execute("SELECT Family, Repeat_Name, Counts FROM `HERV_Groupby_Families`")
	HERV_GbF_rows = cur.fetchall()
	wisit="Family"
	return render_template("family_table_HERV.html", data="", firstcolumn="")

@app.route("/family_table_HERV_Families")
def family_table_HERV_Families():
	cur.execute("SELECT Family, Repeat_Name, Counts FROM `HERV_Groupby_Families`")
	HERV_GbF_rows = cur.fetchall()
	wisit="Family"
	return render_template("family_table_HERV.html", data=HERV_GbF_rows, firstcolumn=wisit)

@app.route("/family_table_HERV_Superfamilies")
def family_table_HERV_Superfamilies():
	cur.execute("SELECT Superfamily, Repeat_Name, Counts FROM `HERV_Groupby_Superfamilies`")
	HERV_GbS_rows = cur.fetchall()
	wisit="Superfamily"
	return render_template("family_table_HERV.html", data=HERV_GbS_rows, firstcolumn=wisit)

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
	global recordID
	global result_seq_one
	global data11
	result_seq_one=[]
	result_seq_multi=[]
	result_seq=""
	rows_count = ""
	error_empty2 = "This is an empty file! Please upload a populated FASTA file"
	no_match= "No Match was found"


	# If data has been submitted to the page i.e uploaded, then the POST method engages
	if request.method == "POST":
		# If the Textbox has been filled with peptide sequences, DB is searched for matching sequence and FAMILY + sequence returned
		# Otherwise if no match found, displays no match
		if request.form["fasta_content"] != "":
			fastaseq = request.form["fasta_content"]
			rows_count = cur.execute("SELECT Family, Sequence FROM prelim2 WHERE Sequence = %s", fastaseq)
			cur.execute("SELECT Family, Sequence FROM prelim2 WHERE Sequence = %s", fastaseq)
			result_seq =  cur.fetchall()
			result_seq_one.append(result_seq)
			DF_PD=pd.DataFrame(result_seq_one)
			result_seq_df=DF_PD.to_html()
			if not cur.rowcount:
			  return render_template("peptide_seq_ident.html", result_family=no_match)
			else:
			  return render_template("peptide_seq_ident.html", result_family=result_seq_df)
		# If a file has been uploaded (file2 - name of upload form), this if statement occurs
		elif 'file2' in request.files:
			# Creates a path to the specified folder
			target = os.path.join(APP_ROOT, "sequence_ident/")
			# Creates directory if it doesnt already exist
			if not os.path.isdir(target):
				os.mkdir(target)

			  # Saves the file to the target folder explicitly mentioned earlier
			for file in request.files.getlist("file2"):
				filename = file.filename
				destination = "/".join([target, filename])
				file.save(destination)

			if os.getcwd() == APP_ROOT:
				os.chdir("sequence_ident")
			elif os.getcwd() == APP_ROOT+"\sequence_ident":
				pass
			elif os.getcwd() == APP_ROOT+"\uploaded":
				os.chdir("..\sequence_ident")

			#fpath = os.path.join(direct, filename)
			# Goes through fasta file and checks whether if its empty
			seqfile = SeqIO.parse(filename, "fasta")
			if os.stat(filename).st_size == 0:
				return render_template("peptide_seq_ident.html", empty = error_empty2)
			else:
				# if file contains information, retrieve data from the DB,
				# if request from DB does not have data returned (rowcount=0), show no match,
				# otherwise display Family + sequence
				record_list = list(SeqIO.parse(filename, "fasta"))
				if len(record_list) == 1:
					#If only 1 fasta sequence in file
					for record in SeqIO.parse(filename, "fasta"):
						recordID = record.seq
						rows_count = cur.execute("SELECT Family, Sequence FROM prelim2 WHERE Sequence = %s", recordID)
						cur.execute("SELECT Family, Sequence FROM prelim2 WHERE Sequence = %s", recordID)
						result_seq =  cur.fetchall()
						result_seq_one.append(result_seq)
						DF_PD=pd.DataFrame(result_seq_one)
						result_seq_df=DF_PD.to_html()
					if not cur.rowcount:
					  return render_template("peptide_seq_ident.html", result_family=no_match)
					else:
						return render_template("peptide_seq_ident.html", result_family=result_seq_df)#, result_seq1=result_seq_df)#result_seq[0][1])
				else:
					for record in SeqIO.parse(filename, "fasta"):
						# Loop which iterates through ever fasta sequence and appends the results
						recordID = record.seq
						rows_count = cur.execute("SELECT Family, Sequence FROM prelim2 WHERE Sequence = %s", recordID)
						cur.execute("SELECT Family, Sequence FROM prelim2 WHERE Sequence = %s", recordID)
						result_seq =  cur.fetchall()
						if cur.rowcount > 0:
							result_seq_multi.append(result_seq)
					DF_PD=pd.DataFrame(result_seq_multi)
					result_seq_df=DF_PD.to_html()
					if not cur.rowcount:
					  return render_template("peptide_seq_ident.html", result_family=no_match)
					else:
					  return render_template("peptide_seq_ident.html", result_family=result_seq_df, data=result_seq_multi)#result_seq_multi)
		elif request.form["fasta_content"] == "":
			return render_template("peptide_seq_ident.html", empty = error_empty2)
	else:
		return render_template("peptide_seq_ident.html")



@app.route("/upload_peptide", methods=["GET","POST"])
def upload_peptide():
	global filename2
	global pepseq
	global list_of_matches
	global list_of_pep_seqs
	global empty_error
	global rows_count
	global result_seq
	global no_match
	global recordID
	global result_seq_multi
	result_seq_multi=[]
	result_seq=""
	rows_count = ""
	error_empty2 = "This is an empty file! Please upload a populated FASTA file"
	no_match= "No Match was found"
	list_of_matches = []
	list_of_pep_seqs = []

	if request.method == "POST":
		target = os.path.join(APP_ROOT, "uploaded/")

		if not os.path.isdir(target):
			os.mkdir(target)

		for file in request.files.getlist("file"):
			filename2 = file.filename
			destination = "/".join([target, filename2])
			file.save(destination)

		if os.getcwd() == APP_ROOT:
			os.chdir("uploaded")
		elif os.getcwd() == APP_ROOT+"\uploaded":
			pass
		elif os.getcwd() == APP_ROOT+"\sequence_ident":
			os.chdir("..\uploaded")

		if filename2.rsplit('.', 1)[1].lower() in ALLOWED_EX_XML:
			file_mz_seq = open(filename2, "r")
			whole_file = file_mz_seq.read()

			regex = r">[a-zA-z]+"
			matches = re.finditer(regex, whole_file)
			for matchNum, match in enumerate(matches):
			    matchNum = matchNum + 1
			    list_of_matches.append(match.group())

			list_of_pep_seqs = [character.replace('>', '') for character in list_of_matches]

		elif filename2.rsplit('.', 1)[1].lower() in ALLOWED_EX_MZTAB:

			file_mztab_seq = open(filename2, "r")
			whole_file = file_mztab_seq.read()

			regex = r"[A-Z]+\S\B[A-Z]"
			matches = re.finditer(regex, whole_file)
			for matchNum, match in enumerate(matches):
			    matchNum = matchNum + 1
			    list_of_matches.append(match.group())
			list_of_words = ["UNIMOD", "PSM", "COM", "TRAQ", "MTD", "PRIDE"]
			mztab_seq_mixed = [word for word in list_of_matches if word not in list_of_words]
			list_of_pep_seqs = [sequence for sequence in mztab_seq_mixed if len(sequence) > 5]

		else:
			result_seq_multi = "Incorrect filetype uploaded! Please Upload a MzIdent or mzTab formatted file"
			return render_template("upload_peptide.html", result_family=result_seq_multi)

		if len(list_of_pep_seqs) == 1:
			#If only 1 fasta sequence in file
			for seqs in list_of_pep_seqs:
				rows_count = cur.execute("SELECT Family, Sequence FROM prelim2 WHERE Sequence = %s", seqs)
				cur.execute("SELECT Family, Sequence FROM prelim2 WHERE Sequence = %s", seqs)
				result_seq =  cur.fetchall()
			if not cur.rowcount:
			  return render_template("upload_peptide.html", result_family=no_match)
			else:
				return render_template("upload_peptide.html", result_family=result_seq[0][0], result_seq1=result_seq[0][1])
		else:
			for seqs in list_of_pep_seqs:
				rows_count = cur.execute("SELECT Family, Sequence FROM prelim2 WHERE Sequence = %s", seqs)
				cur.execute("SELECT Family, Sequence FROM prelim2 WHERE Sequence = %s", seqs)
				result_seq =  cur.fetchall()
				if cur.rowcount > 0:
					result_seq_multi.append(result_seq)
			DF_PD=pd.DataFrame(result_seq_multi)
			result_seq_df=DF_PD.to_html()
			return render_template("upload_peptide.html", result_family=result_seq_df)#result_seq_multi)
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

from flask import Flask, send_from_directory, render_template, request, url_for, flash
import os
import os.path
import MySQLdb
import pandas as pd
import numpy as py
from Bio import SeqIO
from Bio import Phylo
import re
import csv
import hashlib
import json
import cgi

#Serverside processing used for huge protein tables from DB
#Serverside processing: https://github.com/SergioLlana/datatables-flask-serverside
from mod_tables.models import TableBuilder #Serverside processing - Gets the protein data from DB (mod_tables/models.py)

import matplotlib
import matplotlib.pyplot as plt
import timeit

#Packages for "Create your own tree function"
try:
    import pylab
    import pygraphviz
except:
    print "PLEASE INSTALL GRAPHVIZ PROGRAM FROM THE WEBSITE, NOT PIP !!!!!!!!!!!!!!!!!!!!!"
    print "PLEASE INSTALL pip install pygraphviz !!!!!!!!!!!!!!!!!!!!!!"
    print "PLEASE INSTALL pip install matplotlib !!!!!!!!!!!!!!!!!"

#Serverside processing:
table_builder = TableBuilder()

app = Flask(__name__)
from common.routes import main #Serverside processing - Creates a root for the serverside tables (common/routes.py)
from mod_tables.controllers import tables #Serverside processing - Convert DB data into json obsects
app.register_blueprint(main) #roots defined above
app.register_blueprint(tables) #roots defined above
APP_ROOT = os.path.dirname(os.path.abspath(__file__))

#Safety measures
ALLOWED_EXTENSIONS = set(["xml", "mzid", "mzTab", "mztab" ,"fasta"])
ALLOWED_EX_XML = set(["xml", "mzid"])
ALLOWED_EX_MZTAB = set(["mzTab", "mztab"])
ALLOWED_EX_FASTA = set(["fasta","fa", "faa"])
ALLOWED_EX_NWK = set(["nwk","newick", "ph"])

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
cur2 = MySQLdb.cursors.SSCursor(connection)

@app.route("/")
def indexpage():
    # Gets the time and count data from MySQL table to display in browser

    cur = connection.cursor()
    cur.execute("SELECT COUNT(*) FROM herv_repeats")
    herv_index_count = cur.fetchall()
    cur.execute("SELECT COUNT(*) FROM l1_repeats")
    l1_index_count = cur.fetchall()
    cur.execute("SELECT ((SELECT COUNT(*) FROM l1_repeats)+(SELECT COUNT(*) FROM herv_repeats));")
    total_index_count = cur.fetchall()
    cur.close()
    # Formatting the code to display it neater (creates comma seperators e.g 1,000,000)
    herv_index_count ="{:,.0f}".format(int(herv_index_count[0][0]))
    l1_index_count ="{:,.0f}".format(int(l1_index_count[0][0]))
    total_index_count ="{:,.0f}".format(int(total_index_count[0][0]))
    return render_template("index.html", herv_index = herv_index_count, l1_index = l1_index_count, total_index = total_index_count)

#Family Table - LINE1
@app.route('/family_table_LINE1')
def family_table_LINE1():
    #Get LINE-1 repeats from DB
    cur.execute("SELECT repeat_name, no_repeats, no_proteins_found, proteins_found FROM `l1_groupby_repeat_names` ORDER BY `l1_groupby_repeat_names`.`no_proteins_found` DESC")
    HERV_LbR_rows = cur.fetchall()
    #Return them to the website
    return render_template("family_table_LINE1.html", data=HERV_LbR_rows)

#Family table - HERV Superfamilies
@app.route("/family_table_HERV_Superfamilies")
def family_table_HERV_Superfamilies():
    #Get HERV repeats from DB
    cur.execute("SELECT superfamily, description, no_repeats, no_proteins_found, proteins_found FROM `herv_groupby_superfamilies`")
    HERV_GbS_rows = cur.fetchall()
    protlist = []
    countlist = []

    #Get values for the Bar Chart
    for ele in HERV_GbS_rows:
        countlist.append(int(ele[2])) #Values for the chart
        protlist.append(int(ele[3])) #Legend names

    #Return data for table and chart to the html
    return render_template("family_table_HERV_Superfamilies.html", data=HERV_GbS_rows, arr_sc=countlist, arr_sc_2=protlist)

#Family table - HERV SUperfamilies
@app.route("/family_table_HERV_Families")
def family_table_HERV_Families():
    #Get HERV families from DB
    cur.execute("SELECT * FROM `herv_groupby_families`")
    HERV_GbF_rows = cur.fetchall()
    #Return data for the tables to the html
    return render_template("family_table_HERV_Families.html", data=HERV_GbF_rows)

#Family table - HERV Repeats
@app.route("/family_table_HERV_Repeats")
def family_table_HERV_Repeats():
    #Get HERV Repeats from DB
    cur.execute("SELECT * FROM `herv_groupby_repeat_names`")
    HERV_GbR_rows = cur.fetchall()
    wisit="Repeat Name"
    #Return data for the tables to the html
    return render_template("family_table_HERV_Repeats.html", data=HERV_GbR_rows)

#Distribution
@app.route("/distribution")
def distribution():
    #Get HERV Repeats from DB
    cur.execute("SELECT counts FROM `herv_chromosome_count`") #select the number of HERV repeats for each chromosome
    Herv_count=cur.fetchall()
    H=[element for h in Herv_count for element in h] #removes tuple
    cur.execute("SELECT counts FROM `l1_chromosome_count`") #select the number of LINE1 repeats for each chromosome
    L1_count=cur.fetchall()
    L=[element for l in L1_count for element in l] #removes tuple
    #Return data for the tables to the html
    return render_template("distribution.html", H=H, L=L)

#Protein Table - HERV
@app.route("/protein_table_HERV")
def serverside_table():
    #Opens the url and displays the serverside table
    #This is serverside processing - See above
    return render_template("erv_proteins.html")

#Protein Table - LINE-1
@app.route("/protein_table_L1")
def serverside_table2():
    #Opens the url and displays the serverside table
    #This is serverside processing - See above
    return render_template("l1_proteins.html")

@app.route("/relationship_AA", methods=["GET", "POST"])
def relationship_AA():
    if request.method == "POST":

        #HERV tree image
        if request.form["rv_button"] == "HERV":
            return render_template("herv_rv.html")

        #LINE-1 tree image
        elif request.form["rv_button"] == "LINE1":
            return render_template("line1_rv.html")
    else:
        return render_template("relationship_AA.html")

@app.route("/herv_rv")
def herv_rv():
    return render_template("herv_rv.html")

@app.route("/herv_rv1")
def herv_rv1():
    return render_template("herv_rv1.html")

@app.route("/AA_seq_list")
def AA_seq_list():
    return render_template("AA_seq_list.html")
    
@app.route("/custom_tree2")
def custom_tree2():
    #This function only required to manually access the page for ctrl + f5 refresh, to refresh the page to load new image
    #This is because custom_tree by default returns AA_relationship.html when a GET method is used.
    return render_template("custom_tree.html")

@app.route("/custom_tree", methods=["GET", "POST"])
def custom_tree():
    if request.method == "POST":

        if 'file_rv' in request.files:
            # Creates a path to the specified folder
            target = os.path.join(APP_ROOT, "static/assets/img/customtree/")
            # Creates directory if it doesnt already exist
            if not os.path.isdir(target):
                os.mkdir(target)

              # Saves the file to the target folder explicitly mentioned earlier
            for file in request.files.getlist("file_rv"):
                filename3 = file.filename
                destination = "/".join([target, filename3])
                file.save(destination)

            if "." not in filename3 or filename3.rsplit('.', 1)[1].lower() not in ALLOWED_EX_NWK:
                error_msg = "Incorrect filetype uploaded! Please Upload a Newick formatted file"
                return render_template("relationship_AA.html", error_msg=error_msg)

            if filename3.rsplit('.', 1)[1].lower() in ALLOWED_EX_NWK:
                dpi_type = int(request.form.get('dpi_type'))

                    # Checks the current directory and moves to the correct folder
                if os.getcwd() == APP_ROOT:
                    os.chdir("static/assets/img/customtree/")
                elif os.getcwd() == APP_ROOT+"/static/assets/img/customtree/":
                    pass
                elif os.getcwd() == APP_ROOT+"\uploaded":
                    os.chdir("../static/assets/img/customtree/")
                elif os.getcwd() == APP_ROOT+"\sequence_ident":
                    os.chdir("../static/assets/img/customtree/")

                tree = Phylo.read(filename3, 'newick')
                tree.ladderize()   # Flip branches so deeper clades are displayed at top
                Phylo.draw(tree)

                # Only USE THIS CODE IF USING 32BIT PYTHON AND YOU HAVE PYGRAPHVIZ INSTALLED
                # try:
                #     Phylo.draw_graphviz(tree)
                # except:
                #     pass
                #
                # pylab.savefig("customtree.png", dpi=dpi_type)
                return render_template("custom_tree.html")
    else:
        return render_template("relationship_AA.html")

@app.route("/line1_rv")
def line1_rv():
    return render_template("line1_rv.html")

#Functions needed for FASTA and MZID-MZTAB identifiers

#Split a list into chunks
def chunkify(lst,n):
    return [lst[i::n] for i in xrange(n)]

#Remove identical elements of a list
def py_unique(data):
    return list(set(data))

#Protein sequence identifier
@app.route("/peptide_seq_ident", methods=["GET","POST"])
def peptide_seq_ident():
    result_seq_one = []
    result_seq_multi = []
    result_seq = ""
    rows_count = ""
    error_empty2 = "This is an empty file! Please upload a populated FASTA file"
    err_img = "/static/assets/img/error_img.png"
    no_match = "No Match was found!"
    error_fasta_1 = "Please Enter only 1 fasta sequence in the search bar, for multiple use upload function"
    error_fasta_2 = "Please Remove the Headers and only Search using the peptide Sequences or use the upload function"


    # If data has been submitted to the page i.e uploaded, then the POST method engages
    if request.method == "POST":
        # Change the directory to the root location of the python file
        os.chdir(APP_ROOT)

        # If the Textbox has been filled with peptide sequences, DB is searched for matching sequence and FAMILY + sequence returned
        # Otherwise if no match found, displays no match
        if request.form["fasta_content"] != "":
            fasta_check = request.form["fasta_content"]
            fasta_check = fasta_check.count(">")

            # Checks to see if a header exists and rejects the data if it does
            if fasta_check > 1:
                return render_template("peptide_seq_ident.html", empty = error_fasta_1, err_img = err_img)
            elif fasta_check == 1:
                fastaseq = request.form["fasta_content"]
                return render_template("peptide_seq_ident.html", empty = error_fasta_2, err_img = err_img)

            else:
                # Saves the peptide sequences and removes all spaces, tabs and new lines
                # Does a MYSQL query for the correctly formatted peptide sequence
                # and returns a result if cur.rowcount is not empty
                fastaseq = request.form["fasta_content"]
                fastaseq = fastaseq.replace("\n","").replace("\r","").replace(" ","")
                cur = connection.cursor()
                rows_count = cur.execute("SELECT family, sequence FROM all_prot_seqs WHERE sequence = %s", [fastaseq])
                cur.execute("SELECT family, sequence FROM all_prot_seqs WHERE sequence = %s", [fastaseq])
                result_seq =  cur.fetchall()
                cur.close()
                if not cur.rowcount:
                  return render_template("peptide_seq_ident.html", empty = no_match, err_img = err_img)
                else:
                  return render_template("peptide_seq_ident.html", data1=result_seq)

        # If a file has been uploaded (file2 - ID name of upload form), this if statement occurs
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

            # Checks the current directory and moves to the correct folder
            if os.getcwd() == APP_ROOT:
                os.chdir("sequence_ident")
            elif os.getcwd() == APP_ROOT+"\sequence_ident":
                pass
            elif os.getcwd() == APP_ROOT+"\uploaded":
                os.chdir("..\sequence_ident")


            # Goes through fasta file and checks whether the file is empty
            seqfile = SeqIO.parse(filename, "fasta")
            if os.stat(filename).st_size == 0:
                return render_template("peptide_seq_ident.html", empty = error_empty2, err_img = err_img)


            # Checks the file extension, if it is not in the allowed fasta format, it is rejected
            elif "." not in filename or filename.rsplit('.', 1)[1].lower() not in ALLOWED_EX_FASTA:
                result_seq_multi = "Incorrect filetype uploaded! Please Upload a Fasta formatted file"
                return render_template("peptide_seq_ident.html", empty=result_seq_multi)


            else:
                # if file contains information, retrieve data from the DB,
                # if request from DB does not have data returned (rowcount=0), show no match,
                # otherwise display Family + sequence
                record_list = list(SeqIO.parse(filename, "fasta"))
                temp_list=[]

                for rec in record_list:
                    temp_list.append(str(rec.seq))

                if len(record_list) == 1:
                    #If only 1 fasta sequence in file

                    for record in SeqIO.parse(filename, "fasta"):
                        recordID = record.seq
                        cur = connection.cursor()
                        rows_count = cur.execute("SELECT family, sequence FROM all_prot_seqs WHERE sequence = %s", [recordID])
                        cur.execute("SELECT family, sequence FROM all_prot_seqs WHERE sequence = %s", [recordID])
                        result_seq =  cur.fetchall()
                        cur.close()

                    if not cur.rowcount:
                      return render_template("peptide_seq_ident.html", empty=no_match, err_img = err_img)
                    else:
                        return render_template("peptide_seq_ident.html", data2=result_seq)
                else:
                    # Checks to see the size (number) of the Sequences
                    # If it is below 5000, a large query is sent
                    # otherwise the query is chunked into bits of 5000 and then MySQL query executed
                    record_list=temp_list
                    if len(record_list)<=5000:
                        query = "SELECT family, sequence FROM all_prot_seqs WHERE sequence = "
                        query = query+'"'+str(record_list[0])+'"'+' OR sequence="'
                        sep = '" OR sequence="'
                        query = query+sep.join(record_list[1:])
                        query = query+'"'
                        cur = connection.cursor()
                        cur.execute(query)
                        result_seq=cur.fetchall()
                        cur.close()

                    elif len(record_list)>5000:
                        #Split records into 5000 chunks
                        divide = (len(record_list)/5000)+1
                        record_list = chunkify(record_list, divide)
                        result_seq=[]

                        for ele in record_list:
                            query = "SELECT family, sequence FROM all_prot_seqs WHERE sequence = "
                            query = query+'"'+str(ele[0])+'"'+' OR sequence="'
                            sep = '" OR sequence="'
                            query = query+sep.join(ele[1:])
                            query = query+'"'
                            cur = connection.cursor()
                            cur.execute(query)
                            for ele in cur.fetchall():
                                result_seq.append(ele)
                        cur.close()

                # If no matches are found, the webpage is returned stating no match
                # otherwise if matches are found, data is returned and displayed
                if len(result_seq)==0:
                  return render_template("peptide_seq_ident.html", empty=no_match, err_img = err_img)
                else:
                    return render_template("peptide_seq_ident.html", data=result_seq)
        # if nothing is submitted via text box form, empty error raised
        elif request.form["fasta_content"] == "":
            return render_template("peptide_seq_ident.html", empty = error_empty2, err_img = err_img)
    else:
        # By default, this GET method is returned and displays 1000 sequences and families from the database
        cur = connection.cursor()
        cur.execute("SELECT family, sequence FROM all_prot_seqs LIMIT 0, 1000")
        result_seq =  cur.fetchall()
        cur.close()
        return render_template("peptide_seq_ident.html", data1=result_seq)


@app.route("/upload_peptide", methods=["GET","POST"])
def upload_peptide():
    result_seq_multi=[]
    result_seq_multi2=[]
    result_seq=""
    rows_count = ""
    no_match= "No Match was found"
    list_of_matches = []
    list_of_pep_seqs = []
    rowlist=[]
    rowlist2=[]
    result_seq_write=[]
    tissue_type_write=[]
    disease_type_write=[]
    result_seq2=[]
    hashed=str()
    err_img = "/static/assets/img/error_img.png"

    # Checks for post method (data submitted)
    start_time = timeit.default_timer()
    if request.method == "POST":
        #Move to root directory of project python file
        os.chdir(APP_ROOT)
        target = os.path.join(APP_ROOT, "uploaded/")

        #Creates the folder if it doesnt exist
        if not os.path.isdir(target):
            os.mkdir(target)

        #Saves the file into the specified location above
        for file in request.files.getlist("file"):
            filename2 = file.filename
            destination = "/".join([target, filename2])
            file.save(destination)

        #Checks current folder and opens the correct location
        if os.getcwd() == APP_ROOT:
            os.chdir("uploaded")
        elif os.getcwd() == APP_ROOT+"\uploaded":
            pass
        elif os.getcwd() == APP_ROOT+"\sequence_ident":
            os.chdir("..\uploaded")

        #Checks to see if any file has been uploaded, if it has not then error returned
        if 'file' not in request.files:
            result_seq_multi = "No File uploaded! Please Upload a MzIdent or mzTab formatted file"
            return render_template("upload_peptide.html", empty=result_seq_multi, err_img = err_img)

        # Checks to see for the presence of an extension.
        if "." not in filename2:
            result_seq_multi = "Incorrect filetype uploaded! Please Upload a MzIdent or mzTab formatted file"
            return render_template("upload_peptide.html", empty=result_seq_multi, err_img = err_img)

        # If the filetype is allowed, the file is read in and then regex match used to locate the peptide
        if filename2.rsplit('.', 1)[1].lower() in ALLOWED_EX_XML:
            file_mz_seq = open(filename2, "r")
            whole_file = file_mz_seq.read()

            # Regex used to locate the peptide sequences
            regex = r">[a-zA-z]+"
            matches = re.finditer(regex, whole_file)
            for matchNum, match in enumerate(matches):
                matchNum = matchNum + 1
                list_of_matches.append(match.group())

            list_of_pep_seqs = [character.replace('>', '') for character in list_of_matches]

        # Checks to see filetype for MZTAB format
        elif filename2.rsplit('.', 1)[1].lower() in ALLOWED_EX_MZTAB:

            file_mztab_seq = open(filename2, "r")
            whole_file = file_mztab_seq.read()

            # Regex used to locate the peptide sequences
            # Sequences also filtered out using some keywords typically found and a set size allowed
            regex = r"[A-Z]+\S\B[A-Z]"
            matches = re.finditer(regex, whole_file)
            for matchNum, match in enumerate(matches):
                matchNum = matchNum + 1
                list_of_matches.append(match.group())
            list_of_words = ["UNIMOD", "PSM", "COM", "TRAQ", "MTD", "PRIDE"]
            mztab_seq_mixed = [word for word in list_of_matches if word not in list_of_words]
            list_of_pep_seqs = [sequence for sequence in mztab_seq_mixed if len(sequence) > 5]

        else:

            #If the file extension does not match mztab or mzident, an error stating wrong filetype uploaded relayed back
            result_seq_multi = "Incorrect filetype uploaded! Please Upload a MzIdent or mzTab formatted file"
            return render_template("upload_peptide.html", empty=result_seq_multi, err_img = err_img)

        str_pep_seqs = ''.join(str(i) for i in list_of_pep_seqs) # converts list to string for hash checker
        hashed = str(hashlib.sha224(str_pep_seqs).hexdigest()) # A hash check of the found peptide sequences of the file is conducted providing a unique SHA224 ID

        if len(list_of_pep_seqs) == 1:
            #If only 1 sequence in file
            cur = connection.cursor()
            for seqs in list_of_pep_seqs:
                rows_count = cur.execute("SELECT family, sequence FROM all_prot_seqs2 WHERE sequence LIKE %s", ["%%"+seqs+"%%"])
                cur.execute("SELECT family, sequence FROM all_prot_seqs2 WHERE sequence LIKE %s", ["%%"+seqs+"%%"])
                result_seq2 =  cur.fetchall()
                cur.close()


            #return no match
            if not cur.rowcount:
                cur.close()
                return render_template("upload_peptide.html", empty=no_match, err_img = err_img)

        else:
            #If the inserted file is up to or equal to 5000 sequences

            if len(list_of_pep_seqs)<=5000:

                #A single mysql query for all sequences
                query2 = "SELECT family, sequence FROM all_prot_seqs2 WHERE sequence LIKE "
                query2 = query2+'"%'+str(list_of_pep_seqs[0])+'%"'+' OR sequence LIKE "%'
                sep = '%" OR sequence LIKE "%'
                query2 = query2+sep.join(list_of_pep_seqs[1:])+'%"'
                cur = connection.cursor()
                cur.execute(query2)
                result_seq2=cur.fetchall()
                cur.close()


            #If the inserted file is too big
            elif len(list_of_pep_seqs)>5000:
                divide = (len(list_of_pep_seqs)/5000)+1
                #split the sequences into lists of around 5000 sequences
                list_of_pep_seqs=chunkify(list_of_pep_seqs, divide)

                #query for each ~5000 sequences
                for ele in list_of_pep_seqs:
                    cur = connection.cursor()
                    query2 = "SELECT family, sequence FROM all_prot_seqs2 WHERE sequence LIKE "
                    query2 = query2+'"%'+str(ele[0])+'%"'+' OR sequence LIKE "%'
                    sep = '%" OR sequence LIKE "%'
                    query2 = query2+sep.join(ele[1:])+'%"'
                    cur.execute(query2)
                    for somet in cur.fetchall():
                        result_seq2.append(somet)
                cur.close()

            # If not result is found, display no match
            if len(result_seq2)==0:
                result_seq_multi="No match was found!"
                return render_template("upload_peptide.html", empty = result_seq_multi, err_img = err_img)

        # If a result is found, the data from dropdown boxes are saved into memory
        if len(result_seq2) != 0:
            tissue_type = str(request.form.get('tissue_type'))
            disease_type = str(request.form.get('disease_type'))


            # CHecks to see if the hash checker file is present otehrwise creates it
            if os.path.isfile(APP_ROOT+"/uploaded/hash_checker.csv"):
                pass
            else:
                open("hash_checker.csv", "w")

            #If the unique hash check is not duplicate, the data is saved for the atlas expression
            #otherwise the data is not saved (due to it already existing from the same source file)
            #if no match is found
            # The family matches, and the data from the dropdown boxes are saved into MySQL
            # Using MySQL query. This data can then be viewed in expression atlas.
            with open("hash_checker.csv", "r+b") as f:
                reader = csv.reader(f)
                writer = csv.writer(f)
                for row in reader:
                    # Checks if row is empty, if not, appends rows to memory
                    if not (row):
                        pass
                    else:
                        rowlist2.append(row[0])
                if hashed not in rowlist2:
                    rowlist.append(hashed)
                    writer.writerow(rowlist)
                    result_seq_multi2 = result_seq2
                    cur = connection.cursor()

                    for i in range(0, len(result_seq2)):
                        cur.execute("INSERT INTO exp_atlas(tissue_type, repeat_family, disease_type,sequence) VALUES (%s,%s,%s,%s);", (tissue_type, result_seq2[i][0], disease_type, result_seq2[i][1]))
                        cur.fetchall()
                    query3 = ("UPDATE exp_atlas_count SET counts = (SELECT COUNT(tissue_type) FROM exp_atlas);")
                    cur.execute(query3)
                    try:
                        query4 = ("UPDATE exp_atlas_count SET " + tissue_type + " = (SELECT COUNT(tissue_type) FROM exp_atlas WHERE tissue_type = '" + tissue_type + "');")
                        cur.execute(query4)
                    except:
                        pass
                    try:
                        query5 = ("UPDATE exp_atlas_disease_counts SET " + disease_type + " = (SELECT COUNT(disease_type) FROM exp_atlas WHERE disease_type = '" + disease_type + "');")
                        cur.execute(query5)
                    except:
                        pass
                    connection.commit()
                    cur.close()

            return render_template("upload_peptide.html", data=result_seq2)
    else:
        return render_template("upload_peptide.html")


@app.route("/expression_atlas", methods=["GET","POST"])
def atlas():
    cur = connection.cursor()
    cur.execute("SELECT GROUP_CONCAT(DISTINCT disease_type SEPARATOR ','),GROUP_CONCAT(DISTINCT tissue_type SEPARATOR ','), COUNT(tissue_type) AS RTs_no_found, concat(round(((SELECT COUNT(disease_type))/(SELECT counts from exp_atlas_count)* 100 )),'%'), GROUP_CONCAT(DISTINCT repeat_family SEPARATOR ',') AS RT_found,GROUP_CONCAT(DISTINCT sequence SEPARATOR ',') FROM exp_atlas GROUP BY disease_type;")
    overall_disease_percentage=cur.fetchall()

    cur.execute("SELECT tissue_type, COUNT(tissue_type) AS family_no_found, concat(round(((SELECT COUNT(tissue_type))/(SELECT counts from exp_atlas_count)* 100 )),'%'), GROUP_CONCAT(DISTINCT repeat_family SEPARATOR ',') AS family_found,GROUP_CONCAT(DISTINCT sequence SEPARATOR ',') FROM exp_atlas GROUP BY tissue_type;")
    overall_percentage=cur.fetchall()
    cur.close()
    return render_template("expression_atlas.html",overall_percentage=overall_percentage,overall_disease_percentage=overall_disease_percentage)

@app.route("/documentation")
def documentation():
    return render_template("documentation.html")

@app.route("/about_us")
def about_us():
    return render_template("about_us.html")


@app.route("/expression_atlas_2", methods=["GET","POST"])
def atlas2():
    if request.method == "POST":
        ind_repeat_disease = []
        individual_percentage = []
        cur = connection.cursor()
        try:
            disease_type1 = str(request.form.get('disease_type'))
            cur.execute("SELECT GROUP_CONCAT(DISTINCT disease_type SEPARATOR ','), tissue_type, COUNT(tissue_type) AS RTs_no_found, concat(round(((SELECT COUNT(disease_type))/(SELECT " + disease_type1 + " from exp_atlas_disease_counts)* 100 )),'%'), GROUP_CONCAT(DISTINCT repeat_family SEPARATOR ',') AS RT_found,sequence FROM exp_atlas WHERE disease_type = '" + disease_type1 + "' GROUP BY disease_type,tissue_type, repeat_family;")
            ind_repeat_disease=cur.fetchall()
        except:
            pass

        try:
            tissue_type1 = str(request.form.get('tissue_type'))
            cur.execute("SELECT tissue_type, COUNT(tissue_type) AS family_no_found, concat(round(((SELECT COUNT(repeat_family))/(SELECT " + tissue_type1 + " from exp_atlas_count)* 100 )),'%'), GROUP_CONCAT(DISTINCT repeat_family SEPARATOR ',') AS family_found,GROUP_CONCAT(DISTINCT sequence SEPARATOR ',') FROM exp_atlas WHERE tissue_type = '" + tissue_type1 + "' GROUP BY repeat_family;")
            individual_percentage=cur.fetchall()
        except:
            pass
        cur.close()
        return render_template("expression_atlas_2.html",individual_percentage=individual_percentage,ind_repeat_disease=ind_repeat_disease)
    return render_template("expression_atlas_2.html")

@app.errorhandler(404)
def page_not_found(e):
    # Returns error page if user goes to any page that is not valid
    return render_template('404.html'), 404

#if __name__ == "__main__":
#    app.run(debug=True)

#! /usr/bin/python
"""This script uses retrotransposon downloaded from UCSC Table browser - RepeatMasker track in .fasta file
exported by using the "sequence" option from the "sequence" output format of the Table browser and/or .csv file 
exported by using the "get all fields from selected table" option from the output format field 
to return tables that sort the data.\n 
For L1 repeats only the .fasta format is required.
For HERV repeats the .fasta repeats file, the .csv repeats table and a classification .csv table is needed.
Family table Columns: ID  Repeat_Name  Family  SuperFamily\n
Sequence table Columns: ID   Location   Sequence\n
Family names are based on a classification csv file the is provided by the user."""

from Bio import SeqIO
import pandas as pd
import csv # import the csv module
import re # regex module

def fams_and_seqs(Choice):
    print "Enter the name of the sequence .fasta file: <name>.fasta"
    rmsk_seq = raw_input(">>> ")
    
    if Choice == "HERV":
        print "Enter the name of the .csv table file: <name>.csv"
        rmsk_table = raw_input(">>> ")
        print ""
        print "Classification table format:"
        print ""
        print "Family_Name        Repeat_Names"
        print "HERVG              LTR85, LTR43, HERVG2"
        print "..."
        print ""
        print "Enter the name of the .csv classification table: <name>.csv"
        class_table = raw_input(">>> ")
        print ""



        #The analysis is starting
        #Read the initial csv table
        all_table = pd.read_csv(str(rmsk_table), header = 1)

        #Reading the clasffication file
        class_FILE=open(str(class_table),"r")
        classification=csv.reader(class_FILE)
        classification.next() # skip first line

        famdict = {}
        templist = []
        #Creating the classification dictionary famdict
        for line in classification:
            templist = line[2].split(", ")
            for i in templist:
                famdict[i]=line[1]
                famdict[line[1]]=line[1]

        #Return the columns of the repeatmasker table as lists
        namelist = all_table['repName'].values.tolist()
        superlist = all_table['repFamily'].values.tolist()
        genoName = all_table['genoName'].values.tolist()
        genoStart = all_table['genoStart'].values.tolist()
        genoEnd = all_table['genoEnd'].values.tolist()

        #Name the output lists that will be converted into dataframes
        IDs=[]
        namez=[]
        sfams=[]
        fams=[]
        loc=[]

        #df = pd.dataframe()
        #df[column] = df[column].apply(lambda x: x.replace('-',''))
        #df[family] = df[column].apply(lambda x: famdict.get(x,'unassigned'))

        for i in range(0,len(namelist)):
            fam_name = {}
            unassigned = []
            name_no_int = re.sub('-int', '', namelist[i]) #Bug fixing
            name_no_minus = re.sub('-', '', name_no_int) #Bug fixing
            name_no_end_ = re.sub('_$', '', name_no_minus) #Bug fixing
            IDs.append("ER"+str("%06d" % (i+1,))) #Get the IDs
            namez.append(namelist[i]) #Get the repeatnames
            sfams.append(superlist[i]) #Get the superfamily

            if namelist[i] in famdict.keys(): #Search for families and get the family names
                fams.append(famdict[namelist[i]])
            elif name_no_int in famdict.keys():
                fams.append(famdict[name_no_int])
            elif name_no_minus in famdict.keys():
                fams.append(famdict[name_no_minus])
            elif name_no_end_ in famdict.keys():
                fams.append(famdict[name_no_end_])
            else:
                fams.append("Unassigned")

            #Get the location in format chrXX:NNNNNNNNN-NNNNNNNNNN
            loc.append(str(genoName[i])+":"+str(genoStart[i])+"-"+str(genoEnd[i]))

        #Create and save the output csv table
        print "Please insert the name of the families output file with the .csv"
        as_to_fams = raw_input(">>> ")
        print ""
        erv_assign_to_families = pd.DataFrame({'ID' : IDs, 'Repeat Name' : namez, 'Family' : fams, 'Superfamily' : sfams})
        erv_assign_to_families = erv_assign_to_families[["ID","Repeat Name", "Family", "Superfamily"]]
        erv_assign_to_families.set_index("ID").to_csv(str(as_to_fams))

        #Getting the sequence information:
        #Read the seq file and get description and sequence for each entry:
        description = [seq_record.description for seq_record in SeqIO.parse(str(rmsk_seq), "fasta")]
        sequences = [seq_record.seq for seq_record in SeqIO.parse(str(rmsk_seq), "fasta")]

        #Output lists:
        location = []
        seqs = []

        #Create the output lists based on the seq file:
        for i in range(0,len(description)):
            desc_space = description[i].split(" ")
            chromosome = desc_space[1][6:].split(":")[0]
            region1 = desc_space[1][6:].split(":")[1].split("-")[0]
            region2 = desc_space[1][6:].split(":")[1].split("-")[1]
            whole = str(chromosome+":"+region1+"-"+region2) #Location
            location.append(whole) #Location
            seq = str(sequences[i]) #Sequence
            seqs.append(seq)

        #Create and save the output csv table
        print "Please insert the name of the Sequences output file with the .csv"
        loc_seqs = raw_input(">>> ")
        erv_seq = pd.DataFrame({'ID' : IDs, 'Location' : loc, "Sequence 5'->3'" : seqs })
        erv_seq = erv_seq[["ID", "Location", "Sequence 5'->3'"]]
        erv_seq.set_index("ID").to_csv(str(loc_seqs))
    
    elif Choice == "L1":
        #The analysis is starting
        #Getting the sequence information:
        #Read the seq file and get identifier, description and sequence for each entry:
        identifiers = [seq_record.id for seq_record in SeqIO.parse(str(rmsk_seq), "fasta")]
        description = [seq_record.description for seq_record in SeqIO.parse(str(rmsk_seq), "fasta")]
        sequences = [seq_record.seq for seq_record in SeqIO.parse(str(rmsk_seq), "fasta")]

        #The output lists:
        IDs=[]
        namez=[]
        sfams=[]
        seqs=[]
        fams=[]
        loc=[]

        #Create the output lists based on the seq file for each entry:
        for i in range(0,len(description)):
            desc_space = description[i].split(" ")
            chromosome = desc_space[1][6:].split(":")[0]
            region1 = desc_space[1][6:].split(":")[1].split("-")[0]
            region2 = desc_space[1][6:].split(":")[1].split("-")[1]
            loc.append(str(chromosome)+":"+str(region1)+"-"+str(region2)) #Location
            IDs.append("L1"+str("%06d" % (i+1,))) #ID
            namez.append(identifiers[i][10:]) #Name
            sfams.append("L1") #Family
            seqs.append(str(sequences[i])) #Sequence

        #Create and save the output csv table
        print "Please insert the name of the Sequences output file: <name>.csv"
        l1_seqs = raw_input(">>> ") 
        print "Please insert the name of the Families output file: <name>.csv"
        l1_f = raw_input(">>> ") 
        l1_names_fam = pd.DataFrame({'ID' : IDs, 'Repeat Name' : namez, 'Superfamily' : sfams}) 
        l1_seq = pd.DataFrame({'ID' : IDs, 'Location' : loc, "Sequence 5'->3'": seqs})
        l1_names_fam = l1_names_fam[["ID","Repeat Name", "Superfamily"]]
        l1_seq = l1_seq[["ID", "Location", "Sequence 5'->3'"]]
        l1_names_fam.set_index("ID").to_csv(str(l1_f))
        l1_seq.set_index("ID").to_csv(str(l1_seqs))

while True:
    print "Which repeats do you want to generate/update? L1 or HERV"
    wduw = raw_input(">>> ")

    if wduw == "L1":
        fams_and_seqs("L1")
        break
    elif wduw == "HERV":
        fams_and_seqs("HERV")
        
        break
    else: 
        print "No valid Repeats name. Try again"
        continue    

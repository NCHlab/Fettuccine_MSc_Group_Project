from Bio import SearchIO 
import sys


families = {}
ids=[]



with open(sys.argv[1]) as f:
    list_of_files = f.read().splitlines()



for xmlfile in list_of_files:
    blast_qresult = SearchIO.read(xmlfile, 'blast-xml')



    for i in range(len(blast_qresult)):

        hitID = blast_qresult[i].id
        fam= hitID.split("|")[1]
        if fam not in families:
            families[fam]= []
        if hitID not in families[fam]:
            families[fam].append(hitID)



for key in families.keys():
    thefile = open("_".join([key,sys.argv[2]]), 'w')
    for seqID in families[key]:
        thefile.write("%s\n" % seqID)
    thefile.close()

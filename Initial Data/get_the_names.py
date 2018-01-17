#Handle RepeatMasker Data
#Try to get the names. These are NOT the families! But i'm getting closer
#L1 retrotransposons
import csv
FILE=open("L1-retrotransposons_full.csv","r")
l1data=csv.reader(FILE)
l1data.next()
l1_names=[]
for line in l1data:
    if str(line[10]) not in l1_names:
        l1_names.append(str(line[10]))
    else:
        continue
#HERVs
FILE2=open("ERV-retrotransposons_full.csv","r")
ervdata=csv.reader(FILE2)
ervdata.next()
ervdata.next()
erv_names=[]
for line in ervdata:
    if str(line[10]) not in erv_names:
        erv_names.append(str(line[10]))
    else:
        continue  

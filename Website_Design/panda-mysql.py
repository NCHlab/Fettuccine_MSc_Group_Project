import MySQLdb
import pandas as pd
import numpy as py

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

df = pd.read_csv("data2.csv")


# Get the cursor, which is used to traverse the database, line by line
#cursor = database.cursor()
cur = connection.cursor()
# Create the INSERT INTO sql query
query = """INSERT INTO prelim_data (genoName, genoStart,
        genoEnd, genoLeft, strand, repName,
        repClass, repFamily, repStart, repEnd, repLeft) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"""

for i in range(0, len(df.index)):
        genoName = df.iloc[i,0]
        genoStart = df.iloc[i,1]
        genoEnd	= df.iloc[i,2]
        genoLeft = df.iloc[i,3]
        strand	= df.iloc[i,4]
        repName	= df.iloc[i,5]
        repClass = df.iloc[i,6]
        repFamily = df.iloc[i,7]
        repStart = df.iloc[i,8]
        repEnd = df.iloc[i,9]
        repLeft	= df.iloc[i,-1]

		# Assign values from each row
        values = (genoName, genoStart, genoEnd, genoLeft, strand, repName, repClass, repFamily, repStart, repEnd, repLeft)

		# Execute sql Query
        cur.execute(query, values)
        connection.commit()
#tables2 = cur.fetchall()


#cur.execute(query, values)

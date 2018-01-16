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

book = pd.read_csv("data2.csv")


# Get the cursor, which is used to traverse the database, line by line
#cursor = database.cursor()
cur = connection.cursor()
# Create the INSERT INTO sql query
query = """INSERT INTO fake_data (genoName, genoStart,
        genoEnd, genoLeft, strand, repName,
        repClass, repFamily, repStart, repEnd, repLeft) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)"""

for i in range(0, len(book.index)):
        genoName = book.iloc[i,0]
        genoStart = book.iloc[i,1]
        genoEnd	= book.iloc[i,2]
        genoLeft = book.iloc[i,3]
        strand	= book.iloc[i,4]
        repName	= book.iloc[i,5]
        repClass = book.iloc[i,6]
        repFamily = book.iloc[i,7]
        repStart = book.iloc[i,8]
        repEnd = book.iloc[i,9]
        repLeft	= book.iloc[i,-1]

		# Assign values from each row
        values = (genoName, genoStart, genoEnd, genoLeft, strand, repName, repClass, repFamily, repStart, repEnd, repLeft)

		# Execute sql Query
        cur.execute(query, values)
        connection.commit()
#tables2 = cur.fetchall()


#cur.execute(query, values)

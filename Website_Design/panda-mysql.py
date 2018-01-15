import pandas as pd
import numpy as np
import pymysql
import pymysql.cursors


# Connect to the database
try:
	connection = pymysql.connect(host="localhost",
						port=3306,
						user="root",
						password="root",
						db="genome_data",#genome_data
						charset="utf8mb4",
						cursorclass=pymysql.cursors.DictCursor)
except:
	print "Failed to connect to db!!"
	print "Please ensure you have started your MySQL server and the db name is correct"
	print "and that the port is set to 3306 in your server"
	pass

book = pd.read_csv("data1.csv", index_col=0)
print book.head()

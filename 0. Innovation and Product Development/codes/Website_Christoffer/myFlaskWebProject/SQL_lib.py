import mysql.connector
from mysql.connector import errorcode
import time

#deprecated
#def get_cleaning_list(cleaning_id):
     #query = "SELECT * FROM information_schema.tables WHERE `TABLE_NAME` like 'Cleaning%'"
     #return query

def set_summary_entry(cnx,cursor,summary_id,table_id,user,procedure,pipe_id):
    timestamp=time.strftime('%Y-%m-%d %H:%M:%S')
    query=("INSERT INTO `"+ summary_id +
        "`SET `cleaning_id` = '"+ table_id +"', " +
        "`status` = 'progress', " +
        "`user` = '"+ user +"', " +
        "`datetime` = '"+ timestamp +"', " +
        "`procedure` = '"+ procedure +"', " +
        "`pipe_id` = '"+ pipe_id + "'")
    print query
    try:
        cursor.execute(query)
        cnx.commit()
        print "Operation successfull."
        
    except mysql.connector.Error:
        print(err.msg)
    return query

def update_table_status(cnx,cursor,summary_id,table_id,status):
    query=("UPDATE `"+ summary_id +
        "`SET `status` = '"+ status +
        "' WHERE `"+ summary_id +"`.`cleaning_id` = '"+ table_id + "'")
    print query
    try:
        cursor.execute(query)
        cnx.commit()
        print "Operation successfull."
        
    except mysql.connector.Error:
        print(err.msg)
    return query

def create_cleaning_table(cleaning_id = "CCCCCC"):
    query = ("CREATE TABLE `"+cleaning_id+"` ("
    "  `TS` TIMESTAMP NOT NULL,"     # Time stamp of cleaning
    "  `Temp` FLOAT(6,3)  NOT NULL,"    # Temperature (usually C) NULLABLE
    "  `PH` FLOAT(6,3)  NOT NULL,"    # PH (usually C) NULLABLE
    "  `Pressure` FLOAT(6,3)  NOT NULL,"    # PH (usually C) NULLABLE
    "  `Conductivity` FLOAT(6,3)  NOT NULL,"    # PH (usually C) NULLABLE
    "  PRIMARY KEY (`TS`)"
    ") ENGINE=InnoDB")   #  FLOAT(3,3)    TIMESTAMP() 
    
    return query

def delele_table(table_id = ""):
    query = ("DROP TABLE `"+table_id+"`"                                     
    )   
    return query

def create_DDBB(databasename = ""):
    query = ("CREATE DATABASE "+databasename+";"                                     
    )   
    return query

def excute_query(query,cursor, extra_text = "" ):
    try:
        print(extra_text)
        cursor.execute(query)
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_TABLE_EXISTS_ERROR:
            print("already exists.")
        else:
            print(err.msg)
    else:
        print("OK")

def add_cleaning_data(cleaning_id, data):
    """
    This funcrion parses the data from the sensors into a SQL query.
    """
    query = "INSERT INTO `"+cleaning_id+ "` ( TS, Temp, PH,Pressure,Conductivity) VALUES "
    
    for i in range(len(data)):
        query += "('" + data["TS"][i].strftime('%Y-%m-%d %H:%m:%S') + "','"+ "%.2f"%data["Temp"][i] + "','"+ "%.2f"%data["PH"][i] + "','"+ \
                 "%.2f"%data["Pressure"][i]+ "','"+ "%.2f"%data["Conductivity"][i]+"')"
    
        if(i < len(data)-1):
            query+=","
        else:
             query+=";"
    return query

def get_cleaning_data(cleaning_id):
     query = "SELECT * FROM  `"+cleaning_id+ "`"
     return query

def set_table_summary_entry(cnx,cursor,CleaningTable):
    query=("UPDATE `"+table_id+
        "`SET `availability` = '0', " +
        "`status` = 'empty', " +
        "`owner` = '', " +
        "`assigment-ts` = '0000-00-00 00:00:00', " +
        "`recipient-last-name` = '', `recipient-first-name` = '', " +
        "`recipient-street` = '', `recipient-street-no` = '', " +
        "`recipient-postcode` = '', `recipient-city` = '', `recipient-country` = '', " +
        "`shipment-no` = '' " +
        "WHERE `"+table_id+"`.`box-id` = "+str(box_id))
    try:
        cursor.execute(query)
        cnx.commit()
        print "Operation successfull."
        
    except mysql.connector.Error:
        print(err.msg)  
 

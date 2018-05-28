import mysql.connector
from mysql.connector import errorcode

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

def add_cleanning_data(cleaning_id, data):
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

def get_cleanning_data(cleaning_id):
     query = "SELECT * FROM  `"+cleaning_id+ "`"
     return query
 
TABLES = {}
TABLES['employees'] = (
    "CREATE TABLE `employees` ("
    "  `emp_no` int(11) NOT NULL AUTO_INCREMENT,"
    "  `birth_date` date NOT NULL,"
    "  `first_name` varchar(14) NOT NULL,"
    "  `last_name` varchar(16) NOT NULL,"
    "  `gender` enum('M','F') NOT NULL,"
    "  `hire_date` date NOT NULL,"
    "  PRIMARY KEY (`emp_no`)"
    ") ENGINE=InnoDB")

TABLES['departments'] = (
    "CREATE TABLE `departments` ("
    "  `dept_no` char(4) NOT NULL,"
    "  `dept_name` varchar(40) NOT NULL,"
    "  PRIMARY KEY (`dept_no`), UNIQUE KEY `dept_name` (`dept_name`)"
    ") ENGINE=InnoDB")

TABLES['salaries'] = (
    "CREATE TABLE `salaries` ("
    "  `emp_no` int(11) NOT NULL,"
    "  `salary` int(11) NOT NULL,"
    "  `from_date` date NOT NULL,"
    "  `to_date` date NOT NULL,"
    "  PRIMARY KEY (`emp_no`,`from_date`), KEY `emp_no` (`emp_no`),"
    "  CONSTRAINT `salaries_ibfk_1` FOREIGN KEY (`emp_no`) "
    "     REFERENCES `employees` (`emp_no`) ON DELETE CASCADE"
    ") ENGINE=InnoDB")

TABLES['dept_emp'] = (
    "CREATE TABLE `dept_emp` ("
    "  `emp_no` int(11) NOT NULL,"
    "  `dept_no` char(4) NOT NULL,"
    "  `from_date` date NOT NULL,"
    "  `to_date` date NOT NULL,"
    "  PRIMARY KEY (`emp_no`,`dept_no`), KEY `emp_no` (`emp_no`),"
    "  KEY `dept_no` (`dept_no`),"
    "  CONSTRAINT `dept_emp_ibfk_1` FOREIGN KEY (`emp_no`) "
    "     REFERENCES `employees` (`emp_no`) ON DELETE CASCADE,"
    "  CONSTRAINT `dept_emp_ibfk_2` FOREIGN KEY (`dept_no`) "
    "     REFERENCES `departments` (`dept_no`) ON DELETE CASCADE"
    ") ENGINE=InnoDB")

TABLES['dept_manager'] = (
    "  CREATE TABLE `dept_manager` ("
    "  `dept_no` char(4) NOT NULL,"
    "  `emp_no` int(11) NOT NULL,"
    "  `from_date` date NOT NULL,"
    "  `to_date` date NOT NULL,"
    "  PRIMARY KEY (`emp_no`,`dept_no`),"
    "  KEY `emp_no` (`emp_no`),"
    "  KEY `dept_no` (`dept_no`),"
    "  CONSTRAINT `dept_manager_ibfk_1` FOREIGN KEY (`emp_no`) "
    "     REFERENCES `employees` (`emp_no`) ON DELETE CASCADE,"
    "  CONSTRAINT `dept_manager_ibfk_2` FOREIGN KEY (`dept_no`) "
    "     REFERENCES `departments` (`dept_no`) ON DELETE CASCADE"
    ") ENGINE=InnoDB")

TABLES['titles'] = (
    "CREATE TABLE `titles` ("
    "  `emp_no` int(11) NOT NULL,"
    "  `title` varchar(50) NOT NULL,"
    "  `from_date` date NOT NULL,"
    "  `to_date` date DEFAULT NULL,"
    "  PRIMARY KEY (`emp_no`,`title`,`from_date`), KEY `emp_no` (`emp_no`),"
    "  CONSTRAINT `titles_ibfk_1` FOREIGN KEY (`emp_no`)"
    "     REFERENCES `employees` (`emp_no`) ON DELETE CASCADE"
    ") ENGINE=InnoDB")


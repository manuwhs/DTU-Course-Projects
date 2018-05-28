"""
Routes and views for the flask application.
"""

from forms import CleanTable
import config_mysql
import SQL_lib
# import config_cosmos
import pydocumentdb.document_client as document_client
import mysql.connector
from mysql.connector import errorcode
from datetime import datetime
from flask import render_template
from myFlaskWebProject import app
from flask import Flask, redirect, url_for, request
import numpy as np


def get_data_from_cleaning_id(cleaning_id):
    # Open database connection
    cnx = mysql.connector.connect(user=config_mysql.DB_USER, password=config_mysql.DB_PASSWORD,
                                  host=config_mysql.DB_HOST, port = config_mysql.DB_PORT,
                                  database=config_mysql.DB_NAME)
    # prepare a cursor object using cursor() method
    cursor = cnx.cursor()
    
 ### Get the data:
    query = SQL_lib.get_cleanning_data(cleaning_id)
    SQL_lib.excute_query(query,cursor, extra_text = " Getting data" )
    data = cursor.fetchall()
#    print data
    row=data[0]
    # Get data from colums as list objects
    data_table=np.array(data)
    help_B=np.asmatrix(data_table)
    time_list=help_B[:,0]
    temp_list=help_B[:,1]
    ph_list=help_B[:,2]
    pressure_list=help_B[:,3]
    conduc_list=help_B[:,4]

    data_str = "["
    Npoints, Nvar = help_B.shape
    for i in range(Npoints):
        time = help_B[i,0]
        press = help_B[i,3]
        ph_i = help_B[i,2]
        data_str = data_str + "[Date.UTC(%i,%i,%i,%i,%i,%i), %.2f],"%(time.year,time.month,time.day,time.hour,time.minute,time.second, ph_i)
    data_str = data_str + "]"
    # disconnect from server
    cnx.close()

    result = CleanTable(data_str,temp_list,ph_list,pressure_list,conduc_list)
    result.cleaning_id = cleaning_id
    return result

@app.route('/')
@app.route('/home')
def home():
    cleaning_id = "hello"
    result = get_data_from_cleaning_id(cleaning_id);
    
    print "Informaiton from DDBB fetched"

    return render_template(
        'results_chart.html',
        title='Some charties',
        year=datetime.now().year,
        result = result
    )



#################### Controller Func get Data ##############
@app.route('/get_cleaning_process_by_id',methods = ['POST', 'GET'])
def get_cleaning_process_by_id():
   if request.method == 'POST':
      cleaning_id = request.form['cleaning_id']
      print ("Form cleaning id", cleaning_id)
      result = get_data_from_cleaning_id(cleaning_id);
      return render_template(
            'results_chart.html',
            title='Some charties',
            year=datetime.now().year,
            result = result
        )

   else:
      cleaning_id = request.form['cleaning_id']
      result = get_data_from_cleaning_id(cleaning_id);
      return render_template(
            'results_chart.html',
            title='Some charties',
            year=datetime.now().year,
            result = result
        )
    

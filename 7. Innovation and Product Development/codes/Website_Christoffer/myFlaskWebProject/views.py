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
from datetime import date, datetime
from flask import render_template, make_response
from flask import Flask,jsonify,json
from myFlaskWebProject import app
from random import random
from time import time
import numpy as np


@app.route('/')
@app.route('/home')
def home():

    # Open database connection
    cnx = mysql.connector.connect(user=config_mysql.DB_USER, password=config_mysql.DB_PASSWORD,
                                  host=config_mysql.DB_HOST, port = config_mysql.DB_PORT,
                                  database=config_mysql.DB_NAME)


    # prepare a cursor object using cursor() method
    cursor = cnx.cursor()
    query = SQL_lib.get_cleaning_data("cleaning_summary")
    SQL_lib.excute_query(query,cursor, extra_text = "Getting cleaning summary table" )
    data = cursor.fetchall()
    #print data[0]
    # disconnect from server
    cnx.close()

    # Post-data-processing
    #data_table=np.array(data)
    #help=np.asmatrix(data_table)
    #list=[]
    #list=help[:,0]

    #json_data=[]
   # for result in data:
    #    json_data.append(dict(zip(row_headers,result)))
    #print json.dumps(json_data, default=obj_dict)

    sum=[]
    for row in data:
        print row
        id = row[1]
        datetime = row[2]
        user = row[3]
        status = row[4]
        table = CleanTable(row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8])
        sum.append(table)
    
    jsonStr = json.dumps([e.toJSON() for e in sum])
    print repr(jsonStr)
    #json_string=json.dumps(sum)
    #print json_string

    #return render_template('results.html', title='Home', result=result)
    return render_template(
        'summary.html',
        title='Home',
        year=datetime.now().year,
        sum=jsonStr
    )

    #client = document_client.DocumentClient(config_cosmos.COSMOSDB_HOST, {'masterKey': config_cosmos.COSMOSDB_KEY})
    ## Read databases and take first since id should not be duplicated.
    #db = next((data for data in client.ReadDatabases() if data['id'] == config_cosmos.COSMOSDB_DATABASE))
    ## Read collections and take first since id should not be duplicated.
    #coll = next((coll for coll in client.ReadCollections(db['_self']) if coll['id'] == config_cosmos.COSMOSDB_COLLECTION))
    ## Read documents and take first since id should not be duplicated.
    #doc = next((doc for doc in client.ReadDocuments(coll['_self']) if doc['id'] == config_cosmos.COSMOSDB_DOCUMENT))
    ## Create a model to pass to results.html
    #class VoteObject:
    #    choices = dict()
    #    total_votes = 0
    #vote_object = VoteObject()
    #vote_object.choices = {
    #    "Web Site" : doc['Web Site'],
    #    "Cloud Service" : doc['Cloud Service'],
    #    "Virtual Machine" : doc['Virtual Machine']
    #}
    #vote_object.total_votes = sum(vote_object.choices.values())
    ## Hack to avoid zero detection in empty database
    #if vote_object.total_votes == 0:
    #    vote_object.total_votes = 1


    #return render_template(
    #    'results.html',
    #    title='Some charties',
    #    year=datetime.now().year,
    #    result = result
    #)

#@app.route('/live-data')
#def live_data():
#    ## Open database connection
#    #cnx = mysql.connector.connect(user=config_mysql.DB_USER, password=config_mysql.DB_PASSWORD,
#    #                              host=config_mysql.DB_HOST, port = config_mysql.DB_PORT,
#    #                              database=config_mysql.DB_NAME)


#    ## prepare a cursor object using cursor() method
#    #cursor = cnx.cursor()
#    #query = SQL_lib.get_cleaning_data("hello")
#    #SQL_lib.excute_query(query,cursor, extra_text = "Getting cleaning summary table" )
#    #data = cursor.fetchall()
#    ##print data[0]
#    ## disconnect from server
#    #cnx.close()

#    #help=np.array(data)
#    #data=np.asmatrix(help)
 

#    # Create a PHP array and echo it as JSON
#    #time=datetime.now()
#    #temp=random() * 100
#    #print data

#    #json=time.strftime("%Y-%m-%d %H:%M:%S")
#    #time = strftime("Date.UTC(%Y,%m,%d,%H,%M,%S")
#    #temp_json=json.dumps(json)
#    #chart_data = [[datetime.now(), random() * 100],[datetime.now(), random() * 100]]

#    chart_data = [datetime.now().isoformat(),random() * 100]
#    ptint chart_data
#    response = make_response(json.dumps(chart_data))
#    print json.dumps(chart_data)
#    response.content_type = 'application/json'
#    return chart_data

@app.route('/live-data')
def live_data():
    # Create a PHP array and echo it as JSON
    # data = [time(), random() * 100],[time(), random() * 100]]
    data = [[time(), random() * 100]]
    response = make_response(json.dumps(data))
    response.content_type = 'application/json'
    return response

@app.route('/test')
def test():
    # Open database connection
    cnx = mysql.connector.connect(user=config_mysql.DB_USER, password=config_mysql.DB_PASSWORD,
                                  host=config_mysql.DB_HOST, port = config_mysql.DB_PORT,
                                  database=config_mysql.DB_NAME)


    # prepare a cursor object using cursor() method
    cursor = cnx.cursor()
    # Create a PHP array and echo it as JSON
    # data = [time(), random() * 100],[time(), random() * 100]]
    table_id="hallo2"
    SQL_lib.update_table_status(cnx,cursor,"cleaning_summary",table_id,"fail")


@app.route('/result/<cleaning_id>')
def result(cleaning_id, chartID = 'chart_ID', chart_type = 'line', chart_height = 500):

    print cleaning_id

    # Open database connection
    cnx = mysql.connector.connect(user=config_mysql.DB_USER, password=config_mysql.DB_PASSWORD,
                                  host=config_mysql.DB_HOST, port = config_mysql.DB_PORT,
                                  database=config_mysql.DB_NAME)


    # prepare a cursor object using cursor() method
    cursor = cnx.cursor()

    # cleaning_id1="hello"
    # Get data of MySQL:
    query = SQL_lib.get_cleaning_data(cleaning_id)
    SQL_lib.excute_query(query,cursor, extra_text = " Getting data" )
    data = cursor.fetchall()
    # Disconnect from server
    cnx.close()

    # Get data from colums as list objects
    help=np.array(data)
    data=np.asmatrix(help)
    #timestamp_list=data[:,0]
    #temp_list=data[:,1]
    #ph_list=data[:,2]
    #pressure_list=data[:,3]
    #conduc_list=data[:,4]

    temp_chart_data = "["
    ph_chart_data = "["
    press_chart_data = "["
    conduc_chart_data = "["
    data_str = "["
    Npoints, Nvar = data.shape
    for i in range(Npoints):
        time = data[i,0]
        temp = data[i,1]
        ph = data[i,2]
        press=data[i,3]
        conduc=data[i,4]
        print time
        #data_str = data_str + "[Date.UTC(%i,%i,%i,%i,%i,%i), %.2f],"%(time.year,(time.month-1),time.day,time.hour,time.minute,time.second, temp)
        temp_chart_data = temp_chart_data + "[Date.UTC(%i,%i,%i,%i,%i,%i), %.2f],"%(time.year,(time.month-1),time.day,time.hour,time.minute,time.second, temp)
        ph_chart_data = ph_chart_data + "[Date.UTC(%i,%i,%i,%i,%i,%i), %.2f],"%(time.year,(time.month-1),time.day,time.hour,time.minute,time.second, ph)
        press_chart_data = press_chart_data + "[Date.UTC(%i,%i,%i,%i,%i,%i), %.2f],"%(time.year,(time.month-1),time.day,time.hour,time.minute,time.second, press)
        conduc_chart_data = conduc_chart_data + "[Date.UTC(%i,%i,%i,%i,%i,%i), %.2f],"%(time.year,(time.month-1),time.day,time.hour,time.minute,time.second, conduc)
        print temp_chart_data
    data_str = data_str + "]"
    temp_chart_data = temp_chart_data + "]"
    ph_chart_data = ph_chart_data + "]"
    press_chart_data = press_chart_data + "]"
    conduc_chart_data = conduc_chart_data + "]"

    return render_template(
        'results.html', 
        title='Result', 
        result=data_str,
        temp_data = temp_chart_data,
        ph_data =  ph_chart_data,
        press_data = press_chart_data,
        conduc_data = conduc_chart_data
    )


@app.route('/contact')
def contact():
    """Renders the contact page."""
    return render_template(
        'contact.html',
        title='Contact',
        year=datetime.now().year,
        message='Your contact page.'
    )

@app.route('/help')
def help():
    """Renders the about page."""
    return render_template(
        'help.html',
        title='Help',
        year=datetime.now().year,
        message='Your application description page.'
    )

@app.route('/profile')
def profile():
    """Renders the about page."""
    return render_template(
        'profile.html',
        title='Profile',
        year=datetime.now().year,
        message='Your profile.'
    )

@app.route('/summary')
def summary():
    """Renders the about page."""
    return render_template(
        'summary.html',
        title='Summary',
        year=datetime.now().year
    )

@app.route('/clear')
def clear():
    #"""Renders the contact page."""
    #client = document_client.DocumentClient(config_cosmos.COSMOSDB_HOST, {'masterKey': config_cosmos.COSMOSDB_KEY})
    ## Attempt to delete the database.  This allows this to be used to recreate as well as create
    #try:
    #    db = next((data for data in client.ReadDatabases() if data['id'] == config_cosmos.COSMOSDB_DATABASE))
    #    client.DeleteDatabase(db['_self'])
    #except:
    #    pass
    ## Create database
    #db = client.CreateDatabase({ 'id': config_cosmos.COSMOSDB_DATABASE })
    ## Create collection
    #collection = client.CreateCollection(db['_self'],{ 'id': config_cosmos.COSMOSDB_COLLECTION })
    ## Create document
    #document = client.CreateDocument(collection['_self'],
    #    { 'id': config_cosmos.COSMOSDB_DOCUMENT,
    #      'Web Site': 0,
    #      'Cloud Service': 0,
    #      'Virtual Machine': 0,
    #      'name': config_cosmos.COSMOSDB_DOCUMENT 
    #    })

    return render_template(
       'clear.html',
        title='Clear',
        year=datetime.now().year
    )

@app.route('/create', methods=['GET', 'POST'])
def create(): 
    if form.validate_on_submit(): # is user submitted vote  
        #client = document_client.DocumentClient(config_cosmos.COSMOSDB_HOST, {'masterKey': config_cosmos.COSMOSDB_KEY})
        ## Read databases and take first since id should not be duplicated.
        #db = next((data for data in client.ReadDatabases() if data['id'] == config_cosmos.COSMOSDB_DATABASE))
        ## Read collections and take first since id should not be duplicated.
        #coll = next((coll for coll in client.ReadCollections(db['_self']) if coll['id'] == config_cosmos.COSMOSDB_COLLECTION))
        ## Read documents and take first since id should not be duplicated.
        #doc = next((doc for doc in client.ReadDocuments(coll['_self']) if doc['id'] == config_cosmos.COSMOSDB_DOCUMENT))
        ## Take the data from the deploy_preference and increment our database
        #doc[form.deploy_preference.data] = doc[form.deploy_preference.data] + 1
        #replaced_document = client.ReplaceDocument(doc['_self'], doc)
        ## Create a model to pass to results.html
        #class VoteObject:
        #    choices = dict()
        #    total_votes = 0
        #vote_object = VoteObject()
        #vote_object.choices = {
        #    "Web Site" : doc['Web Site'],
        #    "Cloud Service" : doc['Cloud Service'],
        #    "Virtual Machine" : doc['Virtual Machine']
        #}
        #vote_object.total_votes = sum(vote_object.choices.values())

        result = CleanTable(0,0,0,0,0)

        return render_template(
            'results.html', 
            year=datetime.now().year,
            title= 'Results',
            result = result
        )

    else :
        return render_template(
            'create.html', 
            title = 'Create',
            year=datetime.now().year
        )

def obj_dict(obj):
    if isinstance(obj, (datetime, date)):
        return obj.isoformat()
    else:
        return obj.__dict__


import datetime
from flask import jsonify

class CleanTable(object):

    # Static Cleaning tabke infos
    id=""
    user=""
    start_time=""
    status=""

    # Nominal values or thresholds
    temp_nominal=0
    phVal_nominal=0
    pressure_nominal=0
    conduc_nominal=0

    # List objects of sampled data
    datetime_list = []
    temp_list = []
    phVal_list = []
    pressure_list = []
    conductivity_list = []
    data_list = []

    # The class "constructor" - It's actually an initializer 
    def __init__(self, id, start_time, user, status, temp_nominal, phVal_nominal, pressure_nominal, conduc_nominal):
        self.id = id
        self.user = user
        self.start_time = start_time
        self.status= status
        self.temp_nominal = temp_nominal
        self.phVal_nominal = phVal_nominal
        self.pressure_nominal = pressure_nominal
        self.conduc_nominal = conduc_nominal
       

    def append_data_list(self, datetime, temp, phVal, pressure, conductivity):
        self.datetime_list.append(datetime)
        self.temp_list.append(temp)
        self.phVal_list.append(phVal)
        self.pressure_list.append(pressure)
        self.conductivity_list.append(conductivity)
        data=[datetime, temp, phVal, pressure, conductivity]
        self.data_list.append(data)
        return self

    def toJSON(self):
        start_time = self.start_time.strftime("%Y-%m-%d %H:%M:%S")
        return {"name":self.id,"status":self.status,"user":self.user,"date":start_time,"temp":self.temp_nominal,"ph":self.phVal_nominal,"con":self.conduc_nominal}
class CleanTable(object):
    time = 0
    temperature = 0
    phVal = 0
    pressure = 0
    conductivity = 0


    # The class "constructor" - It's actually an initializer 
    def __init__(self, time, temperature, phVal, pressure, conductivity):
        self.time = time
        self.temperature = temperature
        self.phVal = phVal
        self.pressure = pressure
        self.conductivity = conductivity

def setTable(time, temperature, phVal, pressure, conductivity):
    table = CleanTable(time, temperature, phVal, pressure, conductivity)
    return table

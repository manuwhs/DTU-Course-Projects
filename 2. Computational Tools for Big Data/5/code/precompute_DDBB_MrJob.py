from mrjob.job import MRJob
from mrjob.step import MRStep
import utilfunc as uf
import highfunc as hf
from subprocess import call
import os


## This file is meant to be called from the main program 
# in order to fucking ppreprocess all the DDBB obtaining their
## BoW and storing it into different files.

rootFolder = "./books"

class MRWordVC(MRJob):

    def mapper_VC(self, _, line):
        # This mapper obtains all the paths of all the files in the DDBB
        # and then fucking fucking yields a task for all of them.
        print os.getcwd()
        outFolder = "./csvs"
#        outFolder = "/home/montoya/Desktop/DTU Lec/1st Semester/5. Computational Tools for Big Data/2.Assingments/Final/csvs"

        filedir = line
        BoW = hf.preprocess_file(filedir, typefile = "HTML", MaxWords = 1000)
        file_name = filedir.split("/")[-1]
        filedir = filedir.split(file_name)[0]
#        print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
        uf.create_folder_if_needed(outFolder)
#        print "###################################"
        
        if (outFolder[0] == "."): # If we are given a local path
            outFolder = "/".join(filedir.split("/")[:-3]) + "/csvs"
            print outFolder
        BoW.to_csv(outFolder + "/" + file_name + ".csv", 
                   encoding = "utf-8", 
                   index = False)  # Do not write the index number

#        print "PENEPENEPNEPNEPNEPNEPNE"
        yield "hola", 1
        
    def combiner_VC(self, key, values):
        # For each document we read the file properly and compute their BoW
        
        yield key, values

    def reducer_VC(self, key, values):
        # This bastard could load all the BoWs into a single one for example
        yield key, values
            
    def steps(self):
        return [
            MRStep(mapper=self.mapper_VC)
#                   combiner=self.combiner_VC,
#                   reducer=self.reducer_VC)
        ]
        
if __name__ == '__main__':
    MRWordVC.run()

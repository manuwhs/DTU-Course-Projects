
import numpy as np
import utilfunc as uf
import textprocesslib as tplib
import highfunc as hf
from subprocess import call  # Execute commands
import pandas as pd 

from scipy.sparse import csr_matrix
import pickle_lib as pkl

from graph_lib import gl
#rootFolder = "./GutenbergData"
rootFolder = "./GutenbergData/A"
outFolder = "./books"
csvFolder = "./csvs"

## Read the first wrong library
read_preproc = 0

if (read_preproc == 1):
    # If we need to read and preprocess the raw files
    # and put then in the output folder
    hf.read_and_preprocess(rootFolder, outFolder)

compute_BoW = 0
if (compute_BoW == 1):
    ## Read the preprocessed data, preprocess the text and metadata
    # Creating a bag of word for every file and pickleing it
    # The files are already in the desired HTML format or any
    # other that we modify later
    
    use_MrJob = 1
    if(use_MrJob == 0):
        hf.compute_and_save_BoW(outFolder, csvFolder, MaxWords = 1000)
    ## USE MRJOB INSTEAD !!!

    elif (use_MrJob == 1):
        allPaths = uf.get_allPaths(outFolder)
        pd.DataFrame(allPaths).to_csv("./paths.txt")
    
        call("python ./precompute_DDBB_MrJob.py ./paths.txt", shell=True)
 
load_BoWs = 0
if (load_BoWs == 1):
    # Load all the BoWs and join them into a common format ?
    All_bow, dict_files = hf.load_and_join_BoWs(csvFolder, MaxFiles = 500)

#####################################################################
################### Query Formatting ###########################
################################################################
# Now we can enter a Query, then it will be transformed in the same way as the
# corpus and a similarity measure between it and all the documents will be done
# The system will output the N most similar elements.

query_f = 1
if(query_f == 1):
    query_dirs = ["./books/0/Winchester","Yorkshire Battles", "1342-0.txt", "11-0.txt"]
    query_dir = query_dirs[2]
#    BoW = hf.preprocess_file("https://docs.python.org/2/howto/urllib2.html", 
#                             typefile = "URL",MaxWords = 1000)

    BoW = hf.preprocess_file("./free_ebook.pdf", typefile = "pdf",MaxWords = 1000)
#    BoW = hf.preprocess_file("1342-0.txt", typefile = "text",MaxWords = 1000)
#    BoW = hf.preprocess_file("11-0.txt", typefile = "text",MaxWords = 1000)
    
    BoW = BoW.set_index(['word'])  # Set the index
    BoW.index = BoW.index.str.encode('utf-8') # Changing
    
    eu_dist = hf.get_doc_sims(All_bow, BoW)

    best_common,best_common_i  = uf.sort_and_get_order(eu_dist[:-1].tolist(), 
                                                       reverse = False)
#    print best_common_i[:5]
#    print dict_files
#    plt.close("all")
    gl.set_subplots(1,3)
    for i in range(3):
        BoW_values = All_bow[str(best_common_i[i])].values
        BoW_index = All_bow.index
        Nwords = 30
        b_v, b_vi = uf.sort_and_get_order(BoW_values, reverse = True)
        gl.bar(BoW_index[b_vi[0:Nwords]].tolist(), 
               BoW_values[b_vi[0:Nwords]],
               legend = [dict_files[best_common_i[i]].split("/")[-1]],
                labels = ["", "", "Num"],
                nf = 1)
            
    ## Get the BoW of the most common:
plotting_thins = 0
if(plotting_thins == 1):
    query_dirs = ["./books/0/Winchester","Yorkshire Battles"]
    query_dirs = ["1342-0.txt", "11-0.txt"]
    labels = [ "Pride and Prejudice", "Alice's Adventures in Wonderland", "Trading Book"]    
    
    gl.set_subplots(3,1)
    for i in range(3):

        
        # HTML  text
    #    BoW = hf.preprocess_file(query_dir, typefile = "text",MaxWords = 1000)
    #    BoW = hf.preprocess_file("https://docs.python.org/2/howto/urllib2.html", 
   #     typefile = "URL",MaxWords = 1000)
        if(i < 2):
            query_dir = query_dirs[i]
            BoW = hf.preprocess_file(query_dir, typefile = "text",MaxWords = 1000)
        else:
            BoW = hf.preprocess_file("./free_ebook.pdf", typefile = "pdf",MaxWords = 1000)
        
        BoW = BoW.set_index(['word'])  # Set the index
        BoW.index = BoW.index.str.encode('utf-8') # Changing
        
#        BoW = hf.load_Bow()
        
        Nwords = 30
        BoW = BoW.sort(['num'], ascending=[0])
        print BoW.index[0:Nwords].shape
        print BoW["num"].values[1:Nwords]
        
        
        gl.bar(BoW.index[0:Nwords].tolist(), BoW["num"].values[0:Nwords],
               legend = [labels[i]],
                labels = ["", "", "Num"],
                nf = 1)
    
    caca  = open("./BoW.txt","w+")
    for i in range(len(BoW.index.tolist())):
        for j in range(BoW.iloc[i]["num"]):
            caca.write(BoW.iloc[i].name + " ")
    caca.close()
    
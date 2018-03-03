
import numpy as np
import utilfunc as uf
import textprocesslib as tplib

rootFolder = "./WEB BEST"
rootFolder = "./sex"
outFolder = "./processed"

allPaths = uf.get_allPaths(rootFolder)

MaxFiles = 1000
i = 0
for filedir in allPaths:
    
    print uf.type_file(filedir) 
    # For all files, we check the extension
    if (uf.type_file(filedir) == "HTML document"):
        # Copy the file
        print filedir
        file_type = uf.process_HTML_doc(filedir)
        ## If the HTML file is not a book
        if (type(file_type) == type(None)):
            continue  # Exit this loop execution 
        
        document = file_type
#        processed = tplib.doc_preprocess(document)
#        print processed[:100]
        
        #Copy the file in the new destination
        uf.copy_file(filedir,outFolder)
        
        i = i + 1
        
        if(i >= MaxFiles):
            break
        



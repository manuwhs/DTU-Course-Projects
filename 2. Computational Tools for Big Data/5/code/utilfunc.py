
import os
import magic
import shutil
from bs4 import BeautifulSoup  # For HTML web treatment.
import json
import numpy as np
###############################################
################ Files library ################
###############################################

def create_folder_if_needed (folder):
    # This function creates a path if it does not exist
    if not os.path.exists(folder):
        os.makedirs(folder)
    
def loadJsonFromFile(filename):
    try:
        with open(filename) as data_file:
            jsonData = json.load(data_file)
        # return created json object
        return jsonData

    except IOError:
        print "Error: File does not appear to exist."
        return 0

###############################################
################ HTML treating ################
###############################################

### AUTHOR AND TITLE
#<div class="dochead">
#<h2 class="author">Lawrence Fletcher</h2>
#<h2 class="title">"Zero the Slaver"</h2>
#<hr/>
#</div>

### CHAPTER ARE LIKE
#<div class="bodytext">
#<a href="" name="chap02"></a>
#<h3>Chapter Two.</h3>
#<h4>A Night of Horror.</h4>


def check_cover(filedir):
    # This function checks if the file has the words cover
    # If it does, it is just the data of it, 
    # I couldnt find a proper way to separate them really
    # They do not follow an easy format
    
    filename = filedir.split("/")[-1]

    if (filename.find("(cover)") != -1):
        return 1
    else:
        return 0
        
def process_HTML_doc(filedir):
    ## Read the file with HTML and close it
    fd = open(filedir, 'r') 
    doc_HTML = fd.read()
    fd.close()
    # Use BeutifulSoup to process the HTML, get tittle, plain text...
    soup = BeautifulSoup(doc_HTML)  # Transform plain text HTML into soup structure
    
    ## Estructure of the HTMLS !!

    ## First check that it is a file of a document
    # <meta content="Zero the Slaver" name="Description"/>
    # <meta name="description" content="Project Gutenberg Ebooks." />
#    content_type = soup.find("meta", {"name": "Description"})
#    if(content_type == None):
#        content_type = soup.find("meta", {"name": "description"})
#    
#    if(content_type == None):
#        return None
#    content_type = content_type["content"]
#    if (content_type == "Project Gutenberg Ebooks."):
#        return None

    useful_text= ""
   #### FIND ALL CHAPTER ####

    useful = soup.findAll("p")  # Directori (Dessert > Bannaan )
    for elem in useful:
        useful_text += " " + elem.text
#        keywords += " " + elem.text
        
#    useful = soup.findAll("div", {"class": "bodytext"})  # Directori (Dessert > Bannaan )
#    for elem in useful:
#        useful_text += " " + elem.text
##        keywords += " " + elem.text
             

    return useful_text

def sort_and_get_order (x, reverse = True ):
    # Sorts x in increasing order and also returns the ordered index
    x = np.array(x)
    x = x.flatten()  # Just in case we are given a matrix vector.
    order = range(len(x))
    
    if (reverse == True):
        x = -x
        
    x_ordered, order = zip(*sorted(zip(x, order)))
    
    if (reverse == True):
        x_ordered = -np.array(x_ordered)
        
    return np.array(x_ordered), np.array(order)
    


import os
import magic
import shutil
from bs4 import BeautifulSoup  # For HTML web treatment.
import json
###############################################
################ Files library ################
###############################################

def create_folder_if_needed (folder):
    # This function creates a path if it does not exist
    if not os.path.exists(folder):
        os.makedirs(folder)
        
def get_allPaths(rootFolder):
    ## This function finds all the files in a folder
    ## and its subfolders

    allPaths = []
    for dirName, subdirList, fileList in os.walk(rootFolder):  # FOR EVERY DOCUMENT
        for fname in fileList:
            # Read the file
            path = dirName + '/' + fname;
            allPaths.append(path)
    
    return allPaths

def type_file(filedir):

    mime = magic.Magic()
    filetype = mime.id_filename(filedir)
#    print filetype
    
    # This will be of the kind "image/jpeg" so "type/format"
    filetype = filetype.split(",")[0]
    return filetype

def copy_file(file_source, file_destination, new_name = ""):
    # Copies a file into a new destination.
    # If a name is given, it changes its name

    file_name = "" 
    file_path = ""
    
    file_name = file_source.split("/")[-1]
    file_path = file_source.split("/")[0]
    
    if (len(new_name) == 0): # No new name specified
        file_name = file_source.split("/")[-1]
    else:
        file_name = new_name
    
    create_folder_if_needed(file_destination)
    
    shutil.copy2(file_source, file_destination + "/" + file_name)

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
    content_type = soup.find("meta", {"name": "Description"})
    if(content_type == None):
        content_type = soup.find("meta", {"name": "description"})
        
    content_type = content_type["content"]
    if (content_type == "Project Gutenberg Ebooks."):
        return None

    useful_text= ""
   #### FIND ALL CHAPTER ####

    useful = soup.findAll("div", {"class": "bodytext"})  # Directori (Dessert > Bannaan )
    for elem in useful:
        useful_text += " " + elem.text
#        keywords += " " + elem.text
             

    return useful_text

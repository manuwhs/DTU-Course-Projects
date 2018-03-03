import unicodedata
from nltk.tokenize import word_tokenize, wordpunct_tokenize, sent_tokenize
from nltk.stem import SnowballStemmer
from nltk.corpus import stopwords

## TEXT PREPROCESSING FUNCTIONS 

def doc_preprocess(document, mode = 0):
    # Preprocesees a document:  Tokenization, lemmatization...
    if (mode == 1):print document
    document = strip_accents(document)
    if (mode == 1):print document
    document = doc_tokeniz(document,mode) 
    if (mode == 1):print document
    document = doc_lowercase(document,mode)
    if (mode == 1):print document
    document = doc_rem_punctuation(document,mode) 
    if (mode == 1):print document
    document = doc_rem_stopwords(document,mode) 
    if (mode == 1):print document
    document = doc_stem(document,mode) 
    if (mode == 1):print document
    return document

def strip_accents(s):
   return ''.join(c for c in unicodedata.normalize('NFD', s)
                  if unicodedata.category(c) != 'Mn')

def doc_tokeniz(document, mode):
    tokens = word_tokenize(document) 
    return tokens
    
def doc_lowercase (document, mode):
    low_text = [w.lower() for w in document] 
    return low_text
    
def doc_rem_stopwords(document, mode):
    stopwords_en = stopwords.words('english')
    clean_text = [word for word in document if not word in stopwords_en]
    return clean_text
    
def doc_stem(document, mode):
    stemmer = SnowballStemmer('english')
    steammed_text = [stemmer.stem(word)for word in document]
    return steammed_text
    
def doc_rem_punctuation(document, mode):
    clean_text = [w for w in document if w.isalnum()]
    return clean_text

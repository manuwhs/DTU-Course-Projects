import unicodedata
from nltk.tokenize import word_tokenize, wordpunct_tokenize, sent_tokenize
from nltk.stem import SnowballStemmer
from nltk.corpus import stopwords

import numpy as np
import pandas as pd
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
#   s = s.decode('utf-8').encode('ascii', 'replace')
#   s = s.encode('utf-8', 'replace')
   
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

#####################################################
############ Bag of Words ##########################

from sklearn.feature_extraction.text import CountVectorizer
from sklearn.feature_extraction.text import TfidfVectorizer
# Initialize the "CountVectorizer" object, which is scikit-learn's
# bag of words tool.  


def obtain_BoW(document, MaxWords = 1000):
    # Obtains the BoW and transforms it to a pandas Dataframe
    # Document is a list of words
    vectorizer = CountVectorizer(analyzer = "word",   \
                                 tokenizer = None,    \
                                 preprocessor = None, \
                                 stop_words = None,   \
                                 max_features = MaxWords) 
    
    BoWdoc = vectorizer.fit_transform(document)   
#    print BoWdoc

    # BoW =  Pairs that tell you for each word, 
    # the index file it belongs to
#   ## Create the dictionary with the words
#    vectorizer.fit(document)
#    # Transform the file to Sparse form matrix
#    BoWdoc = vectorizer.transform(document)
    
    # Get the bocabulary
    vocab = vectorizer.get_feature_names()
#    print vocab
    # Transform it to an array form
    BoWdoc_array = BoWdoc.toarray()
    # Sum up the counts of each vocabulary word
    count = np.sum(BoWdoc_array, axis=0)
#    print count
    BoW = dict()
    BoW["word"] = vocab
    BoW["num"] = count
    
    datFr = pd.DataFrame(BoW)
    
    return datFr
    
def obtain_Tfidf(document):
    # Document is a list of words
    vectorizer = TfidfVectorizer(analyzer = "word",   \
                                 tokenizer = None,    \
                                 preprocessor = None, \
                                 stop_words = None,   \
                                 max_features = 50,  \
                                 use_idf= False) 
                                 
    ## Create the dictionary with the words
    vectorizer.fit(document)
    
    # Get the bocabulary
    vocab = vectorizer.get_feature_names()
    print vocab
    
    # Transform the file to Sparse form matrix
    BoWdoc = vectorizer.transform(document)
    # Pairs that tell you for each word, the index file it belongs to
    
    # Transform it to an array form
    BoWdoc_array = BoWdoc.toarray()
    # Sum up the counts of each vocabulary word
    dist = np.sum(BoWdoc_array, axis=0)

    print BoWdoc
    print dist


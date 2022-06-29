import pandas as pd
import re

def word_tokenize(sentence_list):
    s=sentence_list.lower() # lowercase
    s=pd.Series(re.findall('[^\s]+',s)) # get all sequences of non-spaces in a pd.Series
    s=s.apply(lambda x:re.sub('(^[^\w]+|[^\w]+$)','',x)) # remove non-alphanumeric characters at start and end of all items
    c=s.apply(lambda x:bool(re.search('\w',x))) # find items with alphanumeric characters
    s=s[c] # remove items without alphanumeric characters
    return list(s) # return result as list of strings
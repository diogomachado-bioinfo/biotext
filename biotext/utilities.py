#!/usr/bin/python
# -*- coding: utf-8 -*-
import codecs

def import_txt(input_file_name, strip=False):
    if strip:
        func=lambda x:x.strip()
    else:
        func=lambda x:x
    with codecs.open(input_file_name,'r','utf-8') as f:
        string_list = [func(i) for i in f]
    return string_list
#!/usr/bin/python
# -*- coding: utf-8 -*-
import tempfile
import os
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import numpy as np
from scipy import stats
import codecs
from sweep import fas2sweep
import subprocess
    
def run_clustalo(input_file_name, args=""):
    fp = tempfile.TemporaryFile(mode='w',delete=False)
    cmd_string="clustalo -i "+input_file_name+" -o "+fp.name+" --auto --outfmt clu --force"
    cmd_string+=" "+args
    subprocess.call(cmd_string, shell=True)
    alignment = AlignIO.read(fp.name, "clustal")
    fp.close()
    os.unlink(fp.name)
    return alignment

def get_consensus(seqrecord_list):
    seq_list = list(seqrecord_list)
    fastaText = tempfile.TemporaryFile(mode='w',delete=False)
    for i in seq_list:
        fastaText.write('>'+str(i.description)+'\n'+str(i.seq)+'\n')
    fastaText.close()
    
    alignment = []
    if len(seq_list) > 1:
        align1 = run_clustalo(fastaText.name)
        align2 = []
        for i in align1:
            align2.append(list(i.seq))
            alignment.append(str(i.seq))
        align2 = np.array(align2)
        m = stats.mode(align2) # determine mode
        m = m[0][m[1]>=0] # filter characters by minimal occurrence
        consensus = re.sub('\-+','',''.join(m))
    else:
        consensus = str(seq_list[0].seq)
        alignment.append(consensus)
    
    os.unlink(fastaText.name)
    return consensus, alignment
    
def get_header(seqrecord_list):
    """
    Get the header from all items in a list of SeqRecord (Biopython object).
    
    Parameters
    ----------
    seqrecord_list : list of SeqRecord
        List of SeqRecord.

    Output
    --------
    header_list : list of string
        List of all headers extracted from input.
        
    Example
    --------
    Creates seqrecord_list, extract headers and print it.

    >>> import biotext as bt
    >>> seq_list = ['ACTG','GTCA']
    >>> seqrecord_list = bt.fastatools.create_seqrecord_list(seq_list)
    >>> extracted_header_list = bt.fastatools.get_header(seqrecord_list)
    >>> print(extracted_header_list)
    ['0', '1']
    """
    
    seqrecord_list = list(seqrecord_list)
    header_list = []
    for i in seqrecord_list:
        header = i.description
        header_list.append (header)
    return header_list

def get_seq(seqrecord_list):
    """
    Get the sequences from all items in a list of SeqRecord (Biopython object).
    
    Parameters
    ----------
    seqrecord_list : list of SeqRecord
        List of SeqRecord.

    Output
    --------
    seq_list : list of string
        List of all sequences extracted from input.
        
    Example
    --------
    Creates seqrecord_list, extract sequences and print it.

    >>> import biotext as bt
    >>> seq_list = ['ACTG','GTCA']
    >>> seqrecord_list = bt.fastatools.create_seqrecord_list(seq_list)
    >>> extracted_seq_list = bt.fastatools.get_seq(seqrecord_list)
    >>> print(extracted_seq_list)
    ['ACTG', 'GTCA']
    """
    
    seqrecord_list = list(seqrecord_list)
    seq_list = []
    for i in seqrecord_list:
        seq = i.seq
        seq_list.append (str(seq))
    return seq_list

def create_seqrecord_list(seq_list,header_list=None):
    """
    Creates a list of SeqRecord (Biopython object) from a string list.
    
    Parameters
    ----------
    seq_list : list of string
        List of biological sequences in string format.
    header : list of string
        List of headers in string format, if set to 'None' the headers will 
        be automatically defined with numbers in increasing order.

    Output
    --------
    seqrecord_list : list of SeqRecord
        List of SeqRecord.
        
    Example
    --------
    Decode a string.

    >>> import biotext as bt
    >>> seq_list = ['ACTG','GTCA']
    >>> seqrecord_list = bt.fastatools.create_seqrecord_list(seq_list)
    >>> for i in seqrecord_list:
    >>>     print (i)
    ID: 0
    Name: <unknown name>
    Description: 0
    Number of features: 0
    Seq('ACTG')
    ID: 1
    Name: <unknown name>
    Description: 1
    Number of features: 0
    Seq('GTCA')
    """
    if header_list == None:
        header_list = list(range(0,len(seq_list)))
    seqrecord_list = []
    for i in range(0,len(seq_list)):
        record = SeqRecord(Seq(seq_list[i]), description=str(header_list[i]), id=str(i))
        seqrecord_list.append(record)
    return seqrecord_list

def remove_pattern(seq_list,rex):
    seq_list = list(seq_list)
    for i in range(0,len(seq_list)):
        for ii in rex:
            s = re.sub(ii,'',str(seq_list[i].seq)) # find and remove rex
            seq_list[i].seq=Seq(s)
    return seq_list
    
def import_fasta(input_file_name):
    """
    Uses biopython to import a FASTA file.
    
    Parameters
    ----------
    input_file_name : string (valid file name)
        Input fasta file name.

    Output
    --------
    seqrecord_list : list of SeqRecord
        List of SeqRecord imported from file.
        
    Example
    --------
    Creates a FASTA file named 'test.fasta' and import it as a SeqRecord list.

    >>> import biotext as bt
    >>> seq_list = ['ACTG','GTCA']
    >>> seqrecord_list_1 = bt.fastatools.create_seqrecord_list(seq_list)
    >>> bt.fastatools.export_fasta(seqrecord_list,'test.fasta')
    >>> seqrecord_list_2 = bt.fastatools.create_seqrecord_list(seq_list)
    """
    
    with codecs.open(input_file_name,'r','utf-8') as f:
        seqrecord_list = list(SeqIO.parse(f, "fasta"))
    return seqrecord_list
    
def export_fasta(seqrecord_list, output_file_name, header=None):
    """
    Creates a file using a SeqRecord (Biopython object) list.
    
    Parameters
    ----------
    seqrecord_list : list of SeqRecord
        List of SeqRecord.
    output_file_name : string
        Output fasta file name.
        
    Example
    --------
    Export a SeqRecord list as FASTA file named 'test.fasta'.

    >>> import biotext as bt
    >>> seq_list = ['ACTG','GTCA']
    >>> seqrecord_list = bt.fastatools.create_seqrecord_list(seq_list)
    >>> bt.fastatools.export_fasta(seqrecord_list,'test.fasta')
    """
    
    seqrecord_list = list(seqrecord_list)
    
    if header != None: 
        seqrecord_list = create_seqrecord_list([str(i.seq) for i in seqrecord_list],header=header)
    
    outputFile = codecs.open(output_file_name,'w','utf-8')
    for i in seqrecord_list:
        if len(i.seq) > 0:
            outputFile.write('>'+i.description+'\n')
            seq = str(i.seq)
            seq = re.findall('[\w-]{0,'+str(100)+'}',seq)
            seq = '\n'.join(seq)
            outputFile.write(seq)
    outputFile.close()

def fasta_to_mat(seq_list,orth_mat=None,chunk_size=2E+3,mask=[2,1,2],
                 fasta_type='AA',n_jobs=1):
    try:
        seq_list[0].seq
    except AttributeError:
        seq_list=create_seqrecord_list(seq_list)
    fasta_aux=[]
    min_size=[]
    for i in seq_list:
        if len(str(i.seq)) >= sum(mask):
            fasta_aux.append(i)
            min_size.append(True)
        else:
            min_size.append(False)

    mat=np.zeros((len(seq_list),600))
    mat_aux = fas2sweep(fasta_aux,orth_mat=orth_mat,chunk_size=chunk_size,
                        mask=mask,fasta_type=fasta_type,n_jobs=n_jobs)
    mat[min_size] = mat_aux

    return mat
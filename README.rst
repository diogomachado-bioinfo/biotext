Biotext
=======
The biotext library offers resources to support text mining strategy using bioinformatics tool.

Installation
------------
To install aminocode through `pip`::

      pip install biotext


Tested Platforms
----------------
- Python:

 - 3.7.4

- Windows (64bits):

 - 10

- Ubuntu (64bits)

 - 18.04.1 LTS

Required external libraries
---------------------------
- numpy
- pandas
- scipy
- scikit-learn
- matplotlib
- unidecode
- biopython
- sweep

Functions
---------------
.. csv-table::
   :header: "Function Name", "Description", "Input", "Output"
   :stub-columns: 1

   "biotext.aminocode.encode_string", "Encodes a string with AMINOcode.","input_string \: string:
    Natural language text string to be encoded.
   detail \: string:
    Set details in coding. 'd' for details in digits; 'p' for details on the punctuation; 'dp' or 'pd' for both.", "encoded_string \: string
    Encoded text."
   "biotext.aminocode.decode_string", "Decodes a string with AMINOcode reverse.", "input_string \: string:
    Text string encoded with AMINOcode.
   detail \: string:
    Set details in coding. 'd' for details in digits; 'p' for details on the punctuation; 'dp' or 'pd' for both.", "decoded_string \: string
    Decoded text."
   "biotext.aminocode.encode_list ", "Encodes all strings in a list with AMINOcode.", "string_list \: list
    List of string to be encoded.
   detail \: string
     Set details in coding. 'd' for details in digits; 'p' for details on the punctuation; 'dp' or 'pd' for both.
   verbose \: bool
    If True displays progress.", "encoded_list \: list
    List with all encoded text in string format."
   "biotext.aminocode.decode_list", "Decodes all strings in a list with reverse AMINOcode.", "string_list \: list
    List of string encoded with aminocode.
   detail \: string
    Set details in coding. 'd' for details in digits; 'p' for details on the punctuation; 'dp' or 'pd' for both.
   verbose \: bool
    If True displays progress.", "decoded_list \: list of string
     List with all decoded text."
   "biotext.dnabits.encode_string", "Encodes a string with DNAbits.","input_string \: string:
    Natural language text string to be encoded.", "encoded_string \: string
    Encoded text."
   "biotext.dnabits.decode_string", "Decodes a string with DNAbits reverse.", "input_string \: string:
    Text string encoded with AMINOcode.", "decoded_string \: string
    Decoded text."
   "biotext.dnabits.encode_list ", "Encodes all strings in a list with DNAbits.", "string_list \: list
    List of string to be encoded.
   verbose \: bool
    If True displays progress.", "encoded_list \: list
    List with all encoded text in string format."
   "biotext.dnabits.decode_list", "Decodes all strings in a list with reverse DNAbits.", "string_list \: list
    List of string encoded with aminocode.
   verbose \: bool
    If True displays progress.", "decoded_list \: list of string
     List with all decoded text."
   "create_seqrecord_list", "Creates a list of SeqRecord*SeqRecord* from a string list.", "seq_list \: list of string
    List of biological sequences in string format.
   header \: list of string
    List of headers in string format, if set to 'None' the headers will be automatically defined with numbers in increasing order.", "seqrecord_list \: list of SeqRecord*
    List of SeqRecord*."
   "biotext.fastatools.import_fasta", "Uses biopython to import a FASTA file.", "input_file_name \: string (valid file name)
    Input fasta file name.", "seqrecord_list \: list of SeqRecord*
    List of SeqRecord* imported from file."
   "biotext.fastatools.export_fasta", "Creates a file using a SeqRecord*SeqRecord* list.", "seqrecord_list \: list of SeqRecord*
    List of SeqRecord*.
   output_file_name \: string
    Output fasta file name.", "A file is saved with the defined name."
   "biotext.fastatools.get_header", "Get the header from all items in a list of SeqRecord*SeqRecord*.", "seqrecord_list \: list of SeqRecord*
    List of SeqRecord*.", "header_list \: list of string
    List of all headers extracted from input."
   "biotext.fastatools.get_seq", "Get the sequences from all items in a list of SeqRecord*SeqRecord*.", "seqrecord_list \: list of SeqRecord*
    List of SeqRecord*.", "seq_list \: list of string
    List of all sequences extracted from input."
   "biotext.fastatools.remove_pattern", "Removes patterns from a SeqRecord* range based on regular expression.", "seq_list \: list of SeqRecord*
    List of SeqRecord*.
   rex \: string
    regular expression.", "seq_list \: list of SeqRecord*
    List of SeqRecord* with removal applied."
   "biotext.fastatools.run_clustalo", "Uses the Clustal Omega to align the strings in a FASTA file.", "input_file_name \: string (valid file name)
    Input fasta file name.", "alignment \: MultipleSeqAlignment**
    Alignment result."
   "biotext.fastatools.get_consensus", "Applies clustalo and obtain alignment consensus.", "seqrecord_list \: list of SeqRecord*
    List of SeqRecord*.", "consensus \: string
    Alignment consensus.
   alignment \: list of string
    List of sequences with alignment gaps."
   "biotext.fastatools.fasta_to_mat", "Performs a vectorization of a list of SeqRecord* using the SWeeP.", "seq_list \: list of string
    List of strings in FASTA format.", "mat \: ndarray***
    Matrix with the generated vectors."
   "biotext.treetools.mat_to_tree", "Creates a dendrogram in newick format from a matrix.", "mat \: ndarray***
    Matrix.
   ids \: list of string
    List with line identifiers in mat.
   method \: string
    Method to create the dendrogram. Available options are 'complete', scipy library implementation, and 'nj' (neighbor joining), skbio library implementation. The default is the 'complete' method.", "tree \: string
    tree: dendrogram in newick format."

\*Bio.SeqRecord.SeqRecord: Biopython object to store biological sequences and its information, as described in <https://biopython.org/docs/1.76/api/Bio.SeqRecord.html>.

\*\*Bio.Align.MultipleSeqAlignment: Biopython object to store biological multiple sequence alignment, as described in <https://biopython.org/docs/1.76/api/Bio.Align.html>.

\*\*\*numpy.ndarray: Numpy object to represent array, as described in <https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html>.

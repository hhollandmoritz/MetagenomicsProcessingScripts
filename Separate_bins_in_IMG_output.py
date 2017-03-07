#!/usr/bin/env python

__author__ = "Hannah Holland-Moritz"
__email__ = "hannah.hollandmoritz@colorado.edu"
__version__ = "0.0.1"

"""
Parses IMG output files for bin names and separates into separate files. This 
script relies on the SuffixTree package from 
https://github.com/JDonner/SuffixTree to quickly search for strings in long
lists of contigs. And on the SeqIO package from biopython for fasta parsing.

To install SuffixTree with sudo:
git clone https://github.com/JDonner/SuffixTree.git
cd SuffixTree
sudo python setup.py build
sudo python setup.py install

To install SuffixTree without sudo:
git clone https://github.com/JDonner/SuffixTree.git
cd SuffixTree
python setup.py build
python setup.py install --user

Then update $PYTHONPATH to include the directory that contains local packages!
(Usually this looks something like adding these lines to your .bash_profile:
PYTHONPATH=/path/to/home/dir/.local/lib/python2.7/site-packages/:$PYTHONPATH
export PYTHONPATH)
"""

import argparse
import os
from Bio import SeqIO
from collections import defaultdict
import SuffixTree.SubstringDict


def main():
    parser = argparse.ArgumentParser(description=\
        'Parses IMG output files for bin names and separates into separate \
        files. Make sure the python modules SuffixTree and SeqIO are \
        installed before running.')
    req = parser.add_argument_group('required arguments')
    req.add_argument('-i', '--input_fp', required=True,
        help='A file with the file paths of each annotation \
        file you want to parse. One file path per line.')
    req.add_argument('-b', '--bins_list', required=True,
        help='A path to a file with the list of bin names to search for in \
        the first column of the input file. Each bin must have its own line')
    req.add_argument('-c', '--input_contigs', required=True,
        type=str, help='A tab-delimited file with your contig names in the \
        first column and IMG-asseigned contigs in the second column. \
        This is usually the img output file with labeled map. No header.')
    parser.add_argument('-o', '--output_dir',
        default='Separated_annotation_information',
        help='The output directory. Will create if does not exist.')


    ## Main Script ##

    args = parser.parse_args()
    
    out_dir = args.output_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    else:
        raise RuntimeError('Output directory already exists. Remove or specify a \
            new directory.')
    out_dir = out_dir + "/"
    
    # Read the input file and get a list of files to separate
    files = parse_input_fp(args.input_fp)
    
    # Read input bins and extract the contigs associated with each bin
    bin_dictionary = parse_bin_contigs(args.input_contigs, args.bins_list)    
    
    # Separate input files into bins
    for current_file in files:
        separate_file_data(current_file, bin_dictionary, out_dir)
    
    # Inform user script has finished
    print("Yay! It's finished!")


## Function Definitions ##


def parse_input_fp(input_fp):
    files_to_sep = open(input_fp, 'U')
    files_list = []

    for line in files_to_sep:
        file_path = line.strip()
        files_list.append(file_path)

    files_to_sep.close()
    
    return files_list
        
def parse_bin_contigs(contig_list_fp, bin_list_fp):
    bin_list = open(bin_list_fp, 'U')
    # check format of contig list
    with open(contig_list_fp,'U') as contig_list:
        if len(contig_list.readline().split("\t")) != 2:
            raise RuntimeError('Contig list does not appear to be in the\n \
            correct format. Please check it and try again.')  
    
    # read in contig list        
    contig_list = open(contig_list_fp,'U')
    contig_dict = SuffixTree.SubstringDict()
    
    for line in contig_list:
        line_split = line.strip().split('\t')
        contig_dict[line_split[0]] = line_split[1]

    # make dictionary of bins and their contigs    
    bin_dict = defaultdict(list)    
    for current_bin in bin_list:
        current_bin = current_bin.strip()
        current_bin_contig_list = contig_dict[current_bin]
        bin_dict[current_bin] = current_bin_contig_list
    
    bin_list.close()
    contig_list.close()
    
    return bin_dict
    
def separate_file_data(file_to_sep_fp, bin_dict, output_dir):
    # Check if file is fasta format or other annotation format
    annotation_file = open(file_to_sep_fp, 'U')
    with open(file_to_sep_fp,'U') as annotation_file:
        if annotation_file.read(1) == ">":
            fasta = True
        else:
            fasta = False
    annotation_file = open(file_to_sep_fp, 'U')
    file_to_sep_bn = os.path.basename(os.path.normpath(file_to_sep_fp))
    print("Separating %s:" % file_to_sep_bn)
    
    # if fasta create new file name and call fasta parser
    if fasta == True:         
         for current_bin in bin_dict:
             current_bin = current_bin.strip()
             if not os.path.exists(output_dir + current_bin):
                 os.makedirs(output_dir + current_bin)
             new_file_name = (output_dir + current_bin + "/" + current_bin + 
             "." + file_to_sep_bn)
             contigs = bin_dict[current_bin]
             count = read_and_write_fasta(file_to_sep_fp, new_file_name, 
                                          contigs)
             print("%i sequences written for %s." % (count, current_bin))


    # else create new file name and call non-fasta parser
    else:
         for current_bin in bin_dict:
             current_bin = current_bin.strip()
             if not os.path.exists(output_dir + current_bin):
                 os.makedirs(output_dir + current_bin)
             new_file_name = (output_dir + current_bin + "/" + current_bin + 
             "." + file_to_sep_bn)
             contigs = bin_dict[current_bin]
             count = read_and_write_non_fasta(file_to_sep_fp, new_file_name,
                                              contigs)
             print("%i sequences written for %s." % (count, current_bin))
             
             
             
    return new_file_name

def read_and_write_fasta(read_fp, write_fp, contig_list):
    # open files and set up sequence file dictionary
    output_file = open(write_fp, 'w')
    input_file = open(read_fp, 'U')
    fasta_dict = SeqIO.index(read_fp, "fasta")

    # Get a list of fasta headers to check
    record_ids = set()
    for record in SeqIO.parse(read_fp, "fasta"):
        record_ids.add(record.id)
        
    # If contig has an associated fasta header, write it to the file
    count = 0
    for i in contig_list:
        if i in record_ids:
            new_fasta = fasta_dict[i].format("fasta")
            output_file.write(new_fasta)            
            count = count + 1            
    output_file.close()
    input_file.close()
    
    return count
    
def read_and_write_non_fasta(read_fp, write_fp, contig_list):
    # open files and set up suffixtree dictionary
    output_file = open(write_fp, 'w')
    input_file = open(read_fp, 'U')
    input_dict = {}
    
    # Create dictionary of input file
    for input_line in input_file:
        # only split at first column
        input_line_split = input_line.strip().split('\t', 1)
        input_dict[input_line_split[0]] = input_line_split[1]
        
        
    # write the contig_list entries to a file
    count = 0
    for i in contig_list:
        for j in input_dict.keys():
            if i in j and len(input_dict[j]) > 0:
                output_line = j + "\t" + input_dict[j]
                output_file.write(output_line)
                output_file.write("\n")
                count = count + 1
    # close files
    output_file.close()
    input_file.close()
    
    return count
    
if __name__ == "__main__":
    main()
'''
SNID_utils.py

A collection of utility functions for interacting with SNID inputs and outputs
'''

# imports
import numpy as np

def it_line_locate(fl, sstring):
    '''
    iteratively finds the line number in a file that contains the first occurrence of a given string

    Parameters
    ----------
    fl : .txt-like file to search for the string in
    sstring : string to search for in the file

    Returns
    -------
    line_num : integer number of the line the the searched for string first occurs on, or None if string not found
    '''

    line_num = 0
    found = False
    with open(fl, 'r') as f:
        for line in f:
            if sstring in line:
                found = True
                break
            else:
                line_num += 1

    return line_num if found else None

def z_arg(z_host):
    '''
    takes SN host redshift and returns SNID formatted command to force redshift
    returns empty string if no use-able (float or int) redshift is passed

    Parameters
    ----------
    z_host : redshift of the SN host - needs to float or int to be counted

    Returns
    -------
    SNID formatted command to force the redshift unless no use-able redshift passed (then empty string)
    '''

    if type(z_host) is float or type(z_host) is int:
        return 'forcez={}'.format(z_host)
    else:
        return ''

def read_output_file(output_file):
    '''
    parses a snid output ('_snid.output') file and returns type results and template rlap results above the rlap cutoff

    Parameters
    ----------
    output_file : snid output file containing results to read (file ending with '_snid.output')

    Returns
    -------
    type_results : matrix containing type results
                   columns - 'type', 'ntemp', 'fraction', 'slope', 'redshift', 'redshift_error', 'age', 'age_error'
    template_results : matrix of template results (above rlap cutoff)
                       columns - 'no.', 'sn', 'type', 'lap', 'rlap', 'z', 'zerr', 'age', 'age_flag', 'grade'
    '''

    # column names for type results section of file
    tp_col_names = ['type', 'ntemp', 'fraction', 'slope', 'redshift', 'redshift_error', 'age', 'age_error']
    tp_res_start = 39 # lines to skip before reading type results

    # max rows is the stopping point minus the starting point minus one empty line
    max_rows_type = it_line_locate(output_file, '### rlap-ordered template listings ###') - tp_res_start - 5
    
    # read in type results section, attempt to handle and report errors if tp_res_start is incorrect
    while True:

        try:
            type_results = np.genfromtxt(output_file, dtype = None, skip_header = tp_res_start,
                                         max_rows = max_rows_type, names= tp_col_names, encoding = 'utf-8')
            break

        except ValueError:
            print('\nWarning: line number error when reading output file. Attempting to rectify...')
            
            # find the correct results starting point and recalculate number of rows to read
            tp_res_start = it_line_locate(output_file, '#' + ' '.join(tp_col_names)) + 1
            max_rows_type = it_line_locate(output_file, '### rlap-ordered template listings ###') - tp_res_start - 5

            print('Starting point identified as line {}'.format(tp_res_start))
            print('Number of lines identified to be {}'.format(max_rows_type))

    # column names for rlap results section of file
    rlap_col_names = ['no.', 'sn', 'type', 'lap', 'rlap', 'z', 'zerr', 'age', 'age_flag', 'grade']
    rlap_res_start = tp_res_start + max_rows_type + 7 # lines to skip before reading rlap results

    # max rows is the stopping point minus the starting point
    max_rows_rlap = it_line_locate(output_file, '#--- rlap cutoff') - rlap_res_start

    template_results = np.genfromtxt(output_file, dtype = None, skip_header = rlap_res_start,
                                     max_rows = max_rows_rlap, names = rlap_col_names, encoding = 'utf-8')

    return type_results, template_results

def read_lnw(fname):
    '''
    reads a .lnw file as created by the logwave routine provided with SNID and returns its contents

    Parameters
    ----------
    fname : filename of .lnw file containing processed spectral epoch(s) of a SN

    Returns
    -------
    sn_name : name of SN
    sn_type : SN type
    ages : list of floating point age(s) for the spectral epoch(s)
    cols : list of column names (including ages)
    f : matrix of data with column names described by cols
    '''
    
    # open file and read meta data
    with open(fname, 'r') as f:
        
        # extra number of spectra, SN name, and SN type from the first line
        header = f.readline()
        header_list = [item for item in header.split(' ') if item != '']
        spec_num = int(header_list[0])
        sn_name = header_list[5]
        sn_type = header_list[7]
        
        # read lines and store contents until the line containing the ages (phases) is reached
        next_line = f.readline()
        next_line_list = [item for item in next_line.split(' ') if item != '']
        line_count = 2
        while len(next_line_list) != 1 + spec_num:    
            line_count += 1
            next_line = f.readline()
            next_line_list = [item for item in next_line.split(' ') if item != '']
            
    # extract and process ages for use as floats and separately for use in indexing  
    ages = []
    cols = ['wav']
    for age in next_line_list[1:]:
        ages.append(float(age))
        tmp = age.rstrip()
        tmp = ''.join(tmp.split('-'))
        cols.append(''.join(tmp.split('.')))
    
    # read the remainder of the file    
    f = np.genfromtxt(fname, names = cols, dtype = None, skip_header = line_count)

    # return the column names and the matrix of data that they describe
    return sn_name, sn_type, ages, cols, f
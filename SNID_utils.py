'''
SNID_utils.py

A collection of utility functions for interacting with SNID inputs and outputs
'''

# imports
import numpy as np

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

    # extract the SNID output file name from the spectrum file
    output_file = '{}_snid.output'.format(fname.split('.')[0])

    # column names for type results section of file
    tp_col_names = ['type', 'ntemp', 'fraction', 'slope', 'redshift', 'redshift_error', 'age', 'age_error']
    tp_res_start = 39 # lines to skip before reading type results
    max_rows_type = 29 # number of rows to read for type results
    
    # read in type results section, attempt to handle and report errors if tp_res_start is incorrect
    while True:

        try:
            type_results = np.genfromtxt(output_file, dtype = None, skip_header = tp_res_start,
                                         max_rows = max_rows_type, names= tp_col_names, encoding = 'utf-8')
            break

        except ValueError:
            print('Warning: line {} appears NOT to be the line where type results start'.format(tp_res_start))
            
            # attempt to find correct type results starting point
            tp_res_start = 1
            with open(output_file, 'r') as f:
                for line in f:
                    if '#' + ' '.join(tp_col_names) in line:
                        break
                    else:
                        tp_res_start += 1
            print('Starting point identified as line {}'.format(tp_res_start))

    # column names for rlap results section of file
    rlap_col_names = ['no.', 'sn', 'type', 'lap', 'rlap', 'z', 'zerr', 'age', 'age_flag', 'grade']
    rlap_res_start = tp_res_start + max_rows_type + 7 # lines to skip before reading rlap results

    # iteratively find stopping point for rlap results reading
    ctr = 0
    with open(output_file, 'r') as f:
        for line in f:
            if '#-- rlap cutoff' in line:
                break
            else:
                ctr += 1

    # max rows is the stopping point minus the starting point
    max_rows_rlap = ctr - rlap_res_start

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
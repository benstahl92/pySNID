'''
SNID_utils.py

A collection of utility functions for interacting with SNID inputs and outputs
the casual user will have no need to directly interact with the contents of this file
'''

# imports
import os
import shlex
import subprocess
import numpy as np

def _snid_check():
    '''checks that SNID is properly installed and raises an import error if not'''

    import shutil
    s = shutil.which('snid')
    if (s is None) or (s == ''):
        raise ImportError("SNID is not installed as expected --- make sure 'snid' invokes it")

def _it_line_locate(fl, sstring):
    '''
    iteratively find line number in file that containing first occurrence of given string

    Parameters
    ----------
    fl : .txt-like file to search for string in
    sstring : string to search for file

    Returns
    -------
    line_num : integer number of the line the the searched for string first occurs on, or None if string not found

    Notes
    -----
    could be done in less lines and probably much faster with pandas
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

def _z_arg(z):
    '''
    takes SN redshift and returns SNID formatted string to force redshift
    returns empty string if no use-able (float or int) redshift is passed

    Parameters
    ----------
    z : redshift of SN (or its host) -- needs to be float or int to be counted

    Returns
    -------
    SNID formatted argument to force the redshift unless no use-able redshift passed (then empty string)
    '''

    if ((type(z) == type(0.1)) or (type(z) == (1))) and (not np.isnan(z)):
        return 'forcez={}'.format(z)
    else:
        return ''

def _template_arg(template):
    '''
    takes preferred template type (e.g. Ia or Ia-norm) and returns SNID formatted string to force only that type
    returns empty string if all templates are specified

    Parameters
    ----------
    template : string specifying a valid SNID (sub)type to restrict comparison to

    Returns
    -------
    SNID formatted arguments to force use of specific template type unless all templates requested (then empty string)
    '''

    if template == 'all':
        return ''
    else:
        return 'usetype={}'.format(template)

def _rlap_arg(rlap):
    '''
    takes rlap limit and returns SNID formatted string to enforce it
    return empty string if default is specified

    Parameters
    ----------
    rlap : rlap value to set for limit

    Returns
    -------
    NID formatted string to enforce rlap limit unless default is requested (then empty string)
    '''

    if rlap == 'default':
        return ''
    else:
        return 'rlapmin={}'.format(rlap)


def exec_SNID(fname, z = None, template = 'all', rlap = 'default', z_tol = 0.02, print_cmd = False):
    '''
    formulates and executes SNID command

    Parameters
    ----------
    fname : filename of spectrum file to be analyzed
    z : optional, redshift of the SN (or its host) -- needs to float or int to be counted
    template : optional, string specifying a valid SNID (sub)type to restrict comparison to
    rlap : optional, rlap value to run SNID with
    z_tol : optional, allowable redshift tolerance (SNID default is 0.02) -- not typically changed
    print_cmd : optional, print SNID command

    Returns
    -------
    output_file : SNID output file generated by executed run (None if run fails)

    Effects
    -------
    runs SNID on spectrum in fname, subject to optional arguments
    '''

    # formulate arguments
    forcez_arg = _z_arg(z)
    ztol_arg = 'zfilter={}'.format(z_tol)
    rlap_arg = _rlap_arg(rlap)
    template_arg = _template_arg(template)

    # construct command
    snid_command = 'snid {} {} {} {} plot=0 inter=0 verbose=0 {}'.format(rlap_arg, forcez_arg, ztol_arg, template_arg, fname)

    # get output file name (made in working directory) and delete file from previous runs
    output_file = '{}_snid.output'.format(os.path.basename(fname).split('.')[0])
    if os.path.isfile(output_file):
        os.remove(output_file)

    # execute command
    if print_cmd:
        print(snid_command)
    p = subprocess.Popen(shlex.split(snid_command), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = p.communicate() # for now, don't do anything with output

    # return output file name if SNID run is successful
    if os.path.isfile(output_file):
        return output_file
    else:
        None

def read_output_file(output_file):
    '''
    parses a SNID output ('_snid.output') file and returns type results and template rlap results above the rlap cutoff

    Parameters
    ----------
    output_file : SNID output file containing results to read (file ending with '_snid.output')

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
    max_rows_type = _it_line_locate(output_file, '### rlap-ordered template listings ###') - tp_res_start - 5
    
    # read in type results section, attempt to handle and report errors if tp_res_start is incorrect
    while True:

        try:
            type_results = np.genfromtxt(output_file, dtype = None, skip_header = tp_res_start,
                                         max_rows = max_rows_type, names= tp_col_names, encoding = 'utf-8')
            break

        except ValueError:
            print('\nWarning: line number error when reading output file. Attempting to rectify...')
            
            # find the correct results starting point and recalculate number of rows to read
            tp_res_start = _it_line_locate(output_file, '#' + ' '.join(tp_col_names)) + 1
            max_rows_type = _it_line_locate(output_file, '### rlap-ordered template listings ###') - tp_res_start - 5

            print('Starting point identified as line {}'.format(tp_res_start))
            print('Number of lines identified to be {}'.format(max_rows_type))

    # column names for rlap results section of file
    rlap_col_names = ['no.', 'sn', 'type', 'lap', 'rlap', 'z', 'zerr', 'age', 'age_flag', 'grade']
    rlap_res_start = tp_res_start + max_rows_type + 7 # lines to skip before reading rlap results

    # max rows is the stopping point minus the starting point
    max_rows_rlap = _it_line_locate(output_file, '#--- rlap cutoff') - rlap_res_start

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

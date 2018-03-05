'''
run_SNID.py

Primary Function: BSNID

    Runs the SuperNova IDentification code (SNID) [Blondin and Tonry 2007] with the updated 
    templates and prescription presented in BSNIP I [Silverman 2012] to find the type, subtype,
    redshift and error, age and error of a SN from its spectrum
'''

# imports - standard
import os
import numpy as np

# imports - custom
from SNID_utils import z_arg, read_output_file

# helper functions

def SNID_type(fname, path, z_host, rlap, z_tol):
    '''
    runs SNID on spectrum with appropriate args, reads SNID output file, determines SN type
    SN type is determined if type of best matching ('good') template is same as type with
        highest fraction of 'good' matches (must be over 50 percent)

    Parameters
    ----------
    fname : filename of spectrum file to be analyzed
    path : full path to spectrum file
    z_host : redshift of the SN host - needs to float or int to be counted
    rlap : rlap value to run SNID with
    z_tol : allowable redshift tolerance (SNID default is 0.02) - not typically changed

    Returns
    -------
    SN type if determined, None if type is not able to be determined
    '''

    # formulate optional arguments to pass to SNID
    forcez_arg = z_arg(z_host)
    ztol_arg = 'zfilter={}'.format(z_tol)
    rlap_arg = 'rlapmin={}'.format(rlap)

    # combine arguments to formulate snid command and execute it
    snid_command = 'snid {} {} {} plot=0 inter=0 verbose=0 {}'.format(rlap_arg, forcez_arg, ztol_arg, path + fname)
    os.system(snid_command)

    # extract the SNID output file name from the spectrum file
    output_file = '{}_snid.output'.format(fname.split('.')[0])

    # if output_file does not exist return None
    if os.path.isfile(output_file) is False:
        return None

    # if output_file does exist, proceed
    else:
        # read output file from run
        type_results, template_results = read_output_file(output_file)

        # get favored type as type which has highest fraction from results
        # note: this should be robust against detecting subtypes b/c of way snid tabulates results
        favored_type = type_results['type'][type_results['fraction'] == np.max(type_results['fraction'])][0]
        fav_tp_frac = type_results['fraction'][type_results['fraction'] == np.max(type_results['fraction'])][0]

        # get type of best template (in terms of rlap) that is 'good' fit (last indexing strips off subtype information)
        best_template_type = template_results[template_results['grade'] == 'good'][0]['type'][:2]

        # if favored type has fraction over 50 percent and is same as type of best fitting template return it - otherwise None
        return favored_type if fav_tp_frac >= 0.5 and favored_type == best_template_type else None

def SNID_subtype(fname, path, z_host, rlap, z_tol, template_type):
    '''
    runs SNID on spectrum with appropriate args, reads SNID output file, determines SN subtype
    SN subtype is determined if subtype of best matching ('good') template is same as subtype with
        highest fraction of 'good' matches (must be over 50 percent)

    Parameters
    ----------
    fname : filename of spectrum file to be analyzed
    path : full path to spectrum file
    z_host : redshift of the SN host - needs to float or int to be counted
    rlap : rlap value to run SNID with
    z_tol : allowable redshift tolerance (SNID default is 0.02) - not typically changed
    template_type : SN type to force SNID to use with search

    Returns
    -------
    SN subtype if determined, None if type is not able to be determined
    '''

    # formulate optional arguments to pass to SNID
    forcez_arg = z_arg(z_host)
    ztol_arg = 'zfilter={}'.format(z_tol)
    rlap_arg = 'rlapmin={}'.format(rlap)
    template_arg = 'usetype={}'.format(template_type)

    # combine arguments to formulate snid command and execute it
    snid_command = 'snid {} {} {} {} plot=0 inter=0 verbose=0 {}'.format(rlap_arg, forcez_arg, ztol_arg, template_arg, path + fname)
    os.system(snid_command)

    # extract the SNID output file name from the spectrum file
    output_file = '{}_snid.output'.format(fname.split('.')[0])

    # if output_file does not exist return None
    if os.path.isfile(output_file) is False:
        return None

    # if output_file does exist, proceed
    else:
        # read output file from run
        type_results, template_results = read_output_file(output_file)

        # create copy of type_results that omits entry for template_type
        # then find the subtype as entry with the highest fraction
        tp_res_masked = type_results[type_results['type'] != template_type]
        favored_subtype = tp_res_masked['type'][tp_res_masked['fraction'] == np.max(tp_res_masked['fraction'])][0]
        fav_stp_frac = tp_res_masked['fraction'][tp_res_masked['fraction'] == np.max(tp_res_masked['fraction'])][0]

        # get subtype of best template (in terms of rlap) that is 'good' fit
        best_template_type = template_results[template_results['grade'] == 'good'][0]['type']

        # if favored subtype has fraction over 50 percent and is same as subtype of best fitting template return it - otherwise None
        return favored_subtype if fav_stp_frac >= 0.5 and favored_subtype == best_template_type else None

def SNID_redshift(fname, path, template_type):
    '''
    runs SNID on spectrum with appropriate args, reads SNID output file, determines SN redshift and error
    SN redshift is median of redshifts from 'good' template matches
    SN redshift error is the standard deviation of redshifts from 'good' template  matches

    Parameters
    ----------
    fname : filename of spectrum file to be analyzed
    path : full path to spectrum file
    template_type : SN type or subtype to force SNID to use with search

    Returns
    -------
    SN redshift
    SN redshift error
    '''

    # formulate optional arguments to pass to SNID
    template_arg = 'usetype={}'.format(template_type)

    # combine arguments to formulate snid command and execute it
    snid_command = 'snid {} plot=0 inter=0 verbose=0 {}'.format(template_arg, path + fname)
    os.system(snid_command)

    # extract the SNID output file name from the spectrum file
    output_file = '{}_snid.output'.format(fname.split('.')[0])

    # if output_file does not exist return None, None
    if os.path.isfile(output_file) is False:
        return (None, None)

    # if output_file does exist, proceed
    else:
        # read output file from run
        type_results, template_results = read_output_file(output_file)

        # extract redshifts from 'good' template matches
        redshifts = template_results[template_results['grade']=='good']['z']

        # return median and std deviation if there is at least one redshift found
        return (np.median(redshifts), np.std(redshifts)) if len(redshifts) > 0 else (None, None)

def SNID_age(fname, path, z, z_tol, template_type, relax_age_restr):
    '''
    runs SNID on spectrum with appropriate args, reads SNID output file, determines SN age and error
    SN age is median of ages from 'good' template matches that have rlap at least 75 pct of maximum rlap
    SN age error is standard deviation of ages
    SN age error req'd to less than 4 days or 20 pct of SN age (whichever greater) for results to be returned

    Parameters
    ----------
    fname : filename of spectrum file to be analyzed
    path : full path to spectrum file
    z : redshift of the SN (either from the host or from SNID)
    rlap : rlap value to run SNID with
    z_tol : allowable redshift tolerance (SNID default is 0.02) - not typically changed
    template_type : SN subtype to force SNID to use with search
    relax_age_restr : bool that selects whether to enforce age restrictions - should be False

    Returns
    -------
    SN age
    SN age error
    '''

    # formulate optional arguments to pass to SNID
    forcez_arg = z_arg(z)
    ztol_arg = 'zfilter={}'.format(z_tol)
    template_arg = 'usetype={}'.format(template_type)

    # combine arguments to formulate snid command and execute it
    snid_command = 'snid {} {} {} plot=0 inter=0 verbose=0 {}'.format(forcez_arg, ztol_arg, template_arg, path + fname)
    os.system(snid_command)

    # extract the SNID output file name from the spectrum file
    output_file = '{}_snid.output'.format(fname.split('.')[0])

    # if output_file does not exist return None, None
    if os.path.isfile(output_file) is False:
        return (None, None)

    # if output_file does exist, proceed
    else:
        # read output file from run
        type_results, template_results = read_output_file(output_file)

        # get the rlap of the best fit template with 'good' grading
        rlap_best = template_results[template_results['grade'] == 'good']['rlap'][0]

        # get ages as the ages from templates that have an rlap at least 75 percent of rlap_best and have 'good' grading
        ages = template_results[np.logical_and(template_results['grade'] == 'good', template_results['rlap'] >= 0.75 * rlap_best)]['age']

        # get the age and error
        age = np.median(ages)
        age_err = np.std(ages)

        # if age restrictions are not relaxed, return age and age error if age error satisfies below reqs, otherwise return None
        if relax_age_restr is False:
            # return age and error if the error is less than 4d or 20 percent of age (whichever greater)
            return (age, age_err) if age_err < 4 or age_err < 0.2 * age else (None, None)

        # if age restrictions are relaxed, then return age and age error regardless
        else:
            return (age, age_err)

def BSNID(data_dict, base_dir, rlaps = (10,5), z_tol = 0.02, relax_age_restr = False):
    '''
    runs the SuperNova IDentification code (SNID) [Blondin and Tonry 2007] with the updated 
    templates and prescription presented in BSNIP I [Silverman 2012] to find the type, subtype,
    redshift and error, age and error of a SN from its spectrum

    Parameters
    ----------
    data_dict : dictionary containing data relevant to the spectrum to analyze
                must have values for keys: 'Filename', 'Filepath', 'Redshift_Gal'
    base_dir : base directory path that all other files and paths are specified relative to
    rlaps : tuple, optional
            tuple containing the SNID rlap values to use, in order of preference
    z_tol : allowable redshift tolerance (SNID default is 0.02) - not typically changed
    relax_age_restr : bool, optional
                      selects whether to enforce age restrictions - should be False

    Returns
    -------
    SN_type : SN type
    SN_subtype : SN subtype
    z_snid : SN redshift
    z_snid_error : SN redshift error
    age : SN age
    age_error : SN age error

    Outputs
    -------
    running this function will lead to a '_snid.output' file being generated (for each spectrum) and left in the working dir
    '''

    # extract needed quantities from data_dict
    fname = data_dict['Filename']
    fpath = data_dict['Filepath']
    z_host = data_dict['Redshift_Gal']

    # construct full path to spectrum file
    path = base_dir + fpath + '/'

    # instantiate variables for results
    SN_type, SN_subtype, z_snid, z_snid_error, age, age_error = None, None, None, None, None, None

    # run SNID_type with higher rlap value
    SN_type = SNID_type(fname, path, z_host, rlaps[0], z_tol)
        
    # if type isn't found, try again with lower rlap value
    if SN_type is None:
        SN_type = SNID_type(fname, path, z_host, rlaps[1], z_tol)
    
    # if type found, try to find subtype (with higher rlap first)
    if SN_type is not None:
        SN_subtype = SNID_subtype(fname, path, z_host, rlaps[0], z_tol, SN_type)

        # if subtype isn't found, try again with lower rlap value
        if SN_subtype is None:
            SN_subtype = SNID_subtype(fname, path, z_host, rlaps[1], z_tol, SN_type)

        # if subtype found (then guaranteed that type is found), try to find redshift
        if SN_subtype is not None:
            z_snid, z_snid_error = SNID_redshift(fname, path, SN_subtype)

            # if snid redshift and subtype found, try to find age
            if z_snid is not None and z_snid_error is not None:

                # set redshift
                if type(z_host) is float or type(z_host) is int:
                    z = z_host
                else:
                    z = z_snid

                age, age_error = SNID_age(fname, path, z, z_tol, SN_subtype, relax_age_restr)

        # if subtype not found try for redshifts
        if SN_subtype is None:
            z_snid, z_snid_error = SNID_redshift(fname, path, SN_type)

    # return calculated quantities
    return SN_type, SN_subtype, z_snid, z_snid_error, age, age_error

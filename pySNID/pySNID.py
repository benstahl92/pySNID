# imports - standard
import os
import numpy as np

# imports - custom
from .pySNIDutils import exec_SNID, read_output_file

def SNID_type(fname, z = None, rlap = 10, z_tol = 0.02, zmin = None, zmax = None):
    '''
    runs SNID on spectrum with appropriate args, reads SNID output file, determines SN type
    SN type is determined if type of best matching ('good') template is same as type with
        highest fraction of 'good' matches (must be over 50 percent)

    Parameters
    ----------
    fname : filename of spectrum file to be analyzed
    z : optional, redshift of SN (or its host) -- needs to be float or int to be counted
    rlap : optional, rlap value to run SNID with
    z_tol : optional, allowable redshift tolerance (SNID default is 0.02) - not typically changed

    Returns
    -------
    SN type, best matching template info, and number of good matches if determined, None if type is not able to be determined
    '''

    # execute SNID and retrieve output
    output_file = exec_SNID(fname, z = z, rlap = rlap, z_tol = z_tol, zmin = zmin, zmax = zmax)

    # if output_file does not exist return None
    if output_file is None:
        return None, None, None

    # if output_file does exist, proceed
    else:
        # read output file from run
        type_results, template_results = read_output_file(output_file)

        # get favored type as type which has highest fraction from results, also get number of 'good' matches
        # note: this should be robust against detecting subtypes b/c of way snid tabulates results
        favored_type = type_results['type'][type_results['fraction'] == np.max(type_results['fraction'])][0]
        fav_tp_frac = type_results['fraction'][type_results['fraction'] == np.max(type_results['fraction'])][0]
        fav_tp_num = type_results['ntemp'][type_results['fraction'] == np.max(type_results['fraction'])][0]

        # retrieve information on the best template (in terms of rlap) that is 'good' fit
        best_template = template_results[template_results['grade'] == 'good'][0]

        # find the type of the best template
        if best_template['type'] not in ['NotSN','AGN','Gal','LBV','M-star','C-Star','QSO']:
            # if not one of the above, then type is two letters long
            best_template_type = best_template['type'][:2]
        else:
            best_template_type = best_template['type']

        # if favored type has fraction over 50 percent and is same as type of best fitting template return it (and info) - otherwise None
        return (favored_type, best_template, fav_tp_num) if ((fav_tp_frac >= 0.5) and (favored_type == best_template_type)) else (None, None, None)

def SNID_subtype(fname, z = None, template_type = 'all', rlap = 10, z_tol = 0.02, zmin = None, zmax = None):
    '''
    runs SNID on spectrum with appropriate args, reads SNID output file, determines SN subtype
    SN subtype is determined if subtype of best matching ('good') template is same as subtype with
        highest fraction of 'good' matches (must be over 50 percent)

    Parameters
    ----------
    fname : filename of spectrum file to be analyzed
    z : optional, redshift of the SN host - needs to float or int to be counted
    template_type : optional, SN type to force SNID to use with search
    rlap : optional, rlap value to run SNID with
    z_tol : allowable redshift tolerance (SNID default is 0.02) - not typically changed

    Returns
    -------
    SN subtype if determined, None if type is not able to be determined
    '''

    # execute SNID and retrieve output
    output_file = exec_SNID(fname, z = z, template = template_type, rlap = rlap, z_tol = z_tol, zmin = zmin, zmax = zmax)

    # if output_file does not exist return None
    if output_file is None:
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
        return favored_subtype if ((fav_stp_frac >= 0.5) and (favored_subtype == best_template_type)) else None

def SNID_redshift(fname, template_type = 'all', rlap = 'default', z_tol = 0.02, zmin = None, zmax = None):
    '''
    runs SNID on spectrum with appropriate args, reads SNID output file, determines SN redshift and error
    SN redshift is median of redshifts from 'good' template matches
    SN redshift error is the standard deviation of redshifts from 'good' template  matches

    Parameters
    ----------
    fname : filename of spectrum file to be analyzed
    template_type : optional, SN type or subtype to force SNID to use with search
    zmin, zmax : redshift bounds

    Returns
    -------
    SN redshift
    SN redshift error
    '''

    # execute SNID and retrieve output
    output_file = exec_SNID(fname, template = template_type, rlap = rlap, z_tol = z_tol, zmin = zmin, zmax = zmax)

    # if output_file does not exist return None
    if output_file is None:
        return (None, None)

    # if output_file does exist, proceed
    else:
        # read output file from run
        type_results, template_results = read_output_file(output_file)

        # extract redshifts from 'good' template matches
        redshifts = template_results[template_results['grade'] == 'good']['z']

        # return median and std deviation if there is at least one redshift found
        return (np.median(redshifts), np.std(redshifts)) if (len(redshifts) > 0) else (None, None)

def SNID_age(fname, z, template_type = 'all', rlap = 'default', relax_age_restr = False, z_tol = 0.02):
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
    fl_ext : file extension (i.e. .flm) of spectrum file

    Returns
    -------
    SN age
    SN age error
    '''

    # execute SNID and retrieve output
    output_file = exec_SNID(fname, z = z, template = template_type, rlap = rlap, z_tol = z_tol)

    # if output_file does not exist return None
    if output_file is None:
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
            return (age, age_err) if ((age_err < 4) or (age_err < np.abs(0.2 * age))) else (None, None)

        # if age restrictions are relaxed, then return age and age error regardless
        else:
            return (age, age_err)

def pySNID(fname, z, rlaps = (10, 5), z_tol = 0.02, relax_age_restr = False, zd_zmin = 0.0, zd_zmax = 0.5):
    '''
    runs the SuperNova IDentification code (SNID) [Blondin and Tonry 2007] with the updated 
    templates and prescription presented in BSNIP I [Silverman 2012] and used by Stahl et al. 2019
    to find the type, subtype, redshift and error, age and error of a SN from its spectrum

    Parameters
    ----------
    fname : filename of spectrum file to be analyzed (including path)
    z : 
    rlaps : tuple, optional
            tuple containing the SNID rlap values to use, in order of preference
    z_tol : allowable redshift tolerance (SNID default is 0.02) - not typically changed
    relax_age_restr : bool, optional
                      selects whether to enforce age restrictions - should be False
    zd_zmin, zd_zmax : bounds to use when determining redshift

    Returns
    -------
    SN_type : SN type
    SN_subtype : SN subtype
    z_snid : SN redshift
    z_snid_error : SN redshift error
    age : SN age
    age_error : SN age error

    Effects
    -------
    running this function will lead to a '_snid.output' file being generated (for each spectrum) and left in the working dir
    '''

    # instantiate variables for results (want default to be None)
    SN_type, best_templ, good_num, SN_subtype, z_snid, z_snid_error, age, age_error = None, None, None, None, None, None, None, None

    # run SNID_type with higher rlap value
    SN_type, best_templ, good_num = SNID_type(fname, z, rlap = rlaps[0], z_tol = z_tol, zmin = zd_zmin, zmax = zd_zmax)
        
    # if type isn't found, try again with lower rlap value
    if SN_type is None:
        SN_type, best_templ, good_num = SNID_type(fname, z, rlap = rlaps[1], z_tol = z_tol, zmin = zd_zmin, zmax = zd_zmax)
    
    # if type found, try to find subtype (with higher rlap first)
    if SN_type is not None:
        SN_subtype = SNID_subtype(fname, z, template_type = SN_type, rlap = rlaps[0], z_tol = z_tol, zmin = zd_zmin, zmax = zd_zmax)

        # if subtype isn't found, try again with lower rlap value
        if SN_subtype is None:
            SN_subtype = SNID_subtype(fname, z, template_type = SN_type, rlap = rlaps[1], z_tol = z_tol, zmin = zd_zmin, zmax = zd_zmax)

        # if subtype found (then guaranteed that type is found), try to find redshift
        if SN_subtype is not None:
            z_snid, z_snid_error = SNID_redshift(fname, template_type = SN_subtype, zmin = zd_zmin, zmax = zd_zmax)

            # if snid redshift and subtype found, try to find age
            if (z_snid is not None) and (z_snid_error is not None):

                # set redshift
                if (type(z) == type(0.1)) or (type(z) == type(1)):
                    z = z
                else:
                    z = z_snid

                age, age_error = SNID_age(fname, z, z_tol = z_tol, template_type = SN_subtype, relax_age_restr = relax_age_restr)

        # if subtype not found try for redshifts
        if SN_subtype is None:
            z_snid, z_snid_error = SNID_redshift(fname, template_type = SN_type, zmin = zd_zmin, zmax = zd_zmax)

    # return calculated quantities
    return SN_type, best_templ, good_num, SN_subtype, z_snid, z_snid_error, age, age_error

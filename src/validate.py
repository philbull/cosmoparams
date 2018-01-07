#!/usr/bin/env python

def compulsory_parameters_exist(param_list):
    """
    Check if compulsory parameters exist in a set of parameter names.
    """
    # List of compulsory parameters
    compulsory = ['Omega_c', 'Omega_b', 'Omega_l', 'h', 'n_s', 
                  ('A_s', 'sigma_8')]
    all_compulsory = True
    missing = []
    
    # Loop over compulsory parameters
    for c in compulsory:
        # Sets of parameters for which at least one must be specified
        if isinstance(c, tuple):
            any_found = False
            for p in c:
                if p in param_list: any_found = True
            if not any_found:
                all_compulsory = False
                missing.append(c)
            continue
            
        # Individual parameters
        if c not in param_list:
            all_compulsory = False
            missing.append(c)
            
    return all_compulsory, missing


def fill_missing_base_params(params):
    """
    Take a parameter dictionary that may have unspecified base parameters and 
    infer/insert them using consistency relations and default values.
    
    Parameters
    ----------
    params : dict
        Dictionary containing cosmological parameters and their values.
    
    Returns
    -------
    params : dict
        Parameter dictionary with all base parameters filled in according to 
        consistency relations or default values.
    """
    p = params.keys()
    
    # Set parameters to default values if they are missing
    if 'w0' not in p: params['w0'] = -1.
    if 'wa' not in p: params['wa'] = 0.
    if 'Omega_n_mass' not in p: params['Omega_n_mass'] = 0.
    if 'Omega_n_rel' not in p: params['Omega_n_rel'] = 0.
    if 'Omega_g' not in p: params['Omega_g'] = 0.
    if 'mnu' not in p: params['mnu'] = 0.
    if 'N_nu_mass' not in p: params['N_nu_mass'] = 0.
    if 'N_nu_rel' not in p: params['N_nu_rel'] = 3.046
    
    # Fill missing normalization parameter with 'None' only if *one* of these 
    # parameters is missing (if both are missing, the parameter set was invalid)
    if ('sigma_8' in p) != ('A_s' in p):
        if 'sigma_8' in p and 'A_s' not in p: params['A_s'] = None
        if 'A_s' in p and 'sigma_8' not in p: params['sigma_8'] = None
    
    # Set parameters using consistency relations if they are missing
    if 'Omega_m' not in p:
        params['Omega_m'] = params['Omega_c'] + params['Omega_b'] \
                          + params['Omega_n_mass']
    if 'Omega_k' not in p:
        params['Omega_k'] = 1. - (  params['Omega_m'] + params['Omega_l'] 
                                  + params['Omega_g'] + params['Omega_n_rel'] )
    
    # Return updated parameter dict
    return params
    

def validate_cosmo_desc(params):
    """
    Validate a set of cosmological parameters (according to Section 2 of the 
    specification).
    
    Parameters
    ----------
    params : dict
        Dictionary containing cosmological parameters and their values.
    
    Returns
    -------
    valid : bool
        Whether the parameter set was valid (True) or invalid (False).
    """
    valid = True
    
    # Check that compulsory parameters exist
    all_compulsory, missing = compulsory_parameters_exist(params.keys())
    if not all_compulsory: valid = False
    
    # Make sure all base parameters are specified before performing more tests
    params = fill_missing_base_params(params)
    
    # Make sure that consistency conditions are satisfied
    if params['Omega_m'] != (  params['Omega_c'] + params['Omega_b'] \
                             + params['Omega_n_mass'] ):
        valid = False
    
    if params['Omega_k'] != 1. - (  params['Omega_m'] + params['Omega_l'] 
                                  + params['Omega_g'] + params['Omega_n_rel'] ):
        valid = False
    
    # Return value
    return valid


def validate_file_desc(filename):
    """
    Load a file and check whether it validates against the interchange format 
    standard in Section 3 of the specification. The validity of the 
    cosmological parameter set is not tested.
    
    Parameters
    ----------
    filename : str
        Path to the file that will be validated.
    
    Returns
    -------
    valid : bool
        Whether the file was valid (True) or invalid (False).
    
    errors : list of str
        List of validation errors, including line numbers of where errors 
        were found.
    """
    # Load file as a sequence of strings
    f = open(filename, 'r')
    lines = f.readlines() # Effectively enforces (3.6)
    f.close()
    
    params = {}
    errors = []
    
    # Loop over lines and check validity
    for i, l in enumerate(lines):
        
        # Check for empty line and skip
        if len(l.strip()) == 0: continue
        
        # Check for tab characters (3.7)
        if "\t" in l:
            errors.append("Line %d: Tab characters are not allowed." % i)
        
        # Check whether disallowed unicode characters were used (3.1)
        try:
            l.decode('ascii')
        except UnicodeDecodeError:
            errors.append("Line %d: Unicode characters are not allowed." % i)
        
        # Strip whitespace and newline from beginning/end of lines (3.7)
        l = l.strip()
        
        # Test whether line is a comment (3.2)
        if l[0] == '#': continue
        
        # Check for number of colons (only one is allowed) (3.3)
        if l.count(":") == 0:
            errors.append("Line %d: Missing separator ':'." % i)
        if l.count(":") > 1:
            errors.append("Line %d: Too many ':' separators." % i)
        
        # Split into parameter : value pair and trim whitespace from the ends
        pname, pval = l.split(":")
        pname = pname.strip()
        pval = pval.strip()
        
        # Check for validity of parameter name
        invalid_chars = '!"#$%&\'()*+,-./;<=>?@[\\]^{|}~'
        for c in invalid_chars:
            if c in pname:
                errors.append("Line %d: "
                              "Invalid character '%s' in parameter name." 
                              % (i, c))
        
        # Add parameter to dictionary
        try:
            params[pname] = float(pval)
        except:
            params[pname] = str(pval)
        
    # Check for compulsory parameters (3.9)
    all_compulsory, missing = compulsory_parameters_exist(params.keys())
    
    # Error message if compulsory parameters not found
    if not all_compulsory:
        errors.append("Line %d: Compulsory parameters missing (%s)." \
                      % (i, missing))
    
    if len(errors) == 0:
        return True, errors
    return False, errors


def read_cosmo_desc(filename):
    """
    Read parameters and their values from a standardized interchange file.
    
    Parameters
    ----------
    filename : str
        Path of the file from which the parameters should be read.
    
    Returns
    -------
    params : dict
        Dictionary of parameters and values read from the file.
    """
    # Load file as a sequence of strings
    f = open(filename, 'r')
    lines = f.readlines() # Effectively enforces (3.6)
    f.close()
    
    params = {}
    # Loop over lines and add to dict
    for i, l in enumerate(lines):
        
        # Check for empty line and skip
        if len(l.strip()) == 0: continue
        
        # Strip whitespace and newline from beginning/end of lines (3.7)
        l = l.strip()
        
        # Test whether line is a comment (3.2)
        if l[0] == '#': continue
        
        # Check for number of colons (only one is allowed) (3.3)
        if l.count(":") != 1:
            raise SyntaxError("More than one ':' found on line %d." % i)
        
        # Split into parameter : value pair and trim whitespace from the ends
        pname, pval = l.split(":")
        pname = pname.strip()
        pval = pval.strip()
        
        # Add sequence parameter to dictionary
        if "," in pval:
            try:
                vals = [float(v) for v in pval.split(",")] # List of floats
            except:
                vals = [v.strip() for v in pval.split(",")] # List of strings
            continue
            
        # Add string or float parameter to dictionary
        try:
            params[pname] = float(pval)
        except:
            params[pname] = str(pval)
    
    # Validate 
    if not validate_cosmo_desc(params):
        raise ValueError("Cosmological parameter set failed validity checks.")
    return params


def write_cosmo_desc(filename, params):
    """
    Write parameters out to a file in a standardized interchange format.
    
    Parameters
    ----------
    filename : str
        Path where the file should be saved.
    
    params : dict
        Dictionary containing parameters and their values.
    """
    # Loop over (sorted) parameters, adding a parameter per line
    lines = []
    pnames = params.keys()
    pnames.sort()
    for p in pnames:
        assert(" " not in p)
        
        # Handle lists/arrays/tuples
        if hasattr(params[p], '__iter__'):
            val = ",".join(["%10.10e" % v for v in params[p]])
            lines.append( "%s:%s\n" % (p, val) )
            continue
            
        # Handle float vs. string parameters
        try:
            lines.append( "%20s : %10.10e\n" % (p, float(params[p])) )
        except:
            lines.append( "%20s : %s\n" % (p, str(params[p])) )
    
    lines = "".join(lines)
    f = open(filename, 'w')
    f.write(lines)
    f.close()


if __name__ == '__main__':
    
    # Filename for example parameter file
    fname = "test.params"
    
    # Example cosmological parameter dictionary
    params = {
        'Omega_c':  0.26,
        'Omega_b':  0.05,
        'Omega_l':  0.69,
        'h':        0.69,
        'sigma_8':  0.8,
        'n_s':      0.96,
    }
    
    # Validate cosmological parameters
    valid = validate_cosmo_desc(params)
    if valid:
        print("Cosmological parameter set passed validity checks.")
    else:
        print("Cosmological parameter set failed validity checks.")
    
    # Write parameters to file
    write_cosmo_desc(fname, params)
    print("Wrote parameters to '%s'." % fname)
    
    # Test that the output file is valid
    valid, errors = validate_file_desc(fname)
    if valid:
        print("File '%s' passed all validity checks." % fname)
    else:
        print("File '%s' failed some validity checks:" % fname)
        for err in errors: print("\t%s" % err)
    
    # Read parameter file back in
    params_in = read_cosmo_desc(fname)
    valid = validate_cosmo_desc(params_in)
    if valid:
        print("Cosmological parameter set passed validity checks.")
    else:
        print("Cosmological parameter set failed validity checks.")
    

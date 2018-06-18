from hashlib import md5
import numpy as np

# Version of the CosmoParams spec that this code conforms to
__SPEC_VERSION = 0.1

class CosmoParamsBase(object):
    
    def __init__(self, Omega_c=None, Omega_b=None, Omega_m=None, Omega_k=None, 
                 Omega_l=None, Omega_n_mass=0., Omega_n_rel=0., Omega_g=0., 
                 w0=-1., wa=0., h=None, N_nu_mass=0., N_nu_rel=3.046, mnu=0., 
                 A_s=None, n_s=None, sigma_8=None):
        """
        Base class for cosmological parameter container classes, based on the 
        CosmoParams specification.
        
        Parameters
        ----------
        Omega_c : float, optional
            Fractional density of CDM relative to critical density.
        Omega_b : float, optional
            Fractional density of baryons relative to critical density.
        Omega_m : float, optional
            Fractional density of all matter relative to critical density.
        Omega_k : float, optional
            Fractional density of curvature relative to critical density.
        Omega_l : float, optional
            Fractional density of Lambda/DE relative to critical density.
        Omega_n_mass : float, optional
            Omega_nu for massive neutrinos.
        Omega_n_rel : float, optional
            Omega_nu for massless neutrinos.
        Omega_g : float, optional
            Fractional density of photons relative to critical density.
        w0 : float, optional
            Dark energy equation of state (value at a=1).
        wa : float, optional 
            Dark energy equation of state (controls asymptotic value as a -> 0).
        h : float, optional
            Dimensionless Hubble parameter (H_0 = 100 h km/s/Mpc).
        N_nu_mass : float, optional
            Number of different species of massive neutrinos.
        N_nu_rel : float, optional
            Effective number of massless neutrinos.
        mnu : float, optional
            Sum of neutrino masses, in eV.
        A_s : float, optional
            Primordial scalar power spectrum amplitude (dimensionless).
        n_s : float, optional
            Spectral index of primordial scalar power spectrum.
        sigma_8 : float, optional
            Power spectrum amplitude, sigma_8.
        """
        # Define the set of parameters managed by this object
        self.parameters = {
          'Omega_c' : 'Fractional density of CDM relative to critical density',
          'Omega_b' : 'Fractional density of baryons relative to critical '
                      'density',
          'Omega_m' : 'Fractional density of all matter relative to critical ' 
                      'density',
          'Omega_k' : 'Fractional density of curvature relative to critical '
                      'density',
          'Omega_l' : 'Fractional density of Lambda/DE relative to critical '
                      'density',
          'Omega_n_mass' : 'Omega_nu for massive neutrinos',
          'Omega_n_rel' : 'Omega_nu for massless neutrinos',
          'Omega_g' : 'Fractional density of photons relative to critical '
                      'density',
          'w0' : 'Dark energy equation of state (value at a=1)',
          'wa' : 'Dark energy equation of state (controls asymptotic value '
                 'as a -> 0)',
          'h' : 'Dimensionless Hubble parameter',
          'N_nu_mass' : 'Number of different species of massive neutrinos',
          'N_nu_rel' : 'Effective number of massless neutrinos',
          'mnu' : 'Sum of neutrino masses, in eV',
          'A_s' : 'Primordial scalar power spectrum amplitude (dimensionless)',
          'n_s' : 'Spectral index of primordial scalar power spectrum',
          'sigma_8' : 'Power spectrum amplitude, sigma_8',
        }
        
        # Initialize parameter dict
        self.base = {}
        
        # Set all parameters to the specified values
        args = locals()
        for p in self.parameters.keys(): self.base[p] = args[p]
        
        # Define which parameters go into each type of hash (see Sect. 4 of the 
        # CosmoParams spec). These should generally be pre-sorted alphabetically.
        self.hash_types = {
            'expansion':    ['Omega_g', 'Omega_k', 'Omega_l', 'Omega_m', 
                             'Omega_n_rel', 'h', 'w0', 'wa'],
            'growth':       ['Omega_g', 'Omega_k', 'Omega_l', 'Omega_m', 
                             'Omega_n_rel', 'h', 'w0', 'wa'],
            'powspec':      ['A_s', 'N_nu_mass', 'Omega_b', 'Omega_c', 
                             'Omega_k', 'Omega_l', 'Omega_n_mass', 'h', 'mnu', 
                             'n_s', 'sigma_8', 'w0', 'wa']
            'lensing':      [],
            'cmb':          [],
            'cmblens':      [],
        }
        
    
    def __getitem__(self, item):
        """
        Retrieve parameter value.
        """
        return self.base[item]
    
    def __setitem__(self, item, val):
        """
        Perform basic checks and then set parameter value.
        """
        assert(item in self.base.keys(), "Parameter '%s' not recognized." % item)
        self.base[item] = val    
    
    def set_defaults(self):
        """
        Set all base cosmological parameters to their defaults, as stated in 
        Section 2 of the CosmoParams spec.
        """
        # Set all parameters to None initially
        for p in self.parameters.keys(): self.base[p] = None
        
        # Set parameters with explicit defaults
        self.base['w0'] = -1.
        self.base['N_nu_rel'] = 3.046
        zero_vals = ['Omega_n_mass', 'Omega_n_rel', 'Omega_g', 'wa', 
                     'N_nu_mass', 'mnu']
        for p in zero_vals: self.base[p] = 0.
    
    
    def apply_consistency(self):
        """
        Apply consistency conditions to subsets of the parameters, as defined 
        in Section 2 of the CosmoParams spec.
        
        Any parameters that are undefined, but can be inferred from a 
        consistency relation, will be filled in.
        """
        raise NotImplementedError()
        
        # MUST satisfy: Omega_m = Omega_c + Omega_b + Omega_n_mass
        pnames = ['Omega_m', 'Omega_c', 'Omega_b', 'Omega_n_mass']
        vals = [self.base[p] for p in pnames]
        if vals.count(None) > 1: raise ValueError("")
        
        # MUST satisfy: 1 - Omega_k = Omega_m + Omega_l + Omega_g + Omega_n_rel
        
    
    def update_consistent(self, param, val, modify=None):
        """
        Update the value of a parameter, but change other parameter values if 
        needed to maintain self-consistency of the parameter set.
        
        This is useful for changing parameters that are related to other 
        parameters by a constraint.
        
        Parameters
        ----------
        param : str
            Name of the parameter to update.
        val : float
            New value for the specified parameter.
        modify : str, optional
            Specify which other parameter should be updated to maintain 
            consistency of the parameter set. If not specified, a default will 
            be used. Default: None.
        """
        # Check that parameters exist
        assert(param in self.base.keys(), "Parameter '%s' not recognized." % param)
        if modify is not None:
            assert(modify in self.base.keys(), 
                   "Parameter '%s' not recognized." % modify)
        
        # List of parameters related by constraints:
        # (1) MUST satisfy: Omega_m = Omega_c + Omega_b + Omega_n_mass
        constr_1 = ['Omega_m', 'Omega_c', 'Omega_b', 'Omega_n_mass']
        
        # (2) MUST satisfy: 1 - Omega_k = Omega_m + Omega_l + Omega_g + Omega_n_rel
        constr_2 = ['Omega_k', 'Omega_m', 'Omega_l', 'Omega_g', 'Omega_n_rel']
        
        # Make sure that param and modify don't have the same value
        if modify == param: 
            raise ValueError("'param' and 'modify' can't be the same.")
        
        # Perform a specific action for each consistent set of params
        if param in constr_1:
            # Check for invalid modify parameter
            if modify is not None:
                assert(modify in constr_1, 
                       "Parameter '%s' can't be used to satisfy a consistency "
                       "relation that affects parameter '%s'." % (modify, param))
            
            # Update parameter
            self.base[param] = value
            
            # Apply consistency relation
            if param == 'Omega_m':
                if modify is None: modify = 'Omega_c'
                constr_1.remove('Omega_m') # Remove Omega_m from list
                self.base[modify] += sum([self[p] for p in constr_1]) - value
            else:
                if modify is None: modify = 'Omega_m'
                constr_1.remove('Omega_m') # Remove Omega_m from list
                if modify == 'Omega_m':
                    self.base['Omega_m'] = sum([self[p] for p in constr_1])
                else:
                    constr_1.remove(modify)
                    self.base[modify] += sum([self[p] for p in constr_1]) \
                                       - self.base['Omega_m']
        
        elif param in constr_2:
            # Check for invalid modify parameter
            if modify is not None:
                assert(modify in constr_2, 
                       "Parameter '%s' can't be used to satisfy a consistency "
                       "relation that affects parameter '%s'." % (modify, param))
            else:
                if param == 'Omega_m':
                    modify = 'Omega_k'
                else:
                    modify = 'Omega_m'
            
            # Update parameter
            self.base[param] = value
            
            # Apply consistency relation
            constr_2.remove(modify)
            self.base[modify] = 1. - sum([self[p] for p in constr_2])
            
        else:
            # Parameter not affected by consistency conditions; just change it
            self.base[param] = val
    
    
    def descriptions(self):
        """
        Return dictionary containing descriptions of the base cosmological 
        parameters.
        
        Returns
        -------
        descs : dict
            Dictionary of parameter descriptions. Keys are parameter names, and 
            values are description strings.
        """
        return self.parameters
    
    def _hash_params(self, params):
        """
        Generate a hash for a specified set of parameters, following the method 
        defined in Section 4 of the CosmoParams spec.
        
        Parameters
        ----------
        params : list of str
            List of names of the parameters to include in the hash. The 
            parameter names will be sorted alphabetically before 
        
        Returns
        -------
        hash : str
            String containing MD5 hex digest for the set of parameters.
        """
        # Sort names alphabetically
        params.sort()
        
        # Loop over parameters, adding to a string with param names and values
        s = ""
        for p in params:
            try:
                val = self[p]
            except:
                raise KeyError("Parameter '%s' not recognized." % p)
            
            # Replace NaN/None with zero
            if val is None or np.isnan(val): val = 0.
            
            # Form string as per Section 4 of the CosmoParams spec
            string = "%s:%10.10e" % (p, val)
            
            # Trim whitespace
            string = string.strip()
            string = string.replace(" ", "")
            string = string.replace("\n", "")
            string = string.replace("\t", "")
            
            # Append
            s += string
        
        # Calculate md5sum on string
        hsh = md5(string).hexdigest()
        return hsh
            
    
    def hash(self, types=None):
        """
        Calculate hashes for 
        
        Parameters
        ----------
        types : list of str, optional
            List of hash types to calculate. These types must be present in the 
            self.hash_types dict. If None, all available hash types will be 
            calculated. Default: None.
        
        Returns
        -------
        hashes : dict
            Dictionary of hashes that were calculated. The key is the name of 
            the hash type and the value is the hash (a string containing an MD5 
            hex digest).
        """
        # For convenience, convert str to list
        if isinstance(types, str): types = [str,]
        
        # Check that requested kind of hash is in list
        if types is not None:
            for t in types:
                assert(t in self.hash_types.keys(), 
                       "Unknown kind of hash '%s' requested." % kind)
        
        # If types=None, get all available hashes
        if types is None: types = self.hash_types.keys()
        
        # Populate dict with the requested hashes
        hashes = {}
        for t in types:
            hashes[t] = self._hash_params(self.hash_types[t])
        return hashes
        

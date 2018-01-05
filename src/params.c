#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "params.h"

/* ------- ROUTINE: cosmo_params_validate ------
   INPUTS: a cosmo_params struct
   TASK: ensure that parameter values are valid and consistent. This is a 
         basic physical validity check; this will not check if the parameters 
         are within valid ranges for the programs that will use them.
*/
void cosmo_params_validate(cosmo_base params){
    
    // Positivity checks
    assert(params.Omega_c >= 0.);
    assert(params.Omega_b >= 0.);
    assert(params.Omega_l >= 0.);
    assert(params.Omega_m >= 0.);
    assert(params.Omega_n_mass >= 0.);
    assert(params.Omega_n_rel >= 0.);
    assert(params.Omega_g >= 0.);
    assert(params.h >= 0.);
    assert(params.A_s >= 0.);
    assert(params.n_s >= 0.);
    assert(params.N_nu_mass >= 0.);
    assert(params.N_nu_rel >= 0.);
    if(!isnan(params.sigma_8)) assert(params.sigma_8 >= 0.);
    
    // Density parameters: Consistency relations
    assert(params.Omega_m == (params.Omega_b + params.Omega_c 
                            + params.Omega_n_mass));
    assert(params.Omega_k == 1. - (params.Omega_m + params.Omega_l 
                                 + params.Omega_g + params.Omega_n_rel));

}

/* ------- ROUTINE: cosmo_params_create_default ------
   INPUTS: non-optional cosmological parameters
   RETURNS: a cosmo_base struct
   TASK: create a new cosmo_params struct for the minimal set of input 
         parameters, with all other parameters set to their defaults.
*/
cosmo_base cosmo_params_create_default(double Omega_c, double Omega_b, 
                                       double Omega_l, double h, double A_s, 
                                       double n_s)
{
    // Instantiate new cosmo_base parameter struct
    cosmo_base params;
    
    // Density parameters
    params.Omega_c = Omega_c;
    params.Omega_b = Omega_b;
    params.Omega_l = Omega_l;
    params.Omega_g = 0.;
    params.Omega_n_mass = 0.;
    params.Omega_n_rel = 0.; // FIXME
    params.Omega_m = params.Omega_c + params.Omega_b + params.Omega_n_mass;
    params.Omega_k = 1. - (  params.Omega_m + params.Omega_l 
                           + params.Omega_g + params.Omega_n_rel );
    
    // Dark energy equation of state parameters
    params.w0 = -1.;
    params.wa = 0.;

    // Neutrino properties
    params.N_nu_mass = 0;
    params.N_nu_rel = 3.046;
    params.mnu = 0.0;
  
    // Primordial power spectrum parameters
    params.A_s = A_s;
    params.n_s = n_s;
    params.sigma_8 = nan("");
    
    // Return
    return params;
}

/* ------- ROUTINE: cosmo_params_print ------
   INPUTS: a cosmo_base struct
   TASK: print values of cosmological parameters
*/
void cosmo_params_print(cosmo_base params)
{
    // Density parameters
    printf("Omega_m:      %4.4f\n", params.Omega_m);
    printf("Omega_c:      %4.4f\n", params.Omega_c);
    printf("Omega_b:      %4.4f\n", params.Omega_b);
    printf("Omega_l:      %4.4f\n", params.Omega_l);
    printf("Omega_k:      %4.4e\n", params.Omega_k);
    printf("Omega_n_mass: %3.3e\n", params.Omega_n_mass);
    printf("Omega_n_rel:  %3.3e\n", params.Omega_n_rel);
    printf("Omega_g:      %3.3e\n", params.Omega_g);
    
    // Dark energy equation of state parameters
    printf("w0:           %3.3f\n", params.w0);
    printf("wa:           %3.3f\n", params.wa);

    // Neutrino properties
    printf("N_nu_mass:    %3.3f\n", params.N_nu_mass);
    printf("N_nu_rel:     %3.3f\n", params.N_nu_rel);
    printf("mnu:          %4.4f\n", params.mnu);
  
    // Primordial power spectrum parameters
    printf("A_s:          %4.4e\n", params.A_s);
    printf("n_s:          %4.4f\n", params.n_s);
    printf("sigma_8:      %4.4f\n", params.sigma_8);
}


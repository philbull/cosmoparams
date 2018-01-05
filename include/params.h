/** @file */
#ifdef __cplusplus
extern "C" {
#endif
#pragma once

// Struct containing the parameters defining a basic cosmology
typedef struct cosmo_base {

  // Densities: CDM, baryons, total matter, neutrinos, curvature
  double Omega_c; // Fractional density of CDM relative to critical density
  double Omega_b; // Fractional density of baryons relative to critical density
  double Omega_m; // Fractional density of all matter relative to critical density
  double Omega_k; // Fractional density of curvature relative to critical density
  double Omega_l; // Fractional density of Lambda/DE relative to critical density
  double Omega_n_mass;  // Omega_nu for MASSIVE neutrinos
  double Omega_n_rel;   // Omega_nu for MASSLESS neutrinos
  double Omega_g; // Fractional density of photons relative to critical density

  // Dark Energy equation of state parameters, for w(z) ~ w_0 + (a - 1) w_a
  double w0;
  double wa;
  
  // Hubble parameter
  double h;

  // Neutrino properties
  // (Base cosmology only supports equal-mass massive neutrinos)
  double N_nu_mass;     // Number of different species of massive neutrinos
  double N_nu_rel;      // Effective number of massless neutrinos
  double mnu;           // Sum of neutrino masses, in eV
  
  // Primordial power spectrum parameters
  double A_s;       // Primordial scalar power spectrum amplitude (dimensionless)
  double n_s;       // Spectral index of primordial scalar power spectrum
  double sigma_8;   // Power spectrum amplitude, sigma_8
  
} cosmo_base;


void cosmo_params_print(cosmo_base params);
void cosmo_params_validate(cosmo_base params);
cosmo_base cosmo_params_create_default(double Omega_c, double Omega_b, 
                                       double Omega_l, double h, double A_s, 
                                       double n_s);

#ifdef __cplusplus
}
#endif

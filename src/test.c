#include "params.h"
#include <stdio.h>

int main(int argc, const char* argv[])
{
	
	double Omega_c = 0.27;
	double Omega_b = 0.04;
    double Omega_l = 0.69;
    double h = 0.67;
    double A_s = 1e-9;
    double n_s = 0.96;
	
	// Create new cosmo_base params struct
	cosmo_base params = cosmo_params_create_default(Omega_c, Omega_b, Omega_l, 
	                                                h, A_s, n_s);
	
	// Validate parameters
	cosmo_params_validate(params);
	
	// Output parameter values
	cosmo_params_print(params);
	
	return 0;
}


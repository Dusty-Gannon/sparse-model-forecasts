////////////////////////////////////////////////
// This Stan program defines a function
// that can be used to define an approximate
// marginal variance of an AR(p) process based
// on the vector of AR coefficients, phi, and
// the variance of the perturbations, sigma2.
//
// @param max_lag The order of the MA process
// used to represent the AR process, from which
// variance calculations are easier.
///////////////////////////////////////////////


functions {

  real ar_var(vector phi, real sigma2, int max_lag){

    // order of AR process
    int p = size(phi);

    // initiate vector of impulse response coefficients
    vector[max_lag] psi = rep_vector(0.0, max_lag);

    // vector is padded with zeros
    psi[p] = 1;

    // conversion to MA-infinity representation
    for(i in (p + 1):max_lag){
      vector[p] rev_psi_sub;
      for(j in 1:p){
        rev_psi_sub[j] = psi[i - j];
      }
      psi[i] = sum(phi .* rev_psi_sub);
    }

    // compute variance based on expression
    // for variane of MA process.
    return(sum(pow(psi, 2)) * sigma2);

  }

}

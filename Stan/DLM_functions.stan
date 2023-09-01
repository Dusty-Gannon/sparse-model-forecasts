///////////////////////////////////////////////
// This stan program fits an ARMA model using a
// Gaussian DLM representation
///////////////////////////////////////////////

functions {

  // function to create observation matrix
  row_vector make_F(int r) {
    row_vector[r] F = rep_row_vector(0, r);
    F[1] = 1;
    return F;
  }

  // function to construct transition matrix
  matrix make_G(vector phi_aug) {
    int r = rows(phi_aug);
    matrix[r, r] G;
    matrix[r - 1, r - 1] D = diag_matrix(rep_vector(1, r - 1));
    row_vector[r - 1] b = rep_row_vector(0, r - 1);

    G[, 1] = phi_aug;
    G[1:(r-1), 2:r] = D;
    G[r, 2:r] = b;

    return G;
  }

  // function to construct R
  vector make_R(vector theta_aug) {
    int r = rows(theta_aug) + 1;
    vector[r] R;

    R[1] = 1;
    R[2:r] = theta_aug;

    return R;
  }

  // function to augment phi vector for DLM
  vector augment_phi(int r, vector phi) {
    int p = rows(phi);
    vector[r] phi_prime = rep_vector(0, r);

    phi_prime[1:p] = phi;

    return phi_prime;
  }

  // function to augment theta vector for DLM
  vector augment_theta(int r, vector theta) {
    int q = rows(theta);
    vector[r - 1] theta_prime = rep_vector(0, r - 1);

    theta_prime[1:q] = theta;

    return theta_prime;
  }

}


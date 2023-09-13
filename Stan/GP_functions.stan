

functions {

  // fill the off-diagonal elements of a matrix with values from a vector //
  matrix fill_off_diag(vector d, vector v){

    int n = rows(d);

    matrix[n, n] M = diag_matrix(d);

    for(j in 1:n){
      for(i in 1:n){
        if(i != j){
          M[i, j] = v[n * (j - 1) + i - (j - 1) - ((i > j) ? 1 : 0)];
        }
      }
    }

    return(M);

  }

  // kronecker product of two matrices //
  matrix kronecker_prod(matrix M, matrix W){

    int n_m = rows(M);
    int p_m = cols(M);
    int n_w = rows(W);
    int p_w = cols(W);

    int N = n_m * n_w;
    int P = p_m * p_w;
    matrix[N, P] K;

    for(i in 1:n_m){
      for(j in 1:p_m){
        for(i2 in 1:n_w){
          for(j2 in 1:p_w){
            K[(i - 1) * n_w + i2, (j - 1) * p_w + j2] = M[i, j] * W[i2, j2];
          }
        }
      }
    }

    return(K);

  }


}


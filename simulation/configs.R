if (scenario == 1) {
  n = 200
  p = 500
  Sigma_X = diag(1, 500)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 0.5
  n_test = 5000
}

if (scenario == 2) {
  n = 200
  p = 500
  Sigma_X = diag(1, 500)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 1
  n_test = 5000
}

if (scenario == 3) {
  n = 200
  p = 500
  Sigma_X = diag(1, 500)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 3
  n_test = 5000
}

if (scenario == 4) {
  n = 200
  p = 100
  Sigma_X = diag(1, 100)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 0.5
  n_test = 5000
}

if (scenario == 5) {
  n = 200
  p = 100
  Sigma_X = diag(1, 100)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 1
  n_test = 5000
}

if (scenario == 6) {
  n = 200
  p = 100
  Sigma_X = diag(1, 100)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 3
  n_test = 5000
}

if (scenario == 7) {
  n = 200
  p = 500
  Sigma_X = diag(1, 500)
  q = 20
  sparsity = 50
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 0.5
  n_test = 5000
}

if (scenario == 8) {
  n = 200
  p = 500
  Sigma_X = diag(1, 500)
  q = 20
  sparsity = 50
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 1
  n_test = 5000
}

if (scenario == 9) {
  n = 200
  p = 500
  Sigma_X = diag(1, 500)
  q = 20
  sparsity = 50
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 3
  n_test = 5000
}

if (scenario == 10) {
  n = 200
  p = 500
  Sigma_X = diag(1, 500)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0.5
  Sigma_E = diag(1, 20)
  SNR = 0.5
  n_test = 5000
}

if (scenario == 11) {
  n = 200
  p = 500
  Sigma_X = diag(1, 500)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0.5
  Sigma_E = diag(1, 20)
  SNR = 1
  n_test = 5000
}

if (scenario == 12) {
  n = 200
  p = 500
  Sigma_X = diag(1, 500)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0.5
  Sigma_E = diag(1, 20)
  SNR = 3
  n_test = 5000
}

if (scenario == 13) {
  n = 200
  p = 100
  Sigma_X = diag(1, 100)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0.5
  Sigma_E = diag(1, 20)
  SNR = 0.5
  n_test = 5000
}

if (scenario == 14) {
  n = 200
  p = 100
  Sigma_X = diag(1, 100)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0.5
  Sigma_E = diag(1, 20)
  SNR = 1
  n_test = 5000
}

if (scenario == 15) {
  n = 200
  p = 100
  Sigma_X = diag(1, 100)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0.5
  Sigma_E = diag(1, 20)
  SNR = 3
  n_test = 5000
}

if (scenario == 16) {
  n = 200
  p = 500
  Sigma_X = diag(1, 500)
  q = 20
  sparsity = 50
  rank = 3
  C_noise = 0.5
  Sigma_E = diag(1, 20)
  SNR = 0.5
  n_test = 5000
}

if (scenario == 17) {
  n = 200
  p = 500
  Sigma_X = diag(1, 500)
  q = 20
  sparsity = 50
  rank = 3
  C_noise = 0.5
  Sigma_E = diag(1, 20)
  SNR = 1
  n_test = 5000
}

if (scenario == 18) {
  n = 200
  p = 500
  Sigma_X = diag(1, 500)
  q = 20
  sparsity = 50
  rank = 3
  C_noise = 0.5
  Sigma_E = diag(1, 20)
  SNR = 3
  n_test = 5000
}

# correlated X

if (scenario == 19) {
  n = 200
  p = 500
  Sigma_X = diag(0.5, 500) + matrix(0.5, 500, 500)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 0.5
  n_test = 5000
}

if (scenario == 20) {
  n = 200
  p = 500
  Sigma_X = diag(0.5, 500) + matrix(0.5, 500, 500)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 1
  n_test = 5000
}

if (scenario == 21) {
  n = 200
  p = 500
  Sigma_X = diag(0.5, 500) + matrix(0.5, 500, 500)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 3
  n_test = 5000
}

if (scenario == 22) {
  n = 200
  p = 100
  Sigma_X = diag(0.5, 100) + matrix(0.5, 100, 100)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 0.5
  n_test = 5000
}

if (scenario == 23) {
  n = 200
  p = 100
  Sigma_X = diag(0.5, 100) + matrix(0.5, 100, 100)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 1
  n_test = 5000
}

if (scenario == 24) {
  n = 200
  p = 100
  Sigma_X = diag(0.5, 100) + matrix(0.5, 100, 100)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 3
  n_test = 5000
}

if (scenario == 25) {
  n = 200
  p = 500
  Sigma_X = diag(0.5, 500) + matrix(0.5, 500, 500)
  q = 20
  sparsity = 50
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 0.5
  n_test = 5000
}

if (scenario == 26) {
  n = 200
  p = 500
  Sigma_X = diag(0.5, 500) + matrix(0.5, 500, 500)
  q = 20
  sparsity = 50
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 1
  n_test = 5000
}

if (scenario == 27) {
  n = 200
  p = 500
  Sigma_X = diag(0.5, 500) + matrix(0.5, 500, 500)
  q = 20
  sparsity = 50
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 3
  n_test = 5000
}

if (scenario == 28) {
  n = 200
  p = 500
  Sigma_X = diag(0.5, 500) + matrix(0.5, 500, 500)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0.5
  Sigma_E = diag(1, 20)
  SNR = 0.5
  n_test = 5000
}

if (scenario == 29) {
  n = 200
  p = 500
  Sigma_X = diag(0.5, 500) + matrix(0.5, 500, 500)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0.5
  Sigma_E = diag(1, 20)
  SNR = 1
  n_test = 5000
}

if (scenario == 30) {
  n = 200
  p = 500
  Sigma_X = diag(0.5, 500) + matrix(0.5, 500, 500)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0.5
  Sigma_E = diag(1, 20)
  SNR = 3
  n_test = 5000
}

if (scenario == 31) {
  n = 200
  p = 100
  Sigma_X = diag(0.5, 100) + matrix(0.5, 100, 100)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0.5
  Sigma_E = diag(1, 20)
  SNR = 0.5
  n_test = 5000
}

if (scenario == 32) {
  n = 200
  p = 100
  Sigma_X = diag(0.5, 100) + matrix(0.5, 100, 100)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0.5
  Sigma_E = diag(1, 20)
  SNR = 1
  n_test = 5000
}

if (scenario == 33) {
  n = 200
  p = 100
  Sigma_X = diag(0.5, 100) + matrix(0.5, 100, 100)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0.5
  Sigma_E = diag(1, 20)
  SNR = 3
  n_test = 5000
}

if (scenario == 34) {
  n = 200
  p = 500
  Sigma_X = diag(0.5, 500) + matrix(0.5, 500, 500)
  q = 20
  sparsity = 50
  rank = 3
  C_noise = 0.5
  Sigma_E = diag(1, 20)
  SNR = 0.5
  n_test = 5000
}

if (scenario == 35) {
  n = 200
  p = 500
  Sigma_X = diag(0.5, 500) + matrix(0.5, 500, 500)
  q = 20
  sparsity = 50
  rank = 3
  C_noise = 0.5
  Sigma_E = diag(1, 20)
  SNR = 1
  n_test = 5000
}

if (scenario == 36) {
  n = 200
  p = 500
  Sigma_X = diag(0.5, 500) + matrix(0.5, 500, 500)
  q = 20
  sparsity = 50
  rank = 3
  C_noise = 0.5
  Sigma_E = diag(1, 20)
  SNR = 3
  n_test = 5000
}

## big application ##

if (scenario == 37) {
  n = 10000
  p = 20000
  Sigma_X = 0.5
  q = 20
  sparsity = 500
  rank = 4
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 1
  n_test = 2000
}

# missing value # --- simulation_missingness file
if (scenario == 38) {
  n = 200
  p = 500
  Sigma_X = diag(1, 500)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 0.5
  n_test = 5000
  rank_choice = 2
}

# missing value #
if (scenario == 39) {
  n = 200
  p = 500
  Sigma_X = diag(1, 500)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 0.5
  n_test = 5000
  rank_choice = 3
}

# missing value #
if (scenario == 40) {
  n = 200
  p = 500
  Sigma_X = diag(1, 500)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 0.5
  n_test = 5000
  rank_choice = 4
}

# missing value #
if (scenario == 41) {
  n = 200
  p = 500
  Sigma_X = diag(1, 500)
  q = 20
  sparsity = 20
  rank = 3
  C_noise = 0
  Sigma_E = diag(1, 20)
  SNR = 0.5
  n_test = 5000
  rank_choice = 5
}

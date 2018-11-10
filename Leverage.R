# HAT MATRIX CALCS

X_mat <- cbind(rep(1,30000),seq(1,30000,1))
G_mat <- crossprod(X_mat)
G_inv <- solve(G_mat)
H_mat <- (X_mat%*%G_inv)%*%t(X_mat)
H_values <- diag(H_mat)
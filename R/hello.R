x <- matrix(seq(-pi,pi,length.out = 100))
y <- sin(x)


# Testing the GP-BART
bart_test <- bart(x_train = x,y_train = y,x_test = x,n_tree = 200,n_mcmc = 1000,
                  n_burn = 0,tau = 1,mu = 1,
                  tau_mu = 1,naive_sigma = 1,alpha = 0.95,
                  beta = 0.1,a_tau = 1,d_tau = 1,nsigma = 1)


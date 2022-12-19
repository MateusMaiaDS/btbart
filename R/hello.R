x <- matrix(seq(-pi,pi,length.out = 100))
y <- sin(x)

# Test functions
tree_loglike <- function(x,y, tau = 1, tau_mu=1){
     n_obs <- length(y)
     return(-0.5*tau*sum(y^2)-0.5*log(tau_mu+tau*n_obs)+(0.5*tau^2*(sum(y)^2))/(tau_mu+tau*n_obs))
}

tree_loglike_rule <- function(x,y, tau = 1, tau_mu=1,var_split,var_split_rule){
     n_obs <- length(y)
     node_one <- y[x[,var_split]<=var_split_rule]
     node_two <- y[x[,var_split]>var_split_rule]
     n_obs_one <- length(node_one)
     n_obs_two <- length(node_two)
     log_node_one <- -0.5*tau*sum(node_one^2)-0.5*log(tau_mu+tau*n_obs_one)+(0.5*tau^2*(sum(node_one)^2))/(tau_mu+tau*n_obs_one)
     log_node_two <- -0.5*tau*sum(node_two^2)-0.5*log(tau_mu+tau*n_obs_two)+(0.5*tau^2*(sum(node_two)^2))/(tau_mu+tau*n_obs_two)

     print(log_node_one)
     print(log_node_two)
     return(log_node_one + log_node_two)
}

microbenchmark::microbenchmark(tree_loglike(x = x,y = y),test_logtree(X = x,y = y),times = 10000)


test_logtree(X = x,y = y)

testingDisplay()
var_split <- 1
var_split_rule <- 0
tree_loglike_rule(x = x,y = y,tau = 1.0,tau_mu = 1,var_split = 1,var_split_rule = 0.67)

# Doing more tests
createTree(X = x,y = y,n_tree = 1)

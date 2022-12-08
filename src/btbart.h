#include<RcppArmadillo.h>

// Creating the struct
struct Node;
struct modelParam;

struct modelParam {

        arma::mat X;
        arma::vec y;

        // BART prior param specification
        double alpha;
        double beta;
        double tau_mu;
        double tau;
        double a_tau;
        double d_tau;

        // MCMC spec.
        int n_mcmc;
        int n_burn;



};

// Creating a forest
class Forest {

public:
        std::vector<Node*> trees;
        modelParam modelParam_;

        Forest(const arma::mat &X, int n_tree);
        ~Forest();
};



// Creating the node struct
struct Node {

     bool isRoot;
     bool isLeaf;
     Node* left;
     Node* right;
     Node* parent;

     // Branch parameters
     int var_split;
     double var_split_rule;
     double lower;
     double upper;
     double curr_weight; // indicates if the observation is within terminal node or not


     // Leaf parameters
     double mu;

     // Storing sufficient statistics over the nodes
     double r_sq_sum = 0;
     double r_sum = 0;
     double log_likelihood = 0;
     double n_leaf = 0.0;

     // Displaying and check nodes
     void displayNode();
     void displayCurrNode();

     // Creating the methods
     void addingLeaves();
     void deletingLeaves();
     void Stump();
     void updateWeight(const arma::mat X, int i);
     void getLimits(); // This function will get previous limit for the current var
     void sampleSplitVar(int p);
     bool isLeft();
     bool isRight();
     void grow();

     Node();
     ~Node();
};

// Creating a function to get the leaves
void leaves(Node* x, std::vector<Node*>& leaves); // This function gonna modify by address the vector of leaves
std::vector<Node*> leaves(Node*x);
// [[Rcpp::export]]
double rand_unif(){
        double rand_d = std::rand();
        return rand_d/RAND_MAX;
};

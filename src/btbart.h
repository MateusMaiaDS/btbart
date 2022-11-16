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

        // modelParam(Rcpp::List& modelParam_obj);


};

// Creating a forest
class Forest {

public:
        std::vector<Node*> trees;
        modelParam modelParam_;

        Forest(arma::mat &X, int n_tree, Rcpp::List modelParam_obj);
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


     // Displaying and check nodes
     void displayNode();

     // Creating the methods
     void addingLeaves();
     void deletingLeaves();
     void Stump();
     void updateWeight(const arma::mat X, int i);

     Node();
     ~Node();
};

// Creating a function to get the leaves
void leaves(Node* x, std::vector<Node*>& leaves); // This function gonna modify by address the vector of leaves
std::vector<Node*> leaves(Node*x);

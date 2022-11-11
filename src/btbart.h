#include<RcppArmadillo.h>

// Creating the struct
struct Node;

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


     // Leaf parameters
     double mu;

     // Creating the methods
     void addingLeaves();
     void deletingLeaves();

};

// Creating a function to get the leaves
void leaves(Node* x, std::vector<Node*>& leaves); // This function gonna modify by address the vector of leaves
std::vector<Node*> leaves(Node*x);

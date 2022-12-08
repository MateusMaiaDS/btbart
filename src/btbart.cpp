#include "btbart.h"

using namespace std;
// Initialising a node
Node::Node(){
        isLeaf = true;
        isRoot = true;
        left = NULL;
        right = NULL;
        parent = NULL;

        var_split = 0;
        var_split_rule = 0.0;
        lower = 0.0;
        upper = 1.0;
        curr_weight = 0.0;
        mu = 0.0;
        n_leaf = 0.0;
        r_sq_sum = 0.0;
        r_sum = 0.0;
        log_likelihood = 0.0;

}

Node::~Node() {
        if(!isLeaf) {
                delete left;
                delete right;
        }
}

// Initializing a stump
void Node::Stump(){

        // Changing the left parent and right nodes;
        left = this;
        right = this;
        parent = this;

}

void Node::addingLeaves(){

     // Create the two new nodes
     left = new Node; // Creating a new vector object to the
     right = new Node;
     isLeaf = false;

     // Modifying the left node
     left -> isRoot = false;
     left -> isLeaf = true;
     left -> left = left;
     left -> right = left;
     left -> parent = this;
     left -> var_split = 0;
     left -> var_split_rule = 0.0;
     left -> lower = 0.0;
     left -> upper = 1.0;
     left -> mu = 0.0;
     left -> r_sq_sum = 0.0;
     left -> r_sum = 0.0;
     left -> log_likelihood = 0.0;
     left -> n_leaf = 0.0;

     right -> isRoot = false;
     right -> isLeaf = true;
     right -> left = right; // Recall that you are saving the address of the right node.
     right -> right = right;
     right -> parent = this;
     right -> var_split = 0;
     right -> var_split_rule = 0.0;
     right -> lower = 0.0;
     right -> upper = 1.0;
     right -> mu = 0.0;
     right -> r_sq_sum = 0.0;
     right -> r_sum = 0.0;
     right -> log_likelihood = 0.0;
     right -> n_leaf = 0.0;


     return;

}

// Creating boolean to check if the vector is left or right
bool Node::isLeft(){
        return (this == this->parent->left);
}

bool Node::isRight(){
        return (this == this->parent->right);
}

// Sample var
void Node::sampleSplitVar(int p){

        // Sampling one index from 0:(p-1)
        var_split = std::rand()%p;

}
// This functions will get and update the current limits for this current variable
void Node::getLimits(){

        // Creating  a new pointer for the current node
        Node* x = this;
        // Already defined this -- no?
        lower = 0.0;
        upper = 1.0;
        // First we gonna check if the current node is a root or not
        bool tree_iter = x->isRoot ? false: true;
        while(tree_iter){
                bool is_left = x->isLeft(); // This gonna check if the current node is left or not
                x = x->parent; // Always getting the parent of the parent
                tree_iter = x->isRoot ? false : true; // To stop the while
                if(x->var_split == var_split){
                        tree_iter = false ; // This stop is necessary otherwise we would go up til the root, since we are always update there is no prob.
                        if(is_left){
                                upper = x->var_split_rule;
                                lower = x->lower;
                        } else {
                                upper = x->upper;
                                lower = x->var_split_rule;
                        }
                }
        }
}

void Node::displayCurrNode(){

                std::cout << "Node address: " << this << std::endl;
                std::cout << "Node parent: " << parent << std::endl;

                std::cout << "Cur Node is leaf: " << isLeaf << std::endl;
                std::cout << "Cur Node is root: " << isRoot << std::endl;
                std::cout << "Cur The split_var is: " << var_split << std::endl;
                std::cout << "Cur The split_var_rule is: " << var_split_rule << std::endl;

                return;
}

void Node::displayNode(){

        if(isLeaf){
                std::cout << "Node is leaf: " << isLeaf << std::endl;
                std::cout << "Node is root: " << isRoot << std::endl;
                std::cout << "The split_var is: " << var_split << std::endl;
                std::cout << "The split_var_rule is: " << var_split << std::endl;

        } else {
                left->displayNode();
                right->displayNode();
        }
        return;
}

void Node::deletingLeaves(){

     // Should I create some warn to avoid memoery leak
     //something like it will only delete from a nog?
     // Deleting
     delete left; // This release the memory from the left point
     delete right; // This release the memory from the right point
     left = this;  // The new pointer for the left become the node itself
     right = this; // The new pointer for the right become the node itself
     isLeaf = true;

     return;

}
// Getting the leaves (this is the function that gonna do the recursion the
//                      function below is the one that gonna initialise it)
void get_leaves(Node* x,  std::vector<Node*> &leaves_vec) {

        if(x->isLeaf){
                leaves_vec.push_back(x);
        } else {
                get_leaves(x->left, leaves_vec);
                get_leaves(x->right,leaves_vec);
        }

        return;

}



// Initialising a vector of nodes in a standard way
std::vector<Node*> leaves(Node* x) {
        std::vector<Node*> leaves_init(0); // Initialising a vector of a vector of pointers of nodes of size zero
        get_leaves(x,leaves_init);
        return(leaves_init);
}

// Initializing the forest
Forest::Forest(const arma::mat& X, int n_tree){

        // Creatina vector of size of number of trees
        trees.resize(n_tree);
        for(int  i=0;i<n_tree;i++){
                // Creating the stump for each tree
                trees[i] = new Node();
                // Filling up each stump for each tree
                trees[i]->Stump();
        }
}

// Function to delete one tree
Forest::~Forest(){
        for(int  i=0;i<trees.size();i++){
                delete trees[i];
        }
}

// Selecting a random node
Node* sample_node(std::vector<Node*> leaves_){

        // Getting the number of leaves
        int n_leaves = leaves_.size();
        return(leaves_[std::rand()%n_leaves]);
}

// Grow a tree for a given rule
void Node::grow(){

        if(isLeaf){
                addingLeaves();
                // Selecting the var
                sampleSplitVar(1);
                // Updating the limits
                getLimits();
                // Selecting a rule
                var_split_rule = (upper-lower)*rand_unif()+lower;

                if(isRoot){
                        isRoot = false;
                }
        }
}

// Calculating the LLT for a tree
double treeLogLike(Node* tree, const arma::vec& y,
                   const arma::mat& X, double tau, double tau_mu){

        // Getting the number of terminal nodes
        std::vector<Node*> t_nodes = leaves(tree) ;


        // Summing-up
        double tree_log_likelihood = 0;
        // Getting number of levaes
        int n_leaves = t_nodes.size();
        arma::vec w_t_nodes = arma::zeros<arma::vec>(n_leaves);
        arma::vec n_t_nodes = arma::zeros<arma::vec>(n_leaves);
        arma::vec t_node_r = arma::zeros<arma::vec>(n_leaves);
        arma::vec t_node_r_sq = arma::zeros<arma::vec>(n_leaves);

        // If it's a root case
        if(tree->isRoot){
                return -0.5*tau*dot(y,y) - 0.5*log(tau_mu + (y.size()*tau)) + (0.5*(tau*tau)*(sum(y)*sum(y)))/(tau_mu + (y.size()*tau));
        }

        cout << "Number of leaves " <<  n_leaves << endl;

        // Selecting it observation
        for(int i=0;i<X.n_rows;i++){
                tree->updateWeight(X,i); // Getting the weights for each terminal node
                w_t_nodes = arma::zeros<arma::vec>(n_leaves); // Putting zero on the weights terminal nodes
                        for(int l = 0; l<n_leaves;l++){
                                if(t_nodes[l]->curr_weight==1){ // This if is just because you only have 1 value per t_node
                                        w_t_nodes(l) = 1; // I might do not need go through all terminal nodes
                                        break;
                                }
                        }

                // Calculating the sufficient statistics
                n_t_nodes = n_t_nodes + w_t_nodes;
                t_node_r = t_node_r + y(i)*w_t_nodes;
                t_node_r_sq = t_node_r_sq + y(i)*y(i)*w_t_nodes;
        }

        cout << "ALL NODES " <<sum(n_t_nodes) << endl;
        // Getting the sum and the sum_sq
        for(int l = 0; l<n_leaves;l++){
                t_nodes[l]->r_sum = t_node_r(l);
                t_nodes[l]->r_sq_sum = t_node_r_sq(l);
                t_nodes[l]->n_leaf = n_t_nodes(l);
                tree_log_likelihood = tree_log_likelihood - 0.5*tau*t_node_r_sq(l) - 0.5*log(tau_mu + (n_t_nodes(l)*tau)) + (0.5*(tau*tau)*(t_node_r(l)*t_node_r(l)))/( (tau*n_t_nodes(l))+tau_mu);
        }


        return tree_log_likelihood;
}

// Update Weight
// (this function gonna run all terminal nodes and verify if a observation is
// inside this terminal node or not)
void Node::updateWeight(const arma::mat X, int i){

        if(!isLeaf){
                if(X(i,var_split)<=var_split_rule){
                        left->curr_weight = 1;
                        right->curr_weight = 0;
                } else {
                        left->curr_weight = 0;
                        right->curr_weight = 1;
                }
        }
        else {
                left->updateWeight(X,i);
                right->updateWeight(X,i);
        }

        return;
}

// TESTING FUNCTIONS

// Creating a function to create a tree and select terminal nodes to be grown
// [[Rcpp::export]]
void createTree(const arma::mat X,
                const arma::vec y,
                int n_tree){

        Forest all_trees = Forest(X,n_tree);

        // Getting the leaves of my first tree
        for(int k=0;k<10;k++){
                std::vector<Node*> tree_zero_leafs = leaves(all_trees.trees[0]);
                Node* g_leaf = sample_node(tree_zero_leafs);

                cout << "Verifying the first tree: " << all_trees.trees[0] <<endl;
                cout << "Verifying the first tree: " << g_leaf <<endl;

                // Displaying this node before grow
                g_leaf->displayCurrNode();
                // Now I gonna grow the chosen node
                g_leaf->grow();
                g_leaf->displayCurrNode();

                cout << " ================= " << endl;
                cout << " ======= " << k <<" ======" << endl;
                cout << " ================= " << endl;

        }

}

// Testing likelihood calculation

// [[Rcpp::export]]
double test_logtree(const arma::mat X,
                    const arma::vec y,
                    double tau,
                    double tau_mu,
                    int new_split_var,
                    double new_split_var_rule){

        // Creating a tree
        Node node_init ;
        std::cout << "Root loglike" << treeLogLike(&node_init,y,X,tau,tau_mu) <<std::endl;

        // Adding one node
        // node_init.addingLeaves();
        // node_init.isRoot = false;
        // node_init.var_split = new_split_var;
        // node_init.var_split_rule = new_split_var_rule;

        node_init.grow();
        cout<< " ===== " << endl;

        node_init.displayCurrNode();
        Node* point_test = &node_init;
        std::cout << "First split: " << treeLogLike(point_test,y,X,tau,tau_mu) <<std::endl;


        return 0.0;
}

// [[Rcpp::export]]
void adding_two_vec(arma::vec x, arma::vec y) {

        arma::vec z = x+y;
        // Adding those two values
        for(int i = 0; i<x.size(); i++) {
                std::cout << "Z position" << i << "is given by " << z[i] << std::endl;
        }

        return;
}


// [[Rcpp::export]]
void testingDisplay(){
        Node node_init;
        node_init.displayNode();
        node_init.addingLeaves();
        std::vector<Node*> return_leaves_init  = leaves(&node_init);
        std::cout << "Final try first " << return_leaves_init.size()<< std::endl;

        node_init.left->addingLeaves();
        return_leaves_init  = leaves(&node_init);
        std::cout << "Final try second " << return_leaves_init.size()<< std::endl;

        node_init.displayNode();

        Node* point_tree = &node_init;
        // Getting all leaves
        std::vector<Node*> return_leaves = leaves(point_tree);

        for(int i=0;i<return_leaves.size();i++){
                std::cout << "Terminal leaf left: " << return_leaves[i]->left<<std::endl;
                std::cout << "Terminal leaf right: " << return_leaves[i]->left<<std::endl;

        }
        return;
}


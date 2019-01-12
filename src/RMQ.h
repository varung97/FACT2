#ifndef RMQ_H_
#define RMQ_H_

#include <sdsl/bit_vectors.hpp>

class RMQ {
public:
    // Tree must be ordered such that for any node, the first child is the one with largest leafset
    RMQ(Tree* tree);
    
    // Max weight along path from v to w
    // w must be an ancestor of v
    int rmq(Tree::Node* v, Tree::Node* w);

private:
    // Pointer from each node to the root of its centroid path
    std::vector<Tree::Node*> cp_roots;
    
    // Pointers to rmq for each of the centroid paths
    std::vector<gen_rmq_t*> cp_rmqs;

    // Index of the centroid path each node belongs to
    std::vector<int> cp_indices;
    
    // Depth of each centroid path (defined as number of centroid path roots that are proper ancestors of the root of this path)
    std::vector<int> cp_depths;
    
    // Depth of each node in the tree
    std::vector<int> depths;
    
    // Pointers from each node to the rmq for some leaf in its leafset
    std::vector<gen_rmq_t*> leaf_rmq_for_nodes;
    
    // Numbering leaves left to right,
    // For each centroid path, for each node, store the minimum leaf index in all side trees as integers
    // Numbering leaves left to right, minimum leaf index in all side trees of a node
    // For leaves, this is just their own index
    std::vector<std::vector<int>> min_leaf_idx_in_side_trees_cps_int;
    
    // For each centroid path, store a rank structure containing min leaf indices
    // Note that these are offset to the index of the smallest leaf in the centroid path to save space
    std::vector<sdsl::rank_support_v<1>> min_leaf_idx_in_side_trees_cp_ranks;
    
    void preprocess_depths(Tree* tree);
    void preprocess_centroid_paths(Tree* tree);
    void preprocess_leaf_rmqs(Tree* tree);
    int preprocess_min_leaf_idx_in_side_trees_helper(Tree::Node* node, int first_available_index);
    void preprocess_min_leaf_idx_in_side_trees(Tree* tree);
    void populate_min_leaf_idx_in_side_trees_cps(Tree* tree);
    void preprocess_leaf_rank_structures(Tree* tree);
    
    int rmq_q1(Tree::Node* v, Tree::Node* w);
    int rmq_q2_qg1(Tree::Node* v, Tree::Node*w);
    int rmq_qg(Tree::Node* v, Tree::Node* w);
};

#endif /* RMQ_hpp */

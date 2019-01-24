#include "Tree.h"
#include "lca_preprocessing.h"
#include <sdsl/bit_vectors.hpp>

#include "RMQ.h"

// Build depths array
void RMQ::preprocess_depths(Tree* tree) {
    depths = std::vector<int>(tree->get_nodes_num());
    depths[0] = 0;
    for (int i = 1; i < tree->get_nodes_num(); i++) {
        depths[i] = depths[tree->get_node(i)->parent->id] + 1;
    }
}

// Build and preprocess centroid path related information
void RMQ::preprocess_centroid_paths(Tree* tree) {
    cp_roots = std::vector<Tree::Node*>(tree->get_nodes_num());
    cp_rmqs = std::vector<gen_rmq_t*>();
    cp_indices = std::vector<int>(tree->get_nodes_num());
    cp_depths = std::vector<int>(tree->get_nodes_num());
    
    // Build and preprocess centroid path related information
    cp_roots[0] = tree->get_root();
    cp_rmqs.push_back(new gen_rmq_t);
    cp_rmqs[0]->v.push_back(-tree->get_root()->weight);
    cp_indices[0] = 0;
    cp_depths[0] = 0;
    
    // Assumes that nodes are iterated top-down
    for (int i = 1; i < tree->get_nodes_num(); i++) {
        Tree::Node* node = tree->get_node(i);
        if (node->pos_in_parent == 0) {
            cp_roots[i] = cp_roots[node->parent->id];
            cp_indices[i] = cp_indices[node->parent->id];
            cp_rmqs[cp_indices[i]]->v.push_back(-node->weight);
            cp_depths[i] = cp_depths[node->parent->id];
        } else {
            // Create new centroid path
            cp_roots[i] = node;
            
            gen_rmq_t* cp_rmq = new gen_rmq_t;
            cp_rmq->v.push_back(-node->weight);
            cp_rmqs.push_back(cp_rmq);
            
            cp_indices[i] = cp_rmqs.size() - 1;
            cp_depths[i] = cp_depths[node->parent->id] + 1;
        }
    }
    
    for (gen_rmq_t* cp_rmq : cp_rmqs) {
        if (cp_rmq->v.size() > 1) {
            general_rmq_preprocess(cp_rmq);
        }
    }
}

// Build centroid subpath rmq for each leaf
void RMQ::preprocess_leaf_rmqs(Tree* tree) {
    preprocess_centroid_paths(tree);
    preprocess_depths(tree);
    
    // Assumes all taxa are represented in tree
    std::vector<gen_rmq_t*> leaf_rmqs = std::vector<gen_rmq_t*>(Tree::get_taxas_num());
    
    for (int taxon = 0; taxon < Tree::get_taxas_num(); taxon++) {
        gen_rmq_t* leaf_rmq = new gen_rmq_t;

        Tree::Node* curr_node = tree->get_leaf(taxon);
        
        // Get the max weight of each centroid subpath on path from leaf to root
        while (true) {
            Tree::Node* cp_root = cp_roots[curr_node->id];

            int depth_in_cp = depths[curr_node->id] - depths[cp_root->id];
            gen_rmq_t* curr_rmq = cp_rmqs[cp_indices[curr_node->id]];
            
            // If node is already the root of centroid path, then no need to query the rmq
            int max_weight = depth_in_cp == 0 ? curr_rmq->v[0] : curr_rmq->v[general_rmq(curr_rmq, 0, depth_in_cp)];
            
            leaf_rmq->v.push_back(max_weight);
            
            if (cp_root == tree->get_root()) {
                break;
            } else {
                curr_node = cp_root->parent;
            }
        }
        
        if (leaf_rmq->v.size() > 1) {
            general_rmq_preprocess(leaf_rmq);
        }
        
        leaf_rmqs[taxon] = leaf_rmq;
    }
    
    // Point each node to a leaf rmq (for any leaf in its leafset)
    leaf_rmq_for_nodes = std::vector<gen_rmq_t*>(tree->get_nodes_num());

    // Assumes nodes are iterated bottom-up
    for (int i = tree->get_nodes_num() - 1; i >= 0; i--) {
        Tree::Node* curr_node = tree->get_node(i);
        if (curr_node->is_leaf()) {
            leaf_rmq_for_nodes[i] = leaf_rmqs[curr_node->taxa];
        } else {
            leaf_rmq_for_nodes[i] = leaf_rmq_for_nodes[curr_node->children[0]->id];
        }
    }
}

// Numbers leaves in subtree of node from left to right, starting at first_available_index
// Populates taxa_l2r_indices, left_most_leaf, child_for_leaf_vectors and heaviest_child
int RMQ::preprocess_child_for_leaf_helper(Tree::Node* node, int first_available_index) {
    if (node->is_leaf()) {
        taxa_l2r_indices[node->taxa] = first_available_index;
        left_most_leaf[node->id] = node;
        return first_available_index + 1;
    }
    
    // Preprocess heaviest child first and populate related arrays
    first_available_index = preprocess_child_for_leaf_helper(node->children[0], first_available_index);
    left_most_leaf[node->id] = left_most_leaf[node->children[0]->id];
    child_for_leaf_offset[node->id] = first_available_index;
    heaviest_child[node->id] = node->children[0];
    
    std::vector<Tree::Node*> child_for_leaf = std::vector<Tree::Node*>();
    int current_leaf_index = 0;
    for (int i = 1; i < node->get_children_num(); i++) {
        int new_first_available_index = preprocess_child_for_leaf_helper(node->children[i], first_available_index);
        
        // Point each leaf to the child it is a descendant of
        // Offset the leaf index by the smallest index in the side trees
        for (; current_leaf_index < new_first_available_index - first_available_index; current_leaf_index++) {
            child_for_leaf.push_back(node->children[i]);
        }
        
        first_available_index = new_first_available_index;
    }
    
    child_for_leaf_vectors[node->id] = child_for_leaf;
    
    return first_available_index;
}

// Numbers leaves from left to right, starting at first_available index
// Populates taxa_l2r_indices, left_most_leaf, child_for_leaf_vectors and heaviest_child
void RMQ::preprocess_child_for_leaf(Tree* tree) {
    taxa_l2r_indices = std::vector<int>(Tree::get_taxas_num());
    left_most_leaf = std::vector<Tree::Node*>(tree->get_nodes_num());
    child_for_leaf_vectors = std::vector<std::vector<Tree::Node*>>(tree->get_nodes_num());
    child_for_leaf_offset = std::vector<int>(tree->get_nodes_num());
    heaviest_child = std::vector<Tree::Node*>(tree->get_nodes_num());
    preprocess_child_for_leaf_helper(tree->get_root(), 0);
}

RMQ::RMQ(Tree* tree, lca_t* lca_prep_in) : lca_prep(lca_prep_in) {
    preprocess_leaf_rmqs(tree);
    preprocess_child_for_leaf(tree);
}

// Max weight from v to the root of its centroid path
int RMQ::rmq_q1(Tree::Node* v, Tree::Node* w) {
    Tree::Node* cp_root = cp_roots[v->id];
    int depth_in_cp_v = depths[v->id] - depths[cp_root->id];
    gen_rmq_t* curr_rmq = cp_rmqs[cp_indices[v->id]];
    
    if (depth_in_cp_v == 0) {
        // If node is already the root of centroid path, then no need to query the rmq
        // Negate answer since values in rmqs are negative
        return -curr_rmq->v[0];
    } else if (cp_indices[v->id] == cp_indices[w->id]) {
        // If v and w are in the same centroid path, then query for the path between them
        int depth_in_cp_w = depths[w->id] - depths[cp_root->id];
        return -curr_rmq->v[general_rmq(curr_rmq, depth_in_cp_w, depth_in_cp_v)];
    } else {
        // The path from v to root of centroid path must be between v and w
        return -curr_rmq->v[general_rmq(curr_rmq, 0, depth_in_cp_v)];
    }
}

// Max weight from parent of cp root of v to shallowest cp root that is a proper descendant of w
int RMQ::rmq_q2_qg1(Tree::Node* v, Tree::Node* w) {
    if (cp_depths[v->id] - cp_depths[w->id] < 2) {
        // No centroid subpaths between v and w
        return 0;
    }
    
    gen_rmq_t* curr_rmq = leaf_rmq_for_nodes[v->id];
    // The leaf rmq must have at least 3 entries, so no need to check for size > 1
    return -curr_rmq->v[general_rmq(curr_rmq, cp_depths[w->id] + 1, cp_depths[v->id] - 1)];
}

// Max weight from w to deepest node in centroid path of w that is an ancestor of v
int RMQ::rmq_qg(Tree::Node* v, Tree::Node* w) {
    if (cp_depths[v->id] - cp_depths[w->id] == 0) {
        // v and w are on the same centroid path
        return 0;
    }
    
    // Take lca of v and the leftmost leaf of w
    int deepest_node = lca(lca_prep, v->id, left_most_leaf[w->id]->id);
    
    gen_rmq_t* curr_rmq = cp_rmqs[cp_indices[w->id]];
    return -curr_rmq->v[general_rmq(curr_rmq, depths[deepest_node] - depths[cp_roots[deepest_node]->id], depths[w->id] - depths[cp_roots[w->id]->id])];
}

// Max weight along path from v to w
// w must be an ancestor of v
int RMQ::rmq(Tree::Node* v, Tree::Node* w) {
    return std::max({rmq_q1(v, w), rmq_q2_qg1(v, w), rmq_qg(v, w)});
}

// Returns the child of w that is an ancestor of v
// w must be a proper ancestor of v
Tree::Node* RMQ::child_for_descendant(Tree::Node* v, Tree::Node* w) {
    // If heaviest child is ancestor of v, return that
    Tree::Node* heaviest_child_w = heaviest_child[w->id];
    if (lca(lca_prep, v->id, heaviest_child_w->id) == heaviest_child_w->id) {
        return heaviest_child_w;
    }
    
    // Otherwise, query children rank structure
    int left_most_leaf_v = taxa_l2r_indices[left_most_leaf[v->id]->taxa];
    int child_for_leaf_offset_w = child_for_leaf_offset[w->id];
    return child_for_leaf_vectors[w->id][left_most_leaf_v - child_for_leaf_offset_w];
}

int RMQ::max_weight_path(Tree::Node* v, Tree::Node* w) {
    return rmq(v, child_for_descendant(v, w));
}

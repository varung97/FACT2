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

// Numbers leaves is subtree of node from left to right, starting at first_available index
// Populates left_most_leaf and right_most_leaf
// Returns the new first_available_index
int RMQ::preprocess_left_right_most_leaf_helper(Tree::Node* node, int first_available_index) {
    if (node->is_leaf()) {
        left_most_leaf[node->id] = right_most_leaf[node->id] = first_available_index;
        return first_available_index + 1;
    }
    
    left_most_leaf[node->id] = first_available_index;
    
    for (Tree::Node* child : node->children) {
        first_available_index = preprocess_left_right_most_leaf_helper(child, first_available_index);
    }
    
    right_most_leaf[node->id] = first_available_index - 1;
    
    return first_available_index;
}

// Numbers leaves from left to right, starting at first_available index
// Populates left_most_leaf and right_most_leaf
void RMQ::preprocess_left_right_most_leaf(Tree* tree) {
    left_most_leaf = std::vector<int>(tree->get_nodes_num());
    right_most_leaf = std::vector<int>(tree->get_nodes_num());
    preprocess_left_right_most_leaf_helper(tree->get_root(), 0);
}

// Populates the rank structures for each centroid path, containing leftmost leaf in the side trees of each node in that path
void RMQ::populate_min_leaf_idx_in_side_trees_cps(Tree* tree) {
    // Number of centroid paths is the same as number of leaves
    std::vector<sdsl::bit_vector*> min_leaf_idx_in_side_trees_cps_bv = std::vector<sdsl::bit_vector*>(tree->get_leaves_num());
    min_leaf_idx_in_side_trees_cps_ranks = std::vector<sdsl::rank_support_v<1>>(tree->get_leaves_num());
    
    // Assumes all taxa are represented in tree
    for (int i = 0; i < Tree::get_taxas_num(); i++) {
        int cp_root_id = cp_roots[tree->get_leaf(i)->id]->id;
        
        // Each entry is offset by the smallest one to reduce space taken by the bit vector
        int bit_vector_size = right_most_leaf[cp_root_id] - left_most_leaf[cp_root_id] + 1;
        min_leaf_idx_in_side_trees_cps_bv[cp_indices[cp_root_id]] = new sdsl::bit_vector(bit_vector_size, 0);
    }
    
    // Populate the bit vectors
    for (int i = 0; i < tree->get_nodes_num(); i++) {
        Tree::Node* node = tree->get_node(i);
        
        if (node->is_leaf()) {
            continue;
        }
        
        // We also subtract 1 from each entry. The reason is illustrated below:
        // Say the entries are 0, 1, 4
        // Then the leafs in the side trees of the second node in the centroid path are 1, 2, 3
        // All of these should give us the same answer when a rank query is made; however 1 gives a different result
        // To account for this we subtract 1 from the entries
        (*min_leaf_idx_in_side_trees_cps_bv[cp_indices[i]])[left_most_leaf[node->children[1]->id] - left_most_leaf[i] - 1] = 1;
    }
    
    // Create rank support structures for the bit vectors
    for (int i = 0; i < min_leaf_idx_in_side_trees_cps_bv.size(); i++) {
        min_leaf_idx_in_side_trees_cps_ranks[i] = *new sdsl::rank_support_v<1>(min_leaf_idx_in_side_trees_cps_bv[i]);
    }
}

// Populates rank structures for each internal node, containing min leaf indices for each child (excluding heaviest child)
void RMQ::populate_min_leaf_idx_children(Tree* tree) {
    min_leaf_idx_children_ranks = std::vector<sdsl::rank_support_v<1>>(tree->get_nodes_num());
    
    for (int i = 0; i < tree->get_nodes_num(); i++) {
        Tree::Node* node = tree->get_node(i);
        
        if (node->is_leaf()) {
            continue;
        }
        
        // Each entry is offset by leftmost leaf of second child to reduce space taken by the bit vector
        int bit_vector_size = right_most_leaf[node->children.back()->id] - left_most_leaf[node->children[1]->id] + 1;
        sdsl::bit_vector* min_leaf_idx = new sdsl::bit_vector(bit_vector_size, 0);
        
        // Skip second child since its entry will become negative after subtracting 1
        for (int j = 2; j < node->get_children_num(); j++) {
            // Similarly as in populate_min_leaf_idx_in_side_trees_cps, we subtract 1 from each entry
            (*min_leaf_idx)[left_most_leaf[node->children[j]->id] - left_most_leaf[node->children[1]->id] - 1] = 1;
        }
        
        // Create a rank support structure for the bit vector
        min_leaf_idx_children_ranks[i] = *new sdsl::rank_support_v<1>(min_leaf_idx);
    }
}

void RMQ::preprocess_leaf_rank_structures(Tree* tree) {
    preprocess_left_right_most_leaf(tree);
    populate_min_leaf_idx_in_side_trees_cps(tree);
    populate_min_leaf_idx_children(tree);
}

RMQ::RMQ(Tree* tree, lca_t* lca_prep_in) : lca_prep(lca_prep_in) {
    preprocess_leaf_rmqs(tree);
    preprocess_leaf_rank_structures(tree);
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
    
    int w_cp_idx = cp_indices[w->id];

    // Rank query gives index from bottom of centroid path, subtract from length of centroid path - 1 to get real index
    int size_of_centroid_path = depths[left_most_leaf[w->id]] - depths[cp_roots[w->id]->id];
    int deepest_node = size_of_centroid_path - 1 -
                       min_leaf_idx_in_side_trees_cps_ranks[w_cp_idx](left_most_leaf[v->id] - left_most_leaf[w->id]);

    gen_rmq_t* curr_rmq = cp_rmqs[w_cp_idx];
    return -curr_rmq->v[general_rmq(curr_rmq, deepest_node, depths[w->id] - depths[cp_roots[w->id]->id])];
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
    if (lca(lca_prep, v->id, w->children[0]->id) == w->children[0]->id) {
        return w->children[0];
    }
    
    // Otherwise, query children rank structure
    int child_rank = min_leaf_idx_children_ranks[w->id](left_most_leaf[v->id] - left_most_leaf[w->children[1]->id]);
    // Add 1 to rank to account for the fact that heaviest child is excluded
    return w->children[child_rank + 1];
}

int RMQ::max_weight_path(Tree::Node* v, Tree::Node* w) {
    return rmq(v, child_for_descendant(v, w));
}

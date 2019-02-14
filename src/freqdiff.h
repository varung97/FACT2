/*
 * freqdiff.h
 *
 *  Created on: 21 Oct 2015
 *      Author: Mesh
 */

#ifndef FREQDIFF_H_
#define FREQDIFF_H_

#include <iostream>
#include <queue>
#include <cassert>
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_set.hpp>
#include <boost/heap/fibonacci_heap.hpp>

#include "taxas_ranges.h"
#include "lca_preprocessing.h"
#include "utils.h"

#include "Tree.h"
#include "RMQ.h"

struct node_bitvec_t {
	Tree::Node* node;
	boost::dynamic_bitset<>* bitvector;

	node_bitvec_t() : node(NULL), bitvector(NULL) {}
	node_bitvec_t(Tree::Node* node, boost::dynamic_bitset<>* bitvector) : node(node), bitvector(bitvector) {}
	//~node_bitvec_t() { delete bitvector; } FIXME
};
bool operator < (const node_bitvec_t& bv1, const node_bitvec_t& bv2) {
	return *(bv1.bitvector) < *(bv2.bitvector);
}

int* start,* stop;
int* e,* m;
std::vector<Tree::Node*>* rsort_lists;

bool* marked;

Tree::Node** left,** right;
size_t* orig_pos_in_parent;

size_t* counter;
bool* BT;
int* leaf_p_index;

int* vleft,* vright;
int* pointer;
int* levels,* ids;
int* parent;
bool* exists;
Tree::Node** tree_nodes;
Tree::Node empty_node = Tree::Node(-1, -1, 0);


// calculate cluster weights using kn^2 method
void calc_w_kn2(std::vector<Tree*>& trees) {

	size_t n = Tree::get_taxas_num();

	int tot_int_nodes = 0;
	for (Tree* tree : trees) {
		tot_int_nodes += tree->get_nodes_num() - n - 1;
	}

	// Generate bit vectors
	node_bitvec_t* node_bitvecs = new node_bitvec_t[tot_int_nodes];
	int bitvc = 0;
	for (Tree* tree : trees) {
		taxas_ranges_t* tr = build_taxas_ranges(tree);

		size_t nodes_num = tree->get_nodes_num();
		for (size_t i = 1; i < nodes_num; i++) {
			Tree::Node* node = tree->get_node(i);
			if (!node->is_leaf()) {
				// fill bit vector for current node
				node_bitvec_t node_bitvec(node, new boost::dynamic_bitset<>(Tree::get_taxas_num()));
				for (int j = tr->intervals[node->id].start; j <= tr->intervals[node->id].end; j++) {
					node_bitvec.bitvector->set(tr->taxas[j]);
				}
				node_bitvecs[bitvc] = node_bitvec;
				bitvc++;
			}
		}

		delete tr;
	}

	// Sort bit vectors
	std::sort(node_bitvecs, node_bitvecs+tot_int_nodes);

	// count adjacent equal bitvectors and set weights
	for (int i = 0; i < bitvc; ) {
		int w = 1;
		while (i+w < bitvc && *(node_bitvecs[i].bitvector) == *(node_bitvecs[i+w].bitvector)) w++;
		for (int j = i; j < i+w; j++) {
			node_bitvecs[j].node->weight = w;
		}
		i += w;
	}

	for (int i = 0; i < bitvc; i++) {
		delete node_bitvecs[i].bitvector;
	}
	delete[] node_bitvecs;
}


void calc_w_k2n(std::vector<Tree*>& trees) {

	size_t n = Tree::get_taxas_num();

	int tot_int_nodes = 0;
	for (Tree* tree : trees) {
		tot_int_nodes += tree->get_nodes_num();
	}

	int* tree_start = new int[trees.size()];
	int* cache_taxa = new int[tot_int_nodes];
	int* cache_parent_ids = new int[tot_int_nodes];
	int* cache_size = new int[tot_int_nodes];
	int* cache_weight = new int[tot_int_nodes];
	std::fill(cache_weight, cache_weight+tot_int_nodes, 0);

	int pos = 0;
	for (size_t i = 0; i < trees.size(); i++) {
		Tree* tree = trees[i];
		tree_start[i] = pos;
		for (size_t j = 0; j < tree->get_nodes_num(); j++) {
			cache_taxa[pos] = tree->get_node(j)->taxa;
			if (j > 0) {
				cache_parent_ids[pos] = tree->get_node(j)->parent->id;
			}
			cache_size[pos] = tree->get_node(j)->size;
			pos++;
		}
	}

	int* e = new int[n];
	std::pair<int,int>* clusters = new std::pair<int,int>[n*2];
	std::pair<int,int>* H = new std::pair<int,int>[n*2];
	int* H2node = new int[n*2];
	std::pair<int,int>* minmax = new std::pair<int,int>[n*2];
	for (size_t i = 0; i < trees.size(); i++) {
		std::fill(H, H+Tree::get_taxas_num()*2, std::pair<int,int>(INT32_MAX, 0));

		Tree* tree1 = trees[i];

		// calculate e
		int count = 0;
		for (size_t j = 0; j < tree1->get_nodes_num(); j++) {
			if (tree1->get_node(j)->is_leaf()) {
				e[tree1->get_node(j)->taxa] = count++;
			}
			tree1->get_node(j)->weight++;
		}

		// calculate L and R for each cluster
		for (int j = tree1->get_nodes_num()-1; j > 0; j--) {
			Tree::Node* node = tree1->get_node(j);
			if (node->is_leaf()) {
				clusters[j].first = clusters[j].second = e[node->taxa];
			} else {
				clusters[j].first = clusters[node->children[0]->id].first;
				clusters[j].second = clusters[(*node->children.rbegin())->id].second;

				int pos;
				pos = node->pos_in_parent == node->parent->get_children_num()-1 ? clusters[j].first : clusters[j].second;
				H[pos] = clusters[j];
				H2node[pos] = tree_start[i]+j;
			}
		}

		// for each tree, see if cluster is in our ref tree, and if so increase its count by 1
		for (size_t j = i+1; j < trees.size(); j++) {
			Tree* tree2 = trees[j];

			std::fill(minmax, minmax+tree2->get_nodes_num(), std::pair<int,int>(INT32_MAX, 0));
			for (int k = tree2->get_nodes_num()-1; k > 0; k--) {
				int taxa = cache_taxa[tree_start[j]+k];
				int parent_id = cache_parent_ids[tree_start[j]+k];
				int size = cache_size[tree_start[j]+k];

				if (taxa != Tree::Node::NONE) {
					minmax[k].first = minmax[k].second = e[taxa];
				}

				if (minmax[parent_id].first > minmax[k].first) {
					minmax[parent_id].first = minmax[k].first;
				}
				if (minmax[parent_id].second < minmax[k].second) {
					minmax[parent_id].second = minmax[k].second;
				}

				if (taxa == Tree::Node::NONE && minmax[k].second-minmax[k].first+1 == size) {
					if (H[minmax[k].first] == minmax[k]) {
						cache_weight[H2node[minmax[k].first]]++;
						cache_weight[tree_start[j]+k]++;
					} else if (H[minmax[k].second] == minmax[k]) {
						cache_weight[H2node[minmax[k].second]]++;
						cache_weight[tree_start[j]+k]++;
					}
				}
			}
		}
	}

	for (size_t i = 0; i < trees.size(); i++) {
		Tree* tree = trees[i];
		for (size_t j = 0; j < tree->get_nodes_num(); j++) {
			tree->get_node(j)->weight += cache_weight[tree_start[i]+j];
		}
	}

	delete[] tree_start;
	delete[] cache_taxa;
	delete[] cache_parent_ids;
	delete[] cache_size;
	delete[] cache_weight;

	delete[] e;
	delete[] clusters;
	delete[] H;
	delete[] H2node;
	delete[] minmax;
}


// for node_id in T1:
// start[node_id] = smallest rank (i.e. leftmost leaf) in T2 among \Lambda(T1[node_id])
// stop[node_id] = biggest rank (i.e. rightmost leaf) in T2 among \Lambda(T1[node_id])
void compute_start_stop(Tree* tree1, Tree* tree2, int* t2_leaves_ranks) {
	// for each cluster in t1, find its leftmost and rightmost leaves in t2
	for (int i = tree1->get_nodes_num()-1; i >= 0; i--) {
		Tree::Node* node = tree1->get_node(i);
		if (node->is_leaf()) {
			start[i] = stop[i] = t2_leaves_ranks[node->taxa];
		} else {
			start[i] = INT32_MAX;
			stop[i] = 0;
			for (Tree::Node* child : node->children) {
				if (start[i] > start[child->id])
					start[i] = start[child->id];
				if (stop[i] < stop[child->id])
					stop[i] = stop[child->id];
			}
		}
	}
}


void filter_clusters_n2(Tree* tree1, Tree* tree2, taxas_ranges_t* t1_tr, taxas_ranges_t* t2_tr, lca_t* t2_lcas,
		bool* to_del) {

	int* t2_ranks = get_taxas_ranks(t2_tr);
	compute_start_stop(tree1, tree2, t2_ranks);
	delete t2_ranks;

	// mark clusters in t1 to be deleted if a heavier incompatible cluster is in t2
	for (size_t i = 1; i < tree1->get_nodes_num(); i++) {
		Tree::Node* node = tree1->get_node(i);
		if (node->is_leaf())
			continue;

		Tree::Node* leftmost_leaf = tree2->get_leaf(t2_tr->taxas[start[i]]);
		Tree::Node* rightmost_leaf = tree2->get_leaf(t2_tr->taxas[stop[i]]);
		int lcaX = lca(t2_lcas, leftmost_leaf->id, rightmost_leaf->id);
		marked[lcaX] = true;
		for (int j = t1_tr->intervals[i].start; j <= t1_tr->intervals[i].end; j++) {
			//mark all ancestors
			Tree::Node* curr = tree2->get_leaf(t1_tr->taxas[j]);
			while (!marked[curr->id]) {
				marked[curr->id] = true;
				curr = curr->parent;
			}
		}

		for (int j = t2_tr->intervals[lcaX].start; j <= t2_tr->intervals[lcaX].end; j++) {
			Tree::Node* leaf = tree2->get_leaf(t2_tr->taxas[j]);
			if (!marked[leaf->id]) { // leaf not in X, mark all ancestors up to lcaX
				Tree::Node* curr = leaf->parent;
				while (curr->id != lcaX) {
					if (marked[curr->id] && node->weight <= curr->weight) {
						// weight of node i in t1 is lower than an incompatible cluster in t2
						// mark for deletion
						to_del[i] = true;
						break;
					}
					curr = curr->parent;
				}
			}
			if (to_del[i])
				break;
		}

		marked[lcaX] = false;
		for (int j = t1_tr->intervals[i].start; j <= t1_tr->intervals[i].end; j++) {
			//mark all ancestors
			Tree::Node* curr = tree2->get_leaf(t1_tr->taxas[j]);
			while (marked[curr->id]) {
				marked[curr->id] = false;
				curr = curr->parent;
			}
		}
	}
}

void compute_m(Tree::Node* node, int* e, int* m, std::vector<Tree::Node*>* rsort_lists) {
	if (node->is_leaf()) {
		m[node->id] = e[node->taxa];
	}
	for (Tree::Node* child : node->children) {
		compute_m(child, e, m, rsort_lists);
		if (m[node->id] > m[child->id]) {
			m[node->id] = m[child->id];
		}
	}
	if (!node->is_root()) {
		rsort_lists[m[node->id]].push_back(node);
	}
}


// See Section 2.4 of paper [XXX]
void merge_trees(Tree* tree1, Tree* tree2, taxas_ranges_t* t1_tr, lca_t* t2_lcas) {
	//compute e
	for (size_t j = 0; j < Tree::get_taxas_num(); j++) {
		e[t1_tr->taxas[j]] = j;
	}
	std::fill(m, m+tree2->get_nodes_num(), INT32_MAX);
	for (size_t i = 0; i < Tree::get_taxas_num(); i++) rsort_lists[i].clear();

	compute_m(tree2->get_root(), e, m, rsort_lists);

	// sort tree2
	for (size_t i = 0; i < tree2->get_nodes_num(); i++) {
		tree2->get_node(i)->clear_children();
	}
	for (size_t i = 0; i < Tree::get_taxas_num(); i++) {
		for (auto it = rsort_lists[i].begin(); it != rsort_lists[i].end(); it++) {
			(*it)->parent->add_child(*it);
		}
	}

	taxas_ranges_t* t2_tr = build_taxas_ranges(tree2);
	int* t2_ranks = get_taxas_ranks(t2_tr);
	compute_start_stop(tree1, tree2, t2_ranks);
	delete[] t2_ranks;

	// calc x_left and x_right
	for (size_t i = 0; i < Tree::get_taxas_num(); i++) {
		Tree::Node* curr = tree2->get_leaf(i);
		Tree::Node* parent = curr->parent;
		while (parent != NULL && *(parent->children.begin()) == curr) {
			curr = parent;
			parent = curr->parent;
		}
		left[i] = curr;

		curr = tree2->get_leaf(i);
		parent = curr->parent;
		while (parent != NULL && *(parent->children.rbegin()) == curr) {
			curr = parent;
			parent = curr->parent;
		}
		right[i] = curr;
	}

	for (size_t i = 0; i < tree2->get_nodes_num(); i++) {
		orig_pos_in_parent[i] = tree2->get_node(i)->pos_in_parent;
	}

	for (int i = tree1->get_nodes_num()-1; i >= 1; i--) {
		Tree::Node* a = tree2->get_leaf(t2_tr->taxas[start[i]]);
		Tree::Node* b = tree2->get_leaf(t2_tr->taxas[stop[i]]);
		if (a == b) continue; // tree1->node i is a leaf

		Tree::Node* ru = tree2->get_node(lca(t2_lcas, a->id, b->id));

		Tree::Node* a_left = left[a->taxa];
		Tree::Node* b_right = right[b->taxa];

		size_t du_pos = (a_left->depth > ru->depth) ? orig_pos_in_parent[a_left->id] : 0;
		size_t eu_pos = (b_right->depth > ru->depth) ? orig_pos_in_parent[b_right->id] : ru->get_children_num()-1;
		if (du_pos == 0 && eu_pos == ru->get_children_num()-1) continue;

		Tree::Node* newnode = tree2->add_node();
		newnode->weight = tree1->get_node(i)->weight;
		newnode->size = tree1->get_node(i)->size;
		for (size_t j = du_pos; j <= eu_pos; j++) {
			if (ru->children[j] != NULL) {
				newnode->add_child(ru->children[j]);
				// lazy child deletion
				ru->null_child(j);
			}
		}
		ru->set_child(newnode, du_pos);
	}
	tree2->fix_tree();

	delete t2_tr;
}

// Returns true if node has a marked leaf in subtree
// Populates is_lca, true if node is lca of marked leaves
bool mark_lcas(Tree::Node* node, std::vector<bool>* is_marked_leaf, std::vector<bool>* is_lca) {
	if (is_marked_leaf->at(node->id)) {
		is_lca->at(node->id) = true;
		return true;
	}

	int num_children_with_marked_leaves = 0;

	for (Tree::Node* child : node->children) {
		if (mark_lcas(child, is_marked_leaf, is_lca)) {
			num_children_with_marked_leaves++;
		}
	}

	if (num_children_with_marked_leaves >= 2) {
		is_lca->at(node->id) = true;
	}

	return num_children_with_marked_leaves > 0;
}

// Builds the subtree of tree restricted to some leaves given marked lcas
Tree::Node* build_restricted_subtree(Tree* tree, Tree::Node* root, std::vector<bool>* is_lca, std::vector<Tree::Node*>* assoc_nodes) {
	bool root_is_lca = is_lca->at(root->id);

	if (root->is_leaf()) {
		if (root_is_lca) {
			// If root is a marked leaf then create new node
			// Leaves are labelled with 1, this is only useful for trees with a single leaf
			return assoc_nodes->at(root->id) = tree->add_node(root->taxa, 1);
		} else {
			// Else no subtree
			return NULL;
		}
	}

	// If there is only one child with non-NULL subtree, then we can just return that subtree
	// If more than one child have non-NULL subtrees, then root is an lca and we need to connect each of these subtrees to the node
	Tree::Node* node = root_is_lca ? tree->add_node() : NULL;
	for (Tree::Node* child : root->children) {
		Tree::Node* child_subtree = build_restricted_subtree(tree, child, is_lca, assoc_nodes);

		if (child_subtree != NULL) {
			if (root_is_lca) {
				node->add_child(child_subtree);
			} else {
				assoc_nodes->at(root->id) = child_subtree;
				return child_subtree;
			}
		}
	}

	// At this point, either root is an lca or has no subtree
	if (root_is_lca) {
		return assoc_nodes->at(root->id) = node;
	} else {
		return NULL;
	}
}

// Builds the subtree of tree restricted to marked_leaves
// Also populates the associated node for each node in tree, defined as shallowest descendant of the node in the restricted subtree
Tree* restricted_subtree(Tree* tree, std::vector<int>* marked_taxa, std::vector<Tree::Node*>* assoc_nodes) {
	if (marked_taxa->empty()) return NULL;

	std::vector<bool> is_lca(tree->get_nodes_num(), false);
	std::vector<bool> is_marked_leaf(tree->get_nodes_num(), false);

	for (int taxon : *marked_taxa) {
		is_marked_leaf[tree->get_leaf(taxon)->id] = true;
	}

	mark_lcas(tree->get_root(), &is_marked_leaf, &is_lca);

	// Build new tree with only lcas
	Tree* new_tree = new Tree(std::count(is_lca.begin(), is_lca.end(), true));
	Tree::Node* root = build_restricted_subtree(new_tree, tree->get_root(), &is_lca, assoc_nodes);
	new_tree->fix_tree(root);

	return new_tree;
}

// Stably counting sorts elements of a vector along axis given by index
std::vector<std::vector<int>> counting_sort(std::vector<std::vector<int>>& vector, int index) {
	// Asummed to be non-negative
	int max_label = 0;
	for (int i = 0; i < vector.size(); i++) {
		max_label = max_label > vector[i][index] ? max_label : vector[i][index];
	}

	// Each bucket will hold a list of elements
	std::vector<std::vector<std::vector<int>>> buckets(max_label + 1);

	// Bucketize
	for (int i = 0; i < vector.size(); i++) {
		buckets[vector[i][index]].push_back(vector[i]);
	}

	std::vector<std::vector<int>> sorted_vector;

	// Reassemble
	for (int i = 0; i < buckets.size(); i++) {
		for (int j = 0; j < buckets[i].size(); j++) {
			sorted_vector.push_back(buckets[i][j]);
		}
	}

	return sorted_vector;
}

// Labels tree using the labels of the associated nodes in the restricted subtrees
void label_trees_using_assoc_nodes(std::vector<Tree*> trees, std::vector<std::vector<Tree::Node*>>& assoc_nodes_vector1, std::vector<std::vector<Tree::Node*>>& assoc_nodes_vector2) {
	// A label pair looks like (left_label, right_label, tree_index, node_id)
	std::vector<std::vector<int>> label_pairs;

	for (int tree_index = 0; tree_index < trees.size(); tree_index++) {
		for (int node_id = 0; node_id < trees[tree_index]->get_nodes_num(); node_id++) {
			label_pairs.push_back({assoc_nodes_vector1[tree_index][node_id]->label, assoc_nodes_vector2[tree_index][node_id]->label, tree_index, node_id});
		}
	}

	// Stable counting sort along each of the labels
	label_pairs = counting_sort(label_pairs, 1);
	label_pairs = counting_sort(label_pairs, 0);

	// Identical nodes stores (tree_index, node_id) of identical nodes
	std::vector<std::vector<int>> identical_nodes;
	int label = 1;

	for (int i = 0; i <= label_pairs.size(); i++) {
		// If identical nodes is not empty and this pair is different from the previous one then label the identical nodes
		// Also do so if i == label_pairs.size() since there are no more pairs to be checked
		if (i == label_pairs.size() ||
			(!identical_nodes.empty() &&
			 (label_pairs[i - 1][0] != label_pairs[i][0] || label_pairs[i - 1][1] != label_pairs[i][1])
			 )
			) {

			for (int j = 0; j < identical_nodes.size(); j++) {
				trees[identical_nodes[j][0]]->get_node(identical_nodes[j][1])->label = label;
			}

			label++;
			identical_nodes.clear();
		}

		if (i < label_pairs.size()) {
			identical_nodes.push_back({label_pairs[i][2], label_pairs[i][3]});
		}
	}
}

// Divide marked_leaves into 2 parts and recursively label
void label_trees_helper(std::vector<Tree*>& trees, std::vector<int>& marked_taxa) {
	if (marked_taxa.size() == 1) {
		return;
	}

	std::vector<int> marked_taxa1, marked_taxa2;

	int i = 0;
	for (; i < marked_taxa.size() / 2; i++) {
		marked_taxa1.push_back(marked_taxa[i]);
	}
	for (; i < marked_taxa.size(); i++) {
		marked_taxa2.push_back(marked_taxa[i]);
	}

	std::vector<Tree*> trees1, trees2;
	std::vector<std::vector<Tree::Node*>> assoc_nodes_vector1, assoc_nodes_vector2;

	for (Tree* tree : trees) {
		// Default assoc node is the empty node, representing no marked taxa in leafset
		std::vector<Tree::Node*> assoc_nodes1(tree->get_nodes_num(), &empty_node);
		std::vector<Tree::Node*> assoc_nodes2(tree->get_nodes_num(), &empty_node);

		trees1.push_back(restricted_subtree(tree, &marked_taxa1, &assoc_nodes1));
		trees2.push_back(restricted_subtree(tree, &marked_taxa2, &assoc_nodes2));

		assoc_nodes_vector1.push_back(assoc_nodes1);
		assoc_nodes_vector2.push_back(assoc_nodes2);
	}

	label_trees_helper(trees1, marked_taxa1);
	label_trees_helper(trees2, marked_taxa2);

	label_trees_using_assoc_nodes(trees, assoc_nodes_vector1, assoc_nodes_vector2);
}

// Labels nodes in the trees, giving identical labels to nodes with identical leafsets
void label_trees(std::vector<Tree*>& trees) {
	std::vector<int> marked_taxa;
	for (int i = 0; i < Tree::get_taxas_num(); i++) {
		marked_taxa.push_back(i);
	}

	label_trees_helper(trees, marked_taxa);
}

void calc_w_knlogn(std::vector<Tree*>& trees) {
	label_trees(trees);

	// A label looks like (label, tree_index, node_id)
	std::vector<std::vector<int>> labels;

	for (int tree_index = 0; tree_index < trees.size(); tree_index++) {
		for (int node_id = 0; node_id < trees[tree_index]->get_nodes_num(); node_id++) {
			labels.push_back({trees[tree_index]->get_node(node_id)->label, tree_index, node_id});
		}
	}

	labels = counting_sort(labels, 0);

	// Identical nodes stores (tree_index, node_id) of identical nodes
	std::vector<std::vector<int>> identical_nodes;

	for (int i = 0; i <= labels.size(); i++) {
		// If identical nodes is not empty and this label is different from the previous one then compute freq for the identical nodes
		// Also do so if i == labels.size() since there are no more labels to be checked
		if (i == labels.size() ||
			(!identical_nodes.empty() && labels[i - 1][0] != labels[i][0])
			) {

			for (int j = 0; j < identical_nodes.size(); j++) {
				trees[identical_nodes[j][0]]->get_node(identical_nodes[j][1])->weight = identical_nodes.size();
			}

			identical_nodes.clear();
		}

		if (i < labels.size()) {
			identical_nodes.push_back({labels[i][1], labels[i][2]});
		}
	}
}

// Removes a node from the heap and erases its heap handle
void remove_from_heap_if_inside(boost::heap::fibonacci_heap<int>* heap, std::unordered_map<int, boost::heap::fibonacci_heap<int>::handle_type>* heap_handles, Tree::Node* v) {
	if (heap_handles->find(v->id) == heap_handles->end()) {
		return;
	}

	// Increase-key then delete-max
	heap->increase(heap_handles->at(v->id), heap->top() + 1);
	heap->pop();
	heap_handles->erase(v->id);
}

// Puts the maximum weight on the path from v to w into the heap, where the path excludes v and w
// If v was already in the heap, deletes old entry
void push_max_weight_into_heap(boost::heap::fibonacci_heap<int>* heap, std::unordered_map<int, boost::heap::fibonacci_heap<int>::handle_type>* heap_handles, Tree::Node* v, Tree::Node* w, RMQ* rmq) {
	if (v->parent == w) {
		// Path is empty
		return;
	}

	remove_from_heap_if_inside(heap, heap_handles, v);
	heap_handles->insert({v->id, heap->push(rmq->max_weight_path(v->parent, w))});
}

// Returns rootsOfSubtrees(\leafset(T_A))
boost::unordered_set<Tree::Node*>* filter_clusters_nlogn_helper(Tree::Node* root_T_A, Tree* T_B, RMQ* rmq_T_B, bool* to_del_T_A) {
	Tree::Node* curr = root_T_A;
	std::vector<boost::unordered_set<Tree::Node*>> lower_boundaries = std::vector<boost::unordered_set<Tree::Node*>>();

	// recursively call filter_clusters on side trees and build up lower boundaries
	while (!curr->is_leaf()) {
		boost::unordered_set<Tree::Node*> lower_boundary = boost::unordered_set<Tree::Node*>();

		for (int i = 1; i < curr->get_children_num(); i++) {
			lower_boundary.merge(*filter_clusters_nlogn_helper(curr->children[i], T_B, rmq_T_B, to_del_T_A));
		}

		lower_boundaries.push_back(lower_boundary);
		curr = curr->children[0];
	}

	// Was built top to bottom, but more convenient in the other direction
	std::reverse(lower_boundaries.begin(), lower_boundaries.end());

	// First lca, leaf in T_B labelled by leaf in T_A
	Tree::Node* l_i = T_B->get_leaf(curr->taxa);
	l_i->counter = 0;
	l_i->parent->counter = 1;
	int lower_boundary_counter = -1; // Start at -1 since no lower boundary for leaf

	// Stores rootsOfSubtrees
	boost::unordered_set<Tree::Node*>* roots = new boost::unordered_set<Tree::Node*>({l_i});

	// Stores incompatible nodes
	boost::heap::fibonacci_heap<int>* incompatible = new boost::heap::fibonacci_heap<int>();

	// Stores handles to nodes in the heap
	std::unordered_map<int, boost::heap::fibonacci_heap<int>::handle_type>* heap_handles = new std::unordered_map<int, boost::heap::fibonacci_heap<int>::handle_type>;

	while (curr != root_T_A) {
		curr = curr->parent;
		lower_boundary_counter++;

		Tree::Node* prev_l_i = l_i;
		for (Tree::Node* v : lower_boundaries[lower_boundary_counter]) {
			l_i = T_B->get_node(lca(rmq_T_B->lca_prep, l_i->id, v->id));
		}

		// Add nodes from prev_l_i to l_i
		// Include prev_l_i if incompatible with previous centroid path node
		if (prev_l_i->counter == prev_l_i->get_children_num()) {
			push_max_weight_into_heap(incompatible, heap_handles, prev_l_i, l_i, rmq_T_B);
		} else {
			push_max_weight_into_heap(incompatible, heap_handles, prev_l_i->children[0], l_i, rmq_T_B);
		}

		for (Tree::Node* v : lower_boundaries[lower_boundary_counter]) {
			v->parent->counter++;
			roots->insert(v);

			Tree::Node* prev_v = v;
			v = v->parent;

			while (v->counter == v->get_children_num() && v != T_B->get_root()) {
				v->parent->counter++;

				roots->insert(v);
				for (Tree::Node* child : v->children) {
					roots->erase(child);
					remove_from_heap_if_inside(incompatible, heap_handles, child);
				}

				prev_v = v;
				v = v->parent;
			}

			// Add max weight from v until l_i to incompatible
			if (v != l_i->parent) {
				push_max_weight_into_heap(incompatible, heap_handles, prev_v, l_i, rmq_T_B);
			}
		}

		int max_weight = heap_handles->size() == 0 ? 0 : incompatible->top();
		if (curr->weight <= max_weight) {
			to_del_T_A[curr->id] = true;
		}
	}

	// Need to erase counters that would affect other calls to filter_clusters_nlogn_helper
	for (Tree::Node* root : *roots) {
		root->parent->counter = 0;
	}

	return roots;
}

void filter_clusters_nlogn(Tree* T_A, Tree* T_B, RMQ* rmq_T_B, bool* to_del_T_A) {
	filter_clusters_nlogn_helper(T_A->get_root(), T_B, rmq_T_B, to_del_T_A);

	// Need to zero out counters since they might be left over from previous calls
	for (int i = 0; i < T_A->get_nodes_num(); i++) {
		T_A->get_node(i)->counter = 0;
	}

	for (int i = 0; i < T_B->get_nodes_num(); i++) {
		T_B->get_node(i)->counter = 0;
	}
}

Tree* freqdiff(std::vector<Tree*>& trees, bool centroid_paths, bool onlyw) {
	start = new int[Tree::get_taxas_num()*2];
	stop = new int[Tree::get_taxas_num()*2];
	e = new int[Tree::get_taxas_num()];
	m = new int[Tree::get_taxas_num()*2];
	rsort_lists = new std::vector<Tree::Node*>[Tree::get_taxas_num()];

	marked = new bool[Tree::get_taxas_num()*2];
	std::fill(marked, marked+Tree::get_taxas_num()*2, false);
	bool* to_del_t = new bool[Tree::get_taxas_num()*2];
	bool* to_del_ti = new bool[Tree::get_taxas_num()*2];

	left = new Tree::Node*[Tree::get_taxas_num()];
	right = new Tree::Node*[Tree::get_taxas_num()];
	orig_pos_in_parent = new size_t[Tree::get_taxas_num()*2];

	BT = new bool[Tree::get_taxas_num()*2];
	counter = new size_t[Tree::get_taxas_num()*2];
	leaf_p_index = new int[Tree::get_taxas_num()];

	// tree contract
	vleft = new int[Tree::get_taxas_num()*2];
	vright = new int[Tree::get_taxas_num()*2];
	pointer = new int[Tree::get_taxas_num()*2];
	levels = new int[Tree::get_taxas_num()*2];
	ids = new int[Tree::get_taxas_num()*2];
	parent = new int[Tree::get_taxas_num()*2];
	exists = new bool[Tree::get_taxas_num()*2];
	tree_nodes = new Tree::Node*[Tree::get_taxas_num()*2];

	for (size_t i = 0; i < trees.size(); i++) {
		trees[i]->reorder();
	}

	calc_w_knlogn(trees);

	if (onlyw) return NULL;

	lca_t** lca_preps = new lca_t*[trees.size()];
	for (size_t i = 0; i < trees.size(); i++) {
		lca_preps[i] = lca_preprocess(trees[i]);
	}

	RMQ** rmqs = NULL;
	if (centroid_paths) {
		rmqs = new RMQ*[trees.size()];
		for (size_t i = 0; i < trees.size(); i++) {
			rmqs[i] = new RMQ(trees[i], lca_preps[i]);
		}
	}

	Tree* T = new Tree(trees[0]);
	for (size_t i = 1; i < trees.size(); i++) {
		Tree* Ti = new Tree(trees[i]);
		taxas_ranges_t* tr_Ti = build_taxas_ranges(Ti);

		taxas_ranges_t* tr_T = build_taxas_ranges(T);
		lca_t* lca_T = lca_preprocess(T);

		if (centroid_paths) {
			T->reorder();
			Ti->reorder();
		}

		// filter clusters
		std::fill(to_del_ti, to_del_ti+Ti->get_nodes_num(), false);
		if (centroid_paths) {
			RMQ* rmq_T = new RMQ(T, lca_T);
			filter_clusters_nlogn(Ti, T, rmq_T, to_del_ti);
		} else {
			filter_clusters_n2(Ti, T, tr_Ti, tr_T, lca_T, to_del_ti);
		}

		std::fill(to_del_t, to_del_t+T->get_nodes_num(), false);
		if (centroid_paths) {
			filter_clusters_nlogn(T, Ti, rmqs[i], to_del_t);
		} else {
			filter_clusters_n2(T, Ti, tr_T, tr_Ti, lca_preps[i], to_del_t);
		}

		Ti->delete_nodes(to_del_ti);
		T->delete_nodes(to_del_t);

		delete lca_T;
		delete tr_T;
		delete tr_Ti;

		lca_T = lca_preprocess(T);
		tr_Ti = build_taxas_ranges(Ti);

		merge_trees(Ti, T, tr_Ti, lca_T);

		delete lca_T;
		delete tr_Ti;
		delete Ti;
	}

	std::fill(to_del_t, to_del_t+T->get_nodes_num(), false);
	taxas_ranges_t* tr_T = build_taxas_ranges(T);
	for (size_t i = 0; i < trees.size(); i++) {
		if (centroid_paths) {
			T->reorder();
			trees[i]->reorder();
			filter_clusters_nlogn(T, trees[i], rmqs[i], to_del_t);
		} else {
			taxas_ranges_t* tr_Ti = build_taxas_ranges(trees[i]);
			filter_clusters_n2(T, trees[i], tr_T, tr_Ti, lca_preps[i], to_del_t);
			delete tr_Ti;
		}

		delete trees[i];
	}
	T->delete_nodes(to_del_t);
	delete tr_T;

	delete[] start;
	delete[] stop;
	delete[] e;
	delete[] m;
	delete[] rsort_lists;

	delete[] marked;
	delete[] to_del_ti;
	delete[] to_del_t;

	delete[] orig_pos_in_parent;
	delete[] left;
	delete[] right;

	for (size_t i = 0; i < trees.size(); i++) {
		delete lca_preps[i];
	}
	delete[] lca_preps;

	if (centroid_paths) {
		for (size_t i = 0; i < trees.size(); i++) {
			delete rmqs[i];
		}
		delete[] rmqs;
	}

	return T;
}


#endif /* FREQDIFF_H_ */

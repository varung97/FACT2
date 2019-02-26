#include <vector>
#include <cassert>

#include "SimpleTree.h"

SimpleTree::SimpleTree(size_t nodes_num_hint) {
	if (nodes_num_hint > 0) {
		nodes.reserve(nodes_num_hint);
	}
}

SimpleTree::~SimpleTree() {
	for (SimpleNode* node : nodes) {
		delete node;
	}
}

size_t SimpleTree::get_nodes_num() {
	return nodes.size();
}

SimpleTree::SimpleNode* SimpleTree::get_node(int i) {
	return nodes[i];
}

SimpleTree::SimpleNode* SimpleTree::get_root() {
    return root;
}

void SimpleTree::set_root(SimpleTree::SimpleNode *in_root) {
    root = in_root;
}

SimpleTree::SimpleNode* SimpleTree::add_node(int taxa, int label) {
	SimpleTree::SimpleNode* newnode = new SimpleTree::SimpleNode(get_nodes_num(), taxa, label);
	nodes.push_back(newnode);
	return newnode;
}

SimpleTree::SimpleNode::SimpleNode(int id, int taxa, int label) : parent(NULL), id(id), taxa(taxa), label(label) {}

void SimpleTree::SimpleNode::add_child(SimpleNode* child) {
	child->parent = this;
	children.push_back(child);
}

bool SimpleTree::SimpleNode::is_leaf() {
	return taxa != NONE;
}

bool SimpleTree::SimpleNode::is_root() {
	return parent == NULL;
}

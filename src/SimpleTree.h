#ifndef SIMPLETREE_H_
#define SIMPLETREE_H_

#include <vector>

class SimpleTree {
public:
	SimpleTree(size_t nodes_num_hint = 0);
	virtual ~SimpleTree();

	class SimpleNode;

	SimpleNode* get_node(int i);
	size_t get_nodes_num();
	SimpleNode* add_node(int taxa = -1, int label = 1);

	SimpleNode* root;

private:
	std::vector<SimpleNode*> nodes;
};

class SimpleTree::SimpleNode {
public:
	static const int NONE = -1;

	std::vector<SimpleNode*> children;

	int id;
	int taxa;
    int label;

	SimpleNode(int id, int taxa, int label);

	void add_child(SimpleNode* child);
	bool is_leaf();
};

#endif

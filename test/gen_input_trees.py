import random, sys

class Node:
    def __init__(self, taxa=None):
        self.parent = None
        self.children = set()
        self.taxa = taxa

    def is_leaf(self):
        return self.taxa is not None

    def to_newick(self):
        if self.is_leaf():
            return str(self.taxa)

        return '(' + ','.join([str(child) for child in self.children]) + ')'

    def __repr__(self):
        return self.to_newick()


class Tree:
    def __init__(self, newick=None, root=None):
        self.nodes = set([root])
        self.leaves = set()
        self.edges = set()

        if newick is not None:
            self.from_newick(newick)
        else:
            self.root = root

    def add_child(self, parent, child):
        parent.children.add(child)
        child.parent = parent

        self.edges.add((parent, child))
        self.nodes.add(child)
        if child.is_leaf():
            self.leaves.add(child)

    def remove_edge(self, parent, child):
        parent.children.remove(child)
        child.parent = None

        self.edges.remove((parent, child))
        self.nodes.remove(child)
        if child.is_leaf():
            self.leaves.remove(child)

    def delete_node(self, node):
        """Doesn't keep track of edges"""
        self.nodes.remove(node)

        for child in node.children:
            child.parent = node.parent

        if node != self.root:
            node.parent.children.remove(node)
            node.parent.children |= node.children
        else:
            assert len(node.children) == 1
            self.root = list(node.children)[0]

    def to_newick(self):
        return str(self.root)

    def __repr__(self):
        return self.to_newick()

    def copy_helper(self, root, copy_tree):
        if root.is_leaf():
            copy_leaf = Node(root.taxa)
            copy_tree.nodes.add(copy_leaf)
            copy_tree.leaves.add(copy_leaf)
            return copy_leaf

        copy_root = Node()
        for child in root.children:
            copy_child = self.copy_helper(child, copy_tree)
            copy_child.parent = copy_root
            copy_root.children.add(copy_child)

        copy_tree.nodes.add(copy_root)
        return copy_root

    def copy(self):
        copy_tree = Tree(None)
        copy_tree.nodes.remove(None)
        copy_root = self.copy_helper(self.root, copy_tree)
        copy_tree.root = copy_root
        return copy_tree

    def nodes_in_subtree(self, node):
        nodes = set([node])
        for child in node.children:
            nodes |= self.nodes_in_subtree(child)
        return nodes


def delete_ghost(tree_with_ghost):
    ghost = tree_with_ghost.root
    ghost_child = tuple(ghost.children)[0]

    tree_with_ghost.edges.remove((ghost, ghost_child))
    tree_with_ghost.nodes.remove(ghost)
    tree_with_ghost.root = ghost_child

    return tree_with_ghost


def uniform_model(n):
    tree_with_ghost = Tree(newick=None, root=Node())
    tree_with_ghost.add_child(tree_with_ghost.root, Node(0))

    for i in range(1, n):
        leaf = Node(i)
        leaf_parent = Node()
        edge = random.choice(tuple(tree_with_ghost.edges))
        parent = edge[0]
        other_child = edge[1]

        tree_with_ghost.remove_edge(parent, other_child)
        tree_with_ghost.add_child(parent, leaf_parent)
        tree_with_ghost.add_child(leaf_parent, other_child)
        tree_with_ghost.add_child(leaf_parent, leaf)

    return delete_ghost(tree_with_ghost)


def scenario1(k, n):
    tree = uniform_model(n)

    for node in tree.nodes - set([tree.root]) - tree.leaves:
        if random.random() < 0.2:
            tree.delete_node(node)

    trees = []

    for i in range(k):
        copy = tree.copy()

        for j in range(int(0.05 * n)):
            u = random.choice(tuple(copy.nodes - set([copy.root])))
            v = random.choice(tuple(copy.nodes - copy.leaves - copy.nodes_in_subtree(u)))

            if v == u.parent:
                continue

            u.parent.children.remove(u)
            if len(u.parent.children) == 1:
                # u.parent would only have one child left
                copy.delete_node(u.parent)

            u.parent = v
            v.children.add(u)

        trees.append(copy)

    return trees


def scenario2(k, n):
    trees = []

    for i in range(k):
        tree = uniform_model(n)

        for node in tree.nodes - set([tree.root]) - tree.leaves:
            if random.random() < 0.2:
                tree.delete_node(node)

        trees.append(tree)

    return trees


def gen_nex_file(trees, n, filename):
    with open(filename, 'w') as f:
        f.write('BEGIN TAXA;\n')
        f.write('    TAXLABELS ' + ' '.join(map(str, range(n))) + ';\n')
        f.write('END;\n\n')
        f.write('BEGIN TREES;\n')
        f.write('    necessary_otherwise_doesnt_work;\n')
        f.write('    necessary_otherwise_doesnt_work;\n')

        for i, tree in enumerate(trees):
            f.write('    TREE tree' + str(i) + ' = ' + str(tree) + '\n')

        f.write('END;\n')

if __name__ == '__main__':
    n = int(sys.argv[1])
    k = int(sys.argv[2])
    sys.setrecursionlimit(5000)
    gen_nex_file(scenario1(k, n), n, sys.argv[3])

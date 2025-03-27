import ete3
import numpy as np
import matplotlib.pyplot as plt
import random
from .utils import (NEWICK_TREE, one_primitive_step, is_valid_tree, prune_by_branch_length, prune_by_root_to_tip)


def prune_tree_PSFA(name_of_the_tree, threshold=90, longest_to_average=9, name_of_output="psfa_tree.nwk"):
    if not is_valid_tree(name_of_the_tree):
        raise TypeError("The input file is not a valid phylogenetic tree. Please make sure that name_of_the_tree is a valid .nwk tree.")
    
    tree = ete3.Tree(name_of_the_tree, format=1, quoted_node_names=True)
    leaf_names = tree.get_leaf_names()
    num_leaves = len(leaf_names)
    original_num_leaves = num_leaves

    percent_left = 100 - threshold

    while percent_left > 0:
        tree, PRUNE, percent_left = one_primitive_step(tree, original_num_leaves, percent_left, longest_to_average)
        if not PRUNE:
            break

    tree.write(outfile=name_of_output, format=1)
    leaf_names = tree.get_leaf_names()
    num_leaves = len(leaf_names)
    return tree, 100 * num_leaves / original_num_leaves


def prune_tree_CPA(name_of_the_tree, root_to_node_ratio=0.1, min_num_of_roots=15, M_n=0, threshold=90, beta=20, name_of_output="cpa_tree.nwk", show_plot=False, show_pruned_tips=False, safe_tips=[]):
    if not is_valid_tree(name_of_the_tree):
        raise TypeError("The input file is not a valid phylogenetic tree. Please make sure that name_of_the_tree is a valid .nwk tree.")
    
    tree = ete3.Tree(name_of_the_tree, format=1, quoted_node_names=True)
    leaf_names = tree.get_leaf_names()
    num_leaves = len(leaf_names)
    original_num_leaves = num_leaves

    if M_n == 0:
        M_n_int = num_leaves
    else:
        M_n_int = M_n(num_leaves)

    def hang(tree1, root1):
      dists = []
      for leaf in tree1.iter_leaves():
          distance = tree1.get_distance(root1, leaf)
          dists.append(distance)
      return dists

    non_leaf_nodes = []
    num_non_leaf_nodes = 0
    for node in tree.traverse():
        if not node.is_leaf():
            non_leaf_nodes.append(node)
            num_non_leaf_nodes += 1


    num_roots = min(max(int(root_to_node_ratio * num_non_leaf_nodes), min_num_of_roots), num_non_leaf_nodes)
    roots = random.sample(non_leaf_nodes, num_roots)
    midpoint_root = tree.get_midpoint_outgroup()
    if not midpoint_root.is_leaf():
        roots.append(midpoint_root)

    best_p_v = float("inf")
    for root in roots:
        data = hang(tree, root)
        bin_edges = np.linspace(min(data), max(data), M_n_int + 1)
        hist, _ = np.histogram(data, bins=bin_edges)

        cumulative_freq = np.cumsum(hist)
        total_freq = cumulative_freq[-1]
        threshold_index = np.argmax(cumulative_freq >= (threshold / 100) * total_freq)
        stop = M_n_int - 1

        if threshold_index == M_n_int - 1:
            stop = threshold_index
        elif threshold_index == M_n_int - 2:
            if hist[threshold_index + 1] < (beta / 100) * hist[threshold_index]:
                stop = threshold_index
        else:
            for i in range(threshold_index, M_n_int - 1):
                if hist[i + 1] < (beta / 100) * hist[i]:
                    stop = i
                    break

        p_v = bin_edges[stop+1]
        if p_v < best_p_v:
            best_p_v = p_v
            best_root = root
            best_stop = stop

    if show_plot:
        data = hang(tree, best_root)
        bin_edges = np.linspace(min(data), max(data), M_n_int + 1)
        hist, _ = np.histogram(data, bins=bin_edges)

        plt.hist(data, bins=bin_edges, edgecolor='black', align='mid')
        patches = plt.gca().patches

        for i, patch in enumerate(patches):
            if i > best_stop:
                patch.set_facecolor('red')

        plt.xlabel('Distance from root')
        plt.ylabel('Frequency')
        plt.title("CPA")
        plt.show()

    leaves_to_keep = []
    leaves_to_prune = []
    for leaf in tree.iter_leaves():
        distance_to_node = best_root.get_distance(leaf)
        if (distance_to_node <= best_p_v) or (leaf.name in safe_tips):
            leaves_to_keep.append(leaf)
        else:
            leaves_to_prune.append(leaf.name)

    tree.prune(leaves_to_keep, preserve_branch_length=True)

    tree.write(outfile=name_of_output, format = 1)
    leaf_names = tree.get_leaf_names()
    num_leaves = len(leaf_names)

    if show_pruned_tips:
        return tree, 100 * num_leaves / original_num_leaves, leaves_to_prune
    return tree, 100 * num_leaves / original_num_leaves


def prune_tree_IQR(name_of_the_tree, threshold = 90, name_of_output="iqr_tree.nwk"):
    if not is_valid_tree(name_of_the_tree):
     raise TypeError("The input file is not a valid phylogenetic tree. Please make sure that name_of_the_tree is a valid .nwk tree.")

    tree = ete3.Tree(name_of_the_tree, format=1, quoted_node_names=True)
    original_num_leaves = len(tree.get_leaves())
    tree_1, percent_left = prune_by_branch_length(tree, threshold)
    tree_2 = prune_by_root_to_tip(tree_1, percent_left, original_num_leaves)
    num_leaves = len(tree_2.get_leaves())
    tree_2.standardize()
    tree_2.write(outfile=name_of_output, format=1)
    return tree_2, 100 * num_leaves / original_num_leaves

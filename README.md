# Phylogenetic tree pruner
A package for pruning phylogenetic trees. The final package is not yet available, however you can download the latest test version from TestPyPi: https://test.pypi.org/project/treepruner/ .

## How to use
After installing the package use following import: <br>

**`import treepruner`** <br>
or <br>
 **`from treepruner import (prune_tree_PSFA, prune_tree_CPA, prune_tree_IQR)`**


## The following algorithms can be used to prune exceptionally long branches in phylogenetic trees. 

- **Primitive Straight-Forward Approach (PSFA)**
    - Identifies the longest edge and calculates the sizes of the two components resulting from its removal

    - If the size of the smaller component is within the tolerance range and the edge was excessively long, removes it
    
    - Have to define criteria for considering an edge as "excessively long", denoted by `longest_to_average` in the code

    - Have to specify the tolerance range, denoted by `threshold`

    - You can call this function by `treepruner.prune_tree_PSFA`

    - The function has four arguments in the following order:
        
        - `name_of_the_tree` is the name of your tree file, which should be in the same folder as your code. For example `name_of_the_tree = ST2332.nwk`. This argument must be a valid phylogenetic tree, but it can be .nwk, .newick or .tree format.

        - `threshold` is the above mentioned tolerance range. It is the maximum percentage of leaves that the algorithm is allowed to remove. The default value is `threshold = 10`.

        - `longest_to_average` determines how many times longer should a branch be than the average, to consider it too long. The default value is `longest_to_average = 9`.

        - `name_of_output` is the desired name of the output file, which will be downloaded into the same folder as the input file. The default value is `name_of_output = psfa_tree.nwk`. You can also set this to .newick, if you prefer that format.

- **Circular Pruning Algorithm (CPA)**
    - This is a more complex and (most of the time) better algorithm than the previous one, but the runtime is also longer.

    - The goal is to prune the trees to make them roughly circular

    - Chooses a "root" and removes leaves that are beyond a specified distance

    - Two conditions must be satisfied:
    
      - Do not remove too many leaves: the minimum percent of remaining leaves is denoted by `alpha`
        
      - It should not be possible to significantly reduce the radius by removing only a few nodes (adjust `beta` and `M_n`)

    - Randomly selects some roots and prunes to meet both conditions (adjust `root_to_node_ratio` and `min_num_of_roots`)

    - Among all the pruned trees, picks the best one

    - The function has eight arguments in the following order:

        - `name_of_the_tree` is the same as before.

        - `root_to_node_ratio` determines how many nodes does the algorithm try out as "roots". By increasing this parameter, the precision and the runtime becomes greater too. The default value is `root_to_node_ratio = 0.1` meaning that the algorithm will try out a randomly selected 10% of the non-leaf nodes as "roots". Apart from these random nodes the algorithm will always try the "midpoint root" of the tree, as it usually gives a good result. The parameter must be between 0 and 1, although the algorithm can handle the other cases as well (as long as it is a number).

        - `min_num_of_roots` is the minimum number of nodes that the algorithm tries out as "roots". This parameter only matters if the `root_to_node_ratio` chooses less nodes than `min_num_of_roots`. It must be less than the number of non-leaf nodes, although the algorithm can handle the other cases as well (as long as it is a number). The default value is `min_num_of_nodes = 15`.

        - `M_n` is a function, which you should define before calling the algorithm if you intend to use this parameter. During the CPA the leaves are being placed in brackets based on their distances from the "root". The algorithm will delete the leaves in the last couple of brackets if certain conditions are met. `M_n` determines the number of brackets: the total number of brackets will be `M_n(num_leaves)`, where `num_leaves` is the number of leaves on the tree. The default value is `M_n = 0`, which means that the user did not give any function to the argument. In this case `M_n(num_leaves) = num_leaves`. Assigning anything to `M_n`other than a function will result in an error message.

        - `alpha` is the minimum percentage of the leaves that the algorithm must keep. The default value is `alpha = 90`.

        - `beta` determines where to prune the tree. After placing the leaves in the brackets, the first couple of brackets will be safe, until the total number of leaves in them reaches `alpha*num_leaves`. When the CPA finds this threshold, it is ready for pruning. The algorithm will start the pruning if a bracket contains `beta` percent less leaves than the previous one. After this point the CPA will delete every leaf. The default value is `beta = 25`.

        - `name_of_output` is the same as before. The default value is `name_of_output=cpa_tree.nwk`.

        - `show_plot` is a boolean parameter. By default it is set to `show_plot = False`, meaning that after running the code it will not show any plot. By setting this value to `True`, the brackets will be shown to the user. These brackets correspond to the "best root" (the center of the outscribed circle of the tree). On the x axis you can see the distances from the root, while on the y axis you can see how many leaves were placed in each bracket. Blue indicates the preserved brackets, while red indicates the ones which have been deleted.

- **IQR algorithm**
  - Abbreviations:

    - IQR: the upper fence of the Inter Quartile Range of a vector: Q3 + 3 \* (Q3 - Q1) = extreme outlier threshold

    - R2T: root-to-tip distance on a midpoint rooted tree

  - 1<sup>st</sup> part: Pruning the unrooted tree based on the upper fence of IQR branch lengths.

    - 1.  Calculating the minimum number of tips for each branch (unrooted bifurcating tree, both directions are looked up, than the least tip number is chosen to represent the branch).

    - 2.  Calculating the IQR for branch lengths.

    - 3.  Excluding those extreme outlier branches containing less than 0.05 (or a given) proportion of the tips.

  - 2<sup>nd</sup> part: Pruning the midpoint rooted tree based on root-to-tip distances. Excluding tips based on the R2T IQRs, while iteratively midpoint rooting and excluding the top gretaest extreme outlier.

    - 4.  Midpoint rooting the tree (after pruning with method described in the 1<sup>st</sup> part).

    - 5.  Calculating the IQR for root-to-tip distances.

    - 6.  Excluding the most extreme outlier tip based on the IQR for root-to-tip distances.

    - 7.  Repeating point 4. to 6. till there are no more extreme outlier IQR tip is found.

    - 8.  Unrooting the pruned tree. 

  - The function has three arguments in the following order:

    - `name_of_the_tree` is the same as before

    - `tipprop` is a number between 0 and 1. The algorithm only deletes a branch, if the proportion of the leaves on one side of that branch is at most this value. By default, `tipprop = 1`, meaning that the algorithm will always delete the branch, if the other two requirements are satisfied (that is we do not cut more than the threshold and the branch is indeed an outlier) .

    - `threshold` is the proportion of leaves that must be kept on the tree. The default value is `threshold = 0.9` meaning that the algorithm will not cut more than 10% of the leaves.

    - `name_of_output` is the same as before.

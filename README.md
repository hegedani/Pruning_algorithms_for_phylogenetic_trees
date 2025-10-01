# Phylogenetic tree pruner
A package for pruning phylogenetic trees. The package is also available at [PyPI](https://pypi.org/project/treepruner/).

## What are phylogenetic trees and why should we prune them?
Phylogenetic trees describe evolutionary relationships among taxa based on genetic data. However, errors due to low sequencing quality, contamination, or misclassification can distort these trees.

### Why prune phylogenetic trees?

- To remove erroneous taxa that might skew evolutionary insights.

- To ensure the tree represents biological reality more accurately.

- To simplify analysis by focusing on reliable data.

# How to use the package

## Prerequisites:

`python3.12` with packages:
  - `ete3` (3.1.3)
  - `numpy` (2.3.3)
  - `matplotlib` (3.10.6)

The script possible works with other versions of these packages, too.

## Standalone version

    treepruner.py

## In a singularity container

  1. download singularity image from docker: https://hub.docker.com/r/mesti90/treepruner/
     
       singularity pull -F mesti90-treepruner.1.1.sif docker://mesti90/treepruner:1.1
     
  3. ```singularity exec mesti90-treepruner.1.1.sif treepruner @args```

## CLI arguments

	treepruner.py [-h] -i INPUT_TREE [-o OUTPUT_TREE] [--stats-file STATS_FILE]
						 [-m {PSFA,CPA,IQR} [{PSFA,CPA,IQR} ...]] [--log LOG] [--output_dir OUTPUT_DIR]
						 [--psfa-output PSFA_OUTPUT] [--psfa-threshold PSFA_THRESHOLD]
						 [--psfa-longest-to-average PSFA_LONGEST_TO_AVERAGE] [--cpa-output CPA_OUTPUT]
						 [--cpa-threshold CPA_THRESHOLD] [--cpa-root-to-node-ratio CPA_ROOT_TO_NODE_RATIO]
						 [--cpa-min-num-of-roots CPA_MIN_NUM_OF_ROOTS] [--cpa-M-n CPA_M_N]
						 [--cpa-beta CPA_BETA] [--cpa-radius-ratio CPA_RADIUS_RATIO]
						 [--cpa-safe-tips [CPA_SAFE_TIPS [CPA_SAFE_TIPS ...]]] [--cpa-show-plot]
						 [--cpa-show-pruned-tips] [--iqr-output IQR_OUTPUT]
						 [--iqr-threshold IQR_THRESHOLD]

	Prune trees using one or more methods:
		  CPA  - Circular Pruning Algorithm
		  IQR  - Interquartile Range approach
		  PSFA - Primitive Straight-Forward approach

	optional arguments:
	  -h, --help            show this help message and exit
	  -i INPUT_TREE, --input-tree INPUT_TREE
							One or more input tree files (Newick format) (default: None)
	  -o OUTPUT_TREE, --output-tree OUTPUT_TREE
							Output tree file name. Final output will be: <output_dir>/<output_tree> (default: output_tree.nwk)
	  --stats-file STATS_FILE
							Filename for CSV summary of pruning stats (saved in --output_dir) (default: pruning_summary.csv)
	  -m {PSFA,CPA,IQR} [{PSFA,CPA,IQR} ...], --methods {PSFA,CPA,IQR} [{PSFA,CPA,IQR} ...]
							List of methods to apply (space-separated, e.g., 'CPA IQR PSFA') (default: ['PSFA', 'CPA'])
	  --log LOG             Log file path (default: treepruner.log)
	  --output_dir OUTPUT_DIR
							Directory to save all output trees and reports (default: results)

	PSFA options:
	  --psfa-output PSFA_OUTPUT
	  --psfa-threshold PSFA_THRESHOLD
	  --psfa-longest-to-average PSFA_LONGEST_TO_AVERAGE

	CPA options:
	  --cpa-output CPA_OUTPUT
	  --cpa-threshold CPA_THRESHOLD
	  --cpa-root-to-node-ratio CPA_ROOT_TO_NODE_RATIO
	  --cpa-min-num-of-roots CPA_MIN_NUM_OF_ROOTS
	  --cpa-M-n CPA_M_N
	  --cpa-beta CPA_BETA
	  --cpa-radius-ratio CPA_RADIUS_RATIO
	  --cpa-safe-tips [CPA_SAFE_TIPS [CPA_SAFE_TIPS ...]]
	  --cpa-show-plot
	  --cpa-show-pruned-tips

	IQR options:
	  --iqr-output IQR_OUTPUT
	  --iqr-threshold IQR_THRESHOLD


## As a package

Install from pypi

```python
import treepruner 
#or
from treepruner import prune_tree_PSFA, prune_tree_CPA, prune_tree_IQR
```

## The following algorithms can be used to prune exceptionally long branches in phylogenetic trees. 

- **Primitive Straight-Forward Approach (PSFA)**
    - Identifies the longest edge and calculates the sizes of the two components resulting from its removal.

    - If the size of the smaller component is within the tolerance range and the edge was excessively long, removes it.
    
    - You have to define criteria for considering an edge as "excessively long".

    - You also have to specify the tolerance range.

    - You can call this function by `treepruner.prune_tree_PSFA`. The function's output consists of two parts. The first part is the pruned tree in `ete3.Tree` format, while the second part is a number indicating the percentage of leaves that remain in the tree after pruning.

    - The function has four arguments in the following order:
        
        - `name_of_the_tree` (str): the name of your tree file, which should be in the same folder as your code. For example `name_of_the_tree = ST2332.nwk`. This argument must be a valid phylogenetic tree, but it can be .nwk, .newick or .tree format.

        - `threshold` (float, default=90): the minimum percentage of the leaves that the algorithm must keep.

        - `longest_to_average` (float, default=9): determines how many times longer should a branch be than the average, to consider it too long.

        - `name_of_output` (str, default="psfa_tree.nwk"): the desired name of the output file, which will be saved into the same folder as the input file. You can also set this to .newick, if you prefer that format.

- **Circular Pruning Algorithm (CPA)**
    - This is a more complex and (most of the time) better algorithm than the previous one, but the runtime is also longer.

    - The goal is to prune the trees to make them roughly circular.

    - Chooses a "root" and removes leaves that are beyond a specified distance.

    - Two conditions must be satisfied:
    
      - Do not remove too many leaves.
        
      - It should not be possible to significantly reduce the radius by removing only a few nodes.

    - Randomly selects some roots and prunes to meet both conditions.

    - Among all the pruned trees, picks the best one.

    - You can call this function by `treepruner.prune_tree_CPA`. The function's output is the same as above.

    - The function has ten arguments in the following order:

        - `name_of_the_tree` (str): the same as above.

        - `root_to_node_ratio` (float, default=0.1): determines how many nodes does the algorithm try out as "roots". By increasing this parameter, the precision and the runtime becomes greater too. The default value means that the algorithm will try out a randomly selected 10% of the non-leaf nodes as "roots". Apart from these random nodes the algorithm will always try the "midpoint root" of the tree, as it usually gives a good result. The parameter must be between 0 and 1, although the algorithm can handle the other cases as well (as long as it is a number).

        - `min_num_of_roots` (int, default=15): the minimum number of nodes that the algorithm tries out as "roots". This parameter only matters if the `root_to_node_ratio` chooses less nodes than `min_num_of_roots`. It must be less than the number of non-leaf nodes, although the algorithm can handle the other cases as well (as long as it is an integer).

        - `M_n` (function, default=0): a function, which you should define before calling the algorithm if you intend to use this parameter. During the CPA the leaves are being placed in brackets based on their distances from the "root". The algorithm will delete the leaves in the last couple of brackets if certain conditions are met. `M_n` determines the number of brackets: the total number of brackets will be `M_n(num_leaves)`, where `num_leaves` is the number of leaves on the tree. The default value means that the user did not give any function to the argument. In this case `M_n(num_leaves) = num_leaves`. Assigning anything to `M_n`other than a function will result in an error message.

        - `threshold` (float, default=90): the same as above.

        - `beta` (float, default=25): determines where to prune the tree. After placing the leaves in the brackets, the first couple of brackets will be safe, until the total number of leaves in them reaches `alpha*num_leaves`. When the CPA finds this threshold, it is ready for pruning. The algorithm will start the pruning if a bracket contains `beta` percent less leaves than the previous one. After this point the CPA will delete every leaf.
     
        - `radius_ratio` (float, default=0): controls whether pruning occurs based on the tree's radius. If the longest root-to-tip distance in the pruned tree is not at least `radius_ratio` percent smaller than in the original tree (measured from the same root), no pruning is performed.

        - `safe_tips` (list, default=[]): a list where you can input the names of the tips that you do not want to be cut off. After the algorithm chose the best root and started the pruning, it checks which leaves are in the list and leaves them on the tree. The default value means that there are no safe leaves. It is important that the names in the list must be strings.

        - `show_plot` (bool, default=False): the default value means that after running the code it will not show any plot. By setting this value to `True`, the brackets will be shown to the user. These brackets correspond to the "best root" (the center of the outscribed circle of the tree). On the x axis you can see the distances from the root, while on the y axis you can see how many leaves were placed in each bracket. Blue indicates the preserved brackets, while red indicates the ones which have been deleted.

        - `show_pruned_tips` (bool, default=False): responsible for indicating which leaves have been cut from the tree. By setting this value to `True`, the function will return a list that contains the names of the pruned tips, in addition to the default return values. This list will be the third return value.

        - `name_of_output` (str, default="cpa_tree.nwk"): the same as above.

- **IQR algorithm**
  - Abbreviations:

    - IQR: the upper fence of the Inter Quartile Range of a vector: Q3 + 3 \* (Q3 - Q1) = extreme outlier threshold.

    - R2T: root-to-tip distance on a midpoint rooted tree.

  - First part of the algorithm: Pruning the unrooted tree based on the upper fence of IQR branch lengths.

    - Calculating the minimum number of tips for each branch (unrooted bifurcating tree, both directions are looked up, than the least tip number is chosen to represent the branch).

    - Calculating the IQR for branch lengths.

    - Excluding those extreme outlier branches containing less than 0.05 (or a given) proportion of the tips.

  - Second part of the algorithm: Pruning the midpoint rooted tree based on root-to-tip distances. Excluding tips based on the R2T IQRs, while iteratively midpoint rooting and excluding the top greatest extreme outlier.

    - Midpoint rooting the tree (after pruning with method described in the first part).

    - Calculating the IQR for root-to-tip distances.

    - Excluding the most extreme outlier tip based on the IQR for root-to-tip distances.

    - Repeating the previous three steps until there are no more extreme outlier IQR tip is found.

    - Unrooting the pruned tree. 

  - You can call this function by `treepruner.prune_tree_IQR`. The function's output is the same as above.

  - The function has three arguments in the following order:

    - `name_of_the_tree` (str): the same as above.

    - `threshold` (float, default=90): the same as above.

    - `name_of_output` (str, default="iqr_tree.nwk"): the same as above.

## Example usage
The code below shows a combination of our algorithms. In some cases, these combinations provide better results, but be careful when experimenting with them.

```python
import treepruner

def f(x):
    return x**2

# Step 1: Prune using PSFA and save intermediate tree
first_output, first_percentage = treepruner.prune_tree_PSFA("ST2332.nwk", name_of_output="pre-pruned_tree.nwk")

# Step 2: Further prune using CPA
final_output, final_percentage = treepruner.prune_tree_CPA("pre-pruned_tree.nwk", threshold=95, M_n=f)

print(final_percentage)
```

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

The script possibly works with other versions of these packages, too.

## Standalone version

    bin/treepruner.py

## In a singularity container

  1. Download singularity image from docker: https://hub.docker.com/r/mesti90/treepruner/
     
       `singularity pull -F mesti90-treepruner.1.1.sif docker://mesti90/treepruner:1.1`
     
  2. `singularity exec mesti90-treepruner.1.1.sif treepruner @args`

## CLI arguments

	treepruner.py [-h] -i INPUT_TREE [-o OUTPUT_TREE] [--stats-file STATS_FILE]
             [-m {PSFA,CPA,IQR,PSFA_CPA_PSFA} [{PSFA,CPA,IQR,PSFA_CPA_PSFA} ...]] [--log LOG] [--output_dir OUTPUT_DIR]
             [--psfa-output PSFA_OUTPUT] [--psfa-threshold PSFA_THRESHOLD]
             [--psfa-longest-to-average PSFA_LONGEST_TO_AVERAGE] [--cpa-output CPA_OUTPUT]
             [--cpa-threshold CPA_THRESHOLD] [--cpa-root-to-node-ratio CPA_ROOT_TO_NODE_RATIO]
             [--cpa-min-num-of-roots CPA_MIN_NUM_OF_ROOTS] [--cpa-M-n CPA_M_N]
             [--cpa-beta CPA_BETA] [--cpa-radius-ratio CPA_RADIUS_RATIO]
             [--cpa-safe-tips [CPA_SAFE_TIPS [CPA_SAFE_TIPS ...]]] [--cpa-show-plot]
             [--cpa-show-pruned-tips] [--iqr-output IQR_OUTPUT]
             [--iqr-threshold IQR_THRESHOLD] [--psfacpa-psfa1-threshold PSFACPA_PSFA1_THRESHOLD]
             [--psfacpa-cpa-threshold PSFACPA_CPA_THRESHOLD] [--psfacpa-psfa2-threshold PSFACPA_PSFA2_THRESHOLD]
             [--psfacpa-longest-to-average PSFACPA_LONGEST_TO_AVERAGE]
             [--psfacpa-root-to-node-ratio PSFACPA_ROOT_TO_NODE_RATIO]
             [--psfacpa-min-num-of-roots PSFACPA_MIN_NUM_OF_ROOTS] [--psfacpa-M-n PSFACPA_M_N]
             [--psfacpa-beta PSFACPA_BETA] [--psfacpa-radius-ratio PSFACPA_RADIUS_RATIO]
             [--psfacpa-safe-tips [PSFACPA_SAFE_TIPS [PSFACPA_SAFE_TIPS ...]]]
             [--psfacpa-show-plot]

	Prune trees using one or more methods:
      CPA           - Circular Pruning Algorithm
      IQR           - Interquartile Range approach
      PSFA          - Primitive Straight-Forward approach
      PSFA_CPA_PSFA - Sequential pruning (PSFA -> CPA -> PSFA)

	optional arguments:
	    -h, --help            show this help message and exit
	    -i INPUT_TREE, --input-tree INPUT_TREE
	              One or more input tree files (Newick format) (default: None)
	    -o OUTPUT_TREE, --output-tree OUTPUT_TREE
	              Output tree file name. Final output will be: <output_dir>/<output_tree> (default: pruned_tree.nwk)
	    --stats-file STATS_FILE
	              Filename for CSV summary of pruning stats (saved in --output_dir) (default: pruning_summary.csv)
	    -m {PSFA,CPA,IQR,PSFA_CPA_PSFA} [{PSFA,CPA,IQR,PSFA_CPA_PSFA} ...], --methods ...
	              List of methods to apply (space-separated, e.g., 'CPA IQR PSFA_CPA_PSFA') (default: ['PSFA', 'CPA'])
	    --log LOG             Log file path (default: treepruner.log)
	    --output_dir OUTPUT_DIR
	              Directory to save all output trees and reports (default: .)



*(Method-specific CLI options correspond directly to the python arguments documented below.)*

## As a package

Install from pypi

```python
import treepruner 
#or
from treepruner import prune_tree_PSFA, prune_tree_CPA, prune_tree_IQR, prune_tree_PSFA_CPA_PSFA
```

## The following algorithms can be used to prune exceptionally long branches in phylogenetic trees. 

- **Primitive Straight-Forward Approach (PSFA)**
	- Identifies the longest edge and calculates the sizes of the two components resulting from its removal.
	- If the size of the smaller component is within the tolerance range and the edge was excessively long, removes it.
	- You have to define criteria for considering an edge as "excessively long".
	- You also have to specify the tolerance range.
	- You can call this function via `treepruner.prune_tree_PSFA`. The function directly returns a `float` indicating the percentage of leaves that remain in the tree after pruning, and automatically saves the pruned tree to your disk.
	- **Arguments:**
		- `name_of_the_tree` (str): the name of your tree file. This argument must be a valid phylogenetic tree (.nwk, .newick or .tree format).
		- `threshold` (int, default=90): the minimum percentage of the leaves that the algorithm must keep.
		- `longest_to_average` (int, default=9): determines how many times longer should a branch be than the average, to consider it too long.
		- `name_of_output` (str, default="psfa_tree.nwk"): the desired name of the output file.
		- `save_output` (bool, default=True): If set to `True`, the pruned tree will be saved to your disk using `name_of_output`. Set to `False` to only return the statistics.

- **Circular Pruning Algorithm (CPA)**
    - This is a more complex and (most of the time) better algorithm than the previous one, but the runtime is also longer.
	- The goal is to prune the trees to make them roughly circular.
	- Chooses a "root" and removes leaves that are beyond a specified distance.
	- Two conditions must be satisfied:
		- Do not remove too many leaves.
		- It should not be possible to significantly reduce the radius by removing only a few nodes.
	- Randomly selects some roots and prunes to meet both conditions.
	- Among all the pruned trees, picks the best one.
	- You can call this function via `treepruner.prune_tree_CPA`. The function returns a tuple of two values: `(percentage_remaining (float), pruned_tips (list))`, and automatically saves the pruned tree to your disk.
	- **Arguments:**
		- `name_of_the_tree` (str): the same as above.
		- `root_to_node_ratio` (float, default=0.1): determines how many nodes does the algorithm try out as "roots". By increasing this parameter, the precision and the runtime becomes greater too. Must be between 0 and 1.
		- `min_num_of_roots` (int, default=15): the minimum number of nodes that the algorithm tries out as "roots".
		- `M_n` (int or callable, default=0): Determines the number of distance brackets the leaves are placed into. If a function is provided, the total number of brackets will be `M_n(num_leaves)`. If an integer is provided, it acts as a hard limit on the number of brackets. The default value (0) sets the number of brackets equal to the number of leaves.
		- `threshold` (int, default=90): the same as above.
		- `beta` (float, default=20): determines where to prune the tree based on leaf distribution density. 
		- `radius_ratio` (float, default=0): controls whether pruning occurs based on the tree's radius reduction.
		- `safe_tips` (list, default=[]): a list of tip names (strings) that you do not want to be cut off under any circumstances.
		- `show_plot` (bool, default=False): if `True`, displays a histogram showing the leaf brackets, distances, and indicating which brackets have been deleted (in red).
		- `show_pruned_tips` (bool, default=False): *(Deprecated parameter kept for legacy configs; the function now always returns the pruned tips list).*
		- `name_of_output` (str, default="cpa_tree.nwk"): the same as above.
		- `save_output` (bool, default=True): the same as above.


- **IQR algorithm**
	- **Abbreviations:**
		- IQR: the upper fence of the Inter Quartile Range of a vector: Q3 + 3 * (Q3 - Q1) = extreme outlier threshold.
		- R2T: root-to-tip distance on a midpoint rooted tree.
	- **Part 1:** Pruning the unrooted tree based on the upper fence of IQR branch lengths.
	- **Part 2:** Pruning the midpoint rooted tree based on root-to-tip distances. Iteratively midpoint roots the tree and excludes the greatest extreme outlier until no more outliers are found.
	- You can call this function via `treepruner.prune_tree_IQR`. It returns a `float` indicating the percentage of leaves that remain, and automatically saves the pruned tree to your disk.
	- **Arguments:**
		- `name_of_the_tree` (str): the same as above.
		- `threshold` (int, default=90): the same as above.
		- `name_of_output` (str, default="iqr_tree.nwk"): the same as above.
		- `save_output` (bool, default=True): the same as above.


- **Sequential Pruning (PSFA -> CPA -> PSFA)**
	- A specialized combination method that links the three algorithms into a powerful sequential pipeline to deeply clean noisy data.
	- The three threshold arguments dictate the **maximum cumulative percentage of the ORIGINAL tree** that may be removed after each respective stage.
	- You can call this function via `treepruner.prune_tree_PSFA_CPA_PSFA`. It returns a `float` indicating the percentage of leaves remaining in the final tree.
	- **Arguments:**
		- `name_of_the_tree` (str): the same as above.
		- `psfa1_threshold` (int, default=5): max percentage of original tree allowed to be pruned after Step 1.
		- `cpa_threshold` (int, default=10): max percentage of original tree allowed to be pruned after Step 2.
		- `psfa2_threshold` (int, default=15): max percentage of original tree allowed to be pruned after Step 3.
		- *(It also accepts all configuration parameters from PSFA and CPA, such as `longest_to_average`, `M_n`, `beta`, `safe_tips`, etc.)*
		- `name_of_output` (str, default="psfa_cpa_psfa_tree.nwk"): the same as above.
		- `save_output` (bool, default=True): the same as above.


## Example usage
The code below shows how to use our algorithms and how to combine them sequentially for superior results.

```python
import treepruner

# Example custom callable for M_n
def f(num_leaves):
    return num_leaves ** 2

# ---------------------------------------------------------
# Example A: Manual step-by-step combination
# ---------------------------------------------------------

# Step 1: Prune using PSFA (Returns just the float percentage)
first_percentage = treepruner.prune_tree_PSFA(
    "ST2332.nwk", 
    threshold=90,
    name_of_output="pre-pruned_tree.nwk"
)

# Step 2: Further prune the intermediate file using CPA 
# (CPA returns the float percentage AND the list of pruned tips)
final_percentage, pruned_tips = treepruner.prune_tree_CPA(
    "pre-pruned_tree.nwk", 
    threshold=95, 
    M_n=f,
    name_of_output="manual_final.nwk"
)

print(f"Manual Pipeline Retained: {final_percentage}%")
print(f"CPA pruned tips: {len(pruned_tips)}")


# ---------------------------------------------------------
# Example B: Using the built-in Sequential Pipeline
# ---------------------------------------------------------

# This automatically handles the temp files and rolling budgets!
combined_percentage = treepruner.prune_tree_PSFA_CPA_PSFA(
    "ST2332.nwk", 
    psfa1_threshold=3,
    cpa_threshold=10,
    psfa2_threshold=15,
    M_n=f,
    name_of_output="auto_combined_final.nwk"
)

print(f"Auto Combined Pipeline Retained: {combined_percentage}%")
```

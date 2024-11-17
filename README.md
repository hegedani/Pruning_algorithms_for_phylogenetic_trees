The listed algorithms can be used to prune exceptionally long branches in phylogenetic trees. 

- **Primitive Straight-Forward Approach (PSFA)**
    - Identify the longest edge and calculate the sizes of the two components resulting from its removal

    - If the size of the smaller component is within the tolerance range and the edge was excessively long, remove it
    
    - Define criteria for considering an edge as "too long", denoted by `longest_to_average` in the code

    - Specify the tolerance range, denoted by `threshold`

- **Circular Pruning Algorithm (CPA)**
    - The goal is to prune the trees to make them roughly circular

    - Choose a root and remove points that are beyond a specified distance

    - Two conditions must be satisfied:
    
      - Do not remove too many leaves, minimum percent of remaining leaves is denoted by `alpha`
        
      - It should not be possible to significantly reduce the radius by removing only a few nodes (adjust `beta` and `M_n`)

    - Randomly select some roots and prune to meet both conditions (adjust `root_to_node_ratio` and `min_num_of_roots`)

    - Among all the pruned trees, pick the best one

- **CPA + PSFA**
    - After pruning with CPA, excessively long branches may remain
        
    - These are pruned by a post-run PSFA
    
    - The parameters are the same as above

- **PSFA + CPA + PFSA**
    - CPA starts with a solution closer to the optimal

    - Remaining long branches are pruned by a final PSFA
        
    - Shorter runtime
 
    - Define a threshold for bounding the first PSFA (adjust `first_psfa`)
 
    - Apart from `first_psfa`, the parameters are the same as above


## Usage

The New_Algorithms.ipynb file consists of five cells. The first just installs the ete3 package. The other cells are used to run four different variants of the above described algorithms.

In each cell, the usage of the specific algorithm is described in comments. 



## Publications

The algorithms are yet to be published.

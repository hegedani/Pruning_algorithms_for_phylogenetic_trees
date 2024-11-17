TreeShrink is an algorithm for detecting abnormally long branches in one or more phylogenetic trees. 

- **Inputs**: 
    - One or more phylogenetic trees with branch lengths. If more than one, the trees should be on overlapping sets of species (though missing data are allowed). 
    - Optional: a number `k` ≤ the total number of species.
    - Optional: a selection of one of the three implemented algorithms for outlier detection.
    - Optional: a false positive tolerance rate, α
    - Optional: a set of alignments, from which, shrunk sequences will be removed
- **Outputs**:
    - The removing list: the final suggested list of species to be removed from each input tree, computed based on the selected statistical test. 
    - The shrunk trees: the input trees with the suggested leaves removed. 
    - The filtered alignments: the input alignments (if provided) with suggested leaves removed. 
    
Note that the tree diameter is the maximum distance between any two leaves of the tree. When multiple trees are available (e.g., gene trees), the statistical tests can use the information from all genes to decide what branches are too long. 

#### Publications:

The latest version of TreeShrink is described in:

## Usage: 
After installing TreeShrink, you can type 

~~~bash
run_treeshrink.py -h
~~~

to learn about all the options. 

## Examples:
The TreeShrink package comes with several testing trees that can be found in [test_data.zip](test_data.zip). If you downloaded TreeShrink from Github, you should have ```test_data.zip``` in your ```TreeShrink``` folder. If you installed using Anaconda, you should download [test_data.zip](https://github.com/uym2/TreeShrink/blob/master/test_data.zip) to your machine. Unzip ```test_data.zip``` before running the following examples.

### The simplest use case
The following command will produce the shrunk trees and the corresponding list of the species that were removed at false positive error rate `α = 0.05` (default)

~~~bash
run_treeshrink.py  -t test_data/mm10.trees
~~~

After running the command, the program will generate the folder `test_data/mm10_treeshrink/`, inside which you will find the shrunk trees (`output.trees`) and the removed species (`output.txt`). You should see 10 trees in `output.trees` corresponding to 10 trees of the input file `mm10.trees`. Accordingly, there are 10 lines in `output.txt`, each shows the list of species that were removed in the corresponding tree (empty lines indicating that the tree has no species removed). 
If you wish to customize the outputs, use ```-o``` to change the output folder and ```-O``` to change and the output prefix.

### Adjusting α threshold
The α threshold can be adjusted using ```-q``` option. 
You can run TreeShrink with multiple α thresholds, as follow

~~~bash
run_treeshrink.py  -t test_data/mm10.trees -q "0.05 0.10" -o test_data/mm10_treeshrink_multi -O shrunk
~~~
 
The program will generate the folder `test_data/mm10_treeshrink_multi/` inside which there are two sets of shrunk trees and removing sets at α = 0.05 and α = 0.10.
 
As TreeShrink is running, it will output a bunch of messages to the console. We suggest saving these in a text log file:

~~~bash
run_treeshrink.py  -t test_data/mm10.trees > test_data/mm10.trees.treeshrinklog.txt
~~~
 
### Modes


~~~

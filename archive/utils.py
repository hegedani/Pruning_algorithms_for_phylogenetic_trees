import ete3
import numpy as np

NEWICK_TREE = """
((Coli_005251:0.005251993,(((Coli_015270:0.004658564,((Coli_015261:0.001648792,Coli_015268:0.002248997)
0.961:0.000597463,(Coli_015255:0.001646538,Coli_015254:0.001048311)0.261:0.000000005)0.997:0.001496873)0.874:0.000295958
,((Coli_007818:0.028156151,(Coli_007872:0.058757615,Coli_007843:0.099775688)1.000:0.032269015)1.000:0.026010078,(Coli_01
5273:0.001649472,Coli_015274:0.001347749)1.000:0.002261343)0.756:0.000134552)0.819:0.000000005,(((Coli_005900:0.00733119
9,((Coli_015341:0.006776630,(Coli_015330:0.004091904,(Coli_015328:0.003213051,Coli_015331:0.004524901)0.786:0.000733106)
0.864:0.000351166)1.000:0.005292432,(((Coli_015248:0.002078626,(Coli_015246:0.001038116,Coli_015243:0.000443051)0.886:0.
000000005)1.000:0.007584228,(Coli_015266:0.000441077,Coli_015264:0.001640945)1.000:0.003874398)0.956:0.000834697,(((Coli
_015343:0.007647088,(Coli_015334:0.002313442,Coli_015333:0.003193376)1.000:0.006144972)1.000:0.002792010,(Coli_015299:0.0
07244310,(Coli_015344:0.010620236,Coli_015292:0.008859214)0.149:0.000493781)0.952:0.000816629)1.000:0.003401559,(Coli_015
308:0.007731858,((Coli_006020:0.071505078,(Coli_006013:0.000227477,((Coli_006019:0.196202996,(Coli_006017:0.216202649,C
oli_006278:0.082826221)0.853:0.003919287)0.944:0.003375374,Coli_006016:0.159116645)0.925:0.002729205)0.291:0.001287395)1
.000:0.006257664,Coli_015329:0.005284177)0.735:0.000109068)1.000:0.001991166)0.265:0.000280329)0.997:0.001386409)0.702:
0.000146719)1.000:0.001785951,(Coli_015250:0.008551488,(((Coli_015293:0.031270675,((((Coli_015361:0.006170696,(Coli_015
230:0.004818736,((Coli_015225:0.003905217,((Coli_015216:0.000299330,(Coli_015218:0.000299167,Coli_015221:0.000599619)0.
961:0.000599020)0.270:0.000000005,((Coli_015220:0.000598964,Coli_015219:0.000149559)0.963:0.000599387,Coli_015215:0.000
899039)0.760:0.000149657)1.000:0.002102428)0.976:0.000897219,(Coli_015224:0.000449523,(Coli_015223:0.000299295,Coli_015
226:0.001198303)0.909:0.000448912)1.000:0.003299664)0.479:0.000000005)0.785:0.000298563)0.998:0.001950818,(Coli_015232:
0.004808030,((Coli_015238:0.001952556,(Coli_015249:0.003610709,(Coli_015247:0.001651408,(Coli_015245:0.001050274,(Coli_
015241:0.000599427,Coli_015242:0.000149790)0.752:0.000149785)0.910:0.000449468)0.998:0.001499345)0.729:0.000146378)1.00
0:0.003756046,(Coli_015234:0.003604926,(Coli_015235:0.000600498,Coli_015236:0.000000005)1.000:0.003303984)0.967:0.00074
6912)0.948:0.000597877)0.986:0.001197720)1.000:0.011454235,((Coli_001884:0.014956149,(Coli_001885:0.006780479,Coli_0018
89:0.005319577)1.000:0.010304965)1.000:0.005794317,(Coli_002360:0.001484486,Coli_001888:0.003410645)0.998:0.001954848)0
.992:0.001466440)0.994:0.001627639,(Coli_007272:0.017991800,Coli_011224:0.010679395)0.977:0.001612403)0.981:0.001257617
)0.381:0.000440652,(Coli_010951:0.010967406,Coli_015276:0.014149388)0.865:0.000000005)1.000:0.005879011,(Coli_015262:0
.011740996,(Coli_015277:0.005246584,Coli_015285:0.008542483)0.998:0.002511585)1.000:0.004624530)1.000:0.007778136)0.77
3:0.000000005)1.000:0.002545511,(((Coli_015231:0.005252447,(Coli_015214:0.001950115,Coli_015229:0.004802605)0.984:0.00
1044002)0.997:0.001496991,(Coli_015233:0.003908978,((Coli_015239:0.004048347,(Coli_015228:0.000312499,(Coli_015227:0.0
00910960,Coli_015222:0.001483070)0.438:0.000149618)1.000:0.004043955)0.979:0.001045561,(Coli_015237:0.005105379,Coli_0
15244:0.005554506)0.999:0.000000005)0.959:0.000000005)0.994:0.001345453)0.916:0.000447793,(((Coli_015346:0.003770711,C
oli_015340:0.001800838)1.000:0.003006843,(Coli_015325:0.001645346,Coli_015324:0.001347190)0.723:0.000146357)1.000:0.00
2544292,(Coli_015327:0.002098511,(Coli_015301:0.002549540,((Coli_015298:0.026759110,Coli_015304:0.002264138)0.924:0.00
0583770,(Coli_015302:0.001946437,(Coli_015305:0.000299302,Coli_015307:0.000449677)1.000:0.003450489)0.709:0.000000005)
0.571:0.000000005)1.000:0.003900001)0.924:0.000596666)0.746:0.000151563)0.737:0.000148679)0.988:0.000897109)0.986:0.00
1046641)0.751:0.000148944,((Coli_015287:0.003270456,(Coli_015282:0.000714457,Coli_015280:0.000631879)1.000:0.002248196
)0.881:0.000328635,(Coli_015289:0.004954926,(Coli_015269:0.000596638,(Coli_015271:0.000366402,Coli_015272:0.000530833)
0.877:0.000453557)0.982:0.001049091)0.746:0.000149969)0.782:0.000149391,(Coli_015279:0.002097724,(Coli_015284:0.002099
154,(Coli_015275:0.001015440,(Coli_015291:0.000865898,Coli_015290:0.000632221)0.999:0.001801737)0.981:0.000783225)0.96
6:0.000747780)0.946:0.000000005);
"""

def one_primitive_step(tree, original_num_leaves, percent_left, longest_to_average):
    PRUNE = False
    longest_branch_length = 0
    longest_branch = None
    total_branch_length = 0
    branch_count = 0

    for node in tree.traverse():
        current = node.dist
        total_branch_length += current
        branch_count += 1
        if current > longest_branch_length:
            longest_branch_length = current
            longest_branch = node

    average_branch_length = total_branch_length / branch_count
    leaves_on_one_side = longest_branch.get_leaves()

    all_leaves = set(tree.get_leaves())
    num_leaves = len(all_leaves)
    one_side = len(leaves_on_one_side)
    other_side = num_leaves - one_side
    leaf_names = [leaf.name for leaf in leaves_on_one_side]
    leaves_on_other_side = all_leaves - set(leaves_on_one_side)
    other_leaf_names = [leaf.name for leaf in leaves_on_other_side]

    if longest_branch_length > longest_to_average*average_branch_length:
        if (percent_left * original_num_leaves) / 100 > one_side:
            tree.prune(other_leaf_names, preserve_branch_length=True)
            PRUNE = True
            percent_left -= 100 * one_side / original_num_leaves
        elif percent_left*original_num_leaves/100 > other_side:
            tree.prune(leaf_names, preserve_branch_length=True)
            PRUNE = True
            percent_left -= 100 * other_side / original_num_leaves
    return tree, PRUNE, percent_left

def is_valid_tree(file_path):
    try:
        tree_wannabe = ete3.Tree(file_path, format=1)
        return True
    except:
        return False


def calculate_iqr_threshold(values):
    q1, q3 = np.percentile(values, [25, 75])
    iqr = q3 - q1
    upper_fence = q3 + 3 * iqr
    return upper_fence


def prune_by_branch_length(tree, threshold):
    branch_lengths = [node.dist for node in tree.traverse() if not node.is_root()]
    upper_fence = calculate_iqr_threshold(branch_lengths)
    total_tips = len(tree.get_leaves())
    percent_left = 100 - threshold
    for node in tree.traverse():
        if not node.is_root():
            all_leaves = set(tree.get_leaves())
            num_leaves = len(all_leaves)
            leaves_on_one_side = node.get_leaves()
            one_side = len(leaves_on_one_side)
            leaf_names = [leaf.name for leaf in leaves_on_one_side]
            other_side = num_leaves - one_side
            leaves_on_other_side = all_leaves - set(leaves_on_one_side)
            other_side_names = [leaf.name for leaf in leaves_on_other_side]
            min_tips = min(one_side, other_side)
            if node.dist > upper_fence and min_tips < (percent_left * total_tips) / 100 :
                if one_side > other_side:
                    tree.prune(leaves_on_one_side, preserve_branch_length=True)
                    percent_left -= 100 * other_side / total_tips
                else:
                    tree.prune(leaves_on_other_side, preserve_branch_length=True)
                    percent_left -= 100 * one_side / total_tips
    return tree, percent_left


def prune_by_root_to_tip(tree, percent_left, original_num_leaves):
    tree.set_outgroup(tree.get_midpoint_outgroup())
    root_to_tip_distances = {leaf: tree.get_distance(leaf) for leaf in tree.get_leaves()}
    upper_fence = calculate_iqr_threshold(list(root_to_tip_distances.values()))
    outliers = [leaf for leaf, distance in root_to_tip_distances.items() if distance > upper_fence]
    num_outliers = len(outliers)
    if num_outliers > (percent_left * original_num_leaves) / 100:
        return tree
    else:
        removed_counter = 0
        while removed_counter < (percent_left * original_num_leaves) / 100:
            tree.set_outgroup(tree.get_midpoint_outgroup())
            root_to_tip_distances = {leaf: tree.get_distance(leaf) for leaf in tree.get_leaves()}
            upper_fence = calculate_iqr_threshold(list(root_to_tip_distances.values()))

            extreme_outlier = max(root_to_tip_distances, key=root_to_tip_distances.get)
            if root_to_tip_distances[extreme_outlier] > upper_fence:
                extreme_outlier.detach()
                removed_counter += 1
                percent_left -= 100 / original_num_leaves
            else:
                break
    tree.unroot()
    return tree

def validate_cpa_arguments(name_of_the_tree, root_to_node_ratio, min_num_of_roots, M_n, threshold, beta, radius_ratio, safe_tips, show_plot, show_pruned_tips, name_of_output):

    if not isinstance(name_of_the_tree, str):
        raise TypeError("`name_of_the_tree` must be a string.")
    if not is_valid_tree(name_of_the_tree):
        raise TypeError("`name_of_the_tree` must contain a valid .nwk or .newick tree.")

    if not isinstance(root_to_node_ratio, (int, float)):
        raise TypeError("`root_to_node_ratio` must be an integer or a float.")
    if not (0 <= root_to_node_ratio <= 1):
        raise ValueError("`root_to_node_ratio` is outside the expected [0,1] range.")

    if not isinstance(min_num_of_roots, int):
        raise TypeError("`min_num_of_roots` must be an integer.")
    if min_num_of_roots < 1:
        raise ValueError("`min_num_of_roots` must be at least 1.")

    if M_n != 0 and not callable(M_n):
        raise TypeError("`M_n` must be a function or 0.")

    if not isinstance(threshold, (int, float)):
        raise TypeError("`threshold` must be an integer or a float.")
    if not (0 <= threshold <= 100):
        raise ValueError("`threshold` is outside the expected [0,100] range.")

    if not isinstance(beta, (int, float)):
        raise TypeError("`beta` must be an integer or a float.")
    if not (0 <= beta <= 100):
        raise ValueError("`beta` is outside the expected [0,100] range.")

    if not isinstance(radius_ratio, (int, float)):
        raise TypeError("`radius_ratio` must be an integer or a float.")
    if not (0 <= radius_ratio <= 100):
        raise ValueError("`radius_ratio` is outside the expected [0,100] range.")
    
    if not isinstance(safe_tips, list) or not all(isinstance(tip, str) for tip in safe_tips):
        raise TypeError("`safe_tips` must be a list of strings.")

    if not isinstance(show_plot, bool):
        raise TypeError("`show_plot` must be a boolean.")

    if not isinstance(show_pruned_tips, bool):
        raise TypeError("`show_pruned_tips` must be a boolean.")

    if not isinstance(name_of_output, str):
        raise TypeError("`name_of_output` must be a string.")


def validate_psfa_arguments(name_of_the_tree,threshold,longest_to_average,name_of_output):
        
    if not isinstance(name_of_the_tree, str):
        raise TypeError("`name_of_the_tree` must be a string.")
    if not is_valid_tree(name_of_the_tree):
        raise TypeError("`name_of_the_tree` must contain a valid .nwk or .newick tree.")

    if not isinstance(threshold, (int, float)):
        raise TypeError("`threshold` must be an integer or a float.")
    if not (0 <= threshold <= 100):
        raise ValueError("`threshold` is outside the expected [0,100] range.")

    if not isinstance(longest_to_average, (int,float)):
        raise TypeError("`longest_to_average` must be an integer or a float.")
    if longest_to_average < 0:
        raise ValueError("`longest_to_average` must be at least 0.")
    
    if not isinstance(name_of_output, str):
        raise TypeError("`name_of_output` must be a string.")

def validate_iqr_arguments(name_of_the_tree,threshold,name_of_output):

    if not isinstance(name_of_the_tree, str):
        raise TypeError("`name_of_the_tree` must be a string.")
    if not is_valid_tree(name_of_the_tree):
        raise TypeError("`name_of_the_tree` must contain a valid .nwk or .newick tree.")

    if not isinstance(threshold, (int, float)):
        raise TypeError("`threshold` must be an integer or a float.")
    if not (0 <= threshold <= 100):
        raise ValueError("`threshold` is outside the expected [0,100] range.")

    if not isinstance(name_of_output, str):
        raise TypeError("`name_of_output` must be a string.")

# =============
# Imports and global setup
# =============

import ete3
import numpy as np
import matplotlib.pyplot as plt
import random
import argparse
import os
import sys
import tempfile
from dataclasses import dataclass, field, MISSING, fields
from typing import List, get_origin, Optional
import datetime
import logging
import csv
import pdb
import shutil

# =============
# Dataclasses for configuration
# =============

@dataclass
class BaseConfig:
	input_tree: str
	output: str
	threshold: int = 90

@dataclass
class PSFAConfig(BaseConfig):
	longest_to_average: int = 9

@dataclass
class CPAConfig(BaseConfig):
	root_to_node_ratio: float = 0.1
	min_num_of_roots: int = 15
	M_n: int = 0
	beta: int = 20
	radius_ratio: float = 0
	safe_tips: List[str] = field(default_factory=list)
	show_plot: bool = False
	show_pruned_tips: bool = False

@dataclass
class IQRConfig(BaseConfig):
	pass


# =============
# Utility classes and helpers
# =============

class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
	"""Custom formatter to allow argparse to display both defaults and raw text formatting."""
	pass


def init_logging(logfile=None):
	"""
	Initialize logging to console and optionally to a file.

	Args:
		logfile (str, optional): Path to log file. If None, logs only to console.
	"""
	logger = logging.getLogger()
	logger.setLevel(logging.INFO)

	# Remove existing handlers
	if logger.hasHandlers():
		logger.handlers.clear()

	# Console handler
	console_handler = logging.StreamHandler()
	console_handler.setLevel(logging.INFO)
	console_formatter = logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s',
										   datefmt='%Y-%m-%d %H:%M:%S')
	console_handler.setFormatter(console_formatter)
	logger.addHandler(console_handler)

	# File handler
	if logfile:
		os.makedirs(os.path.dirname(os.path.abspath(logfile)), exist_ok=True)
		file_handler = logging.FileHandler(logfile)
		file_handler.setLevel(logging.INFO)
		file_handler.setFormatter(console_formatter)
		logger.addHandler(file_handler)


def msg(message):
	logging.info(message)

def safe_int(value):
	try:
		return int(value)
	except ValueError:
		raise argparse.ArgumentTypeError(f"Invalid int value: {value}")

def safe_float(value):
	try:
		return float(value)
	except ValueError:
		raise argparse.ArgumentTypeError(f"Invalid float value: {value}")

# =============
# File I/O validation helpers
# =============
def check_input_file(file_path):
	"""Check if an input file exists.

	Args:
		file_path (str): Path to file.

	Raises:
		FileNotFoundError: If file does not exist.
	"""
	if not os.path.isfile(file_path):
		raise FileNotFoundError(f"Input tree file does not exist: {file_path}")
	if not is_valid_tree(file_path):
		raise TypeError(f"{file_path} must contain a valid .nwk or .newick tree.")


def ensure_output_path(file_path: str):
	"""Ensure parent directories exist and are writable for the given output file.

	Creates the directory if missing and checks for write permissions.

	Args:
		file_path (str): Full path to output file.

	Raises:
		PermissionError: If directory cannot be created or is not writable.
	"""
	out_dir = os.path.dirname(os.path.abspath(file_path)) or "."
	try:
		os.makedirs(out_dir, exist_ok=True)  # create parents if missing
	except Exception as e:
		raise PermissionError(f"Cannot create output directory {out_dir}: {e}")

	if not os.access(out_dir, os.W_OK):
		raise PermissionError(f"Cannot write to output directory: {out_dir}")

def is_valid_tree(file_path):
	"""
	Check if a file contains a valid Newick tree.

	Args:
		file_path (str): Path to the tree file.

	Returns:
		bool: True if the file is a valid Newick tree, False otherwise.
	"""
	try:
		tree_wannabe = ete3.Tree(file_path, format=1)
		return True
	except:
		return False

# =============
# Parsing commandline arguments
# =============

def add_config_args(parser, config_class, prefix):
	"""
	Dynamically add argparse arguments based on a dataclass definition,
	with safe type conversion and clear error messages.
	"""
	group = parser.add_argument_group(f"{prefix.upper()} options")

	for f in fields(config_class):
		if f.name == "input_tree":
			continue  # skip global input_tree

		arg_name = f"--{prefix}-{f.name.replace('_', '-')}"
		default = None

		if f.default != MISSING:
			default = f.default
		elif f.default_factory != MISSING:
			default = f.default_factory()

		arg_type = f.type

		# Handle typing.List
		if get_origin(arg_type) is list:
			group.add_argument(arg_name, nargs="*", default=default, help=f.metadata.get("help", None))
		# Handle bool
		elif arg_type is bool:
			if default is True:
				group.add_argument(arg_name, dest=f.name, action="store_false", help=f.metadata.get("help", None))
			else:
				group.add_argument(arg_name, dest=f.name, action="store_true", help=f.metadata.get("help", None))
		# Handle int/float safely
		elif arg_type is int:
			group.add_argument(arg_name, type=safe_int, default=default, help=f.metadata.get("help", None))
		elif arg_type is float:
			group.add_argument(arg_name, type=safe_float, default=default, help=f.metadata.get("help", None))
		# Default
		else:
			group.add_argument(arg_name, type=arg_type, default=default, help=f.metadata.get("help", None))

	return parser

def parse_args():
	"""Parse command-line arguments."""
	parser = argparse.ArgumentParser(
		description="Prune trees using one or more methods:\n"
					"      CPA  - Circular Pruning Algorithm\n"
					"      IQR  - Interquartile Range approach\n"
					"      PSFA - Primitive Straight-Forward approach",
					formatter_class=CustomFormatter
	)

	# global arguments
	parser.add_argument(
		"-i", "--input-tree", required=True, 
		help="One or more input tree files (Newick format)"
	)
	parser.add_argument(
		"-o", "--final-output-tree",
		default="pruned_tree.nwk",
		help="Output tree file name. Final output will be: <output_dir>/<output_tree>"
	)
	parser.add_argument(
		"--stats-file",
		default="pruning_summary.csv",
		help="Filename for CSV summary of pruning stats (saved in --output_dir)"
	)
	parser.add_argument(
		"-m", "--methods",
		nargs="+", default=["PSFA","CPA"],
		choices=["PSFA","CPA","IQR"],
		help="List of methods to apply (space-separated, e.g., 'CPA IQR PSFA')"
	)
	parser.add_argument("--log", default="treepruner.log", help="Log file path")
	parser.add_argument("--output_dir", default=".", help="Directory to save all output trees and reports")

	# Add method-specific arguments
	for cfg_class, prefix in [(PSFAConfig, "psfa"), (CPAConfig, "cpa"), (IQRConfig, "iqr")]:
		parser = add_config_args(parser, cfg_class, prefix)

	return parser.parse_args()

# =============
# Config building and validation
# =============

def build_config(config_class, prefix, args):
	"""
	Build a dataclass instance from CLI args.

	Args:
		config_class (dataclass): The dataclass to instantiate.
		prefix (str): Prefix for CLI arguments (e.g., "psfa").
		args (argparse.Namespace): Parsed CLI arguments.

	Returns:
		dataclass instance: Config instance populated from CLI args.
	"""
	kwargs = {}
	for f in fields(config_class):
		if f.name == "input_tree":
			kwargs[f.name] = args.input_tree
		elif f.name == "output":
			# Use user-provided CLI value if present, else default filename
			cli_value = getattr(args, f"{prefix}_output", None)
			if cli_value:
				filename = os.path.basename(cli_value)
			else:
				filename = f"{prefix}_output_tree.nwk"
			out_path = os.path.join(args.output_dir, filename)
			os.makedirs(os.path.dirname(out_path), exist_ok=True)
			kwargs[f.name] = out_path
		else:
			cli_attr = f"{prefix}_{f.name}"
			if hasattr(args, cli_attr):
				kwargs[f.name] = getattr(args, cli_attr)
			else:
				# fallback to default if not present in CLI
				kwargs[f.name] = f.default if f.default != MISSING else None
	return config_class(**kwargs)

def validate_configs(configs, methods):
	"""Validate ranges and constraints for all pruning configs."""
	
	psfa = configs["PSFA"]
	if "PSFA" in methods:
		if not (0 <= psfa.threshold <= 100):
			raise ValueError(f"PSFA threshold {psfa.threshold} must be between 0 and 100")
		if psfa.longest_to_average < 0:
			raise ValueError(f"PSFA longest_to_average {psfa.longest_to_average} must be >= 0")
	
	cpa = configs["CPA"]
	if "CPA" in methods:
		if not (0 <= cpa.threshold <= 100):
			raise ValueError(f"CPA threshold {cpa.threshold} must be between 0 and 100")
		if not (0 <= cpa.beta <= 100):
			raise ValueError(f"CPA beta {cpa.beta} must be between 0 and 100")
		if not (0 <= cpa.radius_ratio <= 100):
			raise ValueError(f"CPA radius_ratio {cpa.radius_ratio} must be between 0 and 100")
		if not (0 <= cpa.root_to_node_ratio <= 1):
			raise ValueError(f"CPA root_to_node_ratio {cpa.root_to_node_ratio} must be between 0 and 1")
		if cpa.min_num_of_roots < 1:
			raise ValueError(f"CPA min_num_of_roots {cpa.min_num_of_roots} must be >= 1")
		if not all(isinstance(tip, str) for tip in cpa.safe_tips):
			raise TypeError("CPA safe_tips must be a list of strings")

	iqr = configs["IQR"]
	if "IQR" in methods:
		if not (0 <= iqr.threshold <= 100):
			raise ValueError(f"IQR threshold {iqr.threshold} must be between 0 and 100")



# =============
# Pruning tools (helper functions)
# =============

def calculate_iqr_threshold(values):
	"""
	Calculate the upper threshold using the interquartile range (IQR).

	Args:
		values (list of float): List of numeric values (e.g., branch lengths or distances).

	Returns:
		float: Upper fence threshold = Q3 + 3 * IQR.
	"""
	q1, q3 = np.percentile(values, [25, 75])
	iqr = q3 - q1
	upper_fence = q3 + 3 * iqr
	return upper_fence


def one_primitive_step(tree, original_num_leaves, percent_left, longest_to_average):
	"""
	Perform one pruning step in the Primitive Straight-Forward Approach (PSFA).

	Args:
		tree (ete3.Tree): The phylogenetic tree to prune.
		original_num_leaves (int): Number of leaves in the original tree.
		percent_left (float): Percentage of leaves that should remain.
		longest_to_average (float): Threshold ratio of longest branch to average branch length for pruning.

	Returns:
		tuple: (pruned tree, PRUNE (bool if pruning occurred), updated percent_left)
	"""
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
		elif percent_left*original_num_leaves/100 > other_side:
			tree.prune(leaf_names, preserve_branch_length=True)
			PRUNE = True
		num_leaves_remaining = len(tree.get_leaves())
		percent_left = 100 * num_leaves_remaining / original_num_leaves # percent_left modified to avoid negative percentages
		
	return tree, PRUNE, percent_left

def prune_by_branch_length(tree, threshold):
	"""
	Prune branches from the tree based on IQR of branch lengths.

	Args:
		tree (ete3.Tree): Tree to prune.
		threshold (float): Percentage of leaves to retain.

	Returns:
		tuple: (pruned tree, percent of leaves remaining)
	"""
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
	"""
	Prune tree leaves based on outliers in root-to-tip distances.

	Args:
		tree (ete3.Tree): The tree to prune.
		percent_left (float): Target percentage of leaves to keep.
		original_num_leaves (int): Number of leaves in the original tree.

	Returns:
		ete3.Tree: Tree after pruning extreme outlier leaves.
	"""
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


# =============
# Pruning methods (core algorithms)
# =============

def _prune_tree_PSFA(config):
	"""Prune a tree using the Primitive Straight-Forward Approach (PSFA)."""
	tree = ete3.Tree(config.input_tree, format=1, quoted_node_names=True)
	leaf_names = tree.get_leaf_names()
	num_leaves = len(leaf_names)
	original_num_leaves = num_leaves
	
	percent_left = 100 - config.threshold

	while percent_left > 0:
		orig_percent_left = percent_left
		tree, PRUNE, percent_left = one_primitive_step(
			tree,
			original_num_leaves,
			percent_left,
			config.longest_to_average
		)
		if not PRUNE or percent_left == orig_percent_left:
			break
	num_leaves = len(tree.get_leaves())
	return tree, 100 * num_leaves / original_num_leaves

def _prune_tree_CPA(config):
	"""Prune a tree using the Circular Pruning Algorithm (CPA)."""
	tree = ete3.Tree(config.input_tree, format=1, quoted_node_names=True)
	original_num_leaves = len(tree.get_leaves())

	if config.M_n == 0:
		M_n_int = original_num_leaves
	elif isinstance(config.M_n, int):
		M_n_int = max(config.M_n, 1)
	else:
		raise TypeError("M_n must be int")

	def hang(tree1, root1):
		return [tree1.get_distance(root1, leaf) for leaf in tree1.iter_leaves()]

	non_leaf_nodes = [node for node in tree.traverse() if not node.is_leaf()]
	num_non_leaf_nodes = len(non_leaf_nodes)

	num_roots = min(max(int(config.root_to_node_ratio * num_non_leaf_nodes), config.min_num_of_roots), num_non_leaf_nodes)
	roots = random.sample(non_leaf_nodes, num_roots)
	
	if not roots:
		raise RuntimeError("No valid roots found in CPA pruning.")

	midpoint_root = tree.get_midpoint_outgroup()
	if not midpoint_root.is_leaf():
		roots.append(midpoint_root)

	best_p_v = float("inf")
	for root in roots:
		data = hang(tree, root)
		original_radius = max(data)
		bin_edges = np.linspace(min(data), original_radius, M_n_int + 1)
		hist, _ = np.histogram(data, bins=bin_edges)

		cumulative_freq = np.cumsum(hist)
		total_freq = cumulative_freq[-1]
		threshold_index = np.argmax(cumulative_freq >= (config.threshold / 100) * total_freq)
		stop = M_n_int - 1

		if threshold_index == M_n_int - 1:
			stop = threshold_index
		elif threshold_index == M_n_int - 2 and hist[threshold_index + 1] < (config.beta / 100) * hist[threshold_index]:
			stop = threshold_index
		else:
			for i in range(threshold_index, M_n_int - 1):
				if hist[i + 1] < (config.beta / 100) * hist[i]:
					stop = i
					break

		p_v = bin_edges[stop + 1]
		if p_v < best_p_v:
			best_p_v = p_v
			best_root = root
			best_stop = stop
			best_original_radius = original_radius
	
	if config.show_plot:
		data = hang(tree, best_root)
		bin_edges = np.linspace(min(data), max(data), M_n_int + 1)
		hist, _ = np.histogram(data, bins=bin_edges)
		plt.hist(data, bins=bin_edges, edgecolor='black', align='mid')
		for i, patch in enumerate(plt.gca().patches):
			if i > best_stop:
				patch.set_facecolor('red')
		plt.xlabel('Distance from root')
		plt.ylabel('Frequency')
		plt.title("CPA")
		plt.show()

	if (100 - 100 * best_p_v / best_original_radius) >= config.radius_ratio:
		leaves_to_keep = []
		leaves_to_prune = []
		for leaf in tree.iter_leaves():
			distance_to_node = best_root.get_distance(leaf)
			if distance_to_node <= best_p_v or leaf.name in config.safe_tips:
				leaves_to_keep.append(leaf)
			else:
				leaves_to_prune.append(leaf.name)

		tree.prune(leaves_to_keep, preserve_branch_length=True)
		num_leaves = len(tree.get_leaves())
		return tree, 100 * num_leaves / original_num_leaves, leaves_to_prune

	return tree, 100, []

def _prune_tree_IQR(config):
	"""Prune a tree using the Interquartile Range (IQR) method."""
	tree = ete3.Tree(config.input_tree, format=1, quoted_node_names=True)
	original_num_leaves = len(tree.get_leaves())
	tree_1, percent_left = prune_by_branch_length(tree, config.threshold)
	tree_2 = prune_by_root_to_tip(tree_1, percent_left, original_num_leaves)
	num_leaves = len(tree_2.get_leaves())
	if num_leaves:
		tree_2.standardize()
	return tree_2, 100 * num_leaves / original_num_leaves


# =============
#Pruning methods - to be called from python script
# =============

def prune_tree_PSFA(
		name_of_the_tree: str,
		threshold: int = 90,
		longest_to_average: int = 9,
		name_of_output: str = "psfa_tree.nwk"
	) -> float:
	"""
	Prune a tree using the PSFA method via Python function call.

	Returns:
		float: Percentage of leaves remaining.
	"""
	config = PSFAConfig(
		input_tree=name_of_the_tree,
		output=name_of_output,
		threshold=threshold,
		longest_to_average=longest_to_average
	)
	tree, percent_left = _prune_tree_PSFA(config)
	return percent_left


def prune_tree_CPA(
		name_of_the_tree: str,
		root_to_node_ratio: float = 0.1,
		min_num_of_roots: int = 15,
		M_n: int = 0,
		threshold: int = 90,
		beta: float = 20,
		radius_ratio: float = 0,
		safe_tips: Optional[List[str]] = None,
		show_plot: bool = False,
		show_pruned_tips: bool = False,
		name_of_output: str = "cpa_tree.nwk"
	) -> float:
	"""
	Prune a tree using the CPA method via Python function call.

	Returns:
		float: Percentage of leaves remaining.
	"""
	if safe_tips is None:
		safe_tips = []

	config = CPAConfig(
		input_tree=name_of_the_tree,
		output=name_of_output,
		root_to_node_ratio=root_to_node_ratio,
		min_num_of_roots=min_num_of_roots,
		M_n=M_n,
		threshold=threshold,
		beta=beta,
		radius_ratio=radius_ratio,
		safe_tips=safe_tips,
		show_plot=show_plot,
		show_pruned_tips=show_pruned_tips
	)
	tree, percent_left, pruned_tips = _prune_tree_CPA(config)
	return percent_left, pruned_tips


def prune_tree_IQR(
		name_of_the_tree: str,
		threshold: int = 90,
		name_of_output: str = "iqr_tree.nwk"
	) -> float:
	"""
	Prune a tree using the IQR method via Python function call.

	Returns:
		float: Percentage of leaves remaining.
	"""
	config = IQRConfig(
		input_tree=name_of_the_tree,
		output=name_of_output,
		threshold=threshold
	)
	tree, percent_left = _prune_tree_IQR(config)
	return percent_left

# =============
# Executing methods
# =============

def run_method_with_stats(method_name, config):
	"""
	Run a pruning method and return stats.
	Returns: pruned_tree, percent_remaining, pruned_tips (list)
	"""
	tree = ete3.Tree(config.input_tree, format=1, quoted_node_names=True)
	pruned_tips = []
	if method_name == "PSFA":
		outtree, percent_remaining = _prune_tree_PSFA(config)
	elif method_name == "CPA":
		outtree, percent_remaining, pruned_tips = _prune_tree_CPA(config)
	elif method_name == "IQR":
		outtree, percent_remaining = _prune_tree_IQR(config)
	else:
		raise RuntimeError(f"Unknown method: {method_name}")
	outtree.write(outfile=config.output, format=1)
	return outtree, percent_remaining, pruned_tips



def run_pruning_methods(tree_file, out_tree_file, methods, configs):
	"""
	Run each method sequentially and collect stats.
	Returns a list of dicts with method stats.
	"""
	stats = []
	prev_file = tree_file
	original_tree = ete3.Tree(tree_file, format=1)
	original_leaves = len(original_tree.get_leaves())

	for method in methods:
		cfg = configs[method]
		logging.info(f"Running {method}")
		cfg.input_tree = prev_file
		
		# Run method
		tree_out, percent_left, pruned_tips = run_method_with_stats(method, cfg)
		
		stats.append({
		 "method": method,
		 "output": cfg.output,
		 "leaves_remaining": len(tree_out.get_leaves()),
		 "percent_remaining": percent_left,
		 "pruned_tips": pruned_tips
		})

		prev_file = cfg.output
	
	shutil.copy(cfg.output, out_tree_file)
	
	return original_leaves, stats


def save_stats_csv(stats, original_leaves, output_path):
	"""
	Save pruning statistics to CSV
	"""
	os.makedirs(os.path.dirname(output_path), exist_ok=True)
	with open(output_path, "w", newline="") as f:
		writer = csv.writer(f)
		writer.writerow(["Method", "Original Leaves", "Leaves Remaining", "Percent Remaining", "Pruned Tips"])
		for s in stats:
			tips_str = ";".join(s["pruned_tips"]) if s["pruned_tips"] else ""
			writer.writerow([
				s["method"],
				original_leaves,
				s["leaves_remaining"],
				f"{s['percent_remaining']:.1f}",
				tips_str
			])
	logging.info(f"Pruning summary saved to {output_path}")


# =============
# Main driver
# =============

def main():
	"""Main entry point for CLI."""
	args = parse_args()
	args.final_output_tree = os.path.join(args.output_dir, args.final_output_tree)
	stats_path = os.path.join(args.output_dir, args.stats_file)
	init_logging(logfile=args.log)
	logging.info(f"Resolved output tree path: {args.final_output_tree}")
	methods = args.methods
	configs = {
		"PSFA": build_config(PSFAConfig, "psfa", args),
		"CPA": build_config(CPAConfig, "cpa", args),
		"IQR": build_config(IQRConfig, "iqr", args),
	}
	
	
	try:
		msg("Validating input file and config")
		check_input_file(args.input_tree)
		ensure_output_path(args.final_output_tree)
		validate_configs(configs, methods)
		msg("Start pruning pipeline")
		original_leaves, stats = run_pruning_methods(args.input_tree, args.final_output_tree, methods, configs)
		save_stats_csv(stats, original_leaves, stats_path)
		msg(f"Pruning completed successfully. Output saved to {args.final_output_tree}")
	except Exception as e:
		# Print only error type and message in red
		tb = e.__traceback__
		while tb.tb_next:
			tb = tb.tb_next
		line_no = tb.tb_lineno
		# Print only error type, message, and line number in red
		logging.error(f"\033[91m[{type(e).__name__}]: line {line_no}: {e}\033[0m", file=sys.stderr)
		sys.exit(1)



if __name__ == "__main__":
	main()

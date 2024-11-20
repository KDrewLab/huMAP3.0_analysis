# Scripts used to identify mutually exclusive and structurally consistent protein pairs from paired predicted dimeric (AlphaFold2) structures 

Structures are sourced from Burke, D. F. et al. Towards a structurally resolved human protein interaction network. Nat. Struct. Mol. Biol. 1–10 (2023) doi:10.1038/s41594-022-00910-8.

<details><summary>Hello</summary>World!</details>

# Installation

<details>
	<summary>Installation information</summary>
	1. **PyMol** is needed to perform the alignments of the paired structures on the common subunit. We recommend using a conda environment with the open source version of PyMol. Version 3.0.0 was used at the time of script development. [Anaconda package for open-source PyMol](https://anaconda.org/conda-forge/pymol-open-source)
	
	2. If the plan is to run thousands of structures, it is recommended to parallelize the process. A straightforward way would be use [GNU Parallel] (https://www.gnu.org/software/parallel/parallel_tutorial.html)
	
	3. After, clone this repo and the user is ready to run the scripts and evaluate whether overlapping interfaces exist between unique protein subunits aligned to a common subunit. 

</details>

# Example commandlines

<details>
	<summary>Examples for running the scripts</summary>
	
	1. To run on a single example of paired structures,you pass the structure files and directly identify which are the unique and common chains for each structure to *evaluate_structure_overlap.py* :
	
	'''bash
	python3 evaluate_structure_overlap.py --pdb1 P62917-P62841.pdb --pdb2 P62917-P62847.pdb --common_ch1 "chain A" --common_ch2 "chain A" --test_ch1 "chain B" --test_ch2 "chain B"	
	'''
	
	2. *process_dimer_pair_lines.py* will take an input of structure file names with two structures per line and feed them to the *evaluate_structure_overlap.py* to perform the alignment and evaluation for many structures. Example of the format can be found in the *example_pairs.txt* 
	
	3. For the structure of this script, it expects protein structures are annotated as protein pairs, where each protein/chain is identified in the PDB filename. For example: structure1 contains protein A and protein B with its name containing these two IDs separated by a delimiter.  
 
	'''text
	proteinA-proteinB
	'''  
 
	This delimiter is passed as an argument (e.g., *--pair_delim*)
	
	4. There should be a delimiter between the two structures (e.g., structure1 contains protein A and protein B, while structure2 contains protein A and protein C)
	'''text
	proteinA-proteinB proteinA-proteinC
	'''
	This delimiter is passed as an argument (e.g., *--dimer_delim*)

	5. Example running from the command-line using a text file input:
	'''bash
	while IFS= read -r line; do python3 ./process_dimer_pair_lines.py --input_line "$line" --pair_delim '-' --dimer_delim ' '; done < example_pairs.txt
	'''

	6. Optionally, you can parallelize this process. Here is a command-line example using GNU parallel:
	'''bash
	cat example_pairs.txt | parallel -j 10 python3 ./process_dimer_pair_lines.py --input_line '{}' --pair_delim '-' --dimer_delim '\ '
	'''

</details>







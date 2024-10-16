import subprocess
import argparse
import os
import sys
from datetime import datetime

# define parser
parser = argparse.ArgumentParser(description="Compare the overlap between heterodimeric pdb structures to test for mutual exclusion. Takes in a .txt file to run many at once.")

# define args
parser.add_argument("--input_line", required=True, help="Line of a file which contains pairs of heterodimers.")
parser.add_argument("--pair_delim", required=True, help="Delimiter between individual proteins within a pair.")
parser.add_argument("--dimer_delim", required=True, help="Delimeter between dimers for each line in the file.")

# parse args
args = parser.parse_args()

# ensure the error log directory exists
error_log_dir = 'error_logs'
os.makedirs(error_log_dir, exist_ok=True)

# generate a unique filename for the error log
timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
error_log_file = os.path.join(error_log_dir, f'error_log_{timestamp}.txt')

def log_error(error_message):
    with open(error_log_file, 'a') as f:
        f.write(error_message)
    print(error_message, file=sys.stderr)

# function to process lines in the textfile to retrieve pdb structures for modeling
def process_dimer_pairs(input_line, pair_delim, dimer_delim):
    """
     This function will process lines read from a .txt file that contains a pair of heterodimers, with one common protein between them. It will identify the common protein and unique proteins between the dimer pairs and then use this info to grab the pdb files and seed the parameters necessary to run the overlap script.
    """ 
    # file processing can be edited for a different use case
    # in this situation, it is assumed that each .pdb file for a dimer is stored in a 
    # directory named after the pair, in order of their chains
    # e.g., /protein1-protein2 protein1-protein2.pdb
    
    try:
        # read file lines
        # get paths for pdb structures
        dimers = input_line.strip().split(dimer_delim)
        pdb1 = f"{dimers[0]}/{dimers[0]}.pdb"
        pdb2 = f"{dimers[1]}/{dimers[1]}.pdb"
           
        # determine common and unique proteins
        prots_pr1 = dimers[0].split(pair_delim)
        print(f"proteins in pair 1: {prots_pr1}")
        set1 = set(prots_pr1)
        prots_pr2 = dimers[1].split(pair_delim)
        print(f"proteins in pair 2: {prots_pr2}")
        set2 = set(prots_pr2)
        common_prot = set1.intersection(set2)
        print(f"The common protein is: {common_prot}")
        unique_prot1 = (set1 - common_prot).pop()
        unique_prot2 = (set2 - common_prot).pop()
        common_prot = common_prot.pop()
            
        # extract chain identifiers
        commn_ch1 = 'chain A' if common_prot == prots_pr1[0] else 'chain B'
        commn_ch2 = 'chain A' if common_prot == prots_pr2[0] else 'chain B'
        u_ch1 = 'chain A' if unique_prot1 == prots_pr1[0] else 'chain B'
        u_ch2 = 'chain A' if unique_prot2 == prots_pr2[0] else 'chain B'
            
        # to get these values iteratively we use yield and not return
        return pdb1, pdb2, commn_ch1, commn_ch2, u_ch1, u_ch2, dimers
    
    except Exception as e:
        error_message = f"Error processing input line '{input_line}': {e}\n"
        log_error(error_message)
        raise

# Implementation
if __name__ == "__main__":
    try:
        pdb1, pdb2, cmn_ch1, cmn_ch2, ch1, ch2, dimers = process_dimer_pairs(args.input_line, args.pair_delim, args.dimer_delim)
        print(f"the pairs are: {dimers}")
        print(f"pdb1 is {pdb1}, with common chain: {cmn_ch1} and unique chain: {ch1}")
        print(f"pdb2 is {pdb2}, with common chain: {cmn_ch2} and unique chain: {ch2}")
        print(f"\n\n")
       
        output_file = f"overlap_results/{dimers[0]}_{dimers[1]}_overlap.txt"

    
        print(f"Running evaluate_structure_overlap.py with the following arguments:")
        print(f"--pdb1 {pdb1} --pdb2 {pdb2} --common_ch1 {cmn_ch1} --common_ch2 {cmn_ch2} --test_ch1 {ch1} --test_ch2 {ch2} --output_file {output_file}")

        # Seed params to evaluate_structure_overlap.py
        subprocess.run([
            "python", "/stor/project/sfisch6/scripts/pymol_scripts/evaluate_structure_overlap.py",
            "--pdb1", pdb1,
            "--pdb2", pdb2,
            "--common_ch1", cmn_ch1,
            "--common_ch2", cmn_ch2,
            "--test_ch1", ch1,
            "--test_ch2", ch2,
            "--output_file", output_file
        ], check=True)
    
    except ValueError as e:
        error_message = f"ValueError processing line: {args.input_line}. Error: {e}\n"
        log_error(error_message)
    except subprocess.CalledProcessError as e:
        error_message = f"Subprocess error processing PDBs {pdb1} and {pdb2} with chains {cmn_ch1}, {cmn_ch2}, {ch1}, {ch2}. Error: {e}\n"
        log_error(error_message)
    except Exception as e:
        error_message = f"General error processing PDBs {pdb1} and {pdb2} with chains {cmn_ch1}, {cmn_ch2}, {ch1}, {ch2}. Error: {e}\n"
        log_error(error_message)

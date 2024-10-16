import os
from pymol import cmd
import argparse
import sys
from datetime import datetime

# set-up parser
parser = argparse.ArgumentParser(description="Compare overlap between two PDB structures.")

# define arguments
parser.add_argument("--pdb1", required=True, help="Path to the first PDB file.")
parser.add_argument("--pdb2", required=True, help="Path to the second PDB file.")
parser.add_argument("--common_ch1", required=True, help="Common chain identifier in the first PDB file. Chain labels are case-sensitive.")
parser.add_argument("--common_ch2", required=True, help="Common chain identifier in the second PDB file. Chain labels are case-sensitive.")
parser.add_argument("--test_ch1", required=True, help="Test chain identifier in the first PDB file. Chain labels are case-sensitive.")
parser.add_argument("--test_ch2", required=True, help="Test chain identifier in the second PDB file. Chain labels are case-sensitive.")
parser.add_argument("--output_file", required=False, help="Path to the output file. If not specified, results will be printed to the console.")

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

# function to evaluate overlap in two pdb structures
def compare_overlap(pdb1, pdb2, common_ch1, common_ch2, test_ch1, test_ch2):
    """
    This function takes in two pdb structures to compare. The goal is to determine the overlap between two unique proteins which are aligned to a common protein. To perform the comparison, the chain for the common protein to be aligned must be defined for each structure (e.g., PSMA2 is 'chain A' in PDB1 and 'chain B' in PDB2). The proteins which will be analyzed for overlap must have their chains defined in each structure as well, which are the test chains. Chains should be defined in the same manner (e.g., 'chain C', 'chain D'). Reminder: chain identifiers are case-sensitive.
    """

    try:
        # load in pdb structures
        cmd.load(pdb1, 'structure1')
        cmd.load(pdb2, 'structure2')

        # align the common protein in both structures
        cmd.align(f'structure1 and {common_ch1}', f'structure2 and {common_ch2}')

        # select the interfaces between common proteins and defined unique proteins in each structure
        cmd.select('struct1_interface', f'byres (structure1 and {common_ch1}) within 4.0 of (structure1 and {test_ch1})')
        cmd.select('struct2_interface', f'byres (structure2 and {common_ch2}) within 4.0 of (structure2 and {test_ch2})')

        # evaluate overlap of interfaces if it exists
        cmd.select('interface_overlap', f'byres struct1_interface within 4.0 of struct2_interface')
        interface_overlap = cmd.get_model('interface_overlap')

        # evaluate overlap of chains if it exists
        cmd.select('chain_overlap', f'byres (structure1 and {test_ch1}) within 4.0 of (structure2 and {test_ch2})')
        chain_overlap = cmd.get_model('chain_overlap')

        # collate overlap information
        interface_overlap_atoms = set([atom.resi for atom in interface_overlap.atom])
        chain_overlap_atoms = set([atom.resi for atom in chain_overlap.atom])

        return {
            'interface_overlap': 'yes' if interface_overlap_atoms else 'no',
            'overlapping_intf_residues': interface_overlap_atoms,
            'num_of_overlapping_interface_res': len(interface_overlap_atoms),
            'chain_overlap': 'yes' if chain_overlap_atoms else 'no',
            'overlapping_chain_res': chain_overlap_atoms,
            'num_overlapping_chain_res': len(chain_overlap_atoms)
            }   
    
    except pymol.CmdException as e:
        error_message = f"PyMOL command exception processing PDBs {pdb1} and {pdb2} with chains {common_ch1}, {common_ch2}, {test_ch1}, {test_ch2}: {e}\n"
        log_error(error_message)
        raise
    except Exception as e:
        error_message = f"General exception processing PDBs {pdb1} and {pdb2} with chains {common_ch1}, {common_ch2}, {test_ch1}, {test_ch2}: {e}\n"
        log_error(error_message)
        raise

# implementation of overlap evaluation
if __name__ == "__main__":
    
    try:
        result = compare_overlap(args.pdb1, args.pdb2, args.common_ch1, args.common_ch2, args.test_ch1, args.test_ch2)

        if args.output_file:
            with open(args.output_file, 'w') as f:
                f.write(f"Interface Overlap: {result['interface_overlap']}\n")
                f.write(f"Overlapping Interface Residues: {result['overlapping_intf_residues']}\n")
                f.write(f"Number of Overlapping Interface Residues: {result['num_of_overlapping_interface_res']}\n")
                f.write(f"Chain Overlap: {result['chain_overlap']}\n")
                f.write(f"Overlapping Chain Residues: {result['overlapping_chain_res']}\n")
                f.write(f"Number of Overlapping Chain Residues: {result['num_overlapping_chain_res']}\n")
        else:
            print(f"Interface Overlap: {result['interface_overlap']}")
            print(f"Overlapping Interface Residues: {result['overlapping_intf_residues']}")
            print(f"Number of Overlapping Interface Residues: {result['num_of_overlapping_interface_res']}")
            print(f"Chain Overlap: {result['chain_overlap']}")
            print(f"Overlapping Chain Residues: {result['overlapping_chain_res']}")
            print(f"Number of Overlapping Chain Residues: {result['num_overlapping_chain_res']}")
    
    except Exception as e:
        error_message = f"Error during overlap evaluation for PDBs {args.pdb1} and {args.pdb2} with chains {args.common_ch1}, {args.common_ch2}, {args.test_ch1}, {args.test_ch2}: {e}\n"
        log_error(error_message)

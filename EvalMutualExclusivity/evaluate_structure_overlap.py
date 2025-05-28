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
        alignment = cmd.align(f'structure1 and {common_ch1}', f'structure2 and {common_ch2}')
        rmsd = alignment[0]

        # select the interfaces between common proteins and defined unique proteins in each structure
        cmd.select('struct1_interface', f'byres (structure1 and {common_ch1}) within 4.0 of (structure1 and {test_ch1})')
        cmd.select('struct2_interface', f'byres (structure2 and {common_ch2}) within 4.0 of (structure2 and {test_ch2})')

        # evaluate overlap of interfaces if it exists
        cmd.select('interface_overlap', f'byres struct1_interface within 4.0 of struct2_interface')
        interface_overlap = cmd.get_model('interface_overlap')

        # evaluate overlap of both test chains if it exists
        cmd.select('chain_overlap_1', f'byres (structure1 and {test_ch1}) within 4.0 of (structure2 and {test_ch2})')
        cmd.select('chain_overlap_2', f'byres (structure2 and {test_ch2}) within 4.0 of (structure1 and {test_ch1})')
        chain_overlap_1 = cmd.get_model('chain_overlap_1')
        chain_overlap_2 = cmd.get_model('chain_overlap_2')

        # collate overlap information
        interface_overlap_atoms = {atom.resi for atom in interface_overlap.atom}
        #chain_overlap_atoms = set([atom.resi for atom in chain_overlap.atom])

        # define parsing functions for pymol to process
        def store_chain_overlap(resi, resn, b):
            chain1_overlap_tpl.append((f"{resi}{resn}", b))
        def store_chain_overlap2(resi, resn, b):
            chain2_overlap_tpl.append((f"{resi}{resn}", b))

        # extract info per chain
        chain1_overlap_tpl = []
        cmd.iterate("chain_overlap_1 and name CA", "store_chain_overlap(resi, resn, b)", space={'store_chain_overlap': store_chain_overlap})
        chain1_overlap_res = [tup[0] for tup in chain1_overlap_tpl]
        chain1_overlap_pLDDT = [tup[1] for tup in chain1_overlap_tpl]
        chain2_overlap_tpl = []
        cmd.iterate("chain_overlap_2 and name CA", "store_chain_overlap2(resi, resn, b)", space={'store_chain_overlap2': store_chain_overlap2})
        chain2_overlap_res = [tup[0] for tup in chain2_overlap_tpl]
        chain2_overlap_pLDDT = [tup[1] for tup in chain2_overlap_tpl]

        # calculate average pLDDT per chain overlap
        av_pLDDT_chain_ov_1 = sum(chain1_overlap_pLDDT)/len(chain1_overlap_res) if chain1_overlap_res else None
        av_pLDDT_chain_ov_2 = sum(chain2_overlap_pLDDT)/len(chain2_overlap_res) if chain2_overlap_res else None

       # print(f"\n\ntest: {chain_overlap_test}")
        return {
            'interface_overlap': 'yes' if interface_overlap_atoms else 'no',
            'overlapping_intf_residues': interface_overlap_atoms,
            'num_of_overlapping_interface_res': len(interface_overlap_atoms),
            'alignment_rmsd': rmsd,
            'chain_overlap': 'yes' if (chain1_overlap_res or chain2_overlap_res) else 'no',
            'overlapping_chain1_res': chain1_overlap_res,
            'overlapping_chain2_res': chain2_overlap_res,
            'num_overlapping_chain1_res': len(chain1_overlap_res),
            'num_overlapping_chain2_res': len(chain2_overlap_res),
            'chain1_overlap_av_pLDDT': av_pLDDT_chain_ov_1,
            'chain2_overlap_av_pLDDT': av_pLDDT_chain_ov_2
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
                if isinstance(result['alignment_rmsd'], (int, float)):
                    f.write(f"Alignment RMSD: {result['alignment_rmsd']:.2f}\n")
                else:
                    f.write("Alignment RMSD: N/A\n")
                f.write(f"Chain Overlap: {result['chain_overlap']}\n")
                f.write(f"Overlapping Chain 1 Residues: {result['overlapping_chain1_res']}\n")
                f.write(f"Number of Overlapping Chain 1 Residues: {result['num_overlapping_chain1_res']}\n")
                if isinstance(result['chain1_overlap_av_pLDDT'], (int, float)):
                    f.write(f"Average pLDDT for Overlapping Chain 1: {result['chain1_overlap_av_pLDDT']:.2f}\n")
                else:
                    f.write("Average pLDDT for Overlapping Chain 1: N/A\n")
                f.write(f"Overlapping Chain 2 Residues: {result['overlapping_chain2_res']}\n")
                f.write(f"Number of Overlapping Chain 2 Residues: {result['num_overlapping_chain2_res']}\n")
                if isinstance(result['chain2_overlap_av_pLDDT'], (int, float)):
                    f.write(f"Average pLDDT for Overlapping Chain 2: {result['chain2_overlap_av_pLDDT']:.2f}\n")
                else:
                    f.write("Average pLDDT for Overlapping Chain 2: N/A\n")
        else:
            print(f"Interface Overlap: {result['interface_overlap']}")
            print(f"Overlapping Interface Residues: {result['overlapping_intf_residues']}")
            print(f"Number of Overlapping Interface Residues: {result['num_of_overlapping_interface_res']}")
            print(f"Alignment RMSD: {result['alignment_rmsd']:.2f}" if result['alignment_rmsd'] else "N/A")
            print(f"Chain Overlap: {result['chain_overlap']}")
            print(f"Overlapping Chain 1 Residues: {result['overlapping_chain1_res']}")
            print(f"Number of Overlapping Chain 1 Residues: {result['num_overlapping_chain1_res']}")
            print(f"Average pLDDT for Overlapping Chain 1: {result['chain1_overlap_av_pLDDT']:.2f}" if result['chain1_overlap_av_pLDDT'] else "N/A")
            print(f"Overlapping Chain 2 Residues: {result['overlapping_chain2_res']}")
            print(f"Number of Overlapping Chain 2 Residues: {result['num_overlapping_chain2_res']}")
            print(f"Average pLDDT for Overlapping Chain 2: {result['chain2_overlap_av_pLDDT']:.2f}" if result['chain2_overlap_av_pLDDT'] else "N/A")

    
    except Exception as e:
        error_message = f"Error during overlap evaluation for PDBs {args.pdb1} and {args.pdb2} with chains {args.common_ch1}, {args.common_ch2}, {args.test_ch1}, {args.test_ch2}: {e}\n"
        log_error(error_message)

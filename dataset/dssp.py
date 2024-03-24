from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import subprocess

# Parse the PDB file
parser = PDBParser()
structure = parser.get_structure('182L', '/Users/allalelhommad/PYT/SBPYT_project/dataset/sample_pocket/182l_pocket.pdb')
model = structure[0]

# Run mkdssp as subprocess to generate DSSP output file
with open('dssp_output.txt', 'w') as output_file:
    dssp_process = subprocess.Popen(['/opt/homebrew/bin/mkdssp', '--calculate-accessibility', '/Users/allalelhommad/PYT/SBPYT_project/dataset/sample_pocket/182l_pocket.pdb'], stdout=output_file, stderr=subprocess.PIPE)
    _, stderr = dssp_process.communicate()

    # Check for errors
    if dssp_process.returncode != 0:
        print("Error running mkdssp:", stderr.decode())
        exit()

# Read DSSP output from file and parse it into a dictionary
with open('dssp_output.txt', 'r') as dssp_file:
    dssp_dict = MMCIF2Dict(dssp_file)

#print(dssp_dict['_dssp_struct_summary.secondary_structure'])
#print(dssp_dict['_dssp_struct_bridge_pairs.acceptor_1_label_seq_id'])
print(dssp_dict['_dssp_struct_bridge_pairs.donor_2_energy'])
#print(dssp_dict['_dssp_struct_summary.helix_3_10'])
#print(dssp_dict['_dssp_struct_summary.helix_alpha'])
#print(dssp_dict['_dssp_struct_summary.helix_pi'])
#print(dssp_dict['_dssp_struct_summary.helix_pp'])
#print(dssp_dict['_dssp_struct_summary.bend'])
#print(dssp_dict['_dssp_struct_summary.chirality'])
#print(dssp_dict['_dssp_struct_summary.TCO'])
#print(dssp_dict['_dssp_struct_summary.kappa'])
#print(dssp_dict['_dssp_struct_summary.alpha'])
#print(dssp_dict['_dssp_struct_summary.phi'])
#print(dssp_dict['_dssp_struct_summary.psi'])
import os
os.remove('dssp_output.txt')

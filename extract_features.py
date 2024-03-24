from Bio.PDB import PDBParser 
import torch
import pandas as pd 
import torch.nn as nn
import subprocess
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import os
import numpy as np



'''

    ProteinFeatures class to extract features from a PDB file, such as amino acid composition, secondary structure, etc. 
    The class will be used to extract features from a list of PDB files and add them to a dataframe for further analysis.
    The class will also be used to encode the amino acid sequence into a word embedding for use in a neural network.
    
'''


class ProteinFeatures:
    
    
    def __init__(self, pdb_file):
        self.pdb_file = pdb_file
        self.amino_acids = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR', 'SEC', 'PYL']
        self.one_hot_encoding = self.create_one_hot_encoding()
        torch.manual_seed(42)  # Set a fixed seed for random number generation
        self.aminoacid_frequencies, self.encoded_sequence = self.extract_sequence()
        self.secondary_structure_codes = ['H', 'B', 'E', 'G', 'I', 'T', 'S', '.']
        self.secondary_structure_encoding = self.create_secondary_structure_encoding()
        self.secondary_structure, self.aminoacid_accesibiliity, self.psi_angles, self.phi_angles, self.total_contact, self.encoded_chirality = self.extract_secondary_structure()


    def create_one_hot_encoding(self):
        one_hot_encoding = {}
        for i, amino_acid in enumerate(self.amino_acids):
            one_hot_vector = torch.zeros(len(self.amino_acids))
            one_hot_vector[i] = 1
            one_hot_encoding[amino_acid] = one_hot_vector
        return one_hot_encoding


    def extract_sequence(self):
        parser = PDBParser()
        structure = parser.get_structure('protein', self.pdb_file)
        model = structure[0]  # Assuming single model
        chain = model['A']  # Assuming chain A, change as needed
        residues = chain.get_residues()

        encoded_sequence = []
        aminoacid_counts = {}
        total_residues = 0

        for residue in residues:
            if residue.get_id()[0] == ' ':
                aminoacid = residue.get_resname()
                encoded_aminoacid = self.one_hot_encoding.get(aminoacid, torch.zeros(len(self.amino_acids)))
                encoded_sequence.append(encoded_aminoacid)
                aminoacid_counts[aminoacid] = aminoacid_counts.get(aminoacid, 0) + 1
                total_residues += 1

        aminoacid_frequencies = {aminoacid: count / total_residues for aminoacid, count in aminoacid_counts.items()}
        aminoacid_frequencies = {self.one_hot_encoding[amino_acid]: frequency for amino_acid, frequency in aminoacid_frequencies.items()}

        
        # Stack the one-hot encoded amino acids
        encoded_sequence = torch.stack(encoded_sequence)

        return aminoacid_frequencies, encoded_sequence
    
    def create_secondary_structure_encoding(self):
        encoding = {}
        for i, code in enumerate(self.secondary_structure_codes):
            one_hot_vector = torch.zeros(len(self.secondary_structure_codes))
            one_hot_vector[i] = 1
            encoding[code] = one_hot_vector
        return encoding
    
    def extract_secondary_structure(self):
        # Run mkdssp as a subprocess
        parser = PDBParser()
        structure = parser.get_structure('182L', 'dataset/182l.pdb')
        model = structure[0]

        # Run mkdssp as subprocess to generate DSSP output file
        with open('dssp_output.txt', 'w') as output_file:
            dssp_process = subprocess.Popen(['/opt/homebrew/bin/mkdssp', '--calculate-accessibility', 'dataset/182l.pdb'], stdout=output_file, stderr=subprocess.PIPE)
            _, stderr = dssp_process.communicate()

            # Check for errors
            if dssp_process.returncode != 0:
                print("Error running mkdssp:", stderr.decode())
                exit()

        # Read DSSP output from file and parse it into a dictionary
        with open('dssp_output.txt', 'r') as dssp_file:
            dssp_dict = MMCIF2Dict(dssp_file)

        secondary_structure = dssp_dict['_dssp_struct_summary.secondary_structure']
        secondary_structure = [self.secondary_structure_encoding[code] for code in secondary_structure]

        aminoacid_accesibiliity = dssp_dict['_dssp_struct_summary.accessibility']
        aminoacid_accesibiliity = [float(accessibility) if accessibility != '.' else np.nan for accessibility in aminoacid_accesibiliity]
        phi_angles = dssp_dict['_dssp_struct_summary.phi']
        psi_angles = dssp_dict['_dssp_struct_summary.psi']
        total_contact = dssp_dict['_dssp_struct_summary.TCO']
        chirality_data = dssp_dict['_dssp_struct_summary.chirality']
        phi_angles = [float(angle) if angle != '.' else np.nan for angle in phi_angles]
        psi_angles = [float(angle) if angle != '.' else np.nan for angle in psi_angles]
        total_contact = [float(angle) if angle != '.' else np.nan for angle in total_contact]
        mapping = {'.': [1, 0, 0], '-': [0, 1, 0], '+': [0, 0, 1]}
        encoded_chirality = np.array([mapping[chirality] for chirality in chirality_data])

        
        
        # Remove the file
        os.remove('dssp_output.txt')

        return secondary_structure, aminoacid_accesibiliity, psi_angles, phi_angles, total_contact,encoded_chirality

        
        

    




# Usage example
pdb_file = 'dataset/sample_pocket/182l_pocket.pdb'
protein = ProteinFeatures(pdb_file)


#print(protein.aminoacid_frequencies)


data = []
for i in range(len(protein.aminoacid_frequencies)):
    aminoacid = list(protein.aminoacid_frequencies.keys())[i]
    frequency = list(protein.aminoacid_frequencies.values())[i]
    data.append({
        'Amino_Acid': aminoacid,
        'Solvent_Accessibility': protein.aminoacid_accesibiliity[i],
        'Secondary_Structure': protein.secondary_structure[i],
        'Relative_Frequency': frequency,
        'Phi_Angle': protein.phi_angles[i],
        'Psi_Angle': protein.psi_angles[i],
        'TCO' : protein.total_contact[i],
        'Chirality': protein.encoded_chirality[i],
        'Is_Binding_Pocket': 1,


        
    })

# Create DataFrame
df = pd.DataFrame(data)

# Display DataFrame
print(df['Solvent_Accessibility'][1:10])


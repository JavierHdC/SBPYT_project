from Bio.PDB import PDBParser, NeighborSearch, Selection, is_aa 
import pandas as pd 
import subprocess
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import os
import numpy as np
import torch
import torch.nn as nn
import os
import numpy as np
from Bio.PDB.DSSP import DSSP
from glob import glob



# Set the working directory
os.chdir('/Users/javierherranzdelcerro/Desktop/PYT_SBI/SBPYT_project')

class ProteinFeatures:
    
    def __init__(self, pdb_file, pocket_pdb_file=None):
        self.dssp_executable = '/opt/homebrew/bin/mkdssp' 
        self.pocket_pdb_file = pocket_pdb_file
        self.pdb_file = pdb_file
        self.parser = PDBParser(QUIET=True)
        self.amino_acids = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 
                            'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
        self.secondary_structure_codes = ['H', 'B', 'E', 'G', 'I', 'T', 'S', ' ', '-', 'P']
        self.pocket_residues = self.load_pocket_residues()
        torch.manual_seed(42)  # Set a fixed seed for random number generation
        self.structure = self.parser.get_structure('protein', self.pdb_file)
        self.encoding_ss, self.encoding_aa = self.one_hot_encoding() 

    
    
    def extract_sequence(self):
        parser = PDBParser()
        structure = parser.get_structure('protein', self.pdb_file)
        model = structure[0]  # Assuming single model
        chain = model['A']  # Assuming chain A, change as needed
        residues = chain.get_residues()

        aminoacid_counts = {}
        total_residues = 0

        for residue in residues:
            if residue.get_id()[0] == ' ':
                aminoacid = residue.get_resname()
                aminoacid_counts[aminoacid] = aminoacid_counts.get(aminoacid, 0) + 1
                total_residues += 1

        aminoacid_frequencies = {aminoacid: count / total_residues for aminoacid, count in aminoacid_counts.items()}

        return aminoacid_frequencies
    
    def calculate_total_contact(self):
        """Calculate the total contact for each residue based on a distance threshold."""
        atom_list = Selection.unfold_entities(self.structure, 'A')  # Use self.structure
        ns = NeighborSearch(atom_list)
        total_contact_dict = {}

        for chain in self.structure.get_chains():  # Use self.structure
            for residue in chain:
                if not is_aa(residue, standard=True):  # Skip non-amino acid entities, ensure to import is_aa
                    continue
                residue_key = (chain.id, residue.resname)
                try:
                    alpha_carbon = residue['CA']  # Assuming the alpha carbon is labeled as 'CA'
                    contacts = ns.search(alpha_carbon.get_coord(), 5.0, level='A')  # Search within 10Å radius
                    total_contact_dict[residue_key] = len(contacts) - 1  # Subtract 1 to exclude the residue itself
                except KeyError:
                    total_contact_dict[residue_key] = 0  # Set total contact to 0 if alpha carbon is missing

        return total_contact_dict


    
    def extract_secondary_structure(self):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', self.pdb_file)
        model = structure[0]

        # Generate DSSP data using Biopython's DSSP wrapper
        dssp = DSSP(model, self.pdb_file, dssp=self.dssp_executable)

        # Structured data to store the features
        structured_results = {}

        # Iterate over DSSP output to populate structured_results
        for key in dssp.keys():
            res_id = (key[0], key[1][1])  # Chain ID and residue number (ignoring insertion code for simplicity)
            _, _, ss, access, phi, psi = dssp[key][:6]  # Extract the necessary data
            structured_results[res_id] = {
                'secondary_structure': ss,
                'solvent_accessibility': access,
                'phi': phi,
                'psi': psi
            }

        return structured_results
    
    def load_pocket_residues(self):
        pocket_residues_set = set()
        if self.pocket_pdb_file:
            with open(self.pocket_pdb_file, 'r') as file:
                for line in file:
                    if line.startswith('ATOM'):
                        chain = line[21]
                        residue = line[22:26].strip()
                        pocket_residues_set.add((chain, residue))
        return pocket_residues_set
    
    def one_hot_encoding(self):
        encoding_ss = {}
        encoding_aa = {}
        
        # Encode secondary structure
        for i, code in enumerate(self.secondary_structure_codes):
            one_hot_vector_ss = torch.zeros(len(self.secondary_structure_codes))
            one_hot_vector_ss[i] = 1
            encoding_ss[code] = one_hot_vector_ss
        
        # Encode amino acids
        for i, amino_acid in enumerate(self.amino_acids):
            one_hot_vector_aa = torch.zeros(len(self.amino_acids))
            one_hot_vector_aa[i] = 1 
            encoding_aa[amino_acid] = one_hot_vector_aa
        
        return encoding_ss, encoding_aa
    
    def extract_features(self):
        structure = self.parser.get_structure('protein', self.pdb_file)
        model = structure[0]  # Asumiendo un único modelo para simplificar
        # Inicializar listas para mantener las características extraídas
        dssp_data = self.extract_secondary_structure()
        all_features = []
        for chain in model.get_chains():
            for residue in chain.get_residues():
                if residue.get_id()[0] == ' ':  # Solo residuos estándar
                    residue_id = residue.get_id()
                    pdb_id = self.pdb_file.split('/')[-1].split('.')[0]
                    residue_name = residue.get_resname()
                    secondary_structure = dssp_data.get((chain.id, residue_id[1]), {}).get('secondary_structure', ' ')
                    solvent_accessibility = dssp_data.get((chain.id, residue_id[1]), {}).get('solvent_accessibility', 0)
                    solvent_accessibility = float(solvent_accessibility) if solvent_accessibility else 0
                    psi_angle = dssp_data.get((chain.id, residue_id[1]), {}).get('psi', np.nan)
                    phi_angle = dssp_data.get((chain.id, residue_id[1]), {}).get('phi', np.nan)
                    In_pocket = int((chain.id, str(residue_id[1])) in self.pocket_residues)
                    total_contact = self.calculate_total_contact().get((chain.id, residue.resname), 0)
                    amino_acid_one_hot = self.encoding_aa[residue_name]
                    secondary_structure_one_hot = [self.encoding_ss[code] for code in secondary_structure]
        
                    feature_dict = {
                        'PDB_ID': pdb_id,
                        'Chain': chain.id,
                        'Residue_ID': amino_acid_one_hot,
                        'Residue_Name': residue_name,
                        'In_Pocket': In_pocket,
                        'Secondary_structure': secondary_structure_one_hot,
                        'Solvent_accesibility': solvent_accessibility,
                        'Psi_angle': psi_angle,
                        'Phi_angle': phi_angle,
                        'Total_contact': total_contact,
                    }
                    
                    all_features.append(feature_dict)
                    
        return all_features

def matched_pdb_files(pdb_dir, pocket_dir):
    # List all .pdb files in the PDB folder, remove the extension and the prefix
    pdb_files = [os.path.basename(f) for f in glob(os.path.join(pdb_dir, '*.pdb'))]
    pdb_ids = [f[3:-4] for f in pdb_files if f.startswith('pdb') and f.endswith('.pdb')]

    # List all pocket files and remove the extension
    pocket_files = [os.path.basename(f) for f in glob(os.path.join(pocket_dir, '*_pocket.pdb'))]
    pocket_ids = [f[:-11] for f in pocket_files if f.endswith('_pocket.pdb')]

    # Find the intersection of the two sets to ensure each PDB has its pocket
    matched_ids = set(pdb_ids).intersection(set(pocket_ids))

    # Construct the full filenames for the matched files
    matched_pdb_files = [os.path.join(pdb_dir, f'pdb{id}.pdb') for id in matched_ids]
    matched_pocket_files = [os.path.join(pocket_dir, f'{id}_pocket.pdb') for id in matched_ids]

    return matched_pdb_files, matched_pocket_files


# def extract_features_for_matched_pdb_files(test_matched_files):
#     all_features_dict = {}
#     for pdb_file, pocket_pdb_file in test_matched_files:
#         pf = ProteinFeatures(pdb_file, pocket_pdb_file)
#         features_list = pf.extract_features()  # This is a list of dictionaries
#         for features in features_list:
#             # Assuming 'PDB_ID' uniquely identifies each feature set, use it as the key
#             key = features['PDB_ID']
#             all_features_dict[key] = features
#     return all_features_dict
def extract_features_for_matched_pdb_files(matched_pdb_files):
    all_features_dict = {}
    for pdb_file, pocket_pdb_file in matched_pdb_files:
        pf = ProteinFeatures(pdb_file, pocket_pdb_file)
        features_list = pf.extract_features()  # This is a list of dictionaries
        
        for feature in features_list:
            # Create a unique identifier for each residue
            unique_residue_id = f"{feature['PDB_ID']}_{feature['Chain']}_{feature['Residue_Name']}"
            
            # Append an index if the residue ID is not unique
            index = 0
            original_unique_residue_id = unique_residue_id
            while unique_residue_id in all_features_dict:
                index += 1
                unique_residue_id = f"{original_unique_residue_id}_{index}"
            
            # Now use this unique_residue_id as the key in your dictionary
            all_features_dict[unique_residue_id] = feature
    
    return all_features_dict


# Example of how to call the matched_pdb_files function
pdb_dir = '/Users/javierherranzdelcerro/Desktop/PYT_SBI/SBPYT_project/raw_pdb'
pocket_dir = '/Users/javierherranzdelcerro/Desktop/PYT_SBI/SBPYT_project/dataset/pocket_pdb'
matched_pdb_files, matched_pocket_files = matched_pdb_files(pdb_dir, pocket_dir)
matched_files = list(zip(matched_pdb_files, matched_pocket_files))
features_dict = extract_features_for_matched_pdb_files(matched_files)
print(features_dict)



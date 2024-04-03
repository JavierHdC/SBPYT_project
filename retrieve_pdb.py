"""
Script for Downloading PDB Files Based on Pocket PDB Filenames.

This script allows users to specify a directory containing pocket PDB files
and a destination directory for downloading the corresponding full PDB files
from the Protein Data Bank. It extracts the PDB IDs from the filenames of
pocket PDB files, downloads the PDB files, and saves them in the specified
destination directory.

Usage:
    python download_pdb.py -i <input_directory> -o <output_directory>

Arguments:
    -i --input_dir    Directory containing pocket PDB files.
    -o --output_dir   Directory where the PDB files will be saved.
"""

from Bio.PDB.PDBList import PDBList
import os
import argparse

def main(input_dir, output_dir):
    """
    Main function to download PDB files based on pocket PDB filenames.

    Parameters:
        input_dir (str): Path to the directory containing pocket PDB files.
        output_dir (str): Path to the directory where PDB files will be saved.
    """
    # Initialize PDBList
    pdbl = PDBList()

    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # List to hold the formatted PDB names
    pdb_names = []

    # Loop through each file in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith("_pocket.pdb"):
            # Extract the PDB identifier from the filename
            pdb_id = filename.split("_pocket.pdb")[0]
            
            # Append the formatted name to the list
            pdb_names.append(pdb_id.upper())  # PDB IDs are typically uppercase

    for pdb_id in pdb_names:
        # Retrieve and save the PDB file
        pdbl.retrieve_pdb_file(pdb_id, pdir=output_dir, file_format='pdb')

        print(f"Downloaded: {pdb_id}")

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Download PDB files based on IDs found in pocket PDB filenames.')
    parser.add_argument('--input_dir', '-i', type=str, required=True,
                        help='Directory containing pocket PDB files.')
    parser.add_argument('--output_dir', '-o', type=str, required=True,
                        help='Directory where the PDB files will be saved.')

    args = parser.parse_args()

    # Execute the main function
    main(args.input_dir, args.output_dir)

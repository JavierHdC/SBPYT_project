import json
import argparse
import torch

def tensor_to_list(obj):
    """
    Recursively convert tensors in nested objects (lists or dictionaries) to standard lists.
    This function facilitates JSON serialization of objects containing PyTorch tensors.

    Parameters:
        obj: A nested structure of dictionaries or lists containing tensors.

    Returns:
        A nested structure identical to `obj` with all tensors converted to lists.
    """
    if isinstance(obj, torch.Tensor):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {k: tensor_to_list(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [tensor_to_list(v) for v in obj]
    else:
        return obj

def save_features_to_json(features_dict, file_path):
    """
    Save the features dictionary, potentially containing tensors, to a JSON file.

    Parameters:
        features_dict (dict): Dictionary containing the features to save.
        file_path (str): Path to the JSON file where the dictionary will be saved.
    """
    with open(file_path, 'w') as f:
        json.dump(features_dict, f, indent=4, default=tensor_to_list)

def main(features_dict_path, json_file_path):
    """
    Main function to convert features from a dictionary, potentially containing tensors,
    to a JSON serializable format and save it to a file.

    Parameters:
        features_dict_path (str): Path to the Python file or method that extracts features into a dictionary.
                                  This is a placeholder parameter and would need to be replaced or extended
                                  in an actual implementation to properly load or compute the features.
        json_file_path (str): Path to the JSON file where the features will be saved.
    """
    # Placeholder for feature extraction method
    # test_features_dict = extract_features_for_matched_pdb_files(test_matched_files)
    # For the purpose of this example, we'll simulate loading or computing the features_dict
    # Here, you would replace the following line with the actual feature extraction, e.g.,
    # test_features_dict = your_feature_extraction_method(features_dict_path)
    
    test_features_dict = {}  # Assuming this is filled with your actual data

    # Save the features dictionary to a JSON file
    save_features_to_json(test_features_dict, json_file_path)

    print(f"Features saved to {json_file_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Save feature dictionary to a JSON file.")
    parser.add_argument('--features_dict_path', required=True, help="Path to the Python file or method that provides the feature dictionary.")
    parser.add_argument('--json_file_path', required=True, help="Path where the JSON file will be saved.")

    args = parser.parse_args()

    main(args.features_dict_path, args.json_file_path)

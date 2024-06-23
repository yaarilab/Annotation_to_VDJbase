import os
import shutil
import json
import pandas as pd
import argparse



def main(metadata_path, target_repo_path):

    try:
        parts = target_repo_path.split('/')
        repo_path = '/'.join(parts[:-3])
        chain = parts[-1]

        with open(metadata_path, 'r') as metadata_file:
            project_metadata = json.load(metadata_file)
        
        # Update repertoire_id by removing '_IGH', '_IGK', '_IGL'
        for repertoire in project_metadata["Repertoire"]:
            original_id = repertoire['repertoire_id']
            updated_id = original_id.replace('_IGH', '').replace('_IGK', '').replace('_IGL', '')
            repertoire['repertoire_id'] = updated_id

        # Write the updated metadata directly to the original file
        with open(metadata_path, 'w') as metadata_file:
            json.dump(project_metadata, metadata_file, indent=4)

    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some inputs.')

    # Add arguments
    parser.add_argument('metadata_filename', type=str, help='Path to the metadata file')
    parser.add_argument('target_repo_path', type=str, help='Path to the target repository')

    # Parse the arguments
    args = parser.parse_args()
    main(args.metadata_filename, args.target_repo_path)
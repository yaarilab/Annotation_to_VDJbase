import os
import shutil
import json
import pandas as pd

SEQUENCE_DATA_STORE = r'/home/bcrlab/malachy7/sequence_data_store_test/'
DIGBY_DATA_SAMPLES_PATH = r'misc/work/VDJbase/digby_data/AIRR-seq/Human/IGH'
REQUIRED_FILES = ['haplotype', 'genotype', 'ogrdb_plots.pdf', 'ogrdb_report.csv']
SPLIT = '/'

# Extracts repertoire, subject, and sample IDs from a JSON file
def get_repertoire_details(file_path):
    
    with open(file_path, 'r') as details:
        data = json.load(details)
        repertoire_id = data['repertoire_id']
        subject_id = data['subject_id']
        sample_id = data['sample_id']
    
    return repertoire_id, subject_id, sample_id

# Merges metadata from various sources into a single JSON file in the destination project
def merge_metadata(metadata_filename, project_dest, tsv_map, pre_processed_map, vdjbase_project_name, repertoire_mapping):
    with open(metadata_filename, 'r') as metadata:
        project_metadata = json.load(metadata)

        for file in tsv_map:
            repertoire_id, subject_id, sample_id = get_repertoire_details(file['repertoire_ids'])
            with open(file['annotation_metadata'], 'r') as annotation_metadata:
                annotation_metadata = json.load(annotation_metadata)
                update_annotated_metadata(project_metadata, repertoire_id, annotation_metadata)
        
        for file in pre_processed_map:
            repertoire_id, subject_id, sample_id = get_repertoire_details(file['repertoire_ids'])
            with open(file['pre_processed_metadata'], 'r') as pre_processed_metadata:
                pre_processed_metadata = json.load(pre_processed_metadata)
                update_pre_processed_metadata(project_metadata, repertoire_id, pre_processed_metadata)

        # Write the updated project_metadata to a new JSON file
        new_metadata_path = os.path.join(project_dest, f'{vdjbase_project_name}.json')
        with open(new_metadata_path, 'w') as new_metadata_file:
            json.dump(project_metadata, new_metadata_file, indent=4)

        copy_required_files(repertoire_mapping, tsv_map, project_dest)

def copy_required_files(repertoire_mapping, tsv_map, project_dest):
    for vdjbase_project in repertoire_mapping:
        vdjbase_project = repertoire_mapping[vdjbase_project]
        vdjbase_project_path = os.path.join(project_dest, 'samples', vdjbase_project['project_number'], vdjbase_project['vdjbase_name'])
        if not os.path.exists(vdjbase_project_path):
            os.makedirs(vdjbase_project_path)
        
        for projcet in tsv_map:
            repertoire_id = None
            with open(projcet['repertoire_ids'], 'r') as repertoire:
                repertoire_data = json.load(repertoire)
                repertoire_id = repertoire_data['repertoire_id']

            if repertoire_id == vdjbase_project['airr_repertoire_id']:
                for file in projcet['required_files']:
                    file_name = file.split(SPLIT)[-1]
                    destination_path = os.path.join(vdjbase_project_path, file_name)
                    copy_file(source_path=file, destination_path= destination_path)
        

    

# Updates project metadata with annotated metadata for a specific repertoire
def update_annotated_metadata(project_metadata, repertoire_id, annotation_metadata):
    new_data = annotation_metadata['sample']['data_processing']
    for repertoire in project_metadata['Repertoire']:
        if repertoire['repertoire_id'] == repertoire_id:
            original_data = repertoire['data_processing'][0]
            repertoire['data_processing'][0] = merge_json_data_recursive(original_data, new_data)

# Updates project metadata with pre-processed metadata for a specific repertoire
def update_pre_processed_metadata(project_metadata, repertoire_id, pre_processed_metadata):
    new_data = pre_processed_metadata['sample']['data_processing']
    for repertoire in project_metadata['Repertoire']:
        if repertoire['repertoire_id'] == repertoire_id:
            original_data = repertoire['data_processing'][0]
            repertoire['data_processing'][0] = merge_json_data_recursive(original_data, new_data)


def merge_json_data_recursive(original_data, new_data):
    """
    Recursively merges new_data into original_data. If a key in new_data already exists in original_data
    and both values are dictionaries, it merges them recursively. If both are lists, it appends the items
    from the new list to the old list. Otherwise, the value in original_data is updated with the value from new_data.
"""
    for key, value in new_data.items():
        if key in original_data:
            if isinstance(original_data[key], dict) and isinstance(value, dict):
                merge_json_data_recursive(original_data[key], value)
            elif isinstance(original_data[key], list) and isinstance(value, list):
                original_data[key].extend(value)
            else:
                original_data[key] = value
        else:
            original_data[key] = value

    return original_data


# Copies content from a source directory to a destination directory and merges metadata
def copy_folder_content(source_folder, target_repo_path, vdjbase_project_name,project_number, metadata_filename, repertoire_mapping):
    # Create the destination directory if it does not exist
    if not os.path.exists(target_repo_path):
        os.makedirs(target_repo_path)
    
    tsv_files_paths, pre_processed_files = find_project_tsv_files(source_folder)

    merge_metadata(metadata_filename, target_repo_path, tsv_files_paths, pre_processed_files, vdjbase_project_name, repertoire_mapping)


def copy_file(source_path, destination_path):
    """
    Copy a file from source_path to destination_path.

    """
    # Check if the source file exists
    if not os.path.isfile(source_path):
        raise FileNotFoundError(f"The source file does not exist: {source_path}")
    
    # Ensure the destination directory exists, if not, create it
    os.makedirs(os.path.dirname(destination_path), exist_ok=True)
    
    # Copy the file
    shutil.copy2(source_path, destination_path)
    print(f"File copied from {source_path} to {destination_path}")


# Finds TSV files and pre-processed files within a project directory
def find_project_tsv_files(project_path):
    pre_processed_folders = []
    try:
        annotated_folder_path = os.path.join(project_path, 'annotated')
        annotated_folders = os.listdir(annotated_folder_path)
        pre_processed_folder_path = os.path.join(project_path, 'pre_processed')
        if os.path.exists(pre_processed_folder_path):
            pre_processed_folders = os.listdir(pre_processed_folder_path)
            
        tsv_files = start_scan(annotated_folder_path,annotated_folders , False)
        pre_processed_files = start_scan(pre_processed_folder_path, pre_processed_folders , True)

    except Exception as e:
        print(e)

    return tsv_files, pre_processed_files

# Initiates scanning of folders for TSV and metadata files
def start_scan(folder_path,folders, pre_processed):
    files = []
    for subject in folders:
        subject_path = os.path.join(folder_path, subject)
        res = scan_subject_folder(subject_path, pre_processed)
        for file in res:
            files.append(file)
    
    return files

# Scans a subject folder for TSV and metadata files
def scan_subject_folder(subject_path, pre_processed):
    tsv_files = []
    samples = os.listdir(subject_path)
    for sample in samples:
        sample_path = os.path.join(subject_path, sample)
        files = scan_run_folder(sample_path, pre_processed)
        for file in files:
            tsv_files.append(file)
    
    return tsv_files

# Scans a run folder for TSV and metadata files
def scan_run_folder(sample_path, pre_processed):
    tsv_files = []
    runs = os.listdir(sample_path)
    for run in runs:
        run_path = os.path.join(sample_path, run)
        run_results = os.listdir(run_path)
        for result in run_results:
            result_path = os.path.join(run_path, result)
            if not pre_processed:
                file = find_tsv_and_metadata_for_annotated(result_path)
            else:
                file = find_metadata_for_pre_processed(result_path)

            if file != None:
                tsv_files.append(file[0])
    
    return tsv_files

# Finds metadata for pre-processed results
def find_metadata_for_pre_processed(result_path):
    res_list = []
    res = {
        'repertoire_ids': None,
        'pre_processed_metadata': None
    }
    result_folders = os.listdir(result_path)
    for folder in result_folders:
        folder_path = os.path.join(result_path, folder)
        folder_files = os.listdir(folder_path)
        
        if 'meta_data' in folder:
            if 'pre_processed_metadata.json' in folder_files:
                res['pre_processed_metadata'] = os.path.join(folder_path, 'pre_processed_metadata.json')
            
            if 'repertoire_id.json' in folder_files:
                res['repertoire_ids'] = os.path.join(folder_path, 'repertoire_id.json')

    check_result_fileds(res, result_path)
    if all(value is not None for value in res.values()):
        res_list.append(res)
        return res_list
    
    return None     
                

# Finds TSV files and their corresponding metadata for annotated results
def find_tsv_and_metadata_for_annotated(result_path):
    res_list = []
    res = {
            'file_path': None,
            'file_name': None,
            'repertoire_ids': None,
            'annotation_metadata': None,
            'required_files' : []
        }
    
    result_folders = os.listdir(result_path)
    for folder in result_folders:
        folder_path = os.path.join(result_path, folder)
        folder_files = os.listdir(folder_path)
        for file in folder_files:
            if 'Finale' in file:
                res['file_path'] = os.path.join(folder_path, file)
                res['file_name'] = file
            
            if file == 'repertoire_id.json':
                res['repertoire_ids'] = os.path.join(folder_path, file)
                
            for required_file in REQUIRED_FILES:
                if required_file in file:
                    res['required_files'].append(os.path.join(folder_path, file))


        if 'meta_data' in folder:
            if 'annotation_metadata.json' in folder_files:
                res['annotation_metadata'] = os.path.join(folder_path, 'annotation_metadata.json')

    check_result_fileds(res, result_path)
    if all(value is not None for value in res.values()):
        res_list.append(res)
        return res_list
    
    return None 

# Checks if all required fields in a result are present
def check_result_fileds(result, folder):
    for key, value in result.items():
        if value == None:
            print(f"{key} was not found in the {folder}")



def verify_directory_exists(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Directory does not exist: {path}")

def verify_and_clear_project_directory(target_repo_path, project_name):
    airr_correspondence_path = os.path.join(target_repo_path, 'airr_correspondence.csv')
    if not os.path.isfile(airr_correspondence_path):
        raise FileNotFoundError(f"airr_correspondence.csv not found in the target repo at {target_repo_path}")

    # Read the CSV file into a DataFrame
    correspondence_df = pd.read_csv(airr_correspondence_path)

    # Check if the project name is in the airr_file column and extract the first matching vdjbase_name
    vdjbase_project_name = None
    for airr_file in correspondence_df['airr_file']:
        if project_name in airr_file:
            vdjbase_project_name = correspondence_df[correspondence_df['airr_file'] == airr_file]['vdjbase_name'].str.split('_').str[0].iloc[0]
            break

    if vdjbase_project_name is None:
        raise ValueError(f"No matching project name {project_name} found in airr_correspondence.csv")

    # Path to the project directory within the target repo
    project_dir_path = os.path.join(target_repo_path, 'samples' ,vdjbase_project_name)
    
    if not os.path.isdir(project_dir_path):
        raise FileNotFoundError(f"Project directory {project_dir_path} does not exist")

    # Clear the project directory
    for filename in os.listdir(project_dir_path):
        file_path = os.path.join(project_dir_path, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print(f"Failed to delete {file_path}. Reason: {e}")

    print(f"Project directory {project_dir_path} has been cleared.")
    return airr_correspondence_path


def derive_vdjbase_project_mapping(airr_correspondence_path, project_name):
    correspondence_df = pd.read_csv(airr_correspondence_path)
    
    # Filter the DataFrame for rows where the vdjbase_name contains the specific project name
    project_specific_df = correspondence_df[correspondence_df['airr_file'].str.contains(f"{project_name}", regex=True)]

    mapping = {}
    for index, row in project_specific_df.iterrows():
        vdjbase_name = row['vdjbase_name']
        airr_repertoire_id = row['airr_repertoire_id']
        
        # Extract individual and sample from vdjbase_name
        project_number, individual, sample = vdjbase_name.split('_')

        mapping[vdjbase_name] = {
            'project_name' : project_name,
            'vdjbase_name' : vdjbase_name, 
            'project_number': project_number,
            'individual': individual,
            'sample': sample,
            'airr_repertoire_id': airr_repertoire_id
        }

    return mapping


def verify_annotations_exist(sequence_data_store_path, airr_correspondence_path, project_name):
    # Read the airr_correspondence.csv to get the expected repertoire IDs and corresponding file patterns
    correspondence_df = pd.read_csv(airr_correspondence_path)
    filtered_df = correspondence_df[correspondence_df['airr_file'].str.contains(project_name)]
    expected_repertoires = filtered_df['airr_repertoire_id'].tolist()
    #expected_repertoires = correspondence_df['airr_repertoire_id'].to_list()
    
    # Dictionary to hold the found 'Final' files
    found_files = []

    # Walk through the annotated directory to find 'Final' files
    annotated_path = os.path.join(sequence_data_store_path, 'annotated')
    if not os.path.isdir(annotated_path):
        raise FileNotFoundError(r"there is no annotated file for {sequence_data_store_path}")
        

    for root, dirs, files in os.walk(annotated_path):
        for file in files:
            if "Final" in file:
                found_files.append(file)

    # Verify that each expected file pattern has at least one matching 'Final' file
    missing_annotations = []
    for repertoire in expected_repertoires:
        if not any(repertoire in f for f in found_files):
            missing_annotations.append(repertoire)

    if missing_annotations:
        raise FileNotFoundError(f"Missing 'Final' annotation files for the following patterns: {', '.join(missing_annotations)}")
    
    return True  # Return True if all checks pass


def copy_files(source, destination, file_names):
    for file_name in file_names:
        source_file = os.path.join(source, file_name)
        if os.path.isfile(source_file):
            shutil.copy(source_file, destination)
        else:
            print(f"Required file not found: {source_file}")
            # If the file is required, you may want to raise an exception here instead of just printing.

def consolidate_metadata(repertoire_metadata_file, additional_metadata_paths):
    with open(repertoire_metadata_file, 'r') as f:
        metadata = json.load(f)
    
    for metadata_path in additional_metadata_paths:
        if os.path.isfile(metadata_path):
            with open(metadata_path, 'r') as f:
                additional_metadata = json.load(f)
            metadata.update(additional_metadata)  # Assumes that additional metadata should overwrite
    
    return metadata

def main(project_name, source_folder, metadata_filename, target_repo_path):
    try:
        # Verify that the source folder and target repo path exist
        verify_directory_exists(source_folder)
        verify_directory_exists(target_repo_path)

        # Verify airr_correspondence.csv file exists and get the mapping
        airr_correspondence_path = verify_and_clear_project_directory(target_repo_path, project_name) #need to add the check of contains references to a file matching the project, in the airr_file column.
        repertoire_mapping = derive_vdjbase_project_mapping(airr_correspondence_path, project_name)
        
        # Verify annotations exist
        verify_annotations_exist(source_folder, airr_correspondence_path, project_name)
        first_key = next(iter(repertoire_mapping.keys()))# Access the first key
        project_number = repertoire_mapping[first_key]['project_number']
        vdjbase_project_name = project_number + ('_' + project_name)
        copy_folder_content(source_folder, target_repo_path, vdjbase_project_name, project_number, metadata_filename, repertoire_mapping)

        print("Data copy completed successfully.")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # This is where you would parse command line arguments to the script
    # For example:
    # project_name = sys.argv[1]
    # source_folder = sys.argv[2]
    # metadata_filename = sys.argv[3]
    # target_repo_path = sys.argv[4]

    # Hardcoded for demonstration purposes
    project_name = r"PRJNA248411"
    source_folder = r"/home/bcrlab/malachy7/sequence_data_store_test/PRJNA248411/runs/current/"
    metadata_filename = r"/home/bcrlab/malachy7/sequence_data_store_test/PRJNA248411/project_metadata/metadata.json"
    target_repo_path = r"/home/bcrlab/malachy7/digby_data/AIRR-seq/Human/IGH/"

    main(project_name, source_folder, metadata_filename, target_repo_path)

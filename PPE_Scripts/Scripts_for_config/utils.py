
def txt_to_csv(input_file_path, output_file_path):
    import pandas as pd
    """
    Converts a space-separated text file into a CSV file.

    Parameters:
    input_file_path (str): Path to the input text file.
    output_file_path (str): Path to the output CSV file.
    """
    # Read the data from the text file
    with open(input_file_path, 'r') as file:
        lines = file.readlines()
    
    # Extract the header and data
    header = lines[0].strip().split()
    data = [line.strip().split() for line in lines[1:]]
    
    # Create a DataFrame
    df = pd.DataFrame(data, columns=header)
    
    # Remove the last three columns
    df = df.iloc[:, :-3]
    
    # Save the DataFrame to a CSV file
    df.to_csv(output_file_path, index=False)
    
    print(f"Data has been successfully converted and saved to {output_file_path}")


def calculate_average_distance(ensemble_data):
    """
    Calculate the average Euclidean distance of each ensemble's parameters to those of other ensembles.
    
    Parameters:
    ensemble_data (numpy.ndarray): A 2D array where each row represents an ensemble and each column represents a parameter.
    
    Returns:
    numpy.ndarray: A 1D array containing the average distance for each ensemble.
    """
    num_ensembles, num_parameters = ensemble_data.shape
    average_distances = np.zeros(num_ensembles)
    
    # Compute pairwise Euclidean distances
    for j in range(num_ensembles):
        distances = []
        for m in range(num_ensembles):
            if j != m:
                # Calculate the Euclidean distance between ensemble j and ensemble m
                dist = np.sqrt(np.sum((ensemble_data[j] - ensemble_data[m]) ** 2))
                distances.append(dist)
        
        # Average distance for ensemble j
        average_distances[j] = np.mean(distances)
    
    # Normalize by the number of parameters and ensembles
    normalized_distances = average_distances / (num_parameters * (num_ensembles))
    
    return normalized_distances

def create_param_ranges(file_path):
    param_ranges = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.split()
            if len(parts) == 3:
                param_name = parts[0]
                param_range = (float(parts[1]), float(parts[2]))
                param_ranges[param_name] = param_range
    return param_ranges

def PPE_values_modify(input_file_path, parameters_file):
    ''' 
    Opens the PPE_values.txt file and modifies the header to make sure the parameters you are
    changing are all (and only them) and in the header for perturbing. 
    '''
    # Read the parameter names from the parameters_file
    with open(parameters_file, 'r') as param_file:
        # Extract the parameter names from the file
        parameter_names = [line.split()[0] for line in param_file.read().strip().split('\n')]
    # Read the first line of the PPE_values.txt file
    with open(input_file_path, 'r') as file:
        lines = file.readlines()
    # Split the first line into a list of strings
    header = lines[0].strip().split()
    header = parameter_names
    # Set last 3 values to 1 
    header.extend(['INSTALL', 'START', 'END'])
    # Include the name of the simulation here
    header.insert(0, '0')
    # Join the modified header back into a single string
    new_header = ' '.join(header)
    
    # Replace the first line with the new header
    lines[0] = new_header + '\n'
    
    # Write the modified lines back to the file
    with open(input_file_path, 'w') as file:
        file.writelines(lines)
    
    print("Header updated successfully.")


def check_column_duplicates(file_path, header_array):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Read the data into a 2D list
    data = [line.strip().split() for line in lines]
    
    # Transpose the data to work with columns
    columns = list(zip(*data))
    #colum = len(line.strip().split())
    num_columns = len(data[0]) if columns else 0

    # Determine indices to skip
    skip_indices = set(range(1)) | set(range(num_columns - 3, num_columns))

    for col_index, column in enumerate(columns):
        if col_index in skip_indices:
            continue  # Skip the first and last 3 columns

        duplicates = set([val for val in column if column.count(val) > 1])
        if duplicates:
            if col_index > num_columns:
                pass
            else:
                print(f"Duplicate values found in {header_array[col_index]}: {duplicates}")
                print(f"STOP AND CHANGE LINE 55 in the bash_script.py to increased decimal points.")

        else:
            print(f"No duplicate values found in {header_array[col_index]} columns")
            pass

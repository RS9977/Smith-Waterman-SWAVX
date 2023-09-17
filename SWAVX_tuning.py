import subprocess
import re
import random
import tempfile
import os
import multiprocessing
import time
import math
import numpy as np
import pickle

def parse_gcc_params(output, best_param, i, numIter, alpha=10, beta=30):
    params = {}
    param_pattern = r"--param=([\w-]+)=<(-?\d+),(-?\d+)>"
    exceptions_param = ['lazy-modules', 'logical-op-non-short-circuit', 'ranger-debug']
    for match in re.finditer(param_pattern, output):
        param_name = match.group(1)
        if param_name not in exceptions_param:
            a = int(match.group(2))
            b = int(match.group(3))
            # Initialize parameter with a random value within the specified range
            if param_name not in best_param or random.randint(0, 100)<alpha:
                params[param_name] = random.randint(a, b)
            else:
                value = random.randint(-1, 1)*math.floor(math.exp(-3*i/numIter)*(b-a)/beta) + best_param[param_name] 
                if value > a and value < b:
                    params[param_name] = value
                else:
                    params[param_name] = best_param[param_name] 

    # For parameters with a list of possible values, choose a random value from the list
    list_pattern = r"--param=([\w-]+)=\[((?:\w+\|)+\w+)\]?"
    for match in re.finditer(list_pattern, output):
        param_name = match.group(1)
        if param_name not in exceptions_param:
            possible_values = match.group(2).split('|')
            # Choose a random value from the list of possible values
            if param_name not in params:
                if param_name not in best_param or random.randint(0, 100)<math.floor(math.exp(-3*i/numIter)*(alpha)):
                    params[param_name] = random.choice(possible_values)
                else:
                    params[param_name] = best_param[param_name]


    # For parameters without specified ranges, assign a random value only if not in param_names
    default_pattern = r"--param=([\w-]+)="
    for match in re.finditer(default_pattern, output):
        param_name = match.group(1)
        if param_name not in exceptions_param:
            # Check if the parameter has already been assigned a value from param_pattern
            if param_name not in params:
                if param_name not in best_param or random.randint(0, 100)<alpha:
                    params[param_name] = random.randint(0, 100000)
                else:
                    if isinstance(best_param[param_name] , str):
                        print(best_param[param_name])
                    value = random.randint(-1, 1)*math.floor(math.exp(-3*i/numIter)*(b-a)/beta) + best_param[param_name] 
                    if value > 0 and value < 100000:
                        params[param_name] = value
                    else:
                        params[param_name] = best_param[param_name] 

    return params

def compile_with_gcc(gcc_params, selected_indices, i=-1, CorGorI=2):
    try:
        # Build the GCC command with parameters
        
        if i>-1:
            param_cnt = 0
            if CorGorI == 1:
                gcc_command = ["gcc", "-O3", "-mavx2", "-D", "SUBMAT", "-D", "L8", "-lpthread", "SWAVX_utils.c", "SWAVX_SubMat.c", "SWAVX_TOP_datasets_MultiThread.c", "SWAVX_256_Func_SeqToSeq_SubMat.c", "-o", "SWAVX"]
            else:
                gcc_command = ["gcc", "-O3", "-mavx2", "-D", "SUBMAT", "-D", "L8", "-lpthread", "SWAVX_utils.c", "SWAVX_SubMat.c", "SWAVX_TOP_datasets_MultiThread.c", "SWAVX_256_Func_SeqToSeq_SubMat.c", "-S"]
            #gcc_command = ["make", "AVX2_PAR", "VAR=L8 SUBMAT", "PAR="]
            for param_name, param_value in gcc_params.items():
                if param_cnt not in selected_indices:
                    continue
                param_cnt += 1
                gcc_command.append(f"--param={param_name}={param_value}")
            
            # Run GCC to compile the program
            subprocess.check_call(gcc_command)
            '''
            if CorGorI == 1:
                print(f"Compilation {i} successful. Binary SWAVX generated.")
            else:
                print(f"Compilation {i} successful. Assemblies generated.")
            '''
            return 1
        else:
            param_cnt = 0
            gcc_command = ["gcc", "-O3", "-mavx2", "-D", "SUBMAT", "-D", "L8", "-lpthread", "SWAVX_utils.c", "SWAVX_SubMat.c", "SWAVX_TOP_datasets_MultiThread.c", "SWAVX_256_Func_SeqToSeq_SubMat.c", "-o", "SWAVX"]
            #gcc_command = ["make", "AVX2_PAR", "VAR=L8 SUBMAT", "PAR="]
            for param_name, param_value in gcc_params.items():
                if param_cnt not in selected_indices:
                    continue
                param_cnt += 1
                gcc_command.append(f"--param={param_name}={param_value}")
            
            # Run GCC to compile the program
            subprocess.check_call(gcc_command)
            param_cnt = 0
            gcc_command = ["gcc", "-O3", "-mavx2", "-D", "SUBMAT", "-D", "L8", "-lpthread", "SWAVX_utils.c", "SWAVX_SubMat.c", "SWAVX_TOP_datasets_MultiThread.c", "SWAVX_256_Func_SeqToSeq_SubMat.c", "-S"]
            #gcc_command = ["make", "AVX2_PAR", "VAR=L8 SUBMAT", "PAR="]
            for param_name, param_value in gcc_params.items():
                if param_cnt not in selected_indices:
                    continue
                param_cnt += 1
                gcc_command.append(f"--param={param_name}={param_value}")
            
            # Run GCC to compile the program
            subprocess.check_call(gcc_command)
            print("Compilation successful. Binary SWAVX and assemblies generated.")
            return 1

    except subprocess.CalledProcessError as e:
        print(f"Error during compilation: {i}")
        print(e.output)
        return 0






def compute_mean_and_filter(data_list):
    # Convert the input list to a NumPy array for efficient operations
    data_array = np.array(data_list)

    # Calculate the mean and standard deviation of the data
    mean = np.mean(data_array)
    std_dev = np.std(data_array)

    # Define the range for values to keep
    lower_bound = mean - std_dev
    upper_bound = mean + std_dev

    # Filter the values that are within the specified range
    filtered_data = data_array[(data_array >= lower_bound) & (data_array <= upper_bound)]

    # Calculate the new mean of the filtered data
    new_mean = np.mean(filtered_data)

    return new_mean


def get_gcups_from_command(command, numIter):
    GCUPS = []
    try:
        for i in range(numIter):
            # Run the command and capture its output
            result = subprocess.run(command, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            # Check if the command was successful (return code 0)
            if result.returncode == 0:
                # Search for the line containing "GCUPS: <number>"
                match = re.search(r'GCUPS:\s+(\d+\.\d+)', result.stdout)
                if match:
                    gcups = float(match.group(1))
                    GCUPS.append(gcups)
                else:
                    return None  # No matching line found
            else:
                # Command failed, return None
                return None
        return compute_mean_and_filter(GCUPS)
    except Exception as e:
        # Handle exceptions, e.g., command not found or other errors
        print(f"An error occurred: {str(e)}")
        return None




def count_instructions_in_directory(directory):
    # Initialize a dictionary to store instruction frequencies
    instruction_count = {}

    # Initialize a variable to store the overall number of instructions
    total_instructions = 0

    try:
        # Iterate through all files in the directory
        for filename in os.listdir(directory):
            if filename.endswith(".s"):
                file_path = os.path.join(directory, filename)

                with open(file_path, 'r') as file:
                    for line in file:
                        # Remove leading and trailing whitespace and split the line into words
                        words = line.strip().split()

                        # Check if the line starts with a dot (.) or is empty
                        if not words or words[0].startswith('.') or words[0][-1]==':':
                            continue  # Ignore this line

                        # Extract the first word (instruction) from the line
                        instruction = words[0]

                        # Update the instruction count dictionary
                        if instruction in instruction_count:
                            instruction_count[instruction] += 1
                        else:
                            instruction_count[instruction] = 1

                        # Increment the overall instruction count
                        total_instructions += 1

        # Return the instruction count dictionary and total number of instructions
        return instruction_count, total_instructions

    except FileNotFoundError:
        print(f"Directory '{directory}' not found.")
        return {}, 0

def get_asm_info(input_directory):
    throughput_dict = {}

    # Find all "*.s" files in the input directory
    #s_files = glob.glob(input_directory + "/*.s")
    s_files = ['./SWAVX_256_Func_SeqToSeq_SubMat.s']
    for s_file in s_files:
        command = f"llvm-mca {s_file}"
        
        try:
            output = subprocess.check_output(command, shell=True, text=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print(f"Error running the command for {s_file}: {e}")
            continue

        pattern = r"Total Cycles:\s+(\d+)"
        match = re.search(pattern, output)

        if match:
            throughput = float(match.group(1))
            throughput_dict[s_file] = throughput
        else:
            print(f"Pattern not found in the command output for {s_file}")

    return throughput_dict[s_files[0]]

def save_dictionary_to_file(dictionary, filename):
    with open(filename, 'wb') as file:
        pickle.dump(dictionary, file)

def load_dictionary_from_file(filename):
    with open(filename, 'rb') as file:
        dictionary = pickle.load(file)
    return dictionary



def main(CorGorI, numPar, numOutIter, numIter, numTest):
    try:
        print('Without Tuning:')        
        try:
            [gcc_params_min,selected_indices_min] = load_dictionary_from_file('Par_Val')
        except:
            selected_indices_min   = []
            gcc_params_min         = {}
        compile_with_gcc(gcc_params_min, selected_indices_min, -1, CorGorI)    
        GCUPS_max = get_gcups_from_command("./SWAVX", numTest)
        cycles_min = get_asm_info('./')
        instruction_count, total_instructions_min = count_instructions_in_directory('./')
        print(f"Initial GCUPS : {GCUPS_max}")
        print(f"Initial total Instructions: {total_instructions_min}")
        print(f"Initial cycles: {cycles_min}")
        #gcc_command = ["gcc", "-O3", "-mavx2", "-D", "SUBMAT", "-D", "L8", "-lpthread", "SWAVX_utils.c", "SWAVX_SubMat.c", "SWAVX_TOP.c", "SWAVX_256_Func_SeqToSeq_SubMat.c", "-o", "SWAVX"]
        #subprocess.check_call(gcc_command)
        #GCUPS_max = get_gcups_from_command('./SWAVX',numTest)
        for j in range(numOutIter):
            print(f"j: {j}")
            for i in range(numIter):
                '''
                if j==numOutIter-1:
                    CorGorI = 1
                '''
                # Run `gcc --help=params` and capture the output
                output = subprocess.check_output(["gcc", "--help=params"], stderr=subprocess.STDOUT, universal_newlines=True)
                # Parse the output to extract parameter names and descriptions
                
                gcc_params = parse_gcc_params(output, gcc_params_min, i, numIter)


                #for i in range(100):
                pool = multiprocessing.Pool(processes=1)
                selected_indices = random.sample(range(0, 270), numPar)
                selected_indices.sort()
                # Use apply_async to run the function in a separate process
                result = pool.apply_async(compile_with_gcc, (gcc_params, selected_indices, i, CorGorI))

                result_value = 0
                try:
                    # Get the result with a timeout of 2 seconds
                    result_value = result.get(timeout=100)
                except multiprocessing.TimeoutError:
                    # Handle the case where the function exceeded the timeout
                    print("Function execution took too long and was terminated.")
                    pool.close()
                    pool.join()
                    continue
                
                # Close the pool
                pool.close()
                pool.join()

                if result_value:
                    if CorGorI == 2:
                        gcups = get_gcups_from_command("./SWAVX", numTest)
                        if GCUPS_max < gcups:
                            print(i, 'of',j, ':= ', GCUPS_max, '->',gcups)
                            GCUPS_max = gcups
                            gcc_params_min = gcc_params
                            selected_indices_min = selected_indices
                            save_dictionary_to_file([gcc_params_min,selected_indices_min],'Par_Val_temp')
                            
                    elif CorGorI == 2:
                        instruction_count, total_instructions = count_instructions_in_directory('./')
                        if total_instructions < total_instructions_min:
                            print(i,'of',j, ':= ', total_instructions_min, '->',total_instructions)
                            total_instructions_min = total_instructions
                            gcc_params_min = gcc_params
                            selected_indices_min = selected_indices
                            save_dictionary_to_file([gcc_params_min,selected_indices_min],'Par_Val_temp')
                    else:
                        cycles = get_asm_info('./')
                        if cycles < cycles_min:
                            print(i,'of',j, ':= ', cycles_min, '->',cycles)
                            cycles_min = cycles
                            gcc_params_min = gcc_params
                            selected_indices_min = selected_indices
                            save_dictionary_to_file([gcc_params_min,selected_indices_min],'Par_Val_temp')
        
        print('With Tuning:')
        compile_with_gcc(gcc_params_min, selected_indices_min, -1, CorGorI)    
        GCUPS_max = get_gcups_from_command("./SWAVX", numTest)
        cycles_min = get_asm_info('./')
        instruction_count, total_instructions_min = count_instructions_in_directory('./')
        print(f"Final GCUPS : {GCUPS_max}")
        print(f"Final total Instructions: {total_instructions_min}")
        print(f"Final cycles: {cycles_min}")
        save_dictionary_to_file([gcc_params_min,selected_indices_min],'Par_Val')
    except subprocess.CalledProcessError as e:
        print("Error running 'gcc --help=params':")
        print(e.output)

if __name__ == "__main__":
    main(0, 265, 5, 150, 100)

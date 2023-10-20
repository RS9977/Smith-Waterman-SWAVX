import subprocess
import random
import ax
from ax import ParameterType, RangeParameter, ChoiceParameter
from ax.modelbridge import Models
import numpy as np
import re
from ax.service.ax_client import AxClient
from ax.service.ax_client import AxClient, ObjectiveProperties
import logging

def parse_gcc_params():
    params = {}
    param_pattern = r"--param=([\w-]+)=<(-?\d+),(-?\d+)>"
    exceptions_param = ['lazy-modules', 'logical-op-non-short-circuit', 'ranger-debug', 'lto-min-partition', 'lto-max-partition', 'vect-max-peeling-for-alignment']
    output = subprocess.check_output(["gcc", "--help=params"], stderr=subprocess.STDOUT, universal_newlines=True)
    for match in re.finditer(param_pattern, output):
        param_name = match.group(1)
        if param_name not in exceptions_param:
            a = int(match.group(2))
            b = int(match.group(3))
            # Initialize parameter with a random value within the specified range
            params[param_name] = [a,b]

    # For parameters with a list of possible values, choose a random value from the list
    list_pattern = r"--param=([\w-]+)=\[((?:\w+\|)+\w+)\]?"
    for match in re.finditer(list_pattern, output):
        param_name = match.group(1)
        if param_name not in exceptions_param:
            possible_values = match.group(2).split('|')
            # Choose a random value from the list of possible values
            if param_name not in params:
                params[param_name] = possible_values


    # For parameters without specified ranges, assign a random value only if not in param_names
    default_pattern = r"--param=([\w-]+)="
    for match in re.finditer(default_pattern, output):
        param_name = match.group(1)
        if param_name not in exceptions_param:
            # Check if the parameter has already been assigned a value from param_pattern
            if param_name not in params:
                params[param_name] = [0,1e5]

    return params

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


def get_gcups_from_command(command, numTest):
    GCUPS = []
    try:
        for i in range(numTest):
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

def get_size_info(binary_name):
    command = f"du -b {binary_name}"
        
    try:
        output = subprocess.check_output(command, shell=True, text=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print(f"Error running the command for {binary_name}: {e}")

    try:
        return int(output.split()[0])
    except:
        print(f"Pattern not found in the command output for {binary_name}")


# Function to compile the code with a set of parameters and return the runtime
def compile_and_evaluate(parameters, output_binary="SWAVX_besys_tuned", numThreads=32, numTest=20):
    #included_params = [param for param, include in parameters if include]
    
    # Build the compiler command with the selected parameters
    compiler_command = "gcc -O3 -mavx2 -lpthread -flto -D L8 SWAVX_utils.c SWAVX_SubMat.c SWAVX_TOP_datasets_MultiThread.c SWAVX_256_Func_SeqToSeq_SubMat.c -o " + output_binary
    exceptions_param = ['lazy-modules', 'logical-op-non-short-circuit', 'ranger-debug']
    for param_name, param_value in parameters.items():
        if param_name not in exceptions_param:
            compiler_command += f" --param={param_name}={param_value}"
    
    # Compile the code and measure runtime
    try:
        process = subprocess.Popen(compiler_command, shell=True)
        process.wait()
    except Exception as e:
        return {"GCUPS": (-1, 0.0)}  # Return negative infinity in case of compilation failure

    # Measure runtime using your method
    #GCUPS = get_gcups_from_command(f"./{output_binary} 'test2.fasta' 'test4.fasta' {numThreads}", numTest)
    size = get_size_info(output_binary)

    return {"size": (size, 0.0)}#, "GCUPS": (GCUPS, 0.0)}


# Function to create random parameter settings with inclusion flags
def define_parameter_space(param_dict):
    parameters = []
    for param_name, param_range in param_dict.items():
        if not isinstance(param_range[0], int):
            parameters.append({"name":param_name, "type": "choice", "values":param_range, "value_type": "str"})
        else:
            parameters.append({"name":param_name, "type": "range", "bounds":param_range, "value_type": "int"})
    return parameters

# Function to optimize the compiler parameters
def optimize_compiler_parameters(numIter=10, output_binary="SWAVX_besys_tuned", numThreads=32, numTest=20, GCUPS_init=22):
    # Define the parameter space
    param_dict = parse_gcc_params()
    parameters = define_parameter_space(param_dict)
    #the = GCUPS_init*0.95

    ax_client = AxClient(verbose_logging = False)
    ax_client.create_experiment(
        name="GCC hyperparameters tuning",
        parameters=parameters,
        objectives={"size": ObjectiveProperties(minimize=True)},
        #outcome_constraints=[f"GCUPS >= {(GCUPS_init*0.95)}"],
    )

    for _ in range(numIter):
        parameters, trial_index = ax_client.get_next_trial()
        ax_client.complete_trial(trial_index=trial_index, raw_data=compile_and_evaluate(parameters, output_binary=output_binary, numThreads=numThreads, numTest=numTest))

    best_parameters, metrics = ax_client.get_best_parameters()

    return best_parameters

# Number of optimization iterations
def Main(numIter=10, output_binary="SWAVX_besys_tuned", numThreads=20, numTest=20):
    GCUPS_init = compile_and_evaluate({}, output_binary=output_binary+"_init", numThreads=numThreads, numTest=numTest)
    print("GCUPS initial: ", GCUPS_init)
    best_parameters = optimize_compiler_parameters(numIter=numIter, output_binary=output_binary, numThreads=numThreads, numTest=numTest, GCUPS_init=GCUPS_init)
    GCUPS_final = compile_and_evaluate(best_parameters, output_binary=output_binary, numThreads=numThreads, numTest=numTest)
    print("GCUPS final: ", GCUPS_final)


if __name__ == "__main__":
    Main(numIter=15, output_binary="SWAVX_besys_tuned", numThreads=20, numTest=1)

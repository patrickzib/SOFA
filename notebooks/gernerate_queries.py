import os
import numpy as np
import pandas as pd

import sys
sys.path.append("../")

def add_gaussian_noise(data, data_type, mean=0, std_dev=0.1):
    noise = np.random.normal(0, np.std(data)* std_dev, size=len(data))
    #print(noise)
    #print(data.astype(np.float32))
    #print(data.astype(np.float32) + noise)
    return (data.astype(np.float64) + noise) #.astype(data_type)

def add_noise(input_file, length, data_type, noise_level=0.1):    
    # File paths
    output_file = input_file +"_noise_"+str(noise_level).replace(".","")
    input_file = path+"/"+input_file+".bin"
    output_file = path+"/generated/"+output_file+".bin"
        
    try:
        with open(input_file, "rb") as f:
            data = np.fromfile(f, dtype=data_type, count=1024*1024)            
    except FileNotFoundError:
        print(f"Input file '{input_file}' not found.")
        return

    
    # Add Gaussian noise
    noisy_data = add_gaussian_noise(data, data_type, std_dev = noise_level)
    noisy_data = np.ceil(noisy_data).astype(data_type)
    # print (noisy_data, noisy_data.dtype)
    # print (noisy_data, noisy_data.dtype)
    
    # Write noisy data to output file
    with open(output_file, "wb") as f:
        noisy_data.tofile(f)

    print(f"Noisy data written to '{output_file}'.")


# Local
# path = "../../queries/"
# convert_files = {
#    "turing_ANNS__head" : [100, np.int8]
#}

# Server
path = "."
convert_files = {
    # "turing_ANNS__head" : [100, np.int8]
    "turingANNs" : [100, np.int8],
    "text-to-image" : [200, np.int8],
    "spacev1B" : [100, np.int8]
}

if __name__ == "__main__":
    for noise_level in [0.1, 1, 2, 5, 10]:
        for input_file in  convert_files:
            length, file_type = convert_files[input_file]
            print (input_file, length, file_type)                        
            add_noise(input_file, length, file_type, noise_level)
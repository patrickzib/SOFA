import os
import struct
import numpy as np
import pandas as pd

import sys
sys.path.append("../")

def add_gaussian_noise(data, mean=0, std_dev=0.1):
    noise = np.random.normal(0, np.std(data)* std_dev, size=len(data))
    return data.astype(np.float64) + noise

def add_noise(input_file, dimensions, data_type, has_header, noise_level=0.1):
    # File paths
    output_file = input_file +"_noise_"+str(noise_level).replace(".","")
    input_file = path+"/"+input_file+".bin"
    output_file = path+"/generated/"+output_file+".bin"
        
    try:
        with open(input_file, "rb") as f:
            if has_header:
                header = f.read(8)
                if len(header) != 8:
                    raise ValueError(f"Input file '{input_file}' has no complete vector header.")
                _, header_dimensions = struct.unpack("<II", header)
                if header_dimensions != dimensions:
                    raise ValueError(
                        f"Input file '{input_file}' declares {header_dimensions} dimensions, "
                        f"expected {dimensions}."
                    )
            data = np.fromfile(f, dtype=data_type, count=1024 * 1024)
    except FileNotFoundError:
        print(f"Input file '{input_file}' not found.")
        return

    
    # Add Gaussian noise
    noisy_data = add_gaussian_noise(data, std_dev=noise_level)
    if np.issubdtype(data_type, np.integer):
        limits = np.iinfo(data_type)
        noisy_data = np.clip(np.rint(noisy_data), limits.min, limits.max).astype(data_type)
    else:
        noisy_data = noisy_data.astype(data_type)
    
    # Write noisy data to output file
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
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
    # Canonical Big ANN types: Turing-ANNS and Text-to-Image are float32;
    # SPACEV is signed int8.  Keep generated queries binary-compatible with
    # their source collection.
    "turingANNs": (100, np.float32, True),
    "text-to-image": (200, np.float32, True),
    "spacev1B": (100, np.int8, False),
}

if __name__ == "__main__":
    for noise_level in [0.1, 1, 2, 5, 10]:
        for input_file in  convert_files:
            dimensions, file_type, has_header = convert_files[input_file]
            print(input_file, dimensions, file_type, "header=" + str(has_header))
            add_noise(input_file, dimensions, file_type, has_header, noise_level)

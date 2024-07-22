import os
import io

path = "/vol/tmp/schaefpa/seismic"
os.environ["SEISBENCH_CACHE_ROOT"] = path
print("Path", os.environ["SEISBENCH_CACHE_ROOT"])

import seisbench as sb
import seisbench.data as sbd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

from tqdm import tqdm
from numba import njit, prange

print("Datasets dir", sb.cache_root)
print("Buffer Size", io.DEFAULT_BUFFER_SIZE)

@njit(cache=True, fastmath=True, parallel=True)
def parallel_znorm(waveform):
    n = waveform.shape[1]
    data = np.zeros((3, n // 256, 256), dtype=np.float32)
    # Write ZNE components
    for d in prange(data.shape[0]):
        for j in prange(data.shape[1]):
            data[d, j] = normalize(
                waveform[d, j * 256:(j + 1) * 256].astype(np.float32)
            )
    return data


@njit(cache=True, fastmath=True)
def normalize(vector):
    mean = np.mean(vector)
    stddev = np.std(vector)
    if stddev <= 1e-8:
        return vector - mean

    vector = (vector - mean) / stddev
    return vector


print("Cache root:", sb.cache_root)

datasets = [
    #"ETHZ",
    #"GEOFON",
    #"InstanceCounts",
    #"Iquique",
    #"ISC_EHB_DepthPhases",
    "LenDB",
    #"LFEStacksCascadiaBostock2015",
    #"LFEStacksMexicoFrank2014",
    #"LFEStacksSanAndreasShelly2017",
    # TODO "MLAAPDE",
    #"NEIC",
    #"OBS",
    #"OBST2024",
    #"PNW",
    #"Meier2019JGR",
    #"Ross2018GPD",
    #"Ross2018JGRFM",
    #"Ross2018JGRPick",
    #"STEAD",
    #"TXED",
]

"""


"""

def switch_dataset(name):
    if name == "ETHZ":
        return sbd.ETHZ(sampling_rate=100, component_order="ZNE")
    elif name == "GEOFON":
        return sbd.GEOFON(sampling_rate=100, component_order="ZNE")
    elif name == "InstanceCounts":
        return sbd.InstanceCounts(sampling_rate=100, component_order="ZNE")
    elif name == "InstanceNoise":
        return sbd.InstanceNoise(sampling_rate=100, component_order="ZNE")
    elif name == "InstanceCountsCombined":
        return sbd.InstanceCountsCombined(sampling_rate=100, component_order="ZNE")
    elif name == "Iquique":
        return sbd.Iquique(sampling_rate=100, component_order="ZNE")
    elif name == "ISC_EHB_DepthPhases":
        return sbd.ISC_EHB_DepthPhases(sampling_rate=100, component_order="ZNE")
    elif name == "LenDB":
        return sbd.LenDB(sampling_rate=100, component_order="ZNE")
    elif name == "LFEStacksCascadiaBostock2015":
        return sbd.LFEStacksCascadiaBostock2015(sampling_rate=100, component_order="ZNE")
    elif name == "LFEStacksMexicoFrank2014":
        return sbd.LFEStacksMexicoFrank2014(sampling_rate=100, component_order="ZNE")
    elif name == "LFEStacksSanAndreasShelly2017":
        return sbd.LFEStacksSanAndreasShelly2017(sampling_rate=100, component_order="ZNE")
    elif name == "MLAAPDE":
        return sbd.MLAAPDE(sampling_rate=100, component_order="ZNE", force=True)
    elif name == "NEIC":
        return sbd.NEIC(sampling_rate=100, component_order="ZNE")
    elif name == "OBS":
        return sbd.OBS(sampling_rate=100, component_order="ZNE")
    elif name == "OBST2024":
        return sbd.OBST2024(sampling_rate=100, component_order="ZNE")
    elif name == "PNW":
        return sbd.PNW(sampling_rate=100, component_order="ZNE")
    elif name == "PNWAccelerometers":
        return sbd.PNWAccelerometers(sampling_rate=100, component_order="ZNE")
    elif name == "PNWExotic":
        return sbd.PNWExotic(sampling_rate=100, component_order="ZNE")
    elif name == "PNWNoise":
        return sbd.PNWNoise(sampling_rate=100, component_order="ZNE")
    elif name == "Meier2019JGR":
        return sbd.Meier2019JGR(sampling_rate=100, component_order="ZNE")
    elif name == "Ross2018GPD":
        return sbd.Ross2018GPD(sampling_rate=100, component_order="ZNE")
    elif name == "Ross2018JGRFM":
        return sbd.Ross2018JGRFM(sampling_rate=100, component_order="ZNE")
    elif name == "Ross2018JGRPick":
        return sbd.Ross2018JGRPick(sampling_rate=100, component_order="ZNE")
    elif name == "STEAD":
        return sbd.STEAD(sampling_rate=100, component_order="ZNE")
    elif name == "TXED":
        return sbd.TXED(sampling_rate=100, component_order="ZNE")

def extact(dataset):
    data = switch_dataset(dataset)
    savePath = path + "/" + data.name + ".bin"
    savePathQueries = path + "/" + data.name + "_queries.bin"

    if True: # not os.path.isfile(savePath):
        # data.preload_waveforms(pbar=True)
        output_file_queries = open(savePathQueries, 'wb')

        num_samples_total = (data.metadata).shape[0]
        print("Processing", data.name, "with", num_samples_total, "samples")
        print("\tPath:", savePath, savePathQueries)

        if "trace_p_arrival_sample" in data.metadata.columns:
            pick = data.metadata["trace_p_arrival_sample"]
        elif "trace_P_arrival_sample" in data.metadata.columns:
            pick = data.metadata["trace_P_arrival_sample"]
        else:
            print("\tNo p-wave arrival times found - skipping", data.name)
            print(data.metadata.columns)
            exit()

        # First 200 waves as queries
        print("\tExtracting up to 200 queries")
        indices = np.argwhere(~np.isnan(pick)).squeeze()

        for d in range(3):
            for a in tqdm(np.arange(min(100, len(indices))), ncols=40, miniters=1):
                p_wave = int(pick[indices[a]])
                # waveform, _ = data.get_sample(indices[a], sampling_rate=100)
                waveform = data.get_waveforms(a, sampling_rate=100)

                # Extract 256 values around p-wave
                if (p_wave + 206) > waveform.shape[1]:
                    p_wave = waveform.shape[1] - 206
                elif p_wave < 50:
                    p_wave = 50

                znormalized = normalize(waveform[d, p_wave - 50:p_wave + 206])

                output_file_queries.write(znormalized.astype(np.float32))

        output_file_queries.flush()
        output_file_queries.close()
    else:
        print("Skipping", data.name, "as it already exists")


for dataset in datasets:
    extact(dataset)

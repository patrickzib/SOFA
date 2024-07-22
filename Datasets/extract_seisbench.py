import os
import io

path = "./datasets/"
# path = "/vol/tmp/schaefpa/seismic"
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
    if stddev <= 1e-4:
        return vector - mean

    vector = (vector - mean) / stddev
    return vector


print("Cache root:", sb.cache_root)

datasets = [
    #"ETHZ",
    #"GEOFON",
    #"InstanceCounts",
    "Iquique",
    "ISC_EHB_DepthPhases",
    "LenDB",
    "LFEStacksCascadiaBostock2015",
    "LFEStacksMexicoFrank2014",
    "LFEStacksSanAndreasShelly2017",
    "MLAAPDE",
    "NEIC",
    "OBS",
    "OBST2024",
    "PNW",
    "Meier2019JGR",
    "Ross2018GPD",
    "Ross2018JGRFM",
    "Ross2018JGRPick",
    "STEAD",
    "TXED",
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
    #elif name == "InstanceNoise":
    #    return sbd.InstanceNoise(sampling_rate=100, component_order="ZNE")
    #elif name == "InstanceCountsCombined":
    #    return sbd.InstanceCountsCombined(sampling_rate=100, component_order="ZNE")
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
        return sbd.MLAAPDE(sampling_rate=100, component_order="ZNE")
    elif name == "NEIC":
        return sbd.NEIC(sampling_rate=100, component_order="ZNE")
    elif name == "OBS":
        return sbd.OBS(sampling_rate=100, component_order="ZNE")
    elif name == "OBST2024":
        return sbd.OBST2024(sampling_rate=100, component_order="ZNE")
    elif name == "PNW":
        return sbd.PNW(sampling_rate=100, component_order="ZNE")
    #elif name == "PNWAccelerometers":
    #    return sbd.PNWAccelerometers(sampling_rate=100, component_order="ZNE")
    #elif name == "PNWExotic":
    #    return sbd.PNWExotic(sampling_rate=100, component_order="ZNE")
    #elif name == "PNWNoise":
    #    return sbd.PNWNoise(sampling_rate=100, component_order="ZNE")
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
    data = sbd.LFEStacksMexicoFrank2014(sampling_rate=100, cache="trace")
    # data = switch_dataset(dataset)
    savePath = path + "/" + data.name + ".bin"

    if not os.path.isfile(savePath):
        data.preload_waveforms(pbar=True)

        num_samples_total = (data.metadata).shape[0]
        print("Processing", data.name, "with", num_samples_total, "samples")

        if "trace_p_arrival_sample" in data.metadata.columns:
            pick = data.metadata["trace_p_arrival_sample"]
        elif "trace_P_arrival_sample" in data.metadata.columns:
            pick = data.metadata["trace_P_arrival_sample"]
        else:
            print("\tNo p-wave arrival times found - skipping", data.name)
            print(data.metadata.columns)
            exit()
        
        # Remaining data
        output_file = open(savePath, 'wb')
        ts_counter = 0
        print("\tExtracting data from", 100, "to", num_samples_total, "samples")
        for idx in tqdm(range(100, num_samples_total), ncols=40, miniters=1):
            # No more than 100 Mio TS
            if ts_counter > 100000000:
                print("Enough time series crawled")
                break

            waveform, _ = data.get_sample(idx, sampling_rate=100)
            # waveform = data.get_waveforms(idx, sampling_rate=100)

            # Standard seismometers will consist of three components, commonly vertical (Z),
            # north-south (N) and east-west (E). Depending on your application, you'll need
            # to arrange the components differently. SeisBench can do this automatically.

            znormed = parallel_znorm(waveform)
            output_file.write(znormed.tobytes())
            ts_counter += znormed.shape[1] * znormed.shape[0]

        print("\tTotal number of time series extracted:", ts_counter, "of length 256")
        print("-----------------")

        output_file.flush()
        output_file.close()
    else:
        print("Skipping", data.name, "as it already exists")


for dataset in datasets:
    extact(dataset)

This is the supporting website for the paper "Fast and Exact Similarity Search in less than a Blink of an Eye".


# To compile MESSI v2.0 (Autotools, from repo root)
```bash
./configure
make
```

# Scripts

See the provided scripts in the `scripts`-folder for examples to run SOFA with SFA summarization.

The SOFA command is 

```bash
FILE_PATH=/vol/tmp/schaefpa/messi_datasets/deep1b.bin
QUERIES_PATH=/vol/tmp/schaefpa/messi_datasets/$QUERY
TS_SIZE=96

COEFF_NUMBER=32
DATASET_SIZE=100000000
SAMPLE_SIZE=1000000
QUERY_SIZE=100

./MESSI --dataset --dataset $FILE_PATH --in-memory --timeseries-size $TS_SIZE  --function-type 4 
--dataset-size $DATASET_SIZE --flush-limit 300000 --read-block 20000 --sax-cardinality 8 
--queries $QUERIES_PATH --queries-size $QUERY_SIZE --queue-number $2 --sample-size $SAMPLE_SIZE 
--sample-type 3 --cpu-type $1 --is-norm --histogram-type 2 --leaf-size 20000 --min-leaf-size 20000 
--initial-lbl-size 20000 --coeff-number $COEFF_NUMBER  --SIMD
```

For help, please type:
```bash
./MESSI --help
```


# Datasets

Instruction for downloading the datasets is in the `datasets` folder. The size of the datasets is too large to provide a direct link.
Some datasets must be downloaded, others generated from seisbench.

## Table: Characteristics of the 17 datasets used

| Dataset Name | Series       | Series Length |
|--------------|--------------|---------------|
| **Astro** [soldi2014long] | 100,000,000   | 256           |
| **BigANN** [simhadri2022results] | 100,000,000   | 100           |
| **Deep1b** [babenko2016efficient] | 100,000,000   | 96            |
| **ETHZ** [woollam2022seisbench] | 4,999,932     | 256           |
| **Iquique** [woollam2019convolutional] | 578,853       | 256           |
| **ISC_EHB_DepthPhases** [munchmeyer2024learning] | 100,000,000   | 256           |
| **LenDB** [magrini2020local] | 37,345,260    | 256           |
| **Meier2019JGR** [woollam2022seisbench] | 6,361,998     | 256           |
| **NEIC** [yeck2021leveraging] | 93,473,541    | 256           |
| **OBS** [bornstein2024pickblue] | 15,508,794    | 256           |
| **OBST2024** [niksejel2024obstransformer] | 4,160,286     | 256           |
| **PNW** [ni2023curated] | 31,982,766    | 256           |
| **SALD** [url:SALD] | 100,000,000   | 128           |
| **SCEDC** [center2013southern] | 100,000,000   | 256           |
| **SIFT1b** [jegou2011searching] | 100,000,000   | 128           |
| **STEAD** [mousavi2019stanford] | 87,323,433    | 256           |
| **TXED** [chen2024txed] | 35,851,641    | 256           |

# Competitors

The competitors are stored within the `competitors` folder.

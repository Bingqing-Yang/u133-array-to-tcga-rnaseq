# Download raw data
## TCGA data
Current PanCanAtlas contains only the RNA-seq data, thus we will be using other sources such as GEO to download the microarray data.

### RNA-seq
Ensure `gdc-client` is installed, otherwise, please refer to [gdc-client github](https://github.com/NCI-GDC/gdc-client) and [GDC Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) for instruction on installation.

Then we will download the raw TCGA dataset from PanCanAtlas.

script: `./pancan/rna-seq/get.sh`

input: `./pancan/rna-seq/manifest/PanCan-General_Open_GDC-Manifest_2.txt` # TODO: we will need to upload the file

outputs: # TODO: add in precise files in each directory
- `./pancan/rna-seq/gdc/00a32f7a-c85f-4f86-850d-be53973cbc4d/`
- `./pancan/rna-seq/gdc/3586c0da-64d0-4b74-a449-5ff4d9136611/`
- `./pancan/rna-seq/gdc/0f4f5701-7b61-41ae-bda9-2805d1ca9781/`
- `./pancan/rna-seq/gdc/4f277128-f793-4354-a13d-30cc7fe9f6b5/`
- `./pancan/rna-seq/gdc/0fc78496-818b-4896-bd83-52db1f533c5c/`
- `./pancan/rna-seq/gdc/55d9bf6f-0712-4315-b588-e6f8e295018e/`
- `./pancan/rna-seq/gdc/1a7d7be8-675d-4e60-a105-19d4121bdebf/`
- `./pancan/rna-seq/gdc/7d4c0344-f018-4ab0-949a-09815f483480/`
- `./pancan/rna-seq/gdc/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81/`
- `./pancan/rna-seq/gdc/99b0c493-9e94-4d99-af9f-151e46bab989/`
- `./pancan/rna-seq/gdc/1c6174d9-8ffb-466e-b5ee-07b204c15cf8/`
- `./pancan/rna-seq/gdc/d82e2c44-89eb-43d9-b6d3-712732bf6a53/`
- `./pancan/rna-seq/gdc/1c8cfe5f-e52d-41ba-94da-f15ea1337efc/`
- `./pancan/rna-seq/gdc/fcbb373e-28d4-4818-92f3-601ede3da5e1/`

### Microarray data

# Preprocess data
## TCGA data
We will extract samples that have both U133 and RNA-seq data.

script: `./pancan/preprocess.R`

# filter genes
Rscript filter.R

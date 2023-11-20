# Download raw data
Ensure `gdc-client` is installed, otherwise, please refer to [gdc-client](https://github.com/NCI-GDC/gdc-client) for instruction on installation.

Then we will download the raw TCGA dataset from PanCanAtlas.

script: `get.sh`
output: 
```
# download raw data
./get.sh
# preprocess data
Rscript preprocess.R
# filter genes
Rscript filter.R
```

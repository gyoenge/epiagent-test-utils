## EpiAgent Finetune for In-silico Cancer, Dataset Preparing 

### Source

- Dataset - Epigenetic regulation during cancer transitions across 11 tumour types.
- GEO:GSE240822 
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240822

### Steps 

- Download GEO raw dataset 

    ```
    $ chmod +x download_geodataset.sh
    $ ./download_geodataset.sh
    ```

- Unzip .tar.gz files 

    ```
    $ chmod +x gunzip_geodataset.sh
    $ ./gunzip_geodataset.sh
    ```

- Delete original .tar.gz files

- Unzip .tsv.gz files 

    ```
    $ chmod +x unzip_tsv.sh
    $ ./unzip_tsv.sh
    ```

- Delete original .tsv.gz files

- Convert .tsv to .bed format

    ```
    $ chmod +x tsv_to_bed.sh
    $ ./tsv_to_bed.sh
    ```

- Make sub directories 

    ```
    $ mkdir temp
    $ mkdir raw_h5ad
    $ mkdir processed_h5ad
    ```

- Prepare Meta Data

    ```
    $ chmod +x prepare_metadata.sh
    $ ./prepare_metadata.sh
    ```

- Preprocess data into .h5ad(annData) with metadata 
    
    requires python environment : including epiagent==0.0.3, pandas

    ```
    $ python preprocessing.py 
    ```

- Integrate dataset for has-meta sub-data-cells

    ```
    $ python hasmeta_integration.py 
    ```

- Downsample integrated dataset 

    ```
    $ python downsampling.py 
    ```

- (final) Convert into cell-sentence, constructing the complete train-test EpiAgent dataset 

    ```
    $ python hasmeta_integration.py 
    ```

### Check data

- With `dataset_check.ipynb`

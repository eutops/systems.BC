# systems.BC

This repository contains code required to reproduce analyses and visualisations in *Systems epigenetic approach towards non-invasive breast cancer detection* (Herzog et al., 2024).

## Data

Datasets have been deposited online and are described in the manuscript, including a full list of accession numbers.

## System requirements

This code is based on the R programming language and has been tested under the R version 4.3.1 (2023-06-16), mac OS and Windows.

This code has no non-standard hardware requirements.

### Dependencies

Code required to run the analyses and generate figures has the following dependencies:

* circlize (0.4.15)
* clusterProfiler (4.8.3)
* devtools (2.4.5)
* DMRcate (2.14.1)
* dplyr (1.1.4)
* EpiDISH (2.16.0)
* epitools (0.5.10.1)
* fs (1.6.3)
* GEOquery (2.68.0)
* glmnet (4.1.8)
* ggforce (0.4.1)
* ggplot2 (3.4.4)
* ggpubr (0.6.0)
* ggsci (3.0.0)
* ggtext (0.1.2)
* grDevices (4.3.1)
* gt (0.10.0)
* gtsummary (1.7.2)
* here (1.0.1)
* IlluminaHumanMethylationEPICanno.ilm10b4.hg19 (0.6.0)
* karyoploteR (1.26.0)
* missMethyl (1.34.0)
* org.Hs.eg.db (3.17.0)
* plotly (4.10.3)
* ReactomePA (1.44.0)
* rms (6.7.1)
* stringr (1.5.1)
* TCGAbiolinks (2.28.4)
* tibble (3.2.1)
* tidyverse (2.0.0)
* tidyr (1.3.1)
* patchwork (1.1.3)
* pROC (1.18.5)

## Installation

Download or fork the repository and open the .Rproj file to open the project. Code to reproduce the analyses can then be found organised by subfolders. 

### Typical install time

The typical install time should be less than 1 minute.

## Instruction for use 

Raw data are deposited in EGA and cannot be provided openly due to data protection laws. 

However, summary data for plotting are provided in the dataframes-plotting folder. These can be used to reproduce the figures.

### Expected output

Depending on the script run, expected output should be raw figures underlying the manuscripts. The run time will depend on the figure but should typically be less than a few seconds.

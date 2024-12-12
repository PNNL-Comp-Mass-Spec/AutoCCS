# AutoCCS #

Automated Collision Cross Section calculation software for ion mobility-mass spectrometry


## Main features ##

AutoCCS supports the following platforms and methods:
- Platforms
  - Drift tube-based ion mobility spectrometry coupled with mass spectrometry (DTIMS-MS) instrument [(Stow,S.M. et al., 2017)](https://doi.org/10.1021/acs.analchem.7b01729), [PMID 28763190](https://pubmed.ncbi.nlm.nih.gov/28763190/)
  - Traveling wave based-IMS MS (TWIMS-MS) such as Structures for Lossless Ion Manipulations (SLIM)-based IMS-MS [(Wojcik,R. et al., 2019)](https://doi.org/10.1021/acs.analchem.9b02808), [PMID 31450886](https://pubmed.ncbi.nlm.nih.gov/31450886/)
- Methods
  - Stepped field method [(Stow,S.M. et al., 2017)](https://doi.org/10.1021/acs.analchem.7b01729), [PMID 28763190](https://pubmed.ncbi.nlm.nih.gov/28763190/)
  - Single field method [(Kurulugama,R.T. et al., 2015)](https://doi.org/10.1039/c5an00991j), [PMID 26178817](https://pubmed.ncbi.nlm.nih.gov/26178817/)
  - Traveling wave-based method 
- Calibration functions for single field and traveling wave-based methods
  - Linear function
  - Polynomial functions (e.g. quadratic or cubic functions)
  - Linearized power function [(Ruotolo,B.T. et al., 2008)](https://doi.org/10.1038/nprot.2008.78), [PMID 18600219](https://pubmed.ncbi.nlm.nih.gov/18600219/)
  

## How to _install_ AutoCCS ##

#### Use `conda` environment (Recommended) ####

Please install [conda](https://docs.anaconda.com/anaconda/install/) and create an environment from an environment.yml file. More details about managing the conda environment can be found on the [Managing environments](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file).
```bash
conda env create -f environment.yml
```
And then activate the new conda environment `autoccs`.
```bash
conda activate autoccs
```

#### Use `pip` ####

Install Python 3.7 (or newer) [[link]](https://www.python.org/downloads/) and use `pip` as follows to install dependencies. 
```bash
pip install -r requirements.txt
```

## Tutorial with demo data ##
In this tutorial, we demonstrated the CCS determination using AutoCCS for the Agilent tune-mix samples in three different platforms: stepped-field DTIMS-MS, single-field RapidFire-DTIMS-MS, SLIM-based IMS.

Demo dataset is publicly available at MassIVE (Accession: [MSV000085979](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000085979))

<i>For readability purposes, the input parameters are split over multiple lines. When using the command line, all parameters must be included as a single line.</i>

#### Stepped-Field DTIMS-MS ####

<!-- 
```bash
python -u autoCCS.py --target_list_file data/SteppedField-DTIMS/TargetList.csv --config_file data/SteppedField-DTIMS/autoCCS_config.xml --framemeta_files "data/SteppedField-DTIMS/ImsMetadata/*.txt" --feature_files "data/SteppedField-DTIMS/Features_csv/*.csv" --output_dir data/SteppedField-DTIMS/Results/ --threshold_n_fields 5 --mode multi &> "data/SteppedField-DTIMS/LogFiles/multi.log"
``` 
-->

```bash
python -u autoCCS.py
  --target_list_file data/SteppedField-DTIMS/TargetList.csv
  --config_file data/SteppedField-DTIMS/autoCCS_config.xml
  --framemeta_files "data/SteppedField-DTIMS/ImsMetadata/*.txt"
  --feature_files "data/SteppedField-DTIMS/Features_csv/*.csv"
  --output_dir data/SteppedField-DTIMS/Results/
  --threshold_n_fields 5
  --mode multi &> "data/SteppedField-DTIMS/LogFiles/multi.log"
```

#### Single-Field RapidFire-DTIMS-MS ####

<!--
```bash
python autoCCS.py --config_file data/SingleField-RapidFire-DTIMS/autoCCS_config.xml --framemeta_files "data/SingleField-RapidFire-DTIMS/ImsMetadata/*.txt" --feature_files "data/SingleField-RapidFire-DTIMS/Features_csv/*.csv" --calibrant_file data/SingleField-RapidFire-DTIMS/TuneMix-CCS.txt --output_dir data/SingleField-RapidFire-DTIMS/Results/ --tunemix_sample_type AgilentTuneMix --sample_meta data/SingleField-RapidFire-DTIMS/Datasets.csv --colname_for_sample_type SampleType --colname_for_filename RawFileName --colname_for_ionization IonPolarity --degree 1 --single_mode batch --mode single &> data/SingleField-RapidFire-DTIMS/LogFiles/single.log
```
-->

```bash
python autoCCS.py
  --config_file data/SingleField-RapidFire-DTIMS/autoCCS_config.xml
  --framemeta_files "data/SingleField-RapidFire-DTIMS/ImsMetadata/*.txt"
  --feature_files "data/SingleField-RapidFire-DTIMS/Features_csv/*.csv"
  --calibrant_file data/SingleField-RapidFire-DTIMS/TuneMix-CCS.txt
  --output_dir data/SingleField-RapidFire-DTIMS/Results/
  --tunemix_sample_type AgilentTuneMix
  --sample_meta data/SingleField-RapidFire-DTIMS/Datasets.csv
  --colname_for_sample_type SampleType
  --colname_for_filename RawFileName
  --colname_for_ionization IonPolarity
  --degree 1
  --single_mode batch
  --mode single &> data/SingleField-RapidFire-DTIMS/LogFiles/single.log
```

#### SLIM-based IMS-MS ####

<!--
```bash
python -u autoCCS.py --config_file data/SLIM-IMS/autoCCS_config.xml --feature_files "data/SLIM-IMS/Features_csv/*.csv" --output_dir data/SLIM-IMS/Results/ --mode single --calibrant_file data/SLIM-IMS/TuneMix-CCS_POS.txt --sample_meta data/SLIM-IMS/Datasets.csv --tunemix_sample_type AgilentTuneMix --colname_for_sample_type SampleType --colname_for_filename RawFileName --colname_for_ionization IonPolarity --single_mode batch --degree 2 --calib_method poly --ppm 150 &> data/SLIM-IMS/LogFiles/slim.log
```
-->

```bash
python -u autoCCS.py
  --config_file data/SLIM-IMS/autoCCS_config.xml
  --feature_files "data/SLIM-IMS/Features_csv/*.csv"
  --output_dir data/SLIM-IMS/Results/
  --mode single
  --calibrant_file data/SLIM-IMS/TuneMix-CCS_POS.txt
  --sample_meta data/SLIM-IMS/Datasets.csv
  --tunemix_sample_type AgilentTuneMix
  --colname_for_sample_type SampleType
  --colname_for_filename RawFileName
  --colname_for_ionization IonPolarity
  --single_mode batch
  --degree 2
  --calib_method poly
  --ppm 150 &> data/SLIM-IMS/LogFiles/slim.log
```
Users are allowed to apply high-order polynomial functions: quadratic (`--degree 2`), cubic (`--degree 3`), quartic (`--degree 4`), and so on.
```bash
  --degree 3 # for cubic
```
Also, it allows users to apply non-linear regression based on the linearized power function.
```bash
  --calib_method power
```

## Citation ##

JY Lee, A Bilbao, CR Conant, KJ Bloodsworth, DJ Orton, M Zhou, ... & TO Metz (2021). AutoCCS: Automated collision cross section calculation software for ion mobility spectrometry-mass spectrometry. <i>Bioinformatics</i>. DOI: [10.1093/bioinformatics/btab429](https://doi.org/10.1093/bioinformatics/btab429), [PMID 28763190](https://pubmed.ncbi.nlm.nih.gov/28763190/)

## Contacts ##

Written by Joon-Yong Lee for the Department of Energy (PNNL, Richland, WA)\
Copyright 2020, Battelle Memorial Institute. All Rights Reserved.\
E-mail: proteomics@pnnl.gov\
Website: https://github.com/PNNL-Comp-Mass-Spec/AutoCCS or https://www.pnnl.gov/integrative-omics


## License ##

AutoCCS is licensed under the BSD 2-Clause License; [License](license.txt)

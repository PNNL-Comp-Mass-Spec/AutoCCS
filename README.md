# AutoCCS #

Automated Collision Cross Section calculation software for ion mobility-mass spectrometry


## Main features ##
AutoCCS supports the following platforms and methods:
- Platforms
  - Drift tube-based ion mobility spectrometry coupled with mass spectrometry (DTIMS-MS) instrument [(Stow,S.M. et al., 2017)](https://doi.org/10.1021/acs.analchem.7b01729)
  - Traveling wave based-IMS MS (TWIMS-MS) such as Structures for Lossless Ion Manipulations (SLIM)-based IMS-MS [(Wojcik,R. et al., 2019)](https://doi.org/10.1021/acs.analchem.9b02808)
- Methods
  - Stepped field method [(Stow,S.M. et al., 2017)](https://doi.org/10.1021/acs.analchem.7b01729)
  - Single field method [(Kurulugama,R.T. et al., 2015)](https://doi.org/10.1039/C5AN00991J)
  - Traveling wave-based method 
- Calibration functions for single field and traveling wave-based methods
  - Linear function
  - Polynomial functions (e.g. quadratic or cubic functions)
  - Linearized power function [(Ruotolo,B.T. et al., 2008)](https://doi.org/10.1038/nprot.2008.78)
  
  
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
Install python(>3.7) [[link]](https://www.python.org/downloads/) and use `pip` as follows to install dependencies. 
```bash
pip install -r requirements.txt
```

## Tutorial with demo data ##
In this tutorial, we demonstrated the CCS determination using AutoCCS for the Agilent tune-mix samples in three different platforms: stepped-field DTIMS-MS, single-field RapidFire-DTIMS-MS, SLIM-based IMS.

<i>For sake of readability, the input parameters are split over multiple lines. When using the command line, however, all parameters should be included as a single line.</i>

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
python autoCCS.py --config_file data/SingleField-RapidFire-DTIMS/autoCCS_config.xml --framemeta_files "data/SingleField-RapidFire-DTIMS/ImsMetadata/*.txt" --feature_files "data/SingleField-RapidFire-DTIMS/Features_csv/*.csv" --calibrant_file data/SingleField-RapidFire-DTIMS/TuneMix-CCS.txt --output_dir data/SingleField-RapidFire-DTIMS/Results/ --tunemix_sample_type AgilentTuneMix --sample_meta data/SingleField-RapidFire-DTIMS/Datasets.csv --colname_for_sample_type SampleType --colname_for_filename RawFileName --colname_for_ionization IonPolarity --degree 1 --single_mode batch --mode single &> data/SingleField-RapidFire-DTIMS/LogFiles/multi.log
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
  --mode single &> data/SingleField-RapidFire-DTIMS/LogFiles/multi.log
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

## Contacts ##
Written by Joon-Yong Lee for the Department of Energy (PNNL, Richland, WA)\
Copyright 2020, Battelle Memorial Institute. All Rights Reserved.\
E-mail: joonyong.lee@pnnl.gov or proteomics@pnnl.gov\
Website: https://omics.pnl.gov/ or https://panomics.pnnl.gov/


## License ##
AutoCCS is licensed under the BSD 2-Clause License; [License](license.txt)





# Multi-Field
python multiCCS.py --target_list_file ../data/IO-Files_SteppedField/TargetList.csv --config_file ../data/autoCCS_config.xml --framemeta_files "../data/IO-Files_SteppedField/IV_ImsMetadata/*.txt" --feature_files "../data/IO-Files_SteppedField/IV_Features_csv/*.csv" --output_dir ../data/IO-Files_SteppedField/IV_Results/ &> "../data/IO-Files_SteppedField/LogFiles/multi.log"

python autoCCS.py --target_list_file ../data/IO-Files_SteppedField/TargetList.csv --config_file ../data/autoCCS_config.xml --framemeta_files "../data/IO-Files_SteppedField/IV_ImsMetadata/*.txt" --feature_files "../data/IO-Files_SteppedField/IV_Features_csv/*.csv" --output_dir ../data/IO-Files_SteppedField/IV_Results/ --mode multi &> "../data/IO-Files_SteppedField/LogFiles/multi.log"

# Single-Field
python singleCCS.py --config_file ../data/IO-Files_SingleField/autoCCS_config.xml --framemeta_files "../data/IO-Files_SingleField/IV_ImsMetadata/*.txt" --feature_files "../data/IO-Files_SingleField/IV_Features_csv/*.csv" --calibrant_file ../data/IO-Files_SingleField/TuneMix-CCS.txt --output_dir ../data/IO-Files_SingleField/IV_Results/ --tunemix_sample_type AgilentTuneMix --sample_meta ../data/IO-Files_SingleField/Datasets.csv --colname_for_sample_type SampleType --colname_for_filename RawFileName --colname_for_ionization IonPolarity --single_mode batch

python autoCCS.py --config_file ../data/IO-Files_SingleField/autoCCS_config.xml --framemeta_files "../data/IO-Files_SingleField/IV_ImsMetadata/*.txt" --feature_files "../data/IO-Files_SingleField/IV_Features_csv/*.csv" --calibrant_file ../data/IO-Files_SingleField/TuneMix-CCS.txt --output_dir ../data/IO-Files_SingleField/IV_Results/ --tunemix_sample_type AgilentTuneMix --sample_meta ../data/IO-Files_SingleField/Datasets.csv --colname_for_sample_type SampleType --colname_for_filename RawFileName --colname_for_ionization IonPolarity --single_mode batch --mode single

# SLIM
python autoCCS.py --config_file ../data/IO-Files_SLIM/autoCCS_config.xml --framemeta_files "../data/IO-Files_SLIM/IV_ImsMetadata/*.txt" --feature_files "../data/IO-Files_SLIM/IV_Features_csv/*.csv" --calibrant_file ../data/IO-Files_SLIM/TuneMix-CCS_POS.txt --output_dir ../data/IO-Files_SLIM/IV_Results/ --tunemix_sample_type AgilentTuneMix --sample_meta ../data/IO-Files_SLIM/Datasets.csv --colname_for_sample_type SampleType --colname_for_filename RawFileName --colname_for_ionization IonPolarity --single_mode batch --degree 3 --ppm 50 --mode single



<!-- CASP -->
python autoCCS.py --config_file /Users/leej324/Documents/papers/2020_AutoCCS/data/CASP/autoCCS_config.xml --feature_files "/Users/leej324/Downloads/CASP_2/All_Features_csv/*.csv" --sample_meta /Users/leej324/Documents/papers/2020_AutoCCS/data/CASP/20200313_CASP_Datasets.csv --tunemix_sample_type AgilentTuneMix --colname_for_sample_type Plate --colname_for_filename RawFileNameNew --colname_for_ionization IonPolarity --calibrant_file /Users/leej324/Documents/papers/2020_AutoCCS/data/CASP/TuneMix-CCS.txt --output_dir /Users/leej324/Downloads/CASP_2/Result_wo_meta/ --single_mode batch --mode single


python autoCCS.py --config_file /Users/leej324/Documents/papers/2020_AutoCCS/data/CASP/autoCCS_config.xml --feature_files "/Users/leej324/Downloads/CASP_2/All_Features_csv/AgilentTuneMix-POS-*.csv" --sample_meta /Users/leej324/Documents/papers/2020_AutoCCS/data/CASP/20200313_CASP_Datasets.csv --tunemix_sample_type AgilentTuneMix --colname_for_sample_type Plate --colname_for_filename RawFileNameNew --colname_for_ionization IonPolarity --framemeta_files "/Users/leej324/Downloads/CASP_2/All_FrameMetadata/AgilentTuneMix-POS-*.txt"  --calibrant_file /Users/leej324/Documents/papers/2020_AutoCCS/data/CASP/TuneMix-POS-CCS.txt --output_dir /Users/leej324/Downloads/CASP_2/Result_pos_calibration_with_meta/ --single_mode fit --mode single



python autoCCS.py --config_file /Users/leej324/Documents/papers/2020_AutoCCS/data/CASP/autoCCS_config.xml --feature_files "/Users/leej324/Downloads/CASP_2/All_Features_csv/AgilentTuneMix-POS-*.csv" --sample_meta /Users/leej324/Documents/papers/2020_AutoCCS/data/CASP/20200313_CASP_Datasets.csv --tunemix_sample_type AgilentTuneMix --colname_for_sample_type Plate --colname_for_filename RawFileNameNew --colname_for_ionization IonPolarity --framemeta_files "/Users/leej324/Downloads/CASP_2/All_FrameMetadata/AgilentTuneMix-POS-*.txt"  --calibrant_file /Users/leej324/Documents/papers/2020_AutoCCS/data/CASP/TuneMix-POS-CCS.txt --output_dir /Users/leej324/Downloads/CASP_2/Result_pos_calibration_with_meta/ --single_mode fit --mode single
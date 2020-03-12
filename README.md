# Multi-Field
python multiCCS.py --target_list_file ../data/IO-Files_SteppedField/TargetList.csv --config_file ../data/autoCCS_config.xml --framemeta_files "../data/IO-Files_SteppedField/IV_ImsMetadata/*.txt" --feature_files "../data/IO-Files_SteppedField/IV_Features_csv/*.csv" --output_dir ../data/IO-Files_SteppedField/IV_Results/ &> "../data/IO-Files_SteppedField/LogFiles/multi.log"

python autoCCS.py --target_list_file ../data/IO-Files_SteppedField/TargetList.csv --config_file ../data/autoCCS_config.xml --framemeta_files "../data/IO-Files_SteppedField/IV_ImsMetadata/*.txt" --feature_files "../data/IO-Files_SteppedField/IV_Features_csv/*.csv" --output_dir ../data/IO-Files_SteppedField/IV_Results/ --mode multi &> "../data/IO-Files_SteppedField/LogFiles/multi.log"

# Single-Field
python singleCCS.py --config_file ../data/IO-Files_SingleField/autoCCS_config.xml --framemeta_files "../data/IO-Files_SingleField/IV_ImsMetadata/*.txt" --feature_files "../data/IO-Files_SingleField/IV_Features_csv/*.csv" --calibrant_file ../data/IO-Files_SingleField/TuneMix-CCS.txt --output_dir ../data/IO-Files_SingleField/IV_Results/ --tunemix_sample_type AgilentTuneMix --sample_meta ../data/IO-Files_SingleField/Datasets.csv --colname_for_sample_type SampleType --colname_for_filename RawFileName --colname_for_ionization IonPolarity --single_mode batch

python autoCCS.py --config_file ../data/IO-Files_SingleField/autoCCS_config.xml --framemeta_files "../data/IO-Files_SingleField/IV_ImsMetadata/*.txt" --feature_files "../data/IO-Files_SingleField/IV_Features_csv/*.csv" --calibrant_file ../data/IO-Files_SingleField/TuneMix-CCS.txt --output_dir ../data/IO-Files_SingleField/IV_Results/ --tunemix_sample_type AgilentTuneMix --sample_meta ../data/IO-Files_SingleField/Datasets.csv --colname_for_sample_type SampleType --colname_for_filename RawFileName --colname_for_ionization IonPolarity --single_mode batch --mode single

# SLIM
python autoCCS.py --config_file ../data/IO-Files_SLIM/autoCCS_config.xml --framemeta_files "../data/IO-Files_SLIM/IV_ImsMetadata/*.txt" --feature_files "../data/IO-Files_SLIM/IV_Features_csv/*.csv" --calibrant_file ../data/IO-Files_SLIM/TuneMix-CCS_POS.txt --output_dir ../data/IO-Files_SLIM/IV_Results/ --tunemix_sample_type AgilentTuneMix --sample_meta ../data/IO-Files_SLIM/Datasets.csv --colname_for_sample_type SampleType --colname_for_filename RawFileName --colname_for_ionization IonPolarity --single_mode batch --degree 3 --ppm 50 --mode single





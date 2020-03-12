import argparse
from multiCCS import multi
from singleCCS import single
from utils import get_config

##########################################################################
# ArgumentParser
##########################################################################
parser = argparse.ArgumentParser()

parser.add_argument(
    '--mode', type=str, required=True, choices=['single','multi'], default="single",
    help='single: single-field mode, multi: multi-field mode')

parser.add_argument(
    '--feature_files', required=True, type=str,
    help='feature files to determine CCS values')

parser.add_argument(
    '--framemeta_files', required=True, type=str,
    help='frame meta info file for samples')

parser.add_argument(
    '--output_dir', type=str, required=True, default='./',
    help='a directory to store output files')

parser.add_argument(
    '--config_file', type=str, required=True, default='config.xml',
    help='Configuration file')

################################################
# args for sigle-field
################################################
parser.add_argument(
    '--calibrant_file', type=str,
    help='calibrant file contains the target m/z and CCS')

parser.add_argument(
    '--tune_mix_regexp', type=str, 
    help='a regular expression for tune mix sample files')

parser.add_argument(
    '--tune_mix_frame_regexp', type=str, default='',
    help='frame feature files to calibrate CCS values')

parser.add_argument(
    '--calibration_curves', type=str,
    help='calibration curves obtained from Tune Mix samples')

parser.add_argument(
    '--sample_meta', type=str,
    help='meta info file for samples')

parser.add_argument(
    '--standard_mass', type=float,
    help='internal standard mass')

parser.add_argument(
    '--colname_for_sample_type', type=str, default="SampleType",
    help='colname for sample type in sample_meta')

parser.add_argument(
    '--colname_for_filename', type=str, default="RawFileNameNew",
    help='colname for filename in sample_meta')

parser.add_argument(
    '--colname_for_ionization', type=str, default="IonPolarity",
    help='colname for ionization in sample_meta')

parser.add_argument(
    '--tunemix_sample_type', type=str, default="AgilentTuneMix",
    help='sample type for tunemix in sample_meta')

parser.add_argument(
    '--skip_calibrated_colname', type=str, default=None,
    help='a column name for the calibrated CCS values \
        when you want to skip some files already having calibrated CCS values')

parser.add_argument(
    '--single_mode', type=str, choices=['fit','ccs','batch'], default="fit",
    help='fit: curve fitting, ccs: CCS calibration')

parser.add_argument(
    '--ppm', type=float,
    help='ppm of mz tolerance. It will override mz_tolerance in config_params')

parser.add_argument(
    '--degree', type=int, default=1,
    help='Degree of the fitting polynomial for CCS calibration curves')

################################################
# args for multi-field
################################################
parser.add_argument(
    '--target_list_file', type=str, default='TargetList.txt',
    help='Target list file (Tab-delimited text format)')

# parser.add_argument(
#     '--data_folder', type=str, default='./',
#     help='Data folder containing all the cef and meta data files')

parser.add_argument(
    '--output', type=str, default='ccs_table.tsv',
    help='Output file to save a output table')

parser.add_argument(
    '--r2_threshold', type=float, default=0.99,
    help='threshold value for r2')

parser.add_argument(
    '--num_isotopes_threshold', type=int, default=1,
    help='threshold value for num_isotopes')

parser.add_argument(
    '--intensity_rank_threshold', type=int, default=3,
    help='threshold value for peak intensity rank in m/z window')

parser.add_argument(
    '--maxint', action='store_true',
    help='select max intensive peaks for ccs computation')

parser.add_argument(
    '--format', type=str, choices=['cef','mzmine'], default='mzmine',
    help='file format for the features, e.g., cef or mzmine')

FLAGS = {}

##########################################################################


def main(FLAGS):
    print("#"*50)
    print("# Parameters:", FLAGS)
    print("#"*50)

    # read a set of configuration parameters
    config_params = get_config(FLAGS.config_file)
    print(config_params)

    if FLAGS.mode=="single":
        print("Single-Field Mode")
        single(FLAGS, config_params)
    elif FLAGS.mode=="multi":
        print("Multi-Field Mode")
        multi(FLAGS, config_params)

if __name__ == '__main__':
    FLAGS = parser.parse_args()
    main(FLAGS)

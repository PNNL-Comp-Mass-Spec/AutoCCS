import pandas as pd
import numpy as np

import xml.etree.ElementTree
import glob
import os

from sys import platform as sys_pf

if sys_pf == 'darwin':
    import matplotlib

    matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import seaborn as sns

from utils import *
from CCSCalibrator import *

import argparse

##########################################################################
# ArgumentParser
##########################################################################
parser = argparse.ArgumentParser()

parser.add_argument(
    '--calibrant_file', type=str,
    help='calibrant file contains the target m/z and CCS')

# @deprecated
parser.add_argument(
    '--tune_mix_regexp', type=str,
    help='a regular expression for tune mix sample files')

# @deprecated
parser.add_argument(
    '--tune_mix_framemeta_regexp', type=str, default='',
    help='frame feature files for tune mix samples to calibrate CCS values')

parser.add_argument(
    '--feature_files', type=str,
    help='feature files to calibrate CCS values')

parser.add_argument(
    '--output_dir', type=str, default='./',
    help='a directory to store output files')

parser.add_argument(
    '--config_file', type=str, default='config.xml',
    help='Configuration file')

parser.add_argument(
    '--single_mode', type=str, required=True, choices=['fit', 'ccs', 'batch'], default="batch",
    help='fit: curve fitting, ccs: CCS calibration, batch: fitting and calibration')

parser.add_argument(
    '--calibration_curves', type=str,
    help='calibration curves obtained from Tune Mix samples')

parser.add_argument(
    '--sample_meta', type=str,
    help='meta info file for samples')

parser.add_argument(
    '--framemeta_files', type=str, default=None,
    help='frame meta info file for samples')

parser.add_argument(
    '--standard_mass', type=float,
    help='internal standard mass')

parser.add_argument(
    '--ppm', type=float,
    help='ppm of mz tolerance. It will override mz_tolerance in config_params')

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
    '--degree', type=int, default=1,
    help='Degree of the fitting polynomial for CCS calibration curves')

FLAGS = {}

# Monoisotopic mass https://fiehnlab.ucdavis.edu/staff/kind/Metabolomics/MS-Adduct-Calculator/
mH = 1.007276
m2H = 2 * mH
mNa = 22.989218
mH2O = 18.010564686
mK = 38.963706861
N2 = 28.006148
He = 4.002603


##########################################################################

# def get_calibrants(calibrant_file, drift_gas_mass=28.006148):
#     '''read calibrants and add reduced mass and ccs values based on the drift gas
#     '''
#     calibrants = pd.read_csv(calibrant_file, sep='\t')
#     # print(calibrants)
#     calibrants['reduced_mass'] = (calibrants['m/z'] * calibrants.z * drift_gas_mass)/ (drift_gas_mass + calibrants['m/z'] * calibrants.z)
#     calibrants['reduced_ccs'] = calibrants.CCS / (calibrants.z / (calibrants.reduced_mass**0.5))
#     return calibrants

# def etree_to_dict(t):
#     d = {t.tag : map(etree_to_dict, t.iter())}
#     d.update(('@' + k, v) for k, v in t.attrib.items())
#     d['text'] = t.text
#     d['tag'] = t.tag
#     return d

def get_frame_meta(metafiles, sep="_FrameMetadata", offset=0):
    '''aggregate multiple metadata files and extract the field information for each frame
        if offset<=0: average
        else: df.iloc[offset-1] offset-th value
    Return
        a pandas dataframe having a field information for each frame
    '''
    metafiles.sort()
    meta_data = {}
    for mf in metafiles:
        print(mf)
        f = mf.split("/")[-1].split(sep)[0]
        mfdf = pd.read_csv(mf, sep='\t')

        # if pressures and temperatures are missing
        if (mfdf.ImsPressure.isna().sum() > 0) | mfdf.ImsTemperature.isna().sum() > 0:
            return None

        if offset <= 0:
            meta_data[f] = [mfdf.ImsPressure.mean(), mfdf.ImsTemperature.mean()]
        else:
            meta_data[f] = [mfdf.ImsPressure.iloc[offset - 1], mfdf.ImsTemperature.iloc[offset - 1]]
    return meta_data


# def get_features_from_cef(cef, max_normalize=True):
#     e = xml.etree.ElementTree.parse(cef).getroot()
#     json = etree_to_dict(e.findall('CompoundList')[0])
#     idx = 0
#     mppid = 0
#     rst = []
#     mspeaks = []
#     in_ms_peaks = False
#     for j in json['CompoundList']:
#         if 'Compound' in j:
#             mppid = j['@mppid']

#         if 'Location' in j:
#             mz = j['@m']
#             rt = j['@rt']
#             intensity = j['@y']
#             dt = j['@dt']
#             rst.append({'mppid':mppid, 'mz':float(mz), 'rt':float(rt), 'intensity':float(intensity), 'dt':float(dt)})

#         if 'MSPeaks' in j:
#             for k in j['MSPeaks']:
#                 if ('p' in k):
#                     mspeaks.append({'mppid':mppid, 'mass':float(k['@x']), 'intensity_org':float(k['@y']), 'z':int(k['@z']), 's':k['@s']})

#     df = pd.DataFrame(rst)
#     mspeaks = pd.DataFrame(mspeaks)
#     if df.shape[0] > 0:
#         df['intensity_z'] = (df.intensity - df.intensity.mean())/df.intensity.std()
#         if max_normalize:
#             df.intensity /= df.intensity.max()
#         mspeaks = mspeaks.drop_duplicates(subset='mppid', keep='first')
#         df = pd.merge(mspeaks, df, left_on="mppid", right_on="mppid", how='inner')
#     return df, mspeaks

# def get_features_from_mzmine_csv(csv, max_normalize=True):
#     df = pd.read_csv(csv)
#     col_id = [c for c in df.columns if c.endswith('row ID')][0]
#     col_area = [c for c in df.columns if c.endswith('Peak area')][0]
#     col_height = [c for c in df.columns if c.endswith('Peak height')][0]
#     col_mz = [c for c in df.columns if c.endswith('Peak m/z')][0]
#     col_dt = [c for c in df.columns if c.endswith('Peak RT')][0]
#     if 'calibrated_ccs' in df.columns:
#         cols = [col_id, col_mz, col_dt, col_area, col_height, 'calibrated_ccs']
#         colnames = ['mppid', 'mass', 'dt', 'intensity', 'height', 'calibrated_ccs']
#     else:
#         cols = [col_id, col_mz, col_dt, col_area, col_height]
#         colnames = ['mppid', 'mass', 'dt', 'intensity', 'height']

#     df = df[cols].copy()
#     df.columns = colnames
#     df['intensity_org'] = df.intensity
#     if df.shape[0] > 0:
#         df['intensity_z'] = (df.intensity - df.intensity.mean())/df.intensity.std()
#         if max_normalize:
#             df.intensity /= df.intensity.max()
#     return df, None

# def is_in_tolerance(x, mass, ppm):
#     delta = mass * ppm * 1.0e-6
#     return (x >= mass - delta) & (x <= mass + delta)

def datetime(r):
    return pd.to_datetime(r)


def nearest_tunemix(row, **kwargs):
    '''find a file name of the tunemix sample with the nearest timestamp
    '''
    tunemix = kwargs['tunemix']
    if row.Ionization != row.Ionization: return None
    temp = tunemix[tunemix.Ionization == row.Ionization]
    minidx = ((temp.time - row.time).abs()).idxmin()
    if minidx != minidx: return None
    return temp.loc[minidx].newFile


# def plot_features(features, selected, ax, f):
#     ax.scatter(features.mass, features.intensity_z, s=3)
#     ax.scatter(selected.mass, selected.intensity_z, s=2, c='red')
#     ax.set_title(f)
#     ax.set_ylabel('intensity(z-score)')

# def plot_multi_features(file_list, calibrants, ppm=20, dformat='cef'):
#     plt.close('all')
#     fig, axis = plt.subplots(len(file_list), sharex=True, sharey=True, figsize=(8,8))

#     if dformat=='cef':
#         get_features = get_features_from_cef
#     elif dformat=='mzmine':
#         get_features = get_features_from_mzmine_csv

#     for i, f in enumerate(file_list):
#         ax = axis[i]
#         features, peaks = get_features(f)
#         selected = select_features(features, calibrants, f, ppm=ppm)
#         plot_features(features, selected, ax, f)
#     plt.xlabel('m/z')
#     plt.show()
#     return selected

# def select_features(features, calibrants, file, ppm=20):
#     selected = pd.DataFrame()
#     for row in calibrants.iterrows():
#         mz = row[1]['m/z']
#         reduced_ccs = row[1]['reduced_ccs']
#         ccs = row[1]['CCS']
#         _ff = features[is_in_tolerance(features.mass, mz, ppm)]
#         _ff = _ff.assign(calibrants_mz=mz, ccs=ccs, reduced_ccs=reduced_ccs, replicate=file)
#         if selected.shape[0] == 0:
#             selected = _ff
#         else:
#             selected = selected.append(_ff)
#     return selected

# def plot_ccs_dt(selected, calibrants, ax, file, color):
#     features = selected.sort_values(by='intensity').drop_duplicates(subset=['calibrants_mz'], keep='last')
#     features = features.sort_values(by='dt')
#     ax.plot(features.dt, features.reduced_ccs, color=color, label="")
#     ax.scatter(features.dt, features.reduced_ccs, color=color, label=os.path.basename(file))
#     ax.set_ylabel('Reduced CCS')

#     # ax.set_title(cef)
#     return features

# def plot_calibrate_ccs(file_list, calibrants, ppm=20, dformat='cef'):
#     '''
#         dformat: file format. cef - Agilent MassProfiler, mzmine - MZmine CSV
#     '''
#     plt.close('all')
#     fig, axis = plt.subplots(len(file_list), sharex=True, sharey=True, figsize=(8,8))
#     colors=sns.color_palette("Paired")
#     num_colors = len(colors)
#     if dformat=='cef':
#         get_features = get_features_from_cef
#     elif dformat=='mzmine':
#         get_features = get_features_from_mzmine_csv
#     for i, f in enumerate(file_list):
#         ax = axis[i]
#         features, _ = get_features(f)
#         _selected = select_features(features, calibrants, f, ppm=ppm)
#         _selected = plot_ccs_dt(_selected, calibrants, ax, f, colors[i%num_colors])
#         if i == 0: selected = _selected
#         else: selected = selected.append(_selected)
#     plt.xlabel('Drift time')

#     plt.show()
#     return selected

# def plot_calibrate_ccs_in_one(file_list, calibrants, ppm=20, dformat='cef', ofile=None):
#     '''
#         dformat: file format. cef - Agilent MassProfiler, mzmine - MZmine CSV
#     '''
#     plt.close('all')
#     fig, ax = plt.subplots(1, sharex=True, sharey=True, figsize=(8,8))
#     colors=sns.color_palette("tab20c")
#     num_colors = len(colors)

#     if dformat=='cef':
#         get_features = get_features_from_cef
#     elif dformat=='mzmine':
#         get_features = get_features_from_mzmine_csv

#     for i, f in enumerate(file_list):
#         features, _ = get_features(f)
#         _selected = select_features(features, calibrants, f, ppm=ppm)
#         _selected = plot_ccs_dt(_selected, calibrants, ax, f, colors[i%num_colors])
#         if i == 0: selected = _selected
#         else: selected = selected.append(_selected)
#     plt.xlabel('Drift time')
#     ax.legend()
#     # ax.legend([os.path.basename(f) for f in file_list])
#     if ofile == None:
#         plt.show();
#     else:
#         plt.tight_layout()
#         plt.savefig(ofile, dpi=600)
#     return selected

# def calibrate_ccs(file, calibrate_function, drift_gas_mass=28.006148, dformat='cef', z=1):
#     if dformat=='cef':
#         get_features = get_features_from_cef
#     elif dformat=='mzmine':
#         get_features = get_features_from_mzmine_csv
#     features, _ = get_features(file)
#     gamma = np.sqrt(features.mass/(features.mass + drift_gas_mass))/z
#     beta = calibrate_functions[1]
#     tfix = calibrate_functions[0]
#     features['calibrated_ccs'] = (features.dt - tfix)/(beta*gamma)
#     return features

# def curve_fit(selected, file_list, drift_gas_mass=28.006148, deg=1, z=1):
#     # curve fitting for each rep
#     curve_fit_fn = []
#     for i, f in enumerate(file_list):
#         data = selected[selected.replicate==f]
#         print(i, f, data.shape)
#         if data.shape[0] == 0: continue;

#         gamma = np.sqrt(data.mass/(data.mass + drift_gas_mass))/z

#         x = gamma*data.ccs
#         y = data.dt

#         p = np.poly1d(np.polyfit(x, y, deg))

#         # r-squared
#         # fit values, and mean
#         yhat = p(x)                      # or [p(z) for z in x]
#         ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
#         ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
#         sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
#         r = ssreg / sstot

#         print(p, r)
#         curve_fit_fn.append({"file":f, "tunemix":f.split('.mzML')[0].split('/')[-1], "tfix":p[0], "beta":p[1]})
#     return curve_fit_fn


def calibrate_ccs(file, calibrate_fn, drift_gas_mass=28.006148, dformat='cef', z=1):
    if dformat == 'cef':
        get_features = get_features_from_cef
    elif dformat == 'mzmine':
        get_features = get_features_from_mzmine_csv
    features, _ = get_features(file)
    gamma = np.sqrt(features.mass / (features.mass + drift_gas_mass)) / z
    beta = calibrate_fn.beta
    tfix = calibrate_fn.tfix
    features['calibrated_ccs'] = (features.dt - tfix) / (beta * gamma)
    return features


def calibrate_ccs_with_framemeta(file, calibrate_fn, frame_info=None, drift_gas_mass=28.006148, dformat='cef', z=1):
    if dformat == 'cef':
        get_features = get_features_from_cef
    elif dformat == 'mzmine':
        get_features = get_features_from_mzmine_csv
    features, _ = get_features(file)

    gamma = np.sqrt(features.mass / (features.mass + drift_gas_mass)) / z
    poly_fn = np.poly1d(calibrate_fn.poly)
    # beta = calibrate_fn.beta
    # tfix = calibrate_fn.tfix

    # print(frame_info)
    if frame_info:
        p_torr, temp = frame_info
        t_k = temp + 273.15
        gamma = gamma * p_torr / np.sqrt(t_k)

    features['calibrated_ccs'] = poly_fn(features.dt) / gamma
    # features['calibrated_ccs'] = (features.dt - tfix)/(beta*gamma)

    return features


def compute_ccs_mzmine(files, meta_info, calibration_curves, frame_meta=None, out_dir="./",
                       drift_gas_mass=28.006148, sep=".mzML", skip_calibrated_colname=None):
    '''
        files: mzmine output csv files
        meta_info: it has acquired time, csv file name, the nearest tunemix
    '''
    sample2tunemix = pd.Series([i.split(sep)[0] for i in meta_info.tunemix.values],
                               index=[i.split(sep)[0] for i in meta_info.newFile]).to_dict()
    # print(sample2tunemix)
    for i, f in enumerate(files):
        filename = os.path.basename(f).split(sep)[0]
        # print('filename', filename, sep)
        ####################################################################################
        ### bypass the empty files (<100 bytes)
        ####################################################################################
        if os.stat(f).st_size < 100:
            print("[ERR] file size is too small (<100 bytes) so that it loooks has no feature.")
            print("See this file: {}".format(filename))
            continue
        ####################################################################################
        ####################################################################################
        ### bypass files already having "calibrated_ccs" (passed by skip_calibrated_colname)
        ####################################################################################
        if skip_calibrated_colname:
            with open(f) as handle:
                line = handle.readline()
                if "," + skip_calibrated_colname in line:
                    # print('skip_calibrated_colname:', filename)
                    continue
        ####################################################################################
        tunemix = sample2tunemix[filename]
        if frame_meta:
            frame_info = frame_meta[filename]
        else:
            frame_info = None
        if calibration_curves[calibration_curves.tunemix == tunemix].shape[0] > 0:
            # print(filename, tunemix, calibration_curves[calibration_curves.tunemix==tunemix])
            calibration_function = calibration_curves[calibration_curves.tunemix == tunemix].iloc[0]
            features = calibrate_ccs_with_framemeta(f, calibration_function, frame_info,
                                                    drift_gas_mass=drift_gas_mass, dformat='mzmine')

            # add CCS to original
            org_df = pd.read_csv(f).dropna(how='all', axis=1)
            org_df['calibrated_ccs'] = features['calibrated_ccs']
            org_df.to_csv(out_dir + "/" + os.path.basename(f), index=False)

        if (i + 1) % 500 == 0: print("[{0}/{1}] {2} - {3}".format(i + 1, len(files), filename, tunemix))


def get_target_ions_all(files, meta_info, mass=395.149, ppm=40, ion_mode='POS'):
    '''
        all ions within mass range
    '''
    internal_standard_features = []
    for i, f in enumerate(files):
        filename = os.path.basename(f).split('.mzML')[0].split('/')[-1] + ".d"
        ####################################################################################
        ### bypass the empty files (<100 bytes)
        ####################################################################################
        if os.stat(f).st_size < 100:
            print("[ERR] file size is too small (<100 bytes) so that it loooks has no feature.")
            print("See this file: {}".format(filename))
            continue

        selected_info = meta_info[meta_info.newFile == filename]
        if selected_info.shape[0] != 1:
            print('[ERR] please check out this file:', f)
            continue
        if selected_info.Ionization.tolist()[0] != ion_mode:  # different ionization
            continue

        # group = selected_info.Ionization.tolist()[0] +'-'+ selected_info.Cartridge.tolist()[0]
        group = selected_info.Group.tolist()[0]
        corrected_time = selected_info.AcquiredTime.tolist()[0]

        features, _ = get_features_from_mzmine_csv(f)
        _ff = features[is_in_tolerance(features.mass, mass=mass, ppm=ppm)].copy()
        if _ff.shape[0] == 0:
            continue
        else:
            _ff['csv'] = filename
            _ff['Group'] = group
            _ff['corrected_time'] = corrected_time
            internal_standard_features.append(_ff)

        if (i + 1) % 500 == 0: print(i + 1, len(internal_standard_features))
    print(i + 1, len(internal_standard_features))
    if len(internal_standard_features) > 0:
        df = pd.concat(internal_standard_features)
        df['time'] = df.corrected_time.apply(datetime)
        return df
    else:
        return pd.DataFrame()


def get_target_ions(files, meta_info, mass=395.149, ppm=40):
    '''
        most intense
    '''
    internal_standard_features = []
    for i, f in enumerate(files):
        filename = os.path.basename(f).split('.mzML')[0].split('/')[-1] + ".d"
        group = meta_info[meta_info.newFile == filename].Group.tolist()[0]
        corrected_time = meta_info[meta_info.newFile == filename].corrected_time.tolist()[0]

        features, _ = get_features_from_mzmine_csv(f)
        _ff = features[is_in_tolerance(features.mass, mass=mass, ppm=ppm)].copy()
        if _ff.shape[0] > 1:
            # print(f, _ff.shape, _ff.intensity_org.tolist(), _ff.dt.tolist())
            _ff = _ff.sort_values(by='intensity_org').iloc[-1, :]
            _ff = _ff.to_dict()
        elif _ff.shape[0] == 0:
            continue
        else:
            _ff = _ff.to_dict(orient='records')[0]
            _ff['csv'] = filename
            _ff['Group'] = group
            _ff['corrected_time'] = corrected_time
            internal_standard_features.append(_ff)

        if (i + 1) % 500 == 0: print(i + 1, len(internal_standard_features))
    print(i + 1, len(internal_standard_features))
    df = pd.DataFrame(internal_standard_features)
    df['time'] = df.corrected_time.apply(datetime)
    return df


def plot_internal_standard(internal_standard_df, y='dt', adduct="[M+H]", mode='POS'):
    df = internal_standard_df.copy()

    ## filtering by intensity
    df = df.sort_values(by='intensity').drop_duplicates(subset=['csv'], keep='last')
    print(df)

    plt.close('all')
    # sns.distplot(df[df.Group==mode+'-NA'].dt, hist=False, label=mode+'-NA')
    sns.distplot(df[df.Group == mode + '-GC'].dt, hist=False, label=mode + '-GC')
    sns.distplot(df[df.Group == mode + '-C18'].dt, hist=False, label=mode + '-C18')
    # sns.distplot(df[df.csv.str.endswith('QC')].dt, hist=False, label='QC')
    plt.savefig('{0}-{1}-internal_standard_dt.pdf'.format(mode, adduct))

    plt.close('all')
    # sns.distplot(df[df.Group==mode+'-NA'].calibrated_ccs, hist=False, label=mode+'-NA')
    sns.distplot(df[df.Group == mode + '-GC'].calibrated_ccs, hist=False, label=mode + '-GC')
    sns.distplot(df[df.Group == mode + '-C18'].calibrated_ccs, hist=False, label=mode + '-C18')

    tdf = df[((df.Group == mode + '-GC') | (df.Group == mode + '-C18'))]
    sd, avg = tdf.calibrated_ccs.std(ddof=1), tdf.calibrated_ccs.mean()
    ccs_info = "mean:{0:.1f}, std:{1:.1f}, %RSD:{2:.1f}%".format(avg, sd, 100 * sd / avg)
    print(ccs_info)

    plt.text(180, 0.2, ccs_info)
    plt.savefig('{0}-{1}-internal_standard_ccs.pdf'.format(mode, adduct))


def find_features_with_calibrator(file_list, calibrator, sep="", ppm=20, dformat='mzmine'):
    if dformat == 'cef':
        get_features = get_features_from_cef
    elif dformat == 'mzmine':
        get_features = get_features_from_mzmine_csv

    # read feature files
    selected = []
    for i, f in enumerate(file_list):
        features, _ = get_features(f)

        ionization = "pos"
        if "NEG" in f.upper(): ionization = "neg"

        tmp = calibrator.find_calibrant_features(features, ppm=ppm, ionization=ionization)
        tmp['filename'] = f.split("/")[-1].split(sep)[0]
        selected.append(tmp)
    selected = pd.concat(selected, ignore_index=True)
    return selected


def curve_fit_with_calibrator(selected, file_list, calibrator,
                              frame_meta=None, sep="", drift_gas_mass=28.006148, deg=1, z=1, fout=None):
    # curve fitting for each rep

    plt.close('all')
    fig, axis = plt.subplots(1, sharex=True, sharey=True, figsize=(10, 8))
    colors = sns.color_palette("Paired")
    num_colors = len(colors)

    curve_fit_fn = []
    file_list.sort()

    print('sep', sep)

    for i, f in enumerate(file_list):
        print("#" * 100)
        fname = f.split('/')[-1].split(sep)[0]

        data = selected[selected.filename == fname]
        if data.shape[0] == 0:
            continue
        print(i, f, data.shape)

        gamma = np.sqrt(data.mass / (data.mass + drift_gas_mass)) / z

        if frame_meta:
            print('frame_meta', frame_meta)
            p_torr, temp = frame_meta[fname]
            print('Pressure (Torr):{}, Temperature (C):{}'.format(p_torr, temp))
            t_k = temp + 273.15
            y = gamma * data.ccs * p_torr / np.sqrt(t_k)
        else:
            y = gamma * data.ccs

        x = data.dt
        p, r = calibrator.fit(x, y, deg=deg)

        print(p, ", r=", r)

        curve_dict = {"file": f, "tunemix": fname, "poly": list(p)}
        curve_fit_fn.append(curve_dict)

        # plot points and lines
        pad = 0.1 * (np.max(x) - np.min(x))
        xp = np.linspace(np.min(x) - pad, np.max(x) + pad, 200)
        axis.scatter(y, x, color=colors[i % num_colors], label=fname)
        axis.plot(p(xp), xp, lw=0.1, color=colors[i % num_colors])
        # axis.scatter(y, x, color=colors[i % num_colors])
        # axis.plot(p(xp), xp, color=colors[i % num_colors], label=fname)

    plt.ylabel('Arrival time (ms)', fontsize=15)
    if frame_meta:
        plt.xlabel('Reduced CCS With Pressure and Temperature', fontsize=15)
    else:
        plt.xlabel('Reduced CCS', fontsize=15)
    plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", frameon=False)

    plt.tight_layout()
    if fout: plt.savefig(fout, dpi=300)
    # print(curve_fit_fn)
    return curve_fit_fn


def fit_calib_curves(FLAGS, config_params):
    drift_gas_mass = config_params['neutral_mass']
    ppm = config_params['mz_tolerance']
    if FLAGS.ppm: ppm = FLAGS.ppm

    sep = config_params['suffix_raw'].split("{")[0]
    prefix_tunemix = config_params['prefix_tunemix']

    # target calibrants to get a calibration curve
    calibrator = CCSCalibrator(FLAGS.calibrant_file, drift_gas_mass)
    # calibrants = get_calibrants(FLAGS.calibrant_file, drift_gas_mass)

    # select features to compute a calibration curve. TODO: to handle multiple replicates
    tunemix_files = [fname for fname in glob.glob(FLAGS.feature_files) if "/" + prefix_tunemix in fname]

    assert len(tunemix_files) > 0, \
        "Feature files for tune-mix samples are not found. Please check your 'prefix_tunemix' in your config file."

    frame_meta = None
    if FLAGS.framemeta_files is None:
        tunemix_framemeta_files = []
    else:
        tunemix_framemeta_files = [fname for fname in glob.glob(FLAGS.framemeta_files) if "/" + prefix_tunemix in fname]

    if len(tunemix_framemeta_files) > 0:
        frame_meta = get_frame_meta(tunemix_framemeta_files,
                                    sep=config_params['suffix_meta'], offset=config_params['frame_offset'])
    if frame_meta is None:
        print('[INFO] No frame meta data file is given.', frame_meta)

    # find the features for calibrants
    selected = find_features_with_calibrator(tunemix_files, calibrator,
                                             sep=sep, ppm=ppm, dformat='mzmine')
    print('[INFO] {} features are selected for determining calibration curves. (m/z tolerance: {}ppm)'.format(
        selected.shape[0], ppm))

    calibrate_functions = curve_fit_with_calibrator(selected, tunemix_files, calibrator,
                                                    frame_meta=frame_meta, sep=sep, drift_gas_mass=drift_gas_mass,
                                                    deg=FLAGS.degree, z=1,
                                                    fout=FLAGS.output_dir + "/calibration_output.pdf")

    if len(calibrate_functions) > 0:
        rst_df = pd.DataFrame(calibrate_functions)[['file', 'tunemix', 'poly']]
        rst_df.to_csv(FLAGS.output_dir + "/calibrate_functions.csv")
        return rst_df
    else:
        print("[ERROR] cannot find any good calibration curves")
        raise Exception


def perform_CCS_computation(FLAGS, config_params, _calib_curves=None):
    drift_gas_mass = config_params['neutral_mass']
    sep = config_params['suffix_raw'].split("{")[0]

    if FLAGS.sample_meta:
        print("#" * 80)
        print("# Collect sample metadata from", FLAGS.sample_meta)
        print("#" * 80)
        sample_meta = pd.read_csv(FLAGS.sample_meta)
        if FLAGS.colname_for_ionization not in sample_meta.columns:
            print(
                "[ERROR] cannot found a column name: {} in {}".format(FLAGS.colname_for_ionization, FLAGS.sample_meta))
        if FLAGS.colname_for_filename not in sample_meta.columns:
            print("[ERROR] cannot found a column name: {} in {}".format(FLAGS.colname_for_filename, FLAGS.sample_meta))
        sample_meta['Ionization'] = sample_meta[FLAGS.colname_for_ionization]
        sample_meta['newFile'] = sample_meta[FLAGS.colname_for_filename]

        if "tunemix" not in sample_meta.columns:
            sample_meta['time'] = sample_meta.AcquiredTime.apply(datetime)
            tunemix = sample_meta[sample_meta[FLAGS.colname_for_sample_type] == FLAGS.tunemix_sample_type].copy()
            sample_meta['tunemix'] = sample_meta.apply(nearest_tunemix, axis=1, tunemix=tunemix)
    else:
        print("[ERROR] --sample_meta is required to perform ccs computation")
        raise Exception

    print("#" * 80)
    if _calib_curves is None:
        print("# Collect calibration curve parameters from", FLAGS.calibration_curves)
        from ast import literal_eval
        calibration_curves = pd.read_csv(FLAGS.calibration_curves)
        calibration_curves[['poly']] = calibration_curves[['poly']].applymap(literal_eval)
    else:
        print("# Collect calibration curve parameters")
        calibration_curves = _calib_curves
    if calibration_curves.shape[0] > 0:
        print("\tOK (size: {})".format(calibration_curves.shape))
    print("#" * 80)

    # collect frame meta information (e.g., pressure and temperature)
    print("#" * 80)
    print("# Collect frame meta information from", FLAGS.framemeta_files)
    if FLAGS.framemeta_files is None:
        frame_meta = None
    else:
        frame_meta = get_frame_meta(glob.glob(FLAGS.framemeta_files), sep=config_params['suffix_meta'],
                                    offset=config_params['frame_offset'])
    if frame_meta is None: print('[INFO] No frame meta data file is given.', frame_meta)
    print("#" * 80)

    print("#" * 80)
    print("# Computing CCS values")
    print("#" * 80)
    compute_ccs_mzmine(glob.glob(FLAGS.feature_files), sample_meta, calibration_curves, frame_meta,
                       out_dir=FLAGS.output_dir,
                       drift_gas_mass=drift_gas_mass, sep=sep, skip_calibrated_colname=FLAGS.skip_calibrated_colname)


def assert_params_enough(FLAGS, config_params):
    '''check if all required parameters are given
    '''
    return True


def single(FLAGS, config_params):
    assert_params_enough(FLAGS, config_params)
    os.makedirs(FLAGS.output_dir, exist_ok=True)

    if FLAGS.single_mode == 'fit':
        calibrate_functions = fit_calib_curves(FLAGS, config_params)
        print(calibrate_functions)

    elif FLAGS.single_mode == 'ccs':
        perform_CCS_computation(FLAGS, config_params)

    elif FLAGS.single_mode == 'batch':
        calibrate_functions = fit_calib_curves(FLAGS, config_params)
        perform_CCS_computation(FLAGS, config_params, calibrate_functions)

    elif FLAGS.single_mode == 'standard':
        ''' 
        python singleCCS.py --sample_meta 20190717_CASP_Datasets.csv --feature_files "All_Features_csv/*.csv" --colname_for_filename RawFileNameNew --colname_for_ionization IonPolarity --single_mode standard --standard_mass 394.1416 --ppm 40
        Rotenone: 394.1416 (CASP)
        '''
        meta_info = pd.read_csv(FLAGS.sample_meta)
        meta_info['time'] = meta_info.AcquiredTime.apply(datetime)
        meta_info['Ionization'] = meta_info[FLAGS.colname_for_ionization]
        meta_info['newFile'] = meta_info[FLAGS.colname_for_filename]

        feature_files = glob.glob(FLAGS.feature_files)

        for ion_mode in ['POS', 'NEG']:
            if ion_mode == 'POS':
                adducts = {
                    '[M+dot]': FLAGS.standard_mass,
                    '[M+H]': FLAGS.standard_mass + mH,
                    '[M+Na]': FLAGS.standard_mass + mNa
                }
            elif ion_mode == 'NEG':
                adducts = {
                    '[M-H]': FLAGS.standard_mass - mH,
                    '[M-H2O-H]': FLAGS.standard_mass - 19.01839,
                    '[M-Na]': FLAGS.standard_mass - mNa,
                    '[M+Na-2H]': FLAGS.standard_mass + 20.974666
                }
            for a in adducts:
                adduct_mass = adducts[a]
                internal_standard_ions = get_target_ions_all(feature_files, meta_info, mass=adduct_mass, ppm=ppm,
                                                             ion_mode=ion_mode)
                if internal_standard_ions.shape[0] > 0:
                    plot_internal_standard(internal_standard_ions, y='calibrated_ccs', adduct=a, mode=ion_mode)


if __name__ == '__main__':
    FLAGS = parser.parse_args()

    # read a set of configuration parameters
    config_params = get_config(FLAGS.config_file)
    print(config_params)

    single(FLAGS, config_params)

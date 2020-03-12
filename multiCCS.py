'''compute CCS in multi-step experiments
'''
import pandas as pd
import numpy as np

import xml.etree.ElementTree
import time
import glob
import os

from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import seaborn as sns

from utils import *

from compute_ccs import SteppedFieldCCS
from shortest_for_ccs import get_possible_ccs_values

import argparse

##########################################################################
# ArgumentParser
##########################################################################
parser = argparse.ArgumentParser()

parser.add_argument(
    '--target_list_file', type=str, default='TargetList.txt',
    help='Target list file (Tab-delimited text format)')

parser.add_argument(
    '--config_file', type=str, default='config.xml',
    help='Configuration file')

# parser.add_argument(
#     '--data_folder', type=str, default='./',
#     help='Data folder containing all the cef and meta data files')
parser.add_argument(
    '--feature_files', type=str,
    help='feature files to calibrate CCS values')

parser.add_argument(
    '--framemeta_files', type=str,
    help='frame meta info file for samples')

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

parser.add_argument(
    '--output_dir', type=str, default='./',
    help='a directory to store output files')

FLAGS = {}

##########################################################################

def get_metadata(mfile, offset, ax=None, label=None):
    '''read metadata file and extract the field information for each frame
    TODO: offset method (choose one frame by offset) or average in a range
    Return
        a pandas dataframe having a field information for each frame
    '''
    try:
        metadata = pd.read_csv(mfile, sep='\t')
        _list = list(metadata.drop_duplicates(subset='FrameMethodId').FrameId+offset-1)
        filtered = metadata[metadata.FrameId.isin(_list)]
        ##################################################
        if ax is not None:
            ax[0].plot(metadata.FrameId, metadata.ImsTemperature, label=label)
            ax[0].scatter(filtered.FrameId, filtered.ImsTemperature, label=None)
            ax[0].set_ylabel('Temperature (C)')
            ax[1].plot(metadata.FrameId, metadata.ImsPressure)
            ax[1].scatter(filtered.FrameId, filtered.ImsPressure)
            ax[1].set_ylabel('Pressure (torr)')
            ax[2].plot(metadata.FrameId, metadata.ImsField)
            ax[2].scatter(filtered.FrameId, filtered.ImsField)
            ax[2].set_ylabel('E (V/cm)')
            ax[2].set_xlabel('Frame ID')
        ##################################################
        return filtered
    except Exception as e:
        return None


def get_target_info(target_list_file):
    '''read the target_list_file
        target_list_file: file path for a config file
    Return
        a pandas dataframe containing target information
    '''
    return pd.read_csv(target_list_file, sep='\t').fillna(method='ffill')


def get_adducts(exact_mass, adducts):
    '''get the adducts mass
        exact_mass: exact mass of the target
        adducts: configuration for adducts in config_file
    Return
        adducts2mass: a dict containing information of positive and negative adducts
    '''
    adducts2mass = {'pos':{}, 'neg':{}}
    for adduct in adducts:
        charges = adduct['charges'].replace(' ','').split(',')
        for c in charges:
            charge = int(c)
            name = '[M'+c+adduct['name']+']' if abs(charge)>1 else '[M'+c[0]+adduct['name']+']'
            mass = (exact_mass + charge * adduct['mass'])/abs(charge)
            if charge > 0:
                adducts2mass['pos'][name] = mass
            elif charge < 0:
                adducts2mass['neg'][name] = mass
    return adducts2mass


def get_features(file, max_normalize=True, fformat='cef'):
    if fformat=='cef': return get_features_from_cef(file, max_normalize)
    elif fformat=='mzmine': return get_features_from_mzmine_csv(file, max_normalize)
    else: print('File format: {0}. This tool doesn\'t support this file format.'.format(fformat))
    return None, None

def get_adducts_colors():
    return {'[M+.]':'m',
            '[M+H]':'b',
            '[M+2H]':'c',
            '[M+Na]':'r',
            '[M+K]':'g',
            '[M-H]':'y'}
        
def is_in_tolerance(x, mass, ppm):
    delta = mass * ppm * 1.0e-6
    #print(mass, delta, mass-delta, mass+delta)
    return (x >= mass - delta) & (x <= mass + delta)

def mass_error(x, mass):
    return abs(x - mass) / mass * 1e6

def find_features_maxint(features, metadata, ion_mz, ppm):
    df = features[is_in_tolerance(features.mass, ion_mz, ppm)]
    if df.shape[0] == 0: return df
    
    #  if 'frame' column in metadata, delete it
    if 'frame' in metadata.columns: del metadata['frame']

    df = df.sort_values(by='intensity_z').drop_duplicates(subset='frame', keep='last')
    df = df.merge(metadata, left_on='frame', right_on='FrameMethodId', how='inner')
    df = df.sort_values(by='frame')
    return df

def find_features(features, metadata, ion_mz, ppm,
                  threshold_num_isotopes=2,
                  threshold_intensity_rank=3):
    if 'num_isotopes' in features.columns:
        df = features[is_in_tolerance(features.mass, ion_mz, ppm) & (features.num_isotopes>=threshold_num_isotopes)]
    else:
        df = features[is_in_tolerance(features.mass, ion_mz, ppm)]
    if df.shape[0] == 0: return df
    
    # filter out small peaks by ranking threshold
    rankings = df.groupby('frame')['intensity_org'].rank(ascending=False)
    df = df[rankings<=threshold_intensity_rank]

    # for f in frames_too_many_features:
    #     filter_by_intensity_rank(df, f, threshold_intensity_rank)

    #  if 'frame' column in metadata, delete it
    if 'frame' in metadata.columns: del metadata['frame']

    # df = df.sort_values(by='intensity_z').drop_duplicates(subset='frame', keep='last')
    df = df.merge(metadata, left_on='frame', right_on='FrameMethodId', how='inner')
    # df = df.sort_values(by='frame')
    # df.to_csv("test_{0:.5f}.txt".format(ion_mz),sep="\t")
    return df

def filter_by_intensity_rank(df, frame, threshold_intensity_rank=3):
    temp = df[df.frame == frame]
    # print(df)
    # print(frame, temp.intensity_org)
    np.argsort(temp.intensity_org)

def ccs_filter(ccs_list):
    # remove the redundant regression lines which share the same start nodes(features)
    first_peaks = []
    last_peaks = []
    for ccs in ccs_list:
        first_peaks.append(int(ccs.mppid[0]))
        last_peaks.append(int(ccs.mppid[-1]))
    
    ufirst_peaks = list(np.unique(first_peaks))
    ulast_peaks = list(np.unique(last_peaks))
    if len(ufirst_peaks) < len(ccs_list):
        print("len(ufirst_peaks) < len(ccs_list)", len(ufirst_peaks),len(ccs_list))
        _ccs_list = []
        for u in ufirst_peaks:
            idx_list = np.where(first_peaks == u)[0]
            if idx_list.shape[0] > 1:
                best_r2 = 0
                best_ccs_u = None
                for ii in idx_list:
                    if (best_r2 < ccs_list[ii].r2):
                        best_ccs_u = ccs_list[ii]
                        best_r2 = ccs_list[ii].r2
                if best_ccs_u != None:
                    _ccs_list.append(best_ccs_u)
            else:
                _ccs_list.append(ccs_list[idx_list[0]])
        return _ccs_list

    elif len(ulast_peaks) < len(ccs_list):
        print("len(ulast_peaks) < len(ccs_list)", len(ulast_peaks),len(ccs_list))
        print("ulast_peaks", ulast_peaks)
        print("last_peaks", last_peaks)
        _ccs_list = []
        for u in ulast_peaks:
            idx_list = np.where(last_peaks == u)[0]
            print('idx_list',u, idx_list)
            if idx_list.shape[0] > 1:
                best_r2 = 0
                best_ccs_u = None
                for ii in idx_list:
                    if (best_r2 < ccs_list[ii].r2):
                        best_ccs_u = ccs_list[ii]
                        best_r2 = ccs_list[ii].r2
                if best_ccs_u != None:
                    _ccs_list.append(best_ccs_u)
            else:
                _ccs_list.append(ccs_list[idx_list[0]])
        return _ccs_list
    else:
        return ccs_list
        
    # find the ccs values of earlist molecules
    pass


def files_not_enough(fname, config_params, fformat='cef'):
    # meta_file = (fname + '{0}.txt').format(config_params['suffix_meta'])
    # if not os.path.isfile(meta_file):
    #     print("[ERROR] a metadata file doesn't exist:", meta_file)
    #     return True
    for step in range(config_params['num_fields']):
        if fformat=='cef': ffile = (fname + '{0}.cef').format(config_params['suffix_raw'].format(step+1))
        else: ffile = (fname + '{0}.csv').format(config_params['suffix_raw'].format(step+1))
        if not os.path.isfile(ffile):
            print("[ERROR] a feature file doesn't exist:", ffile)
            return True
    return False

def get_ccs(FLAGS, comp_id, target_list, config_params):
    '''
    Return
        a list
    '''
    ccs_results = []

    # time_for_feature_finding = 0

    # find the target files by the unique id for a compound
    target_info = target_list[target_list.ID==comp_id]
    if target_info.shape[0]==0: return ccs_results

    # get file names for multiple runs
    rep_files = target_info.RawFileName.tolist()
    rep_files.sort()
    num_reps = len(rep_files)

    # get the unique information for each target
    unique_target_info = target_info.drop(['RawFileName', 'FrameMetaName'], axis=1).drop_duplicates()
    if unique_target_info.shape[0] > 1:
        print("[ERROR] There are more than one targets for this comp_id. comp_id:{}, and unique_target_info:".format(comp_id))
        print(unique_target_info)
    compound_id = unique_target_info.iloc[0].CompID
    exact_mass = unique_target_info.iloc[0].ExactMass
    ionization = unique_target_info.iloc[0].Ionization
    neutral_name = unique_target_info.iloc[0].NeutralName

    print(compound_id, neutral_name, ionization, exact_mass)

    # get adducts
    adducts = get_adducts(target_info.ExactMass.tolist()[0], config_params['adducts'])[target_info.Ionization.tolist()[0]]
    
    # get file informations
    tdf = target_info[['RawFileName', 'FrameMetaName']].dropna()
    if tdf.shape[0] == 0:
        print("[ERROR] cannot find any metadata files for", comp_id)
        return ccs_results
    rawFile2Framemeta = pd.Series(tdf.FrameMetaName.values, index=tdf.RawFileName).to_dict()
    print(rawFile2Framemeta)
    ##################################################
    plt.close('all')
    figs = {}
    is_filled = {}
    axis = {}
    for adduct in adducts:
        figs[adduct], axis[adduct] = plt.subplots(num_reps, sharex=True, sharey=True, figsize=(8,3*num_reps))
        is_filled[adduct] = False
    figs['meta'], axis['meta'] = plt.subplots(3, sharex=True, sharey=False, figsize=(8,8))
    figs['intdist'], axis['intdist'] = plt.subplots(config_params['num_fields'], num_reps, sharex=True, sharey=False, figsize=(6*num_reps, 2*config_params['num_fields']))
    ##################################################

    # compute CCS for each replicate
    try:
        for r, rep_file in enumerate(rep_files):
            if files_not_enough(rep_file, config_params, FLAGS.format):
                ccs_prop = dict()
                tokens = comp_id.rsplit('_', 1)
                ccs_prop['Compound_id'] = compound_id
                ccs_prop['Ionization'] = ionization
                ccs_prop['replicate'] = rep_file
                ccs_prop['name'] = neutral_name
                # ccs_prop['CAS'] = list(target_info.CAS)[0]
                ccs_prop['comments'] = "couldn't find some files to compute CCS"
                ccs_results.append(ccs_prop)
                continue
            
            # meta_file = (fname + '{0}.txt').format(config_params['suffix_meta'])
            meta_file = rawFile2Framemeta[rep_file]
            metadata = get_metadata(meta_file, config_params['frame_offset'], ax=axis['meta'], label=rep_file.split('/')[-1])
            
            # collecting features
            features = []
            for step in range(config_params['num_fields']):
                if FLAGS.format=='cef': ffile = (rep_file + '{0}.cef').format(config_params['suffix_raw'].format(step+1))
                else: ffile = (rep_file + '{0}.csv').format(config_params['suffix_raw'].format(step+1))
                
                _features, _ = get_features(ffile, fformat=FLAGS.format)
                if _features.shape[0] > 0:
                    _features['frame'] = np.ones(_features.shape[0], dtype=np.int32) * (step+1)
                    features.append(_features)
                    
                    ## draw m/z vs intensity
                    if num_reps == 1:
                        ax = axis['intdist'][step]
                    else:
                        ax = axis['intdist'][step, r]
                    plot_intensity_distribution(_features, adducts, ax, config_params['mz_tolerance'])
                else:
                    print("[ERROR] This file has no features: {0}".format(ffile))
            
            if len(features) == 0: continue
            features = pd.concat(features)

            # compute CCS for each adducts
            print("#"*150)
            print("# features")
            print("#"*150)
            print(features)
            print("features size:", features.shape)

            for adduct in adducts:
                adduct_mass = adducts[adduct]
                start_time = time.time()
                
                if (FLAGS.maxint):
                    ccs_features_within_mz = find_features_maxint(features, metadata, adduct_mass, config_params['mz_tolerance'])
                else:
                    ccs_features_within_mz = find_features(features, metadata, adduct_mass, config_params['mz_tolerance'],
                                                       threshold_num_isotopes=FLAGS.num_isotopes_threshold,
                                                       threshold_intensity_rank=FLAGS.intensity_rank_threshold)

                if ccs_features_within_mz.shape[0] > 0:
                    print("#"*150)
                    print("# ccs_features_within_mz")
                    print("#"*150)
                    print(ccs_features_within_mz)
                    print("ccs_features_within_mz size:", ccs_features_within_mz.shape)

                    ccs_list = get_possible_ccs_values(ccs_features_within_mz,
                                                       adduct_mass,
                                                       old_drift_tube_length=config_params['old_drift_tube_length'],
                                                       drift_tube_length=config_params['drift_tube_length'],
                                                       neutral_mass=config_params['neutral_mass'],
                                                       threshold_n_fields=3,
                                                       threshold_r2=FLAGS.r2_threshold)
                    # filtering should be done based on ccs values of across all 3 replicates
                    # Note: i am not sure if r2 is a good metric to do this.
                    ccs_list = ccs_filter(ccs_list)

                    if len(ccs_list) > 0:
                        tokens = comp_id.rsplit('_', 1)
                        for ccs in ccs_list:
                            ccs_prop = ccs.to_dict()
                            print("[{0}] {1} ({2}), CCS: {3}({4})".format(comp_id, adduct, rep_file, ccs_prop['ccs'], ccs_prop['r2']))
                            ccs_prop['Compound_id'] = compound_id
                            ccs_prop['Ionization'] = ionization
                            ccs_prop['adduct'] = adduct
                            ccs_prop['replicate'] = rep_file
                            ccs_prop['name'] = neutral_name
                            # ccs_prop['CAS'] = list(target_info.CAS)[0]
                            ccs_results.append(ccs_prop)

                        if num_reps == 1:
                            _tmp_ax = axis[adduct]
                        else:
                            _tmp_ax = axis[adduct][r]

                        ##################################################
                        plot_ccs_regression_lines2(
                            _tmp_ax,
                            adduct,
                            adduct_mass,
                            ccs_features_within_mz,
                            ccs_list,
                            title=rep_file.rsplit("/", 1)[1],
                            drift_tube_length=config_params['drift_tube_length'])
                        is_filled[adduct] = True
                        ##################################################
        ##################################################                
        for adduct in adducts:
            if is_filled[adduct]:
                figs[adduct].tight_layout()
                figs[adduct].savefig(FLAGS.output_dir+"/"+comp_id+"_"+adduct+".pdf", dpi=300)
        
        axis['meta'][0].legend()
        figs['meta'].tight_layout()
        figs['meta'].savefig(FLAGS.output_dir+"/"+comp_id+"_meta.pdf", dpi=300)
        
        figs['intdist'].tight_layout()
        figs['intdist'].savefig(FLAGS.output_dir+"/"+comp_id+'_intensity_dist.pdf')

        ##################################################
    except Exception as e:
        if hasattr(e, 'strerror'):
            print ("[ERROR]: {0} ({1})".format(e.strerror, rep_file))
        else:
            print ("[ERROR]: ", e)
    # print('Total time for feature finding: {0} sec/compound(e.g., 3 reps and 6 adducts)'.format(time_for_feature_finding))
    return ccs_results
            
def compute(df, ion_mz, config_params):
    '''compute ccs
    '''
    params = {}
    params['temp'] = df.ImsTemperature.tolist()
    params['pressures'] = df.ImsPressure.tolist()
    params['voltages'] = (df.ImsField*config_params['old_drift_tube_length']).tolist()  ## 10.869 * (78.12 / 78.236) = 10.853 for correction
    params['arrival_time'] = df.dt.tolist()
    params['neutral_mass'] = config_params['neutral_mass']
    params['drift_tube_length'] = config_params['drift_tube_length']
    params['mass'] = ion_mz
    # print(params)
    ccs, prop = SteppedFieldCCS(params=params).compute()
    # print("CCS:", ccs)
    return prop

def plot_ccs_regression_lines(axis, adduct, adduct_mass, df, prop, title, drift_tube_length=78.236):
    
    addmass = adduct_mass
    color = get_adducts_colors()[adduct]

    p_v = df.ImsPressure / (df.ImsField * drift_tube_length)
    
    p_vmax = p_v.max()
    p_vmin = p_v.min()
    axis.scatter(p_v, df.dt, c=color)
    axis.text(0.05, 0.8, '{0} {1:.6f}'.format(adduct, addmass),
            verticalalignment='bottom', horizontalalignment='left',
            transform=axis.transAxes,
            color='k', fontsize=15)
    for r in df.itertuples():
        axis.text((r.ImsPressure / (r.ImsField * drift_tube_length) + (p_vmax - p_vmin)/7), r.dt,
                  # '{0:.3f}ppm, {1:.2f}(z_score={2:.3f})'.format(mass_error(r.mass, addmass), r.intensity, r.intensity_z),
                  '{0:.3f}ppm, z_score={1:.2f}'.format(mass_error(r.mass, addmass), r.intensity_z),
                  color='k', fontsize=10)

    axis.plot(p_v, 1000 * (prop['intercept'] + prop['slope']*p_v), 'r', label='fitted line')
    axis.text(0.05, 0.65, 'r-squared:{0:.5f}'.format(prop['r_value']**2),
        verticalalignment='bottom', horizontalalignment='left',
        transform=axis.transAxes,
        color='k', fontsize=15)
    axis.text(0.05, 0.5, 'CCS:{0:.4f}'.format(prop['ccs']),
        verticalalignment='bottom', horizontalalignment='left',
        transform=axis.transAxes,
        color='k', fontsize=15)
    axis.set_title(title)
    axis.set_xlabel('Pressure/Voltages (Torr/V)')
    axis.set_ylabel('Arrival time (ms)')

# def plot_ccs_regression_lines2(axis, adduct, adduct_mass, df, prop, title, drift_tube_length=78.236):
def plot_ccs_regression_lines2(
                    axis,
                    adduct,
                    adduct_mass,
                    df,
                    ccs_list,
                    title,
                    drift_tube_length):
    addmass = adduct_mass
    color = get_adducts_colors()[adduct]

    p_v = df.ImsPressure / (df.ImsField * drift_tube_length)
    
    p_vmax = p_v.max()
    p_vmin = p_v.min()
    pv_width = p_vmax - p_vmin

    for r in df.itertuples():
        axis.scatter(r.ImsPressure / (r.ImsField * drift_tube_length), r.dt,
            c=color, s=1000*r.intensity, alpha=0.2)

    axis.text(0.05, 0.8, '{0} {1:.5f}'.format(adduct, addmass),
            verticalalignment='bottom', horizontalalignment='left',
            transform=axis.transAxes,
            color='k', fontsize=10)

    for ccs in ccs_list:
        prop = ccs.to_dict()
        pv = [ccs.pressures[i] / (ccs.fields[i] * drift_tube_length) for i in range(len(ccs.pressures))]
        dt_diff = [abs(ccs.arrival_time[i-1]-ccs.arrival_time[i]) for i in range(1,len(ccs.arrival_time))]
        for i, f in enumerate(ccs.fields):
            axis.text((pv[i] + (p_vmax - p_vmin)/7), ccs.arrival_time[i],
                      '{0:.3f}ppm, z_score={1:.2f}'.format(ccs.mass_ppm_error[i], ccs.intensity_z[i]),
                      color='k', fontsize=10)
            # axis.scatter(pv[i], ccs.arrival_time[i], s=np.log(ccs.intensity_org[i]), c=color)
            axis.scatter(pv[i], ccs.arrival_time[i], s=1000*ccs.intensity[i], c=color, alpha=0.8)
        
        axis.text(min(pv)-2*(p_vmax - p_vmin)/7, min(ccs.arrival_time)-0.8*min(dt_diff),
            'CCS:{0:.4f}(r2:{1:.5f})'.format(prop['ccs'], prop['r2']),
            color='r', fontsize=10)

        axis.plot(p_v, 1000 * (prop['intercept'] + prop['slope']*p_v), 'r', label='fitted line')
    axis.set_title(title)
    axis.set_xlim(left=p_vmin-pv_width*0.5, right=p_vmax+pv_width)
    axis.set_xlabel('Pressure/Voltages (Torr/V)')
    axis.set_ylabel('Arrival time (ms)')


def plot_intensity_distribution(features, adducts_mass, ax, ppm=50):
    if features.shape[0] > 0:
        colors = get_adducts_colors()
        ddata = np.log(features.intensity_org)
        g = sns.kdeplot(ddata, shade=True, color="b", ax=ax)
        ax.axvline(np.log(np.median(features.intensity_org)), linestyle=':')
        ax.axvline(np.log(10*np.median(features.intensity_org)), linestyle=':')
        ax.axvline(np.log(np.mean(features.intensity_org)+2*np.std(features.intensity_org)), linestyle='-.')
        for adduct in adducts_mass:
            sel = features[is_in_tolerance(features.mass, adducts_mass[adduct], ppm)]
            if sel.shape[0] > 0:
                ax.scatter(np.log(sel['intensity_org']), np.zeros(sel.shape[0]), c=colors[adduct])
        ax.set_xlabel('log(Intensity)')
        ax.set_ylabel('Density')
        ax.set_xlim([np.min(ddata), np.max(ddata)])


def report(FLAGS, ccs_table, target_list):
    if ccs_table.shape[0] == 0:
        print("Unfortunately, we couldn't find any good CCS values.")
        return
    def get_stats_adduct(group):
        return {'ccs_avg_adduct': group.mean(), 'ccs_rsd_adduct': 100*group.std()/group.mean(), 'ccs_count_adduct': group.count()}
    def get_stats_file(group):
        return {'ccs_count_file': group.count()}
    
    ccs_avg = ccs_table.groupby(['Compound_id', 'adduct'])['ccs'].apply(get_stats_adduct).unstack()
    ccs_table = pd.merge(ccs_table, ccs_avg.reset_index(), on=['Compound_id','adduct'], how='left')

    ccs_count_file = ccs_table.groupby(['Compound_id', 'adduct', 'replicate'])['ccs'].apply(get_stats_file).unstack()
    ccs_table = pd.merge(ccs_table, ccs_count_file.reset_index(), on=['Compound_id', 'adduct','replicate'], how='left')
    
    print(ccs_table.head())

    # save to a csv file after reordering the columns
    cols = list(ccs_table.columns)
    if 'ccs_avg_adduct' in cols:
        cols.pop(cols.index('ccs_avg_adduct'))
    else:
        ccs_table['ccs_avg_adduct'] = np.nan
    if 'ccs_rsd_adduct' in cols:
        cols.pop(cols.index('ccs_rsd_adduct'))
    else:
        ccs_table['ccs_rsd_adduct'] = np.nan
    cols.pop(cols.index('Compound_id'))
    cols.pop(cols.index('Ionization'))
    cols.pop(cols.index('adduct'))
    cols.pop(cols.index('ccs'))
    cols.pop(cols.index('adduct_mz'))
    cols.pop(cols.index('name'))
    newcols = ['Compound_id','name','Ionization','adduct','adduct_mz','ccs_avg_adduct','ccs_rsd_adduct','ccs']+cols
    
    df = ccs_table[newcols]
    # df = ccs_table
    df.to_csv(FLAGS.output_dir+'/'+FLAGS.output, sep='\t')
    
def multi(FLAGS, config_params):
    if FLAGS.ppm: config_params['mz_tolerance'] = FLAGS.ppm

    os.makedirs(FLAGS.output_dir, exist_ok=True)

    # read a list of targets
    if FLAGS.target_list_file.endswith('.csv'):
        target_list = pd.read_csv(FLAGS.target_list_file)
    else: target_list = pd.read_csv(FLAGS.target_list_file, sep='\t')
    num_targets = target_list.shape[0]
    
    target_list = pd.concat([target_list]*2,ignore_index=True)
    if "Ionization" not in target_list.columns:
        target_list['Ionization'] = ['pos']*num_targets+['neg']*num_targets
    target_list['ID']= target_list.CompID.str.cat("_"+target_list.Ionization)
    target_list = target_list.fillna(method='ffill')

    # find RawFileName
    import re
    suffix_header = config_params['suffix_raw'].split('{',1)[0]
    print(suffix_header)
    uniqueIDs = set(target_list.UniqueID4DfileNames.drop_duplicates().tolist())
    print(uniqueIDs)
    if ("RawFileName" not in target_list.columns) or ("FrameMetaName" not in target_list.columns):
        feature_files = set(glob.glob(FLAGS.feature_files))
        framemeta_files = set(glob.glob(FLAGS.framemeta_files))
        uniqueIDs_list = []
        for _f in feature_files:
            for uid in uniqueIDs:
                if bool(re.search('[-_]?{}[-_]'.format(uid), _f)):
                    if bool(re.search('[-_]?pos[-_]', _f.lower())):
                        _ion = 'pos'
                    else:
                        _ion = 'neg'
                    # print(_f, uid, _ion)
                    # prefix of file names
                    filename = os.path.basename(_f).split(suffix_header)[0]
                    framemeta_name = ""
                    for framemeta in framemeta_files:
                        if filename in framemeta:
                            framemeta_name = framemeta
                    prefix = _f.split(suffix_header)[0]
                    uniqueIDs_list.append({'RawFileName':prefix, 'FrameMetaName':framemeta_name, 'uid':uid, 'ionizations':_ion})
                    break

        tdf = pd.DataFrame(uniqueIDs_list).drop_duplicates()
        target_list = target_list.merge(tdf, left_on=['Ionization','UniqueID4DfileNames'], right_on=['ionizations','uid'])
        del target_list['ionizations']
        del target_list['uid']
        
    # target_list.to_csv('temp.csv')

    ## e.g., S00001.b if you have a same compound id but different versions.
    # num_comp = list(pd.DataFrame(target_list.CompID.str.split('\.').tolist(), columns = ['CompID','ver']).CompID.drop_duplicates())
    compound_ids = target_list.ID.drop_duplicates().tolist()
    num_pos = (target_list.drop_duplicates(subset='ID').Ionization=='pos').sum()
    num_neg = (target_list.drop_duplicates(subset='ID').Ionization=='neg').sum()
    
    # compounds
    assert len(compound_ids) == num_pos+num_neg,\
        "Please check if there are duplicates in CompID and its Ionization"
    print('Number of compounds: {0} (pos:{1}, neg:{2})'.format(len(compound_ids), num_pos, num_neg))
    print(compound_ids)

    ccs_results = []
    start_time = time.time()
    for cid in compound_ids:
        # compute ccs for this compound
        ccs_results += get_ccs(FLAGS, cid, target_list, config_params)
        print('[{0}] {1:.2f} sec'.format(cid, (time.time()-start_time)))
    print('Total time: {0:.2f} sec/compound(e.g., 3 reps)'.format((time.time()-start_time)/len(compound_ids)))
    ccs_table = pd.DataFrame(ccs_results)
    report(FLAGS, ccs_table, target_list)

if __name__ == '__main__':
    FLAGS = parser.parse_args()
    print("options:", FLAGS)

    # read a set of configuration parameters
    config_params = get_config(FLAGS.config_file)
    print(config_params)

    multi(FLAGS, config_params)

import xml.etree.ElementTree
import pandas as pd
import numpy as np

def is_in_tolerance(x, mass, ppm):
    delta = mass * ppm * 1.0e-6
    return (x >= mass - delta) & (x <= mass + delta)

def etree_to_dict(t):
    '''convert XML tree to json format
    '''
    d = {t.tag : map(etree_to_dict, t.iter())}
    d.update(('@' + k, v) for k, v in t.attrib.items())
    d['text'] = t.text
    d['tag'] = t.tag
    return d

def get_config(config_file):
    '''read the config_file
        config_file: file path for a config file
    Return
        config_params: a dict containing a configuration information
    '''
    e = xml.etree.ElementTree.parse(config_file).getroot()
    json = etree_to_dict(e)

    config_params = {}
    adducts_params=[]
    try:
        for j in json['configuration']:
            if 'parameter' in j:
                if j['@name'] in ['mz_tolerance', 'drift_tube_length', 'old_drift_tube_length', 'neutral_mass']:
                    val = float(j['text'])
                elif j['@name'] in ['frame_offset', 'num_fields']:
                    val = int(j['text'])
                elif j['@name'] != 'adducts':
                    val = j['text']
                else:
                    continue
                config_params[j['@name']] = val
            if 'adduct' in j:
                add = {'name':j['@name'], 'charges':j['@charges'], 'mass':float(j['@mass'])}
                adducts_params.append(add)
        config_params['adducts'] = adducts_params
    except Exception as e:
        print(j, config_params)
        raise e

    return config_params


def get_features_from_cef(cef, max_normalize=True):
    '''get features by reading a cef file
    '''
    e = xml.etree.ElementTree.parse(cef).getroot()
    json = etree_to_dict(e.findall('CompoundList')[0])
    idx = 0
    mppid = 0
    rst = []
    mspeaks = []
    in_ms_peaks = False
    for j in json['CompoundList']:
        if 'Compound' in j:
            mppid = j['@mppid']
        if 'Location' in j:
            mz = j['@m']
            rt = j['@rt']
            intensity = j['@y']
            dt = j['@dt']
            rst.append({'mppid':mppid, 'mz':float(mz), 'rt':float(rt), 'intensity':float(intensity), 'dt':float(dt)})
        if 'MSPeaks' in j:
            for k in j['MSPeaks']:
                if ('p' in k):
                    mspeaks.append({'mppid':mppid, 'mass':float(k['@x']), 'intensity_org':float(k['@y']), 'z':k['@z'], 's':k['@s']})
    df = pd.DataFrame(rst)
    mspeaks = pd.DataFrame(mspeaks)
    if df.shape[0] > 0:
        df['intensity_z'] = (df.intensity - df.intensity.mean())/df.intensity.std()
        if max_normalize:
            df.intensity /= df.intensity.max()
        # num_isotopes = mspeaks.mppid.value_counts()
        mspeaks['num_isotopes'] = mspeaks.groupby('mppid')['mppid'].transform('count')
        # print(num_isotopes)
        mspeaks = mspeaks.drop_duplicates(subset='mppid', keep='first')
        # mspeaks['num_isotopes'] = num_isotopes
        
        df = pd.merge(mspeaks, df, left_on="mppid", right_on="mppid", how='inner')
    return df, mspeaks#, num_isotopes


def get_features_from_mzmine_csv(csv, max_normalize=True):
    df = pd.read_csv(csv)
    if df.shape[0] == 0: return df, None
    col_id = [c for c in df.columns if c.endswith('row ID')][0]
    col_area = [c for c in df.columns if c.endswith('Peak area')][0]
    col_height = [c for c in df.columns if c.endswith('Peak height')][0]
    col_mz = [c for c in df.columns if c.endswith('Peak m/z')][0]
    col_dt = [c for c in df.columns if c.endswith('Peak RT')][0]
    if 'calibrated_ccs' in df.columns:
        cols = [col_id, col_mz, col_dt, col_area, col_height, 'calibrated_ccs']
        colnames = ['mppid', 'mass', 'dt', 'intensity', 'height', 'calibrated_ccs']
    else:
        cols = [col_id, col_mz, col_dt, col_area, col_height]
        colnames = ['mppid', 'mass', 'dt', 'intensity', 'height']
    
    df = df[cols].copy()
    df.columns = colnames
    df['intensity_org'] = df.intensity
    if df.shape[0] > 0:
        df['intensity_z'] = (df.intensity - df.intensity.mean())/df.intensity.std()
        if max_normalize:
            df.intensity /= df.intensity.max()
    return df, None
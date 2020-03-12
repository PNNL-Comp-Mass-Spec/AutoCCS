import pandas as pd
import numpy as np
from utils import is_in_tolerance

__N2__ = 28.01340
__He__ = 4.002602

"""
calibrate CCS from calibrants for ion mobility data (IMS and SLIM)
"""


class CCSCalibrator:
    """Calibrator to estimate CCS values"""

    def __init__(self, calibrant_file, neutral_mass=28.01340):
        self.neutral_mass = neutral_mass
        self.calibrant_file = calibrant_file
        self.calibrate_fn = None
        self.calibrated_features = None

        self.calibrants = self.ready()

    def fit(self, x, y, deg=1):
        """curve fitting to find the calibration function
        """
        poly = np.poly1d(np.polyfit(x, y, deg))
        self.calibrate_fn = poly

        # r-squared
        # fit values, and mean
        yhat = poly(x)  # or [p(z) for z in x]
        ybar = np.sum(y) / len(y)  # or sum(y)/len(y)
        ssreg = np.sum((yhat - ybar) ** 2)  # or sum([ (yihat - ybar)**2 for yihat in yhat])
        sstot = np.sum((y - ybar) ** 2)  # or sum([ (yi - ybar)**2 for yi in y])
        r = ssreg / sstot

        return poly, r

    def ready(self):
        """rread calibrants and add reduced mass and ccs values based on the drift gas
        """
        calibrants = pd.read_csv(self.calibrant_file, sep='\t')
        calibrants['reduced_mass'] = (calibrants['m/z'] * calibrants.z * self.neutral_mass) / (
                    self.neutral_mass + calibrants['m/z'] * calibrants.z)
        calibrants['reduced_ccs'] = calibrants.CCS / (calibrants.z / (calibrants.reduced_mass ** 0.5))
        return calibrants

    def find_calibrant_features(self, features, ppm=20, ionization='pos'):
        """find the features"""
        selected = pd.DataFrame()
        if self.calibrants.shape[0] > 0:
            for row in self.calibrants[self.calibrants.Ionization == ionization].iterrows():
                mz = row[1]['m/z']
                reduced_ccs = row[1]['reduced_ccs']
                ccs = row[1]['CCS']
                _ff = features[is_in_tolerance(features.mass, mz, ppm)]
                _ff = _ff.assign(calibrants_mz=mz, ccs=ccs, reduced_ccs=reduced_ccs)
                if selected.shape[0] == 0:
                    selected = _ff
                else:
                    selected = selected.append(_ff)
            # TODO: now simply select the most intensive peaks
            temp = selected.sort_values(by='intensity_org').drop_duplicates(subset=['calibrants_mz'], keep='last')
            return temp
        else:
            return None

    def calibrate(self, features, calibrate_fn):
        """calibrate CCS
        """
        features['reduced_mass'] = (features['m/z'] * features.z * self.neutral_mass) / (
                    self.neutral_mass + features['m/z'] * features.z)
        features['reduced_ccs'] = calibrate_fn(features.dt)
        features['calibrated_ccs'] = features.reduced_ccs * (features.z / (features.reduced_mass ** 0.5))
        self.calibrated_features = features
        return features

    def export(self, features, fout, export_fn):
        export_fn(features, fout)

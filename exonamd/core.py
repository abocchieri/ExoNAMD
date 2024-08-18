import numpy as np
import pandas as pd


def compute_amdk(m, e, di, a):
    return m * np.sqrt(a) * (1 - np.sqrt(1 - e**2) * np.cos(np.deg2rad(di)))


def compute_namd(amdk, m, sqrt_a):
    return np.sum(amdk) / np.sum(m * sqrt_a)

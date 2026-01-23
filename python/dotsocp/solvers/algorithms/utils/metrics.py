import numpy as np


def norm_l2(x, h):
    return np.sqrt(h) * np.linalg.norm(x)


def fnorm_l2(x, h):
    return np.sqrt(h) * np.linalg.norm(x)

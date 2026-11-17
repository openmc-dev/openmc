"""Miller-style D contour; no unsupported elliptical-torus primitive."""
import numpy as np


def contour(config,minor=None,elongation=None,triangularity=None):
    t=np.linspace(0,2*np.pi,config.contour_points,endpoint=False)
    a=config.minor_radius if minor is None else minor
    k=config.elongation if elongation is None else elongation
    d=config.triangularity if triangularity is None else triangularity
    # Match the installed native TokamakSource's delta convention exactly.
    return np.column_stack([config.major_radius+a*np.cos(t+d*np.sin(t)),k*a*np.sin(t)])

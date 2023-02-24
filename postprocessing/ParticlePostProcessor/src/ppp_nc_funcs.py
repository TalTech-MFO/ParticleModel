#!/usr/bin/env python3

import datetime
import numpy as np
from netCDF4 import Dataset, num2date


def get_var(fname, vname, indices=None):
    with Dataset(fname) as f:
        if (indices is not None):
            return np.array(f.variables[vname][indices])
        return np.array(f.variables[vname][:])


def get_time(fname, vname='time'):
    with Dataset(fname) as f:
        timein = f.variables[vname][:]
        timeunits = f.variables[vname].units
    return np.array([datetime.datetime(val.year, val.month, val.day, val.hour, val.minute, val.second) for val in num2date(timein, timeunits)])


def get_file_attributes(fname, attr_name):
    """Get global file attributes"""
    with Dataset(fname) as f:
        return list(f.__dict__[attr_name])


def get_dimensions(fname):
    with Dataset(fname) as f:
        return f.dimensions.keys()

#!/usr/bin/env python3

from .ppp_nc_funcs import get_var
import numpy as np


PARTICLE_FILLVALUE = 9.e36
PARTICLE_FILLVALUE_INT = -9999


class ParticleFile:
    """
    Particle file class.
    Stores particle data from a single file.
    """

    def __init__(self, filename, id_list=None, ini_file=None):
        self.filename = filename
        self._ini_file = ini_file
        self._lon = None
        self._lat = None
        self._depth = None
        self._state = None
        self._indices = None
        self._id = None

        self._nids = None
        self._ntime = None
        self._nparticles = None

        self._id_list = id_list
        if (id_list is not None):
            self.filter_ids()
        return

    def filter_ids(self):
        id_mask = np.sum(
            [self.id == id for id in self._id_list], axis=0, dtype=bool)
        self._id = self.id[id_mask]
        self._indices = (..., id_mask)

    def positions(self):
        x = self.lon.flatten()
        y = self.lat.flatten()
        return np.concatenate((x[:, np.newaxis], y[:, np.newaxis]), axis=1)

    @property
    def nids(self):
        return len(np.unique(self.id))

    @property
    def nparticles(self):
        if (self._nparticles is None):
            self.shape
        return self._nparticles

    @property
    def ntime(self):
        if (self._ntime is None):
            self.shape
        return self._ntime

    @property
    def shape(self):
        if (self._ntime is None or self._nparticles is None):
            self._ntime, self._nparticles = self.lon.shape
        return (self._ntime, self._nparticles)

    @property
    def id(self):
        if (self._id is None):
            self._id = get_var(self.filename, "id", indices=self._indices)
            self._id[self._id > PARTICLE_FILLVALUE] = np.nan
        return self._id

    @property
    def lon(self):
        if (self._lon is None):
            self._lon = get_var(self.filename, "lon", indices=self._indices)
            self._lon[self._lon > PARTICLE_FILLVALUE] = np.nan
        return self._lon

    @property
    def lat(self):
        if (self._lat is None):
            self._lat = get_var(self.filename, "lat", indices=self._indices)
            self._lat[self._lat > PARTICLE_FILLVALUE] = np.nan
        return self._lat

    @property
    def depth(self):
        if(self._depth is None):
            self._depth = get_var(
                self.filename, "depth", indices=self._indices)
            self._depth[self._depth > PARTICLE_FILLVALUE] = np.nan
        return self._depth

    @property
    def state(self):
        if(self._state is None):
            self._state = get_var(
                self.filename, "state", indices=self._indices)
            self._state[self._state < PARTICLE_FILLVALUE_INT] = -1
        return self._state



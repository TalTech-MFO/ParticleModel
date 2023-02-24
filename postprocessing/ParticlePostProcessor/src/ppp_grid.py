#!/usr/bin/env python3

from .ppp_nc_funcs import get_var, get_dimensions
import numpy as np


class Grid:
    """Base class for grid data"""

    m2deg = 1./(1852.*60.)
    deg2m = 1852.*60.

    def __init__(self, filename):
        self.filename = filename
        self.lon = None
        self.lat = None
        self.depth = None
        self.type = None
        self.read_file()
        return

    @property
    def coords(self):
        pass

    @property
    def dims(self):
        pass

    def create_mesh(self):
        pass

    def set_resolution(self, resolution):
        pass

    def read_file(self):
        dims = get_dimensions(self.filename)
        # Check that the file has "lon" and "lat" dimensions
        # "lon" could also be "lonc" or "longitude"
        # "lat" could also be "latc" or "latitude"
        lon_name = list(set(dims).intersection(
            ["lon", "lonc", "longitude"]))[0]
        lat_name = list(set(dims).intersection(
            ["lat", "latc", "latitude"]))[0]
        depth_name = "bathymetry"
        self.lon = get_var(self.filename, lon_name)
        self.lat = get_var(self.filename, lat_name)
        self.depth = get_var(self.filename, depth_name)

    def counts(self, positions, binx, biny):
        print("Calculating counts...")
        bins = np.array([binx, biny])
        c, _ = np.histogramdd(positions, bins=bins)
        return c.T

    def counts_by_state(self, positions, state, binx, biny):
        print("Calculating counts by state...")
        state = state.flatten()
        states = sorted(list(set(np.abs(state))))
        nstates = len(states)

        nx = len(binx)
        ny = len(biny)
        bins = np.array([binx, biny])

        nanmask = np.isnan(positions)
        nanmask = nanmask[:, 0] + nanmask[:, 1]
        positions = positions[~nanmask]
        state = state[~nanmask]

        counts = np.zeros((nstates, ny-1, nx-1), dtype=int)
        for istate in range(nstates):
            kkstate = np.squeeze(np.where(state == states[istate]))
            c, _ = np.histogramdd(positions[kkstate], bins)
            counts[istate, :, :] = c.T
        return counts

    def counts_by_id(self, positions, id, binx, biny):
        print("Calculating counts by id...")
        id = id.flatten()
        ids = sorted(list(set(id)))
        nids = len(ids)

        nx = len(binx)
        ny = len(biny)
        bins = np.array([binx, biny])

        nanmask = np.isnan(positions)
        nanmask = nanmask[:, 0] + nanmask[:, 1]
        positions = positions[~nanmask]
        id = id[~nanmask]

        counts = np.zeros((nids, ny-1, nx-1), dtype=int)
        for iid in range(nids):
            kkid = np.squeeze(np.where(id == ids[iid]))
            c, _ = np.histogramdd(positions[kkid], bins)
            counts[iid, :, :] = c.T
        return counts


class HorizontalGrid(Grid):
    """Class to store horizontal grid data (maps)"""

    def __init__(self, filename, resolution=None):
        self.filename = filename
        super().__init__(filename)
        self.type = "map"
        self.dlon = None
        self.dlat = None
        self.dx = None
        self.dy = None
        self.cell_area = None
        self.resolution = resolution
        if (resolution is not None):
            self.set_resolution(resolution)
        self.create_mesh()
        return

    @property
    def coords(self):
        # ! np.histogramdd reduces size of each dimension by 1
        return {"lat": self.lat[:-1], "lon": self.lon[:-1]}

    @property
    def dims(self):
        # ! np.histogramdd reduces size of each dimension by 1
        return {"lat": len(self.lat)-1, "lon": len(self.lon)-1}

    def create_mesh(self):
        """Create meshgrid in meters"""
        self.dlon = np.diff(self.lon)
        self.dlat = np.diff(self.lat)
        dlon_m, dlat_m = np.meshgrid(self.dlon, self.dlat)
        self.dx = dlon_m * self.deg2m * \
            np.cos(np.nanmean(self.lat) * np.pi/180.)
        self.dy = dlat_m * self.deg2m
        self.cell_area = self.dx * self.dy

    def set_resolution(self, resolution):
        dlon_bin = resolution*self.m2deg / \
            np.cos(np.nanmean(self.lat) * np.pi/180.)
        self.lon = np.arange(
            np.nanmin(self.lon), np.nanmax(self.lon) + dlon_bin, dlon_bin)
        dlat_bin = resolution*self.m2deg
        self.lat = np.arange(
            np.nanmin(self.lat), np.nanmax(self.lat) + dlat_bin, dlat_bin)

    def get_counts(self, particle_file, sort):
        if (sort == "all"):
            return self.counts(particle_file.positions())
        elif (sort == "id"):
            return self.counts_by_id(particle_file.positions(),
                                     np.tile(particle_file.id, (particle_file.ntime, 1)))
        elif (sort == "state"):
            return self.counts_by_state(particle_file.positions(), particle_file.state)

    def counts(self, positions):
        return super().counts(positions, self.lon, self.lat)

    def counts_by_state(self, positions, state):
        return super().counts_by_state(positions, state, self.lon, self.lat)

    def counts_by_id(self, positions, id):
        return super().counts_by_id(positions, id, self.lon, self.lat)

    def landmask(self):
        pass


class VerticalGrid(Grid):
    """Class to store vertical grid data (profiles)"""
    pass

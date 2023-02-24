#!/usr/bin/env python3

import os
from .ppp_nc_funcs import get_file_attributes, get_dimensions
import numpy as np
from netCDF4 import Dataset


class ResultsFile:
    """
    File class to store results
    """

    def __init__(self, filename, original_file):
        self.filename = filename
        self.original_file = original_file
        self.lon = None
        self.lat = None
        self.depth = None
        self.dims = None
        self.coords = None
        self.type = None
        # Outgoing data
        self.data = None
        self.exists = os.path.exists(filename)
        if (self.exists):
            self.check_dimensions()

    def initialize(self, grid, groups):
        self.type = grid.type
        self.dims = grid.dims
        self.coords = grid.coords
        self.data = {group: [] for group in groups}

        # Check that the dimensions and coordinates are the same
        for dimname, dimsize in self.dims.items():
            if (dimname not in self.coords.keys()):
                raise Exception(
                    f"Results file has dimension {dimname} that is not in grid")
            if (len(self.coords[dimname]) != dimsize):
                raise Exception(
                    f"Results file has dimension {dimname} with wrong size")

        # Create the file
        print(f"Creating results file {self.filename}")
        ncout = Dataset(self.filename, "w", format="NETCDF4")
        for dimname, dimsize in self.dims.items():
            ncout.createDimension(dimname, dimsize)
        for varname, var in self.coords.items():
            ncout.createVariable(varname, "float32", (varname,))
            ncout.variables[varname][:] = var
        # Set global attributes
        ncout.setncatts({"original_file": self.original_file,
                         "type": self.type})
        ncout.close()

    def append(self, data, group):
        """Append data dictionary to outgoing data"""
        self.data[group].append(data)

    def write(self):
        print(f"Writing results to file {self.filename}")
        ncout = Dataset(self.filename, "a", format="NETCDF4")
        for data_group in self.data.keys():
            for data_item in self.data[data_group]:
                # Check that the dimensions and coordinates are the same
                for dimname in self.dims.keys():
                    if (dimname not in data_item["dims"]):
                        raise Exception(
                            f"Results file has dimension {dimname} that is not in data")
                # Create the variable
                for i_layer, (varname, varval) in enumerate(data_item["name_dict"].items()):
                    full_varname = f"{data_group}_{varname}"
                    ncout.createVariable(
                        full_varname, "float32", tuple(data_item["dims"].keys()))
                    # Write the data
                    if (data_item["data"].ndim == 3):
                        ncout.variables[full_varname][:] = data_item["data"][i_layer, :, :]
                    else:
                        ncout.variables[full_varname][:] = data_item["data"]
                    # Add attributes
                    ncout.variables[full_varname].setncatts(data_item["attrs"])
        ncout.close()

    def check_dimensions(self):
        """
        Get dimensions of file
        """
        self.dims = get_dimensions(self.filename)
        if (len(self.dims) != 2):
            raise Exception("Results file has wrong number of dimensions")
        if (len(set(self.dims).intersection(["lon", "lat", "depth"])) != 2):
            raise Exception(f"Results file has wrong dimensions: {self.dims}")

    def read(self):
        """
        Read results from file
        """
        with Dataset(self.filename) as f:
            self.type = get_file_attributes(self.filename, "type")
            if (self.type == "map"):
                self.lon = np.array(f.variables["lon"][:])
                self.lat = np.array(f.variables["lat"][:])
            elif (self.type == "profile"):
                # TODO: Should be able to have a profile that is not just along the GOF thalweg
                self.depth = np.array(f.variables["depth"][:])
                self.lon = np.array(f.variables["lon"][:])
            countvars = []
            concvars = []
            for vname in f.variables.keys():
                if (vname.__contains__("counts")):
                    countvars.append(vname)
                if (vname.__contains__("concentration")):
                    concvars.append(vname)
            self.ncounts = len(countvars)
            self.nconcentration = len(concvars)
            counts = [[]]*self.ncounts
            concentration = [[]]*self.nconcentration
            for ivar, vname in enumerate(countvars):
                counts[ivar] = np.array(f.variables[vname][:])
            for ivar, vname in enumerate(concvars):
                concentration[ivar] = np.array(
                    f.variables[vname][:])
            self.counts = np.array(counts)
            self.concentration = np.array(concentration)

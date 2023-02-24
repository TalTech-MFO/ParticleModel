#!/usr/bin/env python3

class Quantity:
    """
    Base class for quantities
    """

    def __init__(self):
        self.grid = None
        self.units = None
        self.name = None
        self.dims = None
        self.coords = None

    def compute(self, counts):
        raise NotImplementedError("Quantity.compute not implemented")

    def data_to_dict(self, data, name_dict):
        d = {}
        d["data"] = data
        d["name"] = self.name
        # TODO: Here a long name could be added
        d["coords"] = self.coords
        d["dims"] = self.dims
        # Name/value dictionary
        nvd = {f"{self.name}{i}": val for i, val in enumerate(
            name_dict)} if name_dict is not None else {f"{self.name}": "all"}
        d["name_dict"] = nvd
        # if (len(nvd) != data.shape[0]):
        #     raise ValueError(
        #         f"Number of names ({len(nvd)}) does not match number of data ({data.shape[0]})")
        d["attrs"] = {"units": self.units,
                      "name": self.name, "name_dict": str(nvd)}
        return d

    def run(self, counts):
        # TODO: Add some sort of name to the counts
        data = self.compute(counts["counts"])
        return self.data_to_dict(data, counts["name_dict"])


class Concentration(Quantity):
    """
    Class to compute concentration
    """

    def __init__(self, grid):
        super().__init__()
        self.grid = grid
        self.units = "particles/m2"
        self.name = "concentration"
        self.dims = grid.dims
        self.coords = grid.coords

    def compute(self, counts):
        return counts/self.grid.cell_area


class Counts(Quantity):
    """
    Class to compute counts
    """

    def __init__(self, grid):
        super().__init__()
        self.grid = grid
        self.units = "particles/cell"
        self.name = "counts"
        self.dims = grid.dims
        self.coords = grid.coords

    def compute(self, counts):
        return counts

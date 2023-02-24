#!/usr/bin/env python3

import argparse
import numpy as np
from .ppp_particle_file import ParticleFile
from .ppp_result_file import ResultsFile
from .ppp_grid import HorizontalGrid, VerticalGrid
from .ppp_quantity import Concentration, Counts
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)


def plot(args):
    print("Plotting")
    filename = args.source
    result_filename = args.out_file
    raise NotImplementedError("Plotting not implemented yet")


def get_state_name_numbers(particle_file):
    return np.unique(np.abs(particle_file.state))


def get_id_name_numbers(particle_file):
    return sorted(np.unique(particle_file.id))


def process(args):
    print("Processing")
    overwrite = args.overwrite
    filename = args.source
    result_filename = args.out_file
    topo_filename = args.topo
    resolution = args.resolution
    id_list = args.id_list
    # TODO: Default should be a list of all options (loop and compute all)
    sort_by = args.sort
    ini_file = args.ini_file
    last_positions = args.last

    result_file = ResultsFile(result_filename, filename)

    if (result_file.exists and overwrite == False):
        raise Exception(
            f"Results file {result_filename} already exists. Use -O to overwrite")

    particle_file = ParticleFile(filename, id_list, ini_file)
    grid = HorizontalGrid(topo_filename, resolution)
    result_file.initialize(grid, groups=sort_by)

    # Create name dictionary
    name_dict = {"state": get_state_name_numbers(particle_file),
                 "id": get_id_name_numbers(particle_file),
                 "all": None}

    counts = {sort_type: {
        "counts": grid.get_counts(particle_file, sort_type),
        "name_dict": name_dict[sort_type]
    } for sort_type in sort_by}

    # Compute measures
    OrdinaryCounter = Counts(grid)
    ConcentrationCounter = Concentration(grid)
    for sort_type in sort_by:
        data = OrdinaryCounter.run(counts[sort_type])
        result_file.append(data, sort_type)

        data = ConcentrationCounter.run(counts[sort_type])
        result_file.append(data, sort_type)

    # Save results
    result_file.write()

    return


def main():
    # Global parser
    global_parser = argparse.ArgumentParser(add_help=True)
    subparsers = global_parser.add_subparsers(
        help="Subcommands", dest="subcommand")

    # Plot parser
    plot_parser = subparsers.add_parser(
        "plot", help="Plot particle distribution")
    plot_parser.add_argument(
        "-s", "--source", help="Path to data file", type=str, required=True)
    plot_parser.add_argument(
        "-o", "--out-file", help="Results file", type=str, required=False, default="counts.nc")
    plot_parser.set_defaults(func=plot)

    # Process parser
    process_parser = subparsers.add_parser(
        "process", help="Process particle file")
    process_parser.add_argument(
        "-O", "--overwrite", help=f"Overwrite results file", action='store_true')
    process_parser.add_argument(
        "-s", "--source", help="Path to data file", type=str, required=True)
    process_parser.add_argument(
        "-o", "--out-file", help="Results file", type=str, required=False, default="counts.nc")
    process_parser.add_argument(
        "--topo", help="Topography (grid) file", type=str, required=False, default="topo.nc")
    process_parser.add_argument(
        "-dx", "--resolution", help="Resolution of grid (meters)", type=float)
    process_parser.add_argument(
        "--sort", help="Sort particles", choices=["all", "id", "state"], nargs="*", type=str, required=False, default=["all", "id", "state"])
    process_parser.add_argument(
        "--id-list", help="List of particle IDs to use in processing", nargs="+", type=float)
    process_parser.add_argument(
        "--ini-file", help="Initial particle position file", type=str, required=False)
    process_parser.add_argument(
        "--last", help=f"Use only last active position", action='store_true')
    process_parser.set_defaults(func=process)

    args = global_parser.parse_args()
    args.func(args)

    return 0


if __name__ == "__main__":
    main()

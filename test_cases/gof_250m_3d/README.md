## GoF 3D test case

This test case is here to demonstrate working with large datasets.
- Hydrodynamic and biogeochemical data is from the GETM high resolution model of the Gulf of Finland. (250 x 250 m)
- Particle inputs are from rivers in the Gulf of Finland.

To run the test case: 
- Download the data from


    Siht, E. (2024). GoF 3D test case data [Data set].
    Zenodo. https://doi.org/10.5281/zenodo.13935461

    **Make sure you download version v2. The data is approximately 8.9 GB.**

    Unzip the downloaded dataset in the `data` directory.

- Unzip `input.data.zip` and `topo.data.zip` in the `data` directory.

- Run the script:
    ```bash
    ./run
    or
    bash run
    ```

This will run the simulation and postprocessing. The output will be in the `out` directory. Histograms in `out/counts` and figures in `out/figs`. Note that there is not much to see in the figures because the simulation runs for only two days.



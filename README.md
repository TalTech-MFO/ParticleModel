# ParticleModel


## Installing

Download the code:

```
git clone https://github.com/TalTech-MFO/ParticleModel.git
```

Configure and build:

```
cd ParticleModel
cmake -B build
cmake --build build
```

There are no dependencies other than netCDF. You should have `nf-config` available in your path.

Before running `cmake`, have a look at the `include/output.h` file. Select output variables by commenting/uncommenting lines in `include/output.h`. Particle position will always be written.

## Running the model

To run the model, run the `ParticleModel` executable in your run directory:
```bash
/path/to/executable/ParticleModel [-nml general_namelist_file] [-bnml biofouling_namelist_file]
```
The `-nml` and `-bnml` flags are optional and specify the general and biofouling namelist files, respectively. If not specified, the model will look for `input.nml` and `biofouling.nml` in the run directory.
The namelist templates for setting the model parameters are in the `nml` directory, along with a helper script to generate the namelists from environment variables (requires the [f90nml](https://github.com/marshallward/f90nml) package). For setting the namelist parameters, you can copy and modify the run scripts in the `scripts/run` directory.


## TODO

See TODO.txt for known issues/ideas for further development. 

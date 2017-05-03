Lorenz '96 model
================

This repository contains an implementation of the Lorenz '96/'95 model for investigating reduced precision data assimilation. The three level model is integrated to produce a truth run, and observations of this truth are generated by adding random noise. These observations are then assimilated into a two level version using an ensemble square-root filter, where the third level is replaced by a parametrization. The forecast and analysis steps can be run at reduced precision, using the [Reduced Precision Emulator](https://github.com/aopp-pred/rpe) so that the quality of the assimilation can be investigated as precision is reduced.

All questions on this repository should be directed to Sam Hatfield (samuel.hatfield@physics.ox.ac.uk).

## Building

### Install dependencies
- [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) - all output is written using NetCDF. It is assumed that the NetCDF lib files are contained in `/usr/lib` and the include files are contained in `/usr/include`. The location may well change depending on your platform.

- [rpe](https://github.com/aopp-pred/rpe) - in the root of this repository, run

```bash
git clone https://github.com/aopp-pred/rpe.git
cd rpe
make
```

- Jinja2 - assuming Anaconda Python is installed, run

```bash
conda install jinja2
```

The precision of the model/assimilation can then be set when running `make`.
- For full (double) precision, simply run `make`
- For single precision, run `make single`
- For reduced precision, run `make rpe`. By default, half precision variables are emulated.

Also, since I don't really understand makefiles yet, best do a `make clean` before each rebuild.

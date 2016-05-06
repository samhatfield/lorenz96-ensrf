Lorenz '95 model
================

This repository contains various implementations of the Lorenz '95 model (also called the Lorenz '96 model) for investigating data assimilation. Two filters can be used, the ensemble Kalman filter and the ensemble square root filter.

![Lorenz '95 run example](lorenz95.png)

## Building

Full precision version:
```
make
```

Version with reduced precision for forecast ensemble model and assimilation:
```
make rpe
```

Also, since I don't really understand makefiles yet, best do a `make clean` before each rebuild.

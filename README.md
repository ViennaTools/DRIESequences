# DRIE Models

This repository contains geometric models for several different DRIE processes.

## Dependencies (install automatically)

* [ViennaLS](https://github.com/ViennaTools/viennals)

## Usage


```bash
cmake -B build
cmake --build build
```

Execute the models by running:

```bash
cd build
./DEM2D
./DEM3D
./DREM3D
./DREAM
```

Note: The size of the DREM model has been reduced, so it can be executed on most common processors.


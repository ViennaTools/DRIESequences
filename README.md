# DRIE Models

This repository contains geometric models for several different DRIE processes.

## Usage

Compile all models by issuing:

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

Execute the models by running:

```bash
./DEM2D
./DEM3D
./DREM3D
./DREAM
```

Note: The size of the DREM model has been reduced, so it can be executed on most common processors.

## Dependencies

The dependencies are included in the folder ViennaTools, so no additional setup is necessary. However, the dependencies versions are frozen to their respective versions and do not receive updates and fixes.

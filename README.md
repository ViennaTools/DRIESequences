# DRIE Models

This repository contains geometric models for several different DRIE processes.

## Dependencies

* [ViennaLS](https://github.com/ViennaTools/viennals)

## Usage

The CMake configuration automatically checks if ViennaLS is installed. If CMake is unable to find it, the dependency will be built from source with the _buildDependencies_ target.

```bash
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/your/ViennaLS/install/
make buildDependencies # this will install all dependencies and might take a while
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


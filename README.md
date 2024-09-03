# simple-lattice-boltzmann

This example is the "simplest possible" lattice Boltzmann method (LBM)
implementation that uses OpenMP for parallelization on CPU and HDF5
for self-describing file output. The simulation is seeded with random
noise, which quickly diffuses away, leaving large-scale, small-amplitude
features behind.

To build the code, assuming HDF5 is available, you can do the following,
which sets up a build directory and the Ninja build system.

```sh
cmake -G Ninja -B build .
cd build && ninja
```

The executable `simplelb` will produce an output file, `simple_lbm.h5`, which
can be previewed using the provided Python script:

```sh
python preview.py simple_lbm.h5
```

<!-- vim: set ft=markdown: -->

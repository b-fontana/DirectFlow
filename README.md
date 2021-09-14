## DirectFlow

Simulation of charged particles:
- under the influence of a magnetic field
- imposing a sudden deviation to mimic what happens at ALICE's interaction point

#### Prerequisites

- [ROOT](https://root.cern/). This code was tested with ```ROOT```'s [```conda``` installation](https://root.cern/install/#conda).
- [boost](https://www.boost.org/) + ```sudo apt-get install libboost-program-options-dev``` (Ubuntu)

#### Compiling

```bash
git clone git@github.com:b-fontana/DirectFlow.git
cd DirectFlow
make
```

Upon successful compilation you should get

```bash
Executable v1_beam.exe created.
```

To clean the object files and executable:

```bash
make clean
```

#### Running

```
./v1_beam.exe --mode euler --x 0.8 --y 0.8 --energy 1380 --nparticles 10000 --zcutoff 50.
```

or, in parallel:

```
parallel --ungroup --jobs 7 ./v1_beam.exe --mode euler --x 0.0 --y 0.8 --energy 1380 --nparticles 500000 --zcutoff 5000 --mass_interaction 0.139 --npartons {} ::: 1 10 200
```

#### Plotting

```
python python/npartons_distribution.py --mode euler --x 0.0 --y 0.8 --energy 1380 --step_size 10
```

or similar to the other plotting macros stored under ```python/```.



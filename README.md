## DirectFlow

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
./v1_beam.exe --mode euler
```

**Arguments**

- ```mode```: numerical differential equation solver [euler/rk4]


#### Plotting

```
python python/fit_circle.py rk4
```

**Arguments**

- ```mode```: see above




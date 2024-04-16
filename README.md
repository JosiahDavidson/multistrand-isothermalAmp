# About

        ___  ___      _ _   _     _                       _ 
        |  \/  |     | | | (_)   | |                     | |
        | .  . |_   _| | |_ _ ___| |_ _ __ __ _ _ __   __| |
        | |\/| | | | | | __| / __| __| '__/ _` | '_ \ / _` |
        | |  | | |_| | | |_| \__ \ |_| | | (_| | | | | (_| |
        \_|  |_/\__,_|_|\__|_|___/\__|_|  \__,_|_| |_|\__,_|


Multistrand is a nucleic acids kinetic simulator, and is developed by the
Winfree group at the California Institute of Technology in Pasadena, California
(USA). Until 2013, development was lead by Joseph Schaeffer (now Autodesk).
The project is currently maintained by Jake Kaslewicz (Riedel group, University
of Minnesota, Minneapolis-Saint Paul) and Boyan Beronov (Condon group,
University of British Columbia, Vancouver).

[Official website](http://www.multistrand.org/)

[Changelog](CHANGELOG.md)

Please use GitHub issues for technical problems and questions, and direct more
general inquiries to help@multistrand.org.

## Licence

    Multistrand nucleic acid kinetic simulator
    Copyright (c) 2008-2024 California Institute of Technology. All rights reserved.
    The Multistrand Team (help@multistrand.org)

Using this software is permitted for academic non-commercial purposes only. All
copyright is retained by Caltech.

**Disclaimer:** This software is provided "as is", without warrenty of any kind,
express or implied, including but not limited to the warrenties of
merchantability, fitness of a particular purpose and noninfringement. In no
event shall the authors or copyright holders be liable for any claim, damages or
other liability, whether in an action of contract, tort or otherwise, arising
from, out of or in connection with the software or the use or other dealings in
the software.

## Contributors
* Erik Winfree (winfree@caltech.edu)
* Chris Thachuk
* Frits Dannenberg (fdann@caltech.edu)
* Chris Berlind
* Joshua Loving
* Justin Bois
* Joseph Berleant
* Joseph Schaeffer
* Jake Kaslewicz (kasle001@umn.edu)
* Boyan Beronov (beronov@cs.ubc.ca)

# Usage

## Installation
### Requirements
| Dependency                        | Notes              |
|-----------------------------------|--------------------|
| C++11                             | gcc 8+ or clang 8+ |
| Python                            | 3.9+               |
| [NUPACK](https://www.nupack.org/) | 4.0.1+             |
 
The `numpy` and `scipy` Python packages are installed automatically as
dependencies, and `matplotlib` is added if the package extra `[tutorials]` is
specified (see `setup.cfg` for details).
 
### Linux
 - Make sure the requirements above are installed on your host system.
 - `git clone` this repository into your workspace.
 - Run `pip install .` in the Multistrand source directory.

### macOS
 - Make sure the requirements above are installed on your host system.
 - Install `xcode` commandline tools.
 - Install Python through `homebrew`.
 - Follow the Linux installation steps.
 - In `~/.bash_profile`, edit the `$PYTHONPATH` to include
   `/Library/Python/<python version>/site-packages`.
 
### Windows
 - Make sure the requirements above are installed on your host system.
 - Follow the instructions for installing the latest version of the [Microsoft
   C++ Build Tools](https://wiki.python.org/moin/WindowsCompilers).
 - Follow the Linux installation steps.
 
### [Apptainer](https://apptainer.org/) container
 - [Install Apptainer](https://apptainer.org/docs/admin/latest/installation.html).
 - Place `nupack-<version>.zip` into the parent folder of the Multistrand
   source directory, so that the container build script can unpack and install NUPACK.
 - `$> cd tools`
 - [Build the container](
   https://apptainer.org/docs/user/latest/build_a_container.html):
   `$> sudo apptainer build multistrand.sif multistrand.def`
 - In order to enable Multistrand to write its outputs inside a container which
   is read-only by default:
   - [create a container overlay](
     https://apptainer.org/docs/user/latest/persistent_overlays.html)
     (recommended for new users): `$> apptainer overlay create --sparse --size
     1024 multistrand.img`
   - and/or choose a host filesystem path to [bind into the container](
     https://apptainer.org/docs/user/latest/bind_paths_and_mounts.html)
     (recommended for transferring simulation results onto the host system, and
     for development).
 - [Start the container](
   https://apptainer.org/docs/user/latest/quick_start.html#interacting-with-images):
   - with an overlay: `$> apptainer shell --cleanenv --contain
     --pwd /dna/multistrand --overlay multistrand.img multistrand.sif`
   - and/or with bind paths: `$> apptainer shell --cleanenv --contain
     --pwd /dna/multistrand --bind <src>:<dest> multistrand.sif`


## Documentation
For an overview of Multistrand's functionality, see the built-in documentation:

```python
from multistrand import objects, options, system
help(objects)
help(options)
help(system)
```

Further documentation can be found in `doc/` and `tutorials/`, and tutorial
files are organized as follows. The folder `under_the_hood/` contains in-depth
tutorials, and Jupyter versions are located in `under_the_hood_notebooks/`. The
folder `case_hybridization/` contains a case study into hybridization kinetics.
Additional demo files are located in `misc/`.


## Development
### Source tree
The Multistrand library is located under `src/`. `test/` is the test suite, and
`tools/` provides Apptainer container definitions and maintenance scripts.

### Testing
To execute the currently maintained portion of the test suite (including some of
the small tutorials):

 - Install the test dependencies: `$> pip install ".[testing]"`
 - Run: `$> pytest`


## Examples
As a very quick primer, we discuss two small scripts below.

### Hybridization trajectory
A quick test to see if Multistrand is working is to run the following script,
which simulates the hybridization of two complementary strands after their
initial collision, ending the simulation when the two strands either completely
hybridize or seperate. Of the example trajectories below, the first hybridizes
successfully, while the second dissociates within a few elementary steps.

```sh
$> python tutorials/misc/sample_trace.py
========================================================================
trajectory_seed = -3377281620794262823
------------------------------------------------------------------------
                        |   t[us]   | dG[kcal/mol] |     state_seed     
------------------------------------------------------------------------
    GCGTTTCAC+GTGAAACGC
[1] ...(.....+....).... | 0.0000000 |    +1.673    | (25464,12418,37992)
[1] ...(..(..+..).).... | 0.1134000 |    +2.273    | (25330,16011,59063)
    GTGAAACGC+GCGTTTCAC
[1] ..(((....+...)).).. | 0.1503000 |    +1.408    | (61436,61430,25585)
[1] ...((....+...)).... | 0.1703000 |    +0.250    | (11542,54358, 8882)
[1] ..(((....+...)))... | 0.1997000 |    +1.032    | (16576,27017, 9385)
[1] ..(((..(.+.).)))... | 0.2131000 |    +0.475    | (62842,33280, 5339)
[1] ..(((.((.+.)))))... | 0.2661000 |    -0.744    | (22980,50665,23442)
[1] ..(((.(((+))))))... | 0.9345000 |    -2.158    | (24606,22501,64078)
[1] ...((.(((+))))).... | 1.5820000 |    -2.940    | (32520,47683,58666)
[1] ..(((.(((+))))))... | 2.0870000 |    -2.158    | (20738,57299,45486)
[1] ...((.(((+))))).... | 2.5680000 |    -2.940    | (13452,27588,53292)
[1] ...((.((.+.)))).... | 2.7350000 |    -1.526    | (60454,38452,21172)
[1] ...((.(((+))))).... | 3.0790000 |    -2.940    | (15952,34652,51923)
[1] ...(..(((+))).).... | 3.2310000 |    -1.844    | (38282,20850,19636)
[1] ...((.(((+))))).... | 3.2450000 |    -2.940    | (41044,37316, 7361)
[1] ...(..(((+))).).... | 3.7210000 |    -1.844    | (61742,15254,44110)
[1] ...((.(((+))))).... | 4.0800000 |    -2.940    | (40600,46018,62617)
[1] ..(((.(((+))))).).. | 4.3240000 |    -1.782    | (58130,16946,28623)
[1] ...((.(((+))))).... | 4.3460000 |    -2.940    | (48412,34475, 6102)
[1] ..(((.(((+))))))... | 4.6460000 |    -2.158    | (36662, 8321,33119)
[1] ...((.(((+))))).... | 4.7990000 |    -2.940    | (49120,42281,27025)
[1] ...((.((.+.)))).... | 5.7420000 |    -1.526    | (22938,26572,58370)
[1] ..(((.((.+.)))))... | 6.2830000 |    -0.744    | (43748, 1086,59949)
[1] ...((.((.+.)))).... | 6.3640000 |    -1.526    | (58942,10766,24330)
[1] ..(((.((.+.)))))... | 6.6260000 |    -0.744    | (49704,22664,45642)
[1] ...((.((.+.)))).... | 6.9970000 |    -1.526    | ( 6434,10971,57601)
[1] ..(((.((.+.)))))... | 7.1220000 |    -0.744    | (35244, 3014,28541)
[1] ..((..((.+.)).))... | 7.2360000 |    +0.352    | ( 5702,47741,31168)
[1] ...(..((.+.)).).... | 7.2670000 |    -0.430    | (50544,34458,17639)
[1] ...(.(((.+.)))).... | 7.2760000 |    -1.526    | (16810,40799, 6485)
[1] .....(((.+.)))..... | 7.5880000 |    -3.948    | (31092,53904,56764)
[1] ....((((.+.)))).... | 7.8000000 |    -5.215    | (16206,45486,51455)
[1] ..(.((((.+.)))).).. | 7.8810000 |    -5.530    | (59832,55390,30635)
[1] ..((((((.+.)))))).. | 7.9320000 |    -9.325    | (62258, 3388,19229)
[1] .(((((((.+.))))))). | 8.3220000 |    -8.372    | (39484,48749,45001)
[1] ((((((((.+.)))))))) | 8.4360000 |   -10.965    | (33110,64936,36987)
[1] (((((((((+))))))))) | 8.9150000 |   -12.379    | (20224,11160,45541)


========================================================================
trajectory_seed = 970235241
------------------------------------------------------------------------
                        |   t[us]   | dG[kcal/mol] |     state_seed     
------------------------------------------------------------------------
    GCGTTTCAC+GTGAAACGC
[1] ..(......+........) | 0.0000000 |    +0.014    | (25464,15762,41396)
[1] ..((.....+.......)) | 0.3017000 |    +0.506    | (25330,14875,31083)
    GTGAAACGC+GCGTTTCAC
[1] ........(+..)...... | 0.3618000 |    +0.014    | (61436,48646,11671)
[1] .......((+..))..... | 1.2300000 |    +0.506    | (11542,63718,55232)
[1] ........(+..)...... | 1.9950000 |    +0.014    | (16576, 2201,60132)
    GCGTTTCAC GTGAAACGC
[2] ......... ......... | 2.3960000 |    +0.000    | (62842,24464,37661)
```

### Hybridization rates
The following script estimates the hybridization rate for a strand and its
complement. The computation relies on "first step mode" and will only work if
NUPACK is correctly installed. Alternatively, dissociation rates can be computed
by using `dissociation` as the first commandline argument. In that case, the
dissociation rate is computed indirectly by estimating the association rate, and
working out the dissociation rate from the partition function (e.g. `k_forward /
k_backward = exp(-dG / RT)` where `dG` is the partition function for the
complex, `R` is the gas constant and `T` is the temperature).

```sh
$> python tutorials/compute/rate.py hybridization AGCTGA -bootstrap
--------------------------------------------------------------------------------
2024-04-13 19:17:53   Starting Multistrand 2.2  (c) 2008-2024 Caltech

Running first step mode simulations for AGCTGA (with Boltzmann sampling)...

Start states:
Complex:
         Name: 'automatic0'
     Sequence: AGCTGA
    Structure: ......
      Strands: ['top']
    Boltzmann: True
  Supersample: 1

Complex:
         Name: 'automatic1'
     Sequence: TCAGCT
    Structure: ......
      Strands: ['top*']
    Boltzmann: True
  Supersample: 1

Stop conditions:
Stop Condition, tag:SUCCESS
  Sequence  0: AGCTGA+TCAGCT
  Structure 0: ((((((+))))))

Stop Condition, tag:FAILURE
  Sequence  0: AGCTGA
  Structure 0: ......

Using Results Type: FirstStepRate
Computing 1000 trials, using 10 threads ..
 .. and rolling 100 trajectories per thread until 500 successful trials occur.

nForward = 443
nReverse = 357

nForward = 546
nReverse = 454

Found 546 successful trials, terminating.
Done.  0.25239 seconds -- now processing results

The hybridization rate of AGCTGA and the reverse complement is 3.00e+06 /M /s
Bootstrapping FirstStepRate, using 1200 samples.    ..finished in 1.39 sec.

Estimated 95% confidence interval: [2.84e+06,3.16e+06]
Computing took 2.7906 s
```

### Log files
Multistrand automatically creates a logfile that contains some information on
the used model, e.g.:

```
$> cat multistrandRun.log
--------------------------------------------------------------------------------
Multistrand 2.2

sodium        :  1 M
magnesium     :  0 M
temperature   :  298.15 K
rate method   :  3  (1: Metropolis, 2: Kawasaki, 3: Arrhenius)
dangles       :  1  (0: none, 1: some, 2: all)
GT pairing    :  1  (0: disabled, 1: enabled)
concentration :  1 M
Nupack params :  /opt/bitnami/python/lib/python3.11/site-packages/nupack/parameters/dna04-nupack3.json

Kinetic parameters:
  type          End        Loop       Stack  StackStack     LoopEnd    StackEnd   StackLoop
  A      1.2965e+01  1.6424e+01  1.4184e+01  8.0457e+00  2.4224e+00  1.7524e+01  5.8106e+00
  E      3.4980e+00  4.4614e+00  5.2869e+00 -6.2712e-01  8.4934e-02  2.6559e+00 -1.1276e+00
  R      1.3584e+06  5.3080e+07  3.7097e+04  8.0873e+07  9.5396e+01  2.1242e+11  5.0138e+06

        dS_A        dH_A     biScale        kUni
 -0.0000e+00 -0.0000e+00  1.6006e-02  1.0000e+00

Rate matrix [ concentration . k_bi . k_uni(l,r) ]:
  2.1742e+04  1.3591e+05  3.5931e+03  1.6777e+05  1.8221e+02  8.5980e+06  4.1772e+04
  1.3591e+05  8.4962e+05  2.2461e+04  1.0487e+06  1.1390e+03  5.3747e+07  2.6112e+05
  3.5931e+03  2.2461e+04  5.9378e+02  2.7724e+04  3.0111e+01  1.4209e+06  6.9031e+03
  1.6777e+05  1.0487e+06  2.7724e+04  1.2945e+06  1.4059e+03  6.6342e+07  3.2231e+05
  1.8221e+02  1.1390e+03  3.0111e+01  1.4059e+03  1.5269e+00  7.2053e+04  3.5006e+02
  8.5980e+06  5.3747e+07  1.4209e+06  6.6342e+07  7.2053e+04  3.4000e+09  1.6518e+07
  4.1772e+04  2.6112e+05  6.9031e+03  3.2231e+05  3.5006e+02  1.6518e+07  8.0252e+04
```


# Frequently asked questions

## Capabilities
**Q:** Can I simulate leak reactions using Multistrand?

**A:** Yes. We have now added a preliminary tutorial, see `tutorials/leak_casestudy`.

## Troubleshooting
**Q:** How do I adjust the solvent salt concentrations?

**A:** Like so. (units are M = mol / litre)

```python
from multistrand.options import Options
o1 = Options()
o1.sodium = 0.05
o1.magnesium = 0.0125
```

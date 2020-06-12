# LitMod3D 4.0

**LitMod3D** is a software for 3D integrated geophysical-petrological interactive modelling of the lithosphere and underlying upper mantle using a variety of input datasets: potential fields (gravity and magnetic), surface heat flow, elevation (isostasy), seismics, magnetotellurics and geochemical.

**installation**
Inst macOS:
Install https://www.xquartz.org
Install macports https://www.macports.org/install.php

Install gmt4 sudo port install gmt4
Add path to gmt4 in .zshrc


**LitMod3D** is a software for 3D integrated geophysical-petrological interactive modelling of the lithosphere and underlying upper mantle using a variety of input datasets: potential fields (gravity and magnetic), surface heat flow, elevation (isostasy), seismics, magnetotellurics and geochemical.

Ver 3.0 incorporates a highly optimised Python thermal solver (bi-conjugate gradient squared method), crustal petrology features (thermodynamic equilibrium and metastable) and a parallel gravity forward solver. The new version is intended to work with program get-inp (customized interface to Perple_X, http://www.perplex.ethz.ch/) to generate the inputcrustal and mantle compositional files.

Ver 4.0 has a built-in dispersion-curves calculator.



## Installation
Download the distibution archive (not yet available) to your computer and unpack it into a folder (for example `~/litmod`). 

### Linux
You have to have [Python 3](https://www.python.org/), [PGPLOT](http://www.astro.caltech.edu/~tjp/pgplot/), [gfortran](https://gcc.gnu.org/fortran/), [zsh](http://www.zsh.org/) and [GMT4](http://gmt.soest.hawaii.edu) installed on your computer. In Linux, set **zsh** as a main shell. Change directory to the folder with unpacked files. 

### macOS
In **macOS**, first install [iTerm2](https://iterm2.com). Within **iTerm2** terminal install [brew](https://brew.sh) with command:
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
```

Than install **zsh** with command :
```
brew install zsh
```

Make **zsh** a default shell in **iTerm2** by clicking iTerm2 -> Profiles -> Command: /bin/zsh.

Restart **iTerm2** and install [Oh My Zsh](https://ohmyz.sh) with command:
```
sh -c "$(curl -fsSL https://raw.github.com/robbyrussell/oh-my-zsh/master/tools/install.sh)"
```

Now, install [XQuartz](https://www.xquartz.org) from the latest .dmg package:


After, install [pgplot](http://www.astro.caltech.edu/~tjp/pgplot/) with command:
```
brew install pgplot
```


Finally, to compile **LitMod**, use:
```
make
```
Note that if you already had **pgplot** and/or **XQuartz** installed on your Mac, you may have to change the path to libraries. Instructions above were tested on macOS Catalina 10.15.5.


## Usage

### Preparing GMT4
LitMod executables and main calling script `LITMOD_3D.job` use [GMT4](http://gmt.soest.hawaii.edu), and therefore **GMT4** must be set up to be called directly from the the local LitMod folder, such as:
```
minmax -C ./layers/layer1.xyz
```

To install **GMT4** on **macOS**, get [MacPorts](https://www.macports.org/).

Than run:
```
sudo port install gmt4
```

In **Linux**, follow instructions on http://gmt.soest.hawaii.edu to install **GMT4**.


**IN PROGRESS**




### Input files
Examples of input files required to run **LitMod** are included into the archive (https://github.com/javfurchu/litmod/releases/download/v3.0/litmod_3.0.zip)

LitMod requires a number of input files that are managed from the calling script `LITMOD_3D.job`:

- `LITMOD3D.info` model parameters such as mesh size and the thermal boundary conditions.
- `mnt.info` lists the compositional files to be used in the thermodynamic calculations.

In addition other input files are required:

- `layers.info` header file that describes the physical and petrological properties of each layer in the model.
- Geometry of model layers in the folder 'layers' (geographical coordinates),and the folder 'layers_xy' (Cartesian coordinates).
- Compositional files in the folder 'mant_data'.

### LITMOD_3D.job
`LITMOD_3D.job`  contains all the input values required to run LiMod. The script interfaces with [GMT](https://github.com/GenericMappingTools) and generates `LITMOD3D.info` and `mnt.info` files rquired by LitMod. `LITMOD_3D.job`  also processes the files contained in the folder 'layers' (geographical layers) and puts the output (Cartesian layers) in the folder 'layers_xy' as required by LitMod. 

At the beginning of `LITMOD_3D.job` the user can set up the geopgraphical boundaries of the modelling region (lon_min, lon_max, lat_min, lat_max), the number of nodes in the model (N_x, N_y, N_z) and other parameters. 

`LITMOD_3D.job` preprocesses input grids for geophysical data (e.g., gravity anomaly, surface topography, etc.) and reprojects them onto the modelling region (within lon_min, lon_max, lat_min, lat_max) with defined grid spacing (as per N_x, N_y variables). To make a first run, set parameter pre_pro to 1 and run `LITMOD_3D.job`:
```
./LITMOD_3D.job
```
The distibution provides sample observed data grids in the folder `example_obs`. The user can provide customized data grids with the format lon; lat; value that will be preprocessed by `LITMOD_3D.job`. For the gravity gradients a small program (LNOF2MRF) is provided to rotate the tensor from the Local North Oriented Reference Frame (https://earth.esa.int/web/guest/data-access/view-data-product/-/article/goce-gravity-gradients-in-lnof-5775) to the model reference frame in the Cartesian coordinate system used by LitMod. 


After running `LITMOD_3D.job` in preprocessing mode, switch pre_pro variable to 0 and run `LITMOD_3D.job` again. This will run **LitMod** forward modelling code. The ouput geophysical data sets and the 3D lithospheric model can be visualized using **LitMod** graphical interface, which also can be used to interactively modify the 3D model:
```
./litmod_intf
```

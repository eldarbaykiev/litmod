# LitMod3D 4.0

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

### Preparing GMT6
LitMod executables and main calling script `LITMOD_3D_4.0.job` use [GMT6](http://gmt.soest.hawaii.edu). **GMT6** should be installed prior the usage of **LitMod**.

Check the version of **GMT** by running:
```
gmt
```

in the console (for both Mac and Linux).


### Input files
Examples of input files required to run **LitMod** are included into the archive (https://github.com/javfurchu/litmod/releases/download/v3.0/litmod_3.0.zip)

LitMod requires a number of input files that are managed from the calling script `LITMOD_3D.job`:

- `LITMOD3D.info` model parameters such as mesh size and the thermal boundary conditions.
- `mnt.info` lists the compositional files to be used in the thermodynamic calculations.

In addition other input files are required:

- `layers.info` header file that describes the physical and petrological properties of each layer in the model.
- Geometry of model layers in the folder 'layers' (geographical coordinates),and the folder 'layers_xy' (Cartesian coordinates).
- Compositional files in the folder 'mant_data'.

### LITMOD_3D_4.0.job
`LITMOD_3D_4.0.job`  contains all the input values required to run LiMod. The script interfaces with [GMT](https://github.com/GenericMappingTools) and generates `LITMOD3D.info` and `mnt.info` files rquired by LitMod. `LITMOD_3D_4.0.job`  also processes the files contained in the folder 'layers' (geographical layers) and puts the output (Cartesian layers) in the folder 'layers_xy' as required by LitMod. 

At the beginning of `LITMOD_3D_4.0.job` the user can set up the geopgraphical boundaries of the modelling region (lon_min, lon_max, lat_min, lat_max), the number of nodes in the model (N_x, N_y, N_z) and other parameters. 

`LITMOD_3D_4.0.job` preprocesses input grids for geophysical data (e.g., gravity anomaly, surface topography, etc.) and reprojects them onto the modelling region (within lon_min, lon_max, lat_min, lat_max) with defined grid spacing (as per N_x, N_y variables). To make a first run, set parameter pre_pro to 1 and run `LITMOD_3D_4.0.job`:
```
./LITMOD_3D_4.0.job
```
The distibution provides sample observed data grids in the folder `example_obs`. The user can provide customized data grids with the format lon; lat; value that will be preprocessed by `LITMOD_3D_4.0.job`. For the gravity gradients a small program (LNOF2MRF) is provided to rotate the tensor from the Local North Oriented Reference Frame (https://earth.esa.int/web/guest/data-access/view-data-product/-/article/goce-gravity-gradients-in-lnof-5775) to the model reference frame in the Cartesian coordinate system used by LitMod. 


After running `LITMOD_3D_4.0.job` in preprocessing mode, switch pre_pro variable to 0 and run `LITMOD_3D_4.0.job` again. This will run **LitMod** forward modelling code. The ouput geophysical data sets and the 3D lithospheric model can be visualized using **LitMod** graphical interface, which also can be used to interactively modify the 3D model:
```
./litmod_intf
```


### MINEOS grid calculation

To calculate dispersion curves in each grid point, open `LITMOD_3D_4.0.job` and set 
```
mineos_whole_grid=1
```
This would turn on the caluclation of dispersion curves in each node of the model (N_x, N_y). 

You are required to have `ME01.dat` and `ME01b.dat` files with reference reference 1D models with Vp, Vs, rho for z > 400 km for isotropic and anisotropic cases in the folder.

After running `LITMOD_3D_4.0.job`, run

```
python3 periods.py 10.0
```

to create  grids with veolcities and velocity anomalies for 10.0 second period. For other period, enter another number as an argument instead of 10.0. 

Four grids in format LON LAT VALUE would be produced as a result:
`MINEOS_out_grid_love_10.0.xyz`,  `MINEOS_out_grid_ray_10.0.xyz`, `MINEOS_out_grid_love_10.0_anom.xyz`, `MINEOS_out_grid_ray_10.0_anom.xyz`

Grids can be plotted using GMT script   `plot_MINEOS_out_grid.sh`


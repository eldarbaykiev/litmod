#!/bin/sh

unamestr="$(uname)"
#gmt=/usr/local/Cellar/gmt/6.0.0_5/bin/gmt
gmt=gmt



#Region coordinates
lon_min=-12
lon_max=-5
lat_min=51
lat_max=56
step=0.2



layer=MINEOS_out_grid_love_10.0
#projection
proj='X6i/6i'

$gmt surface ${layer}.xyz -Rd$lon_min/$lon_max/$lat_min/$lat_max -I${step}d -G${layer}.grd

$gmt grd2cpt ${layer}.grd -Cpolar -I> ${layer}.cpt



$gmt gmtset PS_MEDIA=a2
$gmt gmtset ANNOT_FONT_SIZE_PRIMARY=25p
$gmt grdimage ${layer}.grd -Rd$lon_min/$lon_max/$lat_min/$lat_max -J${proj} -C${layer}.cpt -P -K > ${layer}.ps
$gmt pscoast -Dl -Rd$lon_min/$lon_max/$lat_min/$lat_max -J${proj} -Ba5g5f1/a5g5f1WeSn -V -Di -Wthick -O -K >> ${layer}.ps
#$gmt grdcontour ${layer}.grd -Rd$lon_min/$lon_max/$lat_min/$lat_max -J${proj} -A2.5+f15p -O -K >> ${layer}.ps
$gmt psscale -D6.3i/3i/15c/1c -C${layer}.cpt -I0  -B5:"vel":/:"": -O >> ${layer}.ps
$gmt ps2raster -A+r ${layer}.ps


layer=MINEOS_out_grid_love_10.0_anom
#projection
proj='X6i/6i'

$gmt surface ${layer}.xyz -Rd$lon_min/$lon_max/$lat_min/$lat_max -I${step}d -G${layer}.grd

$gmt grd2cpt ${layer}.grd -Cpolar -I> ${layer}.cpt



$gmt gmtset PS_MEDIA=a2
$gmt gmtset ANNOT_FONT_SIZE_PRIMARY=25p
$gmt grdimage ${layer}.grd -Rd$lon_min/$lon_max/$lat_min/$lat_max -J${proj} -C${layer}.cpt -P -K > ${layer}.ps
$gmt pscoast -Dl -Rd$lon_min/$lon_max/$lat_min/$lat_max -J${proj} -Ba5g5f1/a5g5f1WeSn -V -Di -Wthick -O -K >> ${layer}.ps
#$gmt grdcontour ${layer}.grd -Rd$lon_min/$lon_max/$lat_min/$lat_max -J${proj} -A2.5+f15p -O -K >> ${layer}.ps
$gmt psscale -D6.3i/3i/15c/1c -C${layer}.cpt -I0 -B5:"Anom":/:"%": -O >> ${layer}.ps
$gmt ps2raster -A+r ${layer}.ps

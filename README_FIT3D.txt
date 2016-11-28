# FIT3D
#
#

FIT3D comprises routines to fit stellar continuum and gas emission
lines from optical spectra in general and IFU data in particular.
It also include some simple visualization tools, and tools to
extract images, slices and individual spectra too.

The routines to fit emission lines are the more simple ones, and are
well described in the 2005-IAC Winter School Tutorial, which can be
found in the following webpage:

http://www.caha.es/sanchez/tutorial_IFUs.html

The repository where to find the different required packages and datasets
is the following one.

ftp://ftp.caha.es/CALIFA/docs/tutorials

#
# Installing the required scripts.
#

To install the different packages (and in particular FIT3D), pleae follow
the instructions in the "INSTALL_ALL.txt" file at the repository.

In addition, you should download the "tutorial.tgz" and "FIT3D_EXAMPLES.tgz"
tar files, where a certain number of configuration files and examples are included.

We will assume that all the required files are installed in the "/SOFT3D" directory
along this document.

#
# DEFINITIONS
#

We name "spectrum" to an ASCII file comprising three columns:

ID   Wavelength    FLUX_INTENSITY

We defined "RSS" file or "Row Stacked Spectrum" a 2D fits file comprising NY-individual
spectra of the same length (NX), with a normalized wavelength step and initial wavelength
defined by the header keywords:
CRPIX1  => Pixel of the starting wavelength  
CRVAL1  => Starting Wavelength in AA
CDELT1  => Wavelength step in AA per pixel

FIT3D assumes that the wavelength starts at the origin of each pixel.


A "RSS" file is normally related to a "position-table" file, as
defined in E3D (sanchez et al. 2005). A position table is an ASCII file 
with NY+1 rows. The first of these rows defines the geometry of the IFU-spaxel
(fiber/lensarray...), while the remaininig ones defined the position (in the
corresponding uints, normally arcsec).

The file has 4 columns. For the 1st row they correspond to:

"SHAPE" (C,H,R)   SIZE_X  SIZE_Y   RUNNING_ID

For the remaining ones they include the position projected in the sky:

J(1...NY)   POS_X  POS_Y   RUNNING_ID

E.g., 

C 0.075 0.075 0
1 0.78 0 1
2 -0.878 -0.169 1
3 0 -0.676 1
4 0.78 -0.676 1
5 -0.78 -0.676 1
...
327 0 0.676 1
328 0.585 0.676 1
329 -0.683 0.507 1
330 -0.683 0.169 1
331 0.683 0.169 1

We define a "datacube" as 3D-fitsfile with the 2 first dimensions related
to the spatial distribution of the pixels (X,Y), and a 3rd direction comprising
the wavelength information of the IFU data (Z). The datacubes, stored in this
way, comprises data with the same spatial step in the X and Y direction (normally
in arcsec), and the same step in wavelength in the Z-direction (normally AA).

The WCS information is stored in the corresponding header entries:

CRPIX1  => Pixel of the starting position in RA (or X)  
CRVAL1  => Starting position in RA (or X)
CDELT1  => Spatial extension of a pixel in the X-direction (in arcsec normally)
CRPIX2  => Pixel of the starting position in DEC (or Y)  
CRVAL2  => Starting position in DEC (or Y)
CDELT2  => Spatial extension of a pixel in the Y-direction (in arcsec normally)
CRPIX3  => Pixel of the starting wavelength  
CRVAL3  => Starting Wavelength in AA
CDELT3  => Wavelength step in AA per pixel

In the case of CALIFA, the datacubes comprises 4 extensions, each one including the
flux-intensity data, the standard-deviation of the noise, a 3D mask (1=good data),
and a weighting cube comprising the information of the co-variance between adjacent
pixel after interpolation.

Other more complicated ways to store the data comprises the Euro3D format (which
comprises a set of FITs-tables, very flexible but complex to handle).
Please read the E3D user guide.

We define a "slice" as an ASCII file comprising the 2D information of a certain parameter.
It has a similar structure as a position table, but the fourth parameter is the value to
be spatially sampled. E.g.,

C 0.075 0.075 0
1 0.78 0 0.5
2 -0.878 -0.169 0.6
3 0 -0.676 0.8

Slices are convinient ways to store the 2D information in non regular grid distribution of
data, like for example the ones usually provided by a fiber-based IFU.

Slices can be interpolated to "maps", ie., 2D images comprising the spatial distribution
information of a particular parameter.


#
# Format of the SSP library templates to be used by FIT3D
#

The template libraries used by FIT3D are RSS, for which in each row we store a single stellar
population. The SSP spectra should have a common reference wavelength calibration, 
with the same step in AA per pixel and the same staring wavelength defined by the header
entries:

CRPIX1
CRVAL1
CDELT1

And for each SSP, it is required a "NAME*" header keyword, with the following format:

NAME0 | spec_ssp_00.0891_z0004.spec 
NAME1 | spec_ssp_00.0891_z019.spec 
NAME2 | spec_ssp_00.0891_z030.spec 
NAME3 | spec_ssp_00.4467_z0004.spec 

i.e., spec_ssp_AGE_zMET.spec

Where AGE should be the SSP age in Gyrs, and MET is the metallicity expressed as the
decimals, ie.,

Z/H = 0.0004  should be 0004
Z/H = 0.02    should be 02

These entries are mandatory to allow the programs to understand the inputs.

EXAMPLE, read the entries of different provided templates:
dump_header.pl miles_20.fits
dump_header.pl miles_6.fits


#
# Cube visualization and data extraction
#

1) Simple visualization of single spectra, 2D RSS or 3D datacubes,
   using command line routines.

Single spectum are defined above. They are interesting intermediate
dataproducts that sometimes requires to be visualized. A simple script
to visualize them is:

spec_plot.pl
USE: spec_plot.pl INPUT_FILE DEVICE  [MIN MAX] [WMIN WMAX] [ONE]

e.g.,
spec_plot.pl spec_NGC2916_peak.txt 66/xs

You can visualize a list of spectra, to compare them:

ls spec_NGC2916*.txt > list_of_spectra.txt
spectra_plot.pl list_of_spectra.txt 66/xs

NOTE: You can scale by a certain value one spectrum in the list by adding the corresponding
multiplicative factor in the file "list_of_spectra.txt".

It is possible to transform an ASCII spectrum to a FITS-file (RSS with a single spectrum), by
using:

spec2img.pl
USE: spec2img.pl spec.txt image.fits ny [crval] [cdelt]

or a list of spectra to a single RSS-file using:

spectra2img.pl
USE: spectra2img.pl list_of_spectra.txt image.fits
spectra2img.pl list_of_spectra.txt test_rss.fits 5000

Or normalizing them at a certain wavelength:

spectra2img_norm.pl
USE: spectra2img_norm.pl list_of_spectra.txt image.fits wave_norm [F_NORM]
spectra2img_norm.pl list_of_spectra.txt test_norm.fits 5000

You can visualize the RSS-files comprises in a certain RSS one-by-one:

spec2D_plot.pl test_norm.fits
USE: spec2D_plot.pl INPUT_FILE.FITS NY DEVICE [MIN MAX] [WMIN WMAX]

spec2D_plot.pl test_norm.fits 0 66/xs

Or you can visualize all the spectra together:
spec2D_plot_all.pl test_rss.fits 0 66/xs -0.1 3

NOTE: Some operations with RSS-files can by found in the R3D user guide.

In a similar way it is possible to visualize single spectrum extracted from a 
cube:
spec_cube_plot.pl
USE: spec_cube_plot.pl INPUT_FILE.CUBE.FITS NX NY DEVICE [MIN MAX] [WMIN WMAX]

e.g.,
spec_cube_plot.pl NGC2916.V500.rscube.fits 25 20 66/xs

This script can be used to extract also a single spectrum from the datacube, since the
considered spectrum is stored in the ASCII file "spec_cube_plot.out".

To extract a 2D slice or integrated map over a certain wavelength range,
you have to use the following scripts:

get_slice.pl
USE: get_slice.pl INPUT_CUBE.fits PREFIX CONF_FILE
CONF_FILE: NAME START_W END_W

This scripts allows to extract multiple "2D maps" from a datacube, by selecting
the initial and final wavelength range, and providing with the mean value along
the covered spectral range spaxel-to-spaxel.

The output files are defined in a considered configuration file including a prefix
together with the wavelength range for the integration. E.g.,

get_slice.pl NGC2916.V500.rscube.fits NGC2916 slice.config
3 slices to cut
Reading cube
NGC2916_map_Ha_SII_6500_7000.fits saved
NGC2916_map_V_4500_5500.fits saved
NGC2916_map_B_3900_4500.fits saved

ds9 NGC2916_map_*.fits &

It is possible to extract all the integrated flux (sum) along the considered wavelength range:
get_slice_sum.pl

Or the convolved cube with a considered filter-response (to reconstruct a certain band image):
get_slice_filter.pl NGC2916.V500.rscube.fits V_Johnson.txt map_V_NGC2916.fits

ds9 NGC2916_map_*.fits  map_V_NGC2916.fits &

NOTE: From filter convolved datacubes it is possible to cross-check the spectrophotometry,
band-by-band, zero-points or possible color effects in the cubes.

It is possible to extract also an ASCII slice from an RSS+PT:

get_slice_filter_rss.pl mos_NGC2916.V500.rss.fits mos_NGC2916.V500.pt.txt V_Johnson.txt slice_NGC2916_V.txt 

Exercise: Test,
get_slice_filter_rss.pl mos_NGC2916.V500.rss.fits mos_NGC2916.V500.pt.txt V_Johnson.txt slice_NGC2916_V.txt none 1

It is possible to visualize the result with the command:

plot_slice.pl slice_NGC2916_V.txt 66/xs C 2.68

You can compare it with a single PPAK pointing:

plot_slice.pl single_slice_NGC2916_V.txt 66/xs C 2.68

To visualize 2D plots extract from a datacube you can use:

plot_maps.pl
USE: plot_maps.pl map.fits min max bright contrast Label factor dev [MASK_MAP]

plot_maps.pl NGC2916_map_Ha_SII_6500_7000.fits -0.1 2 0.5 0.5 'map' 1 66/xs

Or overplot them with a continuum image:
plot_maps_cont.pl

USE: plot_maps_cont.pl map.fits min max bright contrast Label factor dev VAL_TO_WHITE MAP_CONT NLEVELS MIN STEP NWHITE[49] [TABLE REVERSE] [ROT]

plot_maps_cont.pl NGC2916_map_Ha_SII_6500_7000.fits -0.1 1 0.5 0.5 'map' 1 66/xs 0 NGC2916_map_V_4500_5500.fits 10 0.1 0.1 0

Finally, you can directly visualize the content of the datacubes, slice by slice:

plot_maps_cube.pl NGC2916.V500.rscube.fits -0.1 2 0.5 0.5 'map' 1 66/xs 3900 4200


All these command lines are useful to interact with the data, but it is more convinient to use GUI interfaces to explore the data.

   
2) E3D / PINGSoft / IFSview

 The instructions of how to use E3D (tk_e3d.tcl) are included in the
"E3D User Guide" ("E3D.pdf" document), included in the FTP
repository. A complementary information with exercises can be found in
the FTP repository, in the "tutorial.ps.gz" file (comprising the
lessons given in the 2005 IAC Winter School), in section 5.2.


 PINGSoft : Notes by F.F. Rosales-Ortega

 IFSview : It is a visualization tools created in Python, still under development,
which tries to mimic the functionality of E3D. It is not completed or documented.


 2.1) RSS / datacubes.

E3D allows to visualize both datacubes and rss files. For doing so it is required
to use the "Import" options. Let's start open E3D:

tk_e3d.tcl &

Open the "Spaxels Inspector" (Spaxels Menu->Open), and the "Spectral Inspector" (Spectra Menu -> Open).
Then Import a datacube (File Menu -> Import Cube), e.g., NGC2916.V500.rscube.fits. You will see a
complete black screen. Change the Bright and Contract (left side Menu), and the range of shown spectra 
(left side Menu, at the bottom). In the Configuration menu change the colormap to Rainbow. Finally, 
select a different Min-Max range (-0.1, 0.5). To zoom in a particular region, use the left-mouse-button,
clicking a certain wavelength range, and the right one to select it, or move through the spectra using the 
central mouse button clicked. You will see the corresponding "slice-map" in the spaxels inspector.

Using a similar method you can select a particular spectrum in the "spaxels inspector" and you will see
the corresponding spectra in the spectra inspector.

NOTE: We will follow the content in Section 5.2 of the tutorial.ps.gz file, for the subsequent 
examples and exercises.

 2.2) Extraction of a single spectra.

In the "Spaxel Inspector", select Select->Clear. The click on a particular region on the visualized map,
and with the right-mouse-button, you will see the corresponding spectrum in the spectral inspector.

Now on the Spectra Inspector you save it using "File->Save Ascii Table". The corresponding file "spec_1.txt",
can be visualized using "spec_plot.pl":

spec_plot.pl spec_1.txt 66/xs

Now, we will select a "pseudo-slit RSS spectra". For doing so, we go back to the "Spaxel Inspector" and click
the key "s" once. You have changed now from the single spectrum selection mode to the slit selection. Click
once on the "Spaxels Inspector" (left-mouse button), and then move on the screen and click again. You will
the corresponding "simulated slit" overplotted on the "Spaxel Inspector". Clicking the right-mouse-button
will send the selected RSS-spectra to the Spectral Inspector. You can store this RSS-spectra by selecting
the "File->Save FITs Spectral Image" option.

You can visualize these spectra using:
spec2D_plot_all.pl spec_1.fits 0 66/xs

or

ds9 spec_1.fits



 2.3) Extraction of radial apertures or rings.

To perform radial extractions you can use the following commands:

radial_sum_cube.pl
USE: radial_sum_cube.pl CUBE.fits Delta_R X_C Y_C OUTPUT.RSS.fits [PLOT]

radial_sum_rss.pl
USE: radial_sum_rss.pl RSS.fits POS_TABLE.txt Delta_R X_C Y_C OUTPUT.RSS.fits [PLOT]


radial_sum_cube_ring.pl
USE: radial_mean_cube_ring.pl CUBE.fits Delta_R X_C Y_C OUTPUT.RSS.fits [PLOT]

radial_sum_rss_ring.pl
USE: radial_sum_rss.pl RSS.fits POS_TABLE.txt Delta_R X_C Y_C OUTPUT.RSS.fits [PLOT]

radial_mean_cube_ring.pl 
USE: radial_mean_cube_ring.pl CUBE.fits Delta_R X_C Y_C OUTPUT.RSS.fits [PLOT]

radial_mean_rss_ring.pl
USE: radial_sum_rss.pl RSS.fits POS_TABLE.txt Delta_R X_C Y_C OUTPUT.RSS.fits [PLOT]

All of them operates in based on the same principals, i.e., the integrate (sum)  or average (mean),
the spectra within a certain radius or in consecutive annulai (ring).

e.g., 

Not to visualize
radial_sum_cube_ring.pl NGC2916.V500.rscube.fits 3 35 34 rad_NGC2916.rss.fits 1

To visualize in the screen
radial_sum_cube_ring.pl NGC2916.V500.rscube.fits 3 35 34 rad_NGC2916.rss.fits 1

To have a hard-copy of the extracted areas:
radial_sum_cube_ring.pl NGC2916.V500.rscube.fits 3 35 34 rad_NGC2916.rss.fits 2
evince radial_mean_cube_ring_3.ps

This ring-extraction can be used once you have corrected for the velocity structure, to
extract annular rings of spectra and study the radial distributions of certain properties.

NOTE: To correct for the velocity map:

correct_vel_map_cube.pl
USE: correct_vel_map_cube.pl INCUBE.fits VEL_MAP.fits OUTCUBE.fits [FINAL_VEL]

NOTE: It is foreseen elliptical extractions, still not implemented in FIT3D.

2.4) Segmentation extraction.

The segmentation extraction has been instroduced in HIIexplorer. It allows to extract the
spectra from a considered datacube by using a segmentation map similar to the one created
by tools like SExtractor, where the value in each pixel represent the ID of the segmented
regions.

The segmentation maps should have the same X-Y size (NAXIS1,NAXIS2) of the considered cube.
Several segmentation maps can be considered, like those created using a voronoi binning.
We have implemented a segmentation scheme to detect HII regions (explained in Sanchez et al. 2012b),
and an example below, and a segmenetation map forseen to analyze the continuum emission
in a galaxy:

cont_seg_all.pl
USE: reg_spec.pl Ha_map.fits FLUX_LIMIT_PEAK MAX_DIST_PEAK FRAC_PEAK MIN_FLUX SEGFILE.fits DIFFUSE_MASK.fits


This latter segmentation map is created from a considered continuum image, and cluster regions
which peak emission is above a certain limit (FLUX_LIMIT_PEAK), less and a certain distance (MAX_DIST_PEAK), and
which relative flux is larger than a certain fraction of the peak intensity (FRAC_PEAK), all of them above
a limit intensity (MIN_FLUX). The output are a considered segmentation file and a gas for the
pixels below the considered limit (MIN_FLUX) 

cont_seg_all.pl map_V_NGC2916.fits 0.02 7 0.8 0.005 cont_seg.NGC2916.fits DMASK.NGC2916.fits

Explore the segmentation map:
ds9 map_V_NGC2916.fits cont_seg.NGC2916.fits DMASK.NGC2916.fits

We can extract now the corresponding RSS spectra:
spec_extract_cube_mean.pl NGC2916.V500.rscube.fits cont_seg.NGC2916.fits CS.NGC2916.RSS.fits

For CALIFA-data, you can also extract the error-spectra:
spec_extract_cube_error.pl NGC2916.V500.rscube.fits cont_seg.NGC2916.fits e_CS.NGC2916.RSS.fits

In addition to the RSS-spectra, you extract the corresponding PT (so you can visualize the
extracted spectra in E3D). The extracted spectra are the mean values within each considered segmented
region.

Sometimes it is interesting to create a cube using the segmented map and the RSS spectra:
rss_seg2cube.pl CS.NGC2916.RSS.fits cont_seg.NGC2916.fits SEG.cube.fits

NOTE: If you create a certain model of the SSP using a segmentation map, to subtract it to the original
cube, you need to rescale the model to the original cube intensity distribution. For using
so, it is recommended to scale using a certain band distribution.

rss_seg2cube.pl CS.NGC2916.RSS.fits cont_seg.NGC2916.fits SEG.cube.fits
get_slice.pl SEG.cube.fits SEG_img ../legacy//slice_V.conf
cp SEG_img_V_4500_5500.fits SEG.V.fits
clean_nan.pl  SEG.V.fits -1
imarith.pl NGC2916.V.fits '/' SEG.V.fits scale.seg.fits
#*** FIT with Multi-SSP model
imarith.pl SSP_mod.cube.fits '*' scale.seg.fits SSP_mod.NGC2916.cube.fits
imarith.pl NGC2916.V500.rscube.fits '-' SSP_mod.NGC2916.cube.fits GAS.NGC2916.cube.fits

More details are given in the "ana_single.NGC2916.log" file.

#
# Fitting emission lines
#

In this section we describe how to fit the emission lines once you
have subtracted the underlying continuum or assuming that the
underlying continuum it is not relevant for your analysis.

NOTE: A complementary information with exercises can be found in the
FTP repository, in the "tutorial.ps.gz" file (comprising the lessons
given in the 2005 IAC Winter School).

Basic routines:
fit_spec_back
fit_spec_back.pl
fit_spec_back_montecarlo.pl

GENERAL NOTE: All the entries in brakets are optional.

C-routine:

fit_spec_back
USE: fit_spec_back SPECTRUM.TXT CONFIG_FILE [BACKGROUND_LIST NBACK]

In this context "SPECTRUM.TXT" is a 1D spectrum, as described in before. The CONFIG_FILE
comprises the information of the model to fit, in this particular context the most
relevants are "eline" and "poly1d".

Any config file should start with a line with four entries. The first one
is keep unused and should be "0", the second defines the number of "functions"
that comprises the model, the 3rd one comprises the goal "Chi-Sq" and the
4th one the goal minimum variation of the "Chi-Sq" between iteration (ie., 
the converging criteria).

After that line, you sould include one row defining the function
(e.g., "eline" for the emission lines), and 9 rows, each one related
to the possible number or input parameters of the considered
fuction. For each parameter it is needed to define 5 entries: (1) The
input guess of the parameter, (2) a flag indicating if it is fitted
[1] or if it let fixed [0], (3) and (4) the min and maximum values
allowed for the variation of the parameter, and (5) a possible "link"
of this parameter with the same parameter coressponding to other
function.

If there is no "link" the entry (5) is set to "-1", if not, it
indicates the order of the function for which there is a link in the
considered configuration function.

There are two possible "links" considered in FIT3D: (a) additive links
and (b) multiplicative links. To define a link you should set the
entry (5) with the number of the function in the configuration file to
which your actual function is linked in the considered
parameter. Then, entry (4) is set to "0" (additive links) or "1"
(multiplicative links), and entry (3) becomes the value to be "added"
(case a) or "multiplied" (case b), to the linked function parameter to
define the parameter of the current function.

If a function has less than 9 parameters, they have to be included as:
0	 0	 0	 0	 -1
Until filling 9-rows.

For an "eline" the parameters are the following ones:

eline
CENTRAL_WAVELENGTH	 0	 0	 0	 -1
INTEGRATED_INTENSITY	 1	 -0.1	 1e10	 -1
Sigma_of_the_Gaussian	 0	 4.0	 4.5	 -1
Systemic_Velocity	 1	 4200	  7800	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1

For a "poly1d" the parameters are the coefficients of
a polynomial function (with the 1st one being a constant).


The program produces three output files:

i) out.fit_spectra

It is a file contaning the results of the fitting, for the different models in the
configuration file, order by their entry in that file, with one model per row.
The first entry is the number of individual functions included in the modelling:

Each row contains the name of the model, each parameter followed by its estimated errors.
Errors are set to 0 if the parameters are not fitted (set fixed in the configuration file).

eg.,
7 
eline  6562.8169 0.0000 219897.2934 866.6700 4.2135 0.0000 3.0858 0.0330 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
eline  6583.6001 0.0000 41876.9524 866.6700 4.2135 0.0000 3.0858 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
eline  6548.1001 0.0000 13958.8448 866.6700 4.2135 0.0000 3.0858 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
eline  6677.9702 0.0000 2447.2491 866.6700 4.2135 0.0000 3.0858 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
eline  6716.4702 0.0000 2944.0231 866.6700 4.2135 0.0000 3.0858 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
eline  6730.8501 0.0000 4404.5289 866.6700 4.2135 0.0000 3.0858 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
poly1d  592.3288 1593.2570 -61.3372 36.6783 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000

ii) out_mod_res.fit_spectra

A five columns ASCII table. The first column is the wavelength of the fitted data, the second the
original fitted spectra, the third is the model and the 4th is the residual (org-model). The
5th column is reserved.

ii) out_config.fit_spectra

A similar configuration file as the input one, but with the guess changed by the results of 
the fitting. Useful for iterative fitting processes.

#
# EXAMPLES: Fitting the emission lines on a single spectrum using
# "fit_spec_back.pl"
#

Go to the "tutorial" directory and run the program:

fit_spec_back.pl spec_single.txt emission_lines.txt 6400 6800 none 0 none

It will ask you which lines to fit. Select Ha+NII, and answer the corresponding questions:
Guess systemic velocity, guess e-line dispersion, guess intentity.

You can use the output "config" file (tmp.config) for future fits, with similar configurations.

Test the following command:
fit_spec_back.pl spec_single.txt emission_lines.txt 6400 6800 none 0 none Ha_NII_He_SII.config /xs


# Example of a config file defining the a single "system" 
# comprising Halpha+NII doublet and a continuum scale for 
# the background.

# In this configuration file all the emission lines have
# the same systemic velocities and the same width (note the links
# of entries 3 and 4 for 2nd and 3rd emission line).
# The 3rd emission line have the flux intensity fixed a 1/3
# of that of the 2nd line (i.e., the physical link between the NII-doublet lines).
#

0 4 1 0.5
eline
6562.82	 0	 0	 0	 -1
1000	 1	 -0.1	 1e10	 -1
4.1	 0	 4.0	 4.5	 -1
6270.0000     1	 4200	  7800	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
eline
6583.60	 0	 0	 0	 -1
1000	 1	 -0.1 1e10 -1
4.1	 1       0	 0	 1
6000      1	 0	 0	 1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
eline
6548	 0	 0	 0	 -1
1000	 1	 0.333   1 	 2	 
4.1       1	 0	 0	 1
6000      1	 0	 0	 1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
poly1d
0.1	 1	 -1e13	 1e13	 -1
0	 0	 -1e13	 1e13	 -1
0	 0 	 -1e13	 1e13	 -1
0	 0	 -1e13	 1e13	 -1
0	 0	 -1e13	 1e13	 -1
0	 0	 -1e13	 1e13	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1

NOTE: BACK_LIST and NBACK are deprecated.


Perl-Routine (requires the C one to work):

fit_spec_back.pl
USE:fit_spec_back.pl SPECTRUM.TXT LINE_FILE START_WAVELENGTH END_WAVELENGTH BACK_LIST NBACK MASK_LIST [CONFIG] [DEVICE]

"fit_spec_back.pl" is a wrapper of "fit_spec_back" that allows the
user to create the corresponding configuration files in a simple way
or just run a certain fit over a pre-defined configuration file.

The required entries are an ASCII spectrum (defined before), and file
including the relevant emission lines to fit, the starting and ending wavelength
of the considered range, two parameters that are deprecated (BACK_LIST and NBACK),
that should be set to "none" and "0" respectibely, and a mask_list.

The program will run in an interactive way to allow the user to construct the
previously described configuration file, by selecting the guess and range of
each parameter to fit, and the emission lines to including in the fitting process.
The output configuration file will be stored in the "tmp.config" file.

If you already have a pre-defined configuration file, you can use it as an optional
input.

The format of the input files are:

LINE_FILE:
5111.6299    [FeII]         
5158.7798    [FeII]         
5199        SKY_NI
5199.6001    [NI]           
5261.6201    [FeII]
5461        SKY_Hg         
5517.71      [ClIII]
5537.60      [ClIII]
5554.94      OI
5577       SKY_OI
5577.3101    [OI]
5685       SKY_Na   
...

For the mast list each row defines a wavelength range excluded from the fitting process,
in the following format:

MASK_LIST:
6520 6600
6700 6750
6670 6685
5868 5885
5565 5585
4820 5080
4460 4480
4320 4390

The optinal parameter [DEVICE] can be used to visualize and or plot the results of the fitting,
using the PGPLOT standard definion of devices (e.g., /XS for screen plots, file.ps/CPS for
PostScripts).

It is recomended that before fitting a set of spectra, or RSS, first you fit an individual
spectrum to define the configuration files. You can extract a single spectrum from a RSS
in different ways, the most simple is using the script:

>img2spec.pl
USE: img2spec.pl INPUT_FILE.FITS NY OUTPUTFILE.txt

You can use tools like PINGSoft or E3D to extract spectra of a particular region on the galaxy, 
both in the RSS format or individual spectra.


In a similar way it is possible to fit the miession lines of a 2D RSS data set, using the
script: 

kin_back_rss.pl
USE: kin_back_rss.pl rss.fits config_file back_file nback start_wavelength end_wavelength out_file DEVICE [MASK_FILE]

As an example, you can run:

kin_back_rss.pl crun31_01255b.obj.fits Ha_NII_He_SII.config none 0 6500 6760 kin_Ha_NII_He_SII.out /xs

The resulting file "kin_Ha_NII_He_SII.out" is a sequence of the consecutive outputs derived using 
"fit_spec_back.pl", for each individual spectrum:

cat kin_Ha_NII_He_SII.out

...
7 
eline  6562.8169 0.0000 115977.6520 8022.4024 4.2483 0.0000 -21.0699 0.0266 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
eline  6583.6001 0.0000 16711.2134 5144.1249 4.2483 0.0000 -21.0699 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
eline  6548.1001 0.0000 5570.3488 0.0000 4.2483 0.0000 -21.0699 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
eline  6677.9702 0.0000 1338.1127 5163.8604 4.2483 0.0000 -21.0699 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
eline  6716.4702 0.0000 673.0763 5152.3128 4.2483 0.0000 -21.0699 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
eline  6730.8501 0.0000 1256.8093 5161.6761 4.2483 0.0000 -21.0699 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
poly1d  183.0761 586.6401 -19.8004 13.3636 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
7 
eline  6562.8169 0.0000 64535.5779 6125.9096 4.3066 0.0000 -27.1695 0.0673 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
eline  6583.6001 0.0000 8144.4909 3925.7481 4.3066 0.0000 -27.1695 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
eline  6548.1001 0.0000 2714.8032 0.0000 4.3066 0.0000 -27.1695 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
eline  6677.9702 0.0000 626.2983 3941.9549 4.3066 0.0000 -27.1695 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
eline  6716.4702 0.0000 376.4248 3932.0440 4.3066 0.0000 -27.1695 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
eline  6730.8501 0.0000 598.7723 3939.2783 4.3066 0.0000 -27.1695 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
poly1d  796.9905 441.1477 -97.2168 10.0493 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
7 
...

These results are difficult to handle as thy are without the required
spatial information, which is provided in this case with the
corresponding position table "ppak_pt_arc.txt". You can use the two
following scripts to derive 2D maps and or slices:

map_interp_back.pl
USE: map_interp_back.pl slice.txt/POS_TABLE out.model PREFIX_OUT INTERP_MODE GRID_FUNC  GRID_PIX 

Where the three last parameters defines the interpolation procedure, as described
in the E3D User Guide. For this particular example we will adopt the Natural Neighbour
interpolation scheme:

map_interp_back.pl ppak_pt_arc.txt kin_Ha_NII_He_SII.out kin_Ha_NII_He_SII 3 1e-12 1.0

This script creates a map for each of the considered parameters and errors derived for
each emission line (flux intensity, systemic velocity and dispersion). It also provide
the corresponding slice for the flux intensity. The naming of each file is defined
by the "PREFIX_OUT" parameter, followed by an ordering number corresponding to the order
in the configuration file used as input in the "kin_back_rss.pl" file.

For comparison purposes it is possible to interpolate the flux slices with a different scheme,
like the one adopted by CALIFA (v 1.3c):

slice2map_int.pl 
USE: rss2cube_int.pl slice.txt dx OUTPUT.map.fits alpha BOX SCALE
type:
weight = exp(-0.5*(dist/(SCALE*dx))**alpha)

This script is forseen for dithered observations, and therefore its result for single pointings
is patchy:

slice2map_int.pl kin_Ha_NII_He_SII_flux_00.slice 1 kin_Ha_NII_He_SII_flux_00.slice.fits 2 5 1

As you can appreciate in:

ds9 kin_Ha_NII_He_SII_flux_00.fits kin_Ha_NII_He_SII_flux_00.slice.fits

Once derived the maps of the different parameters it is possible to derive high-order products,
like line ratios, analyze the kinematic structure, etc.

It is possible to visualize and work on the slice level, without requiring to perform an interpolation,
e.g.:
plot_slice.pl kin_Ha_NII_He_SII_flux_00.slice /xs C 2.68 -1 1e4 2e6 0.5 -0.5 califaCT 1

The need for one or the other depends on the science case to study.

In a similar way it is feasible to fit the emission lines in particular datacubes, using the
script:

kin_back_cube.pl
USE: kin_back_cube.pl cube.fits config_file back_file nback start_wavelength end_wavelength out_file DEVICE [MASK_FILE] [ASK]

The entries are basically the same as in the previous files, with the main difference that
it add an entry [ASK] in case that you want to run it in a completely automatic way.

We can test the result on three different datacubes included in the tutorial:
- NGC2916.V500.rscube.fits : The V500-setup CALIFA datacube, including the stellar continuum
- NGC2916.V500.zcube.fits  : Kinematic corrected version of the previous cube.
- GAS.NGC2916.cube.fits    : The stellar continuum subtracted version of the previous cube.
- GAS.NGC2916.zcube.fits   : The kinematic corrected (normalized to the systemic velocity of the galaxy)
  			      	version of the previous cube.



To fit this data, know with a different systemic velocity, first it is recommended to
extract a single spectra (e.g., the peak intensity one):
spec_cube_plot.pl NGC2916.V500.rscube.fits.gz 36 34 /xs

The spectra will stored in the temporary file "spec_cube_plot.out" (more spec_cube_plot.out),
that you can copy:
cp spec_cube_plot.out spec_NGC2916_peak.txt

Or you can get it directly:
get_spec_cube.pl NGC2916.V500.rscube.fits 36 34 spec_NGC2916_peak.txt
get_spec_cube.pl NGC2916.V500.rscube.fits[1] 36 34 spec_NGC2916_peak_e.txt

We will use this spectrum to determine the guess parameters for the fitting algorithm.

To fit just the ionized gas without subtracting the continuum it is a mistake, since
it instroduce systematic biases in all the emission lines, and in particular in the Balmer ones.

We extract another spectrum, with more emission line in Ha:

spec_cube_plot.pl NGC2916.V500.rscube.fits.gz 36 14 66/xs 
cp spec_cube_plot.out spec_NGC2916_HII.txt

fit_spec_back.pl spec_NGC2916_HII.txt none 6550 6750 none 0 none Ha_V500_NGC2916.config 66/xs

NOTE: Edit  Ha_V500_NGC2916.config, to determine how sensitive is to the parameters.

Now notice the different with an stellar subtracted version of the same spectrum:

spec_cube_plot.pl GAS.NGC2916.cube.fits.gz 36 14 66/xs 
cp spec_cube_plot.out spec_NGC2916_HII_gas.txt

fit_spec_back.pl spec_NGC2916_HII_gas.txt none 6550 6750 none 0 none Ha_V500_NGC2916.config 66/xs

NOTE: Do the same exercise for Hbeta.

Thus, even in the case of you are interested in the ionized gas, it is important to 
subtract the underlying stellar population.			

#
# Fitting stellar continuum
#

There are different method to fit and subtract the underlying stellar
population. FIT3D have tools to analyze the stellar population based
on three different methods: (1) Assuming a single stellar population
as a good representation of the stellar population, useful when
analyzing E-type galaxies. We used it to analyze the galaxies at the
core of Abell2218 (Sanchez et al. 2007); (2) Analyzing the absorption
indices (lick indices). With a limited use in the presence of emission
lines, and less precise than the SSP-analysis; (3) A multi-component
fitting of the underlying stellar population with a library of SSPs.
The procedure is described in Sanchez et al. (2011), in the analysis
of NGC628, and widely used in different articles.

#
# Fitting single stellar population
#

Fitting single stellar populations it is a pure approximation, considering
the complex SFH of galaxies. However, it is one of the most simple methods
to estimate the properties of the underlying stellar populations, and it may
be still valid. If the SSP-library included is large enough, in particular
if it includes constant starforming templates, it would procude a good
representation of the underlying stellar population.

However, for a precise understanding of the underlying stellar population 
this procedure is too crude (although some purists will find it more
realiable). 

The tools included in FIT3D to derive the best single stellar population
are (they do not fit neither the systemic velocities nor the velocity dispersion
of the stellar population):

compare_back_list.pl
compare_back_list_dust.pl
compare_back_list_3D.pl
compare_back_list_dust_3D.pl

The two first are useful to analyze single spectrum, the latter two are similar
to the previous ones but useful to analyze RSS spectra.

compare_back_list.pl
USE: compare_back_list.pl SPEC1.txt BACK_LIST OUTFILE MASK_LIST REDSHIFT SIGMA_DISP WAVE_SCALE PLOT [min max] [wmin wmax]

It requires a single spectrum, a list of SSP templates, the output file, a mask-table (similar to the one
describe before), the redshift and velocity dispersion of the stellar population, the wavelength region
at which the SSP and the input spectrum are normalized and a plot flag (1 to visualize the result and
0 not to).

compare_back_list_dust.pl
USE: compare_back_list.pl SPEC1.txt BACK_LIST OUTFILE MASK_LIST REDSHIFT SIGMA_DISP WAVE_SCALE PLOT Av_min Av_max delta_Av [min max] [wmin wmax]

All the entries are similar to the ones in the previous script, but adding a range of dust attenuations
to look for.

Example:

compare_back_list.pl spec_NGC2916_peak.txt list.SSP.miles compare_back_NGC2916.out ~/sda2/code/FIT3D_EXAMPLES/mask_elines.txt 0.0118 2.7 5500 1

The output of this analysis is a table including the scaling factor (in case that you want to estimate
the M/L ratio) and the chi-squre of the comparison. This output can be used to plot a Chi-square diagram.
A very crude example would be:

plot_comp_back.pl compare_back_NGC2916.out 66/xs R 0.01 -1 0.8 1.6 0.5 0.5 califaCT 0

For Single-SP, it is normally assumed that the Age and metallicity would correspond to the weighted (by
inverse of Chi^2) of the N (3-4) better fits. For the analysis of ETGs sometimes it is
sufficient, if the template library is large enough.

However, for more complex SFHs this procedure is not good enough.

#
# Fitting stellar continuum with miltiple SSPs
#
# (* Best option to analyze the stellar populations)

The main purpose of fitting the stellar population with multiple SSPs is to reconstruct
the SFHs of galaxies (or sections of galaxies), under the assumption that the stellar continuum
is the sum of different components, each of them corresponding to a particular burst of
starformation, and thus to a particular single stellar population, with its own Age and Metallicity.

Thus the stellar population is a "simple" linear combination of the SSPs comprised in a particular
library. However, in practice this is a bit more complex, since you should include all the possible
SSPs in the library (or you may assume that the inbetweeners are equal to the combination of 
some of them), you have to take into account the dust attenuation (that may or not be different
for each SSP), and the stellar kinematics (systemic velocity and dispersion).

In addition you have the ionized gas emission that may affect the result of the multi-SSP analysis,
since some of the lines affects the stronger features indicative of mostly the Age (e.g., Balmer
lines).

FIT3D takle this somehow complex issue by fitting the stellar
population to a combination of SSPs, linearly. In each step of the
linear fitting the library of SSP are convolved and redshifted to
match the systemic velocity and the stellar dispersion, and attenuated
by a certain dust attenuation (using the Cardielli, Clayton & Mathis,
1989, ApJ, 345, 245, extinction law). Once fitted and subtracted the
stellar continuum it is possible to fit (or not) a set of emission
lines, for which the minimization equation takes into accout the
result of the emission line subtraction. The procedure runs in an
iterative way, that effectively means that the stellar population and
the gas emission are fitted together.

The main scripts to address the multi-SSP isue are the following ones:

auto_ssp_elines_several_Av_log_new.pl
auto_ssp_elines_several_Av_log_rss_new.pl
auto_ssp_elines_several_Av_log_cube_new.pl

In the three cases the inputs are the same, although each of them handles 1D spectrum, 2D RSS spectra 
or 3D datacubes:

USE: auto_ssp.pl SPEC1.txt BACK_LIST.fits OUTFILE MASK_LIST CONFIG_FILE PLOT [min max] [wmin wmax] [redshift_elines_to_mask] [input_redshift delta_redshift min_redshift max_redshift] [input_sigma delta_sigma min_sigma max_sigma] [input_Av delta_Av min_Av max_Av] 
CONFIG_FILE:
redshift delta_redshift min_redshift max_redshift JUNK1 JUNK2 JUNK3 JUNK4 min_wavelength_kin max_wavelength_kin
sigma delta_sigma min_sigma max_sigma
Av delta_Av min_Av max_Av [Same range for all]
N_SYSTEMS
(1) START_W END_W MASK_FILE CONFIG_FILE NPOLY MASK_FILE_POLY
...
(N) START_W END_W MASK_FILE CONFIG_FILE NPOLY MASK_FILE_POLY
MIN_DELTA_CHISQ MAX_NITER CUT_MEDIAN_FLUX
start_w_peak end_w_peak
wavelength_to_norm width_AA new_back_templates.fits


The entries required are:
SPEC1.txt      => The spectrum/RSS/CUBE to be fitted.
BACK_LIST.fits => library template of SSP (described before)
OUTPUT         => Output file.
MASK_LIST => FIT3D standard mask table
CONFIG_FILE => The basic configuration file described before, with the required entries for the 
	    redshift, velocity dispersion and dust attenuation (Guess value, step in each variation,
	    min value, max value), and the entries required to configure the fitting of the emission 
	    lines, the number of systems to fit (or groups of emission lines), followed for the
	    same number of rows with the definitions required to run "fit_spec_back.pl".
	    The two last entries are deprecated.
	    In the first row there are 4 deprecated entries, and to additional ones defining the
	    wavelength range for which the stellar systemic velocity and dispersion are determined.
	    If they are not used, they will be redefined to the min-max wavelength range covered.
PLOT       => A flag that could be set to 0 (no visualization of the fitting), 1 (visualization
	   on the screen) or 2 (hard-copy PS of the fited results).

The remaining entries are optional, but the overide the ones in the configuration file (in this way 
you can use a general configuration file to be run over many different objects or individual spectra).


Lets explore the following example:
(NOTE: Update the directories corresponding to the location of the required files, normally at  FIT3D_EXAMPLES directory):

auto_ssp_elines_several_Av_log_new.pl  spec_NGC2916_peak.txt /disk-b/sanchez/ppak/legacy/miles_12.fits auto_ssp.NGC2916.peak.out /disk-b/sanchez/ppak/legacy//mask_elines.txt /disk-b/sanchez/ppak/legacy/auto_ssp_V500_several.config 1 -3 50 3800 6800 /disk-b/sanchez/ppak/legacy//emission_lines.txt  0.0124609539043173 0.000100069228559446 0.0114602616187229 0.0134616461899118 4.5 0.5 1.2 9.0 0.01 0.3 0.0 2.1

more ~/sda2/code/FIT3D_EXAMPLES/auto_ssp_V500_several.config
0.017 0.002 0.0001 0.027 10 100 0.1 0.5 3800 5500
3.2  0.0    1.9    6.5
0.4  0.1    0.0   2.5
5
6530 6950 mask_elines.txt /disk-b/sanchez/ppak/legacy/Ha_SII_V500.config 9 /disk-b/sanchez/ppak/legacy/mask_elines_H.txt 20 1950
4800 5200 mask_elines.txt /disk-b/sanchez/ppak/legacy/OIII_V500.config 9 /disk-b/sanchez/ppak/legacy/mask_elines_H.txt 20 1950
4310 4500 mask_elines.txt /disk-b/sanchez/ppak/legacy/Hg_V500.config 9 /disk-b/sanchez/ppak/legacy/mask_elines_H.txt 20 1950
4095 4250 mask_elines.txt /disk-b/sanchez/ppak/legacy/Hd_V500.config 9 /disk-b/sanchez/ppak/legacy/mask_elines_H.txt 20 1950
3700 3900 mask_elines.txt /disk-b/sanchez/ppak/legacy/OII_V500.config 9 /disk-b/sanchez/ppak/legacy/mask_elines_H.txt 20 1950
0.0001  1 0.05
6588 6760


The script performs a fitting of the underlying continuum with 12
stellar population from the miles-library, and fit a set of emission
lines distributed in 5 systems, from Ha+[NII]+[SII] to [OII]3727.

The 12 SSPs comprises a grid with 3 metallicities and 4 ages:
NAME0 | spec_ssp_00.0891_z0004.spec 
NAME1 | spec_ssp_00.0891_z019.spec 
NAME2 | spec_ssp_00.0891_z030.spec 
NAME3 | spec_ssp_00.4467_z0004.spec 
NAME4 | spec_ssp_00.4467_z019.spec 
NAME5 | spec_ssp_00.4467_z030.spec 
NAME6 | spec_ssp_01.0000_z0004.spec 
NAME7 | spec_ssp_01.0000_z019.spec 
NAME8 | spec_ssp_01.0000_z030.spec 
NAME9 | spec_ssp_17.7828_z0004.spec 
NAME10 | spec_ssp_17.7828_z019.spec 
NAME11 | spec_ssp_17.7828_z030.spec 

The different output results are the following:

1) The one including the information of the stellar population component:
more  auto_ssp.NGC2916.peak.out

That includes the chi^2, the luminosity weight ages and metallicities,
the average dust attenuation, the redshift and velocity
dispersion,median flux intensity averaged over the entire wavelength
range, the rms of the residual, the mass weigth ages and metallicities
(if the NORMXX keyword in the template fits file comprises the right
M/L ratio at the V-band), and the systemic velocity.

#MIN_CHISQ age_min met_min Av redshift sigma median_FLUX redshift_ssp med_flux rms_residual age_min_mass met_min_mass SYS_VEL /disk-b/sanchez/ppak/legacy/miles_12
.fits 
0.0289987607926226 7.33255645378024 0.0183418361081579 0 0.012350877752902 3.195 1564.29515267909 0.012350877752902 1.12432366609573 0.0246932238270997 7.33255645378024 0.0183418361081579 3702.69999999999

NOTE: More on the interpreation of these quatities will be explained somewhere else.

2) The individual weights of the multi-SSP decomposition, that in essence trace the SFH:

more coeffs.out
ID   AGE     MET    COEFF   Norm.Coeff  M/L   AV
0  0.0891  0.0004  0.1003  0.0976  1.00       0.00
1  0.0891  0.0190  0.0000  0.0000  1.00       0.00
2  0.0891  0.0300  0.0000  0.0000  1.00       0.00
3  0.4467  0.0004  0.0079  0.0077  1.00       0.00
4  0.4467  0.0190  0.0000  0.0000  1.00       0.00
5  0.4467  0.0300  0.0000  0.0000  1.00       0.00
6  1.0000  0.0004  0.0000  0.0000  1.00       0.00
7  1.0000  0.0190  0.0000  0.0000  1.00       0.00
8  1.0000  0.0300  0.1215  0.1183  1.00       0.00
9  17.7828 0.0004  0.0000  0.0000  1.00       0.00
10 17.7828 0.0190  0.7974  0.7764  1.00       0.00
11 17.7828 0.0300  0.0000  0.0000  1.00       0.00

NOTE: FIT3D is focused on the study of the gas emission, and therefore it provides with a minimalistic combination
of SSPs that reproduce the underlying continuum. Tests on the luminosity weight properties indicate that it provides
with similar results (e.g., Sanchez-Blazquez et al. 2013, MNRAS, submitted), but the SFHs may be different.

3) The results from the analysis of the emission lines:

more elines_auto_ssp.NGC2916.peak.out 
5 
eline  6562.6802 0.0000 2.2140 0.3808 4.0864 5.8030 3713.3118 12114.8421 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
eline  6583.4102 0.0000 2.1603 0.3808 4.0864 0.0000 3713.3118 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
eline  6548.0801 0.0000 0.7194 0.3808 4.0864 0.0000 3713.3118 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
eline  6730.7402 0.0000 0.7736 0.3389 3.6369 0.0000 3713.3118 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
eline  6716.3901 0.0000 0.8956 0.3389 3.6369 0.0000 3713.3118 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
3 
eline  5006.8398 0.0000 1.8992 0.1468 3.1087 0.7717 3666.2801 2767.8244 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
eline  4958.9102 0.0000 0.6324 0.1468 3.1087 0.0000 3666.2801 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
eline  4861.3198 0.0000 0.1474 0.1468 3.1087 0.0000 3666.2801 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
1 
eline  4340.4702 0.0000 0.6128 0.4295 5.0000 42.9482 3706.2533 205009.0817 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
1 
eline  4101.7402 0.0000 1.0052 0.3369 5.0000 45.8200 3802.7000 244903.0019 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
1 
eline  3727.3999 0.0000 3.9184 0.4260 3.5103 2.5858 3682.6239 16734.3478 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000

This file is compilation of the different outputs of "fit_spec_back.pl" run over the different systems defined before (look at the
explanation of this output in the descrition of fit_spec_back).

#
# Estimating the Errors.
#

The same program can be run generating a MonteCarlo on the basis of the residuals of each fitting, and providing with the corresponding estimation
of the errors.

auto_ssp_elines_several_Av_log_new_mc.pl


auto_ssp_elines_several_Av_log_new_mc.pl  spec_NGC2916_peak.txt /disk-b/sanchez/ppak/legacy/miles_12.fits auto_ssp_mc.NGC2916.peak.out /disk-b/sanchez/ppak/legacy//mask_elines.txt /disk-b/sanchez/ppak/legacy/auto_ssp_V500_several.config 1 -3 50 3800 6800 /disk-b/sanchez/ppak/legacy//emission_lines.txt  0.0124609539043173 0.000100069228559446 0.0114602616187229 0.0134616461899118 4.5 0.5 1.2 9.0 0.01 0.3 0.0 2.1

The MC simulation affects only to the determination of the errors in the emission lines.

#
# Analysing the physical quantities in the SSP analysis:
#

We can fit the continuum segmentated datacubes to analyze the multi-SSP components:

auto_ssp_elines_several_Av_log_rss_new_error.pl  CS.NGC2916.RSS.fits miles_12.fits auto_ssp.CS.NGC2916.rss.out mask_elines.txt auto_ssp_V500_several.config 0 -2 5 3850 6800 emission_lines.txt  0.0124609539043173 0.000100069228559446 0.0114602616187229 0.0134616461899118 4.5 0.5 1.2 9.0 0.3 0.3 0.0 2.1

Then we can reconstruct the 2D distribution of the different components.

map_auto_ssp.seg.pl elines_auto_ssp.CS.NGC2916.rss.out cont_seg.NGC2916.fits map.CS.NGC2916 

In this way we obtain the luminosity weighted distribution of the Age
(age_ssp), Metallicity (met_ssp), Dust extinction (Av_ssp); the Stellar Kinemacits (vel_ssp) and the SFH, via the 
weights of each SSP in the considered library (NORM_NN_ssp.fits files). 

It is possible to explore the 2D distribution of any of these properties:

ds9 map.CS.NGC2916_flux_ssp.fits map.CS.NGC2916_age_ssp.fits map.CS.NGC2916_met_ssp.fits map.CS.NGC2916_Av_ssp.fits map.CS.NGC2916_vel_ssp.fits &

Or the radial distributions:
plot_radial_pt_rss_no_error.pl CS.NGC2916.RSS.pt.txt auto_ssp.CS.NGC2916.rss.out 34 36 1 1 'R' 'Age (Gyr)' 66/xs
plot_radial_pt_rss_no_error.pl CS.NGC2916.RSS.pt.txt auto_ssp.CS.NGC2916.rss.out 34 36 1 2 'R' 'Z/H' 66/xs
plot_radial_pt_rss_no_error.pl CS.NGC2916.RSS.pt.txt auto_ssp.CS.NGC2916.rss.out 34 36 1 3 'R' 'Av (mag)' 66/xs
plot_radial_pt_rss_no_error.pl CS.NGC2916.RSS.pt.txt auto_ssp.CS.NGC2916.rss.out 34 36 1 12 'R' 'Km/s' 66/xs 0 50 3500 3950


Or to explore the SFH:

ds9 map.CS.NGC2916_NORM_?_ssp.fits map.CS.NGC2916_NORM_1?_ssp.fits &


NOTE: E.Perez will discuss more on the analysis of the SFH.

#
# Analysis of the Lick Indices
#

An alternative procedure to analyze the properties of the stellar component is to analyze the stellar absorption
indices or Lick indices. It is beyond the scope of this document to describe in detail the origin of a classical
method (there is a lot of examples in the literature). In FIT3D the routies to analyze the lick indices were
presented in Sanchez et al. 2007 (Analysis of the core of the Abell2218 cluster) and in Sanchez et al. 2011 (analysis
of the PINGS data on NGC628). We adopted the definitions of rmodel 
(http://pendientedemigracion.ucm.es/info/Astrof/software/rmodel/rmodel.html, by N.Cardiel),
and this program can be used to interpret the results of the indices.

For each index, the central band to measure it, and the blue and red bands to measure the adjacent continuum
are defined  in here:

__DATA__
Hd     4083.500 4122.250 4041.600 4079.750 4128.500 4161.000
Hb     4847.875 4876.625 4827.875 4847.875 4876.625 4891.625
Mgb    5160.125 5192.625 5142.625 5161.375 5191.375 5206.375 
Fe5270 5245.650 5285.650 5233.150 5248.150 5285.650 5318.150
Fe5335 5312.125 5352.125 5304.625 5315.875 5353.375 5363.375
D4000  4050.000 4250.000 3750.000 3950.000 0.000    1.000
Hdmod  4083.500 4122.250 4079     4083     4128.500 4161.000
Hg     4319.75  4363.50  4283.50  4319.75  4367.25  4419.75

The main routines to measure to indices are:

get_index.pl
get_index_rss.pl
get_index_cube.pl

The latter two assumes that the spectra have been corrected and shifted to a common redshift or systemic velocity.

The entries in all the cases are very similar
USE:get_index.pl SPECTRUM.TXT REDSHIFT DEV
USE:get_index_rss.pl SPEC.RSS.fits REDSHIFT DEV
USE:get_index_cube.pl NAME.CUBE.fits REDSHIFT DEV

E.g., 

get_index.pl spec_NGC2916_peak.txt 0.0118 66/xs
spec_NGC2916_peak.txt Hd  -1.08420139455479 Hb  1.14047026815434 Mgb  3.26772162193026 Fe5270  2.14234558475139 Fe5335  1.47296864775813 D4000  1.6584157249783 Hdmod  -1.5040249834666 Hg  -5.57578775655551 


However, in most of the cases the 3D data have its original kinematics. Then you can use:
get_index_cube_velmap.pl
USE:get_index_cube_velmap.pl NAME.CUBE.fits vel_map.fits DEV [OUTPUT]

If you have a good "noise" spectrum, it is possible to perform a Montecarlo and derive the estimated errors of the corresponding parameters:
get_index_montecarlo.pl
get_index_rss_montecarlo.pl

USE:get_index_montecarlo.pl SPEC.txt noise.txt NSIM REDSHIFT DEV

E.g., 

get_index_montecarlo.pl spec_NGC2916_peak.txt spec_NGC2916_peak_e.txt 10 0.0118 66/xs
Hd -1.16358202474294 0.305169365673897 Hb 1.10918573087121 0.182358378302967 Mgb 3.31905278243715 0.135234413451691 Fe5270 2.17964772122336 0.122190218209447 Fe5335 1.43589545571944 0.22501789136417 D4000 1.65503773347337 0.0125152412583075 Hdmod -1.53533825536561 0.74120842633777 Hg -5.46116078308321 0.5013933465301 SN  1.16672539710999  0.0317976689702384


NOTE: No similar tool has been created for the analysis of cubes, but it can be easily adopted.

If you have analyzed the RSS spectra using:
auto_ssp_elines_several_Av_log_rss_new_error.pl 
or
auto_ssp_elines_several_Av_log_rss_new.pl

You can use the output of the fit as the entry of the analysis of the indices, to take into acount the derived kinematics: 

get_index_output_auto_ssp_elines_several_Av_log_rss_outFIT3D.pl

e.g.,

get_index_output_auto_ssp_elines_several_Av_log_rss_outFIT3D.pl output.auto_ssp.CS.NGC2916.rss.out.rss.fits 5 auto_ssp.CS.NGC2916.rss.out /null > indices.CS.NGC2916.rss.out

You can visualize the radial distributions:

plot_radial_pt_rss.pl CS.NGC2916.RSS.pt.txt indices.CS.NGC2916.rss.out 34 36 1 7 'R (arcsec)' 'EW(Mgb) \A' 66/xs
plot_radial_pt_rss.pl CS.NGC2916.RSS.pt.txt indices.CS.NGC2916.rss.out 34 36 1 16 'R (arcsec)' 'D4000 \A' 66/xs 0 50 0 2

Or compare among them:

table_plot_e.pl indices.CS.NGC2916.rss.out 16 7 'D4000 \A' 'EW(Mgb) \A' 66/xs 0 2 -2 11
table_plot_e.pl indices.CS.NGC2916.rss.out 16 4 'D4000 \A' 'EW(H\gb) \A' 66/xs 0.1 2.1 -1 7

To interpret the underlying stellar population on the basis of the indices, it is required to compare with
a template library grid (e.g., like rmodel does).

#
# Fitting 2D distributions (e.g., velocity maps)
# TBW

FIT3D include some tools to fit 2D distributions (slices), or 2D maps. In principal these tools are thought to 
analyze the properties of velocity maps:


fit_2D_map vel_map.txt fit_vel.config

0 1 0.00000002 0.0000000005
rot_cur
0        0       -5     5       -1
0        0       -5     5        -1
387.378  1       64.563 774.756  -1
0        0       -20    20       -1
46.81692     1   44.81692        48.81692        -1
-25.8418723336285   0       -180         180     -1
0.5      1       0.1     10      -1
0        0       0.1     5.49733705340228 -1
0        0       0       0       -1


plot_slice.pl fit_2D_map.mod 55/xs C 1.34

plot_maps_vel_cont_CALIFA.pl ve.NGC2916.vel_map.fits.gz -300 300 0.5 0.5 'NGC2916' 1 66/xs 3750 map_Ha.NGC2916.fits 10 0.05 0.5 0

fit_2D_imag.pl ve.NGC2916.vel_map.fits.gz fit_vel_map.config out_fit2D_imag.fits 66/xs

fit_2D_imag.pl map.CS.NGC2916_vel_ssp.fits fit_vel_map.config out_fit2D_imag.fits 66/xs

#
# kinematrics:
#

table_plot_2_csv.pl kinematics.out 2 4 2 8 'Dist (arcsec)' 'Vel (km/s)' 66/xs 0 80 -300 300
kinematrics.pl ve.NGC2916.vel_map.fits.gz -250 160 0.5 0.5 'NGC2916' 1 66/xs 3730 map_Ha.NGC2916.fits 10 0.2 0.5 0



#
# HIIExplorer
# 

Analysing the properties of the ionized gas, HII regions and diffuse gas.
HIIexplorer provides with tools to identify and segment the HII regions
in a particular datacube. The detail explanation is in Sanchez et al. 2012b,
and the HIIExplorer webpage:
http://www.caha.es/sanchez/HII_explorer/index.html

We will concentrate here in the procedures to analyze the HII-regions in a galaxy.

First we extract a Halpha map (if you have a velocity corrected datacube it will be more accurante).
This procedure extract a narrow-band image of 30AA with around the expected wavelength
of Halpha at the considered redshift:

get_Ha.pl NGC2916.V500.rscube.fits 0.0124609539043173 map_Ha.NGC2916.fits

NOTE: You can provided with any other different Halpha (or other emission line) as input.

Once derived this narrow-band image, we proceed to create the HII-regions segmentation map (see before):
HII_recover.pl map_Ha.NGC2916.fits 0.4 3.5 0.1 0.1  seg_Ha.NGC2916.fits mask_Ha.NGC2916.fits > HII_recover.NGC2916.log

Check the result:
ds9 map_Ha.NGC2916.fits seg_Ha.NGC2916.fits mask_Ha.NGC2916.fits

We extract the corresponding RSS and position table for the HII regions (and the diffuse gas content),
for doing so it is better to use the velocity corrected cube:

spec_extract_cube.pl NGC2916.V500.zcube.fits seg_Ha.NGC2916.fits HII.NGC2916.rss.fits

NOTE: The determination of the stellar kinematics fails in the presence of strong emission lines,
and young stellar populations. The reason for that is that there is a lack of spectral features.
Therefore for HII it can affect the result. It is more convinient to use the stelar

spec_extract_cube_error.pl NGC2916.V500.rscube.fits seg_Ha.NGC2916.fits e_HII.NGC2916.rss.fits


spec_extract_diffuse.pl  NGC2916.V500.rscube.fits  map_Ha.NGC2916.fits diffuse.fits diff.NGC2916.rss.fits 0.01
cp seg.diffuse.fits seg_diff.NGC2916.fits

Then we analyze the HII-regions, using FIT3D, fitting and subtracting the underlying stellar population, and
deriving the emission line intensity for a set of emission lines:

rm -f elines_auto_ssp.HII.NGC2916.rss.out
auto_ssp_elines_several_Av_log_rss_new_error.pl HII.NGC2916.rss.fits miles_12.fits auto_ssp.HII.NGC2916.rss.out mask_elines.txt auto_ssp_V500_several.config 0 -2 5 3850 6800 emission_lines.txt 0.0124609539043173 1.66782047599076e-05 0.0123942410852777 0.012527666723357 4.5 0.5 1.2 9.0 0.5 0.2 0.0 1.6


Then we extract the corresponcing HII slices:

map_interp_back_Av_log_new.pl HII.NGC2916.rss.pt.txt elines_auto_ssp.HII.NGC2916.rss.out FIT.HII.NGC2916 0 0 0

An example of the possible application is to analyse the BTP diagram of the inonized regions:
ls FIT.HII.NGC2916*

NOTE: The last entry in this tables correspond to the all the diffuse gas together, and this has not to be taking into account
for any analysis on the radial distribution of properties.

plot_diag_NAME.pl NGC2916 66/xs 36 34 

plot_slice_out_fit.pl FIT.HII.NGC2916_FLUX.slice 66/xs C 2 3 0 400

We can plot the radial intensity of Halpha.
plot_radial_pt_rss.pl HII.NGC2916.rss.pt.txt FIT.HII.NGC2916_FLUX.slice 34 36 1 3 'R' 'Ha' 66/xs 0 30 -1 400
plot_radial_pt_rss_log.pl HII.NGC2916.rss.pt.txt FIT.HII.NGC2916_FLUX.slice 34 36 1 3 'log(R)' 'Ha' 66/xs 0 30 -1 3

We can plot the radial distribution of a certain emission line.
plot_radial_pt_rss_ratio.pl HII.NGC2916.rss.pt.txt FIT.HII.NGC2916_FLUX.slice 34 36 5 3 'R' '[NII]/Ha' 66/xs
plot_radial_pt_rss_ratio_log.pl HII.NGC2916.rss.pt.txt FIT.HII.NGC2916_FLUX.slice 34 36 5 3 'R (argsec)' 'log([NII]/H\ga)' 66/xs 0 30 -1 1

NOTE: Try to create a gradient based on O3N2

./plot_radial_OH_O3N2.pl HII.NGC2916.rss.pt.txt FIT.HII.NGC2916_FLUX.slice 34 36 5 3 'R (argsec)' '12+log(O/H)' 66/xs 0 30


NOTE: Try to create a diag-map spaxel-by-spaxel



#########################################################################################
# Detailed information in some FIT3D routines
#########################################################################################


Fitting tool.

USE:
velocity
fit_spectra
fit_spec_back
sky_sub (experimental!!)
fit_2D_map


SCRIPTS:
fit_spectra.pl
velocity.pl (deprecated)
map.pl (deprecated)
fit_3D_fg.pl
fit_3D_force.pl 
fit_3D_force_delta.pl 
poly_force.pl 
mean_force.pl 
extract_spec.pl
auto_ssp
auto_ssp_elines_several_Av_log.pl


FIT_SPECTRA (09.02.04)
fit_spec_back (27.09.05)
-----------------------

'fit_spectra' fits a system over a single spectrum.
It requires a config file in the form:

-------------------------
0 2 0.2 0.001
gauss1d
5182.245	 1	 5176.145	 5188.345	 -1
0.4	 1	 0.3	 0.7	 -1
4.6	 1	 4	 7	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
gauss1d
5182.245	 1	 -49.6800000000003	 0	 1
0.13	 1	 0.33	 1	 1
4.6	 1	 0	 0	 1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
eline
4862	 1	 0	 0	 1
0.13	 1	 0.33	 1	 1
4.6	 1	 0	 0	 1
veloc	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
0	 0	 0	 0	 -1
------------------------------------
Where the first 4 numbers are:
1st = Not used (by now)
2nd = number of models to fit
3th = Limit Chi^2
4th = Limit delta Chi^2
Them for each model we have 9 parameters:
NAME-OF-THE-MODEL (gauss1d, poly1d)
CENTRAL_WAVELENGTH   FIT(0/1) MIN  MAX FREE_FLAG
INTESITY	     FIT(0/1) MIN  MAX FREE_FLAG
WIDTH(SIGMA)	     FIT(0/1) MIN  MAX FREE_FLAG
... and 9 more!

eline: 1d gaussian:
CENTRAL_WAVELENGTH   FIT(0/1) MIN  MAX FREE_FLAG
INTESITY	     FIT(0/1) MIN  MAX FREE_FLAG
WIDTH(SIGMA)	     FIT(0/1) MIN  MAX FREE_FLAG
VELOCITY             FIT(0/1) MIN  MAX FREE_FLAG

back: a background input file!
NBACK                FIT(0/1) MIN  MAX FREE_FLAG
INTESITY	     FIT(0/1) MIN  MAX FREE_FLAG
WIDTH(SIGMA)	     FIT(0/1) MIN  MAX FREE_FLAG
VELOCITY             FIT(0/1) MIN  MAX FREE_FLAG

NBACK= Number of file from the list of files

The FREE_FLAG is to link a line with another. 
There are two ways to link: ADD or MULTIPLY a parameter.
In that case we have:
PARAMETER  1 DELTA TYPE(0/1) LINE_TO_LINK

Where LINE_TO_LINK is the line that this one is associated.
TYPE is the kind of link (0=ADD, 1=MULTIPLY).
DELTA is the parameter to ADD or MULTIPLY.

The input spectrum should be in the format:
452 5003 -0.00526132
453 5006 0.000298891
454 5009 -0.0123843
455 5012 0.00534711
456 5015 0.000978749
457 5018 -0.0046241
458 5021 0.00411005
459 5024 0.00710419
460 5027 0.0213663
461 5030 0.0156385

I.e.:
ID  WAVELENGTH FLUX
(It is assumed that SIGMA=SQRT(ABS(FLUX)), on the fit.

It is possible to run a script that helps you to generate the
config file:
"fit_spectra.pl"

For this script you need an emission line file in the format:

3727    [OII]         
5007    [OIII]        
4959    [OIII]        
4861    H\beta        
4340    H\gamma       
5876    HeI           
6548    [NII]         
6563    H\alpha       
6583    [NII]         
6716    [SII]         
6731    [SII]         


fit_spec_back 
include fittings over anykind of background, consisting of
a set of spectra ASCII files, in the format:
id wave flux
-------------------------------
back
1	 0	 0       0       -1
1.0	 1	 0.1      100	  -1
0	 0	 0	  0		  -1
31400	 1	 20000    40000	 -1
0	 0	 0	  0	   -1
0	 0	 0	 0	    -1
0	 0       0	 0	     -1

first value: wich background file to use.
2nd value:  Intensity (multiplicative factor).
3rd value: Velocity of this system.



KINEMATICS.PL
KINEMATICS_RSS.PL -> like kinematics but it does not requires
		  the E3D installation.
---------------
This program uses the FIT_SPECTRA, as a loop, over a certain 
3D FITs file.

kinematics.pl E3D_FITSFILE LINE_FILE START_WAVELENGTH END_WAVELENGTH OUT_FILE

It requires as input the E3D_FILE, the emission-line file, the
star_wavelenght and the end_wavelength of the fitting, and 
the output file.

It produces, appart from the output-file (containing the results from
the fitting) a slice cut (slice.kin), with the needed spatial
information to create a map.


MAP.PL
MAPGRID.PL (only for previosly gridded DATA).
--------
This program creates ASCII MAPS, using the results from the
KINEMATIC.PL.

map.pl slice.txt out.model PREFIX_OUT

For each MODEL included in the "out.model", it creates a file named:
PREFIX_OUT_NN.out
Where NN, is the model number.

Each of this files contains:
ID X Y I e_I wave e_wave sig_gaus e_sig_gaus Quality_flag

E.g.,
6 2.1 7.2 0.5083 0.0508 5185.6140 3.6017 4.9823 0.8095 1
7 2.4 7.2 0.5158 0.0597 5185.0353 3.2409 4.3683 0.7752 1
8 2.7 7.2 0.4471 0.0665 5183.1956 4.6284 4.0000 2.4006 1
9 3 7.2 0.4696 0.0717 5183.0539 5.4320 4.0203 3.5377 1
10 3.3 7.2 0.3824 0.0580 5182.0009 5.5433 4.0000 3.0033 1
11 1.8 7.5 0.6162 0.0793 5186.1130 3.3052 4.4353 1.5653 1
12 2.1 7.5 0.7002 0.0915 5185.9623 2.9712 4.4942 1.5909 1

This maps can be imported by E3D.

TK_E3D.TCL
-------------
Once you have a velocity map, you can upload it to E3D, using
"read_vel_map", and the corresponding velocity map!

e.g.,
read_vel_map test_00.map


MODEL and RESIDUAL
------------------
For having the model and the residual you need the PERL CFITSIO
Module, and use:

ascii2fits.pl res.kinematics res_kin.fits

ascii2fits_rot.pl res.kinematics res_kin.fits

So you have fits file with the RESIDUAL spectra.


**************************************************************
FIT3D:
----------

 3D spectroscopy data can be understood a sequence of images taken at
 different wavelengths. In this case each of these images can be modelled,
 by 2D functions (e.g., Gaussians and Moffat functions for point-like sources,
 and/or Freeman functions, deVaucouleurs functions or Sersic functions for
 extended objects) extracting the integrated flux asociated with the different
 components of the model. Repeating this fitting over the sequence of images
 at different wavelengths it is possible to extract the spectra of the
 different components in the objects. There are several packages that model 2D
 images, like GALFIT or GIMP, or techniques to decouple point-like and
 extended objects in 2D images (2 channel fitting). 

 Different experiments have been done in order to disentangle different
 spectra in a multicomponent-IFS: Becker et al. 2004 (AN, ref?), Wisotzki et
 al. 2003 (A&A), Sanchez et al. 2004 (ApJ), Garcia-Lorenzo et al. (2004).
 These experiments were done addapting procedures used to decouple 2D data to
 IFS (e.g., Becker et al. 2004, Sanchez et al. 2004), or using code created
 ad-hoc for an expecific problem (Wisotski et al. 2003). In all these cases
 the input data were regular gridded IFS. However, in a general situation, 3D
 data are not regular gridded, and more than a sequence of 2D images has to be
 understood as a sequence of 2D distributions of data. No expecific software
 is actually available to disentagle the spectra of a multi-component object
 using IFS in this general situation.

 FIT3D is an intent to solve this problem. FIT3D is a package of tools to fit
 3D data using 2D models, extracting the spectra associated with each of these
 models, creating a 3D model and a residual datacube. The input data does not
 need to be in a regular grid, and therefore, FIT3D handles the most general
 situation of IFS data.

 The different packages/tools that comprises FIT3D are described below (*
 schematic summary * ):

FIT_2D:
'''''''

fit_2D_map SLICE CONFIG
-------------------------
This program fits an ascii 2D table,
to the model described (in the standard way) 
by CONFIG.
The ascii table should be in the format:
ID X Y FLUX

(note: sigma=sqrt(flux) -> latter modifications?).

fit_3D_fg.pl E3D_FITSFILE CONFIG_FILE OUT_FILE
---------------------------
This script goes throughout the datacube
and use fit_2D_map, to fit, for each MAP, 
a model described by CONFIG_FILE.

The output is storaged (as a collection of models), 
at OUT_FILE.

It creates a ORG, MODEL and RES ascii file,
that should be transformed to 2D fitsfiles by
the script:

To create the 2D images:
ascii2fits_rot.pl res.fit_3D_fg res.fit_3D_fg.fits 

fit_3D_force.pl E3D_FITSFILE CONFIG_FILE INPUT_FILE OUT_FILE
fit_3D_force_delta.pl E3D_FITSFILE CONFIG_FILE INPUT_FILE OUT_FILE DELTA BIAS
------------------------------------------------------------
Similar to fit_3D_fg.pl, but it uses as initical parameters,
for each fit, the values storaged in a previous fitting.
It can be used to force the fitting fixing parameter, once
you have run a "first-guess" fitting.

The INPUT_FILE is the output from fit_3D_fg.pl or the output
of this scripts itself.

DELTA = Width of the polichromatic image in spectral pixels
BIAS = Level to add to the noise (Sigma=sqrt(sigma^2+bias))

poly_force.pl PARAM MODEL POLY_ORDER N_MIN N_MAX INPUT_FILE OUT_FILE
--------------------------------------------------------------------
This scripts read the OUTPUT from the 2 previous scripts, and 
fits a polynomical function of order POLY_ORDER to the parameter
PARAM of the model MODEL, in the array range N_MIX,N_MAX.

The result is storaged in OUT_FILE.
This OUT_FILE can be use as input for "fit_3D_force.pl" for producing
a force fitting.

mean_force.pl PARAM MODEL N_SIGMA BOX_FWHM INPUT_FILE OUT_FILE
--------------------------------------------------------------
Does a mean filtering over the data.














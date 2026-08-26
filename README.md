# Public_Synthetic_RSO_CCD_Image_Generator
This repository contains software originally developed for the Geo-Informatics and Space Technology Agency of Thailand (GISTDA) within an internal repository. Scripts, example files, and resources are copied here for reporting purposes.

# garnet_image_simulation

## Description

This image_generator.py can be used to generate simulated images of RSOs with background stars given an inputted TLE and information about the observer, observation time, exposure time, telescope details, and other inputs listed in the **Arguments** section. It is designed for use generating simulated telescope images of RSOs for use training detection algorithms. This generator creates images that simulate excellent seeing conditions (Where the sky is not very bright compared to stars and RSOs), and while some extra noise can be added the star fluxes currently remain as brihgt as possible for their magnitudes. The affect of dust/scratches on lens/mirrors (things that can often be corrected with flat frames) is also not included in these simulated images. 

A .ipynb file used to test the image generator code is included in the examples folder, and examples of generated images and data files are included in the resource folder. Example TLEs to use in image generation are also included in the resource folder, and the format of the single or grouped TLEs should be followed for use with image_generator.py. 

Required modules for this image generator are included in the requirements.txt file. To install:
```
python -m pip install -r "requirements.txt"
```

## Usage

**Prefered Use:** To generate a larger amount of varied images, run the second code example below with your output file, TLE file, and chosen number of images. The simplest use of the script can generate 1 set of files using a TLE of the ISS already encoded, as shown in the first example. 

Pass TLE files and observing information to generate images of other objects, other epochs, images from other observing locations, images with other telescopes or cameras, ect. Arguments that can be passed are listed in the **Arguments** section of this README. 

The script will as a default generate one set of TRACKING mode image data and one set of LEAPFROG mode image data (one set = one .fits file, one .png file, and one labelme .json file). A mode can be specified while running the script to produce only images of that mode. Here "TRACKING" refers to satellite-rate tracking images while "LEAPFROG" refers to sidereal-rate tracking images. These terms are used internally at S-TREC interchangeably. 

**Example Code:**

**1.** Running image_generator in the simplest use, which generates an images of the ISS: 

```
python image_generator.py "path/to/output/file/"
```

**2.** Running image_generator with specified TLEs and a number of images to create. A mode is not specified in this example, therefore for each TLE, 100 TRACKING and 100 LEAPFROG images will be produced. The total amount of images produced will be 200 * (number of TLEs in TLEs.txt): 

```
python image_generator.py "path/to/output/file/" -t "TLEs/example_TLEs.txt" -n 100
```

**3.** Running with TLE file and number of images as well as exposure time, observer, and initial observation date: 

```
python image_generator.py "path/to/output/file/" -t "TLEfile.txt", -w "2026-05-02 23:50:50.4663" -e 0.5 -ob "SKP" -n 100
```

## Arguments


The output file path must always be specifed when running image_generator.py as the first argument. Following that, many optional arguments may be passed to alter the configuration of image(s) produced. All optional arguments accepted by image_generator.py are shown below. 


| Optional Parameter Command | Parameter | Use | Default Value |
| --- | --- | --- | --- |
| -t  | TLE file | Specifies TLE or list of TLEs to use for generation | None (will use preset ISS TLE) |
| -n  | Number of images | Specifies number of images to create for each mode and for each TLE (ex: with 2 TLEs and 2 modes, n = 1 will produce 4 total images, for 1 TLE and TRACKING mode only, n = 1 will produce 1 image) | 1 |
| -m  | Mode Generated | Specifies 'TRACKING' or 'LEAPFROG' mode to generate | None (will create images in both modes) |
| -e  | Exposure Time | Specifies image exposure time in seconds | None (will randomize between 0.3-5 sec) |
| -w  | When Observed | Specifies start of observing run | None (will use next pass after current time) |
| -ob  | Observer | Specifies observer name | None (Will default to SKP settings unless other values input) |
| -lat | Observer Latitude | Specifies observer latitude | None (will use Observer value or SKP value) |
| -lon | Observer Longitude | Specifies observer longitude | None (will use Observer value or SKP value) |
| -alt | Observer Altitude | Specifies observer elevation | None (will use Observer value or SKP value) |
| -rm | RSO magnitude | Specifies RSO magnitude | None (will randomize between 5 and 12) |
| -ml | Magnitude Limmit | Specifies magnitude limit for background star selection | 16 |
| -p | Pixel Size | Specifies size of camera pixel in microns | None (will use Observer value or SKP value) |
| -nc | Number Columns | Specifies number of columns in image array | None (will use Observer value or SKP value) |
| -nr | Number Rows | Specifies number of rows in image array | None (will use Observer value or SKP value) |
| -f | Focal Length | Specifies telescope focal length in mm | None (will use Observer value or SKP value) |
| -bin | Binning Factor | Specifies binning factor for image array | 2 |
| -ofl | Offline | Indicates use of offline catalog | False (Using preferred online catalog) |
| -s | Include Sky | Indicates inclusion of noise from sky brightness | False |
| -b | Incude Bias | Indicates inclusion of bias collumns | False |
| -c | Catalog | Specifies which offline catalog you want to use: Currently supports Gaia ("G"), Hipparchos/Tycho ("H"), or combined NOMAD ("N") | None |
| -cf | Catalog File | Specifies the catalog file to use (necessary for NOMAD or Hipparchos use) | None |


**Note on observer name and location information:** To specify a different observer, use [-lat], [-lon], [-alt] arguments to specify observer latitude, longitude, and altitude. To use a preset obaerver, use [-ob "SKP"] or [-ob "DSC"] to use the observer and telescope information of the SKP or DSC sites. Specifying an observer by string will by defult set the lat, lon, and alt as well as pixel size and focal length, however the preset values will be overwritten if an argument for any is also passed along with the observer string. Ex: [-ob "SKP", -alt 100.0] will use all default SKP site values with an altitude of 100 m instead of the preset 49.8. 


## Note About Astroquery Use

This program uses the astroquery python library to find star information unless offline catalog use is specified. Online astroquery is the default and preferred option as it searches for stars in the large GAIA database. However, a warning message saying the archive is "in evolution" may pop up when image_generator.py is run since the GAIA archive is currently being worked on in preperation for the DR4 release, but can be ignored. Any other astroquery warning messges should not be ignored. An Info message will print out in terminal when the star query is completed.

This generator uses an online asyncronous search of the GAIA database (using astroquery) by defualt. By specifying offline catalog use, the script can run with either an offline gaia catalog (which takes 1+ days to download: https://pypi.org/project/gaiaoffline/), a downloaded NOMAD combined catalog which takes a similar amount of time to download, or a downloaded Hipparchos catalog (https://cdsarc.cds.unistra.fr/ftp/cats/I/239/) which takes minutes to download. The Hipparchos catalog, while able to work offline, contains fewer stars than GAIA so uing GAIA is prefered. Using Hipparchos will result in ~5 stars in an image while using GAIA will result in ~100 stars in an image (and can do more). 


**gaiaoffline downloading instructions:** For information on downloading the full gaiaoffline catalog, see the gaia_offline_creator notebook in the examples folder or go to this webpage: https://pypi.org/project/gaiaoffline/#description. 


## Authors and acknowledgment

Authors: Bruce Ritter


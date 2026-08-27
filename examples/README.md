# Examples


## example.ipynb

The "example" jupyter notebook includes code to test out configurations of image_generator.py.  


## rotate_images.py

Images produced by the image generator will not match the orientation of real images. This was determined not to matter for the purposes of detection algorithm training, however if you wish to verify that the RSO location derived from your input TLE is accurate, this image rotator code can be useful. To check the accuracy of image generation from TLEs, follow these steps: 

1. Find a real image with visible RSO and background stars
2. Find the TLE for the RSO in the image from the nearest time before it was taken
3. Run image_generator with the TLE, with observer location, telescope details, observation time, and exposure time matching the details of the real image. The example code below will generate an image matching a real image taken from the DSC site included in the resource/example_results folder. 

```
python image_generator.py "resource/example_results/" -t "TLEs/NORAD_67683.txt" -ob "DSC" -nr 3194 -nc 4784 -p 3.76 -f 1050.0 -e 0.3 -m 'TRACKING' -ofl True
``` 

4. Run rotate_images.py with the generated TRACKING mode image and a rotation angle that causes the angle of the star streaks to match that of the real image's stars (An angel of ~ 85 degrees works for the above example) and visually compare the patterns of stars to the real image. The RSO position may not exactly match its position in the real image, however they should be close enough that a user can see a matching pattern of background stars.  

## gaia_offline_creator.ipynb

This notebook can be used to create and fill a local database with the full Gaia DR3 catalog of stars. 
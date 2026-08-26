"""
image_tools.py

Version date: 19/6/2026
Author: Bruce Ritter

Contains ImageTools class for use with image_generator.py, a script
written to generate simulated telescope images of RSOs.

List of Functions:
- vignette
- add_noise
- add_psf_object
- line_points_difference
- add_blurred_line
- make_labelme_shape
"""
import numpy as np
from PIL import Image, ImageDraw
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.wcs import WCS
from astropy.io import fits
from astroquery.gaia import Gaia


class ImageTools:

    @staticmethod
    def vignette(data):
        """Function to add a vignette affect (sky brightness affets) to the data array
        Currently this does so with a gaussian distribution from the center, but in futer versions
        it could probably be improved to be slightly off-center or have different gradients.
        Inputs: 
        - data array
        Output: 
        - array of values to add to data to create vignette affect
        """
        nrows, ncols = data.shape       
        flux = np.mean(data)  # aproximating some amount of flux difference at corners, was /10 now is not

        cx, cy = ncols // 2, nrows // 2   # cx = col centre, cy = row centre
        fwhm = nrows
        sigma  = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        rows_idx, cols_idx = np.ogrid[0:nrows, 0:ncols]
        dist_sq = (rows_idx - cy)**2 + (cols_idx - cx)**2
        
        return flux * np.exp(-dist_sq / (2 * sigma**2))

    @classmethod
    def add_noise(cls, data, incl_sky, incl_bias):
        """Function combining all the noise adding functions
        Inputs: 
        - data array
        - incl_bias: boolean to include bias collumns or not (can be removed irl with bias frames)
        - incl_sky: boolean to include sky brightness/telescope aperture affects or not (can be removed irl with flat frames)
        Output: 
        - data array with added noise
        """
        nrows, ncols = data.shape
    
        bias_val = 1100
        bias_im = np.zeros((nrows, ncols), dtype=np.float32) + bias_val 
        if incl_bias:
            number_of_columns = 50
            columns = np.random.randint(0, ncols, size=number_of_columns)
            col_pattern = np.random.randint(0, int(0.8 * bias_val), size=nrows)
            for c in columns:
                bias_im[:, c] = bias_val + col_pattern
        data += bias_im
    
        noise = np.random.randint(100, 800, size=(nrows, ncols), dtype=np.uint16)
        data += noise
    
        if incl_sky:
            vin = cls.vignette(data)
            data += vin # adding vignette affect around edges
    
        # adding random hot pixels:
        for _ in range(100):
            col = np.random.randint(ncols)
            row = np.random.randint(nrows)
            data[row, col] += 65535

    @staticmethod
    def add_psf_object(array, cx, cy, flux, fwhm):
        """Funciton that creates a point source object (Gaussaian PSF) at a given location
        Inputs: 
        - array of data that the point source will be added to
        - cx, cy poition of center of point source
        - flux of object
        - fhwm of object
        Output:
        - creates point source object centered at the given point in the data array
        """
        
        sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        r = int(np.ceil(4 * sigma))
        nrows, ncols = array.shape
        x0, y0 = int(round(cx)), int(round(cy))

        for dy in range(-r, r + 1):
            for dx in range(-r, r + 1):
                xi, yi = int(x0 + dx), int(y0 + dy)
                if 0 <= xi < ncols and 0 <= yi < nrows:
                    array[yi, xi] += flux * np.exp(-(dx**2 + dy**2) / (2 * sigma**2))

    @staticmethod
    def line_points_difference(start_ra, start_dec, end_ra, end_dec, my_wcs):
        """Funciton to compute the x and y pixel difference between starting and ending RA, Dec of an object, 
        For use making streaks (either for the stars or the RSO)
        Inputs: 
        - starting RA and Dec (degrees)
        - ending RA and Dec (degrees)
        - wcs used
        Outpus: 
        - change in x (pixels)
        - change in y (pixels)
        """

        start_coord = SkyCoord(ra = start_ra * u.deg, dec =  start_dec * u.deg)
        end_coord = SkyCoord(ra = end_ra * u.deg, dec =  end_dec * u.deg)

        start_px, start_py = my_wcs.world_to_pixel(start_coord)
        start_px, start_py = float(start_px), float(start_py)
        start_ix, start_iy = int(round(start_px)), int(round(start_py))

        end_px, end_py = my_wcs.world_to_pixel(end_coord)
        end_px, end_py = float(end_px), float(end_py)
        end_ix, end_iy = int(round(end_px)), int(round(end_py))

        return (end_ix - start_ix), (end_iy - start_iy)

    @classmethod
    def add_blurred_line(cls, array: np.array, x1, y1, x2, y2, flux, fwhm):
        """Function to create a Gaussian-blurred line given begining and end points and object parameters
        Inputs: 
        - data array to create the line in
        - x1, y1: first object location
        - x2, y2: last ob ject location
        - flux
        - fwhm
        Output: 
        - creates a line streak in the data array
        """

        nrows, ncols = array.shape
        line = int(np.hypot(x2 - x1, y2 - y1))

        #checking if the line would be less than a pixel (so rounded to zero):
        if line == 0: 
            cls.add_psf_object(array, x1, y1, flux, fwhm)
            return
        else:
            flux_per_line_point = flux / (line  * fwhm)   # adding fwhm multiplier since plotted psf objects will overlap

        # drawing lines using an overlay:
        overlay = Image.new('RGBA', (ncols, nrows), (0, 0, 0, 0)) 
        draw = ImageDraw.Draw(overlay)
        COLOR_LEAPFROG = (255, 255, 255, 200) 

        draw.line(
                [(int(x1), int(y1)), (int(x2), int(y2))],
                fill = COLOR_LEAPFROG,
                width = int(fwhm)
        )
        
        nda = np.array(overlay)
        alpha = np.array(nda[:, :, 3])
        mask = alpha > 0  # making mask from area with streak
        ys, xs = np.where(mask)

        for y, x in zip(ys, xs):
            cls.add_psf_object(array, x, y, flux_per_line_point, fwhm)

    @staticmethod
    def make_labelme_shape(obj_type: str, x1, y1, x2 = None, y2 = None, fwhm = None):
        """Function to make a labelme shape label from object type and location"""

        if 'point' in obj_type:
            shape = {
                "label": obj_type,
                "points": [[x1, y1], [x1 + (3*fwhm), y1 + (3*fwhm)]],
                "group_id": None,
                "description": "",
                "shape_type": "circle",
                "flags": {},
                "mask": None
            }
        else:
            shape = {
                "label": obj_type,
                "points": [[x1, y1], [x2, y2]],
                "group_id": None,
                "description": "",
                "shape_type": "line",
                "flags": {},
                "mask": None
            }
        
        return shape
    


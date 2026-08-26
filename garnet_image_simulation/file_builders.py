"""
file_builders.py

Version date: 19/6/2026
Author: Bruce Ritter

Contains FileBuilders class for use with image_generator.py, a script
written to generate simulated telescope images of RSOs. 

List of Functions: 
- make_labelme_file
- make_rso_streak_image
- make_rso_point_image
- fits_to_png
- generate
"""

import numpy as np
import PIL
from PIL import Image
import json
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.wcs import WCS
from astropy.io import fits
import os
from pathlib import Path

from image_tools import ImageTools
from astro_info import AstroInfo


image_tools = ImageTools()
astro_info = AstroInfo()

class FileBuilders:
    
    @staticmethod
    def make_labelme_file(image_path, file_path, shapes):
        """Function to make a labelme file  
        Inputs: 
        - path to image
        - path to save json file to
        - list of shape dictionaries  (shapes information of image)
        Output: 
        - JSON file with label information"""

        with Image.open(image_path) as im: 
            w, h = im.size
        
        labelme_data = {
            "version": "6.0.0",
            "flags": {},
            "shapes": shapes,
            "imagePath": os.path.basename(image_path),
            "imageData": None,
            "imageHeight": h,
            "imageWidth": w
        }

        with open(f'{file_path}.json', "w", encoding="utf-8") as file:
            json.dump(labelme_data, file, ensure_ascii=False, indent=2)
        
        print(f"Label me file f'{file_path}.json' saved")

    @staticmethod
    def make_rso_streak_image(im_RA, im_Dec, f_im_RA, f_im_Dec, rso_flux, rso_fwhm, 
                            label, star_table, pix_sc, shape, sky_vin, bias, pos_change = 0):
        """Function to make a fits file populated with stars from gaia 
        around a given center coordinate with a RSO streak in the image, 
        simulating LEAPFROG mode where telescope is tracking sky

        Inputs: 
        - im_RA, im_Dec: initial RA and Dec used for star search, comes from RSO TLE
        - f_im_RA, f_im_Dec: final RA and Dec of RSO from TLE
        - flux of rso: rso_flux
        - fwhm of rso: rso_fwhm
        - label (string): image label
        - star_table: table of star query results
        - pix_sc (float): image pixel scale
        - shape: shape of array
        - pos_change: indicates how to change the RSO position to keep in the middle or offset
                either behind or ahead in what would be its motion accross the frame
        Output: 
        - creates (not returns) a .fits file image of the star field at the given coordinates with an included s
            simulated RSO streak
        - returns a list of shape dictionaries for a lableme label
        """
        shapes_list = [] # for labelme labels
        
        #print("Starting to make LEAPFROG image...")
        nrows, ncols = shape
        data = np.zeros(shape, dtype=np.float32)

        # setting WCS 
        myWCS = WCS(naxis=2)
        myWCS.wcs.crpix = [ncols // 2, nrows // 2]
        myWCS.wcs.crval = [im_RA, im_Dec]
        myWCS.wcs.cunit = ['deg', 'deg']
        myWCS.wcs.cdelt = [-pix_sc, -pix_sc]
        myWCS.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        myWCS.wcs.equinox = 2000.0
        myWCS.array_shape = [nrows, ncols]

        image_orientation = 0
        myWCS.wcs.pc = [[np.cos(np.deg2rad(image_orientation)), -np.sin(np.deg2rad(image_orientation))],
                        [np.sin(np.deg2rad(image_orientation)),  np.cos(np.deg2rad(image_orientation))]]

    
        # adding in rso:
        coord = SkyCoord(ra = im_RA * u.deg, dec =  im_Dec * u.deg)
        rso_x, rso_y = myWCS.world_to_pixel(coord)
        rso_x, rso_y = float(rso_x), float(rso_y)

        dx, dy = image_tools.line_points_difference(im_RA, im_Dec, f_im_RA, f_im_Dec, my_wcs = myWCS)
        pos_change = float(pos_change)

        if  0 < (rso_x + dx + (pos_change * dx)) < ncols and 0 < (rso_y + dy) + (pos_change * dy) < nrows: 
            image_tools.add_blurred_line(data, (rso_x + (pos_change * dx)), (rso_y + (pos_change * dy)), 
                            (rso_x + dx + (pos_change * dx)), (rso_y + dy) + (pos_change * dy), 
                            #(rso_flux / rso_fwhm), (rso_fwhm / 2)) 
                            rso_flux, (rso_fwhm / 2))

            shapes_list.append(image_tools.make_labelme_shape("rso_streak", (rso_x + (pos_change * dx)), (rso_y + (pos_change * dy)), 
                            (rso_x + dx + (pos_change * dx)), (rso_y + dy) + (pos_change * dy)))
            
            # populating stars within the frame of the image: 
            for star in star_table:

                coord = SkyCoord(ra = star['ra'] * u.deg, dec =  star['dec'] * u.deg)
                star_flux = star['flux']
                star_fwhm = astro_info.fwhm_from_mag(star['mag'])

                px, py = myWCS.world_to_pixel(coord)
                px, py = float(px), float(py)
                ix, iy = int(px), int(py)

                if 0 <= ix < shape[1] and 0 <= iy < shape[0]:
                    image_tools.add_psf_object(data, float(px), float(py), star_flux, star_fwhm)
                    shapes_list.append(image_tools.make_labelme_shape("star_point", px, py, px + (star_fwhm * 2), py + (star_fwhm * 2), star_fwhm))


            # adding general noise to image: 
            image_tools.add_noise(data, sky_vin, bias)

            # creating and saving .fits file: 
            header = myWCS.to_header()
            primary_hdu = fits.PrimaryHDU(data=data, header=header)
            hdu_list = fits.HDUList([primary_hdu])

            hdu_list.writeto(f'{label}.fits', overwrite=True)
            print(f"generated {label}.fits")

            return shapes_list
        
        else:
            return []

    @staticmethod
    def make_rso_point_image(im_RA, im_Dec, f_im_RA, f_im_Dec, rso_flux, rso_fwhm, 
                            label, star_table, pix_sc, shape, sky_vin, bias):
        """Function to make a fits file populated with streaked stars (from GAIA database) with a point RSO, 
        which simulates TRACKING (telescope is tracking RSO not sky)
        Inputs: 
        - im_RA, im_Dec: initial RA and Dec used for star search, comes from RSO TLE
        - f_im_RA, f_im_Dec: final RA and Dec of RSO from TLE
        - flux of rso: rso_flux
        - fwhm of rso: rso_fwhm
        - label (string): image label
        - star_table: table of star query results
        - pix_sc (float): image pixel scale
        - shape: shape of array
        Output: 
        - creates (not returns) a .fits file image of the star field at the given coordinates with an included s
            simulated RSO point, with RSO in the middle
        - returns a list of shape dictionaries for a lableme label
        """
        shapes_list = [] # for labelme labels
        
        nrows, ncols = shape
        data = np.zeros(shape, dtype=np.float32)

        # setting WCS 
        myWCS = WCS(naxis=2)
        myWCS.wcs.crpix = [ncols // 2, nrows // 2]
        myWCS.wcs.crval = [im_RA, im_Dec]
        myWCS.wcs.cunit = ['deg', 'deg']
        myWCS.wcs.cdelt = [-pix_sc, -pix_sc]
        myWCS.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        myWCS.wcs.equinox = 2000.0
        myWCS.array_shape = [nrows, ncols]

        image_orientation = 0
        myWCS.wcs.pc = [[np.cos(np.deg2rad(image_orientation)), -np.sin(np.deg2rad(image_orientation))],
                        [np.sin(np.deg2rad(image_orientation)),  np.cos(np.deg2rad(image_orientation))]]

        # populating stars within the frame of the image: 
        dx, dy = image_tools.line_points_difference(im_RA, im_Dec, f_im_RA, f_im_Dec, my_wcs = myWCS)

        for star in star_table:

            coord      = SkyCoord(ra=star['ra'] * u.deg, dec=star['dec'] * u.deg)
            star_flux  = star['flux']
            star_fwhm  = astro_info.fwhm_from_mag(star['mag'])
            px, py     = myWCS.world_to_pixel(coord)   # (x=col, y=row)
            px, py     = float(px), float(py)
            ix, iy     = int(round(px)), int(round(py))
            if (0 <= (ix + dx) < shape[1] and 0 <= (iy + dy) < shape[0]):
                image_tools.add_blurred_line(data, px, py, px + dx, py + dy, 
                                             #(star_flux / star_fwhm), star_fwhm)
                                             star_flux, star_fwhm / 2)

                # only adding label to line if it is fully wihtin frame: 
                if  0 <= px  < ncols and 0 <= py < nrows:
                    if  0 <= (px + dx) < ncols and 0 <= (py + dy) < nrows:
                        shapes_list.append(image_tools.make_labelme_shape("star_streak", px, py, px + dx, py + dy))

        #adding rso:
        rso_coord = SkyCoord(ra=im_RA * u.deg, dec=im_Dec * u.deg)
        rso_cx, rso_cy = myWCS.world_to_pixel(rso_coord)
        rso_cx, rso_cy = float(rso_cx), float(rso_cy)
        image_tools.add_psf_object(data, rso_cx, rso_cy, rso_flux, rso_fwhm)
        shapes_list.append(image_tools.make_labelme_shape("rso_point", rso_cx, rso_cy, 
                                            (rso_cx + rso_fwhm * 2),(rso_cy + rso_fwhm * 2), rso_fwhm))
        # adding noise to image: 
        image_tools.add_noise(data, sky_vin, bias)

        # creating and saving .fits file: 
        header = myWCS.to_header()
        primary_hdu = fits.PrimaryHDU(data=data, header=header)
        hdu_list = fits.HDUList([primary_hdu])

        hdu_list.writeto(f'{label}.fits', overwrite=True)
        print(f"generated {label}.fits")

        return shapes_list

    @staticmethod
    def fits_to_png(label, file_path):
        with fits.open(file_path) as hdul:
            data = hdul[0].data.astype(np.float64)

        # Background subtract
        bkg = np.median(data)
        data = data - bkg
        mad = np.median(np.abs(data - np.median(data)))
        sigma = mad / 0.6745

        data = np.clip(data, -3 * sigma, np.percentile(data, 99.99))

        # stretch
        a = 3 * sigma
        data_stretched = np.arcsinh(data / a)

        # Normalize
        lo = np.arcsinh(-3 * sigma / a)
        hi = np.percentile(data_stretched, 99.0) 
        data_norm = np.clip((data_stretched - lo) / (hi - lo) * 255, 0, 255)

        img = Image.fromarray(data_norm.astype(np.uint8))
        img.save(f"{label}.png")
        print(f"saved {label}.png")

    @classmethod
    def generate(cls, coords, rso_mag, star_table, pix_sc, array_shape,
                mode, file_path, exp_time, sky_vin, bias, n_image,
                obs_time=None, sat_id=None, create_time = None):
        """Function to generate an image or images in either
        LEAPFROG mode or TRACKING mode
        Inputs:
        - coords: list of RA and Dec coords for start and end RSO position
        - rso mag: magnitude of RSO
        - array_shape: shape of array to create
        - mode: mode of image created (LEAPFROG or TRACKING). Will also determine number of images created
        - file_path: path to save files too
        - exp_time: float exposure time for image
        - observer: string name of observer
        Output:
        - writes 1 image set (.fits, .png, .json) of the mode specified
        """
        
        # rso flux and fwhm for image:
        rso_flux = astro_info.g_mag_to_flux(rso_mag, exp_time)
        rso_fwhm = astro_info.fwhm_from_mag(rso_mag)

        # label/path for files with mode, obs time, exp time, sat id:
        if obs_time is not None:
            timestamp = obs_time.strftime('%Y_%m_%d_%H_%M_%S')
        else:
            timestamp = f"img_{n_image}"

        sat_tag = sat_id if sat_id is not None else 'sat'

        # folder name format: MODE_CREATETIME_NORAD for each run, observation timestamp used if no creation time specified
        folder_time = create_time if create_time is not None else timestamp
        output_file = Path(file_path) / f'{mode}_{folder_time}_{sat_tag}'
        output_file.mkdir(exist_ok=True, parents=True)
        base_name = f'{timestamp}_exp_{exp_time:.2f}_{mode}_{sat_tag}'
        label = str(output_file / base_name)

        # generating image or images based on mode:
        if mode == 'TRACKING':
            shapes_list = cls.make_rso_point_image(coords[0], coords[1], coords[2], coords[3], rso_flux, rso_fwhm,
                                label, star_table, pix_sc, array_shape, sky_vin, bias)

            cls.fits_to_png(label, f'{label}.fits')
            cls.make_labelme_file(f'{label}.png', label, shapes_list)

        elif mode == 'LEAPFROG':
            # randomly determine if RSO will be centered on image or not, with possible movement ahead or behind calculated RSO location:
            pos_change = np.random.randint(-1, 2, 1) * 3
            shapes_list = cls.make_rso_streak_image(coords[0], coords[1], coords[2], coords[3], rso_flux, rso_fwhm,
                                label, star_table, pix_sc, array_shape, sky_vin, bias, pos_change)
            
            if shapes_list == []:  #if image cannot be made with the position change included, use much smaller pos_change
                pos_change = pos_change * 0.01
                shapes_list = cls.make_rso_streak_image(coords[0], coords[1], coords[2], coords[3], rso_flux, rso_fwhm,
                                label, star_table, pix_sc, array_shape, sky_vin, bias, pos_change)

            if shapes_list == []:  #if image still cannot be made, use pos_change = 0 
                pos_change = 0
                shapes_list = cls.make_rso_streak_image(coords[0], coords[1], coords[2], coords[3], rso_flux, rso_fwhm,
                                label, star_table, pix_sc, array_shape, sky_vin, bias, pos_change)  
            
            cls.fits_to_png(label, f'{label}.fits')
            cls.make_labelme_file(f'{label}.png', label, shapes_list)
                
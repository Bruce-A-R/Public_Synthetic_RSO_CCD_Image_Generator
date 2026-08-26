"""
image_generator.py

Version Date: 19/6/2026
Author: Bruce Ritter

Description: 
Use to generate simulated images of RSOs with background stars in .fits and .png form 
with an included json file of labelme labels for objects within the image frame

Function will generate a number of images each with .fits, .png, and labelme .json 
which will show the "answer" of where in the image there is a: 
    - rso_point
    - rso_streak
    - star_point
    - star_streak

This script is designed for use training GISTDA machine learning algorithms to automate 
RSO detection in images.
"""

#--------imports---------------
import argparse
import os
import numpy as np
import PIL
from PIL import Image, ImageDraw
import json
from datetime import timedelta, datetime
from math import degrees
import warnings
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.wcs import WCS
from astropy.io import fits
from astroquery.gaia import Gaia
import gaiaoffline
import ephem
from glob import glob


from astro_info import AstroInfo as astro_info
from coord_tools import CoordTools as coord_tools
from file_builders import FileBuilders as file_builders


# -------------------

def image_making(file_path, n_generated, mode_generated, exp_time, rso_mag, mag_lim, offline, sky_vin, bias, 
                 catalog, cat_file, line1, line2, start, duration, obs_info, nrow, ncol, pix_sc, sat_id, create_time):
    """Function to do most of the main function tasks, so that it can be repeated if after running through all iterations in one pass length, 
    the number of images generated per mode per TLE still does not meet n_generated"""

    iterations = duration // 5 # duration divided by max exp time amount

    # get pictures for every unit of exposure time across pass duration 
    i = 1
    n_gen = 0
    
    while i < iterations:  # or n_gen <= n_generated:

        if exp_time == None:
            rng = np.random.default_rng()
            exp_time = float(np.round(rng.uniform(0.3, 5.0), decimals=2))

        if rso_mag is None:
            rng = np.random.default_rng()
            frame_rso_mag = float(np.round(rng.uniform(5.0, 12.0), decimals=2))
        else:
            frame_rso_mag = rso_mag

        obs_time = start + timedelta(seconds = exp_time * i)
        coords = coord_tools.coords_from_TLE(line1, line2, exp_time, obs_time, *obs_info)
        shape = (nrow, ncol)

        # check if RSO streak is longer than half the image size, and reduce exp time if so:
        frame_wcs = WCS(naxis=2)
        frame_wcs.wcs.crpix = [ncol // 2, nrow // 2]
        frame_wcs.wcs.crval = [coords[0], coords[1]]
        frame_wcs.wcs.cdelt = [-pix_sc, -pix_sc]
        frame_wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        end_x, end_y = frame_wcs.world_to_pixel(SkyCoord(ra=coords[2] * u.deg, dec=coords[3] * u.deg))

        if not (0 <= end_x < ncol and 0 <= end_y < nrow):
            exp_time = exp_time * 0.1
            coords = coord_tools.coords_from_TLE(line1, line2, exp_time, obs_time, *obs_info)

        # query star catalog for background stars:
        print("Querying for stars...")
        radius_pix = np.hypot((ncol / 2), (nrow / 2))
        radius_deg = radius_pix * pix_sc
        #star_table = astro_info.star_query(coords[0], coords[1], radius_deg, mag_lim, offline, hip_file, exp_time)
        star_table = astro_info.star_query(coords[0], coords[1], radius_deg, mag_lim, offline, catalog, cat_file, exp_time)

        # making the files:
        print("Making files...")
        if mode_generated == 'TRACKING':
            file_builders.generate(coords, frame_rso_mag, star_table, pix_sc, shape, 'TRACKING', file_path, exp_time, 
                sky_vin, bias, i, obs_time = obs_time, sat_id = sat_id, create_time = create_time)
        elif mode_generated == 'LEAPFROG':
            file_builders.generate(coords, frame_rso_mag, star_table, pix_sc, shape, 'LEAPFROG', file_path, exp_time, 
                sky_vin, bias, i, obs_time = obs_time, sat_id = sat_id, create_time = create_time)
        else:
            file_builders.generate(coords, frame_rso_mag, star_table, pix_sc, shape, 'TRACKING', file_path, exp_time, 
                sky_vin, bias, i, obs_time = obs_time, sat_id = sat_id, create_time = create_time)

            file_builders.generate(coords, frame_rso_mag, star_table, pix_sc, shape, 'LEAPFROG', file_path, exp_time, 
                sky_vin, bias, i, obs_time = obs_time, sat_id = sat_id, create_time = create_time)

        n_gen += 1
        if n_gen == n_generated:   #stop loop if number of generated images reaches requested amount
            return n_gen

        i +=1

    return n_gen


def main(file_path, TLE, n_generated, mode_generated, exp_time, when, observer, obs_lat, obs_lon, obs_alt,
         rso_mag, mag_lim, pix_size, nrows, ncols, focal_len, binning, offline, sky_vin, bias, catalog, cat_file):  
    """Main function: generates .fits, .png images with a JSON file of label information 
    to a specified path, based on inputted TLE and other optional parameters
    """

    if TLE == None:  # TLE information
        # example ISS TLE (from https://live.ariss.org/tle/) 
        line1_list = ["1 25544U 98067A   26159.80410962  .00007129  00000-0  13425-3 0  9990"]
        line2_list = ["2 25544  51.6336 341.5878 0006923 148.5365 211.6039 15.49672912570453"]
    else:
        line1_list, line2_list = coord_tools.tle_parser(TLE)

    # single wall-clock timestamp for this run; groups all output folders of one run together
    create_time = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')

    for line1, line2 in zip(line1_list, line2_list):

        # NORAD catalog number from TLE line 1 (columns 3-7); keeps file names unique per satellite
        sat_id = line1[2:7].strip()

        # getting observer information then getting coords from TLE and observer information:
        obs_info, pix_sc, nrow, ncol = coord_tools.observer_info(observer, obs_lat, obs_lon, obs_alt, pix_size, 
                                                                focal_len, nrows, ncols, binning) 

        # finding next pass based on TLE and observer location
        start, end, duration = coord_tools.pass_time_from_TLE(line1, line2, *obs_info, when)

        n_gen_total = 0

        while n_gen_total < n_generated:
            #create images like normal in iterations: 
            n = image_making(file_path, n_generated, mode_generated, exp_time, rso_mag, mag_lim, offline, sky_vin, bias, 
                    catalog, cat_file, line1, line2, start, duration, obs_info, nrow, ncol, pix_sc, sat_id, create_time)
            
            n_gen_total += n
            if n_gen_total <= n_generated:
                start = start + timedelta(seconds = 86400)
                start, end, duration = coord_tools.pass_time_from_TLE(line1, line2, *obs_info, when)
            else: 
                break
            
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''This script will create a simulated .fits and .png image along with labelme
                        label data for stars and an RSO given the RSO TLE and magnitude, observer information, telescope information''')
    parser.add_argument('file_path', type = str, help = 'Output path for created files')
    parser.add_argument('-t', '--TLE', type = str, help = "file path of TLE or TLEs", default = None)
    parser.add_argument('-n', '--n_generated', type = int, help = "# files to make", default = 1)
    parser.add_argument('-m', '--mode_generated', type = str, help = "type of file to make: TRACKING, LEAPGROG, or both", default = None)
    parser.add_argument('-e', '--exp_time', type = float, help = "exposure time", default = None)
    parser.add_argument('-w', '--when', type = str, help = "observation time, format = YYYY-MM-DD HH:MM:SS.SSSS", default = None)
    parser.add_argument('-ob', '--observer', type = str, help = "name of observer (SKP, DSC, or other)",default = None)
    parser.add_argument('-lat', '--obs_lat', type = float, help = "observer latitude",default = None) 
    parser.add_argument('-lon', '--obs_lon', type = float, help = "observer longitude",  default = None)
    parser.add_argument('-alt', '--obs_alt', type = float, help = "observer altitude", default = None)
    parser.add_argument('-rm', '--rso_mag', type = float, help = "RSO magnitude",  default = None)
    parser.add_argument('-ml', '--mag_lim', type = float, help = "stars' limiting magnitude", default = 16)
    parser.add_argument('-p', '--pix_size', type = float, help = "sensor pixel size",  default = None)
    parser.add_argument('-nr', '--nrows', type = int, help = "image rows #",  default = None)
    parser.add_argument('-nc', '--ncols', type = int, help = "image collumns #",  default = None)
    parser.add_argument('-f', '--focal_len', type = float, help = "telescope focal length",  default = None)
    parser.add_argument('-bin', '--binning', type = int, help = 'binning factor for image', default = 2)
    parser.add_argument('-ofl', '--offline', type = bool, help = "using offline gaia catalog True/False", default = False)
    parser.add_argument('-s', '--sky_vin', type = bool, help = "including vignette affect (from instrument reflections, sky brightness) True/False", default = False)
    parser.add_argument('-b', '--bias', type = bool, help = "including bias noise True/False", default = False)
    parser.add_argument('-c', '--catalog', type = str, help = "Catalog to use: G, H, N", default = "G")
    parser.add_argument('-cf', '--cat_file', type = str, help = "File path to catalog used (hipparchos or NOMAD)", default = None)

    args = parser.parse_args()

    main(args.file_path, args.TLE, args.n_generated, args.mode_generated, args.exp_time, args.when, args.observer, 
        args.obs_lat, args.obs_lon, args.obs_alt, args.rso_mag, args.mag_lim, args.pix_size, args.nrows, args.ncols, 
        args.focal_len, args.binning, args.offline, args.sky_vin, args.bias, args.catalog, args.cat_file)


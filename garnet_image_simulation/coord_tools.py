"""
coord_tools.py

Version date: 19/6/2026
Author: Bruce Ritter

Contains CoordTools class for use with image_generator.py, a script
written to generate simulated telescope images of RSOs.

List of Functions:
- get_pix_scale
- observer_info
- coords_from_TLE
- pass_time_from_TLE
- tle_parser
"""

import astropy.units as u
from astropy.coordinates import SkyCoord
import ephem
from math import degrees
from datetime import timedelta, datetime, timezone


class CoordTools:

    @staticmethod
    def get_pix_scale(pix_size , focal_len, binning = 2):
        """Funciton to get the pizel scale using pixel size and focal length
        Inputs: 
            pixeize (float): in um, focallen (float): in mm
        Outputs: 
            pixscale (float) in degrees
        """
        deg = (pix_size * 206.265 / focal_len) / 3600.0
        return deg * binning 

    @classmethod
    def observer_info(cls, observer, lat, lon, alt, pix_size, focal_len, nrow, ncol, binning):
        """Function to get observer information based on either direct inputs or string input of location acronym
        Inputs: 
        - observer (str): observation location, can handle input of SKP or DCS
        Output: 
        - list of observer lat, lon, alt
        - pixel scale (float)
        - number of rows in image (int)
        - number of columns in image (int)
        - pixel size of camera in um (float)
        - telescope focal lenght (float)
        """
        if observer == 'SKP': 
            if not lat: 
                lat = 13.10141
            if not lon: 
                lon = 100.9294 
            if not alt: 
                alt = 49.8

            if not pix_size:
                pix_size = 9.0
            if not focal_len:
                focal_len = 2563.0

            pix_sc = cls.get_pix_scale(pix_size, focal_len, binning)

            if not nrow: 
                nrow = 4096
            if not ncol: 
                ncol = 4096

        elif observer == 'DSC':
            if not lat: 
                lat = -30.526
            if not lon:
                lon = -70.85
            if not alt:
                alt = 1710.0

            if not pix_size:
                pix_size = 3.76
            if not focal_len: 
                focal_len = 1050.0

            pix_sc = cls.get_pix_scale(pix_size, focal_len, binning)
            if not nrow: 
                nrow = 2128
            if not ncol: 
                ncol = 3192

        else:  #handling other observer inputs or no observer input
            if not lat: 
                lat = 13.10141
            if not lon: 
                lon = 100.9294 
            if not alt: 
                alt = 49.8

            if not pix_size:
                pix_size = 3.76
                #print("pix size not specified, set to 3.76 um")
            if not focal_len:
                focal_len = 1050.0
                #print("focal length not specified, set to 1050.0 mm")
            
            pix_sc = cls.get_pix_scale(pix_size, focal_len, binning)

            if not nrow: 
                nrow = 2128
            if not ncol:
                ncol = 3192

        obs_info = [lat, lon, alt]
        return obs_info, pix_sc, nrow, ncol

    @staticmethod
    def coords_from_TLE(line1, line2, exp_time, obs_time, obs_lat, obs_lon, obs_alt):
        """Function to get RA and Dec coordniate information (starting and ending) from an inputted TLE
        Inputs: 
        - Line 1 of TLE (string)
        - Line 2 of TLE (string)
        - exp_time in seconds
        - obs_time as datetime
        - observing latitude, longitude (degrees)
        - observing altitude (meters)
        Outputs: 
        - coords list with: starting RA of sat, starting Dec of sat, ending RA and ending Dec of sat
        """
        observer = ephem.Observer()
        observer.lat = str(obs_lat)
        observer.lon = str(obs_lon)
        observer.elevation = float(obs_alt)
        observer.pressure = 0

        exposure_time = float(exp_time)

        # compute center time:
        sat = ephem.readtle('0 TRACK_OBJ', line1, line2)
        observer.date = obs_time
        sat.compute(observer)
        coord = SkyCoord(ra= degrees(sat.a_ra) * u.deg, dec= degrees(sat.a_dec) * u.deg)
        #print(f"Center Coord of Image: {coord}")

        # compute start time (obs_time - exp/2):
        observer.date = obs_time - timedelta(seconds = exposure_time / 2)
        sat = ephem.readtle('0 TRACK_OBJ', line1, line2)
        sat.compute(observer)
        coord_1 = SkyCoord(ra= degrees(sat.a_ra) * u.deg, dec = degrees(sat.a_dec) * u.deg)

        # compute end time (obs_time + exp/2):
        observer.date = obs_time + timedelta(seconds = exposure_time / 2)
        sat = ephem.readtle('0 TRACK_OBJ', line1, line2)
        sat.compute(observer)
        coord_2 = SkyCoord(ra = degrees(sat.a_ra) * u.deg, dec = degrees(sat.a_dec) * u.deg)

        return [coord_1.ra.degree, coord_1.dec.degree, coord_2.ra.degree, coord_2.dec.degree]
    
    @staticmethod
    def pass_time_from_TLE(line1, line2, obs_lat, obs_lon, obs_alt, when):
        """Function to find a functional observering time from a TLE, where the sat would be above the observer's horizon
        Inputs:
        - TLE lines 1 and 2
        - observer lat
        - observer lon
        - observer alt
        returns: 
        - start time
        - end time
        - total seconds in pass
        """

        
        utc_now = datetime.now(timezone.utc)

        observer = ephem.Observer()
        observer.lat = str(obs_lat)
        observer.lon = str(obs_lon)
        observer.elevation = float(obs_alt)
        observer.pressure = 0

        if when == None: 
            utc_now = datetime.now(timezone.utc)
            observer.date = utc_now
        else: 
            observer.date = when

        sat = ephem.readtle('0 TRACK_OBJ', line1, line2)
        
        try:   # trying to find next pass, will fail for circumpolar sats
            pass_time = observer.next_pass(sat)
            t0 = pass_time[0]
            tf = pass_time[4]

            pass_duration = tf.datetime() - t0.datetime()
            pass_duration_sec = pass_duration.total_seconds()
            return t0.datetime(), tf.datetime(), pass_duration_sec
        except: 
            t0 = observer.date
            t0 = t0.datetime()
            tf = observer.date
            tf = tf.datetime() + timedelta(seconds = 600)
            pass_duration_sec = 600
            return t0, tf, pass_duration_sec


    @staticmethod
    def tle_parser(TLE_file):
        """Function to take .txt file of TLE and return two lines (the first and second ones)
        Input: a text file containing one TLE
        Output: the first and second lines of that TLE (str)"""
        
        #line1, line2 = None, None

        line1_list = []
        line2_list = []

        with open(TLE_file, "r") as file: 
            lines = file.readlines()

            for line in lines: 
                if line[0:2] == '1 ':
                    line1_list.append(line) # = line
                elif line[0:2] == '2 ':
                    line2_list.append(line) # = line
            
        return line1_list, line2_list
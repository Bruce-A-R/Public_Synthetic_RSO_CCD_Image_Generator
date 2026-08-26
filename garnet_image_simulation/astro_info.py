"""
astro_info.py

Version date: 19/6/2026
Author: Bruce Ritter

Contains AstroImage class for use with image_generator.py, a script
written to generate simulated telescope images of RSOs.

List of Functions: 
- load_hipparchos
- flux_from_v_mag
- angular_seperation
- hipparchos_query
- star_query
- g_mag_to_flux
- fwhm_from_mag
"""

import numpy as np
import warnings
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.io import fits
from astroquery.gaia import Gaia
import gaiaoffline


class AstroInfo:
    
    @staticmethod
    def load_hipparcos(fits_path):
        """
        Function to load in the hipparchos catalog from the specified fits file containing it, 
        downloaded from https://cdsarc.cds/unistra.fr/ftp/cats/I/239/

        Inputs:
        - the file to load in the catalog from
        Outputs:
        - dictionary with ra, dec,  v magnitude, and v flux for stars 
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")          # suppress column-name warnings
            with fits.open(fits_path, memmap=True) as hdul:
                raw = hdul[1].data
    
                def to_float(col_name):
                    col = raw[col_name]
                    out = np.empty(len(col), dtype=np.float64)
                    for i, v in enumerate(col):
                        s = str(v).strip()
                        try:
                            out[i] = float(s)
                        except (ValueError, TypeError):
                            out[i] = np.nan
                    return out
    
                ra   = to_float('RAdeg')
                dec  = to_float('DEdeg')
                plx  = to_float('Plx')
                vmag = to_float('Vmag')

        distances = []
        for val in plx: 
            if val > 0:
                distances.append(1000 / val)
            else: 
                distances.append(0)


        dictionary = {
            'ra': ra,
            'dec': dec,
            'plx': plx,
            'vmag': vmag,
            'dist': distances
        }
        return dictionary
    
    @staticmethod
    def flux_from_v_Vega(v_mag):
        """Function to estimate flux from object v magnitude using Vega as assumed mangitude zero-point, 
        Designed for use with offline locall downloaded NOMAD offline data base. May not match with GAIA filters
        """
        #m1 - m2 = -2/15log10(F1/F2)
        #vega flux in V: fn = 3.50 ´ 10 –20 erg s–1 cm–2 Hz–1, need to check units 
        f_vega = 3.50 * 10**-20
        flux = (10 * (v_mag / -2.5) ) * f_vega
        return flux 


    @staticmethod
    def flux_from_v_mag(v_mag, distance):
        """Function to get the flux of a star from exposure time, distance, and magnitude from specifically
        the Hipparchos catalog. ZP will not be correct for other catalogs
        
        Note: Currenlty being worked on, may return inaccurate fluxes"""

        app_mag = v_mag + (2.5 * np.log10((distance / 10)**2)) 
        flux = 10**((app_mag - 250) / -2.5) #V band fzp = vega v band for AT-HYG cited as = 3636 * e^-20 but that's not working

        return flux
    
    @staticmethod
    def angular_seperation(ra1_deg, dec1_deg, ra2_arr, dec2_arr):
        """Function to find the seperation between points and center point
        Inputs: 
        - staring RA and Dec (in degrees)
        - ending RA and dec (in degrees)
        Output:
        - angular seperation between start and end point (deg)
        """
        ra1  = np.deg2rad(ra1_deg)
        dec1 = np.deg2rad(dec1_deg)
        ra2  = np.deg2rad(ra2_arr)
        dec2 = np.deg2rad(dec2_arr)
    
        d_ra  = (ra2  - ra1)  / 2.0
        d_dec = (dec2 - dec1) / 2.0
    
        a = np.sin(d_dec)**2 + np.cos(dec1) * np.cos(dec2) * np.sin(d_ra)**2
        return np.degrees(2.0 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0))))

    @classmethod
    def hipparchos_query(cls, hip_file, c_ra, c_dec, radius, mag_lim):
        """Function to search through a downloaded hipparchos catalog of stars 
        within a certain radius from an ra, dec point and witha magnitude limit on stars returned
        Inputs: 
        - hipparchos catalog file
        - image center RA
        - image center Dec
        - search radius (deg)
        - limit magnitude for selecting stars
        Output: 
        - table of found stars in search area from hipparchos/Tycho catalog
        """
        star_dict = cls.load_hipparcos(hip_file)
        print(f"stars found: {len(star_dict)}")
        
        mag_filter = np.isfinite(star_dict["vmag"]) & (star_dict["vmag"] <= mag_lim)

        #defining area wihtin radius of points, masking stars w/in that area
        dec_lo  = c_dec - radius
        dec_hi  = c_dec + radius
        cos_dec = max(np.cos(np.deg2rad(c_dec)), 1e-6)
        ra_margin = min(radius / cos_dec, 180.0)
        ra_lo   = c_ra - ra_margin
        ra_hi   = c_ra + ra_margin
    
        ra_mod  = star_dict["ra"] % 360.0                         
        if ra_lo < 0:
            ra_wrap = (ra_mod >= (ra_lo % 360)) | (ra_mod <= ra_hi % 360)
        elif ra_hi >= 360:
            ra_wrap = (ra_mod >= ra_lo) | (ra_mod <= (ra_hi % 360))
        else:
            ra_wrap = (ra_mod >= ra_lo) & (ra_mod <= ra_hi)
    
        bbox_mask = mag_filter & (star_dict["dec"] >= dec_lo) & (star_dict["dec"] <= dec_hi) & ra_wrap

        idx_box = np.where(bbox_mask)[0]

        idx_box = np.where(bbox_mask)[0]
        if idx_box.size == 0:
            print("table could be empty")
            return []
    
        sep = cls.angular_seperation(
            c_ra, c_dec,
            star_dict["ra"][idx_box], star_dict["dec"][idx_box]
        )
    
        cone_mask = sep <= radius
        idx_cone = idx_box[cone_mask]
        sep_cone = sep[cone_mask]
    
        if idx_cone.size == 0:
            print("table could be empty")
            return []

        order = np.argsort(sep_cone)
        idx_final = idx_cone[order]
        sep_final = sep_cone[order]
        plx_final = star_dict["plx"][idx_final]
        dist_final = np.where(plx_final > 0, 1000.0 / plx_final, np.nan)

        v_flux = []
        for mag, dist in zip(star_dict["vmag"][idx_final], dist_final):
            v_flux.append(cls.flux_from_v_mag(mag, dist))
    
        #output table:
        result = Table({
            'ra': star_dict["ra"][idx_final],
            'dec': star_dict["dec"][idx_final],
            'vmag': star_dict["vmag"][idx_final],
            'parallax_mas': plx_final,
            'distance_pc': dist_final,
            'sep_deg': sep_final,
            'v_flux': v_flux
        })

        # table making: 
        result['ra'].unit = 'deg'
        result['dec'].unit = 'deg'
        result['vmag'].unit = 'mag'
        result['parallax_mas'].unit = 'mas'
        result['distance_pc'].unit = 'pc'
        result['sep_deg'].unit = 'deg'
    
        result['ra'].format = '.6f'
        result['dec'].format = '.6f'
        result['vmag'].format = '.2f'
        result['parallax_mas'].format = '.2f'
        result['distance_pc'].format  = '.1f'
        result['sep_deg'].format = '.4f'
    
        result.meta['center_ra'] = c_ra
        result.meta['center_dec'] = c_dec
        result.meta['radius_deg'] = radius
        result.meta['mag_limit'] = mag_lim
        result.meta['n_stars'] = len(result)
        result.meta['catalogue'] = 'Hipparcos I/239 (ESA 1997)'
        result.meta['epoch'] = 'ICRS J1991.25'
    
        return result
    

    @staticmethod
    def load_nomad(nomad_file):
        """Function to load in the NOMAD catalog from the specified fits file containing it, 
        downloaded from https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/I/297?format=html&tex=true

        Inputs:
        - the file to load in the catalog from
        Outputs:
        - dictionary with ra, dec, v magnitude, and v flux for all catalog stars 
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")          # suppress column-name warnings
            with fits.open(fits_path, memmap=True) as hdul:
                raw = hdul[1].data
    
                def to_float(col_name):
                    col = raw[col_name]
                    out = np.empty(len(col), dtype=np.float64)
                    for i, v in enumerate(col):
                        s = str(v).strip()
                        try:
                            out[i] = float(s)
                        except (ValueError, TypeError):
                            out[i] = np.nan
                    return out
    
                ra   = to_float('RAdeg')
                dec  = to_float('DEdeg')
                vmag = to_float('Vmag')

        #using vmag to get v flux:
        vflux = []
        for val in vmag: 
            vflux.append(flux_from_v_Vega(val))
        

        dictionary = {
            'ra': ra,
            'dec': dec,
            'vmag': vmag,
            'flux': vflux
        }
        return dictionary

    
    @classmethod
    def nomad_query(cls, RA, Dec, radius, mag_lim, nomad_file):
        """Function to query the NOMAD combined data base, hopefully providing more stars per frame than just Hipparchos
        Inputs: 
        - center RA (deg)
        - center Dec (deg) 
        - search radius (deg)
        - limiting V magnitude
        Outputs: 
        - table of star location(J2000 RA, Dec), v magnitude, v flux
        """

        # load in collumn data from NOMAD file: 
        star_dict = cls.load_nomad(nomad_file)

        #cone search is same as hipparchos query:
        mag_filter = np.isfinite(star_dict["vmag"]) & (star_dict["vmag"] <= mag_lim)

        #defining area wihtin radius of points, masking stars w/in that area
        dec_lo  = Dec - radius
        dec_hi  = Dec + radius
        cos_dec = max(np.cos(np.deg2rad(c_dec)), 1e-6)
        ra_margin = min(radius / cos_dec, 180.0)
        ra_lo   = RA - ra_margin
        ra_hi   = RA + ra_margin
    
        ra_mod  = star_dict["ra"] % 360.0                         
        if ra_lo < 0:
            ra_wrap = (ra_mod >= (ra_lo % 360)) | (ra_mod <= ra_hi % 360)
        elif ra_hi >= 360:
            ra_wrap = (ra_mod >= ra_lo) | (ra_mod <= (ra_hi % 360))
        else:
            ra_wrap = (ra_mod >= ra_lo) & (ra_mod <= ra_hi)
    
        bbox_mask = mag_filter & (star_dict["dec"] >= dec_lo) & (star_dict["dec"] <= dec_hi) & ra_wrap

        idx_box = np.where(bbox_mask)[0]

        idx_box = np.where(bbox_mask)[0]
        if idx_box.size == 0:
            print("table could be empty")
            return []
    
        sep = cls.angular_seperation(
            RA, Dec,
            star_dict["ra"][idx_box], star_dict["dec"][idx_box]
        )
    
        cone_mask = sep <= radius
        idx_cone = idx_box[cone_mask]
        sep_cone = sep[cone_mask]
    
        if idx_cone.size == 0:
            print("table could be empty")
            return []

        order = np.argsort(sep_cone)
        idx_final = idx_cone[order]


        result = Table({
            'ra': star_dict["ra"][idx_final],
            'dec': star_dict["dec"][idx_final],
            'vmag': star_dict["vmag"][idx_final],
            'flux': star_dict["vflux"][idx_final]
        })

        result['ra'].unit = 'deg'
        result['dec'].unit = 'deg'
        result['vmag'].unit = 'mag'

    
        result['ra'].format = '.6f'
        result['dec'].format = '.6f'
        result['vmag'].format = '.2f'
    
        result.meta['center_ra'] = RA
        result.meta['center_dec'] = Dec
        result.meta['radius_deg'] = radius
        result.meta['mag_limit'] = mag_lim
        result.meta['n_stars'] = len(result)
        result.meta['catalogue'] = 'NOMAD I/297 (Zacharias+ 2005)'
        result.meta['epoch'] = 'ICRS J2000'

        return result


    @staticmethod
    def g_mag_from_flux(flux, exp_time):
        """Function to convert GAIA result mean flux to mean mag, 
        using GAIA DR3 reported G filter zero point
        """
        m = -2.5 * np.log10(flux) + 25.6884
        return m


    @classmethod
    def star_query(cls, RA, Dec, radius, mag_lim, offline, catalog, cat_file, exp_time):
        """Function to query databases for stars alist of stars within a given radius from a given RA and Dec.
        Currently it can use astroquery to search stars in the online GAIA catalog, or using an offline catalog depending on 
         what is downloaded on the users device: the gaiaoffline catalog, a Hipparchos/Tycho catalog, 
         or the NOMAD combined catalog (last one untested)

        inputs: 
        - center RA (comes from RSO position)
        - center Dec (also from RSO position)
        - mag_lim: limiting magnitude (dimmest GAIA has is ~20 but dimmest a ground telescope will see is likely brighter (12-16))
        - offline: boolean, should gaiaoffline be used (you have to download the catalog) or astroquery online
        output: 
        - a table of stars within the image frame with magnitude and flux information included
        """
        #print(f'ra: {RA}, dec: {Dec}, radius: {radius}, mag_lim: {mag_lim}, offline: {offline}')
        
        
        if offline == False: # online query, which hopefully doesn't take too long with max 150 sources in table
            ra = RA * u.deg
            dec = Dec * u.deg

            coord = SkyCoord(ra, dec, frame = 'icrs')
            # print(f'{coord.ra.deg}, {coord.dec.deg}')
            query = f"""
            SELECT TOP 100 source_id, ra, dec, phot_g_mean_mag, parallax, pmra, pmdec, phot_g_mean_flux 
            FROM gaiadr3.gaia_source
            WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {coord.ra.deg}, {coord.dec.deg}, {radius})) = 1
            AND phot_g_mean_mag < {mag_lim}
            ORDER BY phot_g_mean_mag ASC
            """

            job = Gaia.launch_job_async(query)
            table = job.get_results()
            table['flux'] = [val * exp_time for val in table['phot_g_mean_flux']]
            table['mag'] = table['phot_g_mean_mag'] 

            #display(table)
            return table
        
        else: # offline query, does not have max table length limmit like astroquery does 
            if catalog == "H":
             # try hipparchos query with hipparchos file: 
                try: 
                    table = cls.hipparchos_query(cat_file, RA, Dec, radius, mag_lim)
                    table['flux'] = [val * exp_time for val in table['v_flux']]
                    table['mag'] = table['vmag']
                    return table
                except Exception as e:
                    print(f"Did you input a hipparchos catalog?")
                    return []
            elif catalog == "G":
                try:
                    with gaiaoffline.Gaia(magnitude_limit = (-3, mag_lim), photometry_output='mag') as gaia: 
                        table = gaia.conesearch(RA, Dec, radius)

                        table['flux'] = [cls.g_mag_to_flux(mag, exp_time) for mag in table['phot_g_mean_mag']]
                        table['mag'] = table['phot_g_mean_mag']

                        #table['flux'] = [val * exp_time for val in table['phot_g_mean_flux']] 
                        # results do not include mag, so convert from flux and known GAIA zero point: 
                        #table['mag'] = [cls.g_mag_from_flux(flux, exp_time) for flux in table['phot_g_mean_flux']]

                        sorted_table = table.sort_values(by = 'mag', ascending=True)
                        astropy_table = Table.from_pandas(sorted_table[:100])

                        #display(astropy_table)
                        
                    return astropy_table
                except Exception as e: 
                    print(f"Do you have a gaia offline catalog downloaded?: expection: {e}")
                    return []
            elif catalog == "N":
                try:
                    table = cls.nomad_query(RA, Dec, radius, mag_lim, cat_file)
                    return table
                except:
                    print(f"Did you input a NOMAD catalog file to use?")
                    return []
            else:
                print("To use an offline catalog, please specify type: G for GAIA, H for Hipparchos, and N for NOMAD can be handled")
                print("For H or N, the downloaded catalog file must also be passed")
                return []

    @staticmethod
    def g_mag_to_flux(mag, exp_time):
        """Function to convert a given magnitude to flux assuming magnitude is in GAIA's G band
        Used to ensure stars and satellite are on the same scale in terms of flux
        Inputs: 
        - object magnitude (either collected from GAIA database or input as RSO value
        Output:
        - flux of object
        """
        flux = 10 ** ((mag - 25.6884) / (-2.5))   # zero point from GAIA DR2 release: 25.6884
        return flux * exp_time 

    @staticmethod
    def fwhm_from_mag(rso_mag):
        """Function to change fwhm for the RSO based on its magnitude, 
        so that stars and RSOs have varying sizes in generated images.
        These give an approximation of the affects of brightness on the size of stars in images
        but will not show stars that oversaturate images, nor show affects of poor seeing conditions. 
        
        Input: 
        - magnitude
        Output: 
        - fwhm in pixels 
        """ 
        if rso_mag >= 12: 
            fwhm = 2
        elif rso_mag >= 9: 
            fwhm = 2.5
        elif rso_mag >= 7: 
            fwhm = 3
        elif rso_mag >= 4:
            fwhm = 3.5
        else: 
            fwhm = 4

        return fwhm

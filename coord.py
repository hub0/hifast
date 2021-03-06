import pandas as pd
import numpy as np
import xlrd

from astropy.coordinates import SkyCoord, AltAz, ICRS
from astropy.time import Time, TimeDelta
import astropy.units as u

from .fast import FAST

def load_ky(ky_file, start_time=None, end_time=None):
    '''
    Load KY .xlsx file

    Parameters
    ----------
    ky_file : str
        The file name of the receiver cabin position file

    start_time : astropy.time.Time
        start time of the observation, to truncate the ky records

    end_time : astropy.time.Time
        end time of the observation, to truncate the ky records

    Returns
    -------
    time : astropy.time.Time
        time sequence of the receiver cabin position records 

    xyz : numpy.ndarray
        measured xyz position of the receiver cabin
    '''

    ky_sheetname = pd.ExcelFile(ky_file).sheet_names[1]
    ky_arr = pd.read_excel(ky_file, sheet_name=ky_sheetname)
    time = Time(list(np.squeeze(ky_arr[['SysTime']].to_numpy())), 
            format='iso', scale='utc') 
    time = time - TimeDelta(60*60*8, format='sec') # convert from GMT+8 to GMT
    xyz = ky_arr[['SwtDPos_X', 'SwtDPos_Y', 'SwtDPos_Z']].to_numpy()
    
    # truncate time and xyz by start_time and end_time
    if start_time != None:
        start_dt = time - start_time
        time = time[start_dt > 0]
        xyz = xyz[start_dt > 0]
    if end_time != None:
        end_dt = end_time - time
        time = time[end_dt > 0]
        xyz = xyz[end_dt > 0]

    return time, xyz
    
def xyz_to_AltAz(time, xyz):
    '''
    Convert XYZ position of the phase centre of the receiver to 
    (Azimuth, Elevation) pointing of the 1st beam.

    Definition of XYZ
    origin: the focus
    X: geographic east
    Y: geographic north
    Z: zenith

    Parameters
    ----------
    time : astropy.time.Time
        Time stamps of xyz records

    xyz : numpy.ndarray
        the 3D position with shape (, 3) of the receiver relative to the 
        focus of the telescope.

    Returns
    -------
    altaz : astropy.coordinates.AltAz
        The AltAz of the 1st beam
    '''
    x = xyz[:, 0]
    y = xyz[:, 1]
    z = xyz[:, 2]

    az = (np.pi * 1.5 - np.arctan2(y, x)) * u.rad 
    alt = (np.arctan((z**2/(x**2+y**2))**0.5)) * u.rad

    return AltAz(az=az, alt=alt, obstime=time, location=FAST.loc())
    
    
def xyz_to_radec(time, xyz):
    '''
    Convert XYZ position of the phase centre of the receiver to 
    the projected ICRS coordinates of the 1st beam.

    Definition of XYZ
    origin: the focus
    X: geographic east
    Y: geographic north
    Z: zenith

    Parameters
    ----------
    time : astropy.time.Time
        Time stamps of xyz records

    xyz : numpy.ndarray
        the 3D position with shape (, 3) of the receiver relative to the 
        focus of the telescope.

    Returns
    -------
    ra : float
        The ICRS R.A. of the phase center in deg

    dec : float
        The ICRS Dec of the phase center in deg
    '''
    altaz = xyz_to_AltAz(time, xyz)
    
    icrs = altaz.transform_to(ICRS)

    return icrs.ra, icrs.dec

from astropy.coordinates import EarthLocation, Angle
import astropy.units as u

class FAST:
    '''
    This class is dedicated for refering FAST parameters.

    loc : astropy.coordinates.EarthLocation
        The position of FAST
    '''
    _lon = Angle('106:51:24.0', unit=u.deg) # degree E
    _lat = Angle('25:39:10.6', unit=u.deg) # degree N
    _elv = 1110.0288 * u.m

    @staticmethod
    def loc():
        return EarthLocation.from_geodetic(
                lon=FAST._lon, lat=FAST._lat, height=FAST._elv)

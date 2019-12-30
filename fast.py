from astropy.coordinates import EarthLocation, Angle
import astropy.units as u

import numpy as np

class FAST:
    '''
    This class is dedicated for refering FAST parameters.

    loc : astropy.coordinates.EarthLocation
        The position of FAST
    '''
    _lon = Angle('106:51:24.0', unit=u.deg) # degree E
    _lat = Angle('25:39:10.6', unit=u.deg) # degree N
    _elv = 1110.0288 * u.m

    _beam_offset = np.array([       # in arcmin
        #  ra    ,  dec     ,          beam_num
        [  0     ,  0      ],       #  1
        [  5.74  ,  0.00811],       #  2
        [  2.88  , -4.97   ],       #  3
        [ -2.86  , -4.98   ],       #  4
        [ -5.74  , -0.0127 ],       #  5
        [ -2.88  ,  4.97   ],       #  6
        [  2.86  ,  4.98   ],       #  7
        [ 11.5   ,  0.0116 ],       #  8
        [  8.61  , -4.96   ],       #  9
        [  5.75  , -9.93   ],       # 10
        [  0.018 , -9.94   ],       # 11
        [ -5.71  , -9.96   ],       # 12
        [ -8.6   , -4.99   ],       # 13
        [-11.5   , -0.03   ],       # 14
        [ -8.63  ,  4.95   ],       # 15
        [ -5.77  ,  9.93   ],       # 16
        [ -0.0181,  9.94   ],       # 17
        [  5.73  ,  9.95   ],       # 18
        [  8.61  ,  4.98   ],       # 19
        # credit: Qian Lei @2019 FAST User Workshop
        ]) * u.arcmin

    @classmethod
    def loc(cls):
        return EarthLocation.from_geodetic(
                lon=cls._lon, lat=cls._lat, height=cls._elv)
    
    @classmethod
    def beam_offset(cls, beam, rot):
        '''
        Get beam offset for 19 beam configuration

               rot<---  N
                
                  18----17----16
                 / \   / \   / \
                /   \ /   \ /   \
               19----07----06----15
              / \   / \   / \   / \
             /   \ /   \ /   \ /   \
         E  08----02----01----05----14
             \   / \   / \   / \   /
         |    \ /   \ /   \ /   \ /
         |     09----03----04----13
         V      \   / \   / \   /
        rot      \ /   \ /   \ /
                  10----11----12

        Parameters
        ----------
        beam : int
            beam number
        
        rot : float
            rotation angle in rad

        Returns
        -------
        off_ra : astropy.units.Quantity 
            ra offset in rad unit 

        off_dec : astropy.units.Quantity
            dec offset in rad unit 
        '''
  
        c, s = np.cos(rot), np.sin(rot)
        x, y = cls._beam_offset[beam-1]
        off_ra = c * x - s * y 
        off_dec = s * x + c * y 
        return off_ra, off_dec

    @classmethod
    def mainbeam_eff(cls, beam, za):
        '''
        Calculate main beam efficiency for zenith angle < 26.4 deg
        
        Parameters
        ----------
        beam : int
            beam number

        za : float or numpy.ndarray
            zenith angle in degree

        Returns
        -------
        mb_eff : float or numpy.ndarray
            main beam efficiency with the same shape of za
        '''
        return
        
        
        


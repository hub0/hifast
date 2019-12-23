import glob
import os
import pickle

import numpy as np

from astropy.io import fits
from astropy.time import Time
import astropy.units as u
from .coord import * 

class Chart:
    '''
    Object indexing FAST HI data for calibration and reduction.

    Each attribute contains a numpy.array, each item of which records 
    a specified meta-data for one spectrum with spectrum stack in time
    series. 

    A Chart contains following meta-data for an observation.
        obj : name of the object
        filename : FITS file name
        index : the sequence number in the FITS of spec
        freq1ch : central frequency of the 1st channel of each spectrum
        chanbw : width of 1 channel
        nchan : number of channels
        time  : time stamp of each spectrum
        coord : sky coordinate of each spectrum

    Parameters
    ----------
    obj : str
        Name of the object

    beam : int
        Beam number. 1~19 for certain beam of the 19 beam receiver,
        0 for a beam-fusion chart.

    filename : numpy.ndarray 
        File name of the FITS which the spectrum is stored.

    index : numpy.ndarray
        The Index number of the spectrum in its FITS.

    freq1ch : astropy.units.Quantity (in frequency unit)
        Central frequency of the 1st channel of each spectrum

    chanbw : astropy.units.Quantity (in frequency unit)
        Width of 1 channel

    nchan : numpy.ndarray
        Number of channels

    npol : numpy.ndarray
        Number of polarizations
    
    time : astropy.time.Time
        Time stamp of each spectrum

    coord : astropy.coordinates.SkyCoord
        Sky coordinates for each spectra
    '''

    def __init__(self, obj, beam, 
            filename=None, index=None, 
            freq1ch, chanbw, nchan, npol,
            time, coord=None, data=None):
        self._obj = obj
        self._beam = beam
        self._filename = filename
        self._index = index.astype('int')
        self._freq1ch = freq1ch.astype('float')
        self._chanbw = chanbw.astype('float')
        self._nchan = nchan.astype('int')
        self._npol = npol.astype('int')
        self._time = time 
        self._coord = coord 
        self._data = data.astype('float') 

    @property
    def obj(self):
        '''Get the name of object'''
        return self._obj
    @obj.setter
    def obj(self, value):
        '''Set coord sequence'''
        self._obj = value

    @property
    def beam(self):
        '''Get the beam number'''
        return self._beam
    @beam.setter
    def beam(self, value):
        '''Set beam number'''
        self._beam = value

    @property
    def filename(self):
        '''Get FITS file name'''
        return self._filename
    @filename.setter
    def filename(self, value):
        '''Set FITS file name'''
        self._filename= value

    @property
    def index(self):
        '''Get the list of index sequence'''
        return self._index
    @index.setter
    def index(self, value):
        '''Set index sequence'''
        self._index = value

    @property
    def time(self):
        '''Get the list of time sequence'''
        return self._time
    @time.setter
    def time(self, value):
        '''Set time sequence'''
        self._time = value

    @property
    def coord(self):
        '''Get the list of ra dec sequence'''
        return self._coord
    @coord.setter
    def coord(self, value):
        '''Set coord sequence'''
        self._coord = value

    @property
    def data(self):
        '''Get the data'''
        return self._data
    @data.setter
    def data(self, value):
        '''Set data'''
        self._data = value

    @property
    def freq1ch(self):
        '''Get the list of the frequency of 1st channel'''
        return self._freq1ch
    @freq1ch.setter
    def freq1ch(self, value):
        '''Set freq1ch sequence'''
        self._freq1ch = value

    @property
    def chanbw(self):
        '''Get the list channel band-width'''
        return self._chanbw
    @chanbw.setter
    def chanbw(self, value):
        '''Set chanbw sequence'''
        self._chanbw = value

    @property
    def nchan(self):
        '''Get the list of channel number'''
        return self._nchan
    @nchan.setter
    def nchan(self, value):
        '''Set nchan sequence'''
        self._nchan = value 

    @property
    def npol(self):
        '''Get the list of pols number'''
        return self._npol
    @npol.setter
    def npol(self, value):
        '''Set npol'''
        self._npol = value

    @property
    def freq(self):
        '''Get the 2D numpy array of freq'''
        chan = np.array([np.arange(x) for x in self.nchan])
        chanbw = np.expand_dims(self.chanbw, axis=1)
        freq1ch = np.expand_dims(self.freq1ch, axis=1)
        return freq1ch + chan * chanbw

    @classmethod
    def create(cls, obj=None, beam=1, rot=0, fits_path=None, ky_file=None):
        '''
        Create a Chart instance from FITS and KY files

        Parameters
        ----------
        obj : str
            Name of the object.

        beam : int
            Beam number.

        rot : float
            Rotation angle in rad.
        
        fits_path : str
            Path of FITS files. FITS in under this path should have 
            names like CygXN_OTF-M01_N_0001.fits.

        ky_file : str
            File name of the KY file. This KY file should cover the
            time span of the whole observation.

        Returns
        -------
        chart : Chart
            A Chart instance
        '''
        filename_arr = np.empty(0, dtype='unicode')
        index_arr = np.array([])
        freq1ch_arr = np.array([])
        chanbw_arr = np.array([])
        nchan_arr = np.array([])
        npol_arr = np.array([])
        coord_list = []
        obs_time_list = []
        data = []
    
        beam_str = '{:02d}'.format(beam)
        fits_list = glob.glob(fits_path + '*M' + beam_str + '_N_*.fits')
        fits_list.sort()

        if len(fits_list) == 0:
            raise ValueError('FITS path contains no matching files')

        for fits_name in fits_list:
            nfits = fits.open(fits_name)

            nspec = nfits[1].header['NAXIS2']

            farr = np.empty(nspec, dtype=f'<U{len(fits_name)}')
            farr[:] = fits_name   
            filename_arr = np.append(filename_arr, farr) 

            index = np.arange(nspec, dtype='int')
            index_arr = np.concatenate((index_arr, index))

            freq1ch = nfits[1].data.field('FREQ')
            freq1ch_arr = np.concatenate((freq1ch_arr, freq1ch))

            chanbw = nfits[1].data.field('CHAN_BW')
            chanbw_arr = np.concatenate((chanbw_arr, chanbw))

            nchan = nfits[1].data.field('NCHAN')
            nchan_arr = np.concatenate((nchan_arr, nchan))

            npol = np.ones(nspec) * nfits[1].data.field('DATA').shape[-1]
            npol_arr = np.concatenate((npol_arr, npol))

            time = list(nfits[1].data.field('DATE-OBS'))
            obs_time_list += time

            for spec in nfits[1].data.field('DATA'):
                data.append(spec)
        
        freq1ch_arr = freq1ch_arr * u.MHz
        chanbw_arr = chanbw_arr * u.MHz
        obs_time_arr = Time(obs_time_list, format='isot', scale='utc')
        data = np.array(data)

        # calculate obs_coord from obs_time and ky_xyz
        #
        if ky_file:
            ky_time, ky_xyz = load_ky(ky_file)
            obs_xyz_list = []
            for obs_time in obs_time_arr:
                obs_dt = ky_time - obs_time
                obs_dt0 = obs_dt[obs_dt<0][-1].sec
                obs_dt1 = obs_dt[obs_dt>0][0].sec
                ky_xyz0 = ky_xyz[obs_dt<0][-1]
                ky_xyz1 = ky_xyz[obs_dt>0][0]
                weight0 = abs(obs_dt0)/(abs(obs_dt0)+abs(obs_dt1))
                weight1 = abs(obs_dt1)/(abs(obs_dt0)+abs(obs_dt1))
                obs_xyz_list.append(ky_xyz0*weight0 + ky_xyz1*weight1)
            obs_xyz_arr = np.array(obs_xyz_list)
            c_ra, c_dec = xyz_to_radec(obs_time_arr, obs_xyz_arr)
            # cal offset for beam and rotation
            off_ra, off_dec = FAST.beam_offset(beam, rot)
            obs_dec = c_dec + off_dec
            obs_ra = c_ra + off_ra / np.cos(obs_dec.rad)
            obs_coord = SkyCoord(ra=obs_ra, dec=obs_dec, 
                    frame='icrs', obstime=obs_time_arr)
        else:
            obs_coord = None

        return cls(obj, beam, 
                filename_arr, index_arr, 
                freq1ch_arr, chanbw_arr, nchan_arr, 
                npol_arr, 
                obs_time_arr, obs_coord,
                data)

    def save(self, pkl_name):
        '''
        save Chart to a binary file
        '''
        if pkl_name[-6:].lower() != '.chart':
            pkl_name += '.chart'

        with open(pkl_name, 'wb') as f:
            pickle.dump(self, f)
        return

    @classmethod
    def load(cls, pkl_name):
        '''
        Load Chart from a pickle file
        '''
        f = open(pkl_name, 'rb')
        return pickle.load(f)

    def __getitem__(self, idx):
        sliced = Chart(self.obj, self.beam, 
                self.filename[idx], self.index[idx],
                self.freq1ch[idx], self.chanbw[idx], self.nchan[idx],
                self.npol[idx], self.time[idx], self.coord[idx])
        return sliced

    def __len__(self):
        return len(self.filename)

    def checkout_rawdata(self):
        raw = np.empty((0, self.nchan[0], self.npol[0]))
        for (fn, ind) in zip(self.filename, self.index):
            hdu_list = fits.open(fn)
            data = hdu_list[1].data.field('DATA')[ind]
            data = np.expand_dims(data, axis=0)
            raw = np.append(raw, data, axis=0)
        self.data = raw
        return 

    def freq_trim(self, freq_range):
        # todo: upgrade to freq_regrid method with freq_range, chan_bw,
        # regrid mode
        freq_lower = freq_range[0]
        freq_upper = freq_range[1]
        
        trmd_data = []
        trmd_freq1ch = []
        trmd_nchan = []
        for (freq1ch, nchan, chanbw, spec) in zip(self.freq1ch, 
                self.nchan, self.chanbw, self.data):
            origin_freq = np.arange(nchan)*chanbw + freq1ch
            mask = (origin_freq > freq_lower) & (origin_freq < freq_upper)
            trmd_freq = origin_freq[mask]
            trmd_freq1ch.append(trmd_freq[0].to(u.MHz).value)
            trmd_nchan.append(len(trmd_freq))
            trmd_data.append(spec[mask])
        trmd_data = np.array(trmd_data)
        trmd_freq1ch = np.array(trmd_freq1ch)
        trmd_nchan = np.array(trmd_nchan)
        return Chart(self.obj, self.beam, self.filename, self.index,
                trmd_freq1ch, self.chanbw, trmd_nchan, self.npol, 
                self.time, self.coord, trmd_data)

    #def stokes_I(self):
    #    '''
    #    convert data from pols to stokes I
    #    '''
    def pols_trim(self):
        trmd = self
        trmd.data = trmd.data[:,:,:2]
        return trmd
       
            

            

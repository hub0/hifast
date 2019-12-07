import glob
import os
import json

import numpy as np

from astropy.io import fits
from astropy.time import Time
import astropy.units as u

class Chart:
    '''
    Object indexing FAST HI data for calibration and reduction.

    Each attribute contains a numpy.array, each item of which records 
    a specified meta-data for one spectrum with spectrum stack in time
    series. 

    A Chart contains following meta-data for an observation.
        obj : name of the object
        index : FITS file name and the sequence number in the FITS
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

    filename : list
        File name of the FITS which the spectrum is stored.

    index : numpy.ndarray
        The Index number of the spectrum in its FITS. 

    freq1ch : astropy.units.Quantity (in frequency unit)
        Central frequency of the 1st channel of each spectrum

    chanbw : astropy.units.Quantity (in frequency unit)
        Width of 1 channel

    nchan : numpy.ndarray
        Number of channels
    
    time : astropy.time.Time
        Time stamp of each spectrum

    coord : astropy.coordinates.SkyCoord
        Sky coordinates for each spectra
    '''

    def __init__(self, obj, beam, filename, index, freq1ch, chanbw, nchan, time, coord=None):
        self._obj = obj
        self._beam = beam
        self._filename = filename
        self._index = index
        self._freq1ch = freq1ch
        self._chanbw = chanbw
        self._nchan = nchan 
        self._time = time 
        self._coord = coord 

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
    def freq(self):
        '''Get the 2D numpy array of freq'''
        chan = np.array([np.arange(x) for x in self.nchan])
        chanbw = np.expand_dims(self.chanbw, axis=1)
        freq1ch = np.expand_dims(self.freq1ch, axis=1)
        return freq1ch + chan * chanbw

    @classmethod
    def create(cls, obj=None, beam=0, fits_path=None, ky_path=None):
        filename_list = []
        index_arr = np.array([])
        freq1ch_arr = np.array([])
        chanbw_arr = np.array([])
        nchan_arr = np.array([])
        coord_list = []
        time_list = []
    
        beam_str = '{:02d}'.format(beam)
        fits_list = glob.glob(fits_path + '*M' + beam_str + '_N_*.fits')
        fits_list.sort()

        if len(fits_list) == 0:
            raise ValueError('FITS path contains no matching files')

        for fits_name in fits_list:
            nfits = fits.open(fits_name)

            nspec = nfits[1].header['NAXIS2']

            filename_list += [fits_name] * nspec

            index = np.arange(nspec, dtype='int') + 1
            index_arr = np.concatenate((index_arr, index))

            freq1ch = nfits[1].data.field('FREQ')
            freq1ch_arr = np.concatenate((freq1ch_arr, freq1ch))

            chanbw = nfits[1].data.field('CHAN_BW')
            chanbw_arr = np.concatenate((chanbw_arr, chanbw))

            nchan = nfits[1].data.field('NCHAN')
            nchan_arr = np.concatenate((nchan_arr, nchan))

            time = list(nfits[1].data.field('DATE-OBS'))
            time_list += time
        
        freq1ch_arr = freq1ch_arr * u.MHz
        chanbw_arr = chanbw_arr * u.MHz
        time_list = Time(time_list, format='isot', scale='utc')

        return cls(obj, beam, filename_list, index_arr, freq1ch_arr, chanbw_arr, 
                nchan_arr, time_list)

    #def save_json(self, json_name):
    #    '''
    #    save Chart to a .json file
    #    '''
    #    if json_name[-5:].lower() != '.json':
    #        json_name += '.json'

    #    with open(json_name, 'w') as f:
    #        json.dump(vars(self), f, cls=NpEncoder, sort_keys=True, indent=4)
    #    return

    #@classmethod
    #def load_json(cls, json_name):
    #    '''
    #    load Chart from a .json file
    #    '''
    #    if json_name[-5:].lower() != '.json':
    #        json_name += '.json'

    #    with open(json_name, 'r') as f:
    #        chart_dict = json.load(f)
    #    return cls(chart_dict['_obj'], chart_dict['_index'], chart_dict['_freq1ch'],
    #            chart_dict['_chanbw'], chart_dict['_nchan'], chart_dict['_time'],
    #            chart_dict['_coord'])


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.int64):
            return int(obj)
        elif isinstance(obj, np.float64):
            return float(obj)
        else:
            return super(NpEncoder, self).default(obj)

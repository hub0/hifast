import glob
import os
import json

import numpy as np

from astropy.io import fits
from astropy.time import Time
import astropy.units as u

class Chart:
    '''
    Create a json chart of OTF observation for an object of one beam.
    A Chart contains following attributes.
        index : [(fits_name, seq_num), ...]
        freq  : [(freq_of_1st_chan, chan_bw, nchan), ...]
        time  : [obs_time_isot, ...]
        coord : [(ra, dec), ...]

    Parameters
    ----------
    obj : str
        Name of the object.

    index : list
        spectral locator in a format of [[fits_name, seqnum], ...]

    freq_setup : list
        frequency setup in a format of [[freq1ch, chan_bw, nchan], ...]

    coord : list
        

    beam_num : int
        The beam number (start from 1).

    fits_path : str
        Path of all fits file for the observation.

    ky_path : str
        Path of KY file.
    '''
    def __init__(self, obj, index, freq1ch, chanbw, nchan, time, coord):
        self._index = index
        self._freq1ch = freq1ch
        self._chanbw = chanbw
        self._nchan = nchan 
        self._time = time 
        self._coord = coord 
        self._obj = obj

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
    def obj(self):
        '''Get the name of object'''
        return self._obj
    @obj.setter
    def obj(self, value):
        '''Set coord sequence'''
        self._obj = value

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
        return np.array(
                [np.arange(x) * self._chanbw + self._freq1ch for x in self._nchan])

    @classmethod
    def create(cls, obj=None, beam_num=1, fits_path=None, ky_path=None):
        index_seq = []
        freq1ch_seq = []
        chanbw_seq = []
        nchan_seq = []
        coord_seq = []
        time_seq = []
    
        beam_str = '{:02d}'.format(beam_num)
        fits_list = glob.glob(fits_path + '*M' + beam_str + '_N_*.fits')
        fits_list.sort()

        if len(fits_list) == 0:
            raise ValueError('FITS path contains no matching files')

        for fits_name in fits_list:
            nfits = fits.open(fits_name)
            nchan =   list(nfits[1].data.field('NCHAN'))
            nchan_seq += nchan
            freq1ch = list(nfits[1].data.field('FREQ'))
            freq1ch_seq += freq1ch
            chanbw = list(nfits[1].data.field('CHAN_BW'))
            chanbw_seq += chanbw
            #freq = [list(x) for x in zip(freq1ch, chan_bw, nchan)]
            #freq_seq += freq
        
            time = list(nfits[1].data.field('DATE-OBS'))
            time_seq += time
        
            nspec = nfits[1].header['NAXIS2']
            index = [[fits_name, x+1] for x in range(nspec)]
            index_seq += index

        return cls(obj, index_seq, freq1ch_seq, chanbw_seq, 
                nchan_seq, time_seq, coord_seq)

    def save_json(self, json_name):
        '''
        save Chart to a .json file
        '''
        if json_name[-5:].lower() != '.json':
            json_name += '.json'

        with open(json_name, 'w') as f:
            json.dump(vars(self), f, cls=NpEncoder, sort_keys=True, indent=4)
        return

    @classmethod
    def load_json(cls, json_name):
        '''
        load Chart from a .json file
        '''
        if json_name[-5:].lower() != '.json':
            json_name += '.json'

        with open(json_name, 'r') as f:
            chart_dict = json.load(f)
        return cls(chart_dict['_obj'], chart_dict['_index'], chart_dict['_freq1ch'],
                chart_dict['_chanbw'], chart_dict['_nchan'], chart_dict['_time'],
                chart_dict['_coord'])


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.int64):
            return int(obj)
        elif isinstance(obj, np.float64):
            return float(obj)
        else:
            return super(NpEncoder, self).default(obj)

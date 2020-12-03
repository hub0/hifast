import numpy as np
import astropy.units as u
from astropy.units import UnitConversionError
from .chart import Chart

def get_tcal(tcal_fn, beam):
    '''
    Get T_cal from measurment files

    Parameters
    ----------
    tcal_name : str
        name of the T_cal measurement file
    bm_num : int
        beam number

    Returns
    -------
    tcal : Chart
    '''
    with open(tcal_fn, 'r') as tcal_f:
        measure = [x.strip() for x in tcal_f if x[0] != '#']

    # initialization
    length = int(len(measure)/8)
    freqs = np.empty(length)
    tcal = np.empty([19, length, 2])

    for i in range(len(measure)):
        index = int(i/8)
        items = measure[i].split()
        if i%8 == 0:
            freqs[index] = float(items[0])
            tcal[0, index, 0] = float(items[1])
            tcal[0, index, 1] = float(items[2])
            tcal[1, index, 0] = float(items[3])
            tcal[1, index, 1] = float(items[4])
        if i%8 == 1:
            tcal[2, index, 0] = float(items[0])
            tcal[2, index, 1] = float(items[1])
            tcal[3, index, 0] = float(items[2])
            tcal[3, index, 1] = float(items[3])
            tcal[4, index, 0] = float(items[4])
        if i%8 == 2:
            tcal[4, index, 1] = float(items[0])
            tcal[5, index, 0] = float(items[1])
            tcal[5, index, 1] = float(items[2])
            tcal[6, index, 0] = float(items[3])
            tcal[6, index, 1] = float(items[4])
        if i%8 == 3:
            tcal[7, index, 0] = float(items[0])
            tcal[7, index, 1] = float(items[1])
            tcal[8, index, 0] = float(items[2])
            tcal[8, index, 1] = float(items[3])
            tcal[9, index, 0] = float(items[4])
        if i%8 == 4:
            tcal[ 9, index, 1] = float(items[0])
            tcal[10, index, 0] = float(items[1])
            tcal[10, index, 1] = float(items[2])
            tcal[11, index, 0] = float(items[3])
            tcal[11, index, 1] = float(items[4])
        if i%8 == 5:
            tcal[12, index, 0] = float(items[0])
            tcal[12, index, 1] = float(items[1])
            tcal[13, index, 0] = float(items[2])
            tcal[13, index, 1] = float(items[3])
            tcal[14, index, 0] = float(items[4])
        if i%8 == 6:
            tcal[14, index, 1] = float(items[0])
            tcal[15, index, 0] = float(items[1])
            tcal[15, index, 1] = float(items[2])
            tcal[16, index, 0] = float(items[3])
            tcal[16, index, 1] = float(items[4])
        if i%8 == 7:
            tcal[17, index, 0] = float(items[0])
            tcal[17, index, 1] = float(items[1])
            tcal[18, index, 0] = float(items[2])
            tcal[18, index, 1] = float(items[3])

    freq1ch = np.array(freqs[0:1]) * u.MHz

    nchan = len(freqs)

    freqs0 = freqs[:-1]
    freqs1 = freqs[1:]
    chanbw_list = [x-y for (x, y) in zip(freqs1, freqs0)]
    chanbw = np.array([sum(chanbw_list) / len(chanbw_list)]) * u.MHz

    data = np.expand_dims(tcal[beam-1, :, :], axis=0) * u.K

    #return freq1ch, nchan, chanbw, data
    return Chart('tcal', beam, 
            freq1ch=freq1ch, chanbw=chanbw, nchan=nchan,
            npol=np.array([2]), data=data)

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string

    credit: https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html
    (modifed to keep the array length)
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if len(x) < window_len:
        raise ValueError("Input array needs to lager than window size on the given axis.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window must be one of 'flat', 'hanning', 'hamming', "
                "'bartlett', 'blackman'")


    s = np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]

    ###  1. transpose the ndarray to move the given axis to the last dimension;
    ###  2. Create left and right reflective boundaries on the last dimentsion
    ###     using ellipsis;
    ###  3. concatenate boundaries with the ndarray
    ###  4. convolve on the last dimension
    ###  5. trim the convolved ndarray on the last dimension
    ###  6. transpose the convolved ndarray to move the last dimension to the 
    ###     axis-th dimension

    ##axes = [n for n in range(x.ndim)]
    ##axes.append(axes.pop(axis))
    ##s = np.transpose(x, axes=axes)
    ##
    ##lb = s[..., window_len-1:0:-1]
    ##rb = s[..., -2:-window_len-1:-1]
    ##
    ##s = np.concatenate((lb, s, rb), axis=-1)

    if window == 'flat': #moving average
        w = np.ones(window_len,'d')
    else:
        w = eval('np.'+window+'(window_len)')

    y = np.convolve(w/w.sum(),s,mode='valid')
    
    y = y[int(window_len/2):-1*int(window_len/2)] # trim to the same size as x
    y = y / np.mean(y) * np.mean(x)  # preserve the total value
    return y

def freq_smooth(c, freq_win, window='hanning'):
    """
    Smooth the data in Chart in frequency regime with given frequency widow.

    Parameters
    ----------
    c: hifast.Chart 
        The chart to be smoothed
    freq_win : astropy.units.quantity.Quantity
        The width of smooth window in frequency unit
    window : str
        The type of smooth window

    Returns
    -------
    smoothed_chart : hifast.Chart
    """

    try:
        freq_win.to(u.Hz)
    except UnitConversionError:
        raise UnitConversionError(
                'freq_win must be with a frequency equivalent unit'
                )

    chan_win = int(freq_win / (c.freq[1] - c.freq[0]))
    
    

    for i in range(c.data.shape[0]):
        for j in range(c.data.shape[2]):
            spec = smooth(c.data[i, :, j], window_len=chan_win, window=window)
            c.data[i, :, j] = spec

    return c

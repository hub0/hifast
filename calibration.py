import numpy as np
import astropy.units as u
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

    nchan = np.array([len(freqs)])

    freqs0 = freqs[:-1]
    freqs1 = freqs[1:]
    chanbw_list = [x-y for (x, y) in zip(freqs1, freqs0)]
    chanbw = np.array([sum(chanbw_list) / len(chanbw_list)]) * u.MHz

    data = np.expand_dims(tcal[beam-1, :, :], axis=0) * u.K

    #return freq1ch, nchan, chanbw, data
    return Chart('tcal', beam, 
            freq1ch=freq1ch, chanbw=chanbw, nchan=nchan,
            npol=np.array([2]), data=data)

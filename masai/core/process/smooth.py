
import numpy as np
import pandas as pd

def smooth1D(p, window_len=11, window='hanning'):
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
    
    np.hanning, np.hamming, np.bartlett, np.blackman, np.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    if isinstance(p, pd.Series):
        x = p.values
        idx = p.index
    else:
        x = p

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len < 3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s = np.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]
    # print(len(s), len(x))
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='valid')
    y = y[(window_len / 2):-(window_len / 2)]

    if isinstance(p, pd.Series):
        return pd.Series(y, index=idx)
    else:
        return y

def smooth(p, window_len=11, window='hanning', axis=-1):
    """
    """

    window_len = int(window_len/2)*2-1

    if isinstance(p, pd.DataFrame):
        x = p.values

    else:
        x = p

    pt = np.apply_along_axis(smooth1D, axis, p, window_len, window)

    if isinstance(p, pd.DataFrame):
        return pd.DataFrame(pt, index=p.index, columns=p.columns)
    else:
        return pt

smooth2D = smooth

if __name__ == '__main__':
    pass
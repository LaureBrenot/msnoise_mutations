"""
This plot shows the final output of MSNoise.


.. include:: ../clickhelp/msnoise-plot-dvv.rst


Example:

``msnoise plot dvv`` will plot all defaults:

.. image:: ../.static/dvv.png

"""

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

from matplotlib.dates import DateFormatter

from ..api import *

import numpy as np
from scipy.stats import linregress
import pandas as pd

def compute_slopes(x, y, step_size: int) :
    """
    Computes the slopes between adjacent points in the given data with a given step size.

    Args:
        x (List[float]): The x values of the data points.
        y (List[float]): The y values of the data points.
        step_size (int): The step size for computing slopes.

    Returns:
        List[float]: A list of slopes computed for adjacent points with the given step size.
    """
    slopes = []
    for i in range(0, len(x)-step_size, step_size):
        slope = (y[i+step_size] - y[i]) / (x[i+step_size] - x[i])
        slopes.append(slope)
    return slopes

def detect_drastic_slopes(x, y, threshold: float, step: int) :
    """
    Detects the indices of the points in the given data where the slope between that point and the previous point
    is greater than the given threshold.

    Args:
        x (List[float]): The x values of the data points.
        y (List[float]): The y values of the data points.
        threshold (float): The threshold for the slope difference between two points.

    Returns:
        List[Tuple[int, int]]: A list of tuples representing the indices of the points where a drastic slope change occurs.
        Each tuple contains two integers: the index of the point where the slope change starts, and the index of the point
        where the slope change ends (inclusive).
    """
    drastic_slopes = []
    prev_slope = None
    start_idx = None
    for i in range(1, len(x)):
        slope = (y[i] - y[i-step]) / (x[i] - x[i-step])
        if prev_slope is not None and abs(slope) > threshold:
            if start_idx is None:
                start_idx = i - step
            end_idx = i
        else:
            if start_idx is not None:
                drastic_slopes.append((start_idx, end_idx))
                start_idx = None
        prev_slope = slope
    if start_idx is not None:
        drastic_slopes.append((start_idx, len(x)-step))

    i = 0
    while i < len(drastic_slopes) - 1:
        if drastic_slopes[i][1] >= drastic_slopes[i+1][0]:
            drastic_slopes[i] = (drastic_slopes[i][0], max(drastic_slopes[i][1], drastic_slopes[i+1][1]))
            del drastic_slopes[i+1]
        else:
            i += 1

    return drastic_slopes

def wavg(group, dttname, errname):
    d = group[dttname]
    group[errname][group[errname] == 0] = 1e-6
    w = 1. / group[errname]
    wavg = (d * w).sum() / w.sum()
    return wavg


def wstd(group, dttname, errname):
    d = group[dttname]
    group[errname][group[errname] == 0] = 1e-6
    w = 1. / group[errname]
    wavg = (d * w).sum() / w.sum()
    N = len(np.nonzero(w)[0])
    wstd = np.sqrt(np.sum(w * (d - wavg) ** 2) / ((N - 1) * np.sum(w) / N))
    return wstd


def get_wavgwstd(data, dttname, errname):
    grouped = data.groupby(level=0)
    g = grouped.apply(wavg, dttname=dttname, errname=errname)
    h = grouped.apply(wstd, dttname=dttname, errname=errname)
    return g, h


def main(mov_stack=None, dttname="M", components='ZZ', filterid=1,
         pairs=[], showALL=False, show=False, outfile=None):
    db = connect()

    start, end, datelist = build_movstack_datelist(db)

    if mov_stack != 0:
        mov_stacks = [mov_stack, ]
    else:
        mov_stack = get_config(db, "mov_stack")
        if mov_stack.count(',') == 0:
            mov_stacks = [int(mov_stack), ]
        else:
            mov_stacks = [int(mi) for mi in mov_stack.split(',')]

    if components.count(","):
        components = components.split(",")
    else:
        components = [components, ]

    low = high = 0.0
    for filterdb in get_filters(db, all=True):
        if filterid == filterdb.ref:
            low = float(filterdb.low)
            high = float(filterdb.high)
            break

    gs = gridspec.GridSpec(len(mov_stacks), 1)
    fig = plt.figure(figsize=(12, 9))
    plt.subplots_adjust(bottom=0.06, hspace=0.3)
    first_plot = True
    for i, mov_stack in enumerate(mov_stacks):
        current = start
        first = True
        alldf = []
        while current <= end:
            for comp in components:
                day = os.path.join('DTT', "%02i" % filterid, "%03i_DAYS" %
                                   mov_stack, comp, '%s.txt' % current)
                if os.path.isfile(day):
                    df = pd.read_csv(day, header=0, index_col=0,
                                     parse_dates=True)
                    alldf.append(df)
            current += datetime.timedelta(days=1)
        if len(alldf) == 0:
            print("No Data for %s m%i f%i" % (components, mov_stack, filterid))
            continue

        alldf = pd.concat(alldf)
        print(mov_stack, alldf.head())
        if 'alldf' in locals():
            errname = "E" + dttname
            alldf.to_csv("tt.csv")
            alldf[dttname] *= -100
            alldf[errname] *= -100

            ALL = alldf[alldf['Pairs'] == 'ALL'].copy()
            allbut = alldf[alldf['Pairs'] != 'ALL'].copy()

            if first_plot == 1:
                ax = plt.subplot(gs[i])
            else:
                plt.subplot(gs[i], sharex=ax)

            for pair in pairs:
                print(pair)
                pair1 = alldf[alldf['Pairs'] == pair].copy()
                print(pair1.head())
                plt.plot(pair1.index, pair1[dttname], label=pair)
                plt.fill_between(pair1.index, pair1[dttname]-pair1[errname],
                                 pair1[dttname]+pair1[errname], zorder=-1,
                                 alpha=0.5)
                pair1.to_csv('%s-m%i-f%i.csv'%(pair, mov_stack, filterid))

            if showALL:
                plt.plot(ALL.index, ALL[dttname], c='r',
                         label='ALL: $\delta v/v$ of the mean network')

            tmp2 = allbut[dttname].resample('D').median()
            #print(tmp2)

            y = tmp2
            x = np.asarray(range(len(y)))
            
            slope_on = 30 #days
            polyfit = False
            if polyfit:
                data = [d for d in y if not pd.isnull(d)]
                x2 = np.asarray(range(len(data)))
                model5 = np.poly1d(np.polyfit(x2, data, 15))
                y = model5(x2)
                x = np.asarray(range(len(y)))

            # Calculate a threshold value for significant slope changes
            slopes_1 = compute_slopes(x, y, slope_on)
            threshold= np.nanmean(slopes_1)+ 3*np.nanstd(slopes_1)
            
            #plt.hist(slopes_1)
            #plt.vlines(threshold,0, 6)
            #plt.vlines(-threshold,0, 6)
            #plt.show()

            drastic_slopes = detect_drastic_slopes(x, y, threshold, slope_on)
            print(drastic_slopes)
            plt.plot(tmp2.index, y, color='blue')

            
     
            for start_idx, end_idx in drastic_slopes:
                plt.axvspan(tmp2.index[start_idx], tmp2.index[end_idx], color='red', alpha=0.5)
                print(datelist[start_idx],'to', datelist[end_idx])
                #plt.show()

                # df= pd.DataFrame(tmp2)
                # df.to_csv('dvvnumber.csv')

                if first_plot == 1:
                    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=4,
                               ncol=2, borderaxespad=0.)
                    left, right = tmp2.index[0], tmp2.index[-1]
                    if mov_stack == 1:
                        plt.title('1 Day')
                    else:
                        plt.title('%i Days Moving Window' % mov_stack)
                    first_plot = False
                else:
                    plt.xlim(left, right)
                    plt.title('%i Days Moving Window' % mov_stack)

                plt.grid(True)
                plt.gca().xaxis.set_major_formatter(DateFormatter("%Y-%m-%d %H:%M"))
                fig.autofmt_xdate()
                title = '%s, Filter %d (%.2f - %.2f Hz)' % \
                        (",".join(components), filterid, low, high)
                plt.suptitle(title)
                #del alldf
        if outfile:
            if outfile.startswith("?"):
                if len(mov_stacks) == 1:
                    outfile = outfile.replace('?', '%s-f%i-m%i-M%s' % (components,
                                                                   filterid,
                                                                   mov_stack,
                                                                   dttname))
                else:
                    outfile = outfile.replace('?', '%s-f%i-M%s' % (components,
                                                                   filterid,
                                                                   dttname))
            outfile = "dvv " + outfile
            print("output to:", outfile)
            plt.savefig(outfile)
    if show:
        plt.show()


if __name__ == "__main__":
    main()

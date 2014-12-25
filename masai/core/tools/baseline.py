import numpy as np
from matplotlib import pyplot as pl
from scipy.interpolate import pchip

from masai.process import smooth


def compute_baseline(pc, regions, scores=None, smoothpoints=3, debug=False):
    """
    Compute the baseline using SVD decomposition
    
    pc is a pd.DataFrame or pd.Series
    
    """
    # initialization
    # the data are PCs
    
    nc, npc = pc.shape
    w = pc.index
    # union of all base_regions set's
    if isinstance(regions[0][0],list):
        # a region for each PC
        multiple=True
        basepoints=[]
        for ipc in range(npc):
            bpoints = set([])
            for item in regions[ipc]:
                l = sorted(item)
                br = w[(w >= l[0]) & (w <= l[-1])]
                bpoints |= set(br.values)
            # be sure that the limits of the spectra are added
            bpoints |= set((w[0], w[1], w[-2], w[-1]))
            basepoints.append(bpoints)
    else:    
        multiple=False
        bpoints = set([])
        for item in regions:
            l = sorted(item)
            br = w[(w >= l[0]) & (w <= l[-1])]
            bpoints |= set(br.values)
        # be sure that the limits of the spectra are added
        bpoints |= set((w[0], w[1], w[-2], w[-1]))
        basepoints = [bpoints for i in range(npc)]
        
    # select only the necessary points for computing the baseline
    pcb = pc.copy()*0.0
    for ipc in xrange(npc):
        wb = np.array(sorted(basepoints[ipc] & set(w.values)))  # intersection
        pcs = pc.iloc[:,ipc].reindex(index=wb)
    
        # compute baseline
        wf = np.array(list(w))
        s = pchip(wb, pcs)
        bf = s(wf)
        #pcb = SVDdenoise(pcb, orda=250, rank=10, threshold=95., remove=1, axis=0, n_iter=1,
        #                 verbose='Full')
        pcb.iloc[:,ipc] = smooth(bf, window_len=2*int((smoothpoints+1)/2)+1, window="hanning")
             
    if debug :
        lim = (w[-1],w[0])
        for ipc in xrange(npc):
            nplot=1
            if scores is not None:
                nplot=2
            ax = pl.figure(figsize=(6.8,1.7)).add_subplot(1,nplot,1)
            pc.iloc[:,ipc].plot(ax=ax, xlim=lim, lw=0.5)
            ax.set_xlabel('wavenumbers (cm$^{-1}$)');
            ax.set_ylabel('a.u.');
            ax.set_title('PC %d'%ipc)
            ax.plot(wf,pcb.iloc[:,ipc])
            if multiple:
                region = regions[ipc]
                for reg in region:
                    ax.axvspan(*reg, color='y', alpha=.7)
            if scores is not None:
                ax2=subplot(1,2,2) 
                scores.iloc[:,ipc].plot(ax=ax2, lw=0.5);
            ax.xaxis.set_label_coords(0.5,-0.15)
            ax.yaxis.set_label_coords(-0.12,0.5)
            pl.subplots_adjust(bottom=0.22, left=0.15, right=0.95)
    
    return pcb
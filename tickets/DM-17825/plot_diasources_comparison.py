# This file can be imported into a jupyter notebook namespace by
# %run -n -i "plot_calexp_template_diffim.py"
# also has a functionality as a standalone script to produce reproducible figures
# =========


import os
import numpy as np
import pandas as pd
import sqlite3
import lsst.daf.persistence as dafPersist
import lsst.geom
from astropy.visualization import (ZScaleInterval, AsinhStretch, ImageNormalize)
import astropy.units as u
import matplotlib.pyplot as plt

from astropy.table import Table
from matplotlib.backends.backend_pdf import PdfPages


def loadAllPpdbObjects(repo, dbName='association.db'):
    """Load select DIAObject columns from a PPDB into a pandas dataframe.

    Parameters
    ----------
    repo : `str`
        Path to an output repository from an ap_pipe run.
    dbName : `str`, optional
        Name of the PPDB, which must reside in (or relative to) repo.

    Returns
    -------
    objTable : `pandas.DataFrame`
        DIA Object Table containing only objects with validityEnd NULL.
        Columns selected are presently hard-wired here.
    """
    connection = sqlite3.connect(os.path.join(repo, dbName))

    # These are the tables available in the ppdb
    tables = {'obj': 'DiaObject', 'src': 'DiaSource', 'ccd': 'CcdVisit'}

    # Only get objects with validityEnd NULL because that means they are still valid
    objTable = pd.read_sql_query('select diaObjectId, ra, decl, nDiaSources, '
    'gPSFluxMean, validityEnd, flags from {0} where validityEnd is NULL;'.format(tables['obj']), connection)
     
    return objTable
    
# ---
def plot_diasources_diffim(repo1, repo2,  srcTable1, srcTable2, diffimType='deepDiff_differenceExp',pdfWriter=None):
    # Upper row repo1 diffims with srcTable1 detections
    # Lower row repo2 diffims with srcTable2 detections
    # Cutout is controlled by repo1
    # Cutout is miniregion
    
    ccdVisitIds = np.unique(srcTable1['ccdVisitId'])
    
    n_plots_per_fig = 4
    n_figs = len(ccdVisitIds) // 4
    if len(ccdVisitIds) % 4 !=0:
        n_figs += 1
    fig_idx = 1 # Figure index for title
    
    fig = None
    ccdVisitIds = list(ccdVisitIds[::-1])
    n_panels = len(ccdVisitIds)

    # Center of the mini region
    centerMini = lsst.geom.SpherePoint(155.25, -5.7, lsst.geom.degrees)
    size = lsst.geom.Extent2I(100, 100)

#    panel_idx = 1 # Panel index (all subplots across all figures that belong to one sci object)

    butler1 = dafPersist.Butler(repo1)
#    butler2 = dafPersist.Butler(repo2)
    while len(ccdVisitIds) > 0:
        if fig is None:
            fig = plt.figure(figsize=(11,6))            
            fig.subplots_adjust(left=0.05, right=0.98, bottom=0.05, wspace=0.15)
#            fig.suptitle('DIAObject {}; {}/{}'.format(obj,fig_idx,n_figs))
            fig_idx += 1  
            splot_idx = 1 # Subplot index within figure (first row only)
        
        ccdVisitId = ccdVisitIds.pop()
        # Upper row calexp
        
        ax = fig.add_subplot(2,n_plots_per_fig,splot_idx)
        T = srcTable1[srcTable1['ccdVisitId']==ccdVisitId]

        visit = ccdVisitId // 100
        ccd = ccdVisitId % 100
        dataId = { 'visit': visit, 'ccdnum': ccd }
        print (visit)
        print (ccd)
        print (T['x'])
        print (T['y'])
        print (T['ra'])
        print (T['decl'])
        
        diffexp = butler1.get(diffimType, dataId)
        ax.set_title('{visit} {ccdnum:02d}'.format(**dataId))
        diffexpCutout = diffexp.getCutout(centerMini, size);
        bbox = diffexpCutout.getBBox()
        extentR = (bbox.getMaxY()+0.5, bbox.getMinY()-0.5, bbox.getMaxX()+0.5, bbox.getMinX()-0.5)
        diffexpArray = diffexpCutout.getMaskedImage().getImage().getArray()
        diffexpNorm = ImageNormalize(diffexpArray, interval=ZScaleInterval(), stretch=AsinhStretch())
        ax.imshow(diffexpArray.T[::-1,::-1], origin='lower', cmap='gray', norm=calexpNorm, extent = extentR)
        ax.grid(True)
        ax.scatter(T['y'],T['x'],s=4, alpha=0.5)
        
        
#        T = srcTable2[srcTable2['ccdVisitId']==ccdVisitId]
#        ax.scatter(T['ra'],T['decl'],s=4, alpha=0.5)
        ax.set_title('{visit} {ccd:02d}; {idx}/{n_panels}'.format(visit=visit, ccd= ccd, idx=panel_idx, n_panels=n_panels))
         
        panel_idx += 1
        splot_idx += 1
        
        if splot_idx > n_plots_per_fig:
            if pdfWriter is not None:
                pdfWriter.savefig(fig)
                plt.close(fig)
            fig = None
             
    if fig is not None and pdfWriter is not None:
        pdfWriter.savefig(fig)
        plt.close(fig)

    
# ------

def plot_diasources_compare(srcTable1, srcTable2, pdfWriter=None):
    """Plot the DiaSources as their RA and DEC on a scatter plot from `srcTable1` and `srcTable2` for comparison
    """
    
    ccdVisitIds = np.unique(srcTable1['ccdVisitId'])
    
    n_plots_per_fig = 4
    n_figs = len(ccdVisitIds) // 4
    if len(ccdVisitIds) % 4 !=0:
        n_figs += 1
    fig_idx = 1 # Figure index for title
    
    fig = None
    ccdVisitIds = list(ccdVisitIds[::-1])
    n_panels = len(ccdVisitIds)
    panel_idx = 1 # Panel index (all subplots across all figures that belong to one sci object)
    
    while len(ccdVisitIds) > 0:
        if fig is None:
            fig = plt.figure(figsize=(11,3.3))            
            fig.subplots_adjust(left=0.05, right=0.98, bottom=0.1, wspace=0.15)
#            fig.suptitle('DIAObject {}; {}/{}'.format(obj,fig_idx,n_figs))
            fig_idx += 1  
            splot_idx = 1 # Subplot index within figure (first row only)
        
        ccdVisitId = ccdVisitIds.pop()
        # Upper row calexp
        ax = fig.add_subplot(1,n_plots_per_fig,splot_idx)
        T = srcTable1[srcTable1['ccdVisitId']==ccdVisitId]
        ax.scatter(T['ra'],T['decl'],s=4, alpha=0.5)
        
        T = srcTable2[srcTable2['ccdVisitId']==ccdVisitId]
        ax.scatter(T['ra'],T['decl'],s=4, alpha=0.5)
        ax.set_title('{visitccd}; {idx}/{n_panels}'.format(visitccd=ccdVisitId, idx=panel_idx, n_panels=n_panels))
         

        panel_idx += 1
        splot_idx += 1
        
        if splot_idx > n_plots_per_fig:
            if pdfWriter is not None:
                pdfWriter.savefig(fig)
                plt.close(fig)
            fig = None
             
    if fig is not None and pdfWriter is not None:
        pdfWriter.savefig(fig)
        plt.close(fig)
        


# ==============

def main():
    import matplotlib
    matplotlib.use('Qt5Agg')

    cwpRepo = '/home/gkovacs/data/repo_DM-17825/ingested/rerun/proc_2019-02-21'
    cwpTemplateRepo = '/home/gkovacs/data/repo_DM-17825/templates'
#    my_dbName = '/home/gkovacs/data/repo_DM-17825/ingested/rerun/proc_2019-02-21/association.db'
#    mrawls_dbName = '/home/gkovacs/data/repo_DM-17825/mrawls_cw_processed2/association.db'

#    butlerCwp = dafPersist.Butler(cwpRepo)
    butlerCwpTemplate = dafPersist.Butler(cwpTemplateRepo)

    patchList = ['10,8', '11,8', '12,8', '13,8',
             '10,7', '11,7', '12,7', '13,7',
             '10,9', '11,9', '12,9', '13,9',
             '10,5', '11,5', '12,5', '13,5',
             '10,6', '11,6', '12,6', '13,6',
             '10,10', '11,10', '12,10', '13,10']
             
    cwpObjTable = loadAllPpdbObjects(cwpRepo)
    cwpMiniRegion = defMiniRegion(cwpObjTable)
    cwpMiniUnflagged = cwpMiniRegion & (cwpObjTable['flags'] == 0)
    cwpObjList = list(cwpObjTable.loc[cwpMiniRegion, 'diaObjectId'])
    cwpObjList.sort()

    # Find the patch that belongs to the mini region
    patch = patchFinder(cwpObjTable.loc[cwpMiniUnflagged,'diaObjectId'].values[0],cwpObjTable,butlerCwpTemplate,patchList)

    with PdfPages('proc_2019-02-21_diffims_mini.pdf') as W:
        for obj in cwpObjList:
            print(obj)
            plot_images(cwpRepo,cwpTemplateRepo,obj,patch,cwpObjTable,plotAllCutouts=True,pdfWriter=W)

    cwpObjList = list(cwpObjTable.loc[cwpMiniUnflagged, 'diaObjectId'])
    cwpObjList.sort()
    with PdfPages('proc_2019-02-21_diffims_mini_unflagged.pdf') as W:
        for obj in cwpObjList:
            print(obj)
            plot_images(cwpRepo,cwpTemplateRepo,obj,patch,cwpObjTable,plotAllCutouts=True,pdfWriter=W)

if __name__ == "__main__":
    main()

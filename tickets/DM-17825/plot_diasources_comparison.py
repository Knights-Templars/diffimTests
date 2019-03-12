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
    
def make_and_plot_one_cutout(ax, diffexp, Tsources, color=None, cmap = 'gray'):
    # We need to determine the center of the miniregion and the exposure as 
    # cutout does not support cutting where the center point is out of the image
    p1 = lsst.geom.SpherePoint(155.3, -5.8, lsst.geom.degrees)
    p2 = lsst.geom.SpherePoint(155.2, -5.6, lsst.geom.degrees)
    wcs = diffexp.getWcs()
    img_point1 = wcs.skyToPixel(p1)
    img_point2 = wcs.skyToPixel(p2)
    img_bbox = lsst.geom.Box2D(img_point1,img_point2)
    plt_bbox = lsst.geom.Box2I(img_bbox)
    img_bbox.clip(lsst.geom.Box2D(diffexp.getBBox()))
    center_point = wcs.pixelToSky(
        0.5 * (img_bbox.getMinX()+img_bbox.getMaxX()), 
        0.5 * (img_bbox.getMinY()+img_bbox.getMaxY()))
    size = lsst.geom.Extent2I(img_bbox.getMaxX() - img_bbox.getMinX() + 1, 
        img_bbox.getMaxY() - img_bbox.getMinY() + 1)
    diffexpCutout = diffexp.getCutout(center_point, size);
    
    # The actual cutout pixel rounding may be a little bit different
    bbox = diffexpCutout.getBBox()
    extentR = (bbox.getMaxY()+0.5, bbox.getMinY()-0.5, bbox.getMaxX()+0.5, bbox.getMinX()-0.5)
    diffexpArray = diffexpCutout.getMaskedImage().getImage().getArray()
    diffexpNorm = ImageNormalize(diffexpArray, interval=ZScaleInterval(), stretch=AsinhStretch())
    ax.imshow(diffexpArray.T[::-1,::-1], origin='lower', cmap=cmap, norm=diffexpNorm, extent = extentR)
    ax.set_xlim(plt_bbox.getMaxY()+1, plt_bbox.getMinY()-1)
    ax.set_ylim(plt_bbox.getMaxX()+1, plt_bbox.getMinX()-1)
    ax.grid(True)
    ax.scatter(Tsources['y'],Tsources['x'],s=8, alpha=0.3, c='red')

# ---
def plot_diasources_diffim(repo1, repo2,  srcTable1, srcTable2, 
        diffimType='deepDiff_differenceExp',pdfWriter=None, Nmaxpanels=None):
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
    if Nmaxpanels is not None:
        ccdVisitIds = ccdVisitIds[:Nmaxpanels]
    n_panels = len(ccdVisitIds)
    panel_idx = 1 # Overall subplot panel index

    butler1 = dafPersist.Butler(repo1)
    butler2 = dafPersist.Butler(repo2)
    while len(ccdVisitIds) > 0:
        if fig is None:
            print(fig_idx)
            fig = plt.figure(figsize=(9,10))            
            fig.subplots_adjust(left=0.05, right=0.98, top=0.94, 
                bottom=0.04, hspace=0.05, wspace=0.1)
            fig.suptitle('{} - {} of {}'.format(panel_idx, min(panel_idx+n_plots_per_fig - 1, n_panels),
                n_panels))
#            fig.suptitle('DIAObject {}; {}/{}'.format(obj,fig_idx,n_figs))
            fig_idx += 1  
            splot_idx = 1 # Subplot index within figure (first row only)
        
        ccdVisitId = ccdVisitIds.pop()

        # Upper row diffexp first repo
        
        ax = fig.add_subplot(2,n_plots_per_fig,splot_idx)
        T = srcTable1[srcTable1['ccdVisitId']==ccdVisitId]

        visit = ccdVisitId // 100
        ccd = ccdVisitId % 100
        dataId = { 'visit': visit, 'ccdnum': ccd }

        diffexp = butler1.get(diffimType, dataId)
        make_and_plot_one_cutout(ax, diffexp, T)
        T2 = srcTable2[srcTable2['ccdVisitId']==ccdVisitId]
        ax.set_title('{visit};{ccd:02d} $N_s$:{s1};{s2}'.format(visit=visit, ccd= ccd, s1=len(T), s2=len(T2)))
        
        T = T2
        ax = fig.add_subplot(2,n_plots_per_fig,splot_idx+n_plots_per_fig)
        diffexp = butler2.get(diffimType, dataId)
        make_and_plot_one_cutout(ax, diffexp, T,cmap='Blues_r')
        
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
#    cwpTemplateRepo = '/home/gkovacs/data/repo_DM-17825/templates'
    my_dbName = '/home/gkovacs/data/repo_DM-17825/ingested/rerun/proc_2019-02-21/association.db'
    mrawls_dbName = '/home/gkovacs/data/repo_DM-17825/mrawls_cw_processed2/association.db'
    mrawls_repo = '/home/gkovacs/data/repo_DM-17825/mrawls_cw_processed2'

    # Our processing
    conn = sqlite3.connect(my_dbName)
#    srcTable = Table.from_pandas(pd.read_sql_query('select * from  DiaSource '
#                                                   'where decl < -5.6 and decl > -5.8'
#                ' and ra > 155.2 and ra < 155.3 and flags = 0', conn))
    srcTable = Table.from_pandas(pd.read_sql_query('select * from  DiaSource '
                                                   'where decl < -5.6 and decl > -5.8'
                ' and ra > 155.2 and ra < 155.3', conn))

    # Meredith's ppdb

    m_conn = sqlite3.connect(mrawls_dbName)
#    m_srcTable = Table.from_pandas(pd.read_sql_query(
#        'select * from DiaSource where decl < -5.6 and decl > -5.8'
#                ' and ra > 155.2 and ra < 155.3 and flags = 0', m_conn))
    m_srcTable = Table.from_pandas(pd.read_sql_query(
        'select * from DiaSource where decl < -5.6 and decl > -5.8'
                ' and ra > 155.2 and ra < 155.3', m_conn))

    with PdfPages('compare_det_sources_all_2019-03-12.pdf') as W:
        plot_diasources_diffim(cwpRepo,mrawls_repo,srcTable,m_srcTable, pdfWriter=W)

if __name__ == "__main__":
    main()

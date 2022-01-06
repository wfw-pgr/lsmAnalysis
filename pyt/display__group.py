import sys, os, glob, re
import numpy                      as np
import nkUtilities.load__config   as lcf
import nkUtilities.cMapTri        as cmt
import nkUtilities.configSettings as cfs


# ========================================================= #
# ===  display__group                                  === #
# ========================================================= #

def display__group():

    # ------------------------------------------------- #
    # --- [1] File name Set                         --- #
    # ------------------------------------------------- #
    outFile  = "png/mshape_lsm.vtu"
    elemFile = "dat/elems.dat"
    nodeFile = "dat/nodes.dat"
    mshpFile = "dat/mshape_lsm.dat"
    
    # ------------------------------------------------- #
    # --- [2] Fetch Data                            --- #
    # ------------------------------------------------- #
    import nkUtilities.load__pointFile as lpf
    elems  = lpf.load__pointFile( inpFile=elemFile, returnType="point" )
    nodes  = lpf.load__pointFile( inpFile=nodeFile, returnType="point" )
    mshape = lpf.load__pointFile( inpFile=mshpFile, returnType="point" )
    elems  = np.array( elems, np.int64 )
    names  = [ "x", "y", "z", "init", "shim", "flag", "group" ]
    
    # ------------------------------------------------- #
    # --- [3] construct uGrid                       --- #
    # ------------------------------------------------- #
    import nkVTKRoutines.construct__uGrid as cug
    cug.construct__uGrid( nodes=nodes, elems=elems, vtkFile=outFile, \
                          cellData=mshape, cellDataName=names )
    return()


# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    display__group()

import os, sys
import numpy as np
import nkUtilities.load__config   as lcf
import nkUtilities.cMapTri        as cmt


# ========================================================= #
# ===  generate element info.                           === #
# ========================================================= #

def generate__element( mshFile="msh/mesh.msh" ):

    x_, y_, z_ = 0, 1, 2
    
    # ------------------------------------------------- #
    # --- [1] load element info.                    --- #
    # ------------------------------------------------- #
    import nkMeshRoutines.load__meshio as lms
    elems, nodes = lms.load__meshio( mshFile=mshFile, elementType="triangle", \
                                     returnType="node-elem" )

    # ------------------------------------------------- #
    # --- [2] load mshape info.                     --- #
    # ------------------------------------------------- #
    import nkUtilities.load__pointFile as lpf
    inpFile = "dat/mshape_svd.dat"
    grids   = lpf.load__pointFile( inpFile=inpFile, returnType="structured" )

    # ------------------------------------------------- #
    # --- [3] calculate element Center              --- #
    # ------------------------------------------------- #
    centers = calculate__elementCenter( elemType="triangle", nodes=nodes, elems=elems )

    # ------------------------------------------------- #
    # --- [4] interpolate on the center of element  --- #
    # ------------------------------------------------- #
    import nkInterpolator.interpolate__grid2point as g2p
    Data    = g2p.interpolate__grid2point( gridData=grids, pointData=centers, \
                                           method="linear" )

    # ------------------------------------------------- #
    # --- [5] add info. about the mshape            --- #
    # ------------------------------------------------- #
    init,flag = np.copy( Data[:,z_] )[:,None], np.ones( (Data.shape[0],1) )
    shim      = np.zeros( (Data.shape[0],1) )
    Data      = np.concatenate( [Data,init,shim,flag], axis=1 )
    names     = [ "xm", "ym", "zm", "mi", "ms", "mf" ]
    
    # ------------------------------------------------- #
    # --- [5] save in a file                        --- #
    # ------------------------------------------------- #
    import nkUtilities.save__pointFile as spf
    outFile   = "dat/mshape_lsm.dat"
    spf.save__pointFile( outFile=outFile, Data=Data, names=names )

    # ------------------------------------------------- #
    # --- [6] save as a png File                    --- #
    # ------------------------------------------------- #
    config                   = lcf.load__config()
    config["FigSize"]        = [6,6]
    config["cmp_position"]   = [0.16,0.16,0.90,0.86]
    config["xTitle"]         = "X (m)"
    config["yTitle"]         = "Y (m)"
    config["cmp_xAutoRange"] = True
    config["cmp_yAutoRange"] = True
    cmt.cMapTri( xAxis=Data[:,x_], yAxis=Data[:,y_], cMap=Data[:,z_], \
    	         pngFile="png/mshape_lsm.png", config=config )


    return()

# ========================================================= #
# ===  calculate__elementCenter                         === #
# ========================================================= #

def calculate__elementCenter( elemType="triangle", nodes=None, elems=None ):

    nd1_, nd2_, nd3_, nd4_ = 0, 1, 2, 3
    
    # ------------------------------------------------- #
    # --- [1] Arguments                             --- #
    # ------------------------------------------------- #
    if ( nodes is None ): sys.exit( "[calculate__elementCenter] nodes == ???" )
    if ( elems is None ): sys.exit( "[calculate__elementCenter] elems == ???" )
    if ( elemType.lower() in [ "triangle"] ):
        pass
    else:
        print( "[calculate__elementCenter] not supported elemType == {0}".format( elemType ) )
        sys.exit()
        
    # ------------------------------------------------- #
    # --- [2] calculate                             --- #
    # ------------------------------------------------- #
    node1 = nodes[ elems[:,nd1_], : ]
    node2 = nodes[ elems[:,nd2_], : ]
    node3 = nodes[ elems[:,nd3_], : ]
    ret   = ( node1 + node2 + node3 ) / 3.0

    # ------------------------------------------------- #
    # --- [3] return                                --- #
    # ------------------------------------------------- #
    return( ret )
    

# ========================================================= #
# ===   Execution of Pragram                            === #
# ========================================================= #

if ( __name__=="__main__" ):
    generate__element()

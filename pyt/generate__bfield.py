import os, sys
import numpy as np
import nkUtilities.load__config as lcf
import nkUtilities.cMapTri      as cmt


# ========================================================= #
# ===  generate__bfield                                 === #
# ========================================================= #

def generate__bfield():

    x_, y_, z_, bz_ = 0, 1, 2, 5
    
    # ------------------------------------------------- #
    # --- [1] load ideal bfields                    --- #
    # ------------------------------------------------- #
    import nkUtilities.load__pointFile as lpf
    inpFile = "dat/pole_ideal_mcoord.dat"
    biaData = lpf.load__pointFile( inpFile=inpFile, returnType="point" )
    
    # ------------------------------------------------- #
    # --- [2] load fem results                      --- #
    # ------------------------------------------------- #
    import nkUtilities.load__pointFile as lpf
    inpFile = "dat/ems_pst.field"
    emsData = lpf.load__pointFile( inpFile=inpFile, returnType="point" )

    # ------------------------------------------------- #
    # --- [3] merge fields                          --- #
    # ------------------------------------------------- #
    coord   = biaData[:,x_:z_+1]
    ideal   = biaData[:,bz_    ]
    backg   = emsData[:,bz_    ]
    Data    = np.concatenate( [coord,ideal[:,None],backg[:,None]], axis=1 )
    
    # ------------------------------------------------- #
    # --- [4] save in a file                        --- #
    # ------------------------------------------------- #
    import nkUtilities.save__pointFile as spf
    outFile   = "dat/bfield_input.dat"
    names     = [ "x_", "y_", "z_", "bi_", "bb_" ]
    spf.save__pointFile( outFile=outFile, Data=Data, names=names )

    # ------------------------------------------------- #
    # --- [5] sample png                            --- #
    # ------------------------------------------------- #
    coord_ = np.reshape( coord, (-1,3) )
    config                   = lcf.load__config()
    config["FigSize"]        = [6,6]
    config["cmp_position"]   = [0.16,0.16,0.90,0.86]
    config["xTitle"]         = "X (m)"
    config["yTitle"]         = "Y (m)"
    config["cmp_xAutoRange"] = True
    config["cmp_yAutoRange"] = True
    bi_, bb_                 = 3, 4
    cmt.cMapTri( xAxis=Data[:,x_], yAxis=Data[:,y_], cMap=Data[:,bb_], \
    	         pngFile="png/bfield_backg.png", config=config )
    cmt.cMapTri( xAxis=Data[:,x_], yAxis=Data[:,y_], cMap=Data[:,bi_], \
    	         pngFile="png/bfield_ideal.png", config=config )
    return()

# ========================================================= #
# ===   Execution of Pragram                            === #
# ========================================================= #
if ( __name__=="__main__" ):
    generate__bfield()

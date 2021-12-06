import sys, os, glob, re
import numpy                      as np
import nkUtilities.load__config   as lcf
import nkUtilities.cMapTri        as cmt
import nkUtilities.configSettings as cfs

# ========================================================= #
# ===  display__bfield                                  === #
# ========================================================= #
def display__bfield( num=None ):

    x_, y_, z_, i_, b_, s_, t_, e_ = 0, 1, 2, 3, 4, 5, 6, 7

    # ------------------------------------------------- #
    # --- [1] latest number check                   --- #
    # ------------------------------------------------- #
    if ( num is None ):
        files = glob.glob( "out/bfield_*.dat" )
        if   ( len( files ) == 0 ):
            print( "[display__bfield.py] no files :: out/bfield_*.dat [ERROR]" )
            sys.exit()
        else:
            nums    = []
            pattern = "out/bfield_([0-9]*).dat"
            matches = [ re.search( pattern, hfile ) for hfile in files ]
            numList = [ int( match.group(1) ) for match in matches if match ]
            num     = max( numList )
            
    # ------------------------------------------------- #
    # --- [2] File name Set                         --- #
    # ------------------------------------------------- #
    datFile = "out/bfield_{0:04}.dat"         .format( num )
    pngFile = "png/bfield_{0}" + "_{0:04}.png".format( num )

    # ------------------------------------------------- #
    # --- [3] Fetch Data                            --- #
    # ------------------------------------------------- #
    import nkUtilities.load__pointFile as lpf
    Data  = lpf.load__pointFile( inpFile=datFile, returnType="point" )
    
    # ------------------------------------------------- #
    # --- [4] config Settings                       --- #
    # ------------------------------------------------- #
    config  = lcf.load__config()
    cfs.configSettings( configType="cMap_def", config=config )
    config["FigSize"]        = (5,5)
    config["cmp_position"]   = [0.16,0.12,0.97,0.88]
    config["xTitle"]         = "X (m)"
    config["yTitle"]         = "Y (m)"
    config["cmp_xAutoRange"] = True
    config["cmp_yAutoRange"] = True
    config["cmp_xRange"]     = [-5.0,+5.0]
    config["cmp_yRange"]     = [-5.0,+5.0]

    # ------------------------------------------------- #
    # --- [5] plot Figure                           --- #
    # ------------------------------------------------- #
    cmt.cMapTri( xAxis=Data[:,x_], yAxis=Data[:,y_], cMap=Data[:,b_], \
    		 pngFile=pngFile.format( "backg" ), config=config )
    cmt.cMapTri( xAxis=Data[:,x_], yAxis=Data[:,y_], cMap=Data[:,i_], \
    		 pngFile=pngFile.format( "ideal" ), config=config )
    cmt.cMapTri( xAxis=Data[:,x_], yAxis=Data[:,y_], cMap=Data[:,e_], \
    		 pngFile=pngFile.format( "error" ), config=config )
    cmt.cMapTri( xAxis=Data[:,x_], yAxis=Data[:,y_], cMap=Data[:,s_], \
    		 pngFile=pngFile.format( "shimf" ), config=config )
    cmt.cMapTri( xAxis=Data[:,x_], yAxis=Data[:,y_], cMap=Data[:,t_], \
    		 pngFile=pngFile.format( "total" ), config=config )

# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    display__bfield()


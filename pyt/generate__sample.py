import numpy as np
import nkUtilities.load__config   as lcf
import nkUtilities.cMapTri        as cmt


# ========================================================= #
# ===  generate bfield sample                           === #
# ========================================================= #

def generate__sampleField():

    x_, y_, z_ = 0, 1, 2
    radius     = 0.9
    
    # ------------------------------------------------- #
    # --- [1]  generate bfield                      --- #
    # ------------------------------------------------- #
    import nkUtilities.equiSpaceGrid as esg
    x1MinMaxNum = [ -0.9, 0.9, 31 ]
    x2MinMaxNum = [ -0.9, 0.9, 31 ]
    x3MinMaxNum = [  0.0, 0.0,  1 ]
    coord       = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                     x3MinMaxNum=x3MinMaxNum, returnType = "point" )
    radii       = np.sqrt( coord[:,x_]**2 + coord[:,y_]**2 )
    index       = np.where( radii < radius )
    coord,radii = coord[index], radii[index]
    bb          = 0.050 * np.exp( ( -0.5 ) * ( radii / radius )**2 * 3.0 ) + 2.500
    db          = 0.050 * np.exp( ( -0.5 ) * ( radii / radius )**2 * 3.0 ) - 0.025
    bi          = bb * 0.0 + 2.5
    Data        = np.concatenate( [coord,bi[:,None],bb[:,None]], axis=1 )

    # ------------------------------------------------- #
    # --- [2] save in a file                        --- #
    # ------------------------------------------------- #
    import nkUtilities.save__pointFile as spf
    outFile   = "dat/bfield_input.dat"
    spf.save__pointFile( outFile=outFile, Data=Data )

    # ------------------------------------------------- #
    # --- [3] save in each file                     --- #
    # ------------------------------------------------- #
    zeros     = np.zeros_like( bb )
    Data1     = np.concatenate( [coord,zeros[:,None],zeros[:,None],bb[:,None]], axis=1 )
    Data2     = np.concatenate( [coord,zeros[:,None],zeros[:,None],bi[:,None]], axis=1 )
    outFile1  = "dat/ems_pst.field"
    outFile2  = "dat/pole_ideal_mcoord.dat"
    spf.save__pointFile( outFile=outFile1, Data=Data1 )
    spf.save__pointFile( outFile=outFile2, Data=Data2 )
    
    # ------------------------------------------------- #
    # --- [4] sample png                            --- #
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
# ===  generate sample shape                            === #
# ========================================================= #

def generate__sampleShape():

    x_, y_, z_  = 0, 1, 2
    zMin,zMax   = 0.020, 0.180
    radius      = 1.050
    
    # ------------------------------------------------- #
    # --- [1]  generate mshape                      --- #
    # ------------------------------------------------- #
    import nkUtilities.equiSpaceGrid as esg
    x1MinMaxNum = [ -1.1, 1.1, 201 ]
    x2MinMaxNum = [ -1.1, 1.1, 201 ]
    x3MinMaxNum = [  0.0, 0.0,   1 ]
    coord       = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                     x3MinMaxNum=x3MinMaxNum, returnType = "structured" )
    shape       = coord.shape[1:]
    coord       = np.reshape( coord, (-1,3) )
    radii       = np.sqrt( coord[:,x_]**2 + coord[:,y_]**2 )
    coord[:,z_] = ( zMax - zMin ) * np.exp( - ( radii/radius )**2 ) + zMin
    coord       = np.reshape( coord, shape )
    
    # ------------------------------------------------- #
    # --- [2] save in a file                        --- #
    # ------------------------------------------------- #
    import nkUtilities.save__pointFile as spf
    outFile   = "dat/mshape_svd.dat"
    spf.save__pointFile( outFile=outFile, Data=coord )

    # ------------------------------------------------- #
    # --- [3] sample png                            --- #
    # ------------------------------------------------- #
    coord_ = np.reshape( coord, (-1,3) )
    config                   = lcf.load__config()
    config["FigSize"]        = [6,6]
    config["cmp_position"]   = [0.16,0.16,0.90,0.86]
    config["xTitle"]         = "X (m)"
    config["yTitle"]         = "Y (m)"
    config["cmp_xAutoRange"] = True
    config["cmp_yAutoRange"] = True
    cmt.cMapTri( xAxis=coord_[:,x_], yAxis=coord_[:,y_], cMap=coord_[:,z_], \
    	         pngFile="png/mshape_svd.png", config=config )
    return()



# ========================================================= #
# ===  generate sample weight                           === #
# ========================================================= #

def generate__sampleWeight():

    x_, y_, z_  = 0, 1, 2
    
    # ------------------------------------------------- #
    # --- [1] load mshape file                      --- #
    # ------------------------------------------------- #
    import nkUtilities.load__pointFile as lpf
    inpFile = "dat/bfield_input.dat"
    coord   = lpf.load__pointFile( inpFile=inpFile, returnType="point" )
    coord   = lpf.load__pointFile( inpFile=inpFile, returnType="point" )

    # ------------------------------------------------- #
    # --- [2] concatenate color & weight info.      --- #
    # ------------------------------------------------- #
    color   = np.ones( ( coord.shape[0],1) )
    weight  = np.ones( ( coord.shape[0],1) )
    Data    = np.concatenate( [coord[:,0:3],color,weight], axis=1 )
    
    # ------------------------------------------------- #
    # --- [3] save in a file                        --- #
    # ------------------------------------------------- #
    import nkUtilities.save__pointFile as spf
    outFile   = "dat/weights.dat"
    spf.save__pointFile( outFile=outFile, Data=Data )
    return()

    

# ========================================================= #
# ===   Execution of Pragram                            === #
# ========================================================= #
if ( __name__=="__main__" ):
    generate__sampleField ()
    generate__sampleShape ()
    generate__sampleWeight()

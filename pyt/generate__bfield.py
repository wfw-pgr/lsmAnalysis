import numpy as np


# ========================================================= #
# ===  generate bfield sample                           === #
# ========================================================= #

def generate__bfield():

    x_, y_, z_ = 0, 1, 2
    
    # ------------------------------------------------- #
    # --- [1]  generate bfield                      --- #
    # ------------------------------------------------- #
    import nkUtilities.equiSpaceGrid as esg
    x1MinMaxNum = [ 0.0, 1.0, 21 ]
    x2MinMaxNum = [ 0.0, 1.0, 21 ]
    x3MinMaxNum = [ 0.0, 0.0,  1 ]
    coord       = esg.equiSpaceGrid( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                     x3MinMaxNum=x3MinMaxNum, returnType = "point" )
    radii       = np.sqrt( coord[:,x_]**2 + coord[:,y_]**2 )
    index       = np.where( radii < 1.0 )
    coord       = coord[index]
    coord[:,z_] = 2.5

    # ------------------------------------------------- #
    # --- [2] save in a file                        --- #
    # ------------------------------------------------- #
    import nkUtilities.save__pointFile as spf
    outFile   = "dat/bfield_input.dat"
    spf.save__pointFile( outFile=outFile, Data=coord )

    return()

# ========================================================= #
# ===   Execution of Pragram                            === #
# ========================================================= #
if ( __name__=="__main__" ):
    generate__bfield()

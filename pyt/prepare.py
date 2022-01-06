import os, sys

# ========================================================= #
# ===  prepare__scripts.py                              === #
# ========================================================= #

def prepare():

    # ------------------------------------------------- #
    # --- [1] file copying                          --- #
    # ------------------------------------------------- #
    import nkUtilities.load__constants as lcn
    cnsFile = "dat/unified.conf"
    const   = lcn.load__constants( inpFile=cnsFile )

    # ------------------------------------------------- #
    # --- [2] prepare parameter                     --- #
    # ------------------------------------------------- #
    import prepare__parameter as prm
    prm.prepare__parameter()

    # ------------------------------------------------- #
    # --- [3] generate mesh to be interpolated      --- #
    # ------------------------------------------------- #
    import generate__circleMesh as gcm
    if   ( const["lsm.pole.meshType"] == "edge" ):
        gcm.generate__circleMesh( side=const["general.side"], lc1=0.050, lc2=0.050 )
    elif ( const["lsm.pole.meshType"] == "direct-math" ):
        gcm.generate__circleMesh( side=const["general.side"], \
                                  directMath=const["lsm.pole.direct-math.mathEval"] )
    gcm.save__elem_and_node()

    # ------------------------------------------------- #
    # --- [4] interpolate height on mesh            --- #
    # ------------------------------------------------- #
    import generate__element as elm
    elm.generate__element()

    # ------------------------------------------------- #
    # --- [5] display grouping                      --- #
    # ------------------------------------------------- #
    import display__group as grp
    grp.display__group()

    return()


# ========================================================= #
# ===   Execution of Pragram                            === #
# ========================================================= #
if ( __name__=="__main__" ):
    prepare()

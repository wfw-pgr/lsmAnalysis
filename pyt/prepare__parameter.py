import os, sys
import nkUtilities.load__constants as lcn
import nkUtilities.save__constants as scn

# ========================================================= #
# ===  prepare__parameter.py                            === #
# ========================================================= #

def prepare__parameter( inpFile="dat/unified.conf", lstFile="dat/parameter.lst" ):

    contentsList = [ "lsm.iterMax" ]
    
    
    # ------------------------------------------------- #
    # --- [1] load config                           --- #
    # ------------------------------------------------- #
    const = lcn.load__constants( inpFile=inpFile )
    
    # ------------------------------------------------- #
    # --- [2] modify contents                       --- #
    # ------------------------------------------------- #

    # ------------------------------------------------- #
    # --- [3] save parameter.lst for Fortran        --- #
    # ------------------------------------------------- #
    #  -- [3-1] rename                              --  #
    namelist = { ( key.replace( "lsm.", "" ) ):const[key] for key in contentsList }
    keys     = list( namelist.keys() )
    
    #  -- [3-2] save in a file                      --  #
    group    = "parameters"
    import nkUtilities.save__namelist as snl
    snl.save__namelist( outFile=lstFile, const=namelist, keys=keys, group=group )

    
# ========================================================= #
# ===   実行部                                          === #
# ========================================================= #
if ( __name__=="__main__" ):
    prepare__parameter()





#     contentsList = [ "svd.iterMax", "svd.nDiv_B", "svd.threshold", "svd.resid_convergence", \
#                      "svd.wLaplace", "svd.wPicard_init", "svd.wPicard_last", \
#                      "svd.iPicard_iter", "svd.MzConstant", "svd.solverType", \
#                      "svd.LSM_Engine", "svd.modelType", "svd.rLim1", "svd.rLim2", \
#                      "svd.pLim1", "svd.pLim2", "svd.zLim1", "svd.zLim2", "svd.LI", \
#                      "svd.LJ", "svd.x1Min", "svd.x1Max", "svd.x2Min", "svd.x2Max", \
#                      "svd.except_sw", "svd.except_radii1", "svd.except_radii2", \
#                      "svd.except_theta1", "svd.except_theta2", "svd.restart_fromSave", \
#                      "svd.restart_filename" ]

#     listkeys     = [ "svd.iterMax", "svd.nDiv_B", "svd.threshold", "svd.wLaplace", \
#                      "svd.wPicard_init", "svd.wPicard_last", "svd.iPicard_iter", \
#                      "svd.MzConstant", "svd.solverType", "svd.LSM_Engine", "svd.modelType", \
#                      "svd.rLim1", "svd.rLim2", "svd.pLim1", "svd.pLim2", \
#                      "svd.zLim1", "svd.zLim2", "svd.resid_convergence" ]

    
#     # ------------------------------------------------- #
#     # --- [1] Arguments                             --- #
#     # ------------------------------------------------- #
#     if ( inpFile is None ): inpFile = "dat/unified.conf"
#     if ( outFile is None ): outFile = "dat/parameter.conf"

#     # ------------------------------------------------- #
#     # --- [2] load config                           --- #
#     # ------------------------------------------------- #
#     const = lcn.load__constants( inpFile=inpFile )
    
#     # ------------------------------------------------- #
#     # --- [3] obtain constants                      --- #
#     # ------------------------------------------------- #
#     #  -- [3-1] set r/z Limit value                 --  #
#     const["svd.rLim1"] = 0.0
#     const["svd.rLim2"] = const["geometry.r_pole"]
#     const["svd.zLim1"] = const["geometry.z_gapMin"]
#     const["svd.zLim2"] = const["geometry.z_gapMax"]

#     #  -- [3-2] check side & set phi Limit          --  #
#     if   ( const["general.side"] == "+" ):
#         const["svd.pLim1"]     = -90.0
#         const["svd.pLim2"]     = +90.0
#         const["svd.modelType"] = "half"
#         if ( const["svd.x1Min"] != 0.0 ):
#             print( "[prepare__parameter.py] svd.x1Min != 0.0, although general.side= + "  )
#             sys.exit()
            
#     elif ( const["general.side"] == "-" ):
#         const["svd.pLim1"]     =  +90.0
#         const["svd.pLim2"]     = +270.0
#         const["svd.modelType"] = "full"
#         if ( const["svd.x1Max"] != 0.0 ):
#             print( "[prepare__parameter.py] svd.x1Max != 0.0, although general.side= - "  )
#             sys.exit()

#     elif ( const["general.side"] in ["+-","-+"] ):
#         const["svd.pLim1"]     = -180.0
#         const["svd.pLim2"]     = +180.0
#         const["svd.modelType"] = "full"
#     else:
#         print( "[prepare__parameter.py] general.side ?? :: {0}".format( const["general.side"] ) )
#         sys.exit()

#     # ------------------------------------------------- #
#     # --- [4] save in parameter.conf                --- #
#     # ------------------------------------------------- #
#     ret = scn.save__constants( const=const, keys=contentsList, outFile=outFile )
#     print( "[prepare__parameter.py] outFile :: {0}".format( outFile ) )

#     # ------------------------------------------------- #
#     # --- [5] save parameter.lst for Fortran        --- #
#     # ------------------------------------------------- #
    
#     #  -- [5-1] rename                              --  #
#     namelist = { ( key.replace( "svd.", "" ) ):const[key] for key in listkeys }
#     keys     = list( namelist.keys() )
    
#     #  -- [5-2] save in a file                      --  #
#     lstFile  = "dat/parameter.lst"
#     group    = "parameters"
#     import nkUtilities.save__namelist as snl
#     snl.save__namelist( outFile=lstFile, const=namelist, keys=keys, group=group )
    


    
# # ========================================================= #
# # ===   実行部                                          === #
# # ========================================================= #

# if ( __name__=="__main__" ):
#     inpFile = "dat/unified.conf"
#     outFile = "dat/parameter.conf"
#     prepare__parameter( inpFile=inpFile, outFile=outFile )






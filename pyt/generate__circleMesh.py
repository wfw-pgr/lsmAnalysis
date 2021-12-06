import numpy   as np
import os, sys, subprocess
import gmsh
# import gmsh_api.gmsh as gmsh


# ========================================================= #
# === save elements / nodes                             === #
# ========================================================= #

def save__elem_and_node( mshFile="msh/mesh.msh" ):
    
    # ------------------------------------------------- #
    # --- [1] load mesh file                        --- #
    # ------------------------------------------------- #
    import nkMeshRoutines.load__meshio as lms
    elems, nodes = lms.load__meshio( mshFile=mshFile, elementType="triangle", \
                                     returnType="node-elem" )

    # ------------------------------------------------- #
    # --- [2] save in a file                        --- #
    # ------------------------------------------------- #
    import nkUtilities.save__pointFile as spf
    outFile   = "dat/elems.dat"
    spf.save__pointFile( outFile=outFile, Data=elems, fmt="%12d" )
    print( np.min( elems ), np.max( elems ) )

    import nkUtilities.save__pointFile as spf
    outFile   = "dat/nodes.dat"
    spf.save__pointFile( outFile=outFile, Data=nodes, fmt="%15.8e" )
    
    return()
    

# ========================================================= #
# ===   generate__circleMesh.py                         === #
# ========================================================= #

def generate__circleMesh( lc1=0.0, lc2=0.0, radius=1.0, side="+" ):

    # ------------------------------------------------- #
    # --- [1] initialization of the gmsh            --- #
    # ------------------------------------------------- #
    gmsh.initialize()
    gmsh.option.setNumber( "General.Terminal", 1 )
    gmsh.model.add( "model" )

    # ------------------------------------------------- #
    # --- [2] Modeling                              --- #
    # ------------------------------------------------- #
    r1, r2  = 0.0, radius
    origin  = [ 0.0, 0.0 ]
    import nkGmshRoutines.generate__sector90 as s90
    if ( side in ["+","+-","-+"] ):
        ret1    = s90.generate__sector90( lc=lc1, r1=r1, r2=r2, quadrant=1, \
                                          origin=origin, defineSurf=True )
        ret2    = s90.generate__sector90( lc=lc2, r1=r1, r2=r2, quadrant=4, \
                                          origin=origin, defineSurf=True )

    if ( side in ["-","+-","-+"] ):
        ret3    = s90.generate__sector90( lc=lc1, r1=r1, r2=r2, quadrant=2, \
                                          origin=origin, defineSurf=True )
        ret4    = s90.generate__sector90( lc=lc2, r1=r1, r2=r2, quadrant=3, \
                                          origin=origin, defineSurf=True )
    gmsh.model.occ.synchronize()
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()
    
    # ------------------------------------------------- #
    # --- [3] Mesh generation & save                --- #
    # ------------------------------------------------- #
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.write( "msh/mesh.msh" )
    gmsh.finalize()


# ========================================================= #
# ===   Execution of Pragram                            === #
# ========================================================= #

if ( __name__=="__main__" ):
    side = "+-"
    generate__circleMesh( side=side )
    save__elem_and_node()

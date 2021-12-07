module allocatorMod
contains
  
  ! ====================================================== !
  ! === allocate__variables                            === !
  ! ====================================================== !
  subroutine allocate__variables
    use variablesMod
    use lstStructMod
    implicit none

    ! ------------------------------------------------------ !
    ! --- [1] preparation                                --- !
    ! ------------------------------------------------------ !
    write(6,"(a)",advance="no") "[allocate__variables] allocating variables..... "
    nMpt = nBpt
    nNpt = nGroups
    nEpt = nElems

    ! ------------------------------------------------------ !
    ! --- [2] allocation                                 --- !
    ! ------------------------------------------------------ !
    allocate( Rmat(nMpt,nEpt), Amat(nMpt,nNpt) )
    allocate( rhs(nMpt), hvec(nEpt), xvec(nNpt), wvec(nMpt) )
    allocate( vertex(dim,nVert,nElems) )

    ! ------------------------------------------------------ !
    ! --- [3] initialization                             --- !
    ! ------------------------------------------------------ !
    Rmat(:,:)     = 0.d0
    Amat(:,:)     = 0.d0
    rhs(:)        = 0.d0
    hvec(:)       = 0.d0
    wvec(:)       = 0.d0
    vertex(:,:,:) = 0.d0
    
    write(6,"(3x,a)") "[Done]"
    return
  end subroutine allocate__variables
  
  
  ! ====================================================== !
  ! === deallocate__variables                          === !
  ! ====================================================== !
  subroutine deallocate__variables
    use variablesMod
    implicit none
    
    ! ------------------------------------------------------ !
    ! --- [1] deallocation                               --- !
    ! ------------------------------------------------------ !
    write(6,"(a)",advance="no") "[deallcate__variables] deallocating variables..... "
    deallocate( Rmat, Amat )
    deallocate( rhs , hvec )
    deallocate( vertex )
    write(6,"(3x,a)") "[Done]"
    return
  end subroutine deallocate__variables
  
  
end module allocatorMod

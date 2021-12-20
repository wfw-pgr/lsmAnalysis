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
    nNpt = nUsed
    nEpt = nElems
    if ( flag__smoothing ) then
       nMpt = nMpt + nUsed
    endif

    ! ------------------------------------------------------ !
    ! --- [2] allocation                                 --- !
    ! ------------------------------------------------------ !
    allocate( Rmat(nBpt,nEpt)         , source=0.d0 )
    allocate( Amat(nMpt,nNpt)         , source=0.d0 )
    allocate( Mmat(nEpt,nEpt)         , source=0.d0 )
    allocate( Lmat(nNpt,nNpt)         , source=0.d0 )
    allocate( rhs (nMpt)              , source=0.d0 )
    allocate( hvec(nEpt)              , source=0.d0 )
    allocate( xvec(nNpt)              , source=0.d0 )
    allocate( wvec(nMpt)              , source=0.d0 )
    allocate( kvec(nEpt,1)            , source=0.d0 )
    allocate( lvec(nUsed)             , source=0.d0 )
    allocate( vertex(dim,nVert,nElems), source=0.d0 )
    
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

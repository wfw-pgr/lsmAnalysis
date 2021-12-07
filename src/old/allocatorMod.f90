module allocatorMod
contains

  ! ===================================================== !
  ! ===  allocation of variables                      === !
  ! ===================================================== !
  subroutine allocate__variables
    use variablesMod
    implicit none
    
    ! ------------------------------------- !
    ! --- [1] size determination        --- !
    ! ------------------------------------- !
    if ( Flag__laplacian ) then
       nMpt = nBpt + nNpt
    else
       nMpt = nBpt
    endif
    ! ------------------------------------- !
    ! --- [2] Matrix                    --- !
    ! ------------------------------------- !
    allocate( Rmat(nBpt,nNpt), Amat(nMpt,nNpt), Lmat(nNpt,nNpt) )
    allocate( rhs(nMpt), shim(nNpt), spectrum(nMpt), singular(nMpt) )
    allocate( lapl(nNpt), surf(nNpt) )
    Amat(:,:)   = 0.d0
    Rmat(:,:)   = 0.d0
    Lmat(:,:)   = 0.d0
    rhs (:)     = 0.d0
    shim(:)     = 0.d0
    spectrum(:) = 0.d0
    singular(:) = 0.d0
    lapl(:)     = 0.d0
    surf(:)     = 0.d0
    
    return
  end subroutine allocate__variables


  ! ===================================================== !
  ! ===  deallocation of variables                      === !
  ! ===================================================== !
  subroutine deallocate__variables
    use variablesMod
    implicit none
    
    deallocate( Rmat, Amat )
    deallocate( rhs, shim )
    deallocate( idtable )
    
    return
  end subroutine deallocate__variables
  
end module allocatorMod

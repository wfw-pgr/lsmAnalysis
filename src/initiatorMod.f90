module initiatorMod
contains

  ! ====================================================== !
  ! === initialize data to be used in the analysis     === !
  ! ====================================================== !
  subroutine initialize__conditions
    use variablesMod
    implicit none
    integer :: ie, iv, irow
    
    ! ------------------------------------------------------ !
    ! --- [1] prepare magnetization vector               --- !
    ! ------------------------------------------------------ !
    mvec(1) = 0.0
    mvec(2) = 0.0
    mvec(3) = MzConst

    ! ------------------------------------------------------ !
    ! --- [2] pre-define iter for save__results          --- !
    ! ------------------------------------------------------ !
    iter    = 0

    ! ------------------------------------------------------ !
    ! --- [3] make vertex for faster access              --- !
    ! ------------------------------------------------------ !
    do ie=1, nElems
       ! -- prepare -- !
       do iv=1, nVert
          vertex(:,iv,ie) = nodes( :, elems(iv,ie) )
       enddo
    enddo

    ! ------------------------------------------------------ !
    ! --- [4] extract weight info from weights           --- !
    ! ------------------------------------------------------ !
    nColor = 0
    do irow=1, nMpt
       wvec(irow) = weights(wt_,irow)
       nColor     = max( nColor, int( weights(cl_,irow) ) )
    enddo
    if ( flag__laplaceRegularization ) then
       nColor = nColor - 1
    endif
    
    return
  end subroutine initialize__conditions

  
end module initiatorMod

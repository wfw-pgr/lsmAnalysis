module invSolverMod
contains


  ! ====================================================== !
  ! === Weighted Least Squared Method Solver           === !
  ! ====================================================== !
  subroutine WLS_MatrixInvertor( Amat, bvec, xvec, wvec, LM, LN, lsm__engine )
    use cglsdenseMod, only : cgls4dense
    use utilitiesMod, only : DisplayClock
    implicit none
    integer                       :: i, j, LdA, LdW, LP, Lwork, info, LdB, nRHS
    double precision              :: winv, rwork(1)
    integer         , intent(in)  :: LM, LN
    double precision, intent(in)  :: Amat(LM,LN), bvec(LM), wvec(LM)
    character(4)    , intent(in)  :: lsm__engine
    double precision, intent(out) :: xvec(LN)
    double precision, allocatable :: yvec(:), Umat(:,:), work(:)
    
    ! ------------------------------------- !
    ! --- [1] Preparation               --- !
    ! ------------------------------------- !
    !  -- [1-1] Nortification           --  !
    write(6,"(a)"          ) "[WLS_MatrixInvertor] weighted Least Squared Method start now..."
    write(6,"(a,i10,a,i10)") "[WLS_MatrixInvertor]    M = ", LM, "  , N = ", LN
    write(6,"(a)"          ) "[WLS_MatrixInvertor] solving Ax=b ..."
    call DisplayClock
    
    !  -- [1-2] Prepare variables       --  !
    LdA   = LM
    LdB   = LM
    nRHS  = 1
    allocate( Umat(LM,LN), yvec(LM) )

    !  -- [1-3] Substitution            --  !
    do i=1, LM
       do j=1, LN
          Umat(i,j) = Amat(i,j) * wvec(i)
       enddo
       yvec(i)      = bvec(i)   * wvec(i)
    enddo

    ! ------------------------------------------------------ !
    ! --- [2] Least Squared Method                       --- !
    ! ------------------------------------------------------ !

    !  -- [2-1] CGLS case                                --  !
    if ( ( lsm__engine.eq."CGLS" ).or.( lsm__engine.eq."cgls" ) ) then
       call cgls4dense( Umat, yvec, xvec, LM, LN )
    endif
    
    !  -- [2-2] Lapack version                           --  !
    if ( ( lsm__engine.eq."QRLS" ).or.( lsm__engine.eq."qrls" ) ) then
       ! -- get optimal size of work -- !
       Lwork = -1
       info  =  0
       call dgels ( 'N', LM, LN, nRHS, Umat, LdA, yvec, LdB, rwork, Lwork, info )
       Lwork = rwork(1)
       allocate( work( Lwork ) )
       work(:) = 0.d0
    
       !  -- call DGELS LSM solver   -- !
       info  =  0
       call dgels ( 'N', LM, LN, nRHS, Umat, LdA, yvec, LdB,  work, Lwork, info )
    
       !  -- erorr detection         -- !
       if ( info.ne.0 ) write(6,'(a)') '[WLS_MatrixInvertor]  info != 0  Error !!!'

       !  -- substitute answer       -- !
       do i=1, LN
          xvec(i) = yvec(i)
       enddo
       deallocate( work )
    end if

    ! ------------------------------------- !
    ! --- [3] Post Process              --- !
    ! ------------------------------------- !
    call DisplayClock
    write(6,"(a)") "[WLS_MatrixInvertor] solving Ax=b ...       [ End ]"
    deallocate( yvec, Umat )
    
    return
  end subroutine WLS_MatrixInvertor

end module invSolverMod

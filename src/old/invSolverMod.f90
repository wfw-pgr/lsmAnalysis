module invSolverMod
contains
  
  ! ===================================================== !
  ! ===  LSM_MatrixInvertor  ( solve [A]{x}={b} )     === !
  ! ===================================================== !
  subroutine LSM_MatrixInvertor( Amat, bvec, xvec, LM, LN, LSM_engine )
    use cglsdenseMod, only : cgls4dense
    use utilitiesMod, only : DisplayClock
    implicit none
    integer                       :: i, j, LdA, nRHS, LdB, Lwork, info
    integer         , intent(in)  :: LM, LN
    double precision, intent(in)  :: Amat(LM,LN), bvec(LM)
    character(4)    , intent(in)  :: LSM_engine
    double precision, intent(out) :: xvec(LN)
    double precision, allocatable :: yvec(:,:), work(:), Umat(:,:)
    double precision              :: rwork(1)

    ! ------------------------------------- !
    ! --- [1] Preparation               --- !
    ! ------------------------------------- !
    !  -- [1-1] Nortification           --  !
    write(6,'(a)') '[LSM_MatrixInvertor] Least Squared Method ( solving Ax=b ) start now...'
    write(6,'(a,i10,a,i10)') '[LSM_MatrixInvertor]  M = ', LM, '  , N = ', LN
    call DisplayClock
    !  -- [1-2] Prepare size of array   --  !
    nRHS  = 1
    LdA   = LM
    LdB   = max( LM, LN )
    allocate( yvec(LdB,nRHS), Umat(LM,LN) )
    yvec  = 0.d0
    !  -- [1-3] Substitution            --  !
    do j=1, LN
       do i=1, LM
          Umat(i,j) = Amat(i,j)
       enddo
    enddo
    do i=1, LM
       yvec(i,1) = bvec(i)
    enddo

    ! ------------------------------------------------------ !
    ! --- [2] Least Squared Method                       --- !
    ! ------------------------------------------------------ !

    !  -- [2-1] CGLS case                                --  !
    if ( trim(LSM_engine).eq.'CGLS' ) then
       call cgls4dense( Umat, yvec, xvec, LM, LN )
       
    endif
    
    !  -- [2-2] Lapack version                           --  !
    if ( trim(LSM_engine).eq.'QRLS' ) then
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
       if ( info.ne.0 ) write(6,'(a)') '[LSM_MatrixInvertor]  info != 0  Error !!!'
       deallocate( work )

       !  -- substitute answer       -- !
       do i=1, LN
          xvec(i) = yvec(i,1)
       enddo
    end if
       
    ! ------------------------------------- !
    ! --- [3] Post Process              --- !
    ! ------------------------------------- !
    call DisplayClock
    write(6,'(a)') '[LSM_MatrixInvertor] Least Squared Method ( solving Ax=b ) End.'
    deallocate( yvec, Umat )
    
    return
  end subroutine LSM_MatrixInvertor


  ! ===================================================== !
  ! ===  WLS_MatrixInvertor  ( solve [A]{x}={b} )     === !
  ! ===================================================== !
  subroutine WLS_MatrixInvertor( Amat, bvec, xvec, wvec, LM, LN, LSM_engine )
    use cglsdenseMod, only : cgls4dense
    use utilitiesMod, only : DisplayClock
    implicit none
    integer                       :: i, j, LdA, LdW, LP, Lwork, info, LdB, nRHS
    double precision              :: winv, rwork(1)
    integer         , intent(in)  :: LM, LN
    double precision, intent(in)  :: Amat(LM,LN), bvec(LM), wvec(LM)
    character(4)    , intent(in)  :: LSM_engine
    double precision, intent(out) :: xvec(LN)
    double precision, allocatable :: yvec(:), Umat(:,:), work(:)
    
    ! ------------------------------------- !
    ! --- [1] Preparation               --- !
    ! ------------------------------------- !
    !  -- [1-1] Nortification           --  !
    write(6,'(a)') '[WLS_MatrixInvertor] Least Squared Method ( solving Ax=b ) start now...'
    write(6,'(a,i10,a,i10)') '[WLS_MatrixInvertor]  M = ', LM, '  , N = ', LN
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
    if ( LSM_engine.eq.'CGLS' ) then
       call cgls4dense( Umat, yvec, xvec, LM, LN )
    endif
    
    !  -- [2-2] Lapack version                           --  !
    if ( LSM_engine.eq.'QRLS' ) then
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
    write(6,'(a)') '[WLS_MatrixInvertor] Least Squared Method ( solving Ax=b ) End.'
    deallocate( yvec, Umat )
    
    return
  end subroutine WLS_MatrixInvertor
  
  
  
  ! ===================================================== !
  ! ===  SVD Matrix Invertor  ( solve [A]{x}={b} )    === !
  ! ===================================================== !
  subroutine SVD_MatrixInvertor( Amat, bvec, xvec, spectrum, singular, LM, LN, threshold )
    use utilitiesMod, only         : DisplayClock, save__matrix
    implicit none
    integer                       :: imode, i, j, LdA, LdU, LdVT, Lwork, info, minMN
    integer         , intent(in)  :: LM, LN
    double precision, intent(in)  :: threshold
    double precision, intent(in)  :: Amat(LM,LN),  bvec(LM)
    double precision, intent(out) :: xvec(LN), spectrum(LM), singular(LM)
    integer         , allocatable :: iwork(:)
    double precision, allocatable :: Wmat(:,:), Umat(:,:), VmatT(:,:)
    double precision, allocatable :: SmatD(:), pvec(:), work(:)
    double precision, parameter   :: alpha = 1.d0, beta  = 0.d0
    integer         , parameter   :: incX  = 1, incY  = 1, inspect_mode = -1
    logical         , parameter   :: Flag__debugMode    = .false.
    logical         , parameter   :: Flag__fasterSolver = .true.
    
    ! ------------------------------------- !
    ! --- [1] Preparation               --- !
    ! ------------------------------------- !
    !  -- [1-1] Nortification           --  !
    write(6,'(a)') '[SVD_MatrixInvertor] Singular Value Decomposition ( solving Ax=b ) start now...'
    write(6,'(a,i10,a,i10)') '[SVD_MatrixInvertor]  M = ', LM, '  , N = ', LN
    call DisplayClock
    !  -- [1-2] size of arrays          --  !
    LdA        = LM
    LdU        = LM
    LdVT       = LN
    if ( LM.lt.LN ) then
       write(6,'(a)') '[SVD_MatrixInvertor] CAUTION !!! ( LM < LN )'
       minMN      = LM
       ! stop
    else
       minMN      = LN
    endif
    !  -- [1-3] allocate arrays         --  !
    allocate( Wmat(LM,LN), Umat(LM,LM), VmatT(LN,LN), SmatD(LN), pvec(LN)  )
    allocate( iwork(8*minMN) )
    Wmat (:,:) = Amat(:,:)
    Umat (:,:) = 0.d0
    VmatT(:,:) = 0.d0
    SmatD(:)   = 0.d0
    pvec (:)   = 0.d0
    iwork(:)   = 0
    
    ! ------------------------------------- !
    ! --- [2] Lapack Execution          --- !
    ! ------------------------------------- !
    !  -- [2-1] inspect optimal Lwork   --  !
    allocate  ( work(    1) )
    work(:) = 0.d0
    if ( Flag__fasterSolver ) then
       ! -- dgesdd :: faster solver than dgesvd  -- !
       call dgesdd( 'A'     , LM, LN, Wmat, LdA, SmatD, Umat, LdU, VmatT, LdVT, work, inspect_mode, iwork, info )
    else
       ! -- dgesvd :: more formal svd algorithm  -- !
       call dgesvd( 'A', 'A', LM, LN, Wmat, LdA, SmatD, Umat, LdU, VmatT, LdVT, work, inspect_mode,        info )
    endif
    !  -- [2-2] allocation of Lwork     --  !
    Lwork   = work(1)
    deallocate( work        )
    allocate  ( work(Lwork) )
    work(:) = 0.d0
    !  -- [2-3] call DGESDD SVD solver  --  !
    if ( Flag__fasterSolver ) then
       ! -- dgesdd :: faster solver than dgesvd  -- !
       call dgesdd( 'A'     , LM, LN, Wmat, LdA, SmatD, Umat, LdU, VmatT, LdVT, work, Lwork, iwork, info )
    else
       ! -- dgesvd :: more formal svd algorithm  -- !
       call dgesvd( 'A', 'A', LM, LN, Wmat, LdA, SmatD, Umat, LdU, VmatT, LdVT, work, Lwork,        info )
    endif
    !  -- [2-4] erorr detection         --  !
    if ( info.ne.0 ) write(6,'(a)') '[SVD_MatrixInvertor]  info != 0  Error !!!'
    if ( Flag__debugMode ) then
       call save__matrix( SmatD, LM,  1, 'chk/Smat.dat' )
       call save__matrix( Umat , LM, LM, 'chk/Umat_new.dat' )
       call save__matrix( VmatT, LN, LN, 'chk/Vmat_new.dat' )
    endif
    
    ! ------------------------------------- !
    ! --- [3] Least Square Solution     --- !
    ! ------------------------------------- !
    !  -- LSM   :: {x} = [V][S+][U*]{b}     !
    !  --      for [A] = [U] [S+] [V*]      !
    !  -----------------------------------  !
    !  -- [3-1] {y} = [U*] {b}          --  !
    spectrum(:) = 0.d0
    call dgemv( 'T', LM, LM, alpha, Umat, LdU, bvec, incX, beta, spectrum, incY )
    !  -- [3-2] {p} = [S+] [U*] {b}     --  !
    !   -  Trancation / division        -   !
    do j=1, minMN
       if ( abs( SmatD(j) ).gt.threshold ) then
          pvec(j)     = spectrum(j) / SmatD(j)
          singular(j) = SmatD(j)
       else
          pvec(j)     = 0.d0
          singular(j) = SmatD(j)
       endif
    enddo
    !  -- [3-3] {x} = [V][S+][U*] {b}   --  !
    call dgemv( 'T', LN, LN, alpha, VmatT, LdVT, pvec, incX, beta, xvec, incY )

    ! ------------------------------------- !
    ! --- [4] post process              --- !
    ! ------------------------------------- !
    deallocate( Wmat, Umat, VmatT, SmatD, pvec  )
    deallocate( iwork, work )
    call DisplayClock
    write(6,'(a)') '[SVD_MatrixInvertor] Singular Value Decomposition ( solving Ax=b ) End.'

    return
  end subroutine SVD_MatrixInvertor


  ! ===================================================== !
  ! ===  GLS_MatrixInvertor  ( solve [A]{x}={b} )     === !
  ! ===================================================== !
  subroutine GLS_MatrixInvertor( Amat, bvec, xvec, wvec, LM, LN )
    use utilitiesMod, only : DisplayClock
    implicit none
    integer                       :: i, j, LdA, LdW, LP, Lwork, info
    double precision              :: rwork(1)
    integer         , intent(in)  :: LM, LN
    double precision, intent(in)  :: Amat(LM,LN), bvec(LM), wvec(LM)
    double precision, intent(out) :: xvec(LN)
    double precision, allocatable :: yvec(:), zvec(:), Umat(:,:), Wmat(:,:), work(:), rhs(:)
    
    ! ------------------------------------- !
    ! --- [1] Preparation               --- !
    ! ------------------------------------- !
    !  -- [1-1] Nortification           --  !
    write(6,'(a)') '[GLS_MatrixInvertor] Least Squared Method ( solving Ax=b ) start now...'
    write(6,'(a,i10,a,i10)') '[GLS_MatrixInvertor]  M = ', LM, '  , N = ', LN
    call DisplayClock
    
    !  -- [1-2] Prepare variables       --  !
    LP    = LM
    LdA   = LM
    LdW   = LM
    allocate( Umat(LM,LN), Wmat(LM,LP), rhs(LM), zvec(LN), yvec(LP) )
    
    !  -- [1-3] Prepare matrix          --  !
    !   - Umat  -   !
    do j=1, LN
       do i=1, LM
          Umat(i,j) = Amat(i,j)
       enddo
    enddo
    !   - Wmat  -   !
    do j=1, LP
       do i=1, LM
          if ( i.eq.j ) then
             Wmat(i,j) = 1.d0 / wvec(i)
          else
             Wmat(i,j) = 0.d0
          endif
       enddo
    enddo
    !   - rhs   -   !
    do i=1, LM
       rhs(i)  = bvec(i)
    enddo
    !   - zvec  -   !
    do i=1, LN
       zvec(i) = 0.d0
    enddo
    !   - yvec  -   !
    do i=1, LP
       yvec(i) = 0.d0
    enddo
    
    ! ------------------------------------------------------ !
    ! --- [2] Lapack Execution                           --- !
    ! ------------------------------------------------------ !
    !  -- [2-1] get optimal work size                    --  !
    Lwork = -1
    info  =  0
    call DGGGLM ( LM, LN, LP, Umat, LdA, Wmat, LdW, rhs, zvec, yvec, rwork, Lwork, info )
    Lwork = rwork(1)
    allocate( work( Lwork ) )
    work(:) = 0.d0
    
    !  -- [2-2] solve weighted Least Square              --  !
    info  =  0
    call DGGGLM ( LM, LN, LP, Umat, LdA, Wmat, LdW, rhs, zvec, yvec, work , Lwork, info )
    
    !  -- [2-3] erorr detection         --  !
    if ( info.ne.0 ) write(6,'(a)') '[GLS_MatrixInvertor]  info != 0  Error !!!'

    ! ------------------------------------- !
    ! --- [3] Post Process              --- !
    ! ------------------------------------- !
    call DisplayClock
    write(6,'(a)') '[GLS_MatrixInvertor] Least Squared Method ( solving Ax=b ) End.'
    do i=1, LN
       xvec(i) = zvec(i)
    enddo
    deallocate( yvec, zvec, work, Umat, Wmat, rhs )
    
    return
  end subroutine GLS_MatrixInvertor

  
end module invSolverMod

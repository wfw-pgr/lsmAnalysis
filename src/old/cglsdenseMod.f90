module cglsdenseMod
contains

  
  ! ====================================================== !
  ! === cgls4dense                                     === !
  ! ====================================================== !
  subroutine cgls4dense( Amat, bvec, xvec, LI, LJ )
    implicit none
    integer         , intent(in)    :: LI, LJ
    double precision, intent(inout) :: Amat(LI,LJ), bvec(LI), xvec(LJ)
    integer                         :: iter
    double precision                :: alpha, beta, rNorm, r0Norm, bNorm
    double precision, allocatable   :: Adp(:), ATAdp(:)
    double precision, allocatable   :: pvec(:), rvec(:), r0vec(:)
    integer                         :: iterMax
    integer         , parameter     :: incX = 1, incY = 1
    integer         , parameter     :: iterMax_Factor = 20
    double precision, parameter     :: convergence    = 1.d-12
    
    ! -- BLAS Reference                                                -- !
    ! call DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
    ! call DAXPY ( n, alpha, x, incx, y, incy )
    ! see... http://www.netlib.org/blas/                               -- !
    ! -- algorithm  ( Bi-CGStab )                                      -- !
    ! see... http://www.jicfus.jp/wiki/index.php?Bi-CGSTAB%20%E6%B3%95 -- !
    ! ------------------------------------------------------------------- !
    
    ! ------------------------------------------------------ !
    ! --- [1] preparation                                --- !
    ! ------------------------------------------------------ !
    iterMax = int( max( LI, LJ ) * iterMax_Factor )
    allocate( Adp (LI), ATAdp(LJ) )
    allocate( pvec(LJ), rvec(LJ) , r0vec(LJ) )
    
    ! ------------------------------------------------------ !
    ! --- [2] initial settings                           --- !
    ! ------------------------------------------------------ !
    !  -- (1) x = 0                 -- !
    xvec(:)  = 0.d0
    rvec(:)  = 0.d0
    !  -- (2) r0 = A^T ( b - Ax0 )  -- !
    Adp(:)   = bvec(:)
    call dgemv(  "N", LI, LJ, -1.d0, Amat, LI, rvec, incX, +1.d0, Adp , incY )
    call dgemv(  "T", LI, LJ, +1.d0, Amat, LI, Adp , incX, +0.d0, rvec, incY )
    !  -- (4) p=r, r0=r             -- !
    pvec(:)  = rvec(:)
    r0vec(:) = rvec(:)
    !  -- (6) rNorm, r0Norm, bNorm  -- !
    bNorm    = inner_product( bvec, bvec, LI )
    rNorm    = inner_product( rvec, rvec, LJ )
    r0Norm   = rNorm

    ! ------------------------------------------------------ !
    ! --- [3] Main Loop                                  --- !
    ! ------------------------------------------------------ !
    !  - rvec :: rk
    !  - svec :: sk, r_k+1
    !  - pvec :: pk
    !  ----------------------------------------------------- !
    
    do iter=1, iterMax

       if ( rNorm.lt.convergence*bNorm ) then
          write(6,*)
          write(6,*) "[cglsdense.f90] CGLS method reached Convergence..."
          write(6,*) "[cglsdense.f90]   residuals             :: ", rNorm
          write(6,*) "[cglsdense.f90]   convergence           :: ", convergence
          write(6,*) "[cglsdense.f90]   criterion             :: ", convergence*bNorm
          write(6,*)
          exit
       endif
       
       ! -- (6) calculate A.p                         -- !
       Adp      = 0.d0
       call dgemv( "N", LI, LJ, +1.d0, Amat, LI, pvec, incX, +0.d0, Adp  , incY )

       ! -- (7) calculation A^T. A.p -- !
       call dgemv( "T", LI, LJ, +1.d0, Amat, LI, Adp , incX, +0.d0, ATAdp, incY )
       
       ! -- (8) alpha = ( r0, r ) / ( r0, A.p )       -- !
       alpha    = rNorm / inner_product( pvec, ATAdp , LJ )

       ! -- (9)  update xvec -- !
       call daxpy( LJ,      alpha,  pvec, incX, xvec, incY )
       ! xvec(:)  = xvec(:) + alpha * pvec(:)
       
       ! -- (10) update rvec -- !
       call daxpy( LJ, -1.0*alpha, ATAdp, incX, rvec, incY )
       ! rvec(:)  = rvec(:) - alpha * ATAdp(:)

       ! -- (11) beta = ( ri+1, ri+1 ) / ( ri, ri )    -- !
       r0Norm   = rNorm
       rNorm    = inner_product( rvec, rvec, LJ )
       beta     = rNorm / r0Norm
       
       ! -- (12) update pvec -- !
       pvec(:)  = rvec(:) + beta * pvec(:)
       
    enddo
    
    return
  end subroutine cgls4dense

  
  ! ====================================================== !
  ! === Function :: inner product :: x.y               === !
  ! ====================================================== !
  Function inner_product( xvec, yvec, Npt )
    implicit none
    integer         , intent(in) :: Npt
    double precision, intent(in) :: xvec(Npt), yvec(Npt)
    integer                      :: i
    double precision             :: inner_product
    
    inner_product = 0.d0
    !$omp parallel default(none) &
    !$omp shared(Npt,inner_product,xvec,yvec) private(i)
    !$omp do reduction(+:inner_product)
    do i=1, Npt
       inner_product = inner_product + xvec(i) * yvec(i)
    enddo
    !$omp end do
    !$omp end parallel
    
    return
  end Function inner_product
  
end module cglsdenseMod

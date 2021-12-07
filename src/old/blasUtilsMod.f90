module blasUtilsMod
contains

  ! ===================================================== !
  ! ===  BLAS Matrix Vector Multiply Wrapper          === !
  ! ===================================================== !
  subroutine MatrixVectorMultiply( Amat, xvec, yvec, LM, LN )
    implicit none
    integer         , intent(in)  :: LM, LN
    double precision, intent(in)  :: Amat(LM,LN), xvec(LN)
    double precision, intent(out) :: yvec(LM)
    integer         , parameter   :: incX=1, incY=1
    double precision, parameter   :: c1=1.d0, c2=0.d0
    
    ! ------------------------------------- !
    ! --- call DGEMV                    --- !
    ! ------------------------------------- !
    yvec(:)   = 0.d0
    call DGEMV( 'N', LM, LN, c1, Amat, LM, xvec, incX, c2, yvec, incY )

    return
  end subroutine MatrixVectorMultiply

  
end module blasUtilsMod

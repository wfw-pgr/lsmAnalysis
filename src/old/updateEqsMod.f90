module updateEqsMod
contains
  
  ! ===================================================== !
  ! ===  update Matrix                                === !
  ! ===================================================== !
  subroutine update__Matrix
    use variablesMod
    implicit none
    
    ! ------------------------------------- !
    ! --- [1] ResponseMatrix            --- !
    ! ------------------------------------- !
    call generate__responseMatrix
    
    ! ------------------------------------- !
    ! --- [2] Laplacian Matrix          --- !
    ! ------------------------------------- !
    if ( flag__Laplacian ) then
       call generate__laplaceMatrix
    endif
    
    return
  end subroutine update__Matrix

  
  ! ===================================================== !
  ! ===  update R.H.S.                                === !
  ! ===================================================== !
  subroutine update__RHS
    use variablesMod
    implicit none
    integer              :: i, j, iB, iN
    integer, parameter   :: z_=3, f_=6, e_=8
    
    ! ------------------------------------- !
    ! --- [1] update error Field        --- !
    ! ------------------------------------- !
    !$omp parallel default(none) shared(rhs,BField,nBpt) private(iB)
    !$omp do
    do iB=1, nBpt
       rhs(iB) = BField(e_,iB)
    enddo
    !$omp end do
    !$omp end parallel
    ! ------------------------------------- !
    ! --- [2] update delta Laplace      --- !
    ! ------------------------------------- !
    if ( flag__Laplacian ) then
       lapl(:)          = 0.d0
       surf(:)          = 0.d0
       !$omp parallel default(none) &
       !$omp shared(LI,LJ,mshape,surf,lapl,idtable) private(i,j)
       !$omp do
       do j=1, LJ
          do i=1, LI
             if ( mshape(f_,i,j).eq.1.d0 ) then
                surf(idtable(i,j)) = mshape(z_,i,j)
                lapl(idtable(i,j)) = 0.d0
             endif
          enddo
       enddo
       !$omp end do
       !$omp end parallel
       call MatrixVectorMultiply( Lmat, surf, lapl, nNpt, nNpt )
       !$omp parallel default(none) shared(rhs,lapl,nNpt,nBpt) private(iN)
       !$omp do
       do iN=1, nNpt
          ! rhs(nBpt+iN) = - lapl(iN)
          rhs(nBpt+iN) = 0.d0
       enddo
       !$omp end do
       !$omp end parallel
    endif
    
    return
  end subroutine update__RHS
  

  ! ===================================================== !
  ! ===  generate responseMatrix                      === !
  ! ===================================================== !
  subroutine generate__responseMatrix
    use variablesMod
    use utilitiesMod
    implicit none
    integer                     :: col, row, i, j, ii, jj
    double precision            :: xm, ym, zm, xm0, ym0, xb, yb, zb, vol, radius, theta
    double precision            :: shmpos(3), bfdpos(3)
    integer         , parameter :: x_=1, y_=2, z_=3, i_=4, f_=6


    ! ------------------------------------- !
    ! --- [1] preparation               --- !
    ! ------------------------------------- !
    Rmat(:,:) = 0.d0
    col       = 0
    zb        = 0.d0
    vol       = ddx * ddy * unit_thick
    ! ------------------------------------- !
    ! --- [2] main Loop                 --- !
    ! ------------------------------------- !
    ! -- for (i,j) [ each shim piece ] loop -- !
    if ( Flag__fullModel ) then
       
       do j=1, LJ
          ym0 = yMin + dble(j-1)*dy - dy*0.5d0 + ddy*0.5d0
          do i=1, LI
             xm0 = xMin + dble(i-1)*dx - dx*0.5d0 + ddx*0.5d0
             zm  = mshape(z_,i,j)

             if ( mshape(f_,i,j).eq.1.d0 ) then
                col = idtable(i,j)
                ! -- for each BField point loop               -- !
                do row = 1, nBpt
                   bfdpos(x_:z_) = BField(x_:z_,row)

                   ! -- for (ii,jj) [ small shim piece ] loop -- !
                   do jj=1, nDiv_B
                      ym = ym0 + (jj-1)*ddy
                      do ii=1, nDiv_B
                         xm     = xm0 + (ii-1)*ddx
                         radius = sqrt( xm**2+ ym**2 )
                         theta  = atan( ym, xm )
                         if ( theta.lt.0.d0 ) theta = theta + pJump

                         if ( ( radius.ge.rLim1 ).and.( radius.le.rLim2 ).and.( theta.ge.pLim1 ).and.( theta.le.pLim2 ) ) then
                            Rmat(row, col) = Rmat(row, col) &
                                 & + vol*( + MagneticMoment(  xm,  ym,  zm, bfdpos, mvec ) &
                                 &         + MagneticMoment(  xm,  ym, -zm, bfdpos, mvec ) )
                         endif

                      end do
                   end do
                   ! -- for (ii,jj) [ small shim piece ] loop -- !
                end do
                ! -- for each BField point loop               -- !

             endif
          enddo
       enddo

    else

       do j=1, LJ
          ym0 = yMin + dble(j-1)*dy - dy*0.5d0 + ddy*0.5d0
          do i=1, LI
             xm0 = xMin + dble(i-1)*dx - dx*0.5d0 + ddx*0.5d0
             zm  = mshape(z_,i,j)

             if ( mshape(f_,i,j).eq.1.d0 ) then
                col = idtable(i,j)
                ! -- for each BField point loop               -- !
                do row = 1, nBpt
                   bfdpos(x_:z_) = BField(x_:z_,row)

                   ! -- for (ii,jj) [ small shim piece ] loop -- !
                   do jj=1, nDiv_B
                      ym = ym0 + (jj-1)*ddy
                      do ii=1, nDiv_B
                         xm     = xm0 + (ii-1)*ddx
                         radius = sqrt( xm**2+ ym**2 )
                         theta  = atan( ym, xm )
                         if ( theta.lt.0.d0 ) theta = theta + pJump

                         if ( ( radius.ge.rLim1 ).and.( radius.le.rLim2 ).and.( theta.ge.pLim1 ).and.( theta.le.pLim2 ) ) then
                            Rmat(row, col) = Rmat(row, col) &
                                 & + vol*( + MagneticMoment(  xm,  ym,  zm, bfdpos, mvec ) &
                                 &         + MagneticMoment(  xm,  ym, -zm, bfdpos, mvec ) &
                                 &         + MagneticMoment( -xm,  ym,  zm, bfdpos, mvec ) &
                                 &         + MagneticMoment( -xm,  ym, -zm, bfdpos, mvec ) )
                         endif

                      end do
                   end do
                   ! -- for (ii,jj) [ small shim piece ] loop -- !
                end do
                ! -- for each BField point loop               -- !

             endif
          enddo
       enddo
       
    end if
    ! -- for (i,j) [ each shim piece ] loop -- !
    ! ------------------------------------- !
    ! --- [3] store in Amat             --- !
    ! ------------------------------------- !
    !$omp parallel default(none) shared(Amat,Rmat,nNpt,nBpt) private(col,row)
    !$omp do
    do col=1, nNpt
       do row=1, nBpt
          Amat(row,col) = Rmat(row,col)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    return
  end subroutine generate__responseMatrix

  
  ! ===================================================== !
  ! ===  generate__laplaceMatrix                      === !
  ! ===================================================== !
  subroutine generate__laplaceMatrix
    use variablesMod
    implicit none
    integer                       :: i, j, col, row, count
    integer         , parameter   :: f_ = 6
    double precision              :: dx2Inv, dy2Inv
    logical                       :: Flag__edge
    double precision              :: xm0, ym0, xm1, xm2, ym1, ym2, rm1, rm2, rm3, rm4
    
    ! ------------------------------------- !
    ! --- [1] Laplace matrix make       --- !
    ! ------------------------------------- !
    dx2Inv    = 1.0 / dx**2
    dy2Inv    = 1.0 / dy**2
    Lmat(:,:) = 0.d0

    if ( Flag__fullModel ) then

       !$omp parallel default(none) &
       !$omp shared(LI,LJ,mshape,idtable,Lmat,dx2Inv,dy2Inv) private(i,j,Flag__edge)
       !$omp do
       do j=1, LJ
          do i=1, LI
             if ( mshape(f_,i,j).eq.1.d0 ) then

                Flag__edge = .false.

                if (   ( idtable(i-1,j).gt.0 ).and.( idtable(i,j-1).gt.0 ).and. &
                     & ( idtable(i+1,j).gt.0 ).and.( idtable(i,j+1).gt.0 ) ) then
                   Lmat(idtable(i,j),idtable(i-1,j)) = + dx2Inv
                   Lmat(idtable(i,j),idtable(i+1,j)) = + dx2Inv
                   Lmat(idtable(i,j),idtable(i  ,j)) = - 2.d0*dx2Inv - 2.d0*dy2Inv
                   Lmat(idtable(i,j),idtable(i,j-1)) = + dy2Inv
                   Lmat(idtable(i,j),idtable(i,j+1)) = + dy2Inv
                else
                   Lmat(idtable(i,j),idtable(i,j)) = + 2.d0*( dx2Inv + dy2Inv )
                endif
                
             endif
          enddo
       enddo
       !$omp end do
       !$omp end parallel
       
    else

       !$omp parallel default(none) &
       !$omp shared(LI,LJ,mshape,idtable,Lmat,dx2Inv,dy2Inv) private(i,j,Flag__edge)
       !$omp do
       do j=1, LJ
          do i=1, LI
             if ( mshape(f_,i,j).eq.1.d0 ) then

                Flag__edge = .false.

                if (   ( idtable(i-1,j).gt.0 ).and.( idtable(i,j-1).gt.0 ).and. &
                     & ( idtable(i+1,j).gt.0 ).and.( idtable(i,j+1).gt.0 ) ) then
                   Lmat(idtable(i,j),idtable(i-1,j)) = + dx2Inv
                   Lmat(idtable(i,j),idtable(i+1,j)) = + dx2Inv
                   Lmat(idtable(i,j),idtable(i  ,j)) = - 2.d0*dx2Inv - 2.d0*dy2Inv
                   Lmat(idtable(i,j),idtable(i,j-1)) = + dy2Inv
                   Lmat(idtable(i,j),idtable(i,j+1)) = + dy2Inv
                else if ( ( i.eq.1  ).and.( idtable(i,j-1).gt.0 ).and. &
                     &    ( idtable(i+1,j).gt.0 ).and.( idtable(i,j+1).gt.0 ) ) then
                   Lmat(idtable(i,j),idtable(i+1,j)) = + 1.d0*dx2Inv
                   Lmat(idtable(i,j),idtable(i  ,j)) = - 1.d0*dx2Inv - 2.d0*dy2Inv
                   Lmat(idtable(i,j),idtable(i,j-1)) = + dy2Inv
                   Lmat(idtable(i,j),idtable(i,j+1)) = + dy2Inv
                else if ( ( i.eq.LI ).and.( idtable(i,j-1).gt.0 ).and. &
                     &    ( idtable(i-1,j).gt.0 ).and.( idtable(i,j+1).gt.0 ) ) then
                   Lmat(idtable(i,j),idtable(i-1,j)) = + 1.d0*dx2Inv
                   Lmat(idtable(i,j),idtable(i  ,j)) = - 1.d0*dx2Inv - 2.d0*dy2Inv
                   Lmat(idtable(i,j),idtable(i,j-1)) = + dy2Inv
                   Lmat(idtable(i,j),idtable(i,j+1)) = + dy2Inv
                else
                   Lmat(idtable(i,j),idtable(i,j)) = + 2.d0*( dx2Inv + dy2Inv )
                endif

             endif
          enddo
       enddo
       !$omp end do
       !$omp end parallel

    endif

    ! ------------------------------------- !
    ! --- [3] store Lmat                --- !
    ! ------------------------------------- !
    !$omp parallel default(none) shared(nNpt,nBpt,wLaplace,Amat,Lmat) private(row,col)
    !$omp do
    do col=1, nNpt
       do row=1, nNpt
          Lmat(     row,col) = wLaplace * 0.25d0 * Lmat(row,col)
          Amat(nBpt+row,col) =                     Lmat(row,col)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    return
  end subroutine generate__laplaceMatrix
  
  
  ! ===================================================== !
  ! ===  calculate Magnet Field                       === !
  ! ===================================================== !
  subroutine update__BField
    !$ use omp_lib
    use variablesMod
    implicit none
    integer                       :: i, j, k, ii, jj, kk, row, ith, nOMP
    double precision              :: xm, ym, zm, xm0, ym0, vol, ddz, radius, theta
    double precision              :: bfdpos(3)
    double precision, allocatable :: BzPara(:,:)
    integer         , parameter   :: x_=1, y_=2, z_=3, i_=4, b_=5, s_=6, t_=7, e_=8
    
    ! ------------------------------------- !
    ! --- [1] Preparation               --- !
    ! ------------------------------------- !
    !$omp parallel default(none) shared(nOMP)
    !$ nOMP = omp_get_num_threads()
    !$omp end parallel
    allocate( BzPara(nBpt,nOMP) )
    BzPara(:,:) = 0.d0

    ! ------------------------------------- !
    ! --- [2] OpenMP calculation        --- !
    ! ------------------------------------- !

    if ( flag__fullModel ) then
       
       !$omp parallel default(none) &
       !$omp shared(BzPara,BField,mshape,LI,LJ,nDiv_B,xMin,yMin,dx,dy,ddx,ddy,nOMP,nBpt,mvec) &
       !$omp shared(theta,rLim1,rLim2,pLim1,pLim2,pJump) &
       !$omp private(i,j,ii,jj,kk,row,xm0,ym0,xm,ym,zm,radius,ddz,vol,ith,bfdpos)
       !$ ith = omp_get_thread_num() + 1
       !$omp do
       do j=1, LJ
          ym0 = yMin + dble(j-1)*dy - dy*0.5d0 + ddy*0.5d0
          do i=1, LI
             xm0    = xMin + dble(i-1)*dx - dx*0.5d0 + ddx*0.5d0
             !  -- for all (i,j) point of the shim piece  --  !

             ddz    = ( mshape(i_,i,j) - mshape(z_,i,j) ) / nDiV_B
             vol    = ddx * ddy * ddz
             !  -- devide shim into 1 / nDiv_B^3  piece   --  !
             do kk = 1, nDiv_B
                zm  = mshape(z_,i,j) + dble(kk-1)*ddz + 0.5d0*ddz
                do jj = 1, nDiv_B
                   ym = ym0 + dble(jj-1)*ddy
                   do ii = 1, nDiv_B

                      xm     = xm0 + dble(ii-1)*ddx                   
                      radius = sqrt( xm**2 + ym**2 )
                      theta  = atan( ym, xm )
                      if ( theta.lt.0.d0 ) theta = theta + pJump

                      if (   ( radius.ge.rLim1 ).and.( radius.le.rLim2 ).and.&
                           & ( theta .ge.pLim1 ).and.( theta .le.pLim2 ) ) then

                         do row=1, nBpt
                            bfdpos(x_:z_)    = BField(x_:z_,row)
                            BzPara(row,ith)  = BzPara(row,ith) &
                                 &           + vol * (   MagneticMoment(  xm,  ym,  zm, bfdpos, mvec ) &
                                 &                     + MagneticMoment(  xm,  ym, -zm, bfdpos, mvec ) )
                         enddo
                      endif

                   enddo
                enddo
             enddo
             !  -- nDiv_B loop ( end )                     --  !
          end do
       end do
       !$omp end do
       !$omp end parallel

    else

       !$omp parallel default(none) &
       !$omp shared(BzPara,BField,mshape,LI,LJ,nDiv_B,xMin,yMin,dx,dy,ddx,ddy,nOMP,nBpt,mvec) &
       !$omp shared(theta,rLim1,rLim2,pLim1,pLim2,pJump) &
       !$omp private(i,j,ii,jj,kk,row,xm0,ym0,xm,ym,zm,radius,ddz,vol,ith,bfdpos)
       !$ ith = omp_get_thread_num() + 1
       !$omp do
       do j=1, LJ
          ym0 = yMin + dble(j-1)*dy - dy*0.5d0 + ddy*0.5d0
          do i=1, LI
             xm0    = xMin + dble(i-1)*dx - dx*0.5d0 + ddx*0.5d0
             !  -- for all (i,j) point of the shim piece  --  !

             ddz    = ( mshape(i_,i,j) - mshape(z_,i,j) ) / nDiV_B
             vol    = ddx * ddy * ddz
             !  -- devide shim into 1 / nDiv_B^3  piece   --  !
             do kk = 1, nDiv_B
                zm  = mshape(z_,i,j) + dble(kk-1)*ddz + 0.5d0*ddz
                do jj = 1, nDiv_B
                   ym = ym0 + dble(jj-1)*ddy
                   do ii = 1, nDiv_B

                      xm     = xm0 + dble(ii-1)*ddx                   
                      radius = sqrt( xm**2 + ym**2 )
                      theta  = atan( ym, xm )
                      if ( theta.lt.0.d0 ) theta = theta + pJump

                      if (   ( radius.ge.rLim1 ).and.( radius.le.rLim2 ).and.&
                           & ( theta .ge.pLim1 ).and.( theta .le.pLim2 ) ) then

                         do row=1, nBpt
                            bfdpos(x_:z_)    = BField(x_:z_,row)
                            BzPara(row,ith)  = BzPara(row,ith) &
                                 &           + vol * (   MagneticMoment(  xm,  ym,  zm, bfdpos, mvec ) &
                                 &                     + MagneticMoment(  xm,  ym, -zm, bfdpos, mvec ) &
                                 &                     + MagneticMoment( -xm,  ym,  zm, bfdpos, mvec ) &
                                 &                     + MagneticMoment( -xm,  ym, -zm, bfdpos, mvec ) )
                         enddo
                      endif

                   enddo
                enddo
             enddo
             !  -- nDiv_B loop ( end )                     --  !
          end do
       end do
       !$omp end do
       !$omp end parallel
       
    endif

    ! ------------------------------------- !
    ! --- [3] Manual Reduction          --- !
    ! ------------------------------------- !
    BField(s_,:) = 0.d0
    do ith=1, nOMP
       do row=1, nBpt
          BField(s_,row) = BField(s_,row) + BzPara(row,ith)
       enddo
    enddo

    ! ------------------------------------- !
    ! --- [4] update BField             --- !
    ! ------------------------------------- !
    !$omp parallel default(none) shared(BField,nBpt) private(row)
    !$omp do
    do row=1, nBpt
       BField(t_,row) = BField(b_,row) + BField(s_,row)
       BField(e_,row) = BField(i_,row) - BField(t_,row)
    enddo
    !$omp end do
    !$omp end parallel
    
    ! ------------------------------------- !
    ! --- [5] PostProcess               --- !
    ! ------------------------------------- !
    deallocate( BzPara )
    return
  end subroutine update__BField

  
  ! ====================================================== !
  ! === calculate MagnetiMoment from position          === !
  ! ====================================================== !
  function MagneticMoment( xm, ym, zm, bfdpos, mvec )
    implicit none
    double precision, intent(in) :: xm, ym, zm, bfdpos(3), mvec(3)
    double precision             :: rabs, mdotr, rvec3Inv, rvec(3), rhat(3)
    double precision             :: MagneticMoment
    double precision, parameter  :: fourpi = 16.d0*atan(1.d0)
    integer         , parameter  :: x_=1, y_=2, z_=3

    rvec(x_)          = bfdpos(x_) - xm
    rvec(y_)          = bfdpos(y_) - ym
    rvec(z_)          = bfdpos(z_) - zm
    rabs              = sqrt( rvec(1)**2 + rvec(2)**2 + rvec(3)**2 )
    rhat(:)           = rvec(:) / rabs
    rvec3Inv          = 1.d0 / ( fourpi * rabs**3 )
    mdotr             = mvec(1)*rhat(1) + mvec(2)*rhat(2) + mvec(3)*rhat(3)
    MagneticMoment    = rvec3Inv * ( 3.d0 * mdotr * rhat(z_) - mvec(z_) )
    return
  end function MagneticMoment


  ! ===================================================== !
  ! ===  update Magnet Shape                          === !
  ! ===================================================== !
  subroutine update__MagnetShape
    use variablesMod
    implicit none
    integer                     :: i, j, iterP
    double precision            :: coef
    integer         , parameter :: z_=3, i_=4, s_=5, f_=6

    ! ------------------------------------- !
    ! --- [1] Add iron(shim) on Magnet  --- !
    ! ------------------------------------- !
    !  -- [1-1] Coefficient             --  !
    iterP   = min( iter-1, iPicard_iter )
    wPicard = + wPicard_init  &
         &    + ( wPicard_last - wPicard_init ) / dble(iPicard_iter) * dble(iterP)
    coef    = - 1.d0 * wPicard
    !  -- [1-2] pileup/dig Magnet       --  !
    !$omp parallel default(none) &
    !$omp shared(LI,LJ,mshape,shim,idtable,coef,zLim1,zLim2) private(i,j)
    !$omp do
    do j=1, LJ
       do i=1, LI
          if ( mshape(f_,i,j).eq.1.d0 ) then
             mshape(s_,i,j) = shim( idtable(i,j) ) * unit_thick
             mshape(z_,i,j) = mshape(z_,i,j) + coef*mshape(s_,i,j)
             mshape(z_,i,j) = max( min( mshape(z_,i,j), zLim2 ), zLim1 )
          endif
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    return
  end subroutine update__MagnetShape


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

  
end module updateEqsMod

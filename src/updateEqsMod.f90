module updateEqsMod
contains
  
  ! ====================================================== !
  ! === update matrix                                  === !
  ! ====================================================== !
  subroutine update__matrix
    implicit none
    
    ! ------------------------------------------------------ !
    ! --- [1] Response Matrix                            --- !
    ! ------------------------------------------------------ !
    call generate__responseMatrix
    
    ! ------------------------------------------------------ !
    ! --- [2] Laplacian Matrix                           --- !
    ! ------------------------------------------------------ !
    ! if ( flag__Laplacian ) then
    !    call generate__laplaceMatrix
    ! endif
    return
  end subroutine update__matrix

  
  ! ====================================================== !
  ! === generate response matrix                       === !
  ! ====================================================== !
  subroutine generate__responseMatrix
    use variablesMod
    use bTriPrismMod
    implicit none
    integer                     :: ib, ie, iv, col, row
    double precision            :: bval, bsum, bfdpos(dim), zp(2), vert(dim,nVert)
    
    ! ------------------------------------------------------ !
    ! --- [1] calculate response matrix                  --- !
    ! ------------------------------------------------------ !
    do ib=1, nBpt
       bfdpos(:)  = bfield(xb_:zb_,ib)
       
       do ie=1, nElems
          zp(lo_)     = mshape(zm_,ie) - 0.5d0*unit__thickness
          zp(hi_)     = mshape(zm_,ie) + 0.5d0*unit__thickness
          
          vert(:,:)   = vertex(:,:,ie)
          call magneticField__triPrism( vert, bfdpos, mvec, zp, bval, nSubdiv, nDiv_z )
          
          Rmat(ib,ie) = bval
       enddo
    enddo
    ! -- single-side ver. -- !
    ! zp(lo_)    = mshape(zm_,ie)
    ! zp(hi_)    = mshape(zm_,ie) + unit__thickness

    ! ------------------------------------------------------ !
    ! --- [2] Store in Amatrix                           --- !
    ! ------------------------------------------------------ !
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
  
  
  ! ====================================================== !
  ! === update rhs                                     === !
  ! ====================================================== !
  subroutine update__rhs
    use variablesMod
    implicit none
    integer :: iB
    
    ! ------------------------------------- !
    ! --- [1] update error Field        --- !
    ! ------------------------------------- !
    !$omp parallel default(none) shared(rhs,bfield,nBpt) private(iB)
    !$omp do
    do iB=1, nBpt
       rhs(iB) = bfield(be_,iB)
    enddo
    !$omp end do
    !$omp end parallel
    return
  end subroutine update__rhs

  
  ! ====================================================== !
  ! === update bfield                                  === !
  ! ====================================================== !
  subroutine update__bfield
    use variablesMod
    use bTriPrismMod
    implicit none
    integer          :: ib, ie, iv
    double precision :: bval, bsum, bfdpos(dim), zp(2), vert(dim,nVert)

    ! ------------------------------------------------------ !
    ! --- [1] summing up bfield from all prisms          --- !
    ! ------------------------------------------------------ !
    do ib=1, nBpt
       bsum      = 0.d0
       bfdpos(:) = bfield(xb_:zb_,ib)
       
       do ie=1, nElems
          zp(lo_)    = mshape(mi_,ie)
          zp(hi_)    = mshape(zm_,ie)
          vert(:,:)  = vertex(:,:,ie)
          call magneticField__triPrism( vert, bfdpos, mvec, zp, bval, nSubdiv, nDiv_z )
          
          bsum       = bsum + bval
       enddo
       bfield(bs_,ib) = bsum
    enddo

    ! ------------------------------------------------------ !
    ! --- [2] update Bfield                              --- !
    ! ------------------------------------------------------ !
    do ib=1, nBpt
       bfield(bt_,ib) = bfield(bb_,ib) + bfield(bs_,ib)
       bfield(be_,ib) = bfield(bt_,ib) - bfield(bi_,ib)
    enddo
    
    return
  end subroutine update__bfield


  ! ====================================================== !
  ! === update Picard coefficient                      === !
  ! ====================================================== !
  subroutine update__Picard
    use variablesMod
    implicit none
    integer          :: ik
    double precision :: dPicard

    ! ------------------------------------------------------ !
    ! --- [1] prepare Error Detection                    --- !
    ! ------------------------------------------------------ !
    coefPicard = -1.d0
    if ( iPicard(1).ne.1 ) then
       write(6,*) "[update__mshape] iPicard is NOT 1 .... [ERROR] "
       stop
    endif
    ! ------------------------------------------------------ !
    ! --- [2] linear interpolation from parameter        --- !
    ! ------------------------------------------------------ !
    do ik=nMaxPicard, 1, -1
       if ( ( iPicard(ik).gt.0 ).and.( iter.ge.iPicard(ik) ) ) then
          ! -- [2-1] determine dPicard                   --  !
          if      ( ik.eq.nMaxPicard ) then
             dPicard = 0.d0
          else if ( iPicard(ik+1).lt.0 ) then
             dPicard = 0.d0
          else
             dPicard = ( wPicard(ik+1) - wPicard(ik) ) / dble( iPicard(ik+1) - iPicard(ik) )
          endif
          ! -- [2-2] linear interpolation                --  !
          coefPicard = wPicard(ik) + dPicard * dble( iter - iPicard(ik) )
          exit
       endif

    enddo
    ! ------------------------------------------------------ !
    ! --- [3] Error Detection                            --- !
    ! ------------------------------------------------------ !
    if ( coefPicard.le.0.d0 ) then
       write(6,*) "[update__Picard] coefPicard is Negative.... [ERROR]", coefPicard
       write(6,*) "[update__Picard] see wPicard & iPicard in parameter.lst...."
       stop
    endif
    write(6,'(a,f10.5)') '[update__Picard] coefPicard        = ', coefPicard

    return
  end subroutine update__Picard
  

  ! ====================================================== !
  ! === update mshape                                  === !
  ! ====================================================== !
  subroutine update__mshape
    use variablesMod
    implicit none
    integer :: ie

    ! ------------------------------------------------------ !
    ! --- [1] pilling-up / digging Magnet pole surface   --- !
    ! ------------------------------------------------------ !
    do ie=1, nElems
       mshape(ms_,ie) = hvec(ie) * unit__thickness
       mshape(zm_,ie) = mshape(zm_,ie) - coefPicard*mshape(ms_,ie)
       mshape(zm_,ie) = max( min( mshape(zm_,ie), zLim2 ), zLim1 )
    enddo

    return
  end subroutine update__mshape

  
end module updateEqsMod




 ! ! ====================================================== !
 !  ! === calculate magnetic field from prism            === !
 !  ! ====================================================== !
 !  subroutine magneticField__triPrism( bfdpos, mvec )
 !    implicit none
 !    integer, parameter            :: nVert = 3
 !    integer                       :: nSubdiv = 2
 !    double precision              :: lnode(3,nVert)
 !    double precision, allocatable :: snode(:,:,:)
 !    double precision, intent(in)  :: mvec(3), bfdpos(3)

 !    allocate( snode(3,nVert,nSubE), CoGs(3,nSubE) )
    
 !    do iE=1, nElems

 !       ! -- fetch node info. from mother -- !
 !       do iv=1, nVert
 !          lnode(:,iv) = nodes( :, elems(iv,iE) )
 !       enddo
       
 !       ! -- subdivide element            -- !
 !       call subdivide__element( lnode, snode, nVert, nSubdiv )
       
 !       ! -- calculate center of gravity  -- !
 !       call calculate__centerOfGravity( snode, CoGs, nSubE, nVert )

 !       ! -- calculate magnetic moment    -- !
 !       ret = 0.d0
 !       do is=1, nSubE
 !          xm   = CoGs(x_,is)
 !          ym   = CoGs(y_,is)
 !          delz = ( zh - zl ) / dble( nDiv_z - 1 )
 !          do iz=1, nDiv_z
 !             zm  =  zl + delz * dble( iz-1 )
 !             ret = ret + MagneticMoment( xm, ym, zm, bfdpos, mvec )
 !          enddo
 !       enddo
 !       write(6,*) ret
 !    end do
    
 !    return
 !  end subroutine magneticField__triPrism


 !  ! ====================================================== !
 !  ! === calculate center of gravity of a triangle      === !
 !  ! ====================================================== !
 !  subroutine calculate__centerOfGravity( snode, CoGs, nSubE, nVert )
 !    implicit none
 !    integer         , intent(in)  :: nSubE, nVert
 !    double precision, intent(in)  :: snode(3,nVert,nSubE)
 !    double precision, intent(out) :: CoGs (3,nSubE)
    
 !    CoGs(:,:) = 0.d0
 !    do is=1, nSubE
 !       do iv=1, nVert
 !          CoGs(x_,is) = CoGs(x_,is) + snode(x_,iv,is)
 !       enddo
 !       CoGs(:,is) = CoGs(:,is) / dble( nVert )
 !    enddo
 !    return
 !  end subroutine calculate__centerOfGravity
  

 !  ! ====================================================== !
 !  ! === subdivision of a triangular element            === !
 !  ! ====================================================== !
 !  subroutine subdivide__element( lnode, snode, nVert, nSubE, nSubdiv )
 !    implicit none
 !    integer         , intent(in)  :: nVert, nSubdiv
 !    double precision, intent(in)  :: lnode(3,nVert)
 !    double precision, intent(out) :: snode(3,nVert,nSubE)

 !    if      ( nSubdiv.eq.0 ) then
 !       do iv=1, nVert
 !          snode(:,iv,1) = lnode(:,iv)
 !       enddo
 !    else
 !       write(6,*) "not supported now."
 !    endif
 !    return
 !  end subroutine subdivide__element

  
 !  ! ====================================================== !
 !  ! === calculate MagnetiMoment from position          === !
 !  ! ====================================================== !
 !  function MagneticMoment( xm, ym, zm, bfdpos, mvec )
 !    implicit none
 !    double precision, intent(in) :: xm, ym, zm, bfdpos(3), mvec(3)
 !    double precision             :: rabs, mdotr, rvec3Inv, rvec(3), rhat(3)
 !    double precision             :: MagneticMoment
 !    double precision, parameter  :: fourpi = 16.d0*atan(1.d0)
 !    integer         , parameter  :: x_=1, y_=2, z_=3

 !    rvec(x_)          = bfdpos(x_) - xm
 !    rvec(y_)          = bfdpos(y_) - ym
 !    rvec(z_)          = bfdpos(z_) - zm
 !    rabs              = sqrt( rvec(1)**2 + rvec(2)**2 + rvec(3)**2 )
 !    rhat(:)           = rvec(:) / rabs
 !    rvec3Inv          = 1.d0 / ( fourpi * rabs**3 )
 !    mdotr             = mvec(1)*rhat(1) + mvec(2)*rhat(2) + mvec(3)*rhat(3)
 !    MagneticMoment    = rvec3Inv * ( 3.d0 * mdotr * rhat(z_) - mvec(z_) )
 !    return
 !  end function MagneticMoment
  

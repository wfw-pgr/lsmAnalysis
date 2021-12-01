module ioUtilityMod
contains
  
  ! ===================================================== !
  ! ===  Load Parameter File                          === !
  ! ===================================================== !
  subroutine load__parameterFile
    use variablesMod
    implicit none
    character(100)              :: name, type, TorF
    double precision, parameter :: unit_convert = 4.d0*atan( 1.d0 ) / 180.d0
    namelist /parameters/ iterMax, nDiv_B, threshold, wLaplace, wPicard_init, wPicard_last, &
         &                iPicard_iter, MzConstant, solverType, LSM_Engine, modelType, &
         &                rLim1, rLim2, pLim1, pLim2, zLim1, zLim2, resid_convergence
    
    ! ------------------------------------------------------ !
    ! --- [1] load namelist                              --- !
    ! ------------------------------------------------------ !
    open( lun, file=trim(listFile),status="old",form="formatted" )
    read( lun, nml=parameters )
    close(lun)
    
    ! ------------------------------------------------------ !
    ! --- [2] store parameters                           --- !
    ! ------------------------------------------------------ !
    iter    = 0
    mvec(1) = 0.d0
    mvec(2) = 0.d0
    mvec(3) = MzConstant
    wPicard = wPicard_init
    if ( pLim1.ge.pLim2  ) then
       write(6,*) "[ERROR] pLim1 > pLim2 :: ", pLim1, pLim2
       stop
    endif
    if ( pLim1.gt.180.d0 ) pLim1 = pLim1 - 360.d0
    if ( pLim2.gt.180.d0 ) pLim2 = pLim2 - 360.d0
    if ( pLim1.ge.pLim2  ) then ! quite import .ge. ( != .gt for 360 deg. )
       pJump = 360.d0
       pLim2 = pLim2 + pJump
    else
       pJump = 0.d0
    endif
    pLim1 = pLim1 * unit_convert
    pLim2 = pLim2 * unit_convert
    pJump = pJump * unit_convert
    write(6,*) 'pLim1, pLim2', pLim1, pLim2
    write(6,*) 'pJump', pJump

    if      ( modelType.eq."full" ) then
       Flag__fullModel = .true.
    else if ( modelType.eq."half" ) then
       Flag__fullModel = .false.
    else
       write(6,*) "[load__parameterFile]  ERROR !! modelType == ???  ", modelType
       stop
    endif
    write(6,*)
    write(6,*) ' solverType == ', solverType
    write(6,*) '  modelType == ', modelType
    write(6,*) ' LSM_engine == ', LSM_engine
    write(6,*)
    
    return
  end subroutine load__parameterFile


  ! ===================================================== !
  ! ===  Load Ideal / Background Bz File              === !
  ! ===================================================== !
  subroutine load__BFieldFile
    use variablesMod
    use utilitiesMod, only    : getNumLines
    implicit none
    integer, parameter       :: x_=1, y_=2, z_=3, i_=4, b_=5, s_=6, t_=7, e_=8
    integer                  :: k, nrow, ncol
    character(cLen)          :: cmt
    double precision         :: xd, yd, zd, bx, by, yCnt

    ! ------------------------------------------------------ !
    ! --- [1] Load Ideal Bz Data                         --- !
    ! ------------------------------------------------------ !
    write(6,'(a,1x,a,1x,a)',advance='no') '[load__BFieldFile]', 'loading ideal      BField File.... :: ', trim(BinpFile)
    open (lun,file=trim(BinpFile),form='formatted',status='old')
    read(lun,*)
    read(lun,*) cmt, nrow, ncol
    read(lun,*)
    nBpt = nrow
    allocate( BField(8,nBpt) )
    BField(:,:) = 0.d0
    do k=1, nBpt
       read(lun,*) BField(x_,k), BField(y_,k), BField(z_,k), BField(i_,k), BField(b_,k)
    enddo
    close(lun)
    write(6,*) '      [Done] '
    
    ! ------------------------------------------------------ !
    ! --- [2] shim, total, error                         --- !
    ! ------------------------------------------------------ !
    do k=1, nBpt
       BField(s_,k) = 0.d0
       BField(t_,k) = BField(b_,k) + BField(s_,k)
       BField(e_,k) = BField(i_,k) - BField(t_,k)
    enddo
    
    return
  end subroutine load__BFieldFile


  ! ====================================================== !
  ! === load weights for weighted Least Square         === !
  ! ====================================================== !
  subroutine load__weights
    use variablesMod
    implicit none
    integer                  :: k, nWeight
    character(cLen)          :: cmt
    integer, parameter       :: x_=1, y_=2, z_=3, cl_=4, wt_=5

    ! ------------------------------------------------------ !
    ! --- [1] Load weights                               --- !
    ! ------------------------------------------------------ !
    write(6,'(a,1x,a,1x,a)',advance='no') '[load__weights]', 'loading weights from :: ', trim(wghtFile)
    open (lun,file=trim(wghtFile),form='formatted',status='old')
    read(lun,*)
    read(lun,*) cmt, nWeight, cmt
    read(lun,*)
    allocate( weight_info(5,nWeight) )
    allocate( weights(nWeight) )
    weight_info(:,:) = 0.d0
    do k=1, nWeight
       read(lun,*) weight_info(:,k)
    enddo
    close(lun)
    do k=1, nWeight
       weights(k) = weight_info(wt_,k)
    enddo
    write(6,*) '      [Done] '
    
    return
  end subroutine load__weights
    

  ! ===================================================== !
  ! ===  Load Initial Magnet Position Data            === !
  ! ===================================================== !
  subroutine load__magnetShapeFile
    use variablesMod
    use utilitiesMod, only : getNumLines
    implicit none
    integer                     :: i, j, count, nComponents
    double precision            :: nodeNum, elemNum
    integer         , parameter :: x_=1, y_=2, z_=3, i_=4, s_=5, f_=6
    double precision, parameter :: onethird = 1.d0 / 3.d0
    double precision            :: v21(3), v31(3)
    character(cLen)             :: cmt

    ! ------------------------------------- !
    ! --- [1] Load Initial Magnet Data  --- !
    ! ------------------------------------- !
    write(6,'(a,1x,a,1x,a)',advance='no') '[load__MagnetShape]', 'loading Magnet Shape File.... :: ', trim(mshpFile)
    open (lun,file=trim(mshpFile),form='formatted',status='old')
    read(lun,*) cmt
    read(lun,*) cmt
    read(lun,*) cmt, LJ, LI, nComponents
    allocate( mshape(nComponents,LI,LJ) )
    nMpt          = LI*LJ
    mshape(:,:,:) = 0.d0
    do j=1, LJ
       do i=1, LI
          read(lun,*) mshape(:,i,j)
       enddo
    enddo
    close(lun)
    write(6,*) '      [Done] '
    ! ------------------------------------- !
    ! --- [2] count up cell to be used  --- !
    ! ------------------------------------- !
    ! idtable :: order of (i,j) cell ( count only active cell )
    allocate( idtable(0:LI+1,0:LJ+1) )
    count                  = 0
    idtable(0:LI+1,0:LJ+1) = 0
    do j=1, LJ
       do i=1, LI
          if ( mshape(f_,i,j).eq.1.d0 ) then
             count         = count + 1
             idtable(i,j)  = count
          else
             idtable(i,j)  = 0
          endif
       enddo
    enddo
    nNpt    = count
    ! ------------------------------------- !
    ! --- [3] xy Range Detection        --- !
    ! ------------------------------------- !
    xMin    = mshape(x_, 1, 1)
    xMax    = mshape(x_,LI, 1)
    yMin    = mshape(y_, 1, 1)
    yMax    = mshape(y_, 1,LJ)
    dx      = ( xMax-xMin ) / LI
    dy      = ( yMax-yMin ) / LJ
    ddx     = dx / nDiv_B
    ddy     = dy / nDiv_B
    write(6,*) ' (xMin,xMax) = ', xMin, xMax 
    write(6,*) ' (yMin,yMax) = ', yMin, yMax
    write(6,*) ' (  dx,  dy) = ',   dx, dy
    write(6,*) ' ( ddx, ddy) = ',  ddx, ddy
    ! ------------------------------------------------------ !
    ! --- [4] import initial shape                       --- !
    ! ------------------------------------------------------ !
    if ( Flag__initShape ) then
       write(6,'(a,1x,a,1x,a)',advance='no') '[load__MagnetShape]', 'loading Magnet Shape File.... :: ', trim(mshpFile)
       open (lun,file=trim(mshpFile),form='formatted',status='old')
       read(lun,*) cmt
       read(lun,*) cmt
       read(lun,*) cmt, LJ, LI, nComponents
       allocate( mshape(nComponents,LI,LJ) )
       nMpt          = LI*LJ
       mshape(:,:,:) = 0.d0
       do j=1, LJ
          do i=1, LI
             read(lun,*) mshape(:,i,j)
          enddo
       enddo
       close(lun)
       write(6,*) '      [Done] '
    endif
    
    return
  end subroutine load__magnetShapeFile


  ! ===================================================== !
  ! ===  write result into File                       === !
  ! ===================================================== !
  subroutine save__result
    use variablesMod 
    implicit none
    integer          :: i, j, k
    double precision :: ratio
    character(  4)   :: citer
    character(300)   :: FileName
    write(citer,'(i4.4)') iter

    ! ------------------------------------- !
    ! --- [1] BField output             --- !
    ! ------------------------------------- !
    FileName = 'out/bfield_' // citer // '.dat'
    open (lun,file=trim(FileName),form='formatted',status='replace')
    write(lun,*) '# xb yb zb ideal background shim total error'
    write(lun,'(a,i12)') '# 8 ', nBpt
    do k=1, nBpt
       write(lun,'(8(e15.8,1x))') BField(1:8,k)
    enddo
    close(lun)
    write(6,'(a,a)') '[save__result] outFile :: ', trim(FileName)
    ! ------------------------------------- !
    ! --- [2] shape output              --- !
    ! ------------------------------------- !
    FileName = 'out/mshape_' // citer // '.dat'
    open (lun,file=trim(FileName),form='formatted',status='replace')
    write(lun,*) "# x_ y_ z_ i_ s_ f_"
    write(lun,*) "#", LI*LJ, 6
    write(lun,*) "#", 1, LJ, LI, 6
    do j=1, LJ
       do i=1, LI
          write(lun,'(6(e15.8,1x))') mshape(1:6,i,j)
       enddo
    enddo
    close(lun)
    write(6,'(a,a)') '[save__result] outFile :: ', trim(FileName)
    ! ------------------------------------- !
    ! --- [3] spectrum output           --- !
    ! ------------------------------------- !
    FileName = 'out/spectr_' // citer // '.dat'
    open (lun,file=trim(FileName),form='formatted',status='replace')
    do k=1, min(nNpt,nBpt)
       if ( singular(k).eq.0.d0 ) then
          ratio = 0.d0
       else
          ratio = spectrum(k)/singular(k)
       endif
       write(lun,'(i8,1x,3(e15.8,1x))') k, singular(k), spectrum(k), ratio
    enddo
    close(lun)
    write(6,'(a,a)') '[save__result] outFile :: ', trim(FileName)

    return
  end subroutine save__result


  ! ===================================================== !
  ! ===  write Error Field Data in File               === !
  ! ===================================================== !
  subroutine save__residual
    use variablesMod
    implicit none
    integer, parameter :: t_=7, e_=8
    double precision   :: avg, std, min, max
    character(cLen)    :: fmt              = '(i8,4(1x,e12.5))'
    logical, save      :: flag__initialize = .true.
    
    ! --------------------------------------------------------- !
    ! --- [1] Calculate Average & Standard Deviation        --- !
    ! --------------------------------------------------------- !
    avg =         sum( BField(e_,:)    )            / dble( nBpt )
    std = sqrt( ( sum( BField(e_,:)**2 ) - avg**2 ) / dble( nBpt ) )
    min =      minval( BField(e_,:) )
    max =      maxval( BField(e_,:) )
    ! --------------------------------------------------------- !
    ! --- [2] write down Index of File                      --- !
    ! --------------------------------------------------------- !
    if ( flag__initialize ) then
       open (lun,file=trim(BresFile),status='replace',form='formatted')
       write(lun,'(a)') '# iter avg(dBz) std(dBz) min(dBz) max(dBz)'
       close(lun)
       flag__initialize = .false.
    end if
    ! --------------------------------------------------------- !
    ! --- [3] Record Residual of B-Field                    --- !
    ! --------------------------------------------------------- !
    minRes(3) = minRes(2)
    maxRes(3) = maxRes(2)
    minRes(2) = minRes(1)
    maxRes(2) = maxRes(1)
    minRes(1) = min
    maxRes(1) = max
    open (lun,file=trim(BresFile),status='old',position='append',form='formatted')
    write(lun,trim(fmt)) iter, avg, std, min, max
    close(lun)
    ! --------------------------------------------------------- !
    ! --- [4] writeResidual                                 --- !
    ! --------------------------------------------------------- !
    write(6,'(a,i8,4(1x,e12.5))') '[save__residual] iter, avg, std, min, max == ', iter, avg, std, min, max
    
    return
  end subroutine save__residual


  ! ===================================================== !
  ! ===  write Error Field Data in File               === !
  ! ===================================================== !
  subroutine save__residual_WLS
    use variablesMod
    implicit none
    integer, parameter :: t_=7, e_=8, cl_=4
    integer            :: i1, i2, i3, i4, ik, ic
    double precision   :: avg1, std1, min1, max1, avg2, std2, min2, max2
    character(cLen)    :: fmt              = '(2(i8,1x),4(e15.8,1x))'
    character(cLen)    :: cmt, resFile
    character(2)       :: cnum
    logical, save      :: flag__initialize = .true.
    integer, save      :: maxcolor = 1
    
    ! --------------------------------------------------------- !
    ! --- [1] write down Index of File / initialize         --- !
    ! --------------------------------------------------------- !
    if ( flag__initialize ) then
       do ik=1, nBpt
          maxcolor = max( maxcolor, int( weight_info(ik,cl_) ) )
       enddo
       allocate( avgs(maxcolor), stds(maxcolor) )
       allocate( sumCnt(maxcolor), sumErr(maxcolor), sumErr2(maxcolor), minErr(maxcolor), maxErr(maxcolor) )
       do ic=1, maxcolor
          write(cnum,"(i2.2)") ic
          resFile = "dat/bfield_resid_" // cnum // ".dat"
          open (lun,file=trim(resFile),status='replace',form='formatted')
          write(lun,'(a)') '# ic iter avg(dBz) std(dBz) min(dBz) max(dBz)'
          close(lun)
       enddo
       flag__initialize = .false.
    end if
    
    ! --------------------------------------------------------- !
    ! --- [2] Calculate Average & Standard Deviation        --- !
    ! --------------------------------------------------------- !

    sumCnt (:) =   0
    sumErr (:) =   0.d0
    sumErr2(:) =   0.d0
    minErr (:) =   1.d10
    maxErr (:) = - 1.d10
    do ik=1, nBpt
       ic          = weight_info( ik, cl_ )
       sumCnt (ic) =      sumCnt (ic) + 1
       sumErr (ic) =      sumErr (ic) + BField(e_,ik)
       sumErr2(ic) =      sumErr2(ic) + BField(e_,ik)**2
       minErr (ic) = min( minErr (ic),  BField(e_,ik) )
       maxErr (ic) = max( maxErr (ic),  BField(e_,ik) )
    enddo
    do ic=1, maxcolor
       avgs(ic) =         sumErr (ic)                 / dble( sumCnt(ic) )
       stds(ic) = sqrt( ( sumErr2(ic) - avgs(ic)**2 ) / dble( sumCnt(ic) ) )
    enddo
    
    ! --------------------------------------------------------- !
    ! --- [3] Record Residual of B-Field                    --- !
    ! --------------------------------------------------------- !
    write(6  ,*)
    write(6  ,*) '[save__residual] pole :: iter, avg, std, min, max == '
    do ic=1, maxcolor
       write(cnum,"(i2.2)") ic
       resFile = "dat/bfield_resid_" // cnum // ".dat"
       open (lun,file=trim(resFile),status='old',position='append',form='formatted')
       write(lun,trim(fmt)) ic, iter, avgs(ic), stds(ic), minErr(ic), maxErr(ic)
       close(lun)
       write(6  ,trim(fmt)) ic, iter, avgs(ic), stds(ic), minErr(ic), maxErr(ic)
    enddo
    write(6  ,*)
    
    return
  end subroutine save__residual_WLS


  ! ====================================================== !
  ! === check too slow conversion                      === !
  ! ====================================================== !
  subroutine check__residConvergence
    use variablesMod
    implicit none
    double precision :: chg12, chg23, chg13, within

    ! ------------------------------------------------------ !
    ! --- [1] check maximum change                       --- !
    ! ------------------------------------------------------ !

    if ( abs( maxRes(1) ).ge.abs( minRes(1) ) ) then
       chg12  = abs( abs( maxRes(1) ) - abs( maxRes(2) ) )
       chg23  = abs( abs( maxRes(2) ) - abs( maxRes(3) ) )
       chg13  = abs( abs( maxRes(1) ) - abs( maxRes(3) ) )
       within = max( max( chg12, chg23 ), chg13 )
    else
       chg12  = abs( abs( minRes(1) ) - abs( minRes(2) ) )
       chg23  = abs( abs( minRes(2) ) - abs( minRes(3) ) )
       chg13  = abs( abs( minRes(1) ) - abs( minRes(3) ) )
       within = max( max( chg12, chg23 ), chg13 )
    endif

    ! ------------------------------------------------------ !
    ! --- [2] judge maximum change                       --- !
    ! ------------------------------------------------------ !
    if ( within.lt.resid_convergence ) then
       flag__exitStatus = .true.
    endif
    
    return
  end subroutine check__residConvergence

  
end module ioUtilityMod



    
    ! avg1 = 0.d0
    ! std1 = 0.d0
    ! min1 = + 1.d10
    ! max1 = - 1.d10
    ! do k=1, nBpt_pole
    !    avg1 = avg1 + BField(e_,k)
    !    std1 = std1 + BField(e_,k)**2
    !    min1 = min( min1, BField(e_,k) )
    !    max1 = max( max1, BField(e_,k) )
    ! enddo
    ! avg2 = 0.d0
    ! std2 = 0.d0
    ! min2 = + 1.d10
    ! max2 = - 1.d10
    ! do k=nBpt_pole, nBpt_peel
    !    avg2 = avg2 + BField(e_,k)
    !    std2 = std2 + BField(e_,k)**2
    !    min2 = min( min2, BField(e_,k) )
    !    max2 = max( max2, BField(e_,k) )
    ! enddo
    ! avg1 = avg1 / dble( nBpt_pole )
    ! avg2 = avg2 / dble( nBpt_pole )
    ! avg1 = avg1 / dble( nBpt_pole )
    ! avg1 = avg1 / dble( nBpt_pole )




  ! ! ====================================================== !
  ! ! === load weights for weighted Least Square         === !
  ! ! ====================================================== !
  ! subroutine load__weights
  !   use variablesMod
  !   implicit none
  !   integer                  :: k, nrow
  !   character(cLen)          :: cmt
  !   integer, parameter       :: cl_=1, wt_=2
  !   double precision, allocatable :: colors(:)

  !   ! ------------------------------------------------------ !
  !   ! --- [1] Load weights                               --- !
  !   ! ------------------------------------------------------ !
  !   write(6,'(a,1x,a,1x,a)',advance='no') '[load__weights]', 'loading weights from :: ', trim(wghtFile)
  !   open (lun,file=trim(wghtFile),form='formatted',status='old')
  !   read(lun,*)
  !   read(lun,*) cmt, nrow, cmt
  !   read(lun,*)
  !   allocate( colors(nrow), weights(nrow) )
  !   colors (:) = 0.d0
  !   weights(:) = 0.d0
  !   do k=1, nrow
  !      read(lun,*) colors(k), weights(k)
  !   enddo
  !   close(lun)
  !   write(6,*) '      [Done] '
    
  !   return
  ! end subroutine load__weights




  ! ! ===================================================== !
  ! ! ===  write Error Field Data in File               === !
  ! ! ===================================================== !
  ! subroutine save__residual_WLS
  !   use variablesMod
  !   implicit none
  !   integer, parameter :: t_=7, e_=8
  !   integer            :: i1, i2, i3, i4
  !   double precision   :: avg1, std1, min1, max1, avg2, std2, min2, max2
  !   character(cLen)    :: fmt              = '(i8,8(1x,e12.5))'
  !   character(cLen)    :: cmt
  !   logical, save      :: flag__initialize = .true.


  !   ! --------------------------------------------------------- !
  !   ! --- [1] write down Index of File / initialize         --- !
  !   ! --------------------------------------------------------- !
  !   if ( flag__initialize ) then
  !      open (lun,file=trim(BresFile),status='replace',form='formatted')
  !      write(lun,'(a)') '# iter pole:avg(dBz) std(dBz) min(dBz) max(dBz) peel:avg(dBz) std(dBz) min(dBz) max(dBz) '
  !      close(lun)
  !      flag__initialize = .false.
  !      open (lun,file=trim(wcnfFile),status='old')
  !      read (lun,*) cmt, cmt, nBpt_pole
  !      read (lun,*) cmt, cmt, nBpt_peel
  !      close(lun)
  !   end if
  !   ! --------------------------------------------------------- !
  !   ! --- [2] Calculate Average & Standard Deviation        --- !
  !   ! --------------------------------------------------------- !
  !   i1   = 1
  !   i2   = nBpt_pole
  !   i3   = nBpt_pole + 1
  !   i4   = nBpt_pole + nBpt_peel
  !   avg1 =         sum( BField(e_,i1:i2)    )             / dble( nBpt )
  !   std1 = sqrt( ( sum( BField(e_,i1:i2)**2 ) - avg1**2 ) / dble( nBpt ) )
  !   min1 =      minval( BField(e_,i1:i2)    )
  !   max1 =      maxval( BField(e_,i1:i2)    )
  !   avg2 =         sum( BField(e_,i3:i4)    )             / dble( nBpt )
  !   std2 = sqrt( ( sum( BField(e_,i3:i4)**2 ) - avg2**2 ) / dble( nBpt ) )
  !   min2 =      minval( BField(e_,i3:i4)    )
  !   max2 =      maxval( BField(e_,i3:i4)    )
  !   ! --------------------------------------------------------- !
  !   ! --- [3] Record Residual of B-Field                    --- !
  !   ! --------------------------------------------------------- !
  !   open (lun,file=trim(BresFile),status='old',position='append',form='formatted')
  !   write(lun,trim(fmt)) iter, avg1, std1, min1, max1, avg2, std2, min2, max2
  !   close(lun)
  !   ! --------------------------------------------------------- !
  !   ! --- [4] writeResidual                                 --- !
  !   ! --------------------------------------------------------- !
  !   write(6,'(a,i8,4(1x,e12.5))') '[save__residual] pole :: iter, avg, std, min, max == ', iter, avg1, std1, min1, max1
  !   write(6,'(a,i8,4(1x,e12.5))') '[save__residual] peel :: iter, avg, std, min, max == ', iter, avg2, std2, min2, max2
    
  !   return
  ! end subroutine save__residual_WLS


    ! ! ------------------------------------------------------ !
    ! ! --- [1] Load parameter.dat                         --- !
    ! ! ------------------------------------------------------ !
    ! write(6,*)
    ! write(6,'(a,1x,a,1x,a)',advance='no') '[load__parameterFile]', 'loading parameter File.... :: ', trim(parmFile)
    ! open (lun,file=trim(parmFile),status='old',form='formatted')
    ! read (lun,*) name, type, iterMax
    ! read (lun,*) name, type, nDiv_B
    ! read (lun,*) name, type, threshold
    ! read (lun,*) name, type, wLaplace
    ! read (lun,*) name, type, wPicard_init
    ! read (lun,*) name, type, wPicard_last
    ! read (lun,*) name, type, iPicard_iter
    ! read (lun,*) name, type, MzConstant
    ! read (lun,*) name, type, solverType
    ! read (lun,*) name, type, LSM_engine
    ! read (lun,*) name, type, rLim1
    ! read (lun,*) name, type, rLim2
    ! read (lun,*) name, type, pLim1
    ! read (lun,*) name, type, pLim2
    ! read (lun,*) name, type, zLim1
    ! read (lun,*) name, type, zLim2
    ! read (lun,*) name, type, modelType
    ! read (lun,*) name, type, resid_convergence
    ! close(lun)
    ! write(6,*) '      [Done] '


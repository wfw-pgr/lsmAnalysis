module ioUtilityMod
contains
  
  ! ===================================================== !
  ! ===  Load Parameter File                          === !
  ! ===================================================== !
  subroutine load__parameterFile
    use variablesMod
    implicit none
    namelist /parameters/ iterMax, nSubdiv, nDiv_z, MzConst, zLim1, zLim2, iPicard, wPicard, &
         &                resid__criterion, lsm__engine, solverType, convergenceType

    ! ------------------------------------------------------ !
    ! --- [1] default value                              --- !
    ! ------------------------------------------------------ !
    iPicard(:) = -1
    wPicard(:) = 0.d0
    
    ! ------------------------------------------------------ !
    ! --- [2] load namelist                              --- !
    ! ------------------------------------------------------ !
    write(6,'(a,1x,a25)',advance='no') '[load__parameterFile] loading :: ', &
         &                             adjustl(trim(listFile))
    open( lun, file=trim(listFile),status="old",form="formatted" )
    read( lun, nml=parameters )
    close(lun)
    write(6,"(3x,a)") "[Done]"

    ! ------------------------------------------------------ !
    ! --- [3] post process                               --- !
    ! ------------------------------------------------------ !
    
    ! ------------------------------------------------------ !
    ! --- [4] display parameters                         --- !
    ! ------------------------------------------------------ !
    
    return
  end subroutine load__parameterFile


  ! ====================================================== !
  ! ===  Load ideal / background Bz File               === !
  ! ====================================================== !
  subroutine load__bfieldFile
    use variablesMod
    implicit none
    integer, parameter       :: x_=1, y_=2, z_=3, i_=4, b_=5, s_=6, t_=7, e_=8
    integer                  :: k, nrow, ncol
    character(cLen)          :: cmt
    double precision         :: xd, yd, zd, bx, by, yCnt

    ! ------------------------------------------------------ !
    ! --- [1] Load ideal / background Bz Data            --- !
    ! ------------------------------------------------------ !
    write(6,'(a,1x,a25)',advance='no') '[load__bfieldFile] loading :: ', &
         &                             adjustl( trim(binpFile) )
    open (lun,file=trim(binpFile),form='formatted',status='old')
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
    write(6,"(3x,a)") "[Done]"
    
    ! ------------------------------------------------------ !
    ! --- [2] shim, total, error                         --- !
    ! ------------------------------------------------------ !
    do k=1, nBpt
       BField(s_,k) = 0.d0
       BField(t_,k) = BField(b_,k) + BField(s_,k)
       BField(e_,k) = BField(i_,k) - BField(t_,k)
    enddo
    
    return
  end subroutine load__bfieldFile

  
  ! ====================================================== !
  ! === load__mshapeFile                               === !
  ! ====================================================== !
  subroutine load__mshapeFile
    use variablesMod
    use lstStructMod
    implicit none
    integer                     :: ie, nrow, ncol
    integer        , parameter  :: nComps = 7
    character(cLen)             :: cmt

    ! ------------------------------------------------------ !
    ! --- [1] load initial magnet pole Data              --- !
    ! ------------------------------------------------------ !
    write(6,'(a,1x,a25)',advance='no') '[load__mshapeFile] loading :: ', &
         &                              adjustl( trim(mshpFile) )
    open (lun,file=trim(mshpFile),form='formatted',status='old')
    read(lun,*) cmt
    read(lun,*) cmt, nrow, ncol
    read(lun,*) cmt
    if ( nrow.ne.nElems ) stop "[load__mshapeFile] incompatible nElems.... [ERROR]"
    if ( ncol.ne.nComps ) stop "[load__mshapeFile] incompatible nComps.... [ERROR]"
    allocate( mshape(nComps,nElems) )
    do ie=1, nElems
       read(lun,*) mshape(:,ie)
    enddo
    close(lun)
    write(6,"(3x,a)") "[Done]"

    ! ------------------------------------------------------ !
    ! --- [2] extract group info.                        --- !
    ! ------------------------------------------------------ !
    nullify( groupList )
    do ie=1, nElems
       call add__elementInList( groupList, elementNum=ie, groupNum=int(mshape(mg_,ie)) )
    enddo
    
    call investigate__listInfo( groupList, max_nCell, nGroups )
    
    allocate( groupNums(nGroups), groupedCells(max_nCell) )
    write(6,*)
    write(6,"(a)"    ) "[load__mshapeFile]  --------  *  grouping  *  ---------"
    write(6,"(a,i10)") "[load__mshapeFile]        nGroup :: ", nGroups
    write(6,"(a,i10)") "[load__mshapeFile]     max_nCell :: ", max_nCell
    write(6,"(a)"    ) "[load__mshapeFile]  --------  *  grouping  *  ---------"
    write(6,*)

    call obtain__groupNumArray( groupList, groupNums, nGroups )
    
    return
  end subroutine load__mshapeFile


  ! ====================================================== !
  ! === load weights for weighted Least Square         === !
  ! ====================================================== !
  subroutine load__weightFile
    use variablesMod
    implicit none
    integer                     :: ik, nrow, ncol
    character(cLen)             :: cmt
    integer         , parameter :: nComps = 5

    ! ------------------------------------------------------ !
    ! --- [1] Load weights                               --- !
    ! ------------------------------------------------------ !
    !  -- [1-1]  header                                  --  !
    write(6,'(a,1x,a25)',advance='no') '[load__weightFile] loading :: ', &
         &                             adjustl( trim(wghtFile) )
    open (lun,file=trim(wghtFile),form='formatted',status='old')
    read(lun,*)
    read(lun,*) cmt, nrow, ncol
    read(lun,*)
    ! if ( nrow.ne.nBpt+nElems ) stop "[load__mshapeFile] incompatible nElems.... [ERROR]"
    if ( ncol.ne.nComps ) stop "[load__mshapeFile] incompatible nComps.... [ERROR]"
    !  -- [1-2]  initialize weights                      --  !
    allocate( weights(nComps,nrow) )
    weights(:,:) = 0.d0
    !  -- [1-3]  read weights from file                  --  !
    do ik=1, nrow
       read(lun,*) weights(:,ik)
    enddo
    close(lun)
    write(6,"(3x,a)") "[Done]"
    
    return
  end subroutine load__weightFile
  

  ! ====================================================== !
  ! === load elems and nodes                           === !
  ! ====================================================== !
  subroutine load__ElemsNodes
    use variablesMod
    implicit none
    character(cLen) :: cmt
    integer         :: ik, nComp

    ! ------------------------------------------------------ !
    ! --- [1] load elements                              --- !
    ! ------------------------------------------------------ !
    write(6,'(a,1x,a25)',advance='no') '[load__elemsNodes] loading :: ', &
         &                             adjustl( trim(elemFile) )
    open (lun,file=trim(elemFile),form="formatted",status="old")
    read(lun,*) cmt
    read(lun,*) cmt, nElems, nComp
    read(lun,*) cmt, nElems, nComp
    allocate( elems(nComp,nElems) )
    elems(:,:) = 0
    do ik=1, nElems
       read(lun,*) elems(:,ik)
    enddo
    close(lun)
    ! -- plus 1 -- !
    elems(:,:) = elems(:,:) + 1
    
    write(6,"(3x,a)") "[Done]"
    
    ! ------------------------------------------------------ !
    ! --- [2] load nodes                                 --- !
    ! ------------------------------------------------------ !
    write(6,'(a,1x,a25)',advance='no') '[load__elemsNodes] loading :: ', &
         &                              adjustl( trim(nodeFile) )
    open (lun,file=trim(nodeFile),form="formatted",status="old")
    read(lun,*) cmt
    read(lun,*) cmt, nNodes, nComp
    read(lun,*) cmt, nNodes, nComp
    allocate( nodes(nComp,nNodes) )
    nodes(:,:) = 0.d0
    do ik=1, nNodes
       read(lun,*) nodes(:,ik)
    enddo
    close(lun)
    write(6,"(3x,a)") "[Done]"
    
    return
  end subroutine load__ElemsNodes


  ! ====================================================== !
  ! === save result in files                           === !
  ! ====================================================== !
  subroutine save__results
    use variablesMod 
    implicit none
    integer          :: ik
    character(  4)   :: citer
    character(cLen)  :: FileName

    ! ------------------------------------------------------ !
    ! --- [0] preparation                                --- !
    ! ------------------------------------------------------ !
    write(citer,'(i4.4)') iter

    ! ------------------------------------------------------ !
    ! --- [1] save bfield                                --- !
    ! ------------------------------------------------------ !
    FileName = 'out/bfield_' // citer // '.dat'
    write(6,'(a,1x,a25)',advance='no') '[save__result] saving :: ', trim(FileName)
    open (lun,file=trim(FileName),form='formatted',status='replace')
    write(lun,*) '# xb yb zb ideal background shim total error'
    write(lun,'(a,1x,i8)') '# 8 ', nBpt
    write(lun,'(a,1x,i8)') '# 8 ', nBpt
    do ik=1, nBpt
       write(lun,'(8(e15.8,1x))') bfield(:,ik)
    enddo
    close(lun)
    write(6,"(3x,a)") "[Done]"

    ! ------------------------------------------------------ !
    ! --- [2] save mshape                                --- !
    ! ------------------------------------------------------ !
    FileName = 'out/mshape_' // citer // '.dat'
    write(6,'(a,1x,a25)',advance='no') '[save__result] saving :: ', trim(FileName)
    open (lun,file=trim(FileName),form='formatted',status='replace')
    write(lun,*) "# xm_ ym_ zm_ mi_ ms_ mf_"
    write(lun,'(a,2(1x,i8))') "# ", nElems, 7
    write(lun,'(a,2(1x,i8))') "# ", nElems, 7
    do ik=1, nElems
       write(lun,'(7(e15.8,1x))') mshape(:,ik)
    enddo
    close(lun)
    write(6,"(3x,a)") "[Done]"

    return
  end subroutine save__results
  
  
  ! ====================================================== !
  ! ===  write Error Field Data in a file              === !
  ! ====================================================== !
  subroutine save__residuals
    use variablesMod
    implicit none
    integer            :: ic, ik
    character(cLen)    :: fmt1             = '(2(a4,1x),5(a12  ,1x))'
    character(cLen)    :: fmt2             = '(2(i4,1x),5(e12.5,1x))'
    character(cLen)    :: fmt3             = '(2(i6,1x),5(e15.8,1x))'
    logical, save      :: flag__initialize = .true.
    character(2)       :: citer
    
    ! --------------------------------------------------------- !
    ! --- [1] write down Index of File / initialize         --- !
    ! --------------------------------------------------------- !
    if ( flag__initialize ) then
       allocate( avgs  (nColor), stds  (nColor), rmse   (nColor) )
       allocate( sumCnt(nColor), sumErr(nColor), sumErr2(nColor) )
       allocate( minErr(nColor), maxErr(nColor) )
       open (lun,file=trim(bresFile),status='replace',form='formatted')
       write(lun,"(a)") '# iter ic avg(dBz) std(dBz) min(dBz) max(dBz) rmse(dBz)'
       close(lun)
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
       ic          = nint( weights( cl_, ik ) )
       sumCnt (ic) =      sumCnt (ic) + 1
       sumErr (ic) =      sumErr (ic) + bfield(be_,ik)
       sumErr2(ic) =      sumErr2(ic) + bfield(be_,ik)**2
       minErr (ic) = min( minErr (ic),  bfield(be_,ik) )
       maxErr (ic) = max( maxErr (ic),  bfield(be_,ik) )
    enddo
    do ic=1, nColor
       avgs(ic)    =         sumErr (ic)                 / dble( sumCnt(ic) )
       stds(ic)    = sqrt( ( sumErr2(ic) - avgs(ic)**2 ) / dble( sumCnt(ic) ) )
       rmse(ic)    = sqrt( ( sumErr2(ic)               ) / dble( sumCnt(ic) ) )
    enddo

    ! ------------------------------------------------------ !
    ! --- [3] record residual of B-Field                 --- !
    ! ------------------------------------------------------ !
    !  -- (2) -> (3)  -- !
    minRes(3) = minRes(2)
    maxRes(3) = maxRes(2)
    rmsRes(3) = rmsRes(2)
    !  -- (1) -> (2)  -- !
    minRes(2) = minRes(1)
    maxRes(2) = maxRes(1)
    rmsRes(2) = rmsRes(1)
    !  -- new -> (1)  -- !
    minRes(1) = minErr(1)
    maxRes(1) = maxErr(1)
    rmsRes(1) = rmse  (1)

    ! ------------------------------------------------------ !
    ! --- [4] display & save in a file                   --- !
    ! ------------------------------------------------------ !
    write(6,*)
    write(6,*         ) "[save__residuals] Residuals :: "
    write(6,trim(fmt1)) "iter", "ic", "avg", "std", "min", "max", "rmse"
    open (lun,file=trim(bresFile),status="old",position="append",form="formatted")
    do ic=1, nColor
       write(  6,trim(fmt2)) iter, ic, avgs(ic), stds(ic), minErr(ic), maxErr(ic), rmse(ic)
       write(lun,trim(fmt3)) iter, ic, avgs(ic), stds(ic), minErr(ic), maxErr(ic), rmse(ic)
    enddo
    close(lun)
    
    return
  end subroutine save__residuals


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
    
    if ( trim(convergenceType).eq."maxmin" ) then
       if ( abs( maxRes(1) ).ge.abs( minRes(1) ) ) then
          ! -- evaluate max -- !
          chg12  = abs( abs( maxRes(1) ) - abs( maxRes(2) ) )
          chg23  = abs( abs( maxRes(2) ) - abs( maxRes(3) ) )
          chg13  = abs( abs( maxRes(1) ) - abs( maxRes(3) ) )
          within = max( max( chg12, chg23 ), chg13 )
       else
          ! -- evaluate min -- !
          chg12  = abs( abs( minRes(1) ) - abs( minRes(2) ) )
          chg23  = abs( abs( minRes(2) ) - abs( minRes(3) ) )
          chg13  = abs( abs( minRes(1) ) - abs( minRes(3) ) )
          within = max( max( chg12, chg23 ), chg13 )
       endif
       
    else if ( trim(convergenceType).eq."min" ) then
       ! -- evaluate min -- !
       chg12  = abs( abs( minRes(1) ) - abs( minRes(2) ) )
       chg23  = abs( abs( minRes(2) ) - abs( minRes(3) ) )
       chg13  = abs( abs( minRes(1) ) - abs( minRes(3) ) )
       within = max( max( chg12, chg23 ), chg13 )

    else if ( trim(convergenceType).eq."max" ) then
       ! -- evaluate max -- !
       chg12  = abs( abs( maxRes(1) ) - abs( maxRes(2) ) )
       chg23  = abs( abs( maxRes(2) ) - abs( maxRes(3) ) )
       chg13  = abs( abs( maxRes(1) ) - abs( maxRes(3) ) )
       within = max( max( chg12, chg23 ), chg13 )

    else if ( trim(convergenceType).eq."rmse" ) then
       ! -- evaluate rms -- !
       chg12  = abs( abs( rmsRes(1) ) - abs( rmsRes(2) ) )
       chg23  = abs( abs( rmsRes(2) ) - abs( rmsRes(3) ) )
       chg13  = abs( abs( rmsRes(1) ) - abs( rmsRes(3) ) )
       within = max( max( chg12, chg23 ), chg13 )

    else if ( trim(convergenceType).eq."none" ) then
       within = 10.d0 * resid__criterion
       
    else
       write(6,*) "[check__residConvergence] unrecognized convergenceType :: ", &
            &     trim( convergenceType ), " [ERROR]"
       stop
    endif
       
    ! ------------------------------------------------------ !
    ! --- [2] judge converged or not                     --- !
    ! ------------------------------------------------------ !
    write(6,*)
    write(6,"(a,a15)"  ) "[check__residConvergence] convergenceType   :: ", trim(convergenceType)
    write(6,"(a,e15.8)") "[check__residConvergence] Max.chg. (within) :: ", within
    write(6,"(a,e15.8)") "[check__residConvergence] resid. criterion  :: ", resid__criterion
    write(6,*)
    if ( ( iter.ge.3 ).and.( within.lt.resid__criterion ) ) then
       flag__exitStatus = .true.
    endif

    return
  end subroutine check__residConvergence


end module ioUtilityMod





    ! ! ------------------------------------------------------ !
    ! ! --- [2] store parameters                           --- !
    ! ! ------------------------------------------------------ !
    ! if ( pLim1.ge.pLim2  ) then
    !    write(6,*) "[ERROR] pLim1 > pLim2 :: ", pLim1, pLim2
    !    stop
    ! endif
    ! if ( pLim1.gt.180.d0 ) pLim1 = pLim1 - 360.d0
    ! if ( pLim2.gt.180.d0 ) pLim2 = pLim2 - 360.d0
    ! if ( pLim1.ge.pLim2  ) then ! quite import .ge. ( != .gt for 360 deg. )
    !    pJump = 360.d0
    !    pLim2 = pLim2 + pJump
    ! else
    !    pJump = 0.d0
    ! endif
    ! pLim1 = pLim1 * unit_convert
    ! pLim2 = pLim2 * unit_convert
    ! pJump = pJump * unit_convert
    ! write(6,*) 'pLim1, pLim2', pLim1, pLim2
    ! write(6,*) 'pJump', pJump

    ! if      ( modelType.eq."full" ) then
    !    Flag__fullModel = .true.
    ! else if ( modelType.eq."half" ) then
    !    Flag__fullModel = .false.
    ! else
    !    write(6,*) "[load__parameterFile]  ERROR !! modelType == ???  ", modelType
    !    stop
    ! endif
    ! write(6,*)
    ! write(6,*) ' solverType == ', solverType
    ! write(6,*) '  modelType == ', modelType
    ! write(6,*) ' LSM_engine == ', LSM_engine
    ! write(6,*)




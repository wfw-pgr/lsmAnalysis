module ioUtilityMod
contains

  ! ===================================================== !
  ! ===  Load Parameter File                          === !
  ! ===================================================== !
  subroutine load__parameterFile
    use variablesMod
    use utilitiesMod
    implicit none
    namelist /parameters/ iterMax

    call print__section( "load__parameterFile", "=", 19, 4, 70 )
    
    ! ------------------------------------------------------ !
    ! --- [1] load namelist                              --- !
    ! ------------------------------------------------------ !
    open( lun, file=trim(listFile),status="old",form="formatted" )
    read( lun, nml=parameters )
    close(lun)

    ! ------------------------------------------------------ !
    ! --- [2] post process                               --- !
    ! ------------------------------------------------------ !

    ! ------------------------------------------------------ !
    ! --- [3] display parameters                         --- !
    ! ------------------------------------------------------ !
    
    
    return
  end subroutine load__parameterFile


  ! ====================================================== !
  ! ===  Load ideal / background Bz File               === !
  ! ====================================================== !
  subroutine load__BFieldFile
    use variablesMod
    use utilitiesMod, only    : getNumLines
    implicit none
    integer, parameter       :: x_=1, y_=2, z_=3, i_=4, b_=5, s_=6, t_=7, e_=8
    integer                  :: k, nrow, ncol
    character(cLen)          :: cmt
    double precision         :: xd, yd, zd, bx, by, yCnt
    
    ! ------------------------------------------------------ !
    ! --- [1] Load ideal / background Bz Data            --- !
    ! ------------------------------------------------------ !
    write(6,"(a,1x,a)",advance="no") "[load__BFieldFile] loading ideal      BField File.... :: ",&
         &                           trim(binpFile)
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
  ! === load__mshapeFile              === !
  ! ====================================================== !
  subroutine load__mshapeFile
    implicit none
    return
  end subroutine load__mshapeFile
  

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
    write(6,"(a)",advance="no") "[load__elemsNodes] loading dat/elems.dat "
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
    write(6,*) "      [Done] "
    
    ! ------------------------------------------------------ !
    ! --- [2] load nodes                                 --- !
    ! ------------------------------------------------------ !
    write(6,"(a)",advance="no") "[load__elemsNodes] loading dat/nodes.dat "
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
    write(6,*) "      [Done] "
    
    return
  end subroutine load__ElemsNodes

  
end module ioUtilityMod





    ! ! ------------------------------------------------------ !
    ! ! --- [2] store parameters                           --- !
    ! ! ------------------------------------------------------ !
    ! iter    = 0
    ! mvec(1) = 0.d0
    ! mvec(2) = 0.d0
    ! mvec(3) = MzConstant
    ! wPicard = wPicard_init
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

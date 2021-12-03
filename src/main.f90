program main
  use variablesMod
  use utilitiesMod
  use ioUtilityMod
  use allocatorMod
  implicit none
  
  ! ------------------------------------------------------ !
  ! --- [1] pre-process                                --- !
  ! ------------------------------------------------------ !
  call TimeMeasure( 0 )
  call show__programLogo
  call load__parameterFile
  call load__BFieldFile
  call load__mshapeFile
  call load__ElemsNodes
  call TimeMeasure( 1 )

  
  
  ! ------------------------------------------------------ !
  ! --- [3] post-process                               --- !
  ! ------------------------------------------------------ !
  call deallocate__variables
  call show__endLogo
  call TimeMeasure( -1 )
  
end program main



  ! use variablesMod
  ! use utilitiesMod
  ! use ioUtilityMod
  ! use allocatorMod
  ! use updateEqsMod
  ! use invSolverMod
  ! use debugToolMod
  ! implicit none

  ! ! ------------------------------------------------------ !
  ! ! --- [1] Pre-Process                                --- !
  ! ! ------------------------------------------------------ !
  ! !  -- [1-1] Initialize clock                         --  !
  ! call TimeMeasure( 0 )
  
  ! !  -- [1-2] Load Settings                            --  !
  ! call load__parameterFile
  ! call load__BFieldFile
  ! call load__mShapeFile
  ! ! call load__weightFile
  ! call TimeMeasure( 1 )
  
  ! !  -- [1-3] Generate Matrix                          --  !
  ! call allocate__variables
  ! call update__Matrix
  ! call update__BField
  ! call update__RHS
  ! call TimeMeasure( 3 )
  
  ! !  -- [1-4] output initial state    --  !
  ! call save__result
  ! if ( solverType.eq.'WLS' ) then
  !    call save__residual_WLS
  ! else
  !    call save__residual
  ! endif
  ! call TimeMeasure( 4 )
  
  ! ! ------------------------------------------------------ !
  ! ! --- [2] Main Loop                                  --- !
  ! ! ------------------------------------------------------ !
  ! do iter=1, iterMax
  !    ! -- [2-1] display loop info.     -- !
  !    write(6,*)
  !    write(6,'(2(a,i10))'  ) '[main] iter / iterMax = ', iter, ' / ', iterMax
  !    write(6,'(2(a,f10.5))') '[main] wPicard        = ', wPicard
     
  !    ! -- [2-2] Inverse Problem Solver -- !
  !    ! if ( solverType.eq.'SVD' ) call SVD_MatrixInvertor( Amat, rhs, shim, spectrum, singular, nMpt, nNpt, threshold )
  !    ! if ( solverType.eq.'LSM' ) call LSM_MatrixInvertor( Amat, rhs, shim,                     nMpt, nNpt, LSM_engine    )
  !    if ( solverType.eq.'WLS' ) call WLS_MatrixInvertor( Amat, rhs, shim, weights ,           nMpt, nNpt, LSM_engine    )
  !    ! if ( Flag__debugMode     ) call check__equations
  !    call TimeMeasure( 2 )
     
  !    ! -- [2-3] update pole & R.H.S.   -- !
  !    ! call update__MagnetShape
  !    ! call update__Matrix
  !    ! call update__BField
  !    ! call update__RHS
  !    call TimeMeasure( 3 )
     
  !    ! -- [2-4] write out Data         -- !
  !    if ( mod(iter,iter__save      ).eq.0 ) call save__result
  !    if ( solverType.eq.'WLS' ) then
  !       call save__residual_WLS
  !    else
  !       call save__residual
  !    endif
  !    if ( ( Flag__residConvergence ).and.( iter.gt.3 ) ) then
  !       call check__residConvergence
  !    end if
  !    call TimeMeasure( 4 )
  !    if ( mod( iter,iter__tMeasure ).eq.0 ) call TimeMeasure( -1 )

  !    ! -- [2-5] exit flag check        -- !
  !    if ( flag__exitStatus ) then
  !       write(6,*) "[main] Enough conversion was obtained.... EXIT main loop. "
  !       exit
  !    endif
     
  ! end do

  ! ! ------------------------------------------------------ !
  ! ! --- [3] post process                               --- !
  ! ! ------------------------------------------------------ !
  ! call deallocate__variables
  ! call TimeMeasure(  5 )
  ! call TimeMeasure( -1 )

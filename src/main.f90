program main
  use variablesMod
  use utilitiesMod
  use ioUtilityMod
  use allocatorMod
  use initiatorMod
  use invSolverMod
  use updateEqsMod
  implicit none
  
  ! ------------------------------------------------------ !
  ! --- [1] load parmeter / Data from File             --- !
  ! ------------------------------------------------------ !
  call TimeMeasure( 0 )
  call show__programLogo
  call print__section( "Load parameters", "=", 15, 4, 70 )
  call load__parameterFile
  
  call print__section( "Load Data File" , "=", 14, 4, 70 )
  call load__ElemsNodes
  call load__bfieldFile
  call load__mshapeFile
  call load__weightFile
  call TimeMeasure( 1 )

  ! ------------------------------------------------------ !
  ! --- [2] initialization                             --- !
  ! ------------------------------------------------------ !
  call print__section( "Initialize variables / Conditions", "=" , 33, 4, 70 )
  call allocate__variables
  call initialize__conditions
  call print__section( "Prepare Matrix & R.H.S.", "=" , 23, 4, 70 )
  call TimeMeasure( 1 )
  
  call update__matrix
  call update__bfield
  call update__rhs
  call TimeMeasure( 3 )

  call save__residuals
  call save__results
  call TimeMeasure( 4 )
  

  call print__section( "Begin Main Loop", "#" , 15, 4, 70 )

  ! ------------------------------------------------------ !
  ! --- [3] iteration of main solver                   --- !
  ! ------------------------------------------------------ !
  do iter=1, iterMax

     ! :::::::::::::::::::::::::::::::::::::::::::::::::::::: !
     ! ::: [3-1] display loop info                        ::: !
     ! :::::::::::::::::::::::::::::::::::::::::::::::::::::: !
     call print__loopTitle( iter, iterMax )
     call update__Picard

     ! :::::::::::::::::::::::::::::::::::::::::::::::::::::: !
     ! ::: [3-2] solve inverse problem                    ::: !
     ! :::::::::::::::::::::::::::::::::::::::::::::::::::::: !
     if ( solverType.eq."wls" ) then
        call WLS_MatrixInvertor( Amat, rhs, hvec, wvec, nMpt, nNpt, lsm__engine )
     endif
     call TimeMeasure( 2 )
     
     ! :::::::::::::::::::::::::::::::::::::::::::::::::::::: !
     ! ::: [3-3] update mshape & R.H.S.                   ::: !
     ! :::::::::::::::::::::::::::::::::::::::::::::::::::::: !
     call update__mshape
     call update__matrix
     call update__bfield
     call update__rhs
     call TimeMeasure( 3 )

     ! :::::::::::::::::::::::::::::::::::::::::::::::::::::: !
     ! ::: [3-4] save data in files                       ::: !
     ! :::::::::::::::::::::::::::::::::::::::::::::::::::::: !
     call save__residuals
     call check__residConvergence
     if ( mod( iter, iter__saveResults  ).eq.0 ) call save__results
     if ( flag__exitStatus ) then
        write(6,*) "[main]  Enough Convergence was attained .... EXIT main loop."
        exit
     endif
     call TimeMeasure( 4 )
     if ( mod( iter, iter__timeDisplay  ).eq.0 ) call TimeMeasure( -1 )

  enddo
  
  call print__section( "End of Main Loop", "#" , 16, 4, 70 )
  
  ! ------------------------------------------------------ !
  ! --- [4] post-process                               --- !
  ! ------------------------------------------------------ !
  call deallocate__variables
  call TimeMeasure(  5 )
  call TimeMeasure( -1 )
  call show__endLogo
  
end program main



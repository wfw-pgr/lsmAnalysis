module debugToolMod
contains

  ! ===================================================== !
  ! ===  check equations to be solved                 === !
  ! ===================================================== !
  subroutine check__equations
    use variablesMod
    implicit none
    integer :: i, j
    
    ! ------------------------------------- !
    ! --- [1] Amat check                --- !
    ! ------------------------------------- !
    write(6,*) '[check__equations] Amat check ::  '
    open(lun,file=trim(AmatFile),status='replace',form='formatted')
    write(lun,*) '# Aij'
    write(lun,*) '# ', nMpt, nNpt
    do i=1, nMpt
       write(lun,*) ( Amat(i,j), j=1, nNpt )
    enddo
    close(lun)
    ! ------------------------------------- !
    ! --- [2] R.H.S. check              --- !
    ! ------------------------------------- !
    write(6,*) '[check__equations] rhs check ::  '
    open(lun,file=trim(rhsbFile),status='replace',form='formatted')
    write(lun,*) '# rhs'
    write(lun,*) '# ', nMpt
    do i=1, nMpt
       write(lun,*) rhs(i)
    enddo
    close(lun)
    ! ------------------------------------- !
    ! --- [3] shim check                --- !
    ! ------------------------------------- !
    write(6,*) '[check__equations] shim check ::  '
    open(lun,file=trim(shimFile),status='replace',form='formatted')
    write(lun,*) '# shim'
    write(lun,*) '# ', nNpt
    do i=1, nNpt
       write(lun,*) shim(i)
    enddo
    close(lun)

    return
  end subroutine check__equations

end module debugToolMod

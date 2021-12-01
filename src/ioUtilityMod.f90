module ioUtilityMod
contains

  ! ====================================================== !
  ! === load elems and nodes                           === !
  ! ====================================================== !
  subroutine load__elemsNodes
    implicit none

    ! ------------------------------------------------------ !
    ! --- [1] load elements                              --- !
    ! ------------------------------------------------------ !
    write(6,"(a)",advance="no") "[load__elemsNodes] loading dat/elems.dat "
    open (lun,file=trim(mshpFile),form="formatted",status="old")
    read(lun,*) cmt
    read(lun,*) cmt, nElem, nComp
    read(lun,*) cmt, nElem, nComp
    allocate( elems(nComp,nElem) )
    elems(:,:) = 0
    do ik=1, nElem
       read(lun,*) elems(:,:)
    enddo
    close(lun)
    write(6,*) "      [Done] "

    ! ------------------------------------------------------ !
    ! --- [2] load nodes                                 --- !
    ! ------------------------------------------------------ !
    write(6,"(a)",advance="no") "[load__elemsNodes] loading dat/nodes.dat "
    open (lun,file=trim(mshpFile),form="formatted",status="old")
    read(lun,*) cmt
    read(lun,*) cmt, nNode, nComp
    read(lun,*) cmt, nNode, nComp
    allocate( nodes(nComp,nNode) )
    nodes(:,:) = 0.d0
    do ik=1, nNode
       read(lun,*) nodes(:,:)
    enddo
    close(lun)
    write(6,*) "      [Done] "
    
    return
  end subroutine load__elemsNodes

  
end module ioUtilityMod

module smoothMatMod
contains

  ! ====================================================== !
  ! === matrix :: simple average matrix                === !
  ! ====================================================== !
  subroutine generate__simpleAverageMatrix( elems, Cmat_e2n2e, nElems, nNodes )
    implicit none
    integer         , intent(in)  :: nElems, nNodes
    integer         , intent(in)  :: elems(3,nElems)
    double precision, intent(out) :: Cmat_e2n2e(nElems,nElems)
    integer                       :: ik
    double precision, allocatable :: Cmat_e2n(:,:), Cmat_n2e(:,:)

    ! ------------------------------------------------------ !
    ! --- [1] generate 2-types of conversion matrix      --- !
    ! ------------------------------------------------------ !
    allocate( Cmat_e2n(nNodes,nElems), Cmat_n2e(nElems,nNodes) )
    call generate__matrix_helem2hnode( elems, Cmat_e2n, nElems, nNodes )
    call generate__matrix_hnode2helem( elems, Cmat_n2e, nElems, nNodes )
    
    ! ------------------------------------------------------ !
    ! --- [2] minus unit matrix :: ( MB-[E] )            --- !
    ! ------------------------------------------------------ !
    call matrix__multiply( Cmat_e2n, Cmat_n2e, Cmat_e2n2e, nElems, nNodes, nElems )
    do ik=1, nElems 
       Cmat_e2n2e(ik,ik) = Cmat_e2n2e(ik,ik) - 1.d0
    enddo
    
    return
  end subroutine generate__simpleAverageMatrix
    
  
  ! ====================================================== !
  ! === matrix :: element height into node height      === !
  ! ====================================================== !
  subroutine generate__matrix_helem2hnode( elems, Cmat, nElems, nNodes )
    implicit none
    integer         , intent(in)  :: nElems, nNodes
    integer         , intent(in)  :: elems(3,nElems)
    double precision, intent(out) :: Cmat(nNodes,nElems)
    integer                       :: iE, iv, iN
    double precision              :: nNeib

    ! ------------------------------------------------------ !
    ! --- [1] geenrate matrix                            --- !
    ! ------------------------------------------------------ !
    do iE=1, nElems
       do iv=1, 3
          iN             = elems( iv, iE )
          Cmat( iN, iE ) = 1.d0
       enddo
    enddo

    ! ------------------------------------------------------ !
    ! --- [2] Normalization                              --- !
    ! ------------------------------------------------------ !
    do iN=1, nNodes
       nNeib      = sum( Cmat(iN,:) )
       Cmat(iN,:) = Cmat(iN,:) / nNeib
    enddo
    
    return
  end subroutine generate__matrix_helem2hnode


  ! ====================================================== !
  ! === matrix :: node height into element height      === !
  ! ====================================================== !
  subroutine generate__matrix_hnode2helem( elems, Cmat, nElems, nNodes )
    implicit none
    integer         , intent(in)  :: nElems, nNodes
    integer         , intent(in)  :: elems(3,nElems)
    double precision, intent(out) :: Cmat(nElems,nNodes)
    integer                       :: iE, iv, iN
    double precision              :: nNeib
    double precision, parameter   :: onethird = 1.d0 / 3.d0

    ! ------------------------------------------------------ !
    ! --- [1] geenrate matrix                            --- !
    ! ------------------------------------------------------ !
    do iE=1, nElems
       do iv=1, 3
          iN             = elems( iv, iE )
          Cmat( iE, iN ) = onethird
       enddo
    enddo

    ! ------------------------------------------------------ !
    ! --- [2] compress by grouping                       --- !
    ! ------------------------------------------------------ !
    ! do iN=1, nNodes
    ! enddo
    
    return
  end subroutine generate__matrix_hnode2helem


  ! ====================================================== !
  ! === matrix multiply subroutine                     === !
  ! ====================================================== !
  subroutine matrix__multiply( Amat, Bmat, Cmat, LI, LJ, LK )
    implicit none
    integer         , intent(in)  :: LI, LJ, LK
    double precision, intent(in)  :: Amat(LI,LJ), Bmat(LJ,LK)
    double precision, intent(out) :: Cmat(LI,LK)

    ! ------------------------------------------------------ !
    ! --- [1] gfortran standard library                  --- !
    ! ------------------------------------------------------ !
    Cmat = matmul( Amat, Bmat )
    
    return
  end subroutine matrix__multiply
    
  
end module smoothMatMod

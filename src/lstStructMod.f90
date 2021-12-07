module lstStructMod

  ! -- how to use  -- !
  ! type(list_type), pointer :: list
  ! nullify( list )   
  ! call add__elementInList( list, elementNum=10, groupNum=1 )
  ! call add__elementInList( list, elementNum=20, groupNum=2 )
  ! call show__nodeInList  ( list )
  ! call investigate__listInfo( list, max_nCell, list_length )
  ! allocate( cells(max_nCell), groupNums(list_length) )
  ! call obtain__cellsInGroup ( list, 3, nCell, cells, max_nCell )
  ! call obtain__groupNumArray( list, groupNum_array, list_length )
  ! ----------------- !
  
  ! ------------------------------------------------------ !
  ! ---  list element                                  --- !
  ! ------------------------------------------------------ !
  type :: list_type
     type(list_type), pointer     :: next
     integer                      :: elementNum
     integer                      :: groupNum
     integer                      :: nCell = 0
     integer        , allocatable :: cells(:)
  end type list_type
  
contains

  
  ! ====================================================== !
  ! === make new node of a list                        === !
  ! ====================================================== !
  subroutine make__nodeInList( list, next, elementNum, groupNum )
    implicit none
    type(list_type), pointer,           intent(out) :: list
    type(list_type), pointer, optional, intent(in ) :: next
    integer        ,          optional, intent(in ) :: elementNum, groupNum

    ! ------------------------------------------------------ !
    ! --- [1] list allocation check                      --- !
    ! ------------------------------------------------------ !
    if ( associated(list) ) then
       write(6,*) "[make__nodeInList] list is allocated.... [ERROR] "
       stop
    endif

    ! ------------------------------------------------------ !
    ! --- [2] allocation of a node                       --- !
    ! ------------------------------------------------------ !
    allocate( list      )
    nullify ( list%next )

    ! ------------------------------------------------------ !
    ! --- [3] store next or elementNum                        --- !
    ! ------------------------------------------------------ !
    if ( present( next  ) ) then
       list%next => next
    endif
    if ( present( groupNum ) ) then
       list%groupNum   = groupNum
    endif
    if ( present( elementNum ) ) then
       list%nCell      = list%nCell + 1
       allocate( list%cells(list%nCell) )
       list%cells(1)   = elementNum
    endif
    
    return
  end subroutine make__nodeInList


  ! ====================================================== !
  ! === delete the list                                === !
  ! ====================================================== !
  subroutine delete__list( list )
    implicit none
    type(list_type), pointer, intent(inout) :: list
    type(list_type), pointer                :: iter, temp

    if ( .not.( associated(list) ) ) then
       return
    endif

    iter => list
    do while( associated( iter%next ) )
       temp => iter
       iter => iter%next
       deallocate( temp )
    enddo

    nullify( list )
    return
  end subroutine delete__list


  ! ====================================================== !
  ! === return the pointer of a list                   === !
  ! ====================================================== !
  function tail__nodeInList(list) result(iter)
    type(list_type), pointer, intent(in) :: list
    type(list_type), pointer             :: iter
    
    iter => list
    do while( associated(iter%next) )
       iter => iter%next
    end do
    return
  end function tail__nodeInList


  ! ====================================================== !
  ! ===  append node to a list                         === !
  ! ====================================================== !
  subroutine append__nodeInList( list, elementNum, groupNum )
    implicit none
    integer        ,          intent(in) :: elementNum, groupNum
    type(list_type), pointer, intent(in) :: list
    type(list_type), pointer             :: iter

    iter => tail__nodeInList( list )
    call make__nodeInList( iter%next, elementNum=elementNum, groupNum=groupNum )
    return
  end subroutine append__nodeInList


  ! ====================================================== !
  ! === add element into list                          === !
  ! ====================================================== !
  subroutine add__elementInList( list, elementNum, groupNum )
    implicit none
    integer        ,          intent(in)    :: elementNum, groupNum
    type(list_type), pointer, intent(inout) :: list
    type(list_type), pointer                :: iter
    integer        , allocatable            :: temp(:)
    integer                                 :: ik
    logical                                 :: flag__groupExist

    if ( associated( list ) ) then
       ! -- list exist      -- !
       iter => list
       flag__groupExist = .false.
       ! -- search groupNum -- !
       do
          if ( iter%groupNum.eq.groupNum ) then
             ! -- copy cells                 -- !
             allocate( temp(iter%nCell) )
             temp(:) = iter%cells(:)
             ! -- resize cells -> nCell+1    -- !
             deallocate( iter%cells )
             iter%nCell = iter%nCell + 1
             allocate( iter%cells(iter%nCell) )
             ! -- append elementNum in cells -- !
             do ik=1, iter%nCell-1
                iter%cells(ik) = temp(ik)
             enddo
             iter%cells(iter%nCell) = elementNum
             flag__groupExist  = .true.
             exit
          endif
          if ( associated(iter%next) ) then
             iter => iter%next
          else
             exit
          endif
       enddo
       if ( .not.( flag__groupExist ) ) then
          ! -- no node for groupNum-- !
          call append__nodeInList( list, elementNum=elementNum, groupNum=groupNum )
       endif
    else
       ! -- first group -- !
       call make__nodeInList( list, elementNum=elementNum, groupNum=groupNum )
    endif

    return
  end subroutine add__elementInList


  ! ====================================================== !
  ! === investigate max. nCell / list length           === !
  ! ====================================================== !
  subroutine investigate__listInfo( list, max_nCell, list_length )
    implicit none
    type(list_type), pointer, intent(in)  :: list
    integer        ,          intent(out) :: max_nCell, list_length
    type(list_type), pointer              :: iter
    
    max_nCell   = 0
    list_length = 0
    if ( associated( list ) ) then
       iter => list
       do
          list_length = list_length + 1
          max_nCell   = max( max_nCell, iter%nCell )
          if ( associated(iter%next) ) then
             iter => iter%next
          else
             exit
          endif
       enddo
    else
       write(6,*) "[investigate__listInfo] list is empty.... "
       return
    endif
    
    return
  end subroutine investigate__listInfo
  
  
  ! ====================================================== !
  ! === obtain cells in a group                        === !
  ! ====================================================== !
  subroutine obtain__cellsInGroup( list, groupNum, nCell, cell_ret, cell_len )
    implicit none
    type(list_type), pointer, intent(in)  :: list
    integer                 , intent(in)  :: groupNum
    integer                 , intent(in)  :: cell_len
    integer                 , intent(out) :: nCell, cell_ret(cell_len)
    type(list_type), pointer              :: iter
    logical                               :: flag__groupExist
    
    flag__groupExist     = .false.
    cell_ret(1:cell_len) = 0
    if ( associated( list ) ) then
       iter => list
       do
          if ( iter%groupNum.eq.groupNum ) then
             cell_ret(1:iter%nCell) = iter%cells(1:iter%nCell)
             nCell                  = iter%nCell
             flag__groupExist       = .true.
             exit
          end if
          if ( associated(iter%next) ) then
             iter => iter%next
          else
             exit
          endif
       enddo
       if ( .not.( flag__groupExist ) ) then
          write(6,*) "[obtain__cellsInGroup] cannot find groupNum in list [ERROR] "
          return
       endif
    else
       write(6,*) "[obtain__cellsInGroup] list is empty.... "
       return
    endif
    
    return
  end subroutine obtain__cellsInGroup


  ! ====================================================== !
  ! === obtain groupNum array                          === !
  ! ====================================================== !
  subroutine obtain__groupNumArray( list, groupNum_array, array_len )
    implicit none
    type(list_type), pointer, intent(in)  :: list
    integer                 , intent(in)  :: array_len
    integer                 , intent(out) :: groupNum_array(array_len)
    type(list_type), pointer              :: iter
    integer                               :: count
    
    if ( associated( list ) ) then
       iter => list
       count = 0
       do
          count                 = count + 1
          groupNum_array(count) = iter%groupNum
          
          if ( associated(iter%next) ) then
             iter => iter%next
          else
             exit
          endif
       enddo
    else
       write(6,*) "[obtain__groupNumArray] list is empty.... "
       return
    endif
    return
  end subroutine obtain__groupNumArray

  
  ! ====================================================== !
  ! === show the contents of the list                  === !
  ! ====================================================== !
  subroutine show__nodeInList( list )
    implicit none
    type(list_type), pointer, intent(in) :: list
    type(list_type), pointer             :: iter
    
    if ( associated( list ) ) then
       iter => list
       do
          write(6,*) "groupNum = "  , iter%groupNum, "nCell = ", iter%nCell, &
               &     "elementNum = ", iter%cells(:)
          
          if ( associated(iter%next) ) then
             iter => iter%next
          else
             exit
          endif
       enddo
    else
       write(6,*) "[show__nodeInList] list is empty.... "
       return
    endif
    return
  end subroutine show__nodeInList

  
end module lstStructMod


  ! ! ====================================================== !
  ! ! === insert of a node at the index                  === !
  ! ! ====================================================== !
  ! subroutine insert__nodeInList( list, index, elementNum )
  !   implicit none
  !   type(list_type), pointer, intent(inout) :: list
  !   integer                 , intent(in)    :: index
  !   integer                 , intent(in)    :: elementNum
  !   integer                                 :: ik
  !   type(list_type), pointer                :: iter, prev, node

  !   ! ------------------------------------------------------ !
  !   ! --- [1] index starts from 1                        --- !
  !   ! ------------------------------------------------------ !
  !   if ( index < 1 ) then
  !      write(6,*) "[insert__nodeInList] index ERROR. no such index. [ERROR]"
  !      write(6,*) "[insert__nodeInList] index == ", index
  !      stop
  !   endif

  !   ! ------------------------------------------------------ !
  !   ! --- [2] skip until index                           --- !
  !   ! ------------------------------------------------------ !
  !   nullify( prev )
  !   iter => list
  !   do ik=1, index-1
  !      if ( .not.( associated( iter%next ) ) ) then
  !         write(6,*) "[insert__nodeInList] index ERROR. no such index. [ERROR]"
  !         write(6,*) "[insert__nodeInList] index == ", index
  !         stop
  !      endif
  !      prev => iter
  !      iter => iter%next
  !   enddo

  !   ! ------------------------------------------------------ !
  !   ! --- [3] create a new node                          --- !
  !   ! ------------------------------------------------------ !
  !   nullify( node )
  !   call make__nodeInList( node, next=iter, elementNum=elementNum )
  !   if ( .not.associated( prev ) ) then ! -- this means first node -- !
  !      list => node
  !   else
  !      prev%next => node
  !   endif
    
  ! end subroutine insert__nodeInList



  ! ! ====================================================== !
  ! ! === remove a node at a given index                 === !
  ! ! ====================================================== !
  ! subroutine remove__nodeInList( list, index )
  !   implicit none
  !   type(list_type), pointer, intent(inout) :: list
  !   integer        ,          intent(in)    :: index
  !   integer                                 :: ik
  !   type(list_type), pointer                :: iter, prev

  !   ! ------------------------------------------------------ !
  !   ! --- [1] index starts from 1                        --- !
  !   ! ------------------------------------------------------ !
  !   if ( index < 1 ) then
  !      write(6,*) "[remove__nodeInList] index ERROR. no such index. [ERROR]"
  !      write(6,*) "[remove__nodeInList] index == ", index
  !      stop
  !   endif

  !   ! ------------------------------------------------------ !
  !   ! --- [2] skip until index                           --- !
  !   ! ------------------------------------------------------ !
  !   nullify( prev )
  !   iter => list
  !   do ik=1, index-1
  !      if ( .not.( associated( iter%next ) ) ) then
  !         write(6,*) "[remove__nodeInList] index ERROR. no such index. [ERROR]"
  !         write(6,*) "[remove__nodeInList] index == ", index
  !         stop
  !      endif
  !      prev => iter
  !      iter => iter%next
  !   enddo

  !   ! ------------------------------------------------------ !
  !   ! --- [3] remove a node to be removed                --- !
  !   ! ------------------------------------------------------ !
  !   if ( .not.associated( prev ) ) then ! -- this means first node -- !
  !      list => iter%next
  !   else
  !      prev%next => iter%next
  !   endif
  !   deallocate(iter)

  ! end subroutine remove__nodeInList
  

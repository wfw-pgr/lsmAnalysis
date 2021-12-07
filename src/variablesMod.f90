module variablesMod
  implicit none
  
  ! ------------------------------------------------------ !
  ! --- [1] Constants                                  --- !
  ! ------------------------------------------------------ !
  integer         , parameter   :: cLen              = 300
  integer         , parameter   :: lun               =  50
  integer         , parameter   :: dim               =   3
  integer         , parameter   :: nVert             =   3
  integer         , parameter   :: nMaxPicard        =  10
  integer         , parameter   :: iter__timeDisplay =   5
  integer         , parameter   :: iter__saveResults =   5

  double precision, parameter   :: unit__thickness   = 1.d-3

  ! ------------------------------------------------------ !
  ! --- [2] suffix designation                         --- !
  ! ------------------------------------------------------ !
  !  -- for bfield, mshape, zp, weights, respectively  --  !
  integer         , parameter   :: xb_=1, yb_=2, zb_=3, bi_=4, bb_=5, bs_=6, bt_=7, be_=8
  integer         , parameter   :: xm_=1, ym_=2, zm_=3, mi_=4, ms_=5, mf_=6, mg_=7
  integer         , parameter   :: lo_=1, hi_=2
  integer         , parameter   :: xw_=1, yw_=2, zw_=3, cl_=4, wt_=5
  
  ! ------------------------------------------------------ !
  ! --- [3] I/O Files                                  --- !
  ! ------------------------------------------------------ !
  character(cLen) , parameter   :: listFile          = "dat/parameter.lst"
  character(cLen) , parameter   :: nodeFile          = "dat/nodes.dat"
  character(cLen) , parameter   :: elemFile          = "dat/elems.dat"
  character(cLen) , parameter   :: binpFile          = "dat/bfield_input.dat"
  character(cLen) , parameter   :: mshpFile          = "dat/mshape_lsm.dat"
  character(cLen) , parameter   :: wghtFile          = "dat/weights.dat"
  character(cLen) , parameter   :: bresFile          = "out/bfield_resid.dat"

  ! ------------------------------------------------------ !
  ! --- [4] flags                                      --- !
  ! ------------------------------------------------------ !
  logical                       :: flag__exitStatus             = .false.
  logical                       :: flag__laplaceRegularization  = .false.
  

  ! ------------------------------------------------------ !
  ! --- [5] parameters to be loaded from list File     --- !
  ! ------------------------------------------------------ !
  integer                       :: nSubdiv, nDiv_z, iterMax
  double precision              :: MzConst, zLim1, zLim2, resid__criterion
  character(4)                  :: lsm__engine       = "cgls"
  character(3)                  :: solverType        = "wls"
  character(10)                 :: convergenceType   = "rmse"

  
  ! ------------------------------------------------------ !
  ! --- [6] system variables                           --- !
  ! ------------------------------------------------------ !
  integer                       :: iter
  integer                       :: nBpt, nMpt, nNpt, nEpt
  integer                       :: nElems, nNodes, nColor
  double precision              :: coefPicard, wPicard(nMaxPicard), iPicard(nMaxPicard)
  
  ! ------------------------------------------------------ !
  ! --- [7] data variables                             --- !
  ! ------------------------------------------------------ !
  double precision              :: mvec(3)
  integer         , allocatable :: elems(:,:)
  double precision, allocatable :: nodes(:,:) , vertex(:,:,:)
  double precision, allocatable :: bfield(:,:), mshape(:,:)
  double precision, allocatable :: Amat(:,:)  , Rmat(:,:)
  double precision, allocatable :: xvec(:)    , hvec(:)
  double precision, allocatable :: wvec(:)    , weights(:,:)
  double precision, allocatable :: rhs (:)

  ! ------------------------------------------------------ !
  ! --- [8] residual statistics                        --- !
  ! ------------------------------------------------------ !
  double precision, allocatable :: avgs  (:), stds  (:), rmse   (:)
  double precision, allocatable :: sumCnt(:), sumErr(:), sumErr2(:)
  double precision, allocatable :: minErr(:), maxErr(:)
  double precision              :: minRes(3), maxRes(3), rmsRes(3)

end module variablesMod

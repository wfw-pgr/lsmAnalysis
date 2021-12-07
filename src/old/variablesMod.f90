module variablesMod
  implicit none
  ! ------------------------------------- !
  ! --- [1]  Constants                --- !
  ! ------------------------------------- !
  double precision, parameter   :: unit_thick      = 1.d-3 ! ( 1 mm : arbitral )
  logical                       :: Flag__Laplacian = .true.
  logical                       :: Flag__debugMode = .false.
  logical                       :: Flag__residConvergence = .true.
  logical                       :: Flag__initShape
  logical                       :: Flag__fullModel
  integer         , parameter   :: iter__tMeasure  =  10
  integer         , parameter   :: iter__save      =   1
  integer         , parameter   :: cLen            = 300
  integer         , parameter   :: lun             =  50
  ! ------------------------------------- !
  ! --- [2]  I/O Files                --- !
  ! ------------------------------------- !
  ! --   input   -- !
  character(cLen) , parameter   :: parmFile   = 'dat/parameter.conf'
  character(cLen) , parameter   :: listFile   = 'dat/parameter.lst'
  character(cLen) , parameter   :: BinpFile   = 'dat/bfield_input.dat'
  character(cLen) , parameter   :: mshpFile   = 'dat/mshape_input.dat'
  character(cLen) , parameter   :: wghtFile   = 'dat/weights.dat'
  character(cLen) , parameter   :: wcnfFile   = 'dat/weights.conf'
  ! --  output   -- !
  character(cLen) , parameter   :: BresFile   = 'out/bfield_resid.dat'
  ! -- for debug -- !
  character(cLen) , parameter   :: AmatFile   = 'chk/Amat.dat'
  character(cLen) , parameter   :: rhsbFile   = 'chk/rhsb.dat'
  character(cLen) , parameter   :: shimFile   = 'chk/shim.dat'
  ! ------------------------------------- !
  ! --- [3] Variables                 --- !
  ! ------------------------------------- !
  integer                       :: iter = 0
  integer                       :: iterMax, LI, LJ
  integer                       :: nBpt, nEpt , nNpt, nMpt, LM, LN, nDiv_B, iPicard_iter, nBpt_pole, nBpt_peel
  integer         , allocatable :: elems(:,:) , idtable(:,:)
  double precision, allocatable :: nodes(:,:) ,   initS(:,:), mshape(:,:,:)
  double precision, allocatable :: MagPos(:,:),  BField(:,:)
  double precision, allocatable :: shim(:), rhs(:), spectrum(:), singular(:), surf(:), lapl(:)
  double precision, allocatable :: Rmat(:,:), Lmat(:,:), Amat(:,:)
  double precision, allocatable :: weights(:), weight_info(:,:)
  double precision, allocatable :: avgs(:), stds(:), sumCnt(:), sumErr(:), sumErr2(:), minErr(:), maxErr(:)
  double precision              :: minRes(3), maxRes(3)
  double precision              :: wLaplace, MzConstant, mvec(3), wPicard, threshold
  double precision              :: wPicard_init, wPicard_last
  double precision              :: rLim1, rLim2, pLim1, pLim2, zLim1, zLim2
  double precision              :: xMin, xMax, yMin, yMax, dx, dy, ddx, ddy
  double precision              :: pJump
  double precision              :: resid_convergence
  character(3)                  :: solverType
  character(4)                  :: modelType
  character(4)                  :: LSM_engine
  
  logical                       :: flag__exitStatus = .false.

end module variablesMod

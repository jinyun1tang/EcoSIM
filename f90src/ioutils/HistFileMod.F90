module HistFileMod

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use ncdio_pio
  use fileUtil, only : iulog
implicit none
  save
  private
  character(len=*),  parameter :: mod_filename =__FILE__
  integer, parameter :: max_chars = 256
  integer :: nfmaster = 0                        ! number of fields in master field list
  integer , public, parameter :: max_tapes = 6          ! max number of history tapes
  integer , public, parameter :: max_flds = 2500        ! max number of history fields
  integer , public, parameter :: max_namlen = 64        ! maximum number of characters for field name
  integer , private, parameter :: hist_dim_name_length = 16 ! lenngth of character strings in dimension names

  type field_info
     character(len=max_namlen) :: name         ! field name
     character(len=max_chars)  :: long_name    ! long name
     character(len=max_chars)  :: standard_name  ! CF standard name
     character(len=max_chars)  :: units        ! units
     character(len=hist_dim_name_length) :: type1d                ! pointer to first dimension type from data type (nameg, etc)
     character(len=hist_dim_name_length) :: type1d_out            ! hbuf first dimension type from data type (nameg, etc)
     character(len=hist_dim_name_length) :: type2d                ! hbuf second dimension type ["levgrnd","levlak","numrad","ltype","natpft","cft","glc_nec","elevclas","subname(n)","month"]
     integer :: beg1d                          ! on-node 1d clm pointer start index
     integer :: end1d                          ! on-node 1d clm pointer end index
     integer :: num1d                          ! size of clm pointer first dimension (all nodes)
     integer :: numdims                        ! the actual number of dimensions, this allows
                                               ! for 2D arrays, where the second dimension is allowed
                                               ! to be 1
     integer :: beg1d_out                      ! on-node 1d hbuf pointer start index
     integer :: end1d_out                      ! on-node 1d hbuf pointer end index
     integer :: num1d_out                      ! size of hbuf first dimension (all nodes)
     integer :: num2d                          ! size of hbuf second dimension (e.g. number of vertical levels)
     integer :: hpindex                        ! history pointer index
!     character(len=8) :: p2c_scale_type        ! scale factor when averaging pft to column
!     character(len=8) :: c2l_scale_type        ! scale factor when averaging column to landunit
!     character(len=8) :: l2g_scale_type        ! scale factor when averaging landunit to gridcell
!     character(len=8) :: t2g_scale_type        ! scale factor when averaging topounit to gridcell
     integer :: no_snow_behavior               ! for multi-layer snow fields, flag saying how to treat times when a given snow layer is absent
  end type field_info


  type history_entry
     type (field_info) :: field                ! field information
     character(len=1)  :: avgflag              ! time averaging flag
     real(r8), pointer :: hbuf(:,:)            ! history buffer (dimensions: dim1d x num2d)
     integer , pointer :: nacs(:,:)            ! accumulation counter (dimensions: dim1d x num2d)
  end type history_entry

  type history_tape
     integer  :: nflds                         ! number of active fields on tape
     integer  :: ntimes                        ! current number of time samples on tape
     integer  :: mfilt                         ! maximum number of time samples per tape
     integer  :: nhtfrq                        ! number of time samples per tape
     integer  :: ncprec                        ! netcdf output precision
!     logical  :: dov2xy                        ! true => do xy average for all fields
     logical  :: is_endhist                    ! true => current time step is end of history interval
     real(r8) :: begtime                       ! time at beginning of history averaging interval
     type (history_entry) :: hlist(max_flds)   ! array of active history tape entries
  end type history_tape


  character(len=max_chars) :: locfnh(max_tapes)  ! local history file names
  ! Master list: an array of master_entry entities
  !
  type master_entry
     type (field_info)  :: field               ! field information
     logical            :: actflag(max_tapes)  ! active/inactive flag
     character(len=1)   :: avgflag(max_tapes)  ! time averaging flag ("X","A","M" or "I",)
  end type master_entry

  type (master_entry) :: masterlist(max_flds)  ! master field list


  character(len=16), parameter :: namec  = 'column'       ! name of columns
  character(len=16), parameter :: namep  = 'pft'          ! name of patches
  type (history_tape), public :: tape(max_tapes)       ! array history tapes

  character(len=max_chars) :: locfnhr(max_tapes) ! local history restart file names
  logical :: htapes_defined = .false.            ! flag indicates history contents have been defined

! ids for netcdf files
  integer :: time_dimid                      ! time dimension id
  integer :: hist_interval_dimid             ! time bounds dimension id
  integer :: strlen_dimid                    ! string dimension id
  type(file_desc_t) :: nfid(max_tapes)       ! file ids
  type(file_desc_t) :: ncid_hist(max_tapes)  ! file ids for history restart files


  public :: hist_addfld1d        ! Add a 1d single-level field to the master field list
  public :: hist_addfld2d        ! Add a 2d multi-level field to the master field list
  public :: hist_addfld3d        ! Add a 2d multi-level field to the master field list
  public :: hist_printflds       ! Print summary of master field list
  public :: hist_htapes_build    ! Initialize history file handler for initial or continue run
  public :: hist_update_hbuf     ! Updates history buffer for all fields and tapes
  public :: hist_htapes_wrapup   ! Write history tape(s)
  public :: hist_restart_ncd     ! Read/write history file restart data
  public :: htapes_fieldlist     ! Define the contents of each history file based on namelist

  contains


  subroutine htape_create(t,histrest)

  !
  !DESCRIPTION
  !create history file for output

  use EcoSIMConfig      , only : case_name
  use GridConsts        , only : JZ,JS,JBR,JC,JP,jpstgs
  use EcoSIMConfig      , only : jcplx=>jcplxc,jsken=>jskenc
  use ElmIDMod          , only : npelms
  implicit none
  integer, intent(in) :: t        ! tape index
  logical, optional, intent(in) :: histrest  !if creating the history restart file

  integer :: ncprec              ! output netCDF write precision
  logical :: lhistrest           ! local history restart flag
  type(file_desc_t) :: lnfid     ! local file id
  character(len=  8) :: curdate  ! current date
  character(len=  8) :: curtime  ! current time
  character(len=256) :: name     ! name of attribute
  character(len=256) :: units    ! units of attribute
  character(len=256) :: str      ! global attribute string
  character(len=  1) :: avgflag  ! time averaging flag
  integer :: numc
  integer :: nump
  integer :: dimid               ! dimension id temporary
  character(len=*),parameter :: subname = 'htape_create'

  if ( present(histrest) )then
    lhistrest = histrest
  else
    lhistrest = .false.
  end if

  call get_grid_info(nc=numc, np=nump)

  ncprec = tape(t)%ncprec

  if ( .not. lhistrest )then
    call ncd_pio_createfile(lnfid, trim(locfnh(t)))
    call ncd_putatt(lnfid, ncd_global, 'title', 'EcoSIM History file information' )
  else
    call ncd_pio_createfile(lnfid, trim(locfnhr(t)))
    call ncd_putatt(lnfid, ncd_global, 'title', &
       'EcoSIM Restart History information, required to continue a simulation' )
    call ncd_putatt(lnfid, ncd_global, 'comment', &
                   "This entire file NOT needed for startup or branch simulations")
  endif

!  call ncd_putatt(lnfid, ncd_global, 'source'  , trim(source))
!  call ncd_putatt(lnfid, ncd_global, 'source_id'  , trim(version))
!  call ncd_putatt(lnfid, ncd_global, 'product'  , 'model-output')
  call ncd_putatt(lnfid, ncd_global, 'case', trim(case_name))
!  call ncd_putatt(lnfid, ncd_global, 'username', trim(username))
!  call ncd_putatt(lnfid, ncd_global, 'hostname', trim(hostname))
!  call ncd_putatt(lnfid, ncd_global, 'git_version' , trim(version))
  call getdatetime(curdate, curtime)

  ! Global compressed dimensions (not including non-land points)
!  call ncd_defdim(lnfid, trim(nameg), numg, dimid)
!  call ncd_defdim(lnfid, trim(namet), numt, dimid)
!  call ncd_defdim(lnfid, trim(namel), numl, dimid)
  call ncd_defdim(lnfid, trim(namec), numc, dimid)
  call ncd_defdim(lnfid, trim(namep), nump, dimid)

  ! "level" dimensions
  call ncd_defdim(lnfid, 'levgsoi', JZ, dimid)
  call ncd_defdim(lnfid, 'levsno',  JS,dimid)
  call ncd_defdim(lnfid, 'levcanopy',JC,dimid)
  call ncd_defdim(lnfid, 'npfts',  JP,dimid)
  call ncd_defdim(lnfid, 'nbranches',JBR,dimid)
  call ncd_defdim(lnfid, 'ngrstages',jpstgs,dimid)
  call ncd_defdim(lnfid, 'elements',npelms,dimid)
  call ncd_defdim(lnfid, 'nkinecomp',jsken,dimid)
  call ncd_defdim(lnfid, 'nomcomplx',jcplx,dimid)

  if ( .not. lhistrest )then
    call ncd_defdim(lnfid, 'hist_interval', 2, hist_interval_dimid)
    call ncd_defdim(lnfid, 'time', ncd_unlimited, time_dimid)
    nfid(t) = lnfid
    write(iulog,*) trim(subname), &
                    ' : Successfully defined netcdf history file ',t
    flush(iulog)

  else
    ncid_hist(t) = lnfid

    write(iulog,*) trim(subname), &
                      ' : Successfully defined netcdf restart history file ',t
    flush(iulog)
  end if

  end subroutine htape_create



!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: getdatetime
!
! !INTERFACE:
  subroutine getdatetime (cdate, ctime)
!
! !DESCRIPTION:
! A generic Date and Time routine
!
! !USES:
! !ARGUMENTS:
  implicit none
  character(len=8), intent(out) :: cdate  !current date
  character(len=8), intent(out) :: ctime  !current time
! !LOCAL VARIABLES:
!EOP
  character(len=8)      :: date       !current date
  character(len=10)     :: time       !current time
  character(len=5)      :: zone       !zone
  integer, dimension(8) :: values !temporary
  integer               :: ier    !MPI error code


  call date_and_time (date, time, zone, values)

  cdate(1:2) = date(5:6)
  cdate(3:3) = '/'
  cdate(4:5) = date(7:8)
  cdate(6:6) = '/'
  cdate(7:8) = date(3:4)

  ctime(1:2) = time(1:2)
  ctime(3:3) = ':'
  ctime(4:5) = time(3:4)
  ctime(6:6) = ':'
  ctime(7:8) = time(5:6)

  end subroutine getdatetime
!-----------------------------------------------------------------------
  subroutine get_grid_info(nc,np)
  use GridConsts, only : bounds
  implicit none
  integer, intent(out) :: nc,np

  nc=bounds%ncols
  np=bounds%npfts
  end subroutine get_grid_info


  !-----------------------------------------------------------------------
  subroutine hist_printflds()
    !
    ! !DESCRIPTION:
    ! Print summary of master field list.
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
  implicit none
  integer :: nf
  character(len=*),parameter :: subname = 'ELM_hist_printflds'
    !-----------------------------------------------------------------------

!    if (masterproc) then
     write(iulog,*) trim(subname),' : number of master fields = ',nfmaster
     write(iulog,*)' ******* MASTER FIELD LIST *******'
     do nf = 1,nfmaster
       write(iulog,9000)nf, masterlist(nf)%field%name, masterlist(nf)%field%units
9000      format (i5,1x,a32,1x,a16)
     end do
     flush(iulog)
!    end if

  end subroutine hist_printflds
end module HistFileMod

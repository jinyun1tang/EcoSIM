module PerturbationMod

! Description
! Code to do ecosystem peturbation. 
! E.g. soil warming, atmospheric warming, snow exclusion
! 
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use MiniMathMod,      only: yearday,             isletter, isLeap
  use ClimForcDataType, only: TKS_ref_vr
  use GridConsts,       only: JY,                  JX
  use GridDataType,     only: CumDepz2LayBottom_vr, NU,       NL
  use EcoSimConst  
  implicit none
  private
  save
  character(len=*), parameter :: mod_filename = &
  __FILE__

  real(r8) :: warm_depz                !lower warming depth [m]
  real(r8) :: warm_dTK                 !warming magnitude, [K]
  integer  :: warm_yearb,warm_yeare    !beginning/end year of warming
  integer  :: warm_doyb,warm_doye      !beginning/end doy of warming for the first day and last day
  integer, target, allocatable  :: warm_LL(:,:)                  !vertical layer number of the warming perturbation
  logical  :: lsoil_warming  =.false.
  character(len=256) :: fname_warming_Tref

  public :: check_Soil_Warming
  public :: get_warming_fname
  public :: config_soil_warming
  public :: set_soil_warming
  public :: InitSoilWarming
  public :: destructSoilWarming
  public :: is_warming_layerL
  contains

  subroutine InitSoilWarming()

  implicit none

  if(lsoil_warming)then
    allocate(warm_LL(JY,JX)); warm_LL=0
  endif


  end subroutine InitSoilWarming
!------------------------------------------------------------------------------------------

  subroutine destructSoilWarming()
  use abortutils, only : destroy
  implicit none

  if(lsoil_warming)then
    call destroy(warm_LL)
  endif
  end subroutine destructSoilWarming

!------------------------------------------------------------------------------------------


  subroutine Extract_year_doy(datestr,year,doy)
  !
  !input format is 
  !YYYY/MM/DD
  implicit none
  character(len=*), intent(in) :: datestr
  integer, intent(out) :: year
  integer, intent(out) :: doy

  integer :: month, day 

  read(datestr(1:4),*)year
  read(datestr(6:7),*)month
  read(datestr(9:10),*)day

  doy=yearday(year,month,day)

  end subroutine Extract_year_doy

!------------------------------------------------------------------------------------------

  pure function check_Soil_Warming(year,I)result(ans)

  implicit none
  integer, intent(in) :: year,I

  logical :: ans

  if(.not.lsoil_warming)then
    ans=.false.
  else
    ans=(year==warm_yearb .and. I>=warm_doyb) .or. (year==warm_yeare .and. I<=warm_doye) &
      .or. (year>warm_yearb .and. year <warm_yeare)
  endif
  end function check_Soil_Warming
!------------------------------------------------------------------------------------------

  subroutine config_soil_warming(warming_exp)

  !!
  ! Description
  ! extract soil warming configuration
  ! warming_exp='4K;1m;Blodget.ctrl.ecosim.h1.xxxx-01-00-00000.nc;2014/01/01:2018/12/31'  

  implicit none
  character(len=*), intent(in) :: warming_exp

  character(len=120) :: str_loc
  integer :: k, k1, it, nl, nl1, jj

  k=1
  k1=1
  it=0
  nl=len(trim(warming_exp))    
  lsoil_warming=nl>20
  !return when soil warming is off
  if(.not.lsoil_warming)return

  do k=1,nl
    if (warming_exp(k:k)==';')then
      it=it+1

      !extract temperature
      if(it==1)then      
        read(warming_exp(k1:k-2),*)warm_dTK
        k1=k+1        
      !extracting depth   
      elseif(it==2)then      
        write(str_loc,*)warming_exp(k1:k-1)
        nl1=len(trim(str_loc))
        do jj = 1,nl1
          if(isletter(str_loc(jj:jj)))then
            read(str_loc(1:jj-1),*)warm_depz
            exit
          endif
        enddo  
        k1=k+1
      !extracting file name
      elseif(it==3)then      
        write(fname_warming_Tref,*)warming_exp(k1:k-1)
        k1=k+1
      endif
      !extracting time info
    elseif(warming_exp(k:k)==':' .and. it==3)then      
      write(str_loc,'(A)')warming_exp(k1:k-1)
      call Extract_year_doy(str_loc,warm_yearb,warm_doyb)      
      
      write(str_loc,'(A)')warming_exp(k+1:k+12)
      call Extract_year_doy(str_loc,warm_yeare,warm_doye)       
    endif
  enddo
  write(*,*)'warming ',warm_dTK, ' K to ',warm_depz, ' m depth from day ',warm_doyb, &
    ' in year ',warm_yearb, 'to day ',warm_doye, ' in year ',warm_yeare

  end subroutine config_soil_warming

!------------------------------------------------------------------------------------------
  subroutine get_warming_fname(year,fname)
  !
  !get warming file name for year 
  implicit none  
  integer, intent(in) :: year
  character(len=256), intent(out) :: fname
  integer :: jj

  write(fname,'(A)')trim(fname_warming_Tref)

  do jj = 1, 256
    if(fname(jj:jj+3)=='xxxx')then
      write(fname(jj:jj+3),'(I4)')year
      exit
    endif
  enddo
  end subroutine get_warming_fname

!------------------------------------------------------------------------------------------

  subroutine set_soil_warming(year,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: year
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NY,NX,L
  integer :: i1,i2

  if(.not.lsoil_warming)return
  if(year==warm_yearb)then    
    i1=(warm_doyb-1)*24+1
    i2=365
    if(isLeap(year))i2=366
    i2=i2*24
  elseif(year==warm_yeare)then    
    i1=1
    i2=warm_doye*24
  else
    i1=1
    i2=365
    if(isLeap(year))i2=366
    i2=i2*24      
  endif

  DO NX=NHW,NHE
    DO NY=NVN,NVS
      DO L=NU(NY,NX),NL(NY,NX)
        IF(CumDepz2LayBottom_vr(L,NY,NX).LE.warm_depz)THEN
          warm_LL(NY,NX)=L
          TKS_ref_vr(i1:i2,L,NY,NX)=TKS_ref_vr(i1:i2,L,NY,NX)+warm_dTK
        ENDIF    
      ENDDO  
    ENDDO
  ENDDO
  end subroutine set_soil_warming

!------------------------------------------------------------------------------------------

  function is_warming_layerL(L,NY,NX)result(ans)

  implicit none 
  integer, intent(in) :: L,NY,NX
  logical :: ans

  ans=warm_LL(NY,NX)>=L
  end function is_warming_layerL
end module PerturbationMod

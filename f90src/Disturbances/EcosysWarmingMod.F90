module EcosysWarmingMod
!
! Description:
! Code to do ecosystem warming. 
! E.g. soil warming, atmospheric warming, snow exclusion
! example warming configurations include
!  warming_exp = 'loc[NY,NX];type[Cable_heating];dT[4K];Depth[1m];hist_ctrl[Blodget.ctrl.ecosim.h1.xxxx-01-00-00000.nc];' &
!    //'beg_time[2014/01/01];end_time[2018/12/31]'

! 
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use MiniMathMod,      only: yearday,   isletter, isLeap,iisLeap
  use ClimForcDataType, only: TKS_ref_vr
  use GridConsts,       only: JY, JX
  use GridDataType,     only: CumDepz2LayBottom_vr, NU_col, NL_col
  use fileUtil,         only: iulog  
  use abortutils,       only : endrun  
  use DebugToolMod,     only : PrintInfo
  use StrToolsMod  
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
  integer  :: ihtime                   !0:all day, 1:day, 2:night
  integer  :: seas_doyb, seas_doye     !beginning and ending ordinal day in a year for warming
  real(r8) :: warm_IR_power            !infrared heating, MJ/m2/hr
  integer  :: igrowth                  !0 for whole year, 1 for growing season only
  real(r8) :: wind_otc                 !wind speed in otc, m/s
  real(r8) :: ems_reduct_otc           !emissivity reduction by otc,%
  real(r8) :: taus_otc                 !solar radiation transmissivity by OTC
  integer, target, allocatable  :: warm_type(:,:)         !warming type for each grid
  integer, target, allocatable  :: warm_LL(:,:)       !vertical layer number of the warming perturbation
  integer :: leap_cor(2),iwtype
  logical  :: lsoil_warming  =.false.
  character(len=256) :: fname_warming_Tref
  integer :: NY_w,NX_w

  public :: check_warming_dates
  public :: get_warming_fname
  public :: config_soil_warming
  public :: apply_soil_cable_warming
  public :: InitSoilWarming
  public :: destructSoilWarming
  public :: is_warming_layerL
  public :: apply_IR_warming
  public :: apply_OTC_warming
  contains

  subroutine InitSoilWarming()

  implicit none

  if(lsoil_warming)then
    allocate(warm_type(JY,JX));warm_type=-1
    allocate(warm_LL(JY,JX)); warm_LL=0
    warm_type(NY_w,NX_w) = iwtype
  endif


  end subroutine InitSoilWarming
!------------------------------------------------------------------------------------------

  subroutine destructSoilWarming()
  use abortutils, only : destroy
  implicit none

  if(lsoil_warming)then
    call destroy(warm_LL)
    call destroy(warm_type)
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
  function Extract_dates2doy(datestr,doy1,doy2)result(date_found)
  implicit none
  character(len=*),intent(in) :: datestr  
  integer, intent(out) :: doy1,doy2
  integer :: lens
  integer :: k1,mon,day,KK
  logical :: date_found

  lens=len_trim(datestr)
  k1=0
  do kk=1,lens
    if(datestr(kk:kk)==':')then
      k1=kk
      exit
    endif
  enddo
  date_found=(k1>0)
  if(.not.date_found)return

  read(datestr(1:2),*)mon; read(datestr(4:5),*)day
  doy1=yearday(1981,mon,day)

  read(datestr(k1+1:k1+2),*)mon; read(datestr(k1+4:k1+5),*)day
  doy2=yearday(1981,mon,day)

  end function Extract_dates2doy
!------------------------------------------------------------------------------------------

  function check_warming_dates(year,I,J)result(ans)
  !
  !Description
  !check if day I is in the warming period
  implicit none
  integer, intent(in) :: year,I
  integer, intent(in) :: J

  logical :: ans

  if(.not.lsoil_warming)then
    ans=.false.
    return
  endif

  ans=(year==warm_yearb .and. I>=warm_doyb) .or. (year==warm_yeare .and. I<=warm_doye) &
    .or. (year>warm_yearb .and. year <warm_yeare)
  
  if(I==1 .and. J==1)then
    
    leap_cor=0
    if(seas_doyb>59)then
      leap_cor(1)=iisleap(year)
    endif

    if(seas_doyb<seas_doye)then
      !within a year
      leap_cor(2)=leap_cor(1)
    else
      !over two years
      if(seas_doye>59)then
        leap_cor(2)=iisleap(year+1)
      endif      
    endif
  endif
  return
  end function check_warming_dates
!------------------------------------------------------------------------------------------

  subroutine config_soil_warming(warming_exp)

  !!
  ! Description
  ! extract soil warming configuration
  
  implicit none
  character(len=*), intent(in) :: warming_exp
  character(len=*), parameter :: subname='config_soil_warming'
  character(len=120) :: str_loc
  integer :: k, k1, it, nl, nl1, jj
  character(len=12) :: variables(10)
  character(len=32) :: values(10)
  integer :: NY,NX,KK
  integer :: num_found

  call PrintInfo('beg '//subname)

  CALL parse_var_val_string(warming_exp, variables, values, num_found)
  
  lsoil_warming=.false.
  iwtype=-1
  
  DO kk=1,num_found  
    
    if(are_strings_equal_icase(variables(kk),'type'))then
      
      if (are_strings_equal_icase(values(kk),'ir_heating'))then
        lsoil_warming=.true.
        iwtype=1
        exit
      elseif(are_strings_equal_icase(values(kk),'cable_heating'))then
        lsoil_warming=.true.
        iwtype=2
        exit
      elseif(are_strings_equal_icase(values(kk),'open_top_chamber'))then
        lsoil_warming=.true.
        iwtype=3
        exit
      endif
    elseif(are_strings_equal_icase(variables(kk),'types'))then
      
      !combined heating, cable heating can be combined with IR or OTC
      if(is_substring_present('ir_heating',values(kk),ignore_case=.true.) .and. &
        is_substring_present('cable_heating',values(kk),ignore_case=.true.)) then
        lsoil_warming=.true.        
        iwtype=21
        exit
      elseif(is_substring_present('open_top_chamber',values(kk),ignore_case=.true.) .and. &      
        is_substring_present('cable_heating',values(kk),ignore_case=.true.)) then
        lsoil_warming=.true.
        iwtype=23  
        exit
      endif  
    endif
  ENDDO

  if(.not.lsoil_warming)return

  NY_w=1;NX_w=1
  !identify warming grid
  DO KK=1,num_found
    if(are_strings_equal_icase(variables(kk),'loc'))then
      read(values(kk),*)NY_w,NX_w
      exit 
    endif
  enddo

  select case(iwtype)
  case (1)
    !IR heating
    call set_IR_heating(variables,values,num_found)
  case (2)
    !cable heating in soil
    call set_cable_heating(variables,values,num_found)
  case (3)
    !open top chamber heating
    call set_open_top_chamber(variables,values,num_found)
  case (21)
    !IR+cable
    call set_IR_heating(variables,values,num_found)  
    call set_cable_heating(variables,values,num_found)
  case (23) 
    !OTC+cable
    call set_open_top_chamber(variables,values,num_found)
    call set_cable_heating(variables,values,num_found)    
  case default
    write(iulog,*)"Requested warming type not supported"
  end select
  call PrintInfo('end '//subname)
  end subroutine config_soil_warming
!------------------------------------------------------------------------------------------
  subroutine set_IR_heating(variables,values,num_found)
  !
  !Description
  !  warming_exp = 'loc[NY,NX];type[IR_heating];Power[60Wm-2];beg_time[2014/01/01];end_time[2018/12/31];htime[day];season[growth]'
  ! season can also be set as season[date1:date2], where date1 and date2 are like 05/03, if date1 is later than date2, warming
  ! is applied over the time crossing two years.
  implicit none
  character(len=*), intent(in) :: variables(:)
  character(len=*), intent(in) :: values(:)
  integer, intent(in) :: num_found
  character(len=*), parameter :: subname='set_IR_heating'
  character(len=5) :: unit_out
  integer :: KK

  call PrintInfo('beg '//subname)
  DO KK=1,num_found
    select case(trim(variables(kk)))
    case ('power')
      if(extract_number_and_unit(values(kk), warm_IR_power, unit_out))then
        if(trim(unit_out)=='Wm-2')then
          !convert from W m-2 to MJ/m2/hr
          warm_IR_power=warm_IR_power*3600._r8*1.e-6_r8 
        endif
      endif
    case('dT')
      if(extract_number_and_unit(values(kk), warm_dTK, unit_out))then
        !if Fahrenheit, scale it to K 
        if (trim(unit_out)=='F')warm_dTK=warm_dTK*0.5556_r8
      else
        call endrun("warming temperature not set in "//mod_filename,__LINE__)
      endif            
    case ('beg_time')
      call Extract_year_doy(trim(values(kk)),warm_yearb,warm_doyb) 
    case ('end_time')
      call Extract_year_doy(trim(values(kk)),warm_yeare,warm_doye)
    case ('htime')
      select case (trim(values(kk)))
      case ('day')
        ihtime=1
      case ('night')
        ihtime=2
      case default
        ihtime=0
      end select
    case ('season')
      if(trim(values(kk))=='growth')then
        !growing season only
        igrowth=1
      elseif(Extract_dates2doy(values(kk),seas_doyb,seas_doye))then
        if(seas_doyb<seas_doye)then
          igrowth=2
        else
          igrowth=3
        endif
      else
        !whole year
        igrowth=0
      endif
    end select
  ENDDO
  write(iulog,*)'IR heating uisng ',warm_dTK, 'K'
  
  call PrintInfo('end '//subname)
  end subroutine set_IR_heating

!------------------------------------------------------------------------------------------
  subroutine set_cable_heating(variables,values,num_found)
  implicit none
  character(len=*), intent(in) :: variables(:)
  character(len=*), intent(in) :: values(:)
  integer, intent(in) :: num_found  
  character(len=*), parameter :: subname='set_cable_heating'
  integer :: KK
  character(len=2) :: unit_out
  !example configuration
  !  warming_exp = 'loc[NY,NX];type[Cable_heating];dT[4K];Depth[1m];hist_ctrl[Blodget.ctrl.ecosim.h1.xxxx-01-00-00000.nc];' &
  !    //'beg_time[2014/01/01];end_time[2018/12/31];htime[all];season[growth]'
  call PrintInfo('beg '//subname)

  DO KK=1,num_found
    select case(trim(variables(kk)))
    case ('dT')
      if(extract_number_and_unit(values(kk), warm_dTK, unit_out))then
        !if Fahrenheit, scale it to K 
        if (trim(unit_out)=='F')warm_dTK=warm_dTK*0.5556_r8
      else
        call endrun("warming temperature not set in "//mod_filename,__LINE__)
      endif      
    case ('hist_ctrl')
      fname_warming_Tref=trim(values(kk))
    case ('beg_time')
      call Extract_year_doy(trim(values(kk)),warm_yearb,warm_doyb) 
    case ('end_time')
      call Extract_year_doy(trim(values(kk)),warm_yeare,warm_doye)
    case ('htime')
      select case (trim(values(kk)))
      case ('day')
        ihtime=1
      case ('night')
        ihtime=2
      case default
        ihtime=0
      end select
    case ('season')
      if(trim(values(kk))=='growth')then
        !growing season only
        igrowth=1
      elseif(Extract_dates2doy(values(kk),seas_doyb,seas_doye))then
        if(seas_doyb<seas_doye)then
          igrowth=2
        else
          igrowth=3
        endif
      else
        !whole year
        igrowth=0
      endif
    case default
    end select
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine set_cable_heating

!------------------------------------------------------------------------------------------
  subroutine set_open_top_chamber(variables,values,num_found)
  !  
  !Description
  ! warming_exp = 'loc[NY,NX];type[Open_top_chamber];wind[0.4m/s];EMSReduct[4%];taus[0.9];beg_time[2014/01/01];end_time[2018/12/31];season[growth]'  
  implicit none
  character(len=*), intent(in) :: variables(:)
  character(len=*), intent(in) :: values(:)
  integer, intent(in) :: num_found
  character(len=*), parameter :: subname='set_open_top_chamber'
  integer :: KK
  character(len=8) :: unit_out

  call PrintInfo('beg '//subname)

  wind_otc       = 0.4_r8 !default to 0.4 m/s
  ems_reduct_otc = 0.07 !default reduce 7%
  DO KK             = 1, num_found
    select case(trim(variables(kk)))
    case ('beg_time')
      call Extract_year_doy(trim(values(kk)),warm_yearb,warm_doyb) 
    case ('end_time')
      call Extract_year_doy(trim(values(kk)),warm_yeare,warm_doye)
    case ('wind')
      if(extract_number_and_unit(values(kk), wind_otc, unit_out))then
      !do unit conversion when necessary
      endif
    case ('taus') 
      !solar transmissivity
      if(extract_number_and_unit(values(kk), taus_otc, unit_out))then
        if(trim(unit_out)=='%')taus_otc=taus_otc*0.01_r8
      endif 
    case ('EMSReduct')  
      write(*,*)values(kk)
      if(extract_number_and_unit(values(kk), ems_reduct_otc, unit_out))then
        if(trim(unit_out)=='%')ems_reduct_otc=ems_reduct_otc*0.01_r8
      endif    
    case ('season')
      if(trim(values(kk))=='growth')then
        !growing season only
        igrowth=1
      elseif(Extract_dates2doy(values(kk),seas_doyb,seas_doye))then
        if(seas_doyb<seas_doye)then
          igrowth=2
        else
          igrowth=3
        endif
      else
        !whole year
        igrowth=0
      endif
    case default
    end select
  ENDDO
  write(iulog,*)'otc warming with wind speed ',wind_otc, 'm/s',' emissivity reduction ',ems_reduct_otc, &
    'and solar radiation transmissivity ',taus_otc 
  
  call PrintInfo('end '//subname)
  end subroutine set_open_top_chamber

!------------------------------------------------------------------------------------------
  logical function get_warming_fname(year,fname)result(ans)
  !
  !get warming file name for year 
  implicit none  
  integer, intent(in) :: year
  character(len=256), intent(out) :: fname
  integer :: jj

  if(iwtype/=2 .and. iwtype<20)then
    ans=.false.
    return
  endif
  ans=.true.
  write(fname,'(A)')trim(fname_warming_Tref)

  do jj = 1, 256
    if(fname(jj:jj+3)=='xxxx')then
      write(fname(jj:jj+3),'(I4)')year
      exit
    endif
  enddo
  end function get_warming_fname

!------------------------------------------------------------------------------------------

  subroutine apply_OTC_warming(I,J,NHW,NHE,NVN,NVS)
  !
  !Description
  !set up cable warming reference tempeature
  use ClimForcDataType,  only: SineSunInclAngle_col,EMS_scalar_col,srad_scalar_col,WindSpeedAtm_col
  use PlantDataRateType, only: GrossCO2Fix_pft
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: NY,NX,L
  integer :: ij
  logical :: daytime_chk
  logical :: season_chk
  
  
  DO NX=NHW,NHE
    DO NY=NVN,NVS
      !make sure it is otc warming
      if(mod(warm_type(NY,NX),10)==3)then
        
        season_chk=check_warming_season(I,NY,NX)
        if(season_chk)then
          !Ref: Flesch and Grant (1992), Wind Flow within Open-Top Growth Chambers and the Gas Exchange Implications
          !When Open Top Chamber (OTC) is applied, the typical wind speed inside the chamber ranges from about 0.3 to 0.8 m/s based on several studies.
          WindSpeedAtm_col(NY,NX)=wind_otc*3600._r8    
          !Open Top Chambers (OTCs) can reduce the effective long-wave radiation emissivity of the plant canopy by about 5% to 10%.
          EMS_scalar_col(NY,NX) = 1._r8-ems_reduct_otc   !reduce canopy emissivity
          srad_scalar_col(NY,NX)= taus_otc               !reduces incoming solar radiation, it tries not to differentiate direct and diffuse radiation
        endif
      endif
    ENDDO
  ENDDO    
  END subroutine apply_OTC_warming
!------------------------------------------------------------------------------------------

  subroutine apply_soil_cable_warming(I,J,NHW,NHE,NVN,NVS)
  !
  !Description
  !set up cable warming reference tempeature
  use ClimForcDataType,  only: SineSunInclAngle_col
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: NY,NX,L
  integer :: ij
  logical :: daytime_chk
  logical :: season_chk
  
  DO NX=NHW,NHE
    DO NY=NVN,NVS
      if(warm_type(NY,NX)==2 .or.warm_type(NY,NX)>20)then
        !it is cable warming, use solar angle to determine day or night 
        daytime_chk=ihtime==0 .or. (ihtime==1 .and. SineSunInclAngle_col(NY,NX) > 1.e-2_r8) &
          .or. (ihtime==2 .and. SineSunInclAngle_col(NY,NX) <= 1.e-2_r8)

        !use plant GPP to test growing season 
        season_chk=check_warming_season(I,NY,NX)
        if(daytime_chk .and. season_chk)then
          ij=(I-1)*24+J
          DO L=NU_col(NY,NX),NL_col(NY,NX)
            IF(CumDepz2LayBottom_vr(L,NY,NX).LE.warm_depz)THEN
              warm_LL(NY,NX)=L
              TKS_ref_vr(ij,L,NY,NX)=TKS_ref_vr(ij,L,NY,NX)+warm_dTK
            ENDIF    
          ENDDO
        endif
      endif  
    ENDDO
  ENDDO
  end subroutine apply_soil_cable_warming
!------------------------------------------------------------------------------------------
  subroutine apply_IR_warming(I,J,NHW,NHE,NVN,NVS)
  !
  !Description
  !Apply infrared warming
  use ClimForcDataType,  only: TairK_col,SineSunInclAngle_col
  use GridDataType,      only: AREA_3D, NU_col
  use PlantDataRateType, only: GrossCO2Fix_pft
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: NY,NX
  logical :: daytime_chk
  logical :: season_chk

  DO NX=NHW,NHE
    DO NY=NVN,NVS
      if(mod(warm_type(NY,NX),10)==1)then      
        !it is ir heating warming, use solar angle to determine day or night 
        daytime_chk= ihtime==0 .or. (ihtime==1 .and. SineSunInclAngle_col(NY,NX) > 1.e-2_r8) &
          .or. (ihtime==2 .and. SineSunInclAngle_col(NY,NX) <= 1.e-2_r8)

        season_chk=check_warming_season(I,NY,NX)
        if(daytime_chk .and. season_chk)then
          TairK_col(NY,NX)=TairK_col(NY,NX)+warm_dTK          
        endif  
      endif
    ENDDO
  ENDDO
  end subroutine apply_IR_warming
!------------------------------------------------------------------------------------------

  function is_warming_layerL(L,NY,NX)result(ans)

  implicit none 
  integer, intent(in) :: L,NY,NX
  logical :: ans

  ans=warm_LL(NY,NX)>=L
  end function is_warming_layerL

!------------------------------------------------------------------------------------------
  pure function check_warming_season(I,NY,NX)result(ans)
  !
  !Description
  !check if day I is within the specified warming dates
  use PlantDataRateType, only: GrossCO2Fix_pft  

  implicit none
  integer, intent(in) :: I,NY,NX
  logical :: ans

  if(igrowth==2)then
    !within a year    
    ans=seas_doyb+leap_cor(1)<=I .and. I<=seas_doye+leap_cor(2)
  elseif(igrowth==3)then
    !over two years
    ans=I>=seas_doyb+leap_cor(1) .or. I<=seas_doye+leap_cor(2)
  else
    !floating growing season, using gpp to check
    ans = igrowth==0 .or. (igrowth==1 .and. MAXVAL(GrossCO2Fix_pft(:,NY,NX))>0._r8)      
  endif
  end function check_warming_season
end module EcosysWarmingMod

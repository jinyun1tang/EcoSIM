program ClimTransformer
  use ClimReadMod
  use data_kind_mod  , only : r8 => DAT_KIND_R8
  use ClimForcDataType
  use abortutils     , only : endrun
  use minimathmod    , only : isLeap
  use fileUtil       , only : file_exists, open_safe
  use EcoSIMHistMod  , only : PREFIX
  USE EcoSIMCtrlMod  , ONLY : LVERB
  use FlagDataType   , only : IWTHR
  use data_const_mod , only : spval  => DAT_CONST_SPVAL
  use ncdio_pio
  use netcdf
  use EcoSIMConfig, only : IFLGW
implicit none
  character(len=*), parameter :: mod_filename = &
  __FILE__

  type(atm_forc_type) :: atmf
  integer :: I,IX,II,JJ
  character(len=1) :: TTYPE
  integer :: NTX,NFX,LYR
  integer :: iyear
  character(len=128) :: clmfile
  character(len=128) :: flistnm
  character(len=128) :: tfilenm
  character(len=128) :: ncfilenm,progname
  character(len=10)   :: clmtype    !D for daily, H for hourly
  character(len=1) :: cc
  type(file_desc_t) :: ncf
  integer :: yearDimID,dayDimID,hourDimID,recordDimID
  integer :: num_args,k
  real(r8) :: data(1,24,366)
  real(r8) :: datat(366)
  real(r8) :: data2(1,366)
  real(r8) :: data1(1)
  integer :: idata(1)

  CALL GETARG(0,progname)
  clmfile=''

  NTX=1
  NFX=0
!  lverb=.false.
  num_args = iargc()
  if (num_args /= 3)then
    call endrun(msg='Please use code as: '//trim(progname)//' inflist ncfilename hour '&
      //' or '//trim(progname)//' inflist ncfilename day ')
  endif
  CALL GETARG(1,flistnm)
  CALL GETARG(2,ncfilenm)
  call getarg(3,clmtype)

  if(.not.file_exists(flistnm))then
    call endrun(msg='Fail to find the list file '//trim(flistnm)//' in ' &
      //mod_filename,line=__LINE__)
  endif
  prefix='./'
  print*,'open_safe',flistnm
  call OPEN_safe(4,prefix,trim(flistnm),'OLD',mod_filename,__LINE__)

  call ncd_pio_createfile(ncf,trim(ncfilenm))
  call ncd_defdim(ncf,'year',ncd_unlimited,yearDimID)
  call ncd_defdim(ncf,'day',366,dayDimID)
  call ncd_defdim(ncf,'ngrid',1,recordDimID)
  
  if (clmtype(1:1)=='H' .or. clmtype(1:1)=='h')then
    !hourly data
    call ncd_defdim(ncf,'hour',24,hourDimID)
    call ncd_defvar(ncf, 'TMPH', ncd_float, dim1name='ngrid', &
      dim2name='hour',dim3name='day',&
      dim4name='year',long_name='hourly air temperature',  &
      units='oC', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'WINDH', ncd_float,dim1name='ngrid', &
      dim2name='hour',dim3name='day',&
      dim4name='year',long_name='horizontal wind speed',  &
      units='m s^-1', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'RAINH', ncd_float, dim1name='ngrid', &
      dim2name='hour',dim3name='day',&
      dim4name='year',long_name='hourly precipitation',  &
      units='mm m^-2 hr^-1', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'DWPTH', ncd_float, dim1name='ngrid', &
      dim2name='hour',dim3name='day',&
      dim4name='year',long_name='atmospheric vapor pressure',  &
      units='kPa', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'SRADH', ncd_float, dim1name='ngrid', &
      dim2name='hour',dim3name='day',&
      dim4name='year',long_name='Incident solar radiation',  &
      units='W m^-2', missing_value=spval, fill_value=spval)
  elseif(clmtype(1:1)=='D' .or. clmtype(1:1) =='d')then
    !daily data
    call ncd_defvar(ncf, 'TMPX', ncd_float, dim1name='ngrid', &
      dim2name='day',dim3name='year',long_name='daily maximum air temperature',  &
      units='oC', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'TMPN', ncd_float, dim1name='ngrid', &
      dim2name='day',dim3name='year',long_name='daily minimum air temperature',  &
      units='oC', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'WIND', ncd_float,dim1name='ngrid', &
      dim2name='day',dim3name='year',long_name='daily mean horizontal wind speed',  &
      units='m hr^-1', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'RAIN', ncd_float, dim1name='ngrid', &
      dim2name='day',dim3name='year',long_name='daily precipitation',  &
      units='m day^-1', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'DWPT', ncd_float, dim1name='ngrid', &
      dim2name='day',dim3name='year',long_name='atmospheric vapor pressure',  &
      units='kPa', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'SRAD', ncd_float, dim1name='ngrid', &
      dim2name='day',dim3name='year',long_name='Incident daily solar radiation',  &
      units='MJ/day', missing_value=spval, fill_value=spval)    
  endif
  call ncd_defvar(ncf, 'year', ncd_int, dim1name='year', &
    long_name='year AD')

  call ncd_defvar(ncf,'Z0G', ncd_float,dim1name='ngrid',dim2name='year', &
    long_name='windspeed measurement height',units='m',&
    missing_value=spval, fill_value=spval)

  call ncd_defvar(ncf,'IFLGW',ncd_int, dim1name='ngrid',dim2name='year', &
    long_name='flag for raising Z0G with vegeation')

  call ncd_defvar(ncf,'ZNOONG',ncd_float,dim1name='ngrid',dim2name='year',&
    long_name='time of solar noon',units='hour', &
    missing_value=spval, fill_value=spval)

  call ncd_defvar(ncf,'PHRG',ncd_float,dim1name='ngrid',dim2name='year', &
    long_name='pH in precipitation')

  call ncd_defvar(ncf,'CN4RIG',ncd_float, dim1name='ngrid',dim2name='year', &
    long_name='NH4 conc in precip',units='gN m^-3')

  call ncd_defvar(ncf,'CNORIG',ncd_float, dim1name='ngrid',dim2name='year', &
    long_name='NO3 conc in precip',units='gN m^-3')

  call ncd_defvar(ncf,'CPORG',ncd_float, dim1name='ngrid',dim2name='year', &
    long_name='H2PO4 conc in precip',units='gP m^-3')

  call ncd_defvar(ncf,'CALRG',ncd_float, dim1name='ngrid',dim2name='year', &
    long_name='Al conc in precip',units='gAl m^-3')

  call ncd_defvar(ncf,'CFERG',ncd_float, dim1name='ngrid',dim2name='year', &
    long_name='Fe conc in precip',units='gFe m^-3')

  call ncd_defvar(ncf,'CCARG',ncd_float, dim1name='ngrid',dim2name='year', &
    long_name='Ca conc in precip',units='gCa m^-3')

  call ncd_defvar(ncf,'CMGRG',ncd_float, dim1name='ngrid',dim2name='year', &
    long_name='Mg conc in precip',units='gMg m^-3')

  call ncd_defvar(ncf,'CNARG',ncd_float, dim1name='ngrid',dim2name='year', &
    long_name='Na conc in precip',units='gNa m^-3')

  call ncd_defvar(ncf,'CKARG',ncd_float, dim1name='ngrid',dim2name='year', &
    long_name='K conc in precip',units='gK m^-3')

  call ncd_defvar(ncf,'CSORG',ncd_float, dim1name='ngrid',dim2name='year', &
    long_name='SO4 conc in precip',units='gS m^-3')

  call ncd_defvar(ncf,'CCLRG',ncd_float, dim1name='ngrid',dim2name='year', &
    long_name='Cl conc in precip',units='gCl m^-3')

  call ncd_enddef(ncf)
  k=0
  read(4,*,END=1000)iyear
  DO while(.true.)
    read(4,'(A)',END=1000)clmfile
    print*,'processing file ', trim(clmfile), iyear,isleap(iyear)
    TMP_hrly=spval;SWRad_hrly=spval;WINDH=spval;RAINH=spval;DWPTH=spval
    call ReadClim(iyear,clmfile,NTX,NFX,I,IX,TTYPE,atmf)
    K=K+1
    LYR=0
    if(IWTHR==1)then
      !daily climate
       datat(1:366)=DWPT(1,1:366)
       if(.not.isleap(iyear))then
         TMPX(366) = tmpx(365)                        !maximum daily air temperature, [oC]
         TMPN(366) = tmpn(365)                        !minimum daily air temperature, [oC]
         SRAD(366) = srad(365)                        !daily solar radiation, [MJ m-2 d-1]
         RAIN(366) = rain(365)                        !daily precipitation, [mm d-1 ]
         WIND(366) = wind(365)                         !daily wind travel, [m d-1]
         datat(366)= DWPT(1,365)            !daily dewpoint temperature, [oC]
       endif
       call reshape2(TMPX,data2); call ncd_putvar(ncf,'TMPX',k,data2)
       call reshape2(TMPN,data2); call ncd_putvar(ncf,'TMPN',k,data2)
       call reshape2(SRAD,data2); call ncd_putvar(ncf,'SRAD',k,data2)
       call reshape2(RAIN,data2); call ncd_putvar(ncf,'RAIN',k,data2)
       call reshape2(WIND,data2); call ncd_putvar(ncf,'WIND',k,data2)
       call reshape2(datat,data2);call ncd_putvar(ncf,'DWPT',k,data2)

    elseif(IWTHR==2)then
      !hourly climate
      if(isleap(iyear))then
        !correct wrong data 
        DO JJ=1,24
          if(any((/TMP_hrly(jj,366),windh(jj,366),dwpth(jj,366),rainh(jj,366),swrad_hrly(jj,366)/)>1.e20))then
            TMP_hrly(jj,366)=TMP_hrly(jj,365)
            windh(jj,366)=windh(jj,365)
            dwpth(jj,366)=dwpth(jj,365)
            rainh(jj,366)=rainh(jj,365)
            swrad_hrly(jj,366)=swrad_hrly(jj,365)
          endif
        enddo    
      endif
      if(isleap(iyear))LYR=1

      DO II=1,365+LYR
        DO JJ=1,24
          SWRad_hrly(JJ,II)=SWRad_hrly(JJ,II)*1.e6/3600._r8  !convert from MJ/m2 into W m^-2
          WINDH(JJ,II)=WINDH(JJ,II)/3600._r8       !convert from m/hour to m s^-1
          RAINH(JJ,II)=RAINH(JJ,II)*1.e3_r8        !convert from m/hr into mm hr^-1
        enddo
      ENDDO

      call reshape1(TMP_hrly,data); call ncd_putvar(ncf,'TMPH',k,data)
      call reshape1(WINDH,data); call ncd_putvar(ncf,'WINDH',k,data)
      call reshape1(DWPTH,data); call ncd_putvar(ncf,'DWPTH',k,data)
      call reshape1(RAINH,data); call ncd_putvar(ncf,'RAINH',k,data)
      call reshape1(SWRad_hrly,data); call ncd_putvar(ncf,'SRADH',k,data)
    endif
    data1(1)=atmf%Z0G; call ncd_putvar(ncf,'Z0G',k,data1)
    data1(1)=atmf%ZNOONG; call ncd_putvar(ncf,'ZNOONG',k,data1)
    data1(1)=atmf%PHRG; call ncd_putvar(ncf,'PHRG',k,data1)
    data1(1)=atmf%CN4RIG; call ncd_putvar(ncf,'CN4RIG',k,data1)
    data1(1)=atmf%CNORIG; call ncd_putvar(ncf,'CNORIG',k,data1)
    data1(1)=atmf%CPORG; call ncd_putvar(ncf,'CPORG',k,data1)
    data1(1)=atmf%CALRG; call ncd_putvar(ncf,'CALRG',k,data1)
    data1(1)=atmf%CFERG; call ncd_putvar(ncf,'CFERG',k,data1)
    data1(1)=atmf%CCARG; call ncd_putvar(ncf,'CCARG',k,data1)
    data1(1)=atmf%CMGRG; call ncd_putvar(ncf,'CMGRG',k,data1)
    data1(1)=atmf%CNARG; call ncd_putvar(ncf,'CNARG',k,data1)
    data1(1)=atmf%CKARG; call ncd_putvar(ncf,'CKARG',k,data1)
    data1(1)=atmf%CSORG; call ncd_putvar(ncf,'CSORG',k,data1)
    data1(1)=atmf%CCLRG; call ncd_putvar(ncf,'CCLRG',k,data1)
    idata(1)=IFLGW; call ncd_putvar(ncf,'IFLGW',k,idata)
    call ncd_putvar(ncf,'year',k,iyear)
    iyear=iyear+1
  ENDDO

1000  close(4)
  call ncd_pio_closefile(ncf)
contains

  subroutine reshape1(arr2d,arr3d)
  implicit none
  real(r8), dimension(:,:), intent(in) :: arr2d
  real(r8), dimension(:,:,:),intent(out) :: arr3d

  DO ii=1,366
    DO jj=1,24
      arr3d(1,jj,ii)=arr2d(jj,ii)
    ENDDO
  ENDDO
  end subroutine reshape1
!----------------------------------------------------------------------

  subroutine reshape2(arr1d,arr2d)
  implicit none
  real(r8), dimension(:), intent(in) :: arr1d
  real(r8), dimension(:,:),intent(out) :: arr2d
  integer :: ii
  DO ii=1,366
    arr2d(1,ii)=arr1d(ii)
  ENDDO
  end subroutine reshape2

end program ClimTransformer

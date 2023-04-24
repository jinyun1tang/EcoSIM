program main
  use ClimReadMod
  use data_kind_mod  , only : r8 => DAT_KIND_R8
  use ClimForcDataType
  use abortutils     , only : endrun
  use minimathmod    , only : isLeap
  use fileUtil       , only : file_exists, open_safe
  use EcoSIMHistMod  , only : PREFIX
  USE EcoSIMCtrlMod  , ONLY : LVERB
  use data_const_mod , only : spval  => DAT_CONST_SPVAL
  use ncdio_pio
  use netcdf
  use EcoSIMConfig, only : IFLGW
implicit none
  character(len=*), parameter :: mod_filename = __FILE__

  type(atm_forc_type) :: atmf
  integer :: I,IX,II,JJ
  character(len=1) :: TTYPE
  integer :: L,NTX,NFX,LYR
  integer :: iyear
  character(len=128) :: clmfile
  character(len=128) :: flistnm
  character(len=128) :: tfilenm
  character(len=128) :: ncfilenm,progname
  character(len=1) :: cc
  type(file_desc_t) :: ncf
  integer :: yearDimID,dayDimID,hourDimID,recordDimID
  integer :: num_args,k
  real(r8) :: data(1,24,366)
  real(r8) :: data1(1)
  integer :: idata(1)

  CALL GETARG(0,progname)
  clmfile=''
  L=1
  NTX=1
  NFX=0
!  lverb=.false.
  num_args = iargc()
  if (num_args/=2)then
    call endrun(msg='Please use code as: '//trim(progname)//' inflist ncfilename')
  endif
  CALL GETARG(1,flistnm)
  CALL GETARG(2,ncfilenm)
  if(.not.file_exists(flistnm))then
    call endrun(msg='Fail to find the list file '//trim(flistnm)//' in ' &
      //mod_filename,line=__LINE__)
  endif
  prefix='./'
  call OPEN_safe(4,prefix,trim(flistnm),'OLD',mod_filename,__LINE__)

  call ncd_pio_createfile(ncf,trim(ncfilenm))
  call ncd_defdim(ncf,'year',ncd_unlimited,yearDimID)
  call ncd_defdim(ncf,'day',366,dayDimID)
  call ncd_defdim(ncf,'hour',24,hourDimID)
  call ncd_defdim(ncf,'ngrid',1,recordDimID)

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
    dim4name='year',long_name='Total precipitation',  &
    units='mm m^-2 hr^-1', missing_value=spval, fill_value=spval)

  call ncd_defvar(ncf, 'DWPTH', ncd_float, dim1name='ngrid', &
    dim2name='hour',dim3name='day',&
    dim4name='year',long_name='atmospheric vapor pressure',  &
    units='kPa', missing_value=spval, fill_value=spval)

  call ncd_defvar(ncf, 'SRADH', ncd_float, dim1name='ngrid', &
    dim2name='hour',dim3name='day',&
    dim4name='year',long_name='Incident solar radiation',  &
    units='W m^-2', missing_value=spval, fill_value=spval)

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
    print*,'processing file ', trim(clmfile), iyear
    TMPH=spval;SRADH=spval;WINDH=spval;RAINH=spval;DWPTH=spval
    call ReadClim(iyear,clmfile,NTX,NFX,L,I,IX,TTYPE,atmf)
    LYR=0
    if(isleap(iyear))LYR=1
    DO II=1,365+LYR
      DO JJ=1,24
        SRADH(JJ,II)=SRADH(JJ,II)*1.e6/3600._r8  !convert from MJ/m2 into W m^-2
        WINDH(JJ,II)=WINDH(JJ,II)/3600._r8       !convert from m/hour to m s^-1
        RAINH(JJ,II)=RAINH(JJ,II)*1.e3_r8        !convert from m/hr into mm hr^-1
      enddo
    ENDDO
    k=k+1
    call reshape1(TMPH,data); call ncd_putvar(ncf,'TMPH',k,data)
    call reshape1(WINDH,data); call ncd_putvar(ncf,'WINDH',k,data)
    call reshape1(DWPTH,data); call ncd_putvar(ncf,'DWPTH',k,data)
    call reshape1(RAINH,data); call ncd_putvar(ncf,'RAINH',k,data)
    call reshape1(SRADH,data); call ncd_putvar(ncf,'SRADH',k,data)
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
end program main

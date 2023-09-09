module ForcWriterMod
  use SoilPhysDataType
  use SOMDataType
  USE MicrobialDataType
  use SoilBGCDataType
  use GridDataType
  use GridConsts
  use AqueChemDatatype
  use SoilPropertyDataType
  use ClimForcDataType
  use SoilWaterDataType
  use ncdio_pio
  use EcoSIMConfig, only : jcplx => jcplxc
  use data_const_mod, only : spval  => DAT_CONST_SPVAL
implicit none

  type, public :: bgc_forc_config_type
    logical :: laddband
    integer :: year
    integer :: doy
    integer :: layer
    character(len=64) :: bgc_fname
  end type bgc_forc_config_type
  type(bgc_forc_config_type) :: bgc_forc_conf

  logical, public :: do_bgcforc_write
  public :: WriteBBGCForc
  contains

!------------------------------------------------------------------------------------------
  subroutine WriteBBGCForc(doy,year)
  !
  ! DESCRIPTION:
  ! Write initial condition for the batch soil bgc model.
  ! Because the single layer bgc does not have band information.
  ! The band nutrients is set to zero, or add to the bulk soil
  implicit none
  integer,intent(in) :: doy
  integer,intent(in) :: year
  integer :: NY,NX,L
  integer :: k,kl,nl
  type(file_desc_t) :: ncf
  integer :: recordDimID
  real(r8),pointer :: data1d(:)
  real(r8), allocatable :: dat1d(:)

  if(year==bgc_forc_conf%year .and. doy==bgc_forc_conf%doy)then
    NY=1;NX=1;L=bgc_forc_conf%Layer
    write(*,*)'write bbgc forc on day ',doy,'year',year
    call ncd_pio_createfile(ncf, trim(bgc_forc_conf%bgc_fname))
    call ncd_defdim(ncf,'jcplx',jcplx,recordDimID)
    call ncd_defdim(ncf,'jsken',jsken,recordDimID)
    call ncd_defdim(ncf,'ndbiomcp',2,recordDimID)
    call ncd_defdim(ncf,'nlbiomcp',3,recordDimID)
    call ncd_defdim(ncf,'NMICBSO',NMICBSO,recordDimID)
    call ncd_defdim(ncf,'NMICBSA',NMICBSA,recordDimID)

    call ncd_defvar(ncf, 'pH', ncd_float, long_name='soil pH',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'VLSoilPoreMicP', ncd_float, long_name='volume of soil layer',  &
            units='m3 d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'ORGC', ncd_float, long_name='total soil organic C',  &
            units='gC d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CFOMC', ncd_float, dim1name='ndbiomcp', &
            long_name='allocation coefficient to humus fractions',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'VLSoilMicP', ncd_float, long_name='micropore volume',  &
            units='m3 d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'BKVL', ncd_float, long_name='mass of soil layer',  &
            units='Mg d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'POROS', ncd_float, long_name='soil porosity',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'FC', ncd_float, long_name='water contents at field capacity',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'SRP', ncd_float, long_name='shape parameter for water desorption',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'WP', ncd_float, long_name='water contents at wilting point',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'EHUM',ncd_float, &
            long_name='partitioning coefficient between humus and microbial residue',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'DLYR3', ncd_float, long_name='thickness of soil layer',  &
            units='m', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CEC', ncd_float, long_name='soil cation exchange capacity',  &
            units='cmol kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'XCEC', ncd_float, long_name='soil cation exchange capacity',  &
            units='mol d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'AEC', ncd_float, long_name='soil anion exchange capacity',  &
            units='cmol kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'XAEC', ncd_float, long_name='soil anion exchange capacity',  &
            units='mol d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CFE', ncd_float, long_name='soil Fe content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CCA', ncd_float, long_name='soil Ca content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CMG', ncd_float, long_name='soil Mg content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CNA', ncd_float, long_name='soil Na content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CKA', ncd_float, long_name='soil K content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CSO4', ncd_float, long_name='soil SO4 content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CCL', ncd_float, long_name='soil Cl content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CAL', ncd_float, long_name='soil Al content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'ZMG', ncd_float, long_name='soil aqueous Mg content micropore',  &
            units='mol d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'ZNA', ncd_float, long_name='soil aqueous Na content micropore',  &
            units='mol d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'ZKA', ncd_float, long_name='soil aqueous K content micropore',  &
            units='mol d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CALPO', ncd_float, long_name='soil AlPO4 content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CFEPO', ncd_float, long_name='soil FePO4 content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CCAPD', ncd_float, long_name='soil CaHPO4 content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'EPOC', ncd_float, &
            long_name='partitioning coefficient between POC and litter',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CCAPH', ncd_float, long_name='soil apatite content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CALOH', ncd_float, long_name='soil AlOH3 content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CFEOH', ncd_float, long_name='soil FeOH3 content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CCACO', ncd_float, long_name='soil CaCO3 content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CCASO', ncd_float, long_name='soil CaSO4 content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'FOSRH', ncd_float, dim1name='jcplx',&
            long_name='fraction of total organic C in complex',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OQC', ncd_float, dim1name='jcplx',&
            long_name='dissolved organic C micropore',  &
            units='gC d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OQN', ncd_float, dim1name='jcplx',&
            long_name='dissolved organic N micropore',  &
            units='gN d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OQP', ncd_float, dim1name='jcplx',&
            long_name='dissolved organic P micropore',  &
            units='gP d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OQA', ncd_float, dim1name='jcplx',&
            long_name='dissolved acetate micropore',  &
            units='gC d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OHA', ncd_float, dim1name='jcplx',&
            long_name='adsorbed acetate',  &
            units='gC d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OHC', ncd_float, dim1name='jcplx',&
            long_name='adsorbed soil C',  &
            units='gC d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OHN', ncd_float, dim1name='jcplx',&
            long_name='adsorbed soil N',  &
            units='gN d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OHP', ncd_float, dim1name='jcplx',&
            long_name='adsorbed soil P',  &
            units='gP d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'ATCS',ncd_float,long_name='Mean annual temperature',  &
            units='oC', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'CNOSC', ncd_float, dim1name='jsken',&
            dim2name='jcplx',long_name='N:C ratios of SOM kinetic components in each complex',  &
            units='gN/gC', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CPOSC', ncd_float, dim1name='jsken',&
            dim2name='jcplx',long_name='P:C ratios of SOM kinetic components in each complex',  &
            units='gN/gC', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OSC', ncd_float, dim1name='jsken',&
            dim2name='jcplx',long_name='humus soil C in each complex',  &
            units='gC d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OSA', ncd_float, dim1name='jsken',&
            dim2name='jcplx',long_name='colonized soil C in each complex',  &
            units='gC d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OSN', ncd_float, dim1name='jsken',&
            dim2name='jcplx',long_name='humus soil N in each complex',  &
            units='gC d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OSP', ncd_float, dim1name='jsken',&
            dim2name='jcplx',long_name='humus soil P in each complex',  &
            units='gP d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'ORC', ncd_float, dim1name='ndbiomcp',&
            dim2name='jcplx',long_name='microbial residue C in each complex',  &
            units='gC d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'ORN', ncd_float, dim1name='ndbiomcp',&
            dim2name='jcplx',long_name='microbial residue N in each complex',  &
            units='gN d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'ORP', ncd_float, dim1name='ndbiomcp',&
            dim2name='jcplx',long_name='microbial residue P in each complex',  &
            units='gP d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OMC', ncd_float, dim1name='nlbiomcp',&
            dim2name='NMICBSO',dim3name='jcplx',&
            long_name='microbial biomass C in each complex',  &
            units='gC d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OMN', ncd_float, dim1name='nlbiomcp',&
            dim2name='NMICBSO',dim3name='jcplx',&
            long_name='microbial biomass N in each complex',  &
            units='gN d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OMP', ncd_float, dim1name='nlbiomcp',&
            dim2name='NMICBSO',dim3name='jcplx',&
            long_name='microbial biomass P in each complex',  &
            units='gP d-2', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'OMCff', ncd_float, dim1name='nlbiomcp',dim2name='NMICBSA',&
            long_name='microbial biomass C in autotrophic complex',  &
            units='gC d-2', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'OMNff', ncd_float, dim1name='nlbiomcp',dim2name='NMICBSA',&
            long_name='microbial biomass N in autotrophic complex',  &
            units='gN d-2', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'OMPff', ncd_float, dim1name='nlbiomcp',dim2name='NMICBSA',&
            long_name='microbial biomass P in autotrophic complex',  &
            units='gP d-2', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'BKDS', ncd_float,long_name='soil bulk density',&
            units='Mg m-3', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'THETY', ncd_float,long_name='air-dry water content',&
            units='m3 m-3', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'PSIMX', ncd_float,long_name='log water potential at field capacity',&
            units='log(MPa)', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'PSIMS', ncd_float,long_name='log water potential at saturation',&
            units='log(MPa)', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'PSIMD', ncd_float, &
            long_name='log water potential at wilting point - log water potential at field capacity',&
            units='log(MPa)', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'PSISD', ncd_float, &
            long_name='log water potential at saturation - log water potential at field capacity',&
            units='log(MPa)', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'PSISE', ncd_float, long_name='water potential at saturation',&
            units='MPa', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'ZNH4S', ncd_float,long_name='NH4 in micropore',  &
            units='gN d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'ZNO3S', ncd_float,long_name='NO3(-) in micropore',  &
            units='gN d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'ZNO2S', ncd_float,long_name='NO2(-) in micropore',  &
            units='gN d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'H2PO4', ncd_float,long_name='H2PO4 in micropore',  &
            units='gP d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'H1PO4', ncd_float,long_name='H1PO4 in micropore',  &
            units='gP d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf,'ATKA', ncd_float, long_name='mean annual air temperature',&
            units='K', missing_value=spval, fill_value=spval)
    call ncd_enddef(ncf)

    call ncd_putvar(ncf,'pH',PH(L,NY,NX))
    call ncd_putvar(ncf,'VLSoilPoreMicP',VLSoilPoreMicP(L,NY,NX))
    call ncd_putvar(ncf,'ORGC',ORGC(L,NY,NX))
    call ncd_putvar(ncf,'CFOMC',CFOMC(:,L,NY,NX))
    call ncd_putvar(ncf,'VLSoilMicP',VLSoilMicP(L,NY,NX))
    call ncd_putvar(ncf,'BKVL',SoilMicPMassLayer(L,NY,NX))
    call ncd_putvar(ncf,'POROS',POROS(L,NY,NX))
    call ncd_putvar(ncf,'FC',FieldCapacity(L,NY,NX))
    call ncd_putvar(ncf,'ATKA',TairKClimMean(NY,NX))
    call ncd_putvar(ncf,'WP',WiltPoint(L,NY,NX))
    call ncd_putvar(ncf,'SRP',SRP(L,NY,NX))
    call ncd_putvar(ncf,'EHUM',EHUM(L,NY,NX))
    call ncd_putvar(ncf,'EPOC',EPOC(L,NY,NX))
    call ncd_putvar(ncf,'THETY',THETY(L,NY,NX))
    call ncd_putvar(ncf,'PSIMX',LOGPSIMX(NY,NX))
    call ncd_putvar(ncf,'PSIMD',LOGPSIMND(NY,NX))
    call ncd_putvar(ncf,'PSIMS',LOGPSIAtSat(NY,NX))
    call ncd_putvar(ncf,'PSISD',LOGPSIMXD(NY,NX))
    call ncd_putvar(ncf,'PSISE',PSISE(L,NY,NX))
    call ncd_putvar(ncf,'DLYR3',DLYR(3,L,NY,NX))
    call ncd_putvar(ncf,'CEC',CEC(L,NY,NX))
    call ncd_putvar(ncf,'XCEC',trcx_solml(idx_CEC,L,NY,NX))
    call ncd_putvar(ncf,'AEC',AEC(L,NY,NX))
    call ncd_putvar(ncf,'XAEC',trcx_solml(idx_AEC,L,NY,NX))
    call ncd_putvar(ncf,'CFE',CFE(L,NY,NX))
    call ncd_putvar(ncf,'CCA',CCA(L,NY,NX))
    call ncd_putvar(ncf,'CMG',CMG(L,NY,NX))
    call ncd_putvar(ncf,'CNA',CNA(L,NY,NX))
    call ncd_putvar(ncf,'CKA',CKA(L,NY,NX))
    call ncd_putvar(ncf,'CSO4',CSO4(L,NY,NX))
    call ncd_putvar(ncf,'CCL',CCL(L,NY,NX))
    call ncd_putvar(ncf,'CAL',CAL(L,NY,NX))
    call ncd_putvar(ncf,'ZMG',trcsa_solml(idsa_Mg,L,NY,NX))
    call ncd_putvar(ncf,'ZNA',trcsa_solml(idsa_Na,L,NY,NX))
    call ncd_putvar(ncf,'ZKA',trcsa_solml(idsa_K,L,NY,NX))
    call ncd_putvar(ncf,'CALPO',CALPO(L,NY,NX))
    call ncd_putvar(ncf,'CFEPO',CFEPO(L,NY,NX))
    call ncd_putvar(ncf,'CCAPD',CCAPD(L,NY,NX))
    call ncd_putvar(ncf,'CCAPH',CCAPH(L,NY,NX))
    call ncd_putvar(ncf,'CALOH',CALOH(L,NY,NX))
    call ncd_putvar(ncf,'CFEOH',CFEOH(L,NY,NX))
    call ncd_putvar(ncf,'CCACO',CCACO(L,NY,NX))
    call ncd_putvar(ncf,'CCASO',CCASO(L,NY,NX))
    call ncd_putvar(ncf,'BKDS',SoiBulkDensity(L,NY,NX))

    call ncd_putvar(ncf,'FOSRH',FOSRH(:,L,NY,NX))
    call ncd_putvar(ncf,'OQC',OQC(:,L,NY,NX))
    call ncd_putvar(ncf,'OQN',OQN(:,L,NY,NX))
    call ncd_putvar(ncf,'OQP',OQP(:,L,NY,NX))
    call ncd_putvar(ncf,'OQA',OQA(:,L,NY,NX))
    call ncd_putvar(ncf,'OHA',OHC(:,L,NY,NX))
    call ncd_putvar(ncf,'OHC',OHC(:,L,NY,NX))
    call ncd_putvar(ncf,'OHN',OHN(:,L,NY,NX))
    call ncd_putvar(ncf,'OHP',OHP(:,L,NY,NX))

    call ncd_putvar(ncf,'CNOSC',CNOSC(:,:,L,NY,NX))
    call ncd_putvar(ncf,'CPOSC',CPOSC(:,:,L,NY,NX))
    call ncd_putvar(ncf,'ATCS',ATCS(NY,NX))
    call ncd_putvar(ncf,'OSC',OSC(:,:,L,NY,NX))
    call ncd_putvar(ncf,'OSA',OSA(:,:,L,NY,NX))
    call ncd_putvar(ncf,'OSN',OSN(:,:,L,NY,NX))
    call ncd_putvar(ncf,'OSP',OSP(:,:,L,NY,NX))
    call ncd_putvar(ncf,'ORC',ORC(:,:,L,NY,NX))
    call ncd_putvar(ncf,'ORN',ORN(:,:,L,NY,NX))
    call ncd_putvar(ncf,'ORP',ORP(:,:,L,NY,NX))

    call ncd_putvar(ncf,'OMC',OMC(:,:,:,L,NY,NX))
    call ncd_putvar(ncf,'OMN',OMN(:,:,:,L,NY,NX))
    call ncd_putvar(ncf,'OMP',OMP(:,:,:,L,NY,NX))

    call ncd_putvar(ncf,'OMCff',OMCff(:,:,L,NY,NX))
    call ncd_putvar(ncf,'OMNff',OMNff(:,:,L,NY,NX))
    call ncd_putvar(ncf,'OMPff',OMPff(:,:,L,NY,NX))

    if(bgc_forc_conf%laddband)then
      call ncd_putvar(ncf,'ZNH4S',trc_solml(ids_NH4,L,NY,NX)+trc_solml(ids_NH4B,L,NY,NX))
      call ncd_putvar(ncf,'ZNO3S',trc_solml(ids_NO3,L,NY,NX)+trc_solml(ids_NO3B,L,NY,NX))
      call ncd_putvar(ncf,'ZNO2S',trc_solml(ids_NO2,L,NY,NX)+trc_solml(ids_NO2B,L,NY,NX))
      call ncd_putvar(ncf,'H2PO4',trc_solml(ids_H2PO4,L,NY,NX)+trc_solml(ids_H2PO4B,L,NY,NX))
      call ncd_putvar(ncf,'H1PO4',trc_solml(ids_H1PO4,L,NY,NX)+trc_solml(ids_H1PO4B,L,NY,NX))
    else
      call ncd_putvar(ncf,'ZNH4S',trc_solml(ids_NH4,L,NY,NX))
      call ncd_putvar(ncf,'ZNO3S',trc_solml(ids_NO3,L,NY,NX))
      call ncd_putvar(ncf,'ZNO2S',trc_solml(ids_NO2,L,NY,NX))
      call ncd_putvar(ncf,'H2PO4',trc_solml(ids_H2PO4,L,NY,NX))
      call ncd_putvar(ncf,'H1PO4',trc_solml(ids_H1PO4,L,NY,NX))
    endif
    call ncd_pio_closefile(ncf)
  endif
  end subroutine WriteBBGCForc

end module ForcWriterMod

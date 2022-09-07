module ForcWriterMod
  use SOMDataType          , only : ORGC,CFOMC
  use SoilPhysDataType     , only : FC
  use SOMDataType
  USE MicrobialDataType
  use SoilBGCDataType
  use GridDataType
  use AqueChemDatatype
  use SoilPropertyDataType
  use ncdio_pio
  use data_const_mod, only : spval  => SHR_CONST_SPVAL
implicit none

!    real(r8) :: CCH4E       !atmospheric CH4 concentration, [g m-3]
!    real(r8) :: COXYE       !atmospheric O2 concentration, [g m-3]
!    real(r8) :: COXQ        !surface irrigation  O2 concentration, [g m-3]
!    real(r8) :: COXR        !precipitation  O2 concentration, [g m-3]
!    real(r8) :: FLQRI       !irrigation flux into surface litter, [m3 d-2 h-1]
!    real(r8) :: FLQRQ       !precipitation flux into surface litter, [m3 d-2 h-1]
!    real(r8) :: OFFSET      !offset for calculating temperature in Arrhenius curves, [oC]

!    real(r8) :: THETY       !air-dry water content, [m3 m-3]
!    real(r8) :: TKS         !temperature in kelvin, [K]
!    real(r8) :: THETW       !volumetric water content [m3 m-3]
!    real(r8) :: pH          !pH value
!    real(r8) :: BKVL        !mass of soil layer	Mg d-2
!    real(r8) :: VOLX        !volume of soil layer	m3 d-2
!    real(r8) :: VOLW        !soil micropore water content [m3 d-2]

!    real(r8) :: VLNOB       !NO3 band volume fracrion, [-]
!    real(r8) :: VLNO3       !NO3 non-band volume fraction,[-] := 1-VLNOB
!    real(r8) :: VLNHB       !NH4 band volume fraction, [-]
!    real(r8) :: VLNH4       !NH4 non-band volume fraction, [-]:=1-VLNHB
!    real(r8) :: VLPOB       !PO4 band volume fracrion, [-]
!    real(r8) :: VLPO4       !PO4 non-band volume fraction,[-]:=1-VLPOB

!    real(r8) :: PSISM       !soil micropore matric water potential, [MPa]
!    real(r8) :: OLSGL       !aqueous O2 diffusivity, [m2 h-1], set in hour1
!    real(r8) :: ORGC        !total soil organic C [g d-2]
!    real(r8) :: CFOMC(1:2)  !allocation coefficient to humus fractions

!    real(r8) :: CNH4B       !NH4 concentration band micropore	[g m-3], derived from ZNH4B
!    real(r8) :: ZNH4B       !NH4 band micropore, [g d-2]
!    real(r8) :: CNH4S       !NH4 concentration non-band micropore	[g m-3], derived from ZNH4S
!    real(r8) :: ZNH4S       !NH4 non-band micropore, [g d-2]
!    real(r8) :: CNO3B       !NO3 concentration band micropore	[g m-3], derived from ZNO3B
!    real(r8) :: ZNO3B       !NO3 band micropore, [g d-2]
!    real(r8) :: CNO3S       !NO3 concentration non-band micropore	[g m-3], derived from ZNO3S
!    real(r8) :: ZNO3S       !NO3 band micropore, [g d-2]
!    real(r8) :: CH2P4       !aqueous PO4 concentration non-band	[g m-3], derived from H2PO4
!    real(r8) :: H2PO4       !PO4 non-band micropore, [g d-2]
!    real(r8) :: CH2P4B      !aqueous PO4 concentration band	[g m-3], derived from H2POB
!    real(r8) :: H2POB       !PO4 band micropore, [g d-2]
!    real(r8) :: CH1P4       !aqueous H1PO4 concentration non-band [g m-3], derived from H1PO4
!    real(r8) :: H1PO4       !soil aqueous HPO4 content micropore non-band, [mol d-2]
!    real(r8) :: CH1P4B      !aqueous H1PO4 concentration band [g m-3], derived from H1POB
!    real(r8) :: H1POB       !soil aqueous HPO4 content micropore band, [mol d-2]
!    real(r8) :: CNO2B       !aqueous HNO2 concentration band [g m-3], derived from ZNO2B
!    real(r8) :: ZNO2B       !NO2  band micropore, [g d-2]
!    real(r8) :: CNO2S       !NO2 concentration non-band micropore	[g m-3], derived from ZNO2S
!    real(r8) :: ZNO2S       !NO2 (nitrate) non-band micropore, [g d-2]
!    real(r8) :: DFGS        !coefficient for dissolution - volatilization, []
!    real(r8) :: FILM        !soil water film thickness , [m]
!    real(r8) :: THETPM      !soil air-filled porosity, [m3 m-3]
!    real(r8) :: TORT        !soil tortuosity, []
!    real(r8) :: VOLPM       !soil air content, [m3 d-2]

!    real(r8) :: EPOC      !partitioning coefficient between POC and litter, [], hour1.f
!    real(r8) :: EHUM      !partitioning coefficient between humus and microbial residue, [], hour1.f
!    real(r8) :: CCO2S     !aqueous CO2 concentration micropore	[g m-3]
!    real(r8) :: CZ2OS     !aqueous N2O concentration micropore	[g m-3]
!    real(r8) :: Z2OS      !aqueous N2O micropore, [g d-2]
!    real(r8) :: COXYS     !aqueous O2 concentration micropore	[g m-3]
!    real(r8) :: OXYS      !aqueous O2  micropore	[g d-2]
!    real(r8) :: COXYG     !gaseous O2 concentration	[g m-3]
!    real(r8) :: CZ2GS     !aqueous N2 concentration micropore	[g m-3]
!    real(r8) :: CH2GS     !aqueous H2 concentration	[g m-3]
!    real(r8) :: H2GS      !aqueous H2 	[g d-2]
!    real(r8) :: CCH4G     !gaseous CH4 concentration	[g m-3]
!    real(r8) :: CH4S      !aqueous CO2  micropore	[g d-2]
!    real(r8) :: SCH4L     !solubility of CH4, [m3 m-3]
!    real(r8) :: SOXYL     !solubility of O2, [m3 m-3]
!    real(r8) :: ZNFN0     !initial nitrification inhibition activity
!    real(r8) :: ZNFNI     !current nitrification inhibition activity

!    real(r8) :: VOLY      !micropore volume, [m3 d-2]
!    real(r8) :: POROS     !soil porosity, [m3 m-3]
!    real(r8) :: FC        !water contents at field capacity, [m3 m-3]
!    real(r8) :: CCLAY
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
    call ncd_defdim(ncf,'JG',JG,recordDimID)
    call ncd_defdim(ncf,'NFGs',NFGs,recordDimID)

    call ncd_defvar(ncf, 'pH', ncd_float, long_name='soil pH',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'VOLX', ncd_float, long_name='volume of soil layer',  &
            units='m3 d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'ORGC', ncd_float, long_name='total soil organic C',  &
            units='gC d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CFOMC', ncd_float, dim1name='ndbiomcp', &
            long_name='allocation coefficient to humus fractions',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'VOLY', ncd_float, long_name='micropore volume',  &
            units='m3 d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'BKVL', ncd_float, long_name='mass of soil layer',  &
            units='Mg d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'POROS', ncd_float, long_name='soil porosity',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'FC', ncd_float, long_name='water contents at field capacity',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'EHUM',ncd_float, &
            long_name='partitioning coefficient between humus and microbial residue',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'DLYR3', ncd_float, long_name='thickness of soil layer',  &
            units='m', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CEC', ncd_float, long_name='soil cation exchange capacity',  &
            units='cmol kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'AEC', ncd_float, long_name='soil anion exchange capacity',  &
            units='cmol kg-1', missing_value=spval, fill_value=spval)
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
    call ncd_defvar(ncf, 'CALPO', ncd_float, long_name='soil AlPO4 content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CFEPO', ncd_float, long_name='soil FePO4 content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CCAPD', ncd_float, long_name='soil CaHPO4 content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
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
    call ncd_defvar(ncf, 'OHC', ncd_float, dim1name='jcplx',&
            long_name='adsorbed soil C',  &
            units='gC d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OHN', ncd_float, dim1name='jcplx',&
            long_name='adsorbed soil N',  &
            units='gN d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OHP', ncd_float, dim1name='jcplx',&
            long_name='adsorbed soil P',  &
            units='gP d-2', missing_value=spval, fill_value=spval)

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
            dim2name='JG',dim3name='NFGs',dim4name='jcplx',&
            long_name='microbial biomass C in each complex',  &
            units='gC d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OMN', ncd_float, dim1name='nlbiomcp',&
            dim2name='JG',dim3name='NFGs',dim4name='jcplx',&
            long_name='microbial biomass N in each complex',  &
            units='gN d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OMP', ncd_float, dim1name='nlbiomcp',&
            dim2name='JG',dim3name='NFGs',dim4name='jcplx',&
            long_name='microbial biomass P in each complex',  &
            units='gP d-2', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'OMCff', ncd_float, dim1name='nlbiomcp',&
            dim2name='JG',dim3name='NFGs',&
            long_name='microbial biomass C in autotrophic complex',  &
            units='gC d-2', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'OMNff', ncd_float, dim1name='nlbiomcp',&
            dim2name='JG',dim3name='NFGs',&
            long_name='microbial biomass N in autotrophic complex',  &
            units='gN d-2', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'OMPff', ncd_float, dim1name='nlbiomcp',&
            dim2name='JG',dim3name='NFGs',&
            long_name='microbial biomass P in autotrophic complex',  &
            units='gP d-2', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'BKDS', ncd_float,long_name='soil bulk density',&
            units='Mg m-3', missing_value=spval, fill_value=spval)

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
    call ncd_enddef(ncf)

    call ncd_putvar(ncf,'pH',PH(L,NY,NX))
    call ncd_putvar(ncf,'VOLX',VOLX(L,NY,NX))
    call ncd_putvar(ncf,'ORGC',ORGC(L,NY,NX))
    call ncd_putvar(ncf,'CFOMC',CFOMC(:,L,NY,NX))
    call ncd_putvar(ncf,'VOLY',VOLY(L,NY,NX))
    call ncd_putvar(ncf,'BKVL',BKVL(L,NY,NX))
    call ncd_putvar(ncf,'POROS',POROS(L,NY,NX))
    call ncd_putvar(ncf,'FC',FC(L,NY,NX))
    call ncd_putvar(ncf,'EHUM',EHUM(L,NY,NX))
    call ncd_putvar(ncf,'DLYR3',DLYR(3,L,NY,NX))
    call ncd_putvar(ncf,'CEC',CEC(L,NY,NX))
    call ncd_putvar(ncf,'AEC',AEC(L,NY,NX))
    call ncd_putvar(ncf,'CFE',CFE(L,NY,NX))
    call ncd_putvar(ncf,'CCA',CCA(L,NY,NX))
    call ncd_putvar(ncf,'CMG',CMG(L,NY,NX))
    call ncd_putvar(ncf,'CNA',CNA(L,NY,NX))
    call ncd_putvar(ncf,'CKA',CKA(L,NY,NX))
    call ncd_putvar(ncf,'CSO4',CSO4(L,NY,NX))
    call ncd_putvar(ncf,'CCL',CCL(L,NY,NX))
    call ncd_putvar(ncf,'CALPO',CALPO(L,NY,NX))
    call ncd_putvar(ncf,'CFEPO',CFEPO(L,NY,NX))
    call ncd_putvar(ncf,'CCAPD',CCAPD(L,NY,NX))
    call ncd_putvar(ncf,'CCAPH',CCAPH(L,NY,NX))
    call ncd_putvar(ncf,'CALOH',CALOH(L,NY,NX))
    call ncd_putvar(ncf,'CFEOH',CFEOH(L,NY,NX))
    call ncd_putvar(ncf,'CCACO',CCACO(L,NY,NX))
    call ncd_putvar(ncf,'CCASO',CCASO(L,NY,NX))
    call ncd_putvar(ncf,'BKDS',BKDS(L,NY,NX))
    call ncd_putvar(ncf,'FOSRH',FOSRH(:,L,NY,NX))
    call ncd_putvar(ncf,'OQC',OQC(:,L,NY,NX))
    call ncd_putvar(ncf,'OQN',OQN(:,L,NY,NX))
    call ncd_putvar(ncf,'OQP',OQP(:,L,NY,NX))
    call ncd_putvar(ncf,'OQA',OQA(:,L,NY,NX))
    call ncd_putvar(ncf,'OHC',OHC(:,L,NY,NX))
    call ncd_putvar(ncf,'OHN',OHN(:,L,NY,NX))
    call ncd_putvar(ncf,'OHP',OHP(:,L,NY,NX))

    call ncd_putvar(ncf,'CNOSC',CNOSC(:,:,L,NY,NX))
    call ncd_putvar(ncf,'CPOSC',CPOSC(:,:,L,NY,NX))

    call ncd_putvar(ncf,'OSC',OSC(:,:,L,NY,NX))
    call ncd_putvar(ncf,'OSA',OSA(:,:,L,NY,NX))
    call ncd_putvar(ncf,'OSN',OSN(:,:,L,NY,NX))
    call ncd_putvar(ncf,'OSP',OSP(:,:,L,NY,NX))
    call ncd_putvar(ncf,'ORC',ORC(:,:,L,NY,NX))
    call ncd_putvar(ncf,'ORN',ORN(:,:,L,NY,NX))
    call ncd_putvar(ncf,'ORP',ORP(:,:,L,NY,NX))

    call ncd_putvar(ncf,'OMC',OMC(:,:,:,:,L,NY,NX))
    call ncd_putvar(ncf,'OMN',OMN(:,:,:,:,L,NY,NX))
    call ncd_putvar(ncf,'OMP',OMP(:,:,:,:,L,NY,NX))

    call ncd_putvar(ncf,'OMCff',OMCff(:,:,:,L,NY,NX))
    call ncd_putvar(ncf,'OMNff',OMNff(:,:,:,L,NY,NX))
    call ncd_putvar(ncf,'OMPff',OMPff(:,:,:,L,NY,NX))

    if(bgc_forc_conf%laddband)then
      call ncd_putvar(ncf,'ZNH4S',ZNH4S(L,NY,NX)+ZNH4B(L,NY,NX))
      call ncd_putvar(ncf,'ZNO3S',ZNO3S(L,NY,NX)+ZNO3B(L,NY,NX))
      call ncd_putvar(ncf,'ZNO2S',ZNO2S(L,NY,NX)+ZNO2B(L,NY,NX))
      call ncd_putvar(ncf,'H2PO4',H2PO4(L,NY,NX)+H2POB(L,NY,NX))
      call ncd_putvar(ncf,'H1PO4',H1PO4(L,NY,NX)+H1POB(L,NY,NX))
    else
      call ncd_putvar(ncf,'ZNH4S',ZNH4S(L,NY,NX))
      call ncd_putvar(ncf,'ZNO3S',ZNO3S(L,NY,NX))
      call ncd_putvar(ncf,'ZNO2S',ZNO2S(L,NY,NX))
      call ncd_putvar(ncf,'H2PO4',H2PO4(L,NY,NX))
      call ncd_putvar(ncf,'H1PO4',H1PO4(L,NY,NX))
    endif
    call ncd_pio_closefile(ncf)
  endif
  end subroutine WriteBBGCForc

end module ForcWriterMod

module MicIDMod
use ChemIDMod
implicit none
  save
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  integer :: cidg_CO2      !gaseous tracer CO2
  integer :: cidg_CH4      !gaseous tracer CH4
  integer :: cid_OXYG      !gaseous tracer O2
  integer :: cid_Z2GG      !gaseous tracer N2
  integer :: cid_Z2OG      !gaseous tracer N2O
  integer :: cid_ZN3G      !gaseous tracer NH3
  integer :: cid_H2GG      !gaseous tracer H2
  integer :: cid_Z2GS      !aqueous N2
  integer :: cid_ZNH3G     !gaseous NH3 IN soil, [gN d-2]
  integer :: cid_ZNH3B     !NH3 band micropore, [gN d-2]
  integer :: cid_ZNH3S     !NH3 non-band micropore, [gN d-2]
  integer :: cid_ZNH4B     !NH4 band micropore, [g d-2]
  integer :: cid_ZNH4S     !NH4 non-band micropore, [g d-2]
  integer :: cid_H1POB     !soil aqueous HPO4 content micropore band, [g d-2]
  integer :: cid_H1PO4     !soil aqueous HPO4 content micropore non-band, [g d-2]
  integer :: cid_ZNO2B     !NO2(-)  band micropore, [g d-2]
  integer :: cid_ZNO2S     !NO2(-)  non-band micropore, [g d-2]
  integer :: cid_H2POB     !H2PO4 band micropore, [g d-2]
  integer :: cid_H2PO4     !H2PO4 non-band micropore, [g d-2]
  integer :: cid_CCO2S     !aqueous CO2 concentration micropore	[g m-3]
  integer :: cid_CNO2S     !NO2(-) concentration non-band micropore	[g m-3]
  integer :: cid_CNO2B     !NO2 concentration band micropore	[g m-3]
  integer :: cid_CZ2OS     !aqueous N2O concentration micropore	[g m-3]
  integer :: cid_Z2OS      !aqueous N2O mass in micropore, [g d-2]
  integer :: cid_COXYS     !aqueous O2 concentration micropore	[g m-3]
  integer :: cid_OXYS      !aqueous O2 mass in micropore	[g m-2]
  integer :: cid_COXYG     !gaseous O2 concentration	[g m-3]
  integer :: cid_CZ2GG     !gaseous N2 concentration [g m-3]
  integer :: cid_CZ2GS     !aqueous N2 concentration [g m-3]
  integer :: cid_CH2GS     !aqueous H2 concentration [g m-3]
  integer :: cid_H2GS      !aqueous H2 mass [g m-3]
  integer :: cid_CCH4G     !gaseous CH4 concentration in micropore [g m-3]
  integer :: cid_CH4S      !aqueous CH4 mass in  micropore	[g d-2]
  integer :: cid_ZNFN0     !initial nitrification inhibition activity
  integer :: cid_ZNFNI     !current nitrification inhibition activity
  integer :: cid_oqc_b,cid_oqc_e     !dissolved mass organic C micropore	[gC d-2]
  integer :: cid_oqn_b,cid_oqn_e     !dissolved mass organic N micropore	[gN d-2]
  integer :: cid_oqp_b,cid_oqp_e     !dissolved mass organic P micropore	[gP d-2]
  integer :: cid_oqa_b,cid_oqa_e     !dissolved mass acetate micropore [gC d-2]
  integer :: cid_ohc_b,cid_ohc_e     !adsorbed mass soil C	[gC d-2]
  integer :: cid_ohn_b,cid_ohn_e     !adsorbed mass soil N	[gN d-2]
  integer :: cid_ohp_b,cid_ohp_e     !adsorbed mass soil P	[gP d-2]
  integer :: cid_oha_b,cid_oha_e     !adsorbed mass soil acetate	[gC d-2]
  integer :: cid_osc_b,cid_osc_e     !humus mass soil C	[gC d-2]
  integer :: cid_osa_b,cid_osa_e     !humus mass soil acetate	[gC d-2]
  integer :: cid_osn_b,cid_osn_e     !humus mass soil N	[gN d-2]
  integer :: cid_osp_b,cid_osp_e     !humus mass soil P	[gP d-2]
  integer :: cid_orc_b,cid_orc_e     !microbial mass residue C [gC d-2]
  integer :: cid_orn_b,cid_orn_e     !microbial mass residue N [gN d-2]
  integer :: cid_orp_b,cid_orp_e     !microbial mass residue P [gP d-2]
  integer :: cid_omc_b,cid_omc_e     !microbial biomass C	[gC d-2]
  integer :: cid_omn_b,cid_omn_e     !microbial biomass N	[gN d-2]
  integer :: cid_omp_b,cid_omp_e     !microbial biomass P	[gP d-2]
  integer :: cid_omcff_b,cid_omcff_e   !autotrophic microbial biomass C	[gC d-2]
  integer :: cid_omnff_b,cid_omnff_e   !autotrophic microbial biomass N	[gN d-2]
  integer :: cid_ompff_b,cid_ompff_e   !autotrophic microbial biomass C	[gC d-2]

  integer :: fid_XCODFS             !CO2 dissolution (+)-volatiziation (-) with respect to atmosphere
  integer :: fid_XCHDFS             !CH4 dissolution (+)-volatiziation (-) with respect to atmosphere
  integer :: fid_XOXDFS             !O2  dissolution (+)-volatiziation (-) with respect to atmosphere
  integer :: fid_XNGDFS             !N2  dissolution (+)-volatiziation (-) with respect to atmosphere
  integer :: fid_XN2DFS             !N2O dissolution (+)-volatiziation (-) with respect to atmosphere
  integer :: fid_XN3DFS             !NH3 dissolution (+)-volatiziation (-) with respect to atmosphere in non-band soil
  integer :: fid_XNBDFS             !NH3 dissolution (+)-volatiziation (-) with respect to atmosphere in band soil
  integer :: fid_XHGDFS             !H2  dissolution (+)-volatiziation (-) with respect to atmosphere
  integer :: fid_XCODFG             !CO2 dissolution (+)-volatiziation (-) in soil
  integer :: fid_XCHDFG             !CH4 dissolution (+)-volatiziation (-) in soil
  integer :: fid_XOXDFG             !O2 dissolution (+)-volatiziation (-) in soil
  integer :: fid_XNGDFG             !N2 dissolution (+)-volatiziation (-) in soil
  integer :: fid_XN2DFG             !N2O dissolution (+)-volatiziation (-) in soil
  integer :: fid_XN3DFG             !NH3 dissolution (+)-volatiziation (-) in soil
  integer :: fid_XNBDFG             !NH3 dissolution (+)-volatiziation (-) in band soil
  integer :: fid_XHGDFG             !H2 dissolution (+)-volatiziation (-) in soil
  integer :: fid_XCOFLG             !CO2 gaseous exchange with atmosphere (-) into atmosphere
  integer :: fid_XCHFLG             !CH4 gaseous exchange with atmosphere (-) into atmosphere
  integer :: fid_XOXFLG             !O2 gaseous exchange with atmosphere (-) into atmosphere
  integer :: fid_XNGFLG             !N2 gaseous exchange with atmosphere (-) into atmosphere
  integer :: fid_XN2FLG             !N2O gaseous exchange with atmosphere (-) into atmosphere
  integer :: fid_XN3FLG             !N3H gaseous exchange with atmosphere (-) into atmosphere
  integer :: fid_XHGFLG             !H2 gaseous exchange with atmosphere (-) into atmosphere


  integer :: fid_ROXYF               !net gaseous O2 flux from previous hour, [g d-2 h-1]
  integer :: fid_ROXYY               !total root + microbial O2 uptake, [g d-2 h-1]
  integer :: fid_RNH4Y               !total root + microbial NH4(+) uptake non-band, [gN d-2 h-1]
  integer :: fid_RNO3Y               !total root + microbial NO3(-) uptake non-band, [gN d-2 h-1]
  integer :: fid_RNO2Y               !total root + microbial NO2(-) uptake non-band, [gN d-2 h-1]
  integer :: fid_RN2OY               !total root + microbial N2O uptake, [g d-2 h-1]
  integer :: fid_RPO4Y               !total root + microbial PO4 uptake non-band, [gP d-2 h-1]
  integer :: fid_RP14Y               !HPO4 demand in non-band by all microbial,root,myco populations [gP d-2 h-1]
  integer :: fid_RNHBY               !total root + microbial NH4 uptake band, [gN d-2 h-1]
  integer :: fid_RN3BY               !total root + microbial NO3(-) uptake band, [gN d-2 h-1]
  integer :: fid_RN2BY               !total root + microbial NO2(-) uptake band, [gN d-2 h-1]
  integer :: fid_RPOBY               !total root + microbial PO4 uptake band, [gP d-2 h-1]
  integer :: fid_RP1BY               !HPO4 demand in band by all microbial,root,myco populations [gP d-2 h-1]
  integer :: fid_ROQCY_b,fid_ROQCY_e !total root + microbial DOC uptake, [gC d-2 h-1]
  integer :: fid_ROQAY_b,fid_ROQAY_e !total root + microbial acetate uptake, [gC d-2 h-1]
  integer :: fid_RINHO_b,fid_RINHO_e
  integer :: fid_RINHB_b,fid_RINHB_e
  integer :: fid_RINOO_b,fid_RINOO_e
  integer :: fid_RINOB_b,fid_RINOB_e
  integer :: fid_RIPOO_b,fid_RIPOO_e
  integer :: fid_RIPBO_b,fid_RIPBO_e
  integer :: fid_RIPO1_b,fid_RIPO1_e
  integer :: fid_RIPB1_b,fid_RIPB1_e
  integer :: fid_ROXYS_b,fid_ROXYS_e
end module MicIDMod

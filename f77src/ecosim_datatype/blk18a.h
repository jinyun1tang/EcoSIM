      real(r8) :: ARLFC(JY,JX),ARSTC(JY,JX),TEVAPP(JY,JX),TEVAPC(JY,JX)
     2,TENGYC(JY,JX),THFLXC(JY,JX),TUPWTR(0:JZ,JY,JX),TUPHT(0:JZ,JY,JX)
     3,TVOLWP(JY,JX),TCOFLA(JZ,JY,JX),TOXFLA(JZ,JY,JX),TCHFLA(JZ,JY,JX)
     4,TN2FLA(JZ,JY,JX),TNHFLA(JZ,JY,JX),TLCO2P(JZ,JY,JX),GPP(JY,JX)
     5,TLOXYP(JZ,JY,JX),TLCH4P(JZ,JY,JX),TLN2OP(JZ,JY,JX),RECO(JY,JX)
     6,TLNH3P(JZ,JY,JX),TCO2S(JZ,JY,JX),TUPOXS(JZ,JY,JX)
     7,TUPCHS(JZ,JY,JX),TUPN2S(JZ,JY,JX),TUPN3S(JZ,JY,JX)
     8,TUPNH4(JZ,JY,JX),TUPNO3(JZ,JY,JX),TUPH2P(JZ,JY,JX)
     9,TUPN3B(JZ,JY,JX),TUPNHB(JZ,JY,JX),TUPNOB(JZ,JY,JX)
     1,TUPH2B(JZ,JY,JX),TUPNF(JZ,JY,JX),CSNT(4,0:1,0:JZ,JY,JX)
     2,ZSNT(4,0:1,0:JZ,JY,JX),PSNT(4,0:1,0:JZ,JY,JX)
     3,TDFOMC(0:4,JZ,JY,JX),TDFOMN(0:4,JZ,JY,JX),TDFOMP(0:4,JZ,JY,JX)
     4,TCO2Z(JY,JX),TOXYZ(JY,JX),TCH4Z(JY,JX),TN2OZ(JY,JX),TNH3Z(JY,JX)
     5,TCO2P(JZ,JY,JX),TUPOXP(JZ,JY,JX),THRMC(JY,JX),TCNET(JY,JX)
     6,ZCSNC(JY,JX),ZZSNC(JY,JX),ZPSNC(JY,JX),WGLFT(JC,JY,JX)
     7,ARLFT(JC,JY,JX),ARSTT(JC,JY,JX),ARLSS(JY,JX),RTDNT(JZ,JY,JX)
     6,TCCAN(JY,JX),TUPH1P(JZ,JY,JX),TUPH1B(JZ,JY,JX)

      COMMON/BLK18/ARLFC,ARSTC,TEVAPP,TEVAPC
     2,TENGYC,THFLXC,TUPWTR,TUPHT
     3,TVOLWP,TCOFLA,TOXFLA,TCHFLA
     4,TN2FLA,TNHFLA,TLCO2P,GPP
     5,TLOXYP,TLCH4P,TLN2OP,RECO
     6,TLNH3P,TCO2S,TUPOXS
     7,TUPCHS,TUPN2S,TUPN3S
     8,TUPNH4,TUPNO3,TUPH2P
     9,TUPN3B,TUPNHB,TUPNOB
     1,TUPH2B,TUPNF,CSNT
     2,ZSNT,PSNT
     3,TDFOMC,TDFOMN,TDFOMP
     4,TCO2Z,TOXYZ,TCH4Z,TN2OZ,TNH3Z
     5,TCO2P,TUPOXP,THRMC,TCNET
     6,ZCSNC,ZZSNC,ZPSNC,WGLFT
     7,ARLFT,ARSTT,ARLSS,RTDNT
     6,TCCAN,TUPH1P,TUPH1B
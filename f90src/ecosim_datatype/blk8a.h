
      real(r8) :: PSIFC(JY,JX),PSIWP(JY,JX),FMPR(0:JZ,JY,JX) &
      ,ALBS(JY,JX),CDPTH(0:JZ,JY,JX),ROCK(JZ,JY,JX),TSED(JY,JX) &
      ,BKDS(0:JZ,JY,JX),FC(0:JZ,JY,JX),WP(0:JZ,JY,JX),SCNV(0:JZ,JY,JX) &
      ,SCNH(JZ,JY,JX),CSAND(JZ,JY,JX),CSILT(JZ,JY,JX),CCLAY(JZ,JY,JX) &
      ,FHOL(JZ,JY,JX),PHOL(JZ,JY,JX),DHOL(JZ,JY,JX),HRAD(JZ,JY,JX) &
      ,BKDSI(JZ,JY,JX),PH(0:JZ,JY,JX),CEC(JZ,JY,JX),AEC(JZ,JY,JX) &
      ,CORGC(0:JZ,JY,JX),CORGR(JZ,JY,JX),CORGN(JZ,JY,JX),CORGP(JZ,JY,JX) &
      ,CNH4(JZ,JY,JX),CNO3(JZ,JY,JX),CPO4(JZ,JY,JX),CAL(JZ,JY,JX) &
      ,CFE(JZ,JY,JX),CCA(JZ,JY,JX),CMG(JZ,JY,JX),CNA(JZ,JY,JX) &
      ,CKA(JZ,JY,JX),CSO4(JZ,JY,JX),CCL(JZ,JY,JX),CALOH(JZ,JY,JX) &
      ,CFEOH(JZ,JY,JX),CCACO(JZ,JY,JX),CCASO(JZ,JY,JX),CALPO(JZ,JY,JX) &
      ,CFEPO(JZ,JY,JX),CCAPD(JZ,JY,JX),CCAPH(JZ,JY,JX),SED(JY,JX) &
      ,GKC4(JZ,JY,JX),GKCH(JZ,JY,JX),GKCA(JZ,JY,JX),GKCM(JZ,JY,JX) &
      ,GKCN(JZ,JY,JX),GKCK(JZ,JY,JX),THW(JZ,JY,JX),THI(JZ,JY,JX) &
      ,RSC(0:2,0:JZ,JY,JX),RSN(0:2,0:JZ,JY,JX),RSP(0:2,0:JZ,JY,JX) &
      ,CNOFC(4,0:2),CPOFC(4,0:2),DETS(JY,JX),DETE(JY,JX),CER(JY,JX) &
      ,XER(JY,JX),SLOPE(0:3,JY,JX),CDPTHI(JY,JX),CORGCI(JZ,JY,JX) &
      ,POROSI(JZ,JY,JX),FHOLI(JZ,JY,JX),PTDSNU(JY,JX),VLS(JY,JX)
      integer :: NU(JY,JX),NUI(JY,JX),NJ(JY,JX),NK(JY,JX) &
      ,NLI(JV,JH),NL(JV,JH),ISOILR(JY,JX),NHOL(JZ,JY,JX),NUM(JY,JX)


      COMMON/BLK8A/PSIFC,PSIWP,FMPR &
      ,ALBS,CDPTH,ROCK,TSED &
      ,BKDS,FC,WP,SCNV &
      ,SCNH,CSAND,CSILT,CCLAY &
      ,FHOL,PHOL,DHOL,HRAD &
      ,BKDSI,PH,CEC,AEC &
      ,CORGC,CORGR,CORGN,CORGP &
      ,CNH4,CNO3,CPO4,CAL &
      ,CFE,CCA,CMG,CNA &
      ,CKA,CSO4,CCL,CALOH &
      ,CFEOH,CCACO,CCASO,CALPO &
      ,CFEPO,CCAPD,CCAPH,SED &
      ,GKC4,GKCH,GKCA,GKCM &
      ,GKCN,GKCK,THW,THI &
      ,RSC,RSN,RSP &
      ,CNOFC,CPOFC,DETS,DETE,CER &
      ,XER,SLOPE,CDPTHI,CORGCI &
      ,POROSI,FHOLI,PTDSNU,VLS &
      ,NU,NUI,NJ,NK,NLI,NL &
      ,ISOILR,NHOL,NUM
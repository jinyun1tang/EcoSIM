
      real(r8) :: PHR(JY,JX),CN4RI(JY,JX),CNORI(JY,JX),CN4R(JY,JX) &
      ,CN3R(JY,JX),CNOR(JY,JX),CPOR(JY,JX),CALR(JY,JX),CFER(JY,JX) &
      ,CHYR(JY,JX),CCAR(JY,JX),CMGR(JY,JX),CNAR(JY,JX),CKAR(JY,JX) &
      ,COHR(JY,JX),CSOR(JY,JX),CCLR(JY,JX),CC3R(JY,JX),CHCR(JY,JX) &
      ,CCHR(JY,JX),CAL1R(JY,JX),CAL2R(JY,JX),CAL3R(JY,JX),CAL4R(JY,JX) &
      ,CALSR(JY,JX),CFE1R(JY,JX),CFE2R(JY,JX),CFE3R(JY,JX),CFE4R(JY,JX) &
      ,CFESR(JY,JX),CCAOR(JY,JX),CCACR(JY,JX),CCAHR(JY,JX),CCASR(JY,JX) &
      ,CMGOR(JY,JX),CMGCR(JY,JX),CMGHR(JY,JX),CMGSR(JY,JX),CNACR(JY,JX) &
      ,CNASR(JY,JX),CKASR(JY,JX),CH0PR(JY,JX),CH1PR(JY,JX),CH3PR(JY,JX) &
      ,CF1PR(JY,JX),CF2PR(JY,JX),CC0PR(JY,JX),CC1PR(JY,JX),CC2PR(JY,JX) &
      ,CM1PR(JY,JX),CCOR(JY,JX),COXR(JY,JX),CNNR(JY,JX),CN2R(JY,JX) &
      ,PHQ(366,JY,JX),CN4Q(366,JY,JX),CN3Q(366,JY,JX),CNOQ(366,JY,JX) &
      ,CPOQ(366,JY,JX),CALQ(366,JY,JX),CFEQ(366,JY,JX),CHYQ(366,JY,JX) &
      ,CCAQ(366,JY,JX),CMGQ(366,JY,JX),CNAQ(366,JY,JX),CKAQ(366,JY,JX) &
      ,COHQ(366,JY,JX),CSOQ(366,JY,JX),CCLQ(366,JY,JX),CC3Q(366,JY,JX) &
      ,CHCQ(366,JY,JX),CAL1Q(366,JY,JX),CAL2Q(366,JY,JX) &
      ,CAL3Q(366,JY,JX),CAL4Q(366,JY,JX),CALSQ(366,JY,JX) &
      ,CFE1Q(366,JY,JX),CFE2Q(366,JY,JX),CFE3Q(366,JY,JX) &
      ,CFE4Q(366,JY,JX),CFESQ(366,JY,JX),CCAOQ(366,JY,JX) &
      ,CCACQ(366,JY,JX),CCAHQ(366,JY,JX),CCASQ(366,JY,JX) &
      ,CMGOQ(366,JY,JX),CMGCQ(366,JY,JX),CMGHQ(366,JY,JX) &
      ,CMGSQ(366,JY,JX),CNACQ(366,JY,JX),CNASQ(366,JY,JX) &
      ,CKASQ(366,JY,JX),CH0PQ(366,JY,JX),CH1PQ(366,JY,JX) &
      ,CH3PQ(366,JY,JX),CF1PQ(366,JY,JX),CF2PQ(366,JY,JX) &
      ,CC0PQ(366,JY,JX),CC1PQ(366,JY,JX),CC2PQ(366,JY,JX) &
      ,CM1PQ(366,JY,JX),CSTRQ(366,JY,JX),DENS0(JY,JX) &
      ,CSTRR(JY,JX)


      COMMON/BLK2B/PHR,CN4RI,CNORI,CN4R &
      ,CN3R,CNOR,CPOR,CALR,CFER &
      ,CHYR,CCAR,CMGR,CNAR,CKAR &
      ,COHR,CSOR,CCLR,CC3R,CHCR &
      ,CCHR,CAL1R,CAL2R,CAL3R,CAL4R &
      ,CALSR,CFE1R,CFE2R,CFE3R,CFE4R &
      ,CFESR,CCAOR,CCACR,CCAHR,CCASR &
      ,CMGOR,CMGCR,CMGHR,CMGSR,CNACR &
      ,CNASR,CKASR,CH0PR,CH1PR,CH3PR &
      ,CF1PR,CF2PR,CC0PR,CC1PR,CC2PR &
      ,CM1PR,CCOR,COXR,CNNR,CN2R &
      ,PHQ,CN4Q,CN3Q,CNOQ &
      ,CPOQ,CALQ,CFEQ,CHYQ &
      ,CCAQ,CMGQ,CNAQ,CKAQ &
      ,COHQ,CSOQ,CCLQ,CC3Q &
      ,CHCQ,CAL1Q,CAL2Q &
      ,CAL3Q,CAL4Q,CALSQ &
      ,CFE1Q,CFE2Q,CFE3Q &
      ,CFE4Q,CFESQ,CCAOQ &
      ,CCACQ,CCAHQ,CCASQ &
      ,CMGOQ,CMGCQ,CMGHQ &
      ,CMGSQ,CNACQ,CNASQ &
      ,CKASQ,CH0PQ,CH1PQ &
      ,CH3PQ,CF1PQ,CF2PQ &
      ,CC0PQ,CC1PQ,CC2PQ &
      ,CM1PQ,CSTRQ,DENS0 &
      ,CSTRR
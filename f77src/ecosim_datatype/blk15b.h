      real(r8) :: XCOFHS(3,JD,JV,JH),XCHFHS(3,JD,JV,JH)
     2,XOXFHS(3,JD,JV,JH),XNGFHS(3,JD,JV,JH),XN2FHS(3,JD,JV,JH)
     3,XHGFHS(3,JD,JV,JH),XN4FHW(3,JD,JV,JH),XN3FHW(3,JD,JV,JH)
     4,XNOFHW(3,JD,JV,JH),XH2PHS(3,JD,JV,JH),XNXFHS(3,JD,JV,JH)
     5,XN4FHB(3,JD,JV,JH),XN3FHB(3,JD,JV,JH),XNOFHB(3,JD,JV,JH)
     6,XNXFHB(3,0:JD,JV,JH),XH2BHB(3,JD,JV,JH),XOCFHS(0:4,3,JD,JV,JH)
     7,XONFHS(0:4,3,JD,JV,JH),XOPFHS(0:4,3,JD,JV,JH)
     8,XOAFHS(0:4,3,JD,JV,JH),FLWX(3,JD,JV,JH)
     9,XH1PHS(3,JD,JV,JH),XH1BHB(3,JD,JV,JH)

      COMMON/BLK15B/XCOFHS,XCHFHS
     2,XOXFHS,XNGFHS,XN2FHS
     3,XHGFHS,XN4FHW,XN3FHW
     4,XNOFHW,XH2PHS,XNXFHS
     5,XN4FHB,XN3FHB,XNOFHB
     6,XNXFHB,XH2BHB,XOCFHS
     7,XONFHS,XOPFHS
     8,XOAFHS,FLWX
     9,XH1PHS,XH1BHB
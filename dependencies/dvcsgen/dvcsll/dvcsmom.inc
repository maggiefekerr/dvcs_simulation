c
c    all moments
c
c
c
       real hc0BH,hc1BH,hc2BH, 
     + hc0dvcs,hc1dvcs,hs1dvcs,hs2dvcs,
     + hc0Iunp,hc1Iunp,hs1Iunp,hc2Iunp,hs2Iunp
       COMMON /MOMUNP/ hc0BH,hc1BH,hc2BH,
     + hc0dvcs,hc1dvcs,hs1dvcs,hs2dvcs
       COMMON /MOMUNPI/hc0Iunp,hc1Iunp,hs1Iunp,hc2Iunp,hs2Iunp
c
c
c lpol
c
c234567890123456789012345678901234567890123456789012345678901234567890123456
c
      real hc0BHlp,hc1BHlp,hs1Ilp,hc1Ilp,hc0Ilp
      common/momentslp/hc0BHlp,hc1BHlp,hc0Ilp,hs1Ilp,hc1Ilp
c
c
c   t-pol
c
c
      real hc0BHtpcos,hc1BHtpcos,hs1BHtpsin
      common/momentstpBH/hc0BHtpcos,hc1BHtpcos,hs1BHtpsin
c
      real hs1Itpcos,hs1Itpsin,hc1Itpcos
     6               ,hc1Itpsin, hc0Itpcos,hc0Itpsin
      common/momentstpI/hc0Itpcos,hc0Itpsin,hs1Itpcos
     6               ,hs1Itpsin,hc1Itpcos,hc1Itpsin







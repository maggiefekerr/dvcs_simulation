*
*      REAL elee,elet,elef,elex,eley,elez,xb,q2,w
*      INTEGER nele,eleevtclass,eleh,elerun,eletype,eleevent
*      COMMON /ELEREC/ nele,eleevtclass,eleh,elerun,eletype,eleevent,elee
*     + ,elet,elef,elex,eley,elez,xb,q2,w
*
*      REAL pi0e,pi0t,pi0f,pipe,pipt,pipf,pime,pimt,pimf,proe,prot,prof
*     + ,game,gamt,gamf,gam2t,gam2e,gam2f,neue,neut,neuf,pim,pipx,pipy
*     + ,pipz,pimx,pimy,pimz,prox,proy,proz
*      INTEGER npi0,npip,npim,npro,ngam,nneu
*      COMMON /HADREC/ npi0,npip,npim,npro,ngam,nneu,pi0e,pi0t,pi0f,pipe
*     + ,pipt,pipf,pime,pimt,pimf,proe,prot,prof,game,gamt,gamf,gam2t
*     + ,gam2e,gam2f,neue,neut,neuf,pim,pipx,pipy,pipz,pimx,pimy,pimz
*     + ,prox,proy,proz
*
      REAL gelee,gelet,gelef,gelex,geley,gelez,gelewgt,gxb,gq2,gw
      INTEGER gnele,gelebits,geleh
      COMMON /ELEGEN/ gnele,gelebits,geleh,gelee,gelet,gelef,gelex,geley
     + ,gelez,gelewgt,gxb,gq2,gw
*
      REAL gpi0e,gpi0t,gpi0f,gproe,gprot,gprof
     + ,ggame,ggamt,ggamf,ggam2t,ggam2e,ggam2f,gprox,gproy,gproz
     + ,gprotg,gprotl,gprotv,gprotmis,gprott,gprofs
     + ,gwbh,gwdvcs,gvint,gsin,gsin2
      INTEGER gnpi0,gnpip,gnpim,gnpro,gngam,gnneu
      COMMON /HADGEN/ gnpi0,gnpro,gngam,gnneu,gpi0e,gpi0t
     + ,gpi0f,gproe,gprot,gprof
     + ,ggame,ggamt,ggamf,ggam2t,ggam2e,ggam2f
     + ,gprox,gproy,gproz
     + ,gprotg,gprotl,gprotv,gprotmis,gprott,gprofs
     + ,gwbh,gwdvcs,gvint,gsin,gsin2
c
       INTEGER gproh
       REAL ghp1,ghp2,ghccb,ghcci,gh0bh,gh1bh,gh2bh,gh0dvcs
     + ,gh1dvcs,ghs1dvcs,gh0iunp,gh1iunp,ghs1iunp,ghs2iunp
     + ,gh0bhlp,gh1bhlp,ghs1ilp,gh1ilp,gh0ilp
     + ,gh0bhtpcos,gh1bhtpcos,ghs1bhtpsin
     + ,ghs1itpcos,ghs1itpsin,gh1itpcos,gh1itpsin,gh0itpcos,gh0itpsin
c
       COMMON /MOM0/ gproh,ghp1,ghp2,ghccb,ghcci,gh0bh,gh1bh,gh2bh
     + ,gh0dvcs,gh1dvcs,ghs1dvcs,gh0iunp,gh1iunp,ghs1iunp,ghs2iunp
     + ,gh0bhlp,gh1bhlp,ghs1ilp,gh1ilp,gh0ilp
     + ,gh0bhtpcos,gh1bhtpcos,ghs1bhtpsin
     + ,ghs1itpcos,ghs1itpsin,gh1itpcos,gh1itpsin,gh0itpcos,gh0itpsin

c
        integer iproctype
        real radxb,radxb_true,radq2,radq2_true,radphi,radt,radnu,bornweight,radweight,grade,gradt,gradf
       COMMON /RADISTUFF/iproctype,radxb,radxb_true,radq2,radq2_true,radphi,radt,radnu,bornweight,radweight,
     + grade,gradt,gradf

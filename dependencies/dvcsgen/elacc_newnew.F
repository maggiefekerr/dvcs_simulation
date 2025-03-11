c
c
c     Peter's fiducial cuts
c
c=====================================================================
c FIDUCIAL CUTS FOR 5.764 GeV outbending data
c Volkers acceptance function modified by Vipuli 24/02/2003
c old cuts icutei=1 or 2 
c new cuts icutei=3
c=====================================================================
      subroutine elacc_newnew(ebeam,
     *     secte,pel,thetaeld,phiel,torcur,qe,
     *     icutei,del_phie,weight_e)

* Version 1.0 / November 30, 1999, modified 01/29/02
* this version is from 4/16/03
* ebeam: beam energy in GeV (real*4)
* secte: sector where electron is detected (real*4)
* pel: electron momentum, after corrections in GeV (real*4)
* thetaeld: scattering angle in degrees at dc layer 1 (real*4)
* phiel: azimuthal angle in deg. at dc layer 1 (ie 0 to 60) (real*4)
* note: the above two are trl1_theta and trl1_phi in dst's.
* torcur: torus current in amps (e.g. -1500.) (real*4) 
* icutei defines tightness of fiducial cuts for inbending electrons
* qe: particle charge: should be -1 or 1 (integer)
* icutei = 1  tight cut (~100% efficiency inside fiducial region)
* icutei = 2  loose cut (includes some of the fringes of the CC)  
* icutei = 3  new cut for eg1b (integer)                          
* del_phie: electron acceptance in phi (degrees) for each event, this
*       value should enter in d(omega) = dcos(thetaeld)*del_phie when 
*       calculating a cross section   (real*4)
* weight_e: =1. if electron is accepted,  0. if not accepted (real*4)
* Modifications make by P. Bosted:
* added agruement qe to subroutine
* use torcur * -q to check for in/out bending
* and use abs(torcur) everywhere except to check for in/out
*   instead of just torcur, which can be negative if positrons.
* also added implicit none, defined all variables, lined up If/endifs
* also had to add check on theta>thcut in places in inbending
*      where a sin(theta-thcut)**expon would give an error
*      (you can't exponentiate a negative number)

      implicit none
      real ebeam, secte, pel,thetaeld,phiel, phie_mod, torcur,del_phie 
      real thcut, pnorm,midphi, expon
      integer qe,icutei,weight_e
      
      phie_mod = phiel

* Define acceptance cuts for electrons
      weight_e=0
        
c************************************************************
c Inbending fiducial cuts
c************************************************************       

      If(qe*torcur.lt.0.) then
        if ( ebeam.gt.5.5 ) then
          weight_e = 0
          
          if(icutei.eq.2) then
            thcut =11.0 + 16.72/((pel+0.050)*3375./abs(torcur))
            expon = 0.28*(pel*3375./abs(torcur))**0.25
            if(thetaeld.gt.thcut) 
     >        del_phie = 25.*sin((thetaeld-thcut)*
     >          3.14159/180.)**expon
          endif
          if(icutei.eq.1) then
            thcut =10.0 + 17.5/((pel+0.050)*3375./abs(torcur))
            expon = 0.33*(pel*3375./abs(torcur))**0.33
            if(thetaeld.gt.thcut) 
     >         del_phie = 25.*sin((thetaeld-thcut)*
     >            3.14159/180.)**expon
          endif
c--------------------------------------------------------------
c     combined fiducial cuts -new
c----------------------------------------------------------------
          if(icutei.eq.3) then
            if(pel.le.3.) then
              thcut =11.0 + 16.72/((pel+0.050)*3375./abs(torcur))
              if(thcut.lt.32.) then
                expon = 0.30*(pel*3375./abs(torcur))**0.30
                if(thetaeld.lt.32.and.thetaeld.gt.thcut) then
                  del_phie = 32.*sin((thetaeld-thcut)*
     >               3.14159/180.)**expon
                elseif(thetaeld.ge.32.0) then
                  del_phie = 32.*sin((32-thcut)*
     >               3.14159/180.)**expon
                endif
              elseif(thcut.ge.32.0.and.thetaeld.gt.thcut) then
                expon = 0.28*(pel*3375./abs(torcur))**0.25
                del_phie = 24.*sin((thetaeld-thcut)*
     >             3.14159/180.)**expon
              endif
            elseif(pel.gt.3) then
              thcut =10.7 + 16.72/((pel+0.050)*3375./abs(torcur))
              expon = 0.22*(pel*3375./abs(torcur))**0.25
              if(thetaeld.gt.thcut) 
     >           del_phie = 25.*sin((thetaeld-thcut)*
     >              3.14159/180.)**expon
            endif
          endif
c-------------------------------------------------------------------
c-------------------------------------------------------------------
c          print *,thetaeld,thcut,phie_mod,del_phie
c
          if(thetaeld.gt.thcut.and.thetaeld.lt.45.) then
            if(phie_mod.lt.(30.+del_phie).and.
     >         phie_mod.gt.(30.-del_phie)) then
              weight_e=1
            else
              weight_e=0
            endif
          else
            del_phie= 0.
            weight_e=0
          endif
          del_phie=2.*del_phie

          if(secte.eq.5.and.icutei.eq.3) then
            if(thetaeld.gt.(17.5+8.5/(pel*3375./abs(torcur))).and.
     *          thetaeld.lt.(21.5+14.5/(pel*3375./abs(torcur)))) then
              del_phie=del_phie/2.
              if(phie_mod.gt.30.) then
                weight_e=0
              endif
            endif
          endif
          
        endif ! end check 1.x or 5.x gev
        
      Endif ! end check inbending or outbending
       
      return
      end
c

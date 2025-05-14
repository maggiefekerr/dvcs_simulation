      integer function accep11(pid,p,t,f,torus)
      implicit none
       real torus 
       real newphi11,e,p,f,t,torcur,pi,phis,thetad,d_phi,weight
       integer pid
       pi=acos(-1.0)
       torcur=torus
       phis=newphi11(f*180.0/pi)
       thetad=t*180.0/pi
c
       call clas_at11g(pid,p,thetad,phis,torcur,d_phi,weight)
c       print *,thetad,gelef,phis,weight
       accep11=weight
       end
*     
C
      real Function newphi11(phi)
      real phi,phinew
      if (phi.gt.330.) then
        phinew = phi-360.
      elseif (phi.ge.0.0.and.phi.le.30.) then
        phinew = phi
      elseif (phi.gt.30.0.and.phi.le.90.) then
        phinew = phi-60.
      elseif (phi.gt.90.0.and.phi.le.150.) then
        phinew = phi-120.
      elseif (phi.gt.150.0.and.phi.le.210.) then
        phinew = phi-180.
      elseif (phi.gt.210.0.and.phi.le.270.) then
        phinew = phi-240.
      elseif (phi.gt.270.0.and.phi.le.330.) then
        phinew = phi-300.
      endif
c       print *,phi,phinew
      newphi11=phinew
      end  


c==================================================
      subroutine clas_at11g(id,p,thetad,phis,torcur,d_phi,weight)
*     
*     Version 1.1 / April 18, 2001
*     Simulates acceptance and smearing for CLAS12
*
*     Nov. 30, 2004 - major changes as I separate functions into
*     separate subroutines - Mestayer
*
* Inputs -
*     id:     Particle ID according to PDG (11-electron, 2212-proton, 
*     211-pi+, -211-pi-, 22-photon ...)
*     p:      momentum, in GeV/c, on output will be smeared
*     phis:   phi in the sector in degrees (at mid plane phi=0.), 
*                 on output will be smeared
*     thetad: scattering angle in degrees, 
*                 on output will be smeared
*     torcur: torus current in Amps (e.g. -1500.)
* Outputs -
*     d_phi:  part of the phi acceptance, full acceptance should be 
*                 delta_phi=2*pi*d_phi
*     weight: =1. if accepted,  0. if not
c     
      implicit none
      real p,thetad,phis,torcur,d_phi,weight
      real pout,thout,phiout
      integer id
c     
      real phie_mod,thcut,pnorm,expon
      logical pcutl
      real ec_r,ec_th,rl,rh
      data ec_r,ec_th,rl,rh/510.32,0.436332,562.22,215.67/
      real ylow,yhi,tgrho,dl,rc0
      data ylow,yhi,tgrho,dl,rc0/-182.974,189.956,1.95325,100.,120./
      real r,xh,yh,zh,xcut
      real uh,vh,wh,xi,yi
      real sinth,costh,sinph,cosph,dth
      real rcoil,hcoil,dcoil
c     
      real dthe,dphi,tmax
      data dthe,dphi,tmax/0.0286,0.0573,3375./
      real d2r,pi,r2d,dgap,rran(3)
      data d2r,pi,r2d,dgap/0.01745,3.1415926,57.299578,2./
      integer lrun
      data lrun/3/
      real charge
c  common block keeping the logical values of the track type;      
      logical inbend,outbend,neutral,fwd,central
      common/track_type/inbend,outbend,neutral,fwd,central
c
c  initialize weight and d_phi     
      weight=0.
      d_phi=0.
c
c  1st decide if particle charge is +1, -1, 0 or unknown
c
      call charge_type(id,charge)
c      print *,'id, charge',id,charge
c  a charge type of 999. is an unknown
      if(charge.eq.999.) then
         print *,'charge is 999. - unknown particle'
         return
      endif
c
c  2nd, determine which track type we are concerned with;
c  the subroutine tracktype fills the logical variables in common
c  block, track_type
c
      call tracktype(charge,torcur,thetad)
c      print *,'in,out,neu,fwd,cent',inbend,outbend,neutral,fwd,central
c
c  3rd, call acceptance function with charge, torus current,
c   and p,th,phi as inputs; return is weight (1, accepted; 0, not)
c   and d_phi which is the phi part of the acceptance.  Note, this
c   is the tracking acceptance only, it does NOT contain the additional
c   fiducial volume cuts of Cerenkov counters imposed on electrons.
c
c  get a 3-dim. random number
      call myrnorml(rran,lrun)
c
      if(abs(charge).eq.1.or.charge.eq.0)then
         call accep_fun(rran,charge,torcur,p,thetad,phis,weight,d_phi)
c         print *,'accep - weight,dphi',weight,d_phi
      endif
c
c  4th, if particle is accepted (weight=1.) call smearing routine
c   with input momenta and output momenta
c
      if(weight.eq.1.) then
      call myrnorml(rran,lrun)
c
         call smear_fun(rran,charge,torcur,p,thetad,phis,pout,thout,
     $        phiout)
c  now set p,th,phi to smeared values
         p=pout
         thetad=thout
         phis=phiout
      endif
c     
      return
      end
C=======================================================================
      Subroutine ec_xyz_duvw(x,y,z,u,v,w,xi,yi)
      implicit none
      REAL x,y,z,u,v,w,xi,yi,zi
      REAL EC_the,ec_phi,phi,ylow,yhi,tgrho,sinrho,cosrho
      data EC_the/0.4363323/
      data ylow,yhi/-182.974,189.956/
      data tgrho,sinrho,cosrho/1.95325,0.8901256,0.455715/
      real rot(3,3)
c
      phi=atan2(y,x)*57.29578
      if(phi.lt.0.)phi=phi+360.
      phi=phi+30.
      if(phi.ge.360.)phi=phi-360. 	
      Ec_phi=int(phi/60.)*1.0471975
c     
      rot(1,1)=cos(Ec_the)*cos(Ec_phi)
      rot(1,2)=-sin(Ec_phi)
      rot(1,3)=sin(Ec_the)*cos(Ec_phi)
      rot(2,1)=cos(Ec_the)*sin(Ec_phi)
      rot(2,2)=cos(Ec_phi)
      rot(2,3)=sin(Ec_the)*sin(Ec_phi)
      rot(3,1)=-sin(Ec_the)
      rot(3,2)=0.
      rot(3,3)=cos(Ec_the)
c     
      yi=x*rot(1,1)+y*rot(2,1)+z*rot(3,1)
      xi=x*rot(1,2)+y*rot(2,2)+z*rot(3,2)
      zi=x*rot(1,3)+y*rot(2,3)+z*rot(3,3)
      zi=zi-510.32
      u=(yi-ylow)/sinrho
      v=(yhi-ylow)/tgrho-xi+(yhi-yi)/tgrho
      w=((yhi-ylow)/tgrho+xi+(yhi-yi)/tgrho)/2./cosrho
      end
C++++++++++++++++++++++++++++++++++++++++
      subroutine charge_type(id,charge)
c subroutine which returns charge of a particle given its PDG id
      implicit none
      real charge
      integer id
c
      if(id.eq.11.or.id.eq.-211)then
         charge=-1.
      elseif(id.eq.2212.or.id.eq.211.or.id.eq.45)then
         charge=1.
      elseif(id.eq.22)then
         charge=0.
      else   
         charge=999.
      endif   
      end
c+++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++
      subroutine tracktype(charge,torcur,thetad)
c  subroutine which figures out the track type
c  common block keeping the logical values of the track type;      
      implicit none
      real charge,torcur,thetad
      logical inbend,outbend,neutral,fwd,central
      common/track_type/inbend,outbend,neutral,fwd,central
c  initialize logicals to false
      inbend=.false.
      outbend=.false.
      neutral=.false.
      fwd=.false.
      central=.false.
c  6 functions: (inbend,outbend,neutral) track X (fwd,central) det.
c
      if(charge.eq.0) neutral=.true.
      if(abs(charge).eq.1)then
         if(charge*torcur.gt.0.)outbend=.true.
         if(charge*torcur.lt.0.)inbend=.true.
      endif
      if(thetad.lt.35.)fwd=.true.
      if(thetad.ge.35.)central=.true.
c
      end
c+++++++++++++++++++++++++++++++++++++++++++++++
      subroutine accep_fun(ran,charge,torcur,p,thetad,phis,weight,
     $     d_phi)
      implicit none
      real charge,torcur,p,thetad,phis,weight,d_phi
      real pnorm,tmax,thcut,expon,d2r,sinth,costh
      real sinph,cosph,r,ec_r,ec_th,dl,xh,yh,zh,uh,vh,wh
      real xi,yi,rcoil,rc0,hcoil,dcoil,ylow,yhi,xcut,tgrho,dth
      real r2d,ran(3)
      logical inbend,outbend,neutral,fwd,central
      logical pcutl
      data tmax/3375./
      data d2r,r2d/0.01745,57.299578/
c common block containing track_type; defined in sub. tracktype
      common/track_type/inbend,outbend,neutral,fwd,central
c
c  since our acceptance functions care about sagitta and not momentum, per se,
c  we will form a "normalized momentum", pnorm
c
      pnorm=p*tmax/abs(torcur)
c
c  initialize values
      weight=0.
      d_phi=0.
c  select which of the 6 functions we want
      if(inbend.and.fwd)then
         thcut =5. + 17./pnorm
         if(thetad.gt.thcut)then
            expon = 0.25*(pnorm)**0.2
            d_phi = 35.*sin((thetad-thcut)*d2r)**expon 
c            if(id.eq.-211)d_phi=d_phi*1.1
            if(abs(phis).lt.d_phi)then
               weight=1.
            else
               weight=0.
            endif
         else
            d_phi= 0.
            weight=0.
         endif
      elseif(inbend.and.central)then
         d_phi=30.
         weight=1.
c
      elseif(outbend.and.fwd)then
c         thcut = 5.0  - 25.*(1.-pnorm*0.089)**24
         thcut=5.0
c  don't scale by torus current twice !!
c         pcutl=pnorm.gt.(.6+4.*(thetad/45.)**6)*abs(torcur)/tmax
         pcutl=pnorm.gt.(.6+4.*(thetad/45.)**6)
         if((thetad.gt.thcut).and.pcutl)then
            d_phi=25.*sin((thetad-5.)*d2r)**0.15* 
     *           (pnorm*0.089)**(1./24.)
c         print *,'thetad,thcut,pcutl,phis,dphi',thetad,thcut,pcutl,
c     $        phis,d_phi
            if(abs(phis).lt.d_phi)then 
               weight=1.
            else
               weight=0.
            endif
         else
            weight=0.
            d_phi=0.
         endif
c
      elseif(outbend.and.central)then
         d_phi=30.
         weight=1.
      elseif(neutral.and.fwd)then
         if(p.lt.0.2)return
         sinth=sin(thetad*d2r)
         costh=cos(thetad*d2r)
         sinph=sin(phis*d2r)
         cosph=cos(phis*d2r)
c     
         r=(ec_r*(sin(ec_th)**2+cos(ec_th))+dl)/
     1        (costh+sinth*cosph*sin(ec_th))
         xh=r*cosph*sinth
         yh=r*sinph*sinth
         zh=r*costh-dl
c     
         call ec_xyz_duvw(xh,yh,zh,uh,vh,wh,xi,yi)
c     
         rcoil=dl*costh+sqrt(rc0**2-(dl*sinth)**2)
         hcoil=rcoil*sinth*cosph
         dcoil=hcoil*sin(0.5236)/cos(0.5236)-abs(rcoil*sinth*sinph)
c     
         if(yi.gt.ylow+10..and.yi.lt.yhi-5.)then
            xcut=(yi-ylow)/tgrho
            if(abs(xi).lt.xcut-5.)then
               weight=1
            elseif(abs(xi).gt.xcut.and.dcoil.lt.6.)then
               weight=1
            endif
            d_phi=(1.-5./(abs(r*sinph*sinth)+8.))*30.
         endif
c      
      elseif(neutral.and.central)then
        dth=atan(2./r/sqrt(p))*r2d
            d_phi=30.
            weight=1.
            p=p+0.05*ran(1)*sqrt(p)
            thetad=thetad+ran(2)*dth
            phis=phis+ran(3)*dth
      endif
      end
c
c
      subroutine smear_fun(ran,charge,torcur,p,th,phi,pout,thout,phiout)
      implicit none
      real ran(3),charge,torcur,p,th,phi,pout,thout,phiout
      logical inbend,outbend,neutral,fwd,central
      real pnorm,tmax
      real dthe,dphi
      data dthe,dphi,tmax/0.0286,0.0573,3375./
      integer lrun
      data lrun/3/
c common block containing track_type; defined in sub. tracktype
      common/track_type/inbend,outbend,neutral,fwd,central
c
c  since our acceptance functions care about sagitta and not momentum, per se,
c  we will form a "normalized momentum", pnorm
c
      pnorm=p*tmax/abs(torcur)
c
c  select which of the 6 functions we want
      if(inbend.and.fwd)then
         pout=p+ran(1)*0.005*(p*tmax/torcur)
         thout=th+ran(2)*dthe
         phiout=phi+ran(3)*dphi
c
      elseif(inbend.and.central)then
         pout=p+ran(1)*p*0.02
         thout=th+ran(2)*dthe*2.
         phiout=phi+ran(3)*dphi
c
      elseif(outbend.and.fwd)then
         pout=p+ran(1)*p*0.007
         thout=th+ran(2)*dthe
         phiout=phi+ran(3)*dphi
c
      elseif(outbend.and.central)then
         pout=p+ran(1)*p*0.02
         thout=th+ran(2)*dthe*2.
         phiout=phi+ran(3)*dphi
c
      elseif(neutral.and.fwd)then
         pout=p+0.05*ran(1)*sqrt(p)
         thout=th+ran(2)*dthe
         phiout=phi+ran(3)*dphi
c
      elseif(neutral.and.central)then
         pout=p+0.05*ran(1)*sqrt(p)
         thout=th+ran(2)*dthe
         phiout=phi+ran(3)*dphi
      endif
      end
c

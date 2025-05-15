       subroutine dvcsycol(d,xb,q2,ycol)
       implicit none
       double precision d,xb,q2,ycol
       ycol=(q2+d)/(q2+xb*d)
       return
       end
c
c

       real function thetagg()
       implicit none
#include "ntupgdvcs.inc"
#include "dvcs.inc"
       integer ip
       real ebeam,pien,vdotm,qmod,gmod,qg
       real pi4(4),qiu4(4),el04(4),elf4(4)
c
        ebeam=cl_beam_energy
        if(gngam.eq.0) print *,'Bad gamma'
        pien=ggame
        pi4(4)=pien
        pi4(1)=pien*cos(ggamf)*sin(ggamt) 
        pi4(2)=pien*sin(ggamf)*sin(ggamt)
        pi4(3)=pien*cos(ggamt)
c
         pien=sqrt(ebeam*ebeam-mele*mele)
        el04(1)=0
        el04(2)=0
        el04(3)=pien
        el04(4)=ebeam
c
         pien=sqrt(gelee*gelee-mele*mele)
        elf4(1)=pien*sin(gelet)*cos(gelef)
        elf4(2)=pien*sin(gelet)*sin(gelef)
        elf4(3)=pien*cos(gelet)
        elf4(4)=gelee
c
         call vdifm(el04,elf4,qiu4,4)
       qmod=vdotm(qiu4,qiu4,3)
       gmod=vdotm(pi4,pi4,3)
       qg=  vdotm(pi4,qiu4,3)
       if(qmod.gt.0.0 .and. gmod.gt.0.0) then
         thetagg=qg/sqrt(qmod)/sqrt(gmod)
       else
         thetagg=-1
       endif
c
       return
       end
         
c

       real function mismas(ip)
       implicit none
#include "ntupgdvcs.inc"
#include "dvcs.inc"
       integer ip
c
c      calculating azimuthal angle from lab variables
c      of hadron  
c
c      input variables: pie,pit,pif-momentum,theta,phi of hadron
c                       eleq2,elee,elef-Q^2,momentum,phi of eleitron                 
c
c     CALL CROSS(A,B,C) C=[AxB]
c     VDOT(A,B,N), CALL VMUL(A,B,X,N) X_i=A_i.B_i  N=3  reals
c     VDOTN(A,B,N)=ab/|a||b|
c     VMOD (A,N)  =|a|
       real tetgam,anu,pien,beam
       real eleq,pitg1,cospff1
       real xm,tetgam1,elety,ebeam
       real pi4(4),qiu4(4),el04(4),elf4(4),tnorm(4)
       real pro4(4),tnorm2(4)
       real vmass,vmass2,vangle,vdotm,phigstar
        ebeam=cl_beam_energy
c
c     define all 4momenta
c
c
       if(ip.eq.0) then
c
c      pi0
c
        if(gnpi0.eq.0) print *,'Bad pi0'
        pien=sqrt(gpi0e*gpi0e-0.135*0.135)
        pi4(4)=gpi0e
        pi4(1)=pien*cos(gpi0f)*sin(gpi0t) 
        pi4(2)=pien*sin(gpi0f)*sin(gpi0t)
        pi4(3)=pien*cos(gpi0t)
c
       elseif(ip.eq.2) then
c      real gamma
c
        if(gngam.eq.0) print *,'Bad gamma'
        pien=ggame
        pi4(4)=pien
        pi4(1)=pien*cos(ggamf)*sin(ggamt) 
        pi4(2)=pien*sin(ggamf)*sin(ggamt)
        pi4(3)=pien*cos(ggamt)
c
       elseif(ip.eq.4) then
c
c      pi0
c
        if(gnpro.eq.0) print *,'Bad pro'
        pien=sqrt(gproe*gproe-Mp*Mp)
        pi4(4)=gproe
        pi4(1)=pien*cos(gprof)*sin(gprot) 
        pi4(2)=pien*sin(gprof)*sin(gprot)
        pi4(3)=pien*cos(gprot)
       endif
c       Now el and gamma
c
c
c       el0
c
         pien=sqrt(ebeam*ebeam-mele*mele)
        el04(1)=0
        el04(2)=0
        el04(3)=pien
        el04(4)=ebeam
c
         pien=sqrt(gelee*gelee-mele*mele)
        elf4(1)=pien*sin(gelet)*cos(gelef)
        elf4(2)=pien*sin(gelet)*sin(gelef)
        elf4(3)=pien*cos(gelet)
        elf4(4)=gelee
c
c     gamma*
c
c proton
c
       pro4(4)=Mp
       pro4(3)=0
       pro4(2)=0
       pro4(1)=0
c
c
         call vdifm(el04,elf4,qiu4,4)
c
         call vsumm(qiu4,pro4,tnorm,4)
         call vdifm(tnorm,pi4,tnorm2,4)
         mismas=vmass2(tnorm2)
         return
         end
c
       real function vmass2(a)
       implicit none
       real vm,a(4), vdotm
        vm= vdotm(a,a,4)
        vmass2=-vm
       return
       end
c 
       real function phistar(ip)
       implicit none
#include "ntupgdvcs.inc"
#include "dvcs.inc"
       integer ip
       real tetgam,anu,pien,beam
       real eleq,pitg1,cospff1
       real xm,tetgam1,elety,ebeam
       real pi4(4),qiu4(4),el04(4),elf4(4),tnorm(4)
       real pro4(4),prf4(4),tnorm2(4)
       real vmass,vmass2,vangle,vdotm,phigstar
        ebeam=cl_beam_energy
c
c     define all 4momenta
c
        el04(1)=0
        el04(2)=0
        el04(3)=pien
        el04(4)=ebeam
c
         pien=sqrt(gelee*gelee-mele*mele)
        elf4(1)=pien*sin(gelet)*cos(gelef)
        elf4(2)=pien*sin(gelet)*sin(gelef)
        elf4(3)=pien*cos(gelet)
        elf4(4)=gelee
c
c proton
c
        pro4(4)=Mp
        pro4(3)=0
        pro4(2)=0
        pro4(1)=0
c
         call vdifm(el04,elf4,qiu4,4)
c
       if(ip.eq.0) then
c
c      pi0
c
        if(gnpi0.eq.0) print *,'Bad pi0'
        pien=sqrt(gpi0e*gpi0e-0.135*0.135)
        pi4(4)=gpi0e
        pi4(1)=pien*cos(gpi0f)*sin(gpi0t) 
        pi4(2)=pien*sin(gpi0f)*sin(gpi0t)
        pi4(3)=pien*cos(gpi0t)
c
       elseif(ip.eq.2) then
c      real gamma
c
        if(gngam.eq.0) print *,'Bad gamma'
        pien=ggame
        pi4(4)=pien
        pi4(1)=pien*cos(ggamf)*sin(ggamt) 
        pi4(2)=pien*sin(ggamf)*sin(ggamt)
        pi4(3)=pien*cos(ggamt)
c
       elseif(ip.eq.4) then
c
c      pi0
c
        if(gnpro.eq.0) print *,'Bad pro'
        pien=sqrt(gproe*gproe-Mp*Mp)
        prf4(4)=gproe
        prf4(1)=pien*cos(gprof)*sin(gprot)
        prf4(2)=pien*sin(gprof)*sin(gprot)
        prf4(3)=pien*cos(gprot)
c
c       el0
c
         pien=sqrt(ebeam*ebeam-mele*mele)
         call vsumm(qiu4,pro4,tnorm,4)
         call vdifm(tnorm,prf4,pi4,4)

       endif
c       Now el and gamma
c
c     gamma*
c
c
c
         call crossm(qiu4,el04,tnorm)
         if (vdotm(tnorm,pi4,3).ge.0) then    ! for pi+
          phistar=vangle(qiu4,el04,qiu4,pi4)
         else
          phistar=2*pi-vangle(qiu4,el04,qiu4,pi4)
         endif
         return
         end

c
c
       real function tchan(ip)
       implicit none
#include "ntupgdvcs.inc"
#include "dvcs.inc"
       integer ip,iel
c
c      calculating azimuthal angle from lab variables
c      of hadron  
c
c      input variables: pie,pit,pif-momentum,theta,phi of hadron
c                       eleq2,elee,elef-Q^2,momentum,phi of eleitron                 
c
c     CALL CROSS(A,B,C) C=[AxB]
c     VDOT(A,B,N), CALL VMUL(A,B,X,N) X_i=A_i.B_i  N=3  reals
c     VDOTN(A,B,N)=ab/|a||b|
c     VMOD (A,N)  =|a|
       real tetgam,anu,pien,beam
       real eleq,pitg1,cospff1
       real xm,tetgam1,elety,ebeam
       real pi4(4),qiu4(4),el04(4),elf4(4),tnorm(4)
       real pro4(4),tnorm2(4)
       real vmass,vangle,vdotm,phigstar
        ebeam=cl_beam_energy
c
c     define all 4momenta
c
       if(ip.eq.0) then
c
c      pi0
c
        if(gnpi0.eq.0) print *,'Bad gpi0'
        pien=sqrt(gpi0e*gpi0e-0.135*0.135)
        pi4(4)=gpi0e
        pi4(1)=pien*cos(gpi0f)*sin(gpi0t) 
        pi4(2)=pien*sin(gpi0f)*sin(gpi0t)
        pi4(3)=pien*cos(gpi0t)
c
       elseif(ip.eq.2) then
c      real gamma
c
        if(gngam.eq.0) print *,'Bad ggamma'
        pien=ggame
        pi4(4)=pien
        pi4(1)=pien*cos(ggamf)*sin(ggamt) 
        pi4(2)=pien*sin(ggamf)*sin(ggamt)
        pi4(3)=pien*cos(ggamt)
c
       elseif(ip.eq.4) then
c
c      gpi0
c
        if(gnpro.eq.0) print *,'Bad pro'
        pien=sqrt(gproe*gproe-Mp*Mp)
        pi4(4)=gproe
        pi4(1)=pien*cos(gprof)*sin(gprot) 
        pi4(2)=pien*sin(gprof)*sin(gprot)
        pi4(3)=pien*cos(gprot)
       endif
c       Now el and ggamma
c
c
c       el0
c
         pien=sqrt(ebeam*ebeam-mele*mele)
        el04(1)=0
        el04(2)=0
        el04(3)=pien
        el04(4)=ebeam
c
         pien=sqrt(gelee*gelee-mele*mele)
        elf4(1)=pien*sin(gelet)*cos(gelef)
        elf4(2)=pien*sin(gelet)*sin(gelef)
        elf4(3)=pien*cos(gelet)
        elf4(4)=gelee
c
c     ggamma*
c
c proton
c
       pro4(4)=Mp
       pro4(3)=0
       pro4(2)=0
       pro4(1)=0
c
         call vdifm(el04,elf4,qiu4,4)
c
         if(ip.lt.3) then
         call vdifm(qiu4,pi4,tnorm,4)
         tchan=vdotm(tnorm,tnorm,4)
         else
         call vdifm(pro4,pi4,tnorm,4)
         tchan=vdotm(tnorm,tnorm,4)
         endif         
         return
         end
       


       real function gamel(ip,iel)
       implicit none
#include "ntupgdvcs.inc"
#include "dvcs.inc"
       integer ip,iel
c
c      calculating azimuthal angle from lab variables
c      of hadron  
c
c      input variables: pie,pit,pif-momentum,theta,phi of hadron
c                       eleq2,elee,elef-Q^2,momentum,phi of eleitron                 
c
c     CALL CROSS(A,B,C) C=[AxB]
c     VDOT(A,B,N), CALL VMUL(A,B,X,N) X_i=A_i.B_i  N=3  reals
c     VDOTN(A,B,N)=ab/|a||b|
c     VMOD (A,N)  =|a|
       real tetgam,anu,pien
       real eleq,pitg1,cospff1
       real xm,tetgam1,elety,ebeam
       real pi4(4),qiu4(4),el04(4),elf4(4),tnorm(4)
       real pro4(4),prf4(4),tgam(4)
       real vmass,vmass2,vangle,vdotm,phigstar,ang
        ebeam=cl_beam_energy
c
c     define all 4momenta
c
c
c       el0
c
         pien=sqrt(ebeam*ebeam-mele*mele)
        el04(1)=0
        el04(2)=0
        el04(3)=pien
        el04(4)=ebeam
c
         pien=sqrt(gelee*gelee-mele*mele)
        elf4(1)=pien*sin(gelet)*cos(gelef)
        elf4(2)=pien*sin(gelet)*sin(gelef)
        elf4(3)=pien*cos(gelet)
        elf4(4)=gelee
c
c     gamma*
c
c proton -0
c
       pro4(4)=Mp
       pro4(3)=0
       pro4(2)=0
       pro4(1)=0

         call vdifm(el04,elf4,qiu4,4)

       if(ip.eq.0) then
         if (gnpro.lt.1) then
           print *,'bad proton'
           return
         endif
c
c       define final proton -f
c
        pien=sqrt(gproe*gproe-Mp*Mp)
        prf4(4)=gproe
        prf4(1)=pien*cos(gprof)*sin(gprot)
        prf4(2)=pien*sin(gprof)*sin(gprot)
        prf4(3)=pien*cos(gprot)
         call vsumm(qiu4,pro4,tnorm,4)
         call vdifm(tnorm,prf4,tgam,4)

       elseif(ip.eq.1) then
         if (gngam.lt.1) then
           print *,'bad ggam1'
           return
         endif
        pien=ggame
        tgam(4)=pien
        tgam(1)=pien*cos(ggamf)*sin(ggamt) 
        tgam(2)=pien*sin(ggamf)*sin(ggamt)
        tgam(3)=pien*cos(ggamt)
c
       elseif(ip.eq.2) then
         if (gngam.lt.2) then
           print *,'bad gam2'
           return
         endif
        pien=ggam2e
        tgam(4)=pien
        tgam(1)=pien*cos(ggam2f)*sin(ggam2t) 
        tgam(2)=pien*sin(ggam2f)*sin(ggam2t)
        tgam(3)=pien*cos(ggam2t)
       elseif(ip.eq.3) then
         if (gnpi0.eq.0) then
           print *,'bad gpi0'
           return
         endif
        pien=sqrt(gpi0e*gpi0e-0.135*0.135)
        tgam(4)=gpi0e
        tgam(1)=pien*cos(gpi0f)*sin(gpi0t) 
        tgam(2)=pien*sin(gpi0f)*sin(gpi0t)
        tgam(3)=pien*cos(gpi0t)

       elseif(ip.eq.4) then
         call vdifm(el04,elf4,tgam,4)
       endif

c
c
         if(iel.eq.1) then
          gamel=ang(tgam,el04)
         elseif (iel.eq.2) then
          gamel=ang(tgam,elf4)
         elseif (iel.eq.3) then 
          gamel=ang(tgam,qiu4)
         else
          print *,'dmboooooooooooo'
         endif
c

         return
         end

c
       real function ang(x,y)
       real x(4),y(4),xm,ym,vcos,vdotm
         xm=vdotm(x,x,3)
         ym=vdotm(y,y,3)
       if(xm.gt.0.0 .and. ym.gt.0.0) then
         vcos=vdotm(x,y,3)/sqrt(xm)/sqrt(ym)
       else
         vcos=1.0
       endif
         if(abs(vcos).lt.1.0) then
            ang=acos(vcos)
         else
c            print *,xm,ym,vcos
            ang=-1
         endif
       return
       end


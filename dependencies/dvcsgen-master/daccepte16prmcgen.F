       real function daccepte16prmcgen()
       implicit none
#include "dvcspi0gen.inc"
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
       real amp,ebeam,phirad,deltz,mismasrho,ptran
       real pi,degrad,kf_mag,phi_edn,the_ed,accvb2,accvb1,proex
       integer isec
       pi=acos(-1.0)
       degrad=pi/180.0
        amp=0.93827
        ebeam=5.754
c     define all 4momenta
        daccepte16prmcgen=0
c
       phirad=gprof/degrad
        call newphi(phirad,phi_edn,isec)
          kf_mag=proex(gproe)
          kf_mag=sqrt(kf_mag*kf_mag-amp*amp)
         the_ed=gprot*180/3.1415  
         daccepte16prmcgen=accvb1(1,kf_mag,the_ed,phi_edn,3375.0)
         return
         end
c
c
        real function proex(e)
        real p,pc,e
        p=sqrt(e*e-0.8803)
        pc=p-0.15518+0.80867*p-1.6268*p*p+1.4419*p*p*p-0.47052*p*p*p*p
        if(pc.lt.1.0) then
         proex=sqrt(pc*pc+0.8803)
        else
         proex=e
        endif
        return
        end


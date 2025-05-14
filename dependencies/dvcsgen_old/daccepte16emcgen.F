       integer function daccepte16emcgen()
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
       real amp,ebeam,phirad,deltz,mismasrho,ptran
       real degrad,kf_mag,phi_edn,the_ed,accvb2,accvb1
       integer isec
       degrad=pi/180.0
       ebeam=cl_beam_energy
c     define all 4momenta
        daccepte16emcgen=0
c
       if((1.0-gelee/ebeam).gt.0.85 ) return
c
       phirad=gelef/degrad
        call newphi(phirad,phi_edn,isec)
          kf_mag=gelee
         the_ed=gelet*180/3.1415  
         daccepte16emcgen=accvb1(0,kf_mag,the_ed,phi_edn,3375.0)
         return
         end

C======================================================================
      FUNCTION accvb1(part_type,p,theta,phi,t_current)
C----------------------------------------------------------------------
C-
C- Purpose and Methods : 
C-
C----------------------------------------------------------------------
      IMPLICIT NONE
C----------------------------------------------------------------------
C
      REAL p,theta,phi,t_current,acc  
      REAL t_max
      REAL phi0_el, phi0_nh, phi0_ph
      REAL theta0_el, theta0_nh, theta0_ph
      REAL thetas_el, thetas_nh, thetas_ph
      REAL p_shift, cel_ex, pel_ex
      REAL ch_ex,theta_cut
      REAL theta_min, theta_max,delta_phi, exp
      INTEGER part_type,electron,pos_hadron, neg_hadron
      REAL d2r
      real accvb1
      
      data t_max/3375./
      data phi0_el/30./
      data phi0_nh/25./
      data phi0_ph/25./
      data theta0_el/12.5/
      data theta0_nh/10./
      data theta0_ph/5./
      data thetas_el/15./
      data thetas_nh/15./
      data thetas_ph/25./
      data theta_max/50./
      data p_shift/0.15/
      data pel_ex/0.333/
      data cel_ex/0.35/
      data theta_cut/75./
      data ch_ex/0.3/
      data electron/0/
      data pos_hadron/1/
      data neg_hadron/-1/
c      data pi/3.1415926/
ccc      data d2r/0.01754533/      
      data d2r/0.02754533/      
     
      Acc = 0.0
      IF(part_type.EQ.electron) THEN
        theta_min = theta0_el+thetas_el/(p*t_max/t_current+p_shift)
        if(theta.gt.theta_min.and.theta.lt.theta_max)then
          exp = cel_ex*(p*t_max/t_current)**pel_ex
          delta_phi = phi0_el*sin((theta-theta_min)*d2r)**exp
          if(abs(phi).lt.delta_phi)Acc=1.
        endif
      ELSEIF(part_type.EQ.pos_hadron) THEN
        theta_min = theta0_ph+thetas_ph*(1.-p*t_max/t_current/5.)**24
        if(theta.gt.theta_min)then
          exp=(p*t_max/t_current/5.)**(1./8.)
          delta_phi = phi0_ph*cos((theta-theta_cut)*d2r)**exp
          if(abs(phi).lt.delta_phi)Acc=1.
        endif
      ELSEIF(part_type.EQ.neg_hadron) THEN
        theta_min = theta0_nh+thetas_nh/(p*t_max/t_current+p_shift)
        if(theta.gt.theta_min)then
          exp = ch_ex*(p*t_max/t_current)**pel_ex
          delta_phi = phi0_nh*sin((theta-theta_min)*d2r)**exp
          if(abs(phi).lt.delta_phi)Acc=1.
        endif
      ENDIF
      
      accvb1 = acc
      
      end

C======================================================================
      FUNCTION accvb2(part_type,p,theta,phi,t_current)
C----------------------------------------------------------------------
C-
C- Purpose and Methods : 
C-
C----------------------------------------------------------------------
      IMPLICIT NONE
C----------------------------------------------------------------------
C
      REAL p,theta,phi,t_current,acc  
      REAL t_max
      REAL phi0_el, phi0_nh, phi0_ph
      REAL theta0_el, theta0_nh, theta0_ph
      REAL thetas_el, thetas_nh, thetas_ph
      REAL p_shift, cel_ex, pel_ex
      REAL ch_ex,theta_cut
      REAL theta_min, theta_max,delta_phi, exp
      INTEGER part_type,electron,pos_hadron, neg_hadron
      REAL d2r
      real accvb2
      
      data t_max/3375./
      data phi0_el/32./
      data phi0_nh/27./
      data phi0_ph/27./
      data theta0_el/10.5/
      data theta0_nh/8./
      data theta0_ph/8./
      data thetas_el/15./
      data thetas_nh/15./
      data thetas_ph/25./
      data theta_max/50./
      data p_shift/0.15/
      data pel_ex/0.333/
      data cel_ex/0.35/
      data theta_cut/75./
      data ch_ex/0.3/
      data electron/0/
      data pos_hadron/1/
      data neg_hadron/-1/
c      data pi/3.1415926/
      data d2r/0.01754533/      
      Acc = 0.0
      
      IF(part_type.EQ.electron) THEN
        theta_min = theta0_el+thetas_el/(p*t_max/t_current+p_shift)
        if(theta.gt.theta_min.and.theta.lt.theta_max)then
          exp = cel_ex*(p*t_max/t_current)**pel_ex
          delta_phi = phi0_el*sin((theta-theta_min)*d2r)**exp
          if(abs(phi).lt.delta_phi)Acc=1.
        endif
      ELSEIF(part_type.EQ.pos_hadron) THEN
        theta_min = theta0_ph+thetas_ph*(1.-p*t_max/t_current/5.)**24
        if(theta.gt.theta_min)then
          exp=(p*t_max/t_current/5.)**(1./8.)
          delta_phi = phi0_ph*cos((theta-theta_cut)*d2r)**exp
          if(abs(phi).lt.delta_phi)Acc=1.
        endif
      ELSEIF(part_type.EQ.neg_hadron) THEN
        theta_min = theta0_nh+thetas_nh/(p*t_max/t_current+p_shift)
        if(theta.gt.theta_min)then
          exp = ch_ex*(p*t_max/t_current)**pel_ex
          delta_phi = phi0_nh*sin((theta-theta_min)*d2r)**exp
          if(abs(phi).lt.delta_phi)Acc=1.
        endif
      ENDIF
      
      accvb2 = acc
      
      end   


C======================================================================
      subroutine newphi(phi,phinew,sec)
C----------------------------------------------------------------------
C-
C- Purpose and Methods : 
C-
C----------------------------------------------------------------------
      IMPLICIT NONE
C----------------------------------------------------------------------
C
C
      real phinew, phi
      integer sec
        phinew = phi
        sec=1
c
      if (phi.gt.330.) then
        phinew = phi-360.
        sec=1
      elseif (phi.gt.30.0.and.phi.le.90.) then
        phinew = phi-60.
        sec=2
      elseif (phi.gt.90.0.and.phi.le.150.) then
        phinew = phi-120.
        sec=3
      elseif (phi.gt.150.0.and.phi.le.210.) then
        phinew = phi-180.
        sec=4
      elseif (phi.gt.210.0.and.phi.le.270.) then
        phinew = phi-240.
        sec=5
      elseif (phi.gt.270.0.and.phi.le.330.) then
        phinew = phi-300.
        sec=6
      endif
      
      end  

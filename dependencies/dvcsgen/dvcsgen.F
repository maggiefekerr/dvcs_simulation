         PROGRAM DVCSGEN
cc
        IMPLICIT NONE
c
C******************************************************************
C
C       Main program for DVCS simulations for CLAS kinematics
C
C***************************************************************** 
#include "dvcs.inc"

C..  Number of events 
      INTEGER       I,NEVENT
c
      call CL_OPTIONS                     ! command line processing
#ifdef  hbook
      if(ntOK) call INIT_PAW()            ! init hbook if requested
#endif
c
C..   number of events   
c
       NEVENT=cl_triggers
c
       CALL DVCSINIT       ! init limits/cuts for dvcs/pi0
c
       if(cl_verblev.gt.0) WRITE(*,*) 'NEVENT=',NEVENT
c
C.. Event generation: 
c
        DO I = 1,NEVENT
c
         CALL CLEANNT
         CALL genEVENT(I)       
         CALL CLASFILL(I)    ! process event
c
         IF (MOD(I,cl_nprint).EQ.0.and.cl_verblev.gt.0)
     6   WRITE(*,*) I,'  from ',NEVENT,cl_nprint,cl_be
c
        ENDDO
c
         if (datfileOK) close(41)
#ifdef  hbook
         if(ntOK) call CLOSE_PAW()
#endif
c
       WRITE(*,*) 'Total ',NEVENT,' in pb ',smax
C.. Saving the last random number configuration 
       END      
c
c
c
#ifdef  hbook
       SUBROUTINE INIT_PAW()
       implicit none
       integer h,nwpawc
       parameter (nwpawc=1000000)
       integer  NTID_OUT

       common /PAWC/ h(nwpawc)
       call HLIMIT(nwpawc)
       NTID_OUT=22
       call hini(NTID_OUT)
       return
       end
c
         subroutine  CLOSE_PAW()
         implicit none
         integer icycle
          call HCDIR('//ZXC',' ')
          call HROUT(0,icycle,'')
          call HREND('ZXC')
          close(50)
         return
         end
c 
      subroutine hini(NTID_OUT)
      implicit none
#include "ntupgdvcs.inc"
#include "dvcs.inc"
#include "dvcsmom.inc"
c
      integer NTID_OUT,NTID_OUTM,IQUEST
      COMMON/QUEST/IQUEST(100) 
      integer lrec,istat
      character*14 hbooknm
      if (cl_pi0) then
       write(hbooknm,'(A,A,A)') 'p',bosout,'.hbook'
      else if (cl_eta) then
       write(hbooknm,'(A,A,A)') 'e',bosout,'.hbook'
      else
       write(hbooknm,'(A,A,A)') 'd',bosout,'.hbook'
      endif
      lrec=8191
c     open koko
c
       IQUEST(10)=65000
       print *,'INIT PAW hbook=',hbooknm
       CALL HROPEN(50,'ZXC',hbooknm,'NQE',LREC,ISTAT)
       call hbset('BSIZE',8176,istat)
      if (istat.ne.0) stop
      call HBNT(NTID_OUT,'semiout',' ')
      call HBNAME(NTID_OUT,'elegen',gnele,
     +'gnele:I,gelebits:I,geleh:I,'//
     +'gelee:R,gelet:R,'//
     +'gelef:R,gelex:R,geley:R,'//
     +'gelez:R,gelewgt:R,gxb:R,gq2:R,gw:R')
c234567890123456789012345678901234567890123456789012345678901234567890123456
      call HBNAME(NTID_OUT,'hadgen',gnpi0,
     +'gnpi0:I,gnpro:I,gngam:I,'//
     +'gnneu:I,gpi0e:R,gpi0t:R,gpi0f:R,'//
     +'gproe:R,gprot:R,gprof:R,ggame:R,'//
     +'ggamt:R,ggamf:R,'//
     +'ggam2t:R,ggam2e:R,ggam2f:R,'//
     +'gprox:R,gproy:R,gproz:R,'//
     +'gprotg:R,gprotl:R,gprotv:R,gprotmis:R,gprott:R,gprofs:R,'//
     +'gwbh:R,gwdvcs:R,gvint:R,gsin:R,gsin2:R')
      if (cl_mom) then    
        if(cl_pol.eq.1) then
       call HBNAME(NTID_OUT,'mom0',gproh,
     +'gproh:I,ghp1:R,ghp2:R,ghccb:R,ghcci:R,gh0bh:R,gh1bh:R,gh2bh:R,'//
     +'gh0dvcs:R,gh1dvcs:R,ghs1dvcs:R,'//
     +'gh0iunp:R,gh1iunp:R,ghs1iunp:R,ghs2iunp:R,'//
     +'gh0bhlp:R,gh1bhlp:R,ghs1ilp:R,gh1ilp:R,gh0ilp:R')
c
       else if(cl_pol.eq.2) then
       call HBNAME(NTID_OUT,'mom0',gproh,
     +'gproh:I,gh0bh:R,gh1bh:R,gh2bh:R,'//
     +'gh0dvcs:R,gh1dvcs:R,ghs1dvcs:R,'//
     +'gh0iunp:R,gh1iunp:R,ghs1iunp:R,ghs2iunp:R,'//
     +'gh0bhtpcos:R,gh1bhtpcos:R,ghs1bhtpsin:R,'//
     +'ghs1itpcos:R,ghs1itpsin:R,gh1itpcos:R,'//
     +'gh1itpsin:R,gh0itpcos:R,gh0itpsin:R') 
c
       else
       call HBNAME(NTID_OUT,'mom0',gproh,
     +'gproh:I,gh0bh:R,gh1bh:R,gh2bh:R,'//
     +'gh0dvcs:R,gh1dvcs:R,ghs1dvcs:R,'//
     +'gh0iunp:R,gh1iunp:R,ghs1iunp:R,ghs2iunp:R')
       endif
      endif
c
      if (cl_radgen) then    
       call HBNAME(NTID_OUT,'radgen',iproctype,
     +'iproctype:I,radxb:R,radxb_true:R,radq2:R,radq2_true:R,radphi:R,radt:R,radnu:R,bornweight:R,radweight:R,'//
     +'grade:R,gradt:R,gradf:R')
      endif
      end
c
#endif
c
       SUBROUTINE CL_OPTIONS
       implicit none
#include "dvcs.inc"
       integer i,iran,irannom,ivalnum
       real VEC1,valnum,rndnmbr
       character*20 cnumber
       integer  numopts,IARGC
       integer*4 now(3)
c
c    command line processing
c
c     defaults
c
      cl_dvcs=.TRUE.         ! default dvcs on
      cl_pi0 =.FALSE.        ! default no pi0
      cl_eta =.FALSE.        ! default no etta
      cl_bh=3                ! full x-section 1-only BH
      cl_gpd=3               ! gpd model
      cl_scale=1.0           ! gpd model
      cl_bpi0=1.3            ! b=1.7
      cl_ycol=0.005          ! default min for P1
      cl_xpos=0              ! x-position
      cl_ypos=0              ! y-position
      cl_zpos=0              ! z-position
      cl_rast=0.0           ! raster diameter in cm
      cl_zwidth=0            ! z-width ( zpos+/-zwidth)
      cl_pol=0               ! unpolarized target
      cl_beam_energy=10.6    ! def clas12
      cl_wmin=4.0            ! def w2min
      cl_ymin=0.05           ! def ymin -> Emax=(1-ymin)*E
      cl_ymax=0.9            ! def ymax
      cl_thmin=0.09           ! def e'tmin
      cl_thmax=1.57            ! def e'tmax
      cl_tmax=1.0             ! def tmax
      cl_tmin=0.0             ! def tmin
      cl_phigd=0.0            ! def phigd for rc factor
      cl_q2max=15.0            ! def Q2 max
      cl_q2min=1.0            ! def Q2 min
      cl_xmax=0.75             ! def Q2 max
      cl_xmin=0.05             ! def Q2 min
      cl_delta=1d-2             ! def delta (minimum rad photon energy)
      cl_vv2cut=0.1     ! def vv2cut (minimum missing mass_ep squared)
      cl_target='proton'     ! (def proton target)
c
      cl_proloss=0            ! no proton loss by default
      cl_smear=0            ! no smearing by default
      cl_sma=0.006          ! A
      cl_smb=0.001         ! B
      cl_smc=0.0008         ! C
      cl_smd=0.001         ! D
      cl_sme=0.11          ! photon E
      cl_smf=0.003         ! F
      cl_smg=0.015         ! G
c
      cl_seed=0
      cl_verblev=0
      cl_writef=0            ! 0-clas12 1-gsim lund format
      cl_nprint=1000         ! print every cl_nprint event
      cl_triggers = 10000  
      cl_nmax = 2000         ! max number of events in the file
      cl_printgpd=.FALSE.
      cl_printrad=.FALSE.
      cl_printdstot=.FALSE.
      cl_mom=.FALSE.
      cl_ktcor=.FALSE.
      cl_radgen=.FALSE.
      cl_radstable=.FALSE.
      cl_mod=0              ! write all, otherwise filter
      datfileOK=.TRUE.      ! write a data file
      weightOK=.FALSE.
      acce16=.FALSE.
      acceg1=.FALSE.
      acc12=.FALSE.
      ntOK=.TRUE.
      cl_docker=.FALSE.
      bosout='dvcsgen'
      numopts  = IARGC()
c      
      if(numopts.gt.0 ) then
        i=0
      else
        goto 2
      endif
 1    i=i+1
        CALL GETARG(i,cnumber)
        if(cnumber.eq.'--trig'.and.i.lt.numopts) then
          i=i+1
          CALL GETARG(i,cnumber)
          cl_triggers=valnum(cnumber)
        elseif(cnumber.eq.'--nmax'.and.i.lt.numopts) then
          i=i+1
          CALL GETARG(i,cnumber)
          cl_nmax=valnum(cnumber)
        elseif(cnumber.eq.'--v'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_verblev=valnum(cnumber)
        elseif(cnumber.eq.'--writef'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_writef=valnum(cnumber)
        elseif(cnumber.eq.'--mod'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_mod=valnum(cnumber)
        elseif(cnumber.eq.'--print'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_nprint=valnum(cnumber)
        elseif(cnumber.eq.'--lpol') then
           cl_pol=1
        elseif(cnumber.eq.'--tpol') then
           cl_pol=2
        elseif(cnumber.eq.'--radgen') then
           cl_radgen=.TRUE.
        elseif(cnumber.eq.'--radstable') then
           cl_radstable=.TRUE.
        elseif(cnumber.eq.'--mom') then
           cl_mom=.TRUE.
        elseif(cnumber.eq.'--pi0') then
           cl_pi0=.TRUE.
        elseif(cnumber.eq.'--eta') then
           cl_eta=.TRUE.
        elseif(cnumber.eq.'--nodvcs') then
           cl_dvcs=.FALSE.
        elseif(cnumber.eq.'--ktcor') then
           cl_ktcor=.TRUE.
        elseif(cnumber.eq.'--beam'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_beam_energy=valnum(cnumber)           
        elseif(cnumber.eq.'--bh'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_bh=valnum(cnumber)           
        elseif(cnumber.eq.'--scale'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_scale=valnum(cnumber)           
        elseif(cnumber.eq.'--t'.and.i.lt.numopts+1) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_tmin=valnum(cnumber)
           i=i+1
           CALL GETARG(i,cnumber)
           cl_tmax=valnum(cnumber)
        elseif(cnumber.eq.'--th'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_thmin=valnum(cnumber)
           i=i+1
           CALL GETARG(i,cnumber)
           cl_thmax=valnum(cnumber)
        elseif(cnumber.eq.'--q2'.and.i.lt.numopts+1) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_q2min=valnum(cnumber)
           i=i+1
           CALL GETARG(i,cnumber)
           cl_q2max=valnum(cnumber)
        elseif(cnumber.eq.'--x'.and.i.lt.numopts+1) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_xmin=valnum(cnumber)
           i=i+1
           CALL GETARG(i,cnumber)
           cl_xmax=valnum(cnumber)
        elseif(cnumber.eq.'--w'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_wmin=valnum(cnumber)
        elseif(cnumber.eq.'--phi'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_phigd=valnum(cnumber)
        elseif(cnumber.eq.'--y'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_ymin=valnum(cnumber)
           i=i+1
           CALL GETARG(i,cnumber)
           cl_ymax=valnum(cnumber)
        elseif(cnumber.eq.'--xpos'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_xpos=valnum(cnumber)
        elseif(cnumber.eq.'--ypos'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_ypos=valnum(cnumber)
        elseif(cnumber.eq.'--zpos'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_zpos=valnum(cnumber)
        elseif(cnumber.eq.'--ycol'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_ycol=valnum(cnumber)           
        elseif(cnumber.eq.'--gpd'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_gpd=valnum(cnumber)
        elseif(cnumber.eq.'--zwidth'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_zwidth=valnum(cnumber)
         elseif(cnumber.eq.'--raster'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_rast=valnum(cnumber)           
        elseif(cnumber.eq.'--file'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,bosout)
        elseif(cnumber.eq.'--delta'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_delta=valnum(cnumber)           
        elseif(cnumber.eq.'--vv2cut'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_vv2cut=valnum(cnumber)           
        elseif(cnumber.eq.'--heli'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_heli=valnum(cnumber)           
        elseif(cnumber.eq.'--targ'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cl_target)
        elseif(cnumber.eq.'--nodat') then
           datfileOK=.FALSE.
        elseif(cnumber.eq.'--printgpd') then
           cl_printgpd=.TRUE.
        elseif(cnumber.eq.'--printrad') then
           cl_printrad=.TRUE.
        elseif(cnumber.eq.'--printdstot') then
           cl_printdstot=.TRUE.
        elseif(cnumber.eq.'--acce16') then
           acce16=.TRUE.
        elseif(cnumber.eq.'--acc12') then
           acc12=.TRUE.
        elseif(cnumber.eq.'--acceg1') then
           acceg1=.TRUE.
        elseif(cnumber.eq.'--weight') then
           weightOK=.TRUE.
        elseif(cnumber.eq.'--docker') then
           cl_docker=.TRUE.
        elseif(cnumber.eq.'--nont') then
           ntOK=.FALSE.
        elseif(cnumber.eq.'--seed'.and.i.lt.numopts) then
          i=i+1
          CALL GETARG(i,cnumber)
          cl_seed=ivalnum(cnumber)
        elseif(cnumber.eq.'--proloss') then
          cl_proloss=1
        elseif(cnumber.eq.'--smear') then
          cl_smear=1
        elseif(cnumber.eq.'--A'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_sma=valnum(cnumber)           
        elseif(cnumber.eq.'--B'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_smb=valnum(cnumber)           
        elseif(cnumber.eq.'--C'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_smc=valnum(cnumber)           
        elseif(cnumber.eq.'--D'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_smd=valnum(cnumber)           
        elseif(cnumber.eq.'--E'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_sme=valnum(cnumber)           
        elseif(cnumber.eq.'--F'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_smf=valnum(cnumber)           
        elseif(cnumber.eq.'--G'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_smf=valnum(cnumber)           
        else
          call printusage
          stop
        endif
        if(i.lt. numopts) goto 1
 2      continue
        call printvalues
c
c randomize
c
        call RNDMSEED()
        return 
        end
c
        SUBROUTINE  printrad
        implicit none
#include "dvcs.inc"
C..
        integer nx,nq,nt,nfe,nfg
        real  xmax,xmin,q2max,q2min,tmax,tmin 
        common /xmax/  xmax,xmin,q2max,q2min,tmax,tmin 
     6                 ,nx,nq,nt,nfe,nfg
C..
C..
      double precision xbd,Q2d,del2d,phield,phigd,dstot,stotint,yb,w2
     6                ,avgdstot,avgprobn,avgprobs,avgprobp,xsec_born,xsec_obs,sum_born,sum_obs,unc_obs,var_obs
     6                ,del_born, del_obs, delu,delp, A_obs, delta_u, delta_p, delta_A
     6                , unc_A_obs, unc_delta_u, unc_delta_p, unc_delta_A, dsigma_born, dsigma_obs, unc_dsigma_born, unc_dsigma_obs
      integer la0,itrial,iloop,numn,nump,nums
c internal variables
c
        logical outoflimits
        real altime,alltime(2)
        real*8 t1lim, t2lim
c
        DO iloop=0,9
        la0=0
        sum_obs = 0
        sum_born = 0
        del_obs = 0
        del_born = 0
        if(cl_pol.ge.1) la0=-1
          DO helpi=la0,1,2
            DO heli=-1,1,2
               nums = 0
               nump = 0
               numn = 0
               avgdstot = 0
               avgprobs = 0
               avgprobp = 0
               avgprobn = 0
          xbd=xmin!+ix*dx
          Q2d=q2min!+iq*dq
           yb=Q2d/(2D0*Mp*xbd*cl_be)
           w2=Mp*Mp+2.0*Mp*cl_be*yb-Q2d
c
         if( .not.outoflimits(xbd,Q2d) ) then
c
              t2lim=-0.5d0*((q2d/xbd-q2d)*(q2d/xbd+sqrt((q2d/xbd)**2+4d0*Mp*Mp*q2d))+2d0*Mp*Mp*q2d)/w2
              t1lim=-0.5d0*((q2d/xbd-q2d)*(q2d/xbd-sqrt((q2d/xbd)**2+4d0*Mp*Mp*q2d))+2d0*Mp*Mp*q2d)/w2
              tmax = min(-cl_tmin, t1lim)
              tmin = max(-cl_tmin, t2lim)
              if (tmax<tmin) then
                 print *, "check -t limit"
                 return
              endif
              del2d=tmin!+it*dt
              phigd=cl_phigd
c
              CALL bmkxsec(xbd, Q2d, del2d, phield,phigd,dstot)
              xsec_born = dstot*2*pi*0.001

              print *, "Loop over 100 iterations for rad. cross sections"

              xff = xbd
              q2ff = Q2d
              tff = del2d
              phiff = pi-phigd
              ebeamff = cl_beam_energy
              ikeygene = 10
              DO itrial=0,99
                 CALL bhradgen(ebeamff, xff, q2ff, tff, phiff,cl_vv2cut,cl_delta,egamma,thetag,phig,ichannel,cl_pol,ikeygene,dstot,probn,probp,probs)
                 dstot = dstot*0.001
                 avgdstot = avgdstot + dstot/100
                 avgprobp = avgprobp + probp/100
                 avgprobs = avgprobs + probs/100
                 avgprobn = avgprobn + probn/100
                 if (ichannel.eq.1) then
                    numn = numn+1
C                     print *, egamma, dstot, ichannel, heli, helpi, xsec_born
                 endif
                 if (ichannel.eq.2) then
                    nums = nums+1
C                     print *, egamma, dstot, ichannel, heli, helpi
                 endif
                 if (ichannel.eq.3) then
                  nump = nump+1
C                   print *, egamma, dstot, ichannel, heli, helpi
                 endif
              ENDDO
              xsec_obs = avgdstot

              sum_obs = sum_obs + xsec_obs
              sum_born = sum_born + xsec_born
              del_obs = del_obs + xsec_obs*heli
              del_born = del_born + xsec_born*heli

              print *, "=======sanity check for reasonable event generation=========loop====",iloop
              print *, "average prob_n, average prob_s, average prob_p, nevents_n, nevents_s, nevents_p"
              print *,  avgprobn, avgprobs,avgprobp, numn, nums, nump
              print *, "=======kinematic variables in ", heli,"helicity=========loop====",iloop
              print *, "helicity, xb, Q2, -t, phi, xsec_born, xsec_obs"
              print *, heli, xbd, Q2d, del2d, phigd, xsec_born/2, xsec_obs/2
c
c
c
            else 
               print *, "check xb and q2"
            endif
           ENDDO
        ENDDO
              print *, "=======RC factors in loop====",iloop
              print *, "A_born, A_obs, delta_u, delta_p, delta_A"
              print *, del_born/sum_born, del_obs/sum_obs, sum_obs/sum_born, del_obs/del_born, (del_obs/del_born)/(sum_obs/sum_born)-1
              A_obs = A_obs + del_obs/sum_obs/10
              dsigma_obs = dsigma_obs + sum_obs/10
              dsigma_born = dsigma_born + sum_born/10
              delta_u = delta_u + sum_obs/sum_born/10
              delta_p = delta_p + del_obs/del_born/10
              delta_A = delta_A + ((del_obs/del_born)/(sum_obs/sum_born)-1)/10

              unc_A_obs = unc_A_obs+(del_obs/sum_obs)**2/10 ! rms
              unc_dsigma_obs = unc_dsigma_obs + sum_obs**2/10 ! rms
              unc_dsigma_born = unc_dsigma_born + sum_born**2/10 ! rms
              unc_delta_u = unc_delta_u+(sum_obs/sum_born)**2/10 ! rms
              unc_delta_p = unc_delta_p+(del_obs/del_born)**2/10 ! rms
              unc_delta_A = unc_delta_A+((del_obs/del_born)/(sum_obs/sum_born)-1)**2/10 ! rms
      ENDDO
              unc_A_obs = sqrt(unc_A_obs - A_obs**2)
              unc_dsigma_obs = sqrt(unc_dsigma_obs - dsigma_obs**2)
              unc_dsigma_born = sqrt(unc_dsigma_born - dsigma_born**2)
              unc_delta_u = sqrt(unc_delta_u - delta_u**2)
              unc_delta_p = sqrt(unc_delta_p - delta_p**2)
              unc_delta_A = sqrt(unc_delta_A - delta_A**2)

              print *, "=======RC factors averaged over 10 loops============"
              print *, "A_born, A_obs, delta_u, delta_p, delta_A"
              print *, "averaged:", del_born/sum_born, A_obs, delta_u, delta_p, delta_A
              print *, "unc_A_born, unc_A_obs, unc_delta_u, unc_delta_p, unc_delta_A"
              print *, 0, unc_A_obs, unc_delta_u, unc_delta_p, unc_delta_A
              print *, "Frad_with_error", delta_u, unc_delta_u
              print *, "differential_x_sec", dsigma_obs/2, unc_dsigma_obs/2, dsigma_born/2, unc_dsigma_born/2
        return
        end


        SUBROUTINE  printdstot 
C module for reweighting        
        implicit none
#include "dvcs.inc"
C..
        integer nx,nq,nt,nfe,nfg
        real  xmax,xmin,q2max,q2min,tmax,tmin 
        common /xmax/  xmax,xmin,q2max,q2min,tmax,tmin 
     6                 ,nx,nq,nt,nfe,nfg
C..
C..
      double precision xbd,Q2d,del2d,phield,phigd,dstot,stotint,yb,w2
     6                ,avgdstot,avgprobn,avgprobs,avgprobp,xsec_born,xsec_obs,sum_born,sum_obs,unc_obs,var_obs
     6                ,del_born, del_obs, delu,delp, A_obs, delta_u, delta_p, delta_A
     6                , unc_A_obs, unc_delta_u, unc_delta_p, unc_delta_A, dsigma_born, dsigma_obs, unc_dsigma_born, unc_dsigma_obs
      integer la0,itrial,iloop,numn,nump,nums
c internal variables
c
        logical outoflimits
        real altime,alltime(2)
        real*8 t1lim, t2lim

        helpi = 0 ! unpol target
c
          xbd=xmin!+ix*dx
          Q2d=q2min!+iq*dq
           yb=Q2d/(2D0*Mp*xbd*cl_be)
           w2=Mp*Mp+2.0*Mp*cl_be*yb-Q2d
c
         if( .not.outoflimits(xbd,Q2d) ) then
c
           t2lim=-0.5d0*((q2d/xbd-q2d)*(q2d/xbd+sqrt((q2d/xbd)**2+4d0*Mp*Mp*q2d))+2d0*Mp*Mp*q2d)/w2
           t1lim=-0.5d0*((q2d/xbd-q2d)*(q2d/xbd-sqrt((q2d/xbd)**2+4d0*Mp*Mp*q2d))+2d0*Mp*Mp*q2d)/w2
           tmax = min(-cl_tmin, t1lim)
           tmin = max(-cl_tmin, t2lim)
           if (tmax<tmin) then
              print *, "check -t limit"
              return
           endif
           del2d=tmin!+it*dt
           phigd=cl_phigd

           if (cl_heli.eq.1) then
              heli = 1 ! polarized beam, (+) heli
              CALL bmkxsec(xbd, Q2d, del2d, phield,phigd,dstot)
              xsec_born = dstot
              xff = xbd
              q2ff = Q2d
              tff = del2d
              phiff = pi-phigd
              ebeamff = cl_beam_energy
              CALL bhradgen(ebeamff, xff, q2ff, tff, phiff,cl_vv2cut,cl_delta,egamma,thetag,phig,ichannel,cl_pol,ikeygene,dstot,probn,probp,probs)
              dstot = dstot/2/pi
              print *, "xsec_obs,xsec_born:", dstot, xsec_born
           else if (cl_heli.eq.-1) then

              heli = -1 ! polarized beam, (-) heli
              CALL bmkxsec(xbd, Q2d, del2d, phield,phigd,dstot)
              xsec_born = dstot
              xff = xbd
              q2ff = Q2d
              tff = del2d
              phiff = pi-phigd
              ebeamff = cl_beam_energy
              CALL bhradgen(ebeamff, xff, q2ff, tff, phiff,cl_vv2cut,cl_delta,egamma,thetag,phig,ichannel,cl_pol,ikeygene,dstot,probn,probp,probs)
              dstot = dstot/2/pi
              print *, "xsec_obs,xsec_born:", dstot,  xsec_born
            else
               print *, "check helicity"
               return
            endif
         else 
            print *, "check xb and q2"
         endif

        return
        end


        subroutine printgpd
        IMPLICIT NONE
#include "dvcs.inc"
#include "dvcsmom.inc"
        integer nx,nq,nt,nfe,nfg
        real  xmax,xmin,q2max,q2min,tmax,tmin 
        common /xmax/  xmax,xmin,q2max,q2min,tmax,tmin 
     6                 ,nx,nq,nt,nfe,nfg
      double precision xbd,Q2d,del2d,phield,phigd,dstot,stotint,yb,w2
      double precision nu,W,qmod,E1cm,P1cm,E2cm,P2cm,del2max,del2minc
c
      real  h1uimag,h1dimag,h1tuimag,h1tdimag
      common /h1imag/h1uimag,h1dimag,h1tuimag,h1tdimag
      double precision xb, ybc, Q2, del2, del2min,phip,phipel,P1,P2,Kfac,Jfac,ds
      common/todvcs/ xb, ybc, Q2, del2, del2minc,phip,phipel,P1,P2,Kfac,Jfac,ds
c
c
        real dx,dq,dt,cpi,alu,rP1,rP2
        integer ix,iq,it
        logical outoflimits
      double precision F1pn,F2pn
      common/formfac/ F1pn(2), F2pn(2)
      double precision skew,H1_RE, H1_IM, H1T_RE, H1T_IM,
     &             E1_RE, E1_IM, E1T_RE, E1T_IM , C_I_im
c
        print *, '            x         Q2        t      A_LU      F1        F2       H1-IM
     6      H1-RE     H1T-IM     E1-IM    E1T-IM, H1-RE,    H1T-RE     E1-RE    E1T-RE h1uim h1dim h1tuim h1tdim c_I_im'
        cpi=-3.1415926
        dx=0.05
        dq=0.5
        dt=0.1
        nx=(xmax-xmin)/dx+1
        nq=(q2max-q2min)/dq+1
        nt=(tmax-tmin)/dt+1
c
         DO ix=1,nx
         xbd=xmin+(ix-1)*dx
          DO iq=1,nq
          Q2d=q2min+(iq-1)*dq
           yb=Q2d/(2D0*Mp*xbd*cl_be)
           w2=Mp*Mp+2.0*Mp*cl_be*yb-Q2d
           nu  = Q2d/(2D0*Mp*xbd)
           W   = sqrt(W2)
           qmod = sqrt(nu**2 + Q2d)
c           print *,'xbd,Q2d',xbd,Q2d,outoflimits(xbd,Q2d),nx,nq,nt,tmin,tmax
         if( .not.outoflimits(xbd,Q2d) ) then
c
             DO it=1,nt
              del2d=tmax-(it-1)*dt
              CALL bhdvcs(xbd, Q2d, del2d, 2.2689D0, 2.2689D0,rP1,rP2)
c
              E1cm = Mp*(Mp + nu)/W
              P1cm = Mp*qmod/W
              E2cm = (W2 + Mp**2)/(2D0*W)
              P2cm = (W2 - Mp**2)/(2D0*W)
              del2max = 2D0*(Mp**2 - E1cm*E2cm - P1cm*P2cm)
              del2min = 2D0*(Mp**2 - E1cm*E2cm + P1cm*P2cm)
c            asymmetries at 90-degree
             alu=hs1Iunp/(hc0BH+hc0dvcs+hc0Iunp)
c           print *,'xbd,Q2d,del2',xbd,Q2d,del2d,del2min,del2max
      if( del2d.le.del2min .and. del2d.ge.del2max ) then  ! in limits
      call nuclFF( del2d )
      skew = xbd/(2D0 - xbd)
      IF(IGPD.lt.100) then
      call amptab(skew, del2d, 
     &             H1_RE, H1_IM, H1T_RE, H1T_IM,
     &             E1_RE, E1_IM, E1T_RE, E1T_IM )
      ELSE
      call amptabgag(skew, del2d, 
     &             H1_RE, H1_IM, H1T_RE, H1T_IM,
     &             E1_RE, E1_IM, E1T_RE, E1T_IM )
      endif
      C_I_im = F1pn(1)*H1_IM + xbd/(2D0-xbd)*(F1pn(1)+F2pn(1))*H1T_IM 
c     &- del2/Mp**2/4D0*F2pn(1)*E1_IM
c
c
        write(6,'(A,18F10.3)') 'BM-GPDs ',xbd,q2d,-del2d,alu,F1pn(1), F2pn(1),
     &             -H1_IM,-H1_RE, -H1T_IM,
     &             -E1_IM, -E1T_IM,H1_RE,H1T_RE,E1_RE,E1T_RE,C_I_im
        write(6,'(A,18F10.5)') 'BM-GPD2 ',del2minc,P1,P2,Kfac,Jfac,ds,
     &hc0BH,hc1BH/hc0BH,hc2BH/hc0BH,hs1Iunp/hc0BH,hc0Iunp/hc0BH,hc0dvcs/hc0BH
        write(6,'(A,18F10.5)') 'BM-MOM  ',ac0int,ac1int,as1int,ac0dvcs,ac1dvcs
        endif ! if( del2d.....
        enddo
        endif  !if( .not.outoflimits
        enddo
        enddo
c
        print *,'STOP-here'
        return
        end
c
        subroutine printusage
        print *,'dvcspi0gen [options]            '
        print *,'      option  value      default    comment  '
        print *,'     --pi0                    exclusive pi-0 on '
        print *,'     --eta                    exclusive eta on '
        print *,'     --nodvcs                 DVCS off '
        print *,'     --v    verbos_level    0  additional printout '
        print *,'     --trig nevents  10      number of triggers '
        print *,'     --seed seed_number 0-use local time '
        print *,'     --nmax nevents   40k     # of triggers per file '
        print *,'     --beam beam_energy   5.754 e- momentum in GeV'
        print *,'     --ycol P_1 cut        0.005      P_1>ycol_cut'
        print *,'     --x  x-min   x-max   0.1 0.65   min max x_Bj '
        print *,'     --q2 Q2-min Q2-max 1.0  10. min max Q2 in GeV^2 '
        print *,'     --y y-min y-max 0.05 0.9    min max for y=nu/E'
        print *,'     --w w^2-min         4.0 min for w^2=M^2+2MyE-Q^2'
        print *,'     --t tmin tmax  0 1.0      t  min/max in GeV^2 '
        print *,'     --th thmin thmax  0.2 1 theta min/max \theta rad'
        print *,'     --xpos x-position      0 beam x position in cm '
        print *,'     --ypos y-position      0 beam y position in cm '
        print *,'     --zpos z-position      0 target z position in cm '
        print *,'     --zwidth z-width 0  width z in cm (zpos+/-zwidth)'
        print *,'      --raster diameter 0.75   raster diameter in cm  '
        print *,'     --weight   flat distributions with weight(part12)'
        print *,'     --printgpd               print gpds and exit '
        print *,'     --nont               do not write out the ntuple '
        print *,'     --file              dvcspi0gen   filename'
        print *,'     --gpd  Igpd 3  GPD model(1-A,2-B,3-C,4-D)>=100VGG'
        print *,'     --scale  scale      1.0   scale the sinphi-mom '
        print *,'     --targ target       proton   deut/neut possible '
        print *,'     --lpol                    Long.pol.target  '
        print *,'     --tpol                    Trans.pol.target  '
        print *,'     --writef    format      0-lund12, 1-lundgsim'
        print *,'     --mod     write-mod      0-all, 1-cut on events'
        print *,'     --mom                include moments in ntuple'
        print *,'     --proloss                  add proton loss'
        print *,'     --ktcor          FALSE   turn on k_t cor for A_LU'
        print *,'     --radgen                   include bhradgen'
        print *,'     --radstable                use born cross sections for rejection sampling'
        print *,'     --nodat               do not write a data file'
        print *,'     --acce16            include e16 acceptance for e-'
        print *,'     --acceg1           include eg1 acceptance for e-'
        print *,'     --acc12         include clas++ acceptance for e-'
        print *,'     --smear                   smear moments'
        print *,'     --A    value         0.006   smear param-A'
        print *,'     --B    value         0.001  smear param-B'
        print *,'     --C    value         0.0008  smear param-C'
        print *,'     --D    value         0.001  smear param-D'
        print *,'     --nmax   value     2000  maximum events per file'
        print *,'     --print nprint     1000   print ev nprint event '
        print *,'     --bh  value      3 BH status:3-All, 1-only BH  ' 
        print *,'     --delta  value      0.01 Minimum rad photon energy (GeV)'
        print *,'     --vv2cut value      0.1 cuts on missing mass ep squared (GeV^2)'
        return
        end
c
        subroutine printvalues
        implicit none
#include "dvcs.inc"
        if(cl_pi0) cl_dvcs=.FALSE.
        if(cl_eta) cl_dvcs=.FALSE.
        print *,'dvcspi0gen input values:'    
        print *,'generators:       DVCS->',cl_dvcs,'    pi-0->',cl_pi0     
        print *,'NEVENTS:             ',cl_triggers    
        print *,'BEAM ENERGY  :         ',cl_beam_energy  
        print *,'BH status:3-All, 1-only BH  :  ',cl_bh  
        print *,'w^2   min    :         ',cl_wmin  
        print *,'y   min/max:         ',cl_ymin,cl_ymax  
        print *,'Eprime thmin/thmax in rad: ',cl_thmin,cl_thmax  
        print *,'-Delta  tmin/tmax:         ',cl_tmin,cl_tmax  
        print *,'Q2  q2min/q2max:         ',cl_q2min,cl_q2max  
        print *,'x  xmin/xmax:         ',cl_xmin,cl_xmax  
        print *,'TARGET:              ', cl_target 
        print *,'BEAM x-position in cm:   ', cl_xpos 
        print *,'BEAM y-position in cm:   ', cl_ypos 
        print *,'TARGET z-position in cm:   ', cl_zpos 
        print *,'TARGET 1/2 width in cm:   ', cl_zwidth 
        print *,'raster diameter(cm)  :   ', cl_rast
        print *,'Tpol (0-U,1-L,2-T):   ', cl_pol 
        print *,'GPD model A,B,C,D:   ', cl_gpd 
        print *,'GPD sinphi scale:   ', cl_scale 
        print *,'BHRADGEN ON        :   ', cl_radgen 
        print *,'Use Born Cross sections when bhradgen on :   ', cl_radstable
        print *,'y_col cut        :   ', cl_ycol 
        print *,'acceptance       :   ',acceg1,acce16,acc12
        print *,'nmax-events       :   ',cl_nmax
        print *,'write form      :   ',cl_writef
        print *,'write mod       :   ',cl_mod
        print *,'Minimum rad photon energy:  ',cl_delta
        print *,'cuts on missing mass ep (GeV^2):  ',cl_vv2cut
        print *,'trento angle to be used at the printing rc factor:  ',cl_phigd
        print *,'write file_name  :   ',bosout
        if(cl_mom) print *,'include moments in the ntuple:   ', cl_mom 
        if(cl_smear.eq.1) 
     6print *,'smearing A,B,C,D:   ', cl_sma,cl_smb,cl_smc,cl_smd
     6,cl_proloss  
        print *,'__________________________________'
        print *,'verbosity level',  cl_verblev
        return
        end

C******************************************************************


        SUBROUTINE DVCSINIT
        IMPLICIT NONE

C*******************************************************************
#include "dvcs.inc"

C.. Set up cuts and parameters which are different from the default
C.. valules 
c
C..
        integer nx,nq,nt,nfe,nfg
        real  xmax,xmin,q2max,q2min,tmax,tmin 
        common /xmax/  xmax,xmin,q2max,q2min,tmax,tmin 
     6                 ,nx,nq,nt,nfe,nfg
C..
c
          CHARACTER*2 CTARGET
          integer lst40
c
C..
        real  elmax,elmin,ttmax,ttmin
C..
        Ifl=1   !dsigma/dx/dq/dt/dphiel/dpigam
        Ivar=cl_bh  ! BH+DVCS if 3
        IGPD=cl_gpd  ! GPD type
        cl_be=cl_beam_energy
        Ich=-1
C.. 
        nx=40
        nq=20
        nt=40
        if(cl_pol.eq.2) then 
            nfe=36
        else
            nfe=1
        endif
        nfg=36 
C.. $x$ min cut  
       xmin=cl_xmin
C.. $x$ max cut  
       xmax=cl_xmax
C.. $Q^2$ max cut
       q2max=cl_q2max
C.. $Q^2$ min cut    
       q2min=cl_q2min
C.. t max cut
       tmax=-cl_tmin
C.. t min cut    
       tmin=-cl_tmax
c
C.. lower limit scat. lep. energy
       elmin=(1.0-cl_ymin)*cl_be
C..  upper limit scat. lep. energy
       elmax=(1.0-cl_ymax)*cl_be
c
C.. lower limit scat lep. angle
       ttmin=cl_thmin !0.52-0.08          ! 29.75 degree+/-4
C..  upper limit scat lep. angle
       ttmax=cl_thmax !0.52+0.08
c
C.. Momentum of incoming lepton (in GeV)
       Ed = cl_be  

C.. Simulation for a proton
      if(cl_target.eq.'neutron') then
        Ipn = 2
        ctarget='P'
      elseif(cl_target.eq.'proton') then
        Ipn = 1
        ctarget='P'
      else               ! proton
          call printusage
          stop       
      endif
C
C       get the maximum 5-dim x-section
C

        if(cl_printdstot) then 
        call printdstot
        stop
        endif

        if(cl_printrad) then 
        call printrad
        stop
        endif

        if(cl_printgpd) then 
        call printgpd
        stop
        endif
        print *,'CALL getSCALE-->',Ipn,xmax,xmin,tmax,tmin
        CALL getSCALE
c
        if(cl_radgen) then
c         lst40=0  ! unpolarized case
c         beam_energy=cl_beam_energy
c         call radgen_init(ctarget,beam_energy,lst40)
         ikeygene = 4


        endif 
c
        RETURN
        END
c
c
        SUBROUTINE  getSCALE
        implicit none
#include "dvcs.inc"
C..
        integer nx,nq,nt,nfe,nfg
        real  xmax,xmin,q2max,q2min,tmax,tmin 
        common /xmax/  xmax,xmin,q2max,q2min,tmax,tmin 
     6                 ,nx,nq,nt,nfe,nfg
C..
C..
      double precision xbd,Q2d,del2d,phield,phigd,dstot,stotint,yb,w2
c internal variables
c
        real dx,dq,dt,dfe,dfg,cix,cit,ciq,cife,cifg
        integer ix,iq,it,ife,ifg,la0
        logical outoflimits
        real altime,alltime(2)
        real*8 t1lim, t2lim
c
        smax=0.0
        la0=0
        dx=(xmax-xmin)/nx
        dq=(q2max-q2min)/nq
        dfe=2*pi/nfe
        dfg=2*pi/nfg
        if(cl_pol.ge.1) la0=-1
       if(cl_verblev.gt.50) print *,'loop limits',nx,xmin,dx,nq,q2min,dq,nt,tmin,dt
          DO helpi=la0,1,2
            DO heli=-1,1,2
         stotint=0
         DO ix=1,nx
         xbd=xmin+ix*dx
         cix=1.0
         if(ix.eq.1.or.ix.eq.nx) cix=0.5
          DO iq=1,nq
          Q2d=q2min+iq*dq
          ciq=1.0
          if(iq.eq.1.or.iq.eq.nx) ciq=0.5
           yb=Q2d/(2D0*Mp*xbd*cl_be)
           w2=Mp*Mp+2.0*Mp*cl_be*yb-Q2d
c
         if( .not.outoflimits(xbd,Q2d) ) then
c
             DO it=1,nt
                 t2lim=-0.5d0*((q2d/xbd-q2d)*(q2d/xbd+sqrt((q2d/xbd)**2+4d0*Mp*Mp*q2d))+2d0*Mp*Mp*q2d)/w2
                 t1lim=-0.5d0*((q2d/xbd-q2d)*(q2d/xbd-sqrt((q2d/xbd)**2+4d0*Mp*Mp*q2d))+2d0*Mp*Mp*q2d)/w2
                 tmax = min(-cl_tmin, t1lim)
                 tmin = max(-cl_tmax, t2lim)
                 if (tmax<tmin) then
                    continue
                 endif
              dt=(tmax-tmin)/nt
              del2d=tmin+it*dt
              cit=1.0
             if(it.eq.1.or.it.eq.nx) cit=0.5
              DO ife=1,nfe
               phield=ife*dfe
               cife=1.0
               if(ife.eq.1.or.ife.eq.nx) cife=0.5
               DO ifg=1,nfg
                phigd=ifg*dfg
                cifg=1.0
                if(ifg.eq.1.or.ifg.eq.nx) cifg=0.5
c
             if(cl_dvcs) then
                if(cl_radgen.and.(.not.cl_radstable)) then
                    xff = xbd
                    q2ff = Q2d
                    tff = del2d
                    phiff = pi-phigd
                    ebeamff = cl_beam_energy
                    ikeygene = 4
                    CALL bhradgen(ebeamff, xff, q2ff, tff, phiff,cl_vv2cut,cl_delta,egamma,thetag,phig,ichannel,cl_pol,ikeygene,dstot,probn,probp,probs)
                    dstot = dstot/2.0/pi
                else
                    CALL bmkxsec(xbd, Q2d, del2d, phield,phigd,dstot)
                endif
             else if(cl_pi0.or.cl_eta) then
              CALL excpi0(xbd, Q2d, del2d, phield,phigd,dstot)
             else
              print *,'ERROR! at least one reaction should be ON'
              stop
             endif
c
             stotint=stotint+dstot*cix*cit*ciq*cife*cifg
c
           if(smax.lt.dstot) then
            smax=dstot
            if(cl_verblev.gt.1) print *,
     6'smax',heli,helpi,xbd,q2d,del2d,stotint,dstot,smax

           endif
c
       if(cl_verblev.gt.50) print *,'loopall',ix,iq,it,ife,ifg,dstot
            ENDDO
           ENDDO
        ENDDO
       altime=etime(alltime)
       if(cl_verblev.gt.50) print *,'loopin getscale',ix,iq,it,ife,ifg,dstot
     6,altime,alltime(1),alltime(2)
       else
c        print *,'bad-y',xmin,xmax,xbd,q2min,q2max,q2d,del2d
       endif
       ENDDO
       ENDDO
       if(cl_verblev.gt.-1) print *,'Maximum found at',SMAX,
     6'x-sec(pb) = ',stotint*dx*dt*dq*dfe*dfg,dx*dt*dq*dfe*dfg
          ENDDO
         ENDDO
        return
        end            

        SUBROUTINE  genEVENT(I)      
        implicit none
#include "dvcs.inc"
#include "ntupgdvcs.inc"
#include "tailcom.inc"
ccc-Mick
	common/radmick/PhRAD
ccc-Mick
ccc- Hyon-Suk
        common/radhs/radtheta, ichannel
        real radtheta
ccc- Hyon-Suk

C..
C..
        integer nx,nq,nt,nfe,nfg
        real  xmax,xmin,q2max,q2min,tmax,tmin 
        common /xmax/  xmax,xmin,q2max,q2min,tmax,tmin 
     6                 ,nx,nq,nt,nfe,nfg
C..     generate a DVCS event
      double precision xbd,Q2d,del2d,phield,phigd,dstot,yb,w2
c
c       radiative stuff
c
        real VPGEN(4),VPRAD(4),PhRAD(4),q2tr,utr,weight
c
        REAL PolarTheta, Azimuthalphi
        REAL VEC(8),e_photon,prob_max,prob_curr,ggam2p
        INTEGER I,ITRY,LEN
        logical outoflimits
        real  beam_energy
        real*8 t1lim, t2lim
c
        beam_energy=cl_beam_energy
        LEN=8
        ITRY=0
 10     CALL MYRANMAR(VEC,LEN)
        ITRY=ITRY+1
        xbd=xmin+(xmax-xmin)*VEC(1)
        Q2d=q2min+(q2max-q2min)*VEC(2)

C         if(cl_radgen) then
C ccc- Hyon-Suk
C          PhRAD(1)=0.
C          PhRAD(2)=0.
C          PhRAD(3)=0.
C          PhRAD(4)=0.
C ccc- Hyon-Suk         
C          radq2=Q2d
C          radnu=Q2d/2.0/Mp/xbd
C          if( outoflimits(xbd,Q2d) ) then
C           if(itry .gt.10000.and.cl_verblev.eq.1000) 
C      6print *,'stuck before radgen:',itry,istatus,xbd,q2d
C           goto 10 
C          endif
C          call getgammastar(xbd,Q2d,phield,VPGEN)
C  !	 call RADGEN(beam_energy,VPGEN,VPRAD,PhRAD,q2tr,utr,weight)
C          Q2d=q2tr
C          xbd=Q2d/2.0/Mp/utr
C          if( outoflimits(xbd,Q2d) ) then
C           if(itry .gt.10000.and.cl_verblev.eq.1000) 
C      6print *,'stuck after radgen:',itry,istatus,xbd,q2d
C           goto 10 
C          endif
C c
C          ggam2p=PhRAD(1)**2+PhRAD(2)**2+PhRAD(3)**2
C          PhRAD(4)=sqrt(ggam2p)
C          if(ggam2p.gt.0) then
C            gradt=PolarTheta(PhRAD(1),PhRAD(2),PhRAD(3))
C            gradf=Azimuthalphi(PhRAD(1),PhRAD(2))
C            grade=PhRAD(4)
C           else 
C            gradt=0
C            gradf=0
C            grade=0
C          endif   
C            radweight=weight !del2d
C            gelebits=ita
C ccc- Hyon-Suk
C          radtheta=gradt
C ccc- Hyon-Suk
C         endif


        phield=2*pi*VEC(4)
        phigd=2*pi*VEC(5)
        yb=Q2d/(2D0*Mp*xbd*cl_beam_energy)
        w2=Mp*Mp+2.0*Mp*cl_beam_energy*yb-Q2d
                 t2lim=-0.5d0*((Q2d/xbd-Q2d)*(Q2d/xbd+sqrt((Q2d/xbd)**2+4d0*Mp*Mp*Q2d))+2d0*Mp*Mp*Q2d)/w2
                 t1lim=-0.5d0*((Q2d/xbd-Q2d)*(Q2d/xbd-sqrt((Q2d/xbd)**2+4d0*Mp*Mp*Q2d))+2d0*Mp*Mp*Q2d)/w2
                 tmax = min(-cl_tmin, t1lim)
                 tmin = max(-cl_tmax, t2lim)

        del2d=tmin+(tmax-tmin)*VEC(3)
ccc- Hyon-Suk
C         if(radtheta.gt.0.15) then
C          yb=Q2d/(2D0*Mp*xbd*cl_beam_energy)
C          w2=Mp*Mp+2.0*Mp*cl_beam_energy*yb-Q2d
C         else
C          yb=Q2d/(2D0*Mp*xbd*(cl_beam_energy-PhRAD(4)))
C          w2=Mp*Mp+2.0*Mp*(cl_beam_energy-PhRAD(4))*yb-Q2d
C         endif
ccc- Hyon-Suk
        if(yb.gt.cl_ymax.or.w2.lt.cl_wmin) then 
             goto 10
         endif
c
c       random helicities of beam and target
c
        heli=1
        if(VEC(7).gt.0.5) heli=-1
        if(cl_pol.gt.0) then
         helpi=1
         if(VEC(8).gt.0.5) helpi=-1
        else
         helpi=0
        endif
c
        if(cl_dvcs) then
C..     Generation of an event in GENDVCS 
            if (cl_radgen) then
                    xff = xbd
                    q2ff = Q2d
                    tff = del2d
                    phiff = pi-phigd
                    ebeamff = cl_beam_energy
                    ikeygene = 4
                    CALL bmkxsec(xbd, Q2d, del2d, phield,phigd,dstot)
                    bornweight = dstot
                    CALL bhradgen(ebeamff, xff, q2ff, tff, phiff,cl_vv2cut,cl_delta,egamma,thetag,phig,ichannel,cl_pol,ikeygene,dstot,probn,probp,probs)
                    dstot = dstot/2.0/pi
                    radweight = dstot
                    radq2 = Q2d
                    radxb = xbd
                    radt = del2d
                    radphi = phigd
                    iproctype = ichannel
                    PhRAD(4) = egamma
            else
                CALL bmkxsec(xbd, Q2d, del2d, phield,phigd,dstot)
            endif
        else
        CALL excpi0(xbd, Q2d, del2d, phield,phigd,dstot)
        endif
c
       if(itry.GE.100000)  print *,'Check the limits',itry,smax
     & ,heli,xbd, Q2d, del2d, phield, phigd, dstot
        if(weightOK) then
          if(dstot.le.1.0E-8.or.istatus.ne.0) goto 10 !put in the weight 
        else
          if(dstot.le.smax*VEC(6)) goto 10 ! use x-section
        endif
        if(cl_dvcs) then
           call getphoton(xbd,Q2d,del2d,phield,phigd)
        else
           call getpi0(xbd,  Q2d, del2d,phield,phigd)
        endif
       if(cl_verblev.gt.100) print *,'hel,x,q2,-t,tot'
     & ,heli,helpi,xbd, Q2d, del2d,phigd,dstot
       gelewgt=dstot/smax
       gprox=xbd
       gproy=phield
       gproz=0
       return
       end

C**********************************************************************
C      random number generator
C**********************************************************************
         real function random_num()
c         CALL RANDOM_SEED()
         CALL RANDOM_NUMBER(random_num)
         return
         end


C***********************************************************************
      SUBROUTINE CLEANNT
C***********************************************************************
#include "ntupgdvcs.inc"
#include "dvcsmom.inc"
#include "ludat1234.inc"
       call VZERO(gnele,12)
       call VZERO(gnpi0,30)
       call resetmom()

       return
       end
c
C***********************************************************************
      SUBROUTINE CLASFILL(IEVENT)
C***********************************************************************
#include "dvcs.inc"
#include "dvcsmom.inc"
#include "ntupgdvcs.inc"
#include "ludat1234.inc"
ccc-Mick
	common/radmick/PhRAD
ccc-Mick
ccc- Hyon-Suk
        common/radhs/radtheta
        real radtheta
ccc- Hyon-Suk

        real PhRAD(4)
	real theta_radgam,phi_radgam
        REAL PolarTheta, Azimuthalphi
c
c
        real ppmod2
        real pmy(4000,5)
        integer kmy(4000,5)
c
        integer ierr,lunout
        INTEGER j,jj
        real pcharge(10)
c
C..     fill lund
c
c    electron
       DO j=1,3
        pmy(1,j)=V3k2(j)
        pmy(2,j)=V3p2(j)
       ENDDO
c 
        pmy(1,5)=mele ! mass of el
        pmy(2,5)=Mp    !mass of proton
        pmy(1,4)=sqrt(mele**2+pmy(1,1)**2+pmy(1,2)**2+pmy(1,3)**2)
        ppmod2 =pmy(2,1)**2+pmy(2,2)**2+pmy(2,3)**2
        if(ppmod2.gt.0.and.ppmod2.lt.1000.0) then
         pmy(2,4)=sqrt(Mp**2+ppmod2)
        else
c             dump this events
         return
        endif
c
        pcharge(1)=-1
        pcharge(2)=1
c
       kmy(1,2)=11  ! electron
       kmy(2,2)=2212 ! proton
       kmy(1,1)=1    ! final
       kmy(2,1)=1    ! final
c       
       do j=1,N-1
        do l=1,5
         pmy(2+j,l)=p(j+1,l)  ! j-is the decaying particle
         kmy(2+j,l)=k(j+1,l)
        enddo
       enddo

       jj=N+1
c
        IF( phrad(4).gt.0.010 ) THEN
           jj=jj+1
           theta_radgam=PolarTheta(PhRAD(1),PhRAD(2),PhRAD(3))
           phi_radgam=Azimuthalphi(PhRAD(1),PhRAD(2))
           pmy(jj,1)=phrad(4)*sin(theta_radgam)*cos(phi_radgam)
           pmy(jj,2)=phrad(4)*sin(theta_radgam)*sin(phi_radgam)
           pmy(jj,3)=phrad(4)*cos(theta_radgam)
           pmy(jj,4)=phrad(4)
           kmy(jj,1)=1
           kmy(jj,2)=22 ! radiative photon
           pcharge(jj)=0
        ENDIF

c
         do  j=1,jj
          call fillevent(j,pmy,kmy)
         enddo
c
          call goodevent()     ! check the event
c
          if(dvcsgenOK) then
c            print *,'write-nt',gelee,ntOK
            if (datfileOK) call  writefile(jj,pmy,kmy)
#ifdef  hbook
            if(ntOK) call HFNT(22)
#endif
          endif
c
       return
       end  

C***********************************************************************
c write compatible with kine_lund format
c
       subroutine writefile(jj,p,k)
       implicit none
#include "ntupgdvcs.inc"
#include "dvcs.inc"
        integer i,ifile,itotal,iparl1,iparl2
        logical lFirst
        DATA lFirst /.TRUE./
        save lFirst
          DATA ifile /0/
          DATA itotal /20000/
          character*21  fname
c
        integer jj,j,jl,k(4000,5)
        real gy,gw2,u,p(4000,5),rhelpi,rheli
        real randnum,random_num
        real PARL(2),plu(4),V(4)                   ! for LUND-compatibility
        data PARL/1.0,1.0/
        data plu/-1.0,1.0,0.0,0.0/
        data v/0.0,0.0,0.0,0.0/
c
          if((itotal.eq.cl_nmax).or.lFirst) then
           close (41)
           lFirst=.TRUE.
           ifile=ifile+1
           itotal=0
          endif
c
          if (lFirst) then   ! open the file
           if(cl_docker) then
             write(fname,'(A)') 'dvcsgen.dat'
             cl_nmax=cl_triggers+100
           else
           if(ifile.le.9) then
             write(fname,'(A,i1,A)') bosout,ifile,'.dat'
             write(6,'(A,i1,A)') bosout,ifile,'.dat'
           else if(ifile.ge.10.and.ifile.le.99) then
            write(fname,'(A,i2,A)') bosout,ifile,'.dat'
            write(6,'(A,i2,A)') bosout,ifile,'.dat'
           else if(ifile.ge.100.and.ifile.le.999) then
            write(fname,'(A,i3,A)') bosout,ifile,'.dat'
            write(6,'(A,i3,A)') bosout,ifile,'.dat'
           else
            write(fname,'(A,i4,A)') bosout,ifile,'.dat'
            write(6,'(A,i4,A)') bosout,ifile,'.dat'
           endif 
           endif ! docker check          
             OPEN(unit=41,file=fname,status='unknown')
             lFirst=.FALSE.
          endif
c
                itotal=itotal+1
c
           gy=1.0-gelee/cl_beam_energy
           gw2=gw*gw
           u=heli+10*helpi
            randnum=random_num()-0.5
 9         v(1)=(cl_rast*(random_num()-0.5)+cl_xpos)
           v(2)=(cl_rast*(random_num()-0.5)+cl_ypos)
      if(sqrt((v(1)-cl_xpos)**2+(v(2)-cl_ypos)**2).gt.cl_rast*0.5) goto 9 
           v(3)=cl_zpos+randnum*cl_zwidth
          gelex=v(1)
          geley=v(2)
          gelez=v(3)
c
          iparl1=PARL(1)
          iparl2=PARL(2)
          rhelpi=helpi
          rheli=heli
c101      FORMAT(2x,I10,2I6,2F11.2,I5,F10.3,2I5,2E14.7,2F6.2)
         if(cl_writef.eq.0) then                  ! lund12 format
                write(41,111) jj,iparl1,iparl2,rhelpi,rheli,11,
     6                   cl_beam_energy,1,1,gsin2
            do  j=1,jj
                 P(j,4)=sqrt(P(j,1)**2+P(j,2)**2+P(j,3)**2+P(j,5)**2)
                 write (41,112) j,plu(j),1,(k(j,jl),jl=2,4),(P(j,jl),jl=1,5)
     6                                                   ,(V(jl),jl=1,3)
            enddo
          
         else if(cl_writef.eq.1) then                  ! gsim-lund format
           write(41,102) jj,gxb,gy,gw2,gq2,u
            do  j=1,jj
              write(41,101) (k(j,jl),jl=1,4),heli
              write(41,103) (p(j,jl),jl=1,5) 
              write(41,103) v(1),v(2),v(3),0.0,0.0 
            enddo

         else if(cl_writef.eq.2) then                  ! variation of lund12 format
                write(41,111) jj,iparl1,iparl2,rhelpi,rheli,11,
     6                   cl_beam_energy,1,1,radweight
                 write (41,112) 1,radxb,1,(k(1,jl),jl=2,3),iproctype,(P(1,jl),jl=1,3),radq2,-radt
     6                                                   ,(V(jl),jl=1,3)
                 write (41,112) 2,radphi,1,(k(2,jl),jl=2,4),(P(2,jl),jl=1,3),radxb_true,radq2_true
     6                                                   ,(V(jl),jl=1,3)
                 write (41,113) 3,plu(3),1,(k(3,jl),jl=2,4),(P(3,jl),jl=1,4),bornweight
     6                                                   ,(V(jl),jl=1,3)

                 if (jj.eq.4.) then
                 P(jj,4)=sqrt(P(jj,1)**2+P(jj,2)**2+P(jj,3)**2+P(jj,5)**2)
                 write (41,112) jj,plu(j),1,(k(jj,jl),jl=2,4),(P(jj,jl),jl=1,5)
     6                                                   ,(V(jl),jl=1,3)
                 endif

         else
            print *,'ERROR: unknown write format'
            stop            
         endif
c
 111      FORMAT(2x,I10,2I6,F3.0,2x,F3.0,I5,F10.3,2I5,E14.7)
c 111      FORMAT(2x,I10,2F4.0,2I3,5F7.3)
 112      FORMAT(2x,I3,2x,F10.4,2I5,2I3,5F10.4,2X,5F10.4)
 113      FORMAT(2x,I3,2x,F10.4,2I5,2I3,4F10.4,E14.7,2X,5F10.4)
 101      FORMAT(2x,10I10)
 102      FORMAT(2x,I5,10F10.4)
 103      FORMAT(2x,10F10.4)
        return
        end

C***********************************************************************
      SUBROUTINE goodevent()
C***********************************************************************
       implicit none
#include "ntupgdvcs.inc"
#include "dvcs.inc"
       real pi4(4),qiu4(4),el04(4),elf4(4),tnorm(4)
       real anu,pien,ebeam
       real pro4(4),tnorm2(4)
       real vmass,vangle,vdotm
       real gt,gmismasrho,gphistar,phistar,mismas,thetagg
       real gamel,tchan
       integer daccepte16emcgen,eg1amcaccept,accep11,iaccept,iOK
       dvcsgenOK=.FALSE.
       ebeam=cl_beam_energy
       anu=ebeam-gelee
           if(acceg1) then
        iaccept=eg1amcaccept(cl_beam_energy,gelee,gelet,gelef)
           if(iaccept.eq.0) return
           endif
c
           if(acc12) then
        iaccept=accep11(11,gelee,gelet,gelef,2250.0)
           if(iaccept.eq.0) return
           endif
c
c        print *,'here1'
       if (gnpro.gt.0) then 
          gprotg=gamel(0,1)   ! prot_gamm e
          gprotl=gamel(0,2)   ! prot_gamm e'
          gprotv=gamel(0,3)   ! prot_gamm virt_gamma angle
          gprotmis=mismas(4)
          gprott=tchan(4)
          gprofs=phistar(4)
        endif
c
        if (gngam.gt.0) then
          gproz=thetagg()
        endif 
        if (gnpi0.gt.0) then
          gprox=mismas(0)
        endif

c
c      check the write mod
c
        if(cl_mod.eq.1) then
         if (ggamt.lt.0.0348.or.ggamt.gt.0.087) then   ! photons from 2 to 5 degree
          return
         endif
        endif
c
        if(cl_mod.eq.2) then     ! FX cuts
        iOK=0
         if (0.15.lt.gxb.and.0.2.gt.gxb.and.
     6       (gelet*180/3.1415.gt.26.0).and.gelee.gt.0.8.and.
     6       gprott.gt.0.2.and.gprott.lt.0.25) iOK=1 ! bin 1
         if (0.3.lt.gxb.and.0.4.gt.gxb.and.
     6       (gelet*180/3.1415.gt.21.0).and.(23.5.gt.gelet*180/3.1415)
     6       .and.gprott.gt.0.2.and.gprott.lt.0.25) iOK=1 ! bin 2
         if (0.3.lt.gxb.and.0.4.gt.gxb.and.
     6       (gelet*180/3.1415.gt.21.0).and.(23.5.gt.gelet*180/3.1415)
     6       .and.gprott.gt.0.5.and.gprott.lt.0.80) iOK=1 ! bin 3
          if (iOK.eq.0) return
        endif
c
        if(acce16) then
          iOK=daccepte16emcgen()
          if(iOK.eq.1) dvcsgenOK=.TRUE.  ! write only e16 acceptance
        else
          dvcsgenOK=.TRUE.    ! write everything
        endif
c
         return
         end
c
        real function  PolarTheta(vx,vy,vz)
      implicit none
        real vx,vy,vz,pmod,theta
        pmod=vx*vx+vy*vy+vz*vz
        if(pmod .gt. 0 ) then
         theta=acos(vz/sqrt(pmod))
        else
         theta=-100
        endif
        PolarTheta=theta
        return
        end

        real function Azimuthalphi(vx,vy )
      implicit none
        real pi
        parameter ( PI=3.1415926)
        real vx,vy,pmod,phi,cosf
        pmod=vx*vx+vy*vy
        if(pmod .gt. 0 ) then
         pmod=sqrt(pmod)
         cosf=vx/pmod
        else
         cosf=1.0
        endif
        if(abs(cosf) .le. 1.0) phi=acos(cosf);
        if(vy .lt. 0.0) phi= 2*PI-phi;
        Azimuthalphi=phi
        return
        end
c
        subroutine fillevent(j,p,k)
        implicit none
#include "dvcs.inc"
#include "dvcsmom.inc"
#include "ntupgdvcs.inc"
        integer j,k(4000,5)
        real p(4000,5)
c
       real el04(4),elf4(4),qiu4(4),pro4(4)
       integer part,part_id,parent_id,lund2geantid
       real cx,cy,cz,theta,phi,pmom,pener,px,py,pz
       real Azimuthalphi,PolarTheta,anu,gw2,pien
       real pi0x,pi0y,pi0z
c
       part_id=lund2geantid(k(j,2))
       px=p(j,1)
       py=p(j,2)
       pz=p(j,3)
       pmom=sqrt(px*px+py*py+pz*pz)         ! momentum
       pener=sqrt(p(j,5)*p(j,5)+pmom*pmom)  ! energy
c
       if(pmom.gt.0) then
        cz=pz/pmom
        cx=px/pmom
        cy=py/pmom
       else
        cx=0
        cy=0
        cz=0
       endif
c
        if(cz.le.1.0) then
           theta=acos(cz)
        else
           theta=0
        endif
        phi=Azimuthalphi(px,py) !atan2(cy,cx)
c
       SELECT CASE (part_id)  
        CASE (3) ! electron
         gnele=gnele+1
         if(gnele.eq.1 .or. gelee.lt.pener ) then
           gelee=pener
           geleh=heli+10*helpi
           gelet=theta
           gelef=phi
c
          if(cl_smear.eq.1) call smear(part_id,gelee,gelet,gelef)  ! smear     
c
          anu=cl_beam_energy-gelee
          gq2=2.*cl_beam_energy*gelee*(1.0-cos(gelet))
          gw2=Mp*Mp+2.0*Mp*anu-gq2
c
          if(gw2.gt.0) then
            gw=sqrt(gw2)
          else
            gw=0
          endif
c
          pien=sqrt(cl_beam_energy*cl_beam_energy-mele*mele)
          el04(1)=0
          el04(2)=0
          el04(3)=pien
          el04(4)=cl_beam_energy
c
          pien=sqrt(gelee*gelee-mele*mele)
          elf4(1)=pien*sin(gelet)*cos(gelef)
          elf4(2)=pien*sin(gelet)*sin(gelef)
          elf4(3)=pien*cos(gelet)
          elf4(4)=gelee
          call vdifm(el04,elf4,qiu4,4)
c         print *,'q2',gq2,vdotm(qiu4,qiu4,4)          
          gxb=gq2/2.0/Mp/anu
        endif 
c
        CASE (14) ! proton
         gnpro=gnpro+1
         if(gnpro.eq.1 .or. gproe.lt.pener ) then
          gproe=pener
          gprot=theta
          gprof=phi
          if(cl_smear.eq.1) call smear(part_id,gproe,gprot,gprof)  ! smear 
         endif 
c
        CASE (7) ! pi0
         gnpi0=gnpi0+1
         if(gnpi0.eq.1 .or. gpi0e.lt.pener ) then
          gpi0e=pener
          gpi0t=theta
          gpi0f=phi
         endif 

        CASE (1) ! gamma
         gngam=gngam+1
         if(gngam.eq.1 .or. ggame.lt.pener ) then  ! high energy  gamma-1
c
           if(gngam.gt.1) then   ! if high energy is the second store the first
            ggam2e=ggame
            ggam2t=ggamt
            ggam2f=ggamf
           endif
c
          ggame=pener
          ggamt=theta
          ggamf=phi
          if(cl_smear.eq.1) call smear(part_id,ggame,ggamt,ggamf)  ! smear 
         else
            ggam2e=pener
            ggam2t=theta
            ggam2f=phi
          if(cl_smear.eq.1) call smear(part_id,ggam2e,ggam2t,ggam2f)  ! smear 
         endif
         if(gngam.gt.1) then
          gnpi0=1
          gpi0e=ggam2e+ggame
          pi0x=ggame*sin(ggamt)*cos(ggamf)+ggam2e*sin(ggam2t)*cos(ggam2f)
          pi0y=ggame*sin(ggamt)*sin(ggamf)+ggam2e*sin(ggam2t)*sin(ggam2f)
          pi0z=ggame*cos(ggamt)+ggam2e*cos(ggam2t)
          gpi0f=Azimuthalphi(pi0x,pi0y)
          gpi0t=PolarTheta(pi0x,pi0y,pi0z)
         endif
        CASE DEFAULT

       END SELECT       
         return
         end

          subroutine smear(pid,e,t,f)
          IMPLICIT NONE
#include "dvcs.inc"
          integer pid
          real p,e,t,f,rn1,rn2,rn3,rn4,sin6
          real dts,ps,es,ts,fs,rn(4)
           es=e
           ts=t
           fs=f
           sin6=9.6                  ! sin6^o
           call MYRNORML(rn,4)       ! replacing normal random generator
c           call rannor(rn1,rn2)
c           call rannor(rn3,rn4)
            rn1=rn(1)
            rn2=rn(2)
            rn3=rn(3)
            rn4=rn(4)
          if(pid.eq.1) then ! photon separately
           e=es+cl_sme*sqrt(es)*rn1
           t=ts+cl_smf*rn2
           f=fs+cl_smg*rn3
           return
          endif
          SELECT CASE (pid)  
          CASE (3) ! electron
           p=e
          CASE (14) ! proton
           p=e*e-0.83827**2
           if(p.gt.0) then 
            p=sqrt(p)
           else
            p=0
           endif
          CASE (8) ! pi+
           p=e*e-0.139**2
           if(p.gt.0) then 
            p=sqrt(p)
           else
            p=0
           endif
          CASE (9) ! pi-
           p=e*e-0.139**2
           if(p.gt.0) then 
            p=sqrt(p)
           else
            p=0
           endif
          CASE (11) ! K+
           p=e*e-0.494**2
           if(p.gt.0) then 
            p=sqrt(p)
           else
            p=0
           endif

          CASE DEFAULT
          END SELECT       

           ps=p+(cl_sma+cl_smb*p)*p*rn1
           dts=(cl_smd+cl_smc/p)*rn2
           ts=t+dts
           fs=f+dts/sin6*rn3
           
          SELECT CASE (pid)  
          CASE (3) ! electron
           es=ps
          CASE (14) ! proton
           es=sqrt(ps*ps+0.83827**2)
           if(cl_proloss.eq.1)  then ! add loss
            if(es.lt.1.0) then
            es=es-0.002+(es-0.93827)/30.865       ! 0 at 1
            endif
           endif
          CASE (8) ! pi+
           es=sqrt(ps*ps+0.139**2)
          CASE (9) ! pi-
           es=sqrt(ps*ps+0.139**2)
          CASE (11) ! K+
           es=sqrt(ps*ps+0.494**2)
          CASE DEFAULT
          END SELECT       
           e=es
           f=fs
           t=ts
          return
          end



      INTEGER FUNCTION LUND2GEANTID(I)
      PARAMETER (NSEL=44)
      INTEGER IPLUND(NSEL),IGE
      DATA IPLUND/
     +    22,   -11,    11,    12,   -13,    13,   111,   211,
     +  -211,   130,   321,  -321,  2112,  2212, -2212,   310,
     +   221,  3122,  3222,  3212,  3112,  3322,  3312,  3334,
     + -2112, -3122, -3112, -3212, -3222, -3322, -3312, -3334,
     +   -15,    15,   411,  -411,   421,  -421,   431,  -431,
     +  4122,    24,   -24,    23/

       DO IGE=1,NSEL
        IF(I.EQ.IPLUND(IGE)) THEN
         LUND2GEANTID=IGE
         RETURN
        ENDIF
       ENDDO
         LUND2GEANTID=0 
       RETURN
       END
c
      INTEGER FUNCTION GEANT2LUNDID(I)
      PARAMETER (NSEL=44)
      INTEGER IPLUND(NSEL),I
      DATA IPLUND/
     +    22,   -11,    11,    12,   -13,    13,   111,   211,
     +  -211,   130,   321,  -321,  2112,  2212, -2212,   310,
     +   221,  3122,  3222,  3212,  3112,  3322,  3312,  3334,
     + -2112, -3122, -3112, -3212, -3222, -3322, -3312, -3334,
     +   -15,    15,   411,  -411,   421,  -421,   431,  -431,
     +  4122,    24,   -24,    23/

      IF(I.GT.0.AND.I.LE.NSEL) THEN
       GEANT2LUNDID=IPLUND(I)
      ELSE
       GEANT2LUNDID=0
      ENDIF 
       RETURN
      END

C********************************************************************* 
c
       subroutine vsumm(a,b,c,n)
       real a(n),b(n),c(n)
       integer i,n
        do i=1,n
          c(i)=a(i)+b(i)
        enddo
       return
       end
c
       subroutine vdifm(a,b,c,n)
       real a(n),b(n),c(n)
       integer i,n
        do i=1,n
          c(i)=a(i)-b(i)
        enddo
       return
       end
c
c
       real function vdotm(a,b,n)
       real a(n),b(n),s
       integer i,n
       s=0.0
       do i=1,3
         s=s+a(i)*b(i)
       enddo
       if(n.eq.4) s=s-a(n)*b(n)
       vdotm=s
       return
       end
c   
       real function vangle(a,b,c,d)
       real a(3),b(3),c(3),d(3),xm,ym,vcos
       real x(3),y(3),pi
       pi=acos(-1.0)
       call crossm(a,b,x)
       call crossm(c,d,y)
       xm=vdotm(x,x,3)
       ym=vdotm(y,y,3)
       if(xm.gt.0.0 .and. ym.gt.0.0) then
         vcos=vdotm(x,y,3)/sqrt(xm)/sqrt(ym)
         if(abs(vcos).lt.1.0) then
            vangle=acos(vcos)
         else
            if(vcos.ge.1.0)  vangle=0
            if(vcos.le.-1.0)  vangle=pi
         endif 
       else
         vangle=0
       endif
       return
       end
c
       subroutine crossm(a,b,c)
       real a(3),b(3),c(3)
       c(1)=a(2)*b(3)-a(3)*b(2)
       c(2)=a(3)*b(1)-a(1)*b(3)
       c(3)=a(1)*b(2)-a(2)*b(1)
       return
       end
c
c
       real function vmass(a)
       real vm,a(4)
        vm= vdotm(a,a,4)
        if (vm.lt.0.0) then
          vmass=sqrt(-vm)
        else
          vmass=-1.0
        endif 
       return
       end

c
      logical function  outoflimits(xb, Q2)
      implicit none
#include "dvcs.inc"
      double precision E, xb, Q2, del2,Phi_e,Phi_g,Phi_s,Phi_gb,dsigma
      double precision nu,W2,W,qmod,E1cm,P1cm,E2cm,P2cm,del2max,del2min
      double precision yb,ymax,ycol
      double precision  xmin1,xmax1
      double precision dsIunp,dsBHlp, dsIlp,dsBHtp,dsItp
c
      istatus=0                    ! kinematic range OK
      Ed=cl_be
      E=Ed
      xmin1 = Q2/(2D0*Mp*E)
      xmax1 = 1D0
      nu  = Q2/(2D0*Mp*xb)
      W2  = Mp**2 + 2D0*Mp*nu - Q2
c
      if( xb.lt.xmin1 .or. xb.gt.xmax1 ) istatus=1  !    x  out of range
      if( Q2.lt.cl_q2min .or. Q2.gt. cl_q2max) istatus=2  !    Q2  out of range
      yb=nu/Ed
      if(yb.gt.cl_ymax .or. yb.lt.cl_ymin) istatus=3        ! ymin<y<ymax
      if(W2.lt.cl_wmin) istatus=4
c$$$      call dvcsycol(del2,xb,Q2,ycol)
c$$$      if( (ycol-yb).le.cl_ycol ) istatus=4                   ! y-too big
      if(istatus.gt.0) then
       outoflimits=.TRUE.
      else
       outoflimits=.FALSE.
      endif
      return
      end


C        Ich  -  positron(electron) = +1(-1)
C int   hel  -  lepton beam polarization
C int   help -  target polarization
C double  dsigma   d^5\sig / dx dQ^2 d|\Delta^2| d\phi_e d\phi_\gamma
C                  (pb GeV^{-4})
c
      subroutine bmkxsec(xb, Q2, del2, Phi_e, Phi_g,dsigma)
      implicit none
#include "dvcsmom.inc"
#include "dvcs.inc"
#include "ntupgdvcs.inc"
      double precision E, xb, Q2, del2,Phi_e,Phi_g,Phi_s,Phi_gb,dsigma
      double precision nu,W2,W,qmod,E1cm,P1cm,E2cm,P2cm,del2max,del2min
      double precision yb,ymax,ycol
      double precision  xmin1,xmax1
      double precision dsIunp,dsBHlp, dsIlp,dsBHtp,dsItp
      real corraul,rP1,rP2
c
c
      istatus=0                    ! kinematic range OK
c      Ed=cl_be
      E=Ed
      xmin1 = Q2/(2D0*Mp*E)
      xmax1 = 1D0
      nu  = Q2/(2D0*Mp*xb)
      W2  = Mp**2 + 2D0*Mp*nu - Q2
      W   = sqrt(W2)
      qmod = sqrt(nu**2 + Q2)
c
      E1cm = Mp*(Mp + nu)/W
      P1cm = Mp*qmod/W
      E2cm = (W2 + Mp**2)/(2D0*W)
      P2cm = (W2 - Mp**2)/(2D0*W)
      del2max = 2D0*(Mp**2 - E1cm*E2cm - P1cm*P2cm)
      del2min = 2D0*(Mp**2 - E1cm*E2cm + P1cm*P2cm)
c
      if( xb.le.xmin1 .or. xb.gt.xmax1 ) istatus=1           !    x  out of range
      if( del2.ge.del2min .or. del2.le.del2max ) istatus=2   ! delta out of range
      yb=nu/Ed
      if(yb.gt.cl_ymax) istatus=3                               ! y<ymax
      call dvcsycol(del2,xb,Q2,ycol)

      Phi_gb=pi - Phi_g
      Phi_s=Phi_e
      call bhdvcs(xb,  Q2,  del2,  Phi_s,Phi_gb,rP1,rP2) 
c
      if(abs(rP1).le.cl_ycol ) istatus=4                   ! y-too big
c      if((ycol-yb).le.cl_ycol ) istatus=4                   ! y-too big

      if (istatus.eq.0) then
c
      dsBH =hc0BH +hc1BH*cos(Phi_gb)+hc2BH*cos(2D0*Phi_gb)
      dsDVCS=hc0dvcs+hc1dvcs*cos(Phi_gb)+hs1dvcs*sin(Phi_gb)+hs2dvcs*sin(2*Phi_gb)
      dsIunp=hc0Iunp +heli*hs1Iunp*sin(Phi_gb)+hc1Iunp*cos(Phi_gb) 
     6     +heli*hs2Iunp*sin(2*Phi_gb)+hc2Iunp*cos(2*Phi_gb)
c      
c
c      scale the sin\phi moment
       dsIunp=dsIunp*cl_scale
c
       if(cl_ktcor) then
         dsIunp=dsIunp*corraul(del2,Q2)     ! correct~20% of A_LU at small t
       endif
c
       dsigma=dsBH
       if(cl_bh.gt.1)  dsigma=dsigma+dsDVCS-Ich*dsIunp 

c
c      L-POL
c
       if(cl_pol.eq.1) then 
         dsBHlp =heli*helpi*(hc0BHlp +hc1BHlp*cos(Phi_gb))
         dsIlp=heli*helpi*hc0Ilp+heli*helpi*hc1Ilp*cos(Phi_gb)
     +         +helpi*hs1Ilp*sin(Phi_gb)
         dsigma=dsigma+dsBHlp -Ich*dsIlp
       endif
c
c      T-POL
c 
       if(cl_pol.eq.2) then 
         dsBHtp =heli*helpi*cos(Phi_s)*(hc0BHtpcos +hc1BHtpcos*cos(Phi_gb)) 
     6       +heli*hs1BHtpsin*sin(Phi_s)*sin(Phi_gb)
         dsItp=heli*hc0Itpcos*cos(Phi_s)+hc0Itpsin*sin(Phi_s)
     6     +hc1Itpsin*cos(Phi_gb)*sin(Phi_s)+heli*hc1Itpcos*cos(Phi_gb)*cos(Phi_s)
     6     +heli*hs1Itpsin*sin(Phi_gb)*sin(Phi_s)+hs1Itpcos*sin(Phi_gb)*cos(Phi_s)
         dsigma=dsigma+dsBHtp -Ich*dsItp
       endif
c
c
c
      else
        dsigma=0
c       print *,'out of limits ',xb,xmin1,xmax1,del2,del2min,del2max,istat      
      endif
       gwbh=dsBH
       gwdvcs=dsDVCS
       gvint=dsIunp
c       gsin=hs1Iunp/(hc0BH+hc0dvcs+hc0Iunp)
       gsin=hs1Iunp/hc0BH
       gsin2=dsigma
        gproh=helpi
        ghp1=hcp1
        ghp2=hcp2
        ghccb=hccb
        ghcci=hcci
        gh0bh=hc0BH
        gh1bh=hc1BH
        gh2bh=hc2BH
        gh0dvcs=hc0dvcs
        gh1dvcs=hc1dvcs
        ghs1dvcs=hs1dvcs
        gh0iunp=hc0Iunp
        gh1iunp=hc1Iunp
        ghs1iunp=hs1Iunp
        ghs2iunp=hs2Iunp
       if(cl_pol.eq.1) then
        gh0bhlp=hc0bhlp
        gh1bhlp=hc1bhlp
        ghs1ilp=hs1ilp
        gh1ilp=hc1ilp
        gh0ilp=hc0ilp
       else if(cl_pol.eq.2) then
        gh0bhtpcos = hc0BHtpcos
        gh1bhtpcos = hc1BHtpcos
        ghs1bhtpsin= hs1BHtpsin
        ghs1itpcos = hs1Itpcos
        ghs1itpsin = hs1Itpsin
        gh1itpcos  = hc1Itpcos
        gh1itpsin  = hc1Itpsin
        gh0itpcos  = hc0Itpcos
        gh0itpsin  = hc0Itpsin
      endif
c      print *,'bhdvcs rP1= ',rP1,rP2,dsigma,cl_ycol,ycol-yb,xb,Q2,del2
      return
      end
c
      subroutine excpi0(xb, Q2, del2, Phi_e, Phi_g,dsigma)
      implicit none
#include "dvcsmom.inc"
#include "dvcs.inc"
#include "ntupgdvcs.inc"
      double precision E, xb, Q2, del2,Phi_e,Phi_g,Phi_s,Phi_gb,dsigma
      double precision nu,W2,W,qmod,E1cm,P1cm,E2cm,P2cm,del2max,del2min
      double precision yb,ymax,ycol,dpi0u,dpi0l,dpi0t,sinpi0
      double precision  xmin1,xmax1,skew
      real HTre2,UFx,PFxq2,Tfac,u,aac2s0,dubbs,ddbbs,rskew
      real tt,tQ2,A1_PHI,A2_PHI,BT_Q2

c
c
      istatus=0                    ! kinematic range OK
      E=Ed
      xmin1 = Q2/(2D0*Mp*E)
      xmax1 = 1D0
      nu  = Q2/(2D0*Mp*xb)
      W2  = Mp**2 + 2D0*Mp*nu - Q2
      W   = sqrt(W2)
      qmod = sqrt(nu**2 + Q2)
      tt=-del2
      tQ2=Q2
c
      E1cm = Mp*(Mp + nu)/W
      P1cm = Mp*qmod/W
      E2cm = (W2 + Mp**2 - mpi0**2)/(2D0*W)     ! vpk noticed
      P2cm = (W2 - Mp**2)/(2D0*W)
      del2max = 2D0*(Mp**2 - E1cm*E2cm - P1cm*P2cm)
      del2min = 2D0*(Mp**2 - E1cm*E2cm + P1cm*P2cm)
c
      if( xb.le.xmin1 .or. xb.gt.xmax1 ) istatus=1           !    x  out of range
      if( del2.ge.del2min .or. del2.le.del2max ) istatus=2   ! delta out of range
      yb=nu/Ed
      if(W2.lt.cl_wmin) print *,'Debiles lriv!!!!!!!!!',W2
      if (istatus.eq.0) then
c
       sinpi0=2*xb*(1-xb)**5
c       PFxq2=(1-xb)**2/(Mp*Mp*xb/Q2+1-xb)**2
       PFxq2=84*xb*xb*xb*(1-xb) !(1-xb)**2/(Mp*Mp*xb/Q2+1-xb)**2
c       Tfac=exp(4*(del2-del2min))
cccc before 1/10/07       Tfac=exp(cl_bpi0*(del2-del2min))
        Tfac=exp(BT_Q2(tQ2)*(del2)) !-del2min))   ! from vpk
c
c
c      aac2s0(xb,Q2)-leading order pdf from http://spin.riken.bnl.gov/aac/ (hep-ph/0001046)
c
       skew= xb/(2.0 - xb)
c       UFx=aac2s0(skew,Q2) !HTre2(xb,del2min) !(du(xi)-dd(xi))**2
       rskew=xb !skew
       UFx=(2*dubbs(rskew)+ddbbs(rskew))**2
       dpi0u=0.3*Tfac*UFx*PFxq2
ccc before 1/10/07       dpi0u=dpi0u/Q2/Q2/(Q2+Mp**2)**2*(1.0+sinpi0*heli*sin(Phi_g)) 
       dpi0u=dpi0u/Q2/Q2/(Q2+Mp**2)**2*    ! vpk
     6(1.0+A1_PHI(tt)*cos(Phi_g)+A2_PHI(tt)*cos(2*Phi_g)
     6    +sinpi0*heli*sin(Phi_g))
       dsigma=dpi0u
c       print *,xb,u(xb)*xb,UFx,(u(xb)*xb)**2,UFx/(u(xb)*xb)**2
c       print *,'exc->',smax,q2,del2,'->',xb,PFxq2,UFx,dsigma
c
c      L-POL
c
c       print *,'------',cl_pol,dsigma
c
       if(cl_pol.eq.1) then 
        dpi0l=6*xb*xb*(1-xb)**5*helpi*sin(Phi_g) !????    ! longpol
        dsigma=dsigma*(1.0+dpi0l)
       endif
c
c      T-POL
c 
       if(cl_pol.eq.2) then 
         dpi0t=0 !???? ! tpol
         dsigma=dsigma*(1.0+dpi0t)
       endif
c
      else
        dsigma=0
c       print *,'out of limits ',xb,xmin1,xmax1,del2,del2min,del2max,istat      
      endif
c
       gwbh=dpi0u
       gwdvcs=dpi0l
       gvint=dpi0t
       gsin=sinpi0
      return
      end

c
       real function umrs(x)
       implicit none
       double precision x
c
c  x-dep MRS(G)
c
       umrs=2.26*x**(0.593-1.0)*(1.0-0.76*sqrt(x)+4.20*x)*(1.0-x)**3.96
c       u=2.26*x**(0.559-1.0)*(1.0-0.54*sqrt(x)+4.65*x)*(1.0-x)**3.96
       return
       end
c
c
       real function dmrs(x)
       implicit none
       double precision x
c
c  x-dep MRS(G)
c
       dmrs=0.279*x**(0.335-1.0)*(1.0+8.63*sqrt(x)+0.32*x)*(1.0-x)**4.41
       return
       end
c
      real function uplus(x)
c
c hep-ph/9401328:
c
c QCD Constrains on the Shape of Polarized Quark and Gluon Distributions
c    Author: Stanley J. Brodsky, Matthias Burkardt, Ivan Schmidt
c    Journal-ref: Nucl.Phys. B441 (1995) 197-214

      implicit none
      real qf
      real x,x1,x3,x4,alpha
      alpha=1.12
      x1=(1.0-x)
      x3=x1**3
      x4=x1*x3
      uplus=1.0/x**alpha*(3.784*x3-3.672*x4)
      return
      end
c
      real function dplus(x)
      implicit none
      real x,x1,x3,x4,alpha
      alpha=1.12
      x1=(1.0-x)
      x3=x1**3
      x4=x1*x3
      dplus=1.0/x**alpha*(0.757*x3-0.645*x4)
      return
      end
c
      real function umnus(x)
      implicit none
      real x,x1,x5,x6,alpha
      alpha=1.12
      x1=(1.0-x)
      x5=x1**5
      x6=x1*x5
      umnus=1.0/x**alpha*(2.004*x5-1.892*x6)
      return
      end
c
      real function dmnus(x)
      implicit none
      real x,x1,x5,x6,alpha
      alpha=1.12
      x1=(1.0-x)
      x5=x1**5
      x6=x1*x5
      dmnus=1.0/x**alpha*(3.23*x5-3.118*x6)
      return
      end

c
      real function ubbs(x)
      implicit none
      real t,x,uplus,umnus,dplus,dmnus
      ubbs=(uplus(x)+umnus(x))
      return
      end
c
c
      real function dubbs(x)
      implicit none
      real x,uplus,umnus,dplus,dmnus
         dubbs=uplus(x)-umnus(x)
      return
      end
c
c
      real function ddbbs(x)
      implicit none
      real x,uplus,umnus,dplus,dmnus
         ddbbs=dplus(x)-dmnus(x)
      return
      end
c
c
      subroutine bhdvcs(xbd, Q2d, del2d, phield, phigd,rP1,rP2)
      implicit double precision (A-H,O-Z)
#include "dvcsmom.inc"
#include "dvcs.inc"
c
C BH+DVCS+INT by V. Korotkov
C Tpol and Lpol added by  H.Avakian
c
C   Diff. cross-section for Bethe-Heitler, DVCS and interf.term 
C
C         d^5\sig / dx dQ^2 d|\Delta^2| d\phi_e d\phi_\gamma
C              (pb GeV^{-4})
C INPUT:
C        Ivar -  1 (BH only), 2 (DVCS + int only), 3 (BH + DVCS + int)
C        IGPD -  GPD variant
C        Ipn  -  proton(neutron) target = 1(2)
C real   E(GeV), xb, Q2(GeV^2), del2(GeV^2)(negative) 
C real   phiel(rad) - scatt. electron azimuthal angle
C real   phig(rad)  - photon azimuthal angle around momentum transfer vector q
C OUTPUT:
C real   hs*,hc* sin and cos moments of 5-fold diff. cross-section for
C                            BH, DVCS and interf. terms
*

      parameter ( alpha = 1D0/137.036D0, hc2  = 0.38938D0)
      parameter ( coeff = 1D+9*hc2*alpha**3 )
c
      real rP1,rP2
      double precision nu, k1pl, k2pl, Kfac, Jfac
      common/todvcs/ xb, yb, Q2, del2, del2min,phip,phipel,P1,P2,Kfac,Jfac,ds
*
      common/formfac/ F1pn(2), F2pn(2)
*
      xb    = xbd
      Q2    = Q2d
      del2  = del2d
*
      nu  = Q2d/(2D0*Mp*xbd)
      qmod = sqrt(nu**2 + Q2d)
      yb = nu/Ed
      Esc = Ed - nu

c
      eps = 2D0*xbd*Mp/sqrt(Q2)
      eps2=eps*eps
      qeps2=1D0 + eps2
      sqeps2=sqrt(qeps2)
      ds = coeff*(xbd*yb**2/(16D0*pi**2*Q2**2))/sqeps2
c
      del2min=-Q2*(2D0*(1D0-xb)*(1D0-sqeps2)+eps2)
      del2min=del2min/(4D0*xb*(1D0-xb)+eps2)
      tau  = del2d/(4D0*Mp**2)
      taum1=1D0-tau
      xtau=xbd*xbd/tau
      del2q2=del2d/Q2
      del2q4=del2q2*del2q2
      del2q2m1=1D0-del2q2
*
*
      phip = phigd 
      phipel = phield 
*
       call nuclFF( del2d )
C       call nuclFF_YAHL( del2d )

*
      y1eps=1D0 - yb - yb*yb*eps2/4D0
      sqy1eps=sqrt(y1eps)
      Kfac = sqrt((-del2q2)*(1D0 - xbd)*y1eps*
     *       (1D0 - del2min/del2d)*(sqrt(1D0 + eps2) + 
     *  ((4D0*xbd*(1D0 - xbd) + eps2)/(4D0*(1D0 - xbd)))*
     *                                       ((del2d - del2min)/Q2)))
      Jfac = (1D0 - yb - yb*eps2/2D0)*(1D0 + del2q2) - 
     *       (1D0 - xbd)*(2D0 - yb)*del2q2
      P1 = -(Jfac + 2D0*Kfac*cos(phip))/(yb*(1D0 + eps2))
      P2 = 1D0 + del2q2 - P1
      rP1=P1
      rP2=P2
c
        if( Ivar.eq.1 .or. Ivar.eq.3 ) then   
       F1 = F1pn(Ipn)
       F2 = F2pn(Ipn)    
       F1_M_F = F1 + tau*F2
       F1_M_F2 = F1**2 - tau*F2**2
       F1_P_F = F1 + F2
       F1_P_F2 = F1_P_F*F1_P_F
       c01_BH = 8D0*Kfac**2*((2D0 + 3D0*eps2)*Q2*F1_M_F2/del2d +
     *                                          2D0*xbd**2*F1_P_F2)
       c02_BH = (2D0 - yb)**2*((2D0 + eps2)*F1_M_F2*
     *       ((2D0*xbd*Mp)**2*(1D0 + del2q2)**2/del2d + 
     *                         4D0*(1D0 - xbd)*(1D0 + xbd*del2q2)) +
     *     4D0*xbd**2*F1_P_F2*(xbd + 
     *                   (1D0 - xbd + eps2/2D0)*(del2q2m1)**2 -
     *                    xbd*(1D0 - 2D0*xbd)*(del2q2)**2))
       c03_BH = 8D0*(1D0 + eps2)*(1D0 - yb - yb*yb*eps2/4D0)*
     *       (2D0*eps2*(1D0 - del2d/(4D0*Mp**2))*F1_M_F2 - 
     *           xbd**2*(del2q2m1)**2*F1_P_F2)
       c0_BH = c01_BH + c02_BH + c03_BH
       c1_BH = 8D0*Kfac*(2D0 - yb)*(
     *      F1_M_F2*(4D0*(xbd*Mp)**2/del2d - 2D0*xbd - eps2) +
     *      F1_P_F2*2D0*xbd**2*(1D0 - (1D0 - 2D0*xbd)*del2q2))
       c2_BH = 8D0*(xbd*Kfac)**2*(F1_M_F2*4D0*Mp**2/del2d + 2D0*F1_P_F2)
c
c     BH-lpol part
c
      bhkin1=8D0*xbd*yb*sqeps2/taum1*F1_P_F
      c01_BHlp=(xbd/2D0*del2q2m1-tau)/2D0
      c02_BHlp=2D0-xbd-2D0*(1D0-xbd)**2*del2q2+eps2*del2q2m1
     6         -xbd*(1D0-2D0*xbd)*del2q4
      c03_BHlp=1D0-(1D0-xbd)*del2q2
      c04_BHlp=xtau/4D0*(1D0+del2q2)**2+(1D0-xbd)*(1D0+xbd*del2q2)
      c0_BHlp=c01_BHlp*c02_BHlp*F1_P_F+c03_BHlp*c04_BHlp*F1_M_F
      c0_BHlp=c0_BHlp*(2D0-yb)*bhkin1
c
      c11_BHlp=(2D0*tau-xbd*del2q2m1)*(1D0-xbd+xbd*del2q2)  
      c12_BHlp=1D0+xbd-(3D0-2D0*xbd)*(1D0+xbd*del2q2)-xtau*(1D0+del2q4)
      c1_BHlp=c11_BHlp*F1_P_F+c12_BHlp*F1_M_F
      c1_BHlp=-c1_BHlp*Kfac*bhkin1
c
c
c     BH-Tpol part
c
      c01_BHtpcos=-8D0*(2D0-yb)*yb*sqrt(q2d)/Mp*sqeps2*Kfac/sqy1eps
c
      c02_BHtpcos=xbd*(xbd*Mp)**2/q2d*del2q2m1*F1_P_F
      c03_BHtpcos=1D0-(1D0-xbd)*del2q2
      c04_BHtpcos=xtau/4D0*del2q2m1*F1+xbd/2D0*F2
      c0_BHtpcos=c01_BHtpcos*F1_P_F*(c02_BHtpcos+c03_BHtpcos*c04_BHtpcos)
c
      c11_BHtpcos=-16D0*xbd*yb*sqy1eps*Mp/sqrt(q2d)
      c11_BHtpcos= c11_BHtpcos*sqeps2*F1_P_F
      c12_BHtpcos=2D0*Kfac*Kfac*q2d/del2/y1eps
      c13_BHtpcos=xbd*del2q2m1*F1+tau*F2
      c14_BHtpcos=qeps2*xbd*del2q2m1*F1_M_F
      c1_BHtpcos=c11_BHtpcos*(c12_BHtpcos*c13_BHtpcos+c14_BHtpcos)
c
      s11_BHtpsin=16D0*yb*xbd*xbd*sqy1eps/sqrt(q2d)*Mp
      s12_BHtpsin=qeps2*sqeps2*del2q2m1*F1_P_F*F1_M_F
      s1_BHtpsin=s11_BHtpsin*s12_BHtpsin
c
       BHfact=ds/((xbd*yb*(1D0 + eps2))**2*del2d*P1*P2)
       hcP1=P1
       hcP2=P2
       hccb=BHfact
       hc0BH=c0_BH*BHfact
       hc1BH=c1_BH*BHfact
       hc2BH=c2_BH*BHfact
c
       ac0bh=c0_BH
       ac1bh=c1_BH
       ac2bh=c2_BH
c
c
       hc0BHlp=c0_BHlp*BHfact
       hc1BHlp=c1_BHlp*BHfact
c
       hc0BHtpcos=c0_BHtpcos*BHfact
       hc1BHtpcos=c1_BHtpcos*BHfact
       hs1BHtpsin=s1_BHtpsin*BHfact
c
        endif
c
        if( Ivar.eq.2 .or. Ivar.eq.3 ) then       
         call dvcsfun()
        endif
c
      return
      end

      subroutine dvcsfun()
      implicit double precision (A-H,O-Z)
#include "dvcs.inc"
#include "dvcsmom.inc"
      double precision Kfac,Jfac,Intfac
      common/todvcs/ x, y,Q2,del2,del2min,phip,phipel,P1,P2,Kfac,Jfac,ds
      common/formfac/ F1pn(2), F2pn(2)
c
      skew = x/(2D0 - x)
      if(IGPD.lt.100) then
      call amptab(skew, del2, 
     &             H1_RE, H1_IM, H1T_RE, H1T_IM,
     &             E1_RE, E1_IM, E1T_RE, E1T_IM )
      else
      call amptabgag(skew, del2, 
     &             H1_RE, H1_IM, H1T_RE, H1T_IM,
     &             E1_RE, E1_IM, E1T_RE, E1T_IM )
      endif
*
C         proton/neutron
      F1 = F1pn(Ipn)
      F2 = F2pn(Ipn)
*
      deldel    = 1D0 - del2min/del2
      deldel_sq = sqrt(deldel)
      del2m2    = -del2/Mp**2
      del2m4    = -del2m2/4D0
      delm2_sq  = sqrt(del2m2)
      cy2   = 2D0 - 2D0*y + y**2
      Intfac=ds/(x*y**3*P1*P2*(-del2))
      hcci=Intfac
*
C  DVCS
C     
      a1 = H1_RE**2 + H1_IM**2 + H1T_RE**2 + H1T_IM**2
      a2 = 2.*( H1_RE*E1_RE  +  H1_IM*E1_IM +
     &         H1T_RE*E1T_RE + H1T_IM*E1T_IM )
      a3 =  E1_RE**2 +  E1_IM**2
      a4 = E1T_RE**2 + E1T_IM**2 
      C_DVCS = ( 4D0*(1D0-x)*a1 - a2*x**2 - (x**2 + (2D0-x)**2*del2m4)*a3 
     &           - x**2*del2m4*a4 )/(2D0 - x)**2
      C_DVCS_eff = -x*C_DVCS
      c0_DVCS = 2D0*cy2*C_DVCS
      c1_DVCS = 8D0*((2D0 - y)/(2D0 - x))*C_DVCS_eff
      T_DVCS  = (c0_DVCS + Kfac*c1_DVCS*cos(phip))/(y**2*Q2)
      DVCSfac=ds/(y**2*Q2)
      hc0dvcs=c0_DVCS*DVCSfac
      hc1dvcs=Kfac*c1_DVCS*DVCSfac
      ac0dvcs=c0_DVCS
      ac1dvcs=Kfac*c1_DVCS
C
C  INTERF
C
      C_I_re = F1*H1_RE + x/(2D0-x)*(F1+F2)*H1T_RE - del2m4*F2*E1_RE
      C_I_im = F1*H1_IM + x/(2D0-x)*(F1+F2)*H1T_IM - del2m4*F2*E1_IM
      RE2    = x/(2D0-x)*(H1_RE + E1_RE) + H1T_RE
      C_I_re_eff = -x*C_I_re
      C_I_im_eff = -x*C_I_im
      b1 = (2D0 - x)*(1D0 - y) 
      b1= b1 + (2D0-y)**2/(1D0-y)*Kfac*Kfac/del2*Q2 ! old - (1D0 - x)*(2D0 - y)**2*deldel
      b2 = (1D0 - y)*x*(F1 + F2)
      c0_I = -8D0*(2D0 - y)*( b1*C_I_re - b2*RE2 )
      c1_I = -8D0*cy2*C_I_re
      s1_I =  8D0*y*(2D0 - y)*C_I_im
      c2_I = -16D0*((2D0 - y)/(2D0 - x))*C_I_re_eff
      s2_I =  16D0*(y/(2D0 - x))*C_I_im_eff
c
c
C     moments
c
      hs2Iunp=Kfac*Kfac*s2_I*Intfac
      hs1Iunp=Kfac*s1_I*Intfac
      hc2Iunp=Kfac*Kfac*c2_I*Intfac
      hc1Iunp=Kfac*c1_I*Intfac
      hc0Iunp=del2/Q2*c0_I*Intfac
c
      ac0int=del2/Q2*c0_I
      ac1int=Kfac*c1_I
      as1int=Kfac*s1_I
c
c      print *,'mysl-dvcs',hs1Iunp,hc1Iunp,hc0Iunp
C
C  LPOL
C
      C_LP_re = (F1+F2)*skew*(H1_RE+x/2D0*E1_RE)+
     6          F1*H1T_RE-skew*(x/2D0*F1+del2m4*F2)*E1T_RE
c
c
      C_LP_im = (F1+F2)*skew*(H1_IM+x/2D0*E1_IM)+
     6          F1*H1T_IM-skew*(x/2D0*F1+del2m4*F2)*E1T_IM
c
      DC_LP_re=-skew*(F1+F2)*(H1_RE+x/2D0*E1_RE+skew
     6*(H1T_RE+x/2D0*E1T_RE))
      DC_LP_im=-skew*(F1+F2)*(H1_IM+x/2D0*E1_IM+skew
     6*(H1T_IM+x/2D0*E1T_IM))
c
      yf2=(2D0-y)**2/(1D0-y)+2
      c0lp_I =-8*y*(Kfac**2*yf2*C_LP_re+(1D0-y)*(2D0-x)*del2/Q2
     6        *(DC_LP_re+C_LP_re))
      s1lp_I = 8*Kfac*cy2*C_LP_im
      c1lp_I = -8*Kfac*y*(2D0-y)*C_LP_re

C
C  TPOL
C
      xb2=x*skew
      C_TPP_re_s = (F1+F2)*(xb2*(H1_RE+x/2D0*E1_RE)+x*del2m4*E1_RE)
     6          -xb2*F1*(H1T_RE+x/2D0*E1T_RE)
      C_TPP_re_b1 =del2m4*4D0*(1D0-x)/(2D0-x)*F2*H1T_RE
      C_TPP_re_b21 =-del2m4*x*F1*E1T_RE
      C_TPP_re_b22 =-del2m4*xb2*F2*E1T_RE
      C_TPP_re_b =C_TPP_re_b1+C_TPP_re_b21+C_TPP_re_b22
c
      C_TPP_im_s = (F1+F2)*(xb2*(H1_IM+x/2D0*E1_IM)+x*del2m4*E1_IM)
     6          -xb2*F1*(H1T_IM+x/2D0*E1T_IM)
      C_TPP_im_b =del2m4*(4D0*(1D0-x)/(2D0-x)*F2*H1T_IM
     6            -(x*F1+xb2*F2)*E1T_IM)
c
c
c
      C_TPM_re_s = 1D0/(2D0-x)*(x*x*F1-(1D0-x)*4D0*del2m4*F2)*H1_RE
     6          -xb2*(F1+F2)*(H1T_RE+del2m4*E1T_RE)
      C_TPM_re_b = (del2m4*((2D0-x)*F1+xb2*F2) +xb2*F1)*E1_RE
c
      C_TPM_im_s = 1D0/(2D0-x)*(x*x*F1-(1D0-x)*4D0*del2m4*F2)*H1_IM
     6          -xb2*(F1+F2)*(H1T_IM+del2m4*E1T_IM)
      C_TPM_im_b = (del2m4*((2D0-x)*F1+xb2*F2) +xb2*F1)*E1_IM
c
c
      C_TPM_re=C_TPM_re_s+C_TPM_re_b
      C_TPP_re=C_TPP_re_s+C_TPP_re_b
      C_TPM_im=C_TPM_im_s+C_TPM_im_b
      C_TPP_im=C_TPP_im_s+C_TPP_im_b
c
c
c
      DC_TPP_re=-4D0*del2m4*(F2*H1T_RE-x/(2D0-x)*(F1+x/2D0*F2)*E1T_RE)
      DC_TPM_re=4D0*del2m4*(F2*H1_RE-F1*E1_RE)
c
      DC_TPP_im=-4D0*del2m4*(F2*H1T_IM-x/(2D0-x)*(F1+x/2D0*F2)*E1T_IM)
      DC_TPM_im=4D0*del2m4*(F2*H1_IM-F1*E1_IM)
c
c
c
      qm8=8D0*Mp*sqrt(1D0-y)/sqrt(Q2)
      c0tpcos_I = -qm8*Kfac*y*(((2D0-y)**2/(1D0-y)+2D0)*C_TPP_re+DC_TPP_re)
      c0tpsin_I = qm8*Kfac*(2D0-y)*((2D0-y)**2/(1D0-y)*C_TPM_im+DC_TPM_im)
c
      c1tpcos_I = -qm8*y*(2D0-y)*C_TPP_re
      c1tpsin_I = qm8*cy2*C_TPM_im
c
      s1tpcos_I = qm8*cy2*C_TPP_im
      s1tpsin_I = -qm8*y*(2D0-y)*C_TPM_re

      hc0Itpcos = c0tpcos_I*Intfac
      hc0Itpsin = c0tpsin_I*Intfac
      hc1Itpcos = c1tpcos_I*Intfac
      hc1Itpsin = c1tpsin_I*Intfac
      hs1Itpcos = s1tpcos_I*Intfac
      hs1Itpsin = s1tpsin_I*Intfac
c
c 
      hs1Ilp = s1lp_I*Intfac
      hc1Ilp = c1lp_I*Intfac
      hc0Ilp = c0lp_I*Intfac
c
      return
      end

      real function HTre2(x,del2)
      implicit double precision (A-H,O-Z)
#include "dvcs.inc"
#include "dvcsmom.inc"
c
      skew = x/(2D0 - x)
      if(IGPD.lt.100) then
      call amptab(skew, del2, 
     &             H1_RE, H1_IM, H1T_RE, H1T_IM,
     &             E1_RE, E1_IM, E1T_RE, E1T_IM )
      else
      call amptabgag(skew, del2, 
     &             H1_RE, H1_IM, H1T_RE, H1T_IM,
     &             E1_RE, E1_IM, E1T_RE, E1T_IM )
      endif
C
      HTre2 = H1T_RE*H1T_RE+ H1T_IM*H1T_IM
      return
      end

      subroutine amptabgag(skew, del2, 
     &                   H1_RE, H1_IM, H1T_RE, H1T_IM,
     &                   E1_RE, E1_IM, E1T_RE, E1T_IM )
#include "dvcs.inc"
c
      REAL xi
      REAL mt
      REAL cff_re
      REAL cff_im
      REAL v_re,v_im
c
      double precision skew, del2, H1_RE, H1_IM, H1T_RE, H1T_IM,
     &                             E1_RE, E1_IM, E1T_RE, E1T_IM 
      double precision F1pn(2), F2pn(2)
      common/formfac/ F1pn, F2pn
      real mpi/0.1396/
      data init/1/
c
      SAVE
c
      if( init .eq. 1 ) then
        init = 0
        call GPD_Init(IGPD-100)
      endif
*
      F1u = 2D0*F1pn(1) + F1pn(2)
      F1d = 2D0*F1pn(2) + F1pn(1)
      F2u = 2D0*F2pn(1) + F2pn(2)
      F2d = 2D0*F2pn(2) + F2pn(1)
c
      xi=skew
      mt=-del2
      v_re = GPD_GET_CFF(xi,mt,0)
      v_im = GPD_GET_CFF(xi,mt,3)
      H1_RE=-v_re
      H1_IM=-v_im
      v_re = GPD_GET_CFF(xi,mt,1)
      v_im = GPD_GET_CFF(xi,mt,4)
      H1T_RE=-v_re
      H1T_IM=-v_im
      v_re = GPD_GET_CFF(xi,mt,2)
      v_im = GPD_GET_CFF(xi,mt,5)
      E1_RE=-v_re
      E1_IM=-v_im
      E1T_IM=0
      E1T_RE=0
c
      return
      end
c
c
c
      subroutine amptab(skew, del2, 
     &                   H1_RE, H1_IM, H1T_RE, H1T_IM,
     &                   E1_RE, E1_IM, E1T_RE, E1T_IM )
#include "dvcs.inc"
      double precision skew, del2, H1_RE, H1_IM, H1T_RE, H1T_IM,
     &                             E1_RE, E1_IM, E1T_RE, E1T_IM 
      double precision F1pn(2), F2pn(2)
      common/formfac/ F1pn, F2pn
*
      common/retbl/ vh1ure(51,5),  vh1dre(51,5),
     &              vh1ture(51,5), vh1tdre(51,5),
     &              ve1ure(51,5),  ve1dre(51,5),
     &              ve1ture(51,5), ve1tdre(51,5)
      common/imtbl/ vh1uim(51,5),  vh1dim(51,5),
     &              vh1tuim(51,5), vh1tdim(51,5),
     &              ve1uim(51,5),  ve1dim(51,5),
     &              ve1tuim(51,5), ve1tdim(51,5)
*
      common/retbl1/ vh1urenf(51,21,2), vh1drenf(51,21,2) 
      common/imtbl1/ vh1uimnf(51,21,2), vh1dimnf(51,21,2)
*
      real mpi/0.1396/
      data init/1/
      real  h1uimag,h1dimag,h1tuimag,h1tdimag
      common /h1imag/h1uimag,h1dimag,h1tuimag,h1tdimag
*
      SAVE
*
      if( init .eq. 1 ) then
        init = 0
        call rtable
        skewmin  = 0.01D0
        skewmax  = 1.00D0
        skewminl = log10(skewmin)
        skewmaxl = log10(skewmax)
        dskewl   = (skewmaxl - skewminl)/51D0
        dlmin    = 0.01D0
        dlmax    = 1.00D0
        dlminl   = log10(dlmin)
        dlmaxl   = log10(dlmax)
        ddll     = (dlmaxl - dlminl)/20D0
      endif
*
      F1u = 2D0*F1pn(1) + F1pn(2)
      F1d = 2D0*F1pn(2) + F1pn(1)
      F2u = 2D0*F2pn(1) + F2pn(2)
      F2d = 2D0*F2pn(2) + F2pn(1)
*
      DS   = (log10(skew) - skewminl)/dskewl + 1
      IS   = int(DS)
      FDS  = DS - real(IS)
      FDS1 = 1. - FDS
*
      DT   = (log10(-del2) - dlminl)/ddll + 1
      IT   = int(DT)
      FDT  = DT - real(IT)
      FDT1 = 1. - FDT
*
      if( IGPD .ge. 6 ) then
        IGPD1 = IGPD - 4
        IGPD2 = IGPD - 5
      else
        IGPD1 = IGPD
      endif
*
        if( IGPD .le. 5 ) then
      h1ure =   vh1ure(IS,IGPD1)*FDS1 + vh1ure(IS+1,IGPD1)*FDS
      h1dre =   vh1dre(IS,IGPD1)*FDS1 + vh1dre(IS+1,IGPD1)*FDS
*
      h1uim =   vh1uim(IS,IGPD1)*FDS1 + vh1uim(IS+1,IGPD1)*FDS
      h1dim =   vh1dim(IS,IGPD1)*FDS1 + vh1dim(IS+1,IGPD1)*FDS
        else
      h1ure = (vh1urenf(IS,IT,IGPD2)*FDS1 + 
     &           vh1urenf(IS+1,IT,IGPD2)*FDS)*FDT1 +
     &        (vh1urenf(IS,IT+1,IGPD2)*FDS1 +
     &           vh1urenf(IS+1,IT+1,IGPD2)*FDS)*FDT
      h1dre = (vh1drenf(IS,IT,IGPD2)*FDS1 + 
     &           vh1drenf(IS+1,IT,IGPD2)*FDS)*FDT1 +
     &        (vh1drenf(IS,IT+1,IGPD2)*FDS1 +
     &           vh1drenf(IS+1,IT+1,IGPD2)*FDS)*FDT
*
      h1uim = (vh1uimnf(IS,IT,IGPD2)*FDS1 + 
     &           vh1uimnf(IS+1,IT,IGPD2)*FDS)*FDT1 +
     &        (vh1uimnf(IS,IT+1,IGPD2)*FDS1 +
     &           vh1uimnf(IS+1,IT+1,IGPD2)*FDS)*FDT
      h1dim = (vh1dimnf(IS,IT,IGPD2)*FDS1 + 
     &           vh1dimnf(IS+1,IT,IGPD2)*FDS)*FDT1 +
     &        (vh1dimnf(IS,IT+1,IGPD2)*FDS1 +
     &           vh1dimnf(IS+1,IT+1,IGPD2)*FDS)*FDT
        endif
*
      h1ture = vh1ture(IS,IGPD1)*FDS1 + vh1ture(IS+1,IGPD1)*FDS
      h1tdre = vh1tdre(IS,IGPD1)*FDS1 + vh1tdre(IS+1,IGPD1)*FDS
*
      e1ure =   ve1ure(IS,IGPD1)*FDS1 + ve1ure(IS+1,IGPD1)*FDS
      e1dre =   ve1dre(IS,IGPD1)*FDS1 + ve1dre(IS+1,IGPD1)*FDS
*
      e1ture = ve1ture(IS,IGPD1)*FDS1 + ve1ture(IS+1,IGPD1)*FDS
      e1tdre = ve1tdre(IS,IGPD1)*FDS1 + ve1tdre(IS+1,IGPD1)*FDS
*
      h1tuim = vh1tuim(IS,IGPD1)*FDS1 + vh1tuim(IS+1,IGPD1)*FDS
      h1tdim = vh1tdim(IS,IGPD1)*FDS1 + vh1tdim(IS+1,IGPD1)*FDS
*
      e1uim =   ve1uim(IS,IGPD1)*FDS1 + ve1uim(IS+1,IGPD1)*FDS
      e1dim =   ve1dim(IS,IGPD1)*FDS1 + ve1dim(IS+1,IGPD1)*FDS
*
      e1tuim = ve1tuim(IS,IGPD1)*FDS1 + ve1tuim(IS+1,IGPD1)*FDS
      e1tdim = ve1tdim(IS,IGPD1)*FDS1 + ve1tdim(IS+1,IGPD1)*FDS
*
      gA  = 1.267D0/(1D0 - del2/0.84D0)**2
      gA0 = 0.6D0*gA
      gAu = 0.5D0*( gA + gA0)/( 0.8D0*1.267D0)
      gAd = 0.5D0*(-gA + gA0)/(-0.2D0*1.267D0)
      ha = 4D0*Mp**2*1.267D0/(mpi**2-del2)
*
        if( Ipn .eq. 1 ) then
          if( IGPD .le. 5 ) then
      h1uimag=(F1u/2D0)*h1uim
      h1dimag=(F1d/2D0)*h1dim
      H1_RE = (4D0*(F1u/2D0)*h1ure + F1d*h1dre)/9D0
      H1_IM = (4D0*(F1u/2D0)*h1uim + F1d*h1dim)/9D0
          else
      h1uimag=h1uim
      h1dimag=h1dim
      H1_RE = (4D0*h1ure + h1dre)/9D0
      H1_IM = (4D0*h1uim + h1dim)/9D0
          endif
      h1tuimag=gAu*h1tuim
      h1tdimag=gAd*h1tdim
      H1T_RE = (4D0*gAu*h1ture + gAd*h1tdre)/9D0
      H1T_IM = (4D0*gAu*h1tuim + gAd*h1tdim)/9D0
      E1_RE = (4D0*(F2u/2D0)*e1ure + F2d*e1dre)/9D0
      E1_IM = (4D0*(F2u/2D0)*e1uim + F2d*e1dim)/9D0
      E1T_RE = ha*(4D0*e1ture + e1tdre)/9D0
      E1T_IM = ha*(4D0*e1tuim + e1tdim)/9D0
        elseif( Ipn .eq. 2 ) then
      H1_RE = (4D0*(F1d/2D0)*h1dre + F1u*h1ure)/9D0
      H1_IM = (4D0*(F1d/2D0)*h1dim + F1u*h1uim)/9D0
      H1T_RE = (4D0*gAd*h1tdre + gAu*h1ture)/9D0
      H1T_IM = (4D0*gAd*h1tdim + gAu*h1tuim)/9D0
      E1_RE = (4D0*(F2d/2D0)*e1dre + F2u*e1ure)/9D0
      E1_IM = (4D0*(F2d/2D0)*e1dim + F2u*e1uim)/9D0
      E1T_RE = ha*(4D0*e1tdre + e1ture)/9D0
      E1T_IM = ha*(4D0*e1tdim + e1tuim)/9D0
        endif
*
      return
      end

      subroutine rtable
*
*      A - 1, B - 2, C - 3, D - 4, E - 5 
*
      common/retbl/ vh1ure(51,5),  vh1dre(51,5),
     &              vh1ture(51,5), vh1tdre(51,5),
     &              ve1ure(51,5),  ve1dre(51,5),
     &              ve1ture(51,5), ve1tdre(51,5)
      common/imtbl/ vh1uim(51,5),  vh1dim(51,5),
     &              vh1tuim(51,5), vh1tdim(51,5),
     &              ve1uim(51,5),  ve1dim(51,5),
     &              ve1tuim(51,5), ve1tdim(51,5)
*
*         6, 7
*
      common/retbl1/ vh1urenf(51,21,2), vh1drenf(51,21,2) 
      common/imtbl1/ vh1uimnf(51,21,2), vh1dimnf(51,21,2)
      character*500  pdffile,clasdvcspdf
* 
      call getenv('CLASDVCS_PDF',clasdvcspdf)
       print *,'READING the gpd.dat from directory $CLASDVCS_PDF='
     6,trim(clasdvcspdf)
       write(pdffile,'(A,A)') trim(clasdvcspdf),'/gpd.dat'
      open(unit=11,file=pdffile,status='old')
      do i = 1,51
        read(11,101) vh1ure(i,1),vh1dre(i,1),vh1ture(i,1),vh1tdre(i,1),
     &               ve1ture(i,1),ve1tdre(i,1) 
      enddo
      call ucopy(vh1ure(1,1),ve1ure(1,1),51)
      call ucopy(vh1dre(1,1),ve1dre(1,1),51)
      do j=2,5
        call ucopy(ve1ture(1,1),ve1ture(1,j),51)
        call ucopy(ve1tdre(1,1),ve1tdre(1,j),51)       
      enddo
*
      do j = 2,3
       do i = 1,51
        read(11,102) vh1ure(i,j),vh1dre(i,j),vh1ture(i,j),vh1tdre(i,j)
       enddo
      enddo
      call ucopy(vh1ure(1,2),ve1ure(1,2),102)
      call ucopy(vh1dre(1,2),ve1dre(1,2),102)
*
      do j = 4,5
       do i = 1,51
        read(11,101) vh1ure(i,j),vh1dre(i,j),vh1ture(i,j),vh1tdre(i,j),
     &               ve1ure(i,j),ve1dre(i,j)
       enddo
      enddo
*
      do i = 1,51
        read(11,101) vh1uim(i,1),vh1dim(i,1),vh1tuim(i,1),vh1tdim(i,1),
     &               ve1tuim(i,1),ve1tdim(i,1) 
      enddo
      call ucopy(vh1uim(1,1),ve1uim(1,1),51)
      call ucopy(vh1dim(1,1),ve1dim(1,1),51)
      do j=2,5
        call ucopy(ve1tuim(1,1),ve1tuim(1,j),51)
        call ucopy(ve1tdim(1,1),ve1tdim(1,j),51)       
      enddo
*
      do j = 2,3
       do i = 1,51
        read(11,102) vh1uim(i,j),vh1dim(i,j),vh1tuim(i,j),vh1tdim(i,j)
       enddo
      enddo
      call ucopy(vh1uim(1,2),ve1uim(1,2),102)
      call ucopy(vh1dim(1,2),ve1dim(1,2),102)
*
      do j = 4,5
       do i = 1,51
        read(11,101) vh1uim(i,j),vh1dim(i,j),vh1tuim(i,j),vh1tdim(i,j),
     &               ve1uim(i,j),ve1dim(i,j)
       enddo
      enddo
*
      do IG = 1,2
        do IT = 1,21
          do IS = 1,50
            read(11,103) vh1urenf(IS,IT,IG), vh1drenf(IS,IT,IG)
          enddo
        enddo
      enddo
*
      do IG = 1,2
        do IT = 1,21
          do IS = 1,50
            read(11,103) vh1uimnf(IS,IT,IG), vh1dimnf(IS,IT,IG)
          enddo
        enddo
      enddo
      close( 11 )
*
 101  format(6e12.4)
 102  format(4e12.4)
 103  format(2e12.4)
      return
      end

      subroutine nuclFF( del2 )
C
C  Elastic nucleon's formfactors
C
      implicit double precision (A-H,O-Z)
#include "dvcs.inc"
      common/formfac/ F1pn(2), F2pn(2)
      double precision Mv, kp, kn
      parameter (Mv = 0.843D0, kp = 1.79285D0, kn = -1.91D0)
C
      dipol = 1D0/(1D0 - del2/Mv**2)**2
*
      GE_p = dipol
      GE_n = 0D0
      GM_p = (1D0 + kp)*dipol
      GM_n =        kn*dipol
*
      delm = del2/(2D0*Mp)**2
*
      F1pn(1) = (GE_p - delm*GM_p)/(1D0-delm)
      F1pn(2) = (GE_n - delm*GM_n)/(1D0-delm)
      F2pn(1) = (GM_p - GE_p)/(1D0-delm)         
      F2pn(2) = (GM_n - GE_n)/(1D0-delm)         
      return
      end

      subroutine V3subd( A, B, C)
      implicit double precision (A-H,O-Z)
      dimension A(3), B(3), C(3)
      do i = 1,3
        C(i) = A(i) - B(i)
      enddo
      end

      double precision function V3dotd( A, B)
      implicit double precision (A-H,O-Z)
      dimension A(3), B(3)
      S = 0D0
      do i = 1,3
        S = S + A(i)*B(i)
      enddo
      V3dotd = S
      end

      subroutine nuclFF_YAHL( del2 )
C from the supplementary material of
C Zhihong Ye, John R. Arrington, Richard J.Hill and Gabriel Lee 
C "Proton and neutron electromagnetic form factors and uncertainties"
C Phys. Lett. B 777 (2018) 8--15
C doi: 10.1016/j.physletb.2017.11.023
C
C added by Sangbaek Lee
C
C
      implicit double precision (A-H,O-Z)
#include "dvcs.inc"
      common/formfac/ F1pn(2), F2pn(2)
      double precision Mv, kp, kn
      double precision z, a0(4), a1(4), a2(4), a3(4), a4(4), a5(4), a6(4), a7(4), a8(4), a9(4), a10(4), a11(4), a12(4)
      parameter (Mv = 0.843D0, kp = 1.79285D0, mu_n = -1.91D0, tcut = 4*0.13957**2, t0 = -0.7)
C
      dipol = 1D0/(1D0 - del2/Mv**2)**2
      z = (sqrt(tcut-del2)-sqrt(tcut-t0))/(sqrt(tcut-del2)+sqrt(tcut-t0)) 
*
      a0 = (/ 0.239163298067, 0.264142994136, 0.048919981379, 0.257758326959 /)
      a1 = (/ -1.109858574410, -1.095306122120, -0.064525053912, -1.079540642058 /)
      a2 = (/ 1.444380813060, 1.218553781780, -0.240825897382, 1.182183812195 /)
      a3 = (/ 0.479569465603, 0.661136493537, 0.392108744873, 0.711015085833 /)
      a4 = (/ -2.286894741870, -1.405678925030, 0.300445258602, -1.348080936796 /)
      a5 = (/ 1.126632984980, -1.356418438880, -0.661888687179, -1.662444025208 /)
      a6 = (/ 1.250619843540, 1.447029155340, -0.175639769687, 2.624354426029 /)
      a7 = (/ -3.631020471590, 4.235669735900, 0.624691724461, 1.751234494568 /)
      a8 = (/ 4.082217023790, -5.334045653410, -0.077684299367, -4.922300878888 /)
      a9 = (/ 0.504097346499, -2.916300520960, -0.236003975259, 3.197892727312 /)
      a10 = (/ -5.085120460510, 8.707403067570, 0.090401973470, -0.712072389946 /)
      a11 = (/ 3.967742543950, -5.706999943750, 0.0, 0.0 /)
      a12 = (/ -0.981529071103, 1.280814375890, 0.0, 0.0 /)

      GE_p = (a0(1)*z**0+a1(1)*z**1+a2(1)*z**2+a3(1)*z**3+a4(1)*z**4+a5(1)*z**5+a6(1)*z**6+a7(1)*z**7+a8(1)*z**8+a9(1)*z**9+a10(1)*z**10+a11(1)*z**11+a12(1)*z**12)
      GM_p = (1D0 + kp)*(a0(2)*z**0+a1(2)*z**1+a2(2)*z**2+a3(2)*z**3+a4(2)*z**4+a5(2)*z**5+a6(2)*z**6+a7(2)*z**7+a8(2)*z**8+a9(2)*z**9+a10(2)*z**10+a11(2)*z**11+a12(2)*z**12)
      GE_n = (a0(3)*z**0+a1(3)*z**1+a2(3)*z**2+a3(3)*z**3+a4(3)*z**4+a5(3)*z**5+a6(3)*z**6+a7(3)*z**7+a8(3)*z**8+a9(3)*z**9+a10(3)*z**10+a11(3)*z**11+a12(3)*z**12)
      GM_n = mu_n*(a0(4)*z**0+a1(4)*z**1+a2(4)*z**2+a3(4)*z**3+a4(4)*z**4+a5(4)*z**5+a6(4)*z**6+a7(4)*z**7+a8(4)*z**8+a9(4)*z**9+a10(4)*z**10+a11(4)*z**11+a12(4)*z**12)
*
      delm = del2/(2D0*Mp)**2
*
      F1pn(1) = (GE_p - delm*GM_p)/(1D0-delm)
      F1pn(2) = (GE_n - delm*GM_n)/(1D0-delm)
      F2pn(1) = (GM_p - GE_p)/(1D0-delm)         
      F2pn(2) = (GM_n - GE_n)/(1D0-delm)         
      return
      end

      subroutine resetmom()
c
C     set to 0 all moments
      implicit none
#include "dvcsmom.inc"
c
c  BH
c
       hc0BH=0
       hc1BH=0.0
       hc2BH=0.0
       hc0BHlp=0.0
       hc1BHlp=0.0
       hc0BHtpcos=0.0
       hc1BHtpcos=0.0
       hs1BHtpsin=0.0
c
c  DVCS
c      
       hc0dvcs=0.0
       hc1dvcs=0.0
       hs1dvcs=0.0
       hs2dvcs=0.0
c
c  INTERF
c       
       hs1Iunp=0.0
       hc1Iunp=0.0
       hc0Iunp=0.0
       hs2Iunp=0.0
       hc2Iunp=0.0
c       
       hs1Ilp=0.0
       hc1Ilp=0.0
       hc0Ilp=0.0
c
       hs1Itpcos=0.0
       hs1Itpsin=0.0
       hc1Itpcos=0.0
       hc1Itpsin=0.0
       hc0Itpcos=0.0
       hc0Itpsin=0.0
c            
         return
         end
c
      subroutine getpi0(xbd,Q2d,del2d,phield,phigd)
c
C     set to 0 all moments
      implicit none
#include "dvcs.inc"
#include "ludat1234.inc"
#include "ntupgdvcs.inc"
c
      dimension V3k1(3), V3q(3)
      dimension  V3p1(3), V3plus(3)
      double precision V3k1,V3q,V3p1,V3plus
      double precision xbd,Q2d,del2d,phield,phigd
      double precision yb,nu,Esc,sintel,costel,costVq,sintVq,qmod
      double precision Ep,Egam,costgg,sintgg,Vgx,Vgy,Vgz
      double precision costeg,teteg,coste1g,tete1g,V3dotd
      double precision cosphe,sinphe
      real pi0mom,PolarTheta,Azimuthalphi,pi0momx,pi0momy,pi0momz
      real pi0the,pi0phi
      integer i
ccc- Hyon-Suk
      common/radmick/PhRAD
      real PhRAD(4)
      common/radhs/radtheta
      real radtheta
ccc- Hyon-Suk
c
      nu  = Q2d/(2D0*Mp*xbd)
      qmod = sqrt(nu**2 + Q2d)
ccc- Hyon-Suk
c      if(Ed.ne.cl_beam_energy) print *,'dans getpi0 : Ed=',Ed
      if((PhRAD(4).gt.0.).and.(radtheta.le.0.15)) Ed=cl_beam_energy-PhRAD(4)
ccc- Hyon-Suk
      yb = nu/Ed
      Esc = Ed - nu
      costel = 1D0 - Q2d/(2D0*Ed*Esc)
      sintel = sqrt(1D0 - costel**2)
c
      Ep   = Mp - del2d/(2D0*Mp)
      Egam = nu + del2d/(2D0*Mp)
c
      V3k1(1) = 0D0
      V3k1(2) = 0D0
      V3k1(3) = Ed
      cosphe=cos(phield)
      sinphe=sin(phield)
      V3k2(1) = Esc*sintel*cosphe
      V3k2(2) = Esc*sintel*sinphe   !0D0
      V3k2(3) = Esc*costel
*
      call V3subd( V3k1, V3k2, V3q)
      costVq = V3q(3)/qmod
      sintVq = sqrt(1D0 - costVq**2)
*
      if(cl_pi0) then
       pi0mom=sqrt(Egam*Egam-mpi0*mpi0)
       costgg = (2D0*Egam*(Mp + nu) + Q2d - 2D0*Mp*nu-mpi0*mpi0)
      else if(cl_eta) then
       pi0mom=Egam*Egam-metta*metta
       if(pi0mom.gt.0.0) then
        pi0mom=sqrt(pi0mom)
        costgg = (2D0*Egam*(Mp + nu) + Q2d - 2D0*Mp*nu-metta*metta)
       else
        print *,'problem',pi0mom,Egam,costgg
        pi0mom=0.0
        costgg=2.0
       endif
      else
        print *,'error: unknown particle'
       stop
      endif
      costgg = costgg/(2D0*pi0mom*qmod)
      if(abs(costgg).gt.1) then
c       print *,'costgg-------->',costgg   ! shit happens
       costgg=1.0
       sintgg=0.0
      else
        sintgg = sqrt(1D0 - costgg**2)
      endif
c
      Vgx = pi0mom*sintgg*cos(phigd)
      Vgy = pi0mom*sintgg*sin(phigd)
      Vgz = pi0mom*costgg
*
      V3gam(1) = Vgx*costVq*cosphe - Vgz*sintVq*cosphe - Vgy*sinphe ! Vgx*costVq - Vgz*sintVq
      V3gam(2) = Vgx*costVq*sinphe - Vgz*sintVq*sinphe + Vgy*cosphe ! Vgy
      V3gam(3) = Vgx*sintVq        + Vgz*costVq                     ! Vgx*sintVq + Vgz*costVq
*
      call V3subd( V3q, V3gam, V3p2)
ccc- Hyon-Suk
      if(cl_radgen) then
       if (PhRAD(4).gt.0.) then
c        print *,'getpi0 : radtheta=',radtheta
c        print *,'!!!!dans getpi0: PhRAD = ',PhRAD
        if(radtheta.gt.0.15) then
         PhRAD(1) = PhRAD(4)*sintel*cosphe
         PhRAD(2) = PhRAD(4)*sintel*sinphe
         PhRAD(3) = PhRAD(4)*costel
         V3k2(1)=V3k2(1)-PhRAD(1)
         V3k2(2)=V3k2(2)-PhRAD(2)
         V3k2(3)=V3k2(3)-PhRAD(3)
        else
         PhRAD(1) = 0.
         PhRAD(2) = 0.
         PhRAD(3) = PhRAD(4)
        endif
       endif
      endif
ccc- Hyon-Suk
c
c     decay pi-0
c
       pi0momx=v3gam(1)
       pi0momy=v3gam(2)
       pi0momz=v3gam(3)
       gpi0e=Egam
       pi0phi=Azimuthalphi(pi0momx,pi0momy)
       pi0the=PolarTheta(pi0momx,pi0momy,pi0momz)
       gpi0f=pi0phi
       gpi0t=pi0the
 10    continue
       if(cl_pi0) then
        CALL LU1ENT(1,111,gpi0e,pi0the,pi0phi)
       else
        CALL LU1ENT(1,221,gpi0e,pi0the,pi0phi)
       endif
c
       CALL LUEXEC   ! decay pi0 
       if(cl_verblev.eq.2) write (6,*),gpi0e,N,(k(i,2),i=1,N)
c
       if(k(3,2).ne.22) then
        dalOK=.TRUE.          ! decay pi0 ->gamma,e+/e- 
        if (k(3,2).eq.11) then ! first electron
         v3el(1)=P(3,1)
         v3el(2)=P(3,2)
         v3el(3)=P(3,3)
         v3po(1)=P(4,1)
         v3po(2)=P(4,2)
         v3po(3)=P(4,3)
        else if (k(3,2).eq.-11) then ! positron
         v3po(1)=P(3,1)
         v3po(2)=P(3,2)
         v3po(3)=P(3,3)
         v3el(1)=P(4,1)
         v3el(2)=P(4,2)
         v3el(3)=P(4,3)
        else
          if(cl_verblev.gt.1) print *,'Unknown PID!!!! ',k(3,2)
          goto 10
        endif
c
       else
        dalOK=.FALSE.          ! decay pi0 -> 2gamma  
        v3gam2(1)=P(3,1)
        v3gam2(2)=P(3,2)
        v3gam2(3)=P(3,3)
       endif
c       call lulist(1) 
       v3gam(1)=P(2,1)
       v3gam(2)=P(2,2)
       v3gam(3)=P(2,3)
c
*
      return
      end
c
c
c
      subroutine getgammastar(xbd,Q2d,phield,VPGEN)
c
C     set to 0 all moments
      implicit none
#include "dvcs.inc"
c
      real vepr(4),ve0(4),VPGEN(4)
c
      double precision xbd,Q2d,del2d,phield,qmod
      double precision yb,nu,Esc,sintel,costel
      double precision cosphe,sinphe
c
      nu  = Q2d/(2D0*Mp*xbd)
      qmod = sqrt(nu**2 + Q2d)
      yb = nu/Ed
      Esc = Ed - nu
      costel = 1D0 - Q2d/(2D0*Ed*Esc)
      sintel = sqrt(1D0 - costel**2)
c
      cosphe=cos(phield)
      sinphe=sin(phield)
      vepr(1) = Esc*sintel*cosphe
      vepr(2) = Esc*sintel*sinphe   !0D0
      vepr(3) = Esc*costel
      vepr(4) =sqrt(vepr(1)*vepr(1)+vepr(2)*vepr(2)+vepr(3)*vepr(3))
      ve0(1)=0
      ve0(2)=0
      ve0(3)=cl_beam_energy
      ve0(4)=cl_beam_energy
c
         call vdifm(ve0,vepr,VPGEN,4)
c
      return
      end
c
      subroutine getphoton(xbd,Q2d,del2d,phield,phigd)
c
C     set to 0 all moments
      implicit none
#include "dvcs.inc"
#include "ludat1234.inc"
#include "ntupgdvcs.inc"
*
      dimension V3k1(3), V3q(3)
      dimension  V3p1(3), V3plus(3)
      double precision V3k1,V3q,V3p1,V3plus
      double precision xbd,Q2d,del2d,phield,phigd
      double precision yb,nu,Esc,sintel,costel,costVq,sintVq,qmod
      double precision Ep,Egam,Eele_inner,costgg,sintgg,Vgx,Vgy,Vgz
      double precision costeg,teteg,coste1g,tete1g,V3dotd
      double precision cosphe,sinphe
      integer j
ccc- Hyon-Suk
      common/radmick/PhRAD
      real PhRAD(4)
      common/radhs/radtheta, ichannel
      real radtheta
ccc- Hyon-Suk
c
      nu  = Q2d/(2D0*Mp*xbd)
      qmod = sqrt(nu**2 + Q2d)
ccc- Hyon-Suk
c      if(Ed.ne.cl_beam_energy) print *,'dans getphoton : Ed=',Ed
c      if((PhRAD(4).gt.0.).and.(ichannel.eq.2)) Ed=cl_beam_energy-PhRAD(4)
ccc- Hyon-Suk
      yb = nu/Ed
      Esc = Ed - nu
      costel = 1D0 - Q2d/(2D0*Ed*Esc)
      sintel = sqrt(1D0 - costel**2)
      cosphe=cos(phield)
      sinphe=sin(phield)
c
      if(cl_radgen) then
         if (PhRAD(4).gt.0.) then
            if (ichannel.gt.2) then
               PhRAD(1) = PhRAD(4)*sintel*cosphe
               PhRAD(2) = PhRAD(4)*sintel*sinphe
               PhRAD(3) = PhRAD(4)*costel
            else
               PhRAD(1) = 0.
               PhRAD(2) = 0.
               PhRAD(3) = PhRAD(4)
               Ed=cl_beam_energy-PhRAD(4)
            endif
         endif
      endif

c
      V3k1(1) = 0D0
      V3k1(2) = 0D0
      V3k1(3) = Ed

      V3k2(1) = Esc*sintel*cosphe
      V3k2(2) = Esc*sintel*sinphe   !0D0
      V3k2(3) = Esc*costel

      if(cl_radgen) then
         if (PhRAD(4).gt.0.) then
            if (ichannel.gt.2) then
               Eele_inner = Esc + PhRad(4)
               V3plus(1) = Eele_inner*sintel*cosphe
               V3plus(2) = Eele_inner*sintel*sinphe
               V3plus(3) = Eele_inner*costel
               call V3subd( V3k1, V3plus, V3q)
               Q2d = -(Ed-Esc-PhRad(4))**2 + V3q(1)**2 + V3q(2)**2 + V3q(3)**2
               nu = Ed-Eele_inner
               qmod =  sqrt(V3q(1)**2 + V3q(2)**2 + V3q(3)**2)
            else
               call V3subd( V3k1, V3k2, V3q)
               Q2d = -(Ed-Esc)**2 + V3q(1)**2 + V3q(2)**2 + V3q(3)**2
               nu = Ed-Esc
               qmod =  sqrt(V3q(1)**2 + V3q(2)**2 + V3q(3)**2)
            endif
         else
            call V3subd( V3k1, V3k2, V3q)
         endif
      else
            call V3subd( V3k1, V3k2, V3q)
      endif

      radxb_true = Q2d/2/nu/Mp
      radq2_true = Q2d

      Ep   = Mp - del2d/(2D0*Mp)
      Egam = nu + del2d/(2D0*Mp)

      costVq = V3q(3)/qmod
      sintVq = sqrt(1D0 - costVq**2)
*
      costgg = (2D0*Egam*(Mp + nu) + Q2d - 2D0*Mp*nu)/(2D0*Egam*qmod)
      sintgg = sqrt(1D0 - costgg**2)
      Vgx = Egam*sintgg*cos(phigd)
      Vgy = Egam*sintgg*sin(phigd)
      Vgz = Egam*costgg
*
      V3gam(1) = Vgx*costVq*cosphe - Vgz*sintVq*cosphe - Vgy*sinphe ! Vgx*costVq - Vgz*sintVq
      V3gam(2) = Vgx*costVq*sinphe - Vgz*sintVq*sinphe + Vgy*cosphe ! Vgy
      V3gam(3) = Vgx*sintVq        + Vgz*costVq                     ! Vgx*sintVq + Vgz*costVq
*
      call V3subd( V3q, V3gam, V3p2)
ccc- Hyon-Suk
ccc- Hyon-Suk
      N=2
c
c
        p(2,4)=0.0      ! gamma_E
        p(2,5)=0.0      ! gamma_M
       DO j=1,3
        p(2,j)=V3gam(j)
        p(2,4)=p(2,4)+V3gam(j)*V3gam(j)
       enddo
        if( p(2,4).gt.0) p(2,4)=sqrt(p(2,4))
c
        k(2,2)=22   ! photon
        k(2,1)=1    ! final
      
*
      return
      end

c
      real function aac2s0(xx,QQ2)
      double precision  xx,QQ2
      aac2s0=(duaac(xx,QQ2)-ddaac(xx,QQ2))**2
        return
	end
c
c
      real function duaac(xx,QQ2)
      double precision xx,QQ2
      real*8 x,Q2
      REAL*8 POLPDF(-3:3),STRUCT(3)
      integer iset
      iset=1
      Q2=QQ2
      x=xx
      CALL AACPDF(ISET,Q2,X,POLPDF,STRUCT)      
      duaac=POLPDF(1)/xx
        return
	end
c
      real function ddaac(xx,QQ2)
      double precision  xx,QQ2
      real*8 x,Q2
      REAL*8 POLPDF(-3:3),STRUCT(3)
      integer iset
      iset=1
      Q2=QQ2
      x=xx
      CALL AACPDF(ISET,Q2,X,POLPDF,STRUCT)      
      ddaac=POLPDF(2)/xx
        return
	end
c
c
      integer function length(str)
c     return the string length without the blanks characters

      implicit integer (k-l)
      character *(*) str    
      
      lmax=len(str)      
       
c     search the last non blank character
      do i=lmax,1,-1
      if(str(i:i).ne.' ')then
      length=i
      return
      end if
      end do

      length=lmax

      return
      end     

      real function valnum(str)
c     return the real value contained into a string

      implicit integer (k-l)
      character *(*) str
      logical segno
      

      segno=.false.
      valnum=0.00
      lu=length(str)
      
c     check the number sign
      if(str(1:1).eq.'-')then
      segno=.true.
      str=str(2:lu)
      lu=lu-1
      end if

c     check if number is float or integer
      if(index(str,'.').ne.0)then
      iin=index(str,'.')-1
      else
      iin=lu
      end if

      ifr=lu-(iin+1)

c     translate the integer portion
      doi=1,iin  
      k=ichar(str(i:i))-48
      valnum=valnum+float(k)*10.00**float(iin-i)
      end do         

      if(iin.eq.lu)goto 10
      str=str(iin+2:lu)

c     translate the decimal portion
      doi=1,ifr  
      k=ichar(str(i:i))-48
      valnum=valnum+float(k)/10.00**float(i)
      end do

10    if(segno)valnum=-valnum

      return
      end     
c
      integer function ivalnum(str)
c     return the real value contained into a string

      implicit integer (k-l)
      character *(*) str
      logical segno
      double precision valnum
      

      segno=.false.
      valnum=0.00
      lu=length(str)
      
c     check the number sign
      if(str(1:1).eq.'-')then
      segno=.true.
      str=str(2:lu)
      lu=lu-1
      end if

c     check if number is float or integer
      if(index(str,'.').ne.0)then
      iin=index(str,'.')-1
      else
      iin=lu
      end if

      ifr=lu-(iin+1)

c     translate the integer portion
      doi=1,iin  
      k=ichar(str(i:i))-48
      valnum=valnum+float(k)*10.00**float(iin-i)
      end do         

      if(iin.eq.lu)goto 10
      str=str(iin+2:lu)

c     translate the decimal portion
      doi=1,ifr  
      k=ichar(str(i:i))-48
      valnum=valnum+float(k)/10.00**float(i)
      end do

10    if(segno)valnum=-valnum
      ivalnum=int(valnum)
      return
      end     
c
c     k_T correction to A_LU from Marc Vdh program
c

       real function corraul(delt,dQ2)
       double precision delt,dQ2
       real t,Q2,qfac,tfac
       t=-delt
       Q2=dQ2
       if(Q2.lt.1.5) Q2=1.5
       if(Q2.gt.2.6) Q2=2.6
       qfac=(0.69*t-1.04*t*t)*(q2-1.5)/1.1+1.0
       tfac=(0.908-alog(t)/9.5)
       corraul=1.0/qfac*tfac
       return 
       end
c
       subroutine RNDMSEED()
       implicit none
#include "dvcs.inc"
       integer*4 now(3),ifilenum
       integer iran,irannom
       real rand,rnd
       integer :: n
       integer,allocatable :: seeds(:)

c
c
      if(cl_seed.eq.0) then
      call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
  	  irannom=now(3)+now(1)*3600+now(2)*60
  	  else 
  	  irannom=cl_seed
  	  endif
      call random_seed(size=n) ! get size of the seed used by the OS
      allocate(seeds(n))
      seeds = irannom 
      call random_seed(put=seeds)
      deallocate(seeds)
      CALL SRAND(irannom)
      call random_number(rnd)  ! random number in [0,1]
      print *,'Rndm numb->',irannom,now(1),now(2),now(3),rnd
c           DO 100 iran=1,irannom
c 100       RND=rand()
c
	   return
	   end

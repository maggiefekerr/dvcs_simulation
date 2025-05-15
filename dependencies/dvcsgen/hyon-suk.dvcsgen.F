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
      call CL_OPTIONS          ! command line processing
      call INIT_PAW            ! init hbook
c
C..   number of events   
c
       NEVENT=cl_triggers


C.. Constats and Parameter
C.. Initialization of the lepto program 
c
       CALL  CLASBOSINIT('MCEVENT')
c
      CALL DVCSINIT
      if(cl_verblev.gt.0) WRITE(*,*) 'NEVENT=',NEVENT

C.. Event generation: 
       DO I = 1,NEVENT
          print *,'Event nb: ',I
       
C.. Output of the number of generated events 
         IF (MOD(I,cl_nprint).EQ.0.and.cl_verblev.gt.0)
     6   WRITE(*,*) I,'  from ',NEVENT
C..
C..
c         print *,'juste avant cleannt'
         CALL CLEANNT
         CALL genEVENT(I)       
C..
C..      fill in BOS and ntuples
C.. 
        CALL CLASBOSFILL(I)
C..
       ENDDO

         if (datfileOK) close(41)
         if(ntOK) call CLOSE_PAW()

       CALL CLASBOSEND('MCEVENT')
       WRITE(*,*) 'Total ',NEVENT,' in pb ',smax
C.. Saving the last random number configuration 
       END      
c
c
c
       SUBROUTINE INIT_PAW
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
       hbooknm='exclpi0.hbook'
      else
       hbooknm='dvcsgen.hbook'
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
     +'gproh:I,gh0bh:R,gh1bh:R,gh2bh:R,'//
     +'gh0dvcs:R,gh1dvcs:R,ghs1dvcs:R,'//
     +'gh0iunp:R,gh1iunp:R,ghs1iunp:R,ghs2iunp:R,'//
     +'gh0bhlp:R,gh1bhlp:R,ghs1ilp:R,gh1ilp:R')
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

      if (cl_radgen) then    
       call HBNAME(NTID_OUT,'radgen',iproctype,
     +'iproctype:I,radq2:R,radnu:R,radweight:R,'//
     +'grade:R,gradt:R,gradf:R')
      endif
      end
c
c
       SUBROUTINE CL_OPTIONS
       implicit none
#include "dvcs.inc"
       integer i,iran,irannom
       real VEC1,valnum
       character*20 cnumber
       integer  numopts,NCDECI
       integer*4 now(3)
c
c    command line processing
c
c     defaults
c
      cl_dvcs=.TRUE.         ! default dvcs on
      cl_pi0 =.FALSE.        ! default no pi0
      cl_eta =.FALSE.        ! default no etta
      cl_gpd=4               ! gpd model
      cl_ycol=0.05           ! default ycol-y >0.05 def
      cl_zpos=0              ! z-position
      cl_zwidth=1.25         ! z-width ( zpos+/-zwidth)
      cl_pol=0               ! unpolarized target
      cl_beam_energy=5.75    ! def e1_dvcs
      cl_bh=0.0              ! def BH const C dN/Dk=1/(k+C)
      cl_wmin=1.85           ! def w2min
      cl_ymin=0.2            ! def ymin -> Emax=(1-ymin)*E
      cl_ymax=0.94           ! def ymax
      cl_thmin=0.0            ! def e'tmin
      cl_tmax=6.0             ! def tmax
      cl_tmin=0.0             ! def tmin
      cl_q2max=6.0            ! def Q2 max
      cl_q2min=0.8            ! def Q2 min
      cl_xmax=0.95             ! def x max
      cl_xmin=0.             ! def x min
      cl_thmax=1.57            ! def e'tmax
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
      cl_verblev=1
      cl_writef=1            ! gsim lund format
      cl_nprint=1000         ! print every cl_nprint event
      cl_triggers = 10  
      bosout='dvcsgen.evt'
      cl_mom=.FALSE.
      cl_radgen=.FALSE.
      bosOK=.TRUE.
      bossmearOK=.FALSE.
      datfileOK=.TRUE.      ! write a data file
      weightOK=.FALSE.
      acce16=.FALSE.
      acceg1=.FALSE.
      acc12=.FALSE.
      ntOK=.TRUE.
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
          cl_triggers=NCDECI(cnumber)
        elseif(cnumber.eq.'--v'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_verblev=NCDECI(cnumber)
        elseif(cnumber.eq.'--writef'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_writef=NCDECI(cnumber)
        elseif(cnumber.eq.'--print'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_nprint=NCDECI(cnumber)
        elseif(cnumber.eq.'--lpol') then
           cl_pol=1
        elseif(cnumber.eq.'--tpol') then
           cl_pol=2
        elseif(cnumber.eq.'--radgen') then
           cl_radgen=.TRUE.
        elseif(cnumber.eq.'--mom') then
           cl_mom=.TRUE.
        elseif(cnumber.eq.'--pi0') then
           cl_pi0=.TRUE.
        elseif(cnumber.eq.'--eta') then
           cl_eta=.TRUE.
        elseif(cnumber.eq.'--nodvcs') then
           cl_dvcs=.FALSE.
        elseif(cnumber.eq.'--beam'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_beam_energy=valnum(cnumber)           
        elseif(cnumber.eq.'--bh'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_bh=valnum(cnumber)           
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
           cl_thmin=NCDECI(cnumber)*0.001
           i=i+1
           CALL GETARG(i,cnumber)
           cl_thmax=NCDECI(cnumber)*0.001
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
        elseif(cnumber.eq.'--y'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_ymin=valnum(cnumber)
           i=i+1
           CALL GETARG(i,cnumber)
           cl_ymax=valnum(cnumber)
        elseif(cnumber.eq.'--zpos'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_zpos=NCDECI(cnumber)*0.1           
        elseif(cnumber.eq.'--ycol'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_ycol=valnum(cnumber)           
        elseif(cnumber.eq.'--gpd'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_gpd=NCDECI(cnumber)           
        elseif(cnumber.eq.'--zwidth'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           cl_zwidth=NCDECI(cnumber)*0.1           
        elseif(cnumber.eq.'--targ'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cl_target)
        elseif(cnumber.eq.'--bos') then
           bosOK=.TRUE.
        elseif(cnumber.eq.'--bossmear') then
           bossmearOK=.TRUE.
        elseif(cnumber.eq.'--weight') then
           weightOK=.TRUE.
        elseif(cnumber.eq.'--nodat') then
           datfileOK=.FALSE.
        elseif(cnumber.eq.'--acce16') then
           acce16=.TRUE.
        elseif(cnumber.eq.'--acc12') then
           acc12=.TRUE.
        elseif(cnumber.eq.'--acceg1') then
           acceg1=.TRUE.
        elseif(cnumber.eq.'--weight') then
           weightOK=.TRUE.
        elseif(cnumber.eq.'--nont') then
           ntOK=.FALSE.
        elseif(cnumber.eq.'--rndm'.and.i.lt.numopts) then
           i=i+1
           CALL GETARG(i,cnumber)
           irannom=valnum(cnumber)
           call RLUXGO(4,irannom,0,0)
c        call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
c           irannom=10000*now(3)+100*now(2)+now(1)
c        print *,'Random number for this set->',irannom,now(1),now(2),now(3)
c           DO 100 iran=1,irannom
c 100       CALL RANMAR(VEC1,1)
        elseif(cnumber.eq.'--file'.and.i.lt.numopts) then
          i=i+1
          CALL GETARG(i,bosout)
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
        return 
        end

        subroutine printusage
        print *,'dvcspi0gen [options]            '
        print *,'      option  value      default    comment  '
        print *,'     --pi0                    exclusive pi-0 on '
        print *,'     --eta                    exclusive eta on '
        print *,'     --nodvcs                 DVCS off '
        print *,'     --v    verbos_level    0       additional printout '
        print *,'     --trig nevents         10      number of triggers '
        print *,'     --beam beam_energy     5.776    e- momentum in GeV '
        print *,'     --ycol ycol_cut        0.05      ycol-y>ycol_cut '
        print *,'     --x  x-min   x-max   0.15 0.65   min max x_Bj '
        print *,'     --q2 Q2-min Q2-max   1.2  8.0    min max Q2 in GeV^2 '
        print *,'     --y y-min y-max      0.05 0.9    min max for y=nu/E      '
        print *,'     --w2 w^2-min         4.0         min for w^2=M^2+2MyE-Q^2      '
        print *,'     --t tmin tmax         0 1.0      t  min/max in GeV^2 '
        print *,'     --th thmin thmax            theta min/max for \theta in mrad'
        print *,'     --zpos z-position      0    target z position in mm '
        print *,'     --zwidth z-width       0    width in z in mm (zpos+/-zwidth) '
        print *,'     --file bos_file     dvcsgen.evt  bos outfile '
        print *,'     --weight                    flat distributions with weight(part12)'
        print *,'     --bos                write out the bos file '
        print *,'     --nont               do not write out the ntuple '
        print *,'     --rndm               randomize '
        print *,'     --gpd  Igpd          1    GPD model (1-A,2-B,3-C,4-D) '
        print *,'     --targ target       proton   deuteron/neutron possible '
        print *,'     --lpol                    Long.pol.target  '
        print *,'     --tpol                    Trans.pol.target  '
        print *,'     --mom                     include moments in ntuple  '
        print *,'     --proloss                  add proton loss'
        print *,'     --radgen                   include radgen'
        print *,'     --nodat                    do not write a data 
     +	file'
        print *,'     --acce16                   include e16 acceptance 
     +	for e-'
        print *,'     --acceg1                   include eg1 acceptance 
     +	for e-'
        print *,'     --acc12                   include clas++ acceptance 
     +	for e-'
        print *,'     --smear                   smear moments'
        print *,'     --bossmear                smear also bos output'
        print *,'     --A    value         0.006   smear param-A'
        print *,'     --B    value         0.001  smear param-B'
        print *,'     --C    value         0.0008  smear param-C'
        print *,'     --D    value         0.001  smear param-D'
        print *,'     --print nprint       1000   print every nprint event '
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
        print *,'BH C->dN/dk=1/(k+C)  :  ',cl_bh  
        print *,'w^2   min    :         ',cl_wmin  
        print *,'y   min/max:         ',cl_ymin,cl_ymax  
        print *,'Eprime thmin/thmax:         ',cl_thmin,cl_thmax  
        print *,'-Delta  tmin/tmax:         ',cl_tmin,cl_tmax  
        print *,'Q2  q2min/q2max:         ',cl_q2min,cl_q2max  
        print *,'x  xmin/xmax:         ',cl_xmin,cl_xmax  
        print *,'TARGET:              ', cl_target 
        print *,'TARGET z-position:   ', cl_zpos 
        print *,'TARGET 1/2 width :   ', cl_zwidth 
        print *,'Tpol (0-U,1-L,2-T):   ', cl_pol 
        print *,'GPD model A,B,C,D:   ', cl_gpd 
        print *,'RADGEN ON        :   ', cl_radgen 
        print *,'y_col cut        :   ', cl_ycol 
        print *,'acceptance       :   ',acceg1,acce16,acc12
        print *,'write form      :   ',cl_writef
        if(cl_mom) print *,'include moments in the ntuple:   ', cl_mom 
        if(cl_smear.eq.1) 
     6print *,'smearing A,B,C,D:   ', cl_sma,cl_smb,cl_smc,cl_smd
     6,cl_proloss,'in bos ',bossmearOK  
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
        integer nx,nq,nt2,nfe,nfg
        real  xmax,xmin,q2max,q2min,tmax,tmin 
        common /xmax/  xmax,xmin,q2max,q2min,tmax,tmin 
     6                 ,nx,nq,nt2,nfe,nfg
C..
c
          CHARACTER*2 CTARGET
          integer lst40
          real  beam_energy
c
C..
        real  elmax,elmin,ttmax,ttmin 
C..
        Ifl=1   !dsigma/dx/dq/dt/dphiel/dpigam
        Ivar=3  ! BH+DVCS
        IGPD=cl_gpd  ! GPD type
        cl_be=cl_beam_energy
        Ich=-1
C.. 
        nx=40
        nq=20
        nt2=40
        nfe=36
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
       elmin=(1.0-cl_ymin)*cl_beam_energy
C..  upper limit scat. lep. energy
       elmax=(1.0-cl_ymax)*cl_beam_energy
c
C.. lower limit scat lep. angle
       ttmin=cl_thmin !0.52-0.08          ! 29.75 degree+/-4
C..  upper limit scat lep. angle
       ttmax=cl_thmax !0.52+0.08
c
C.. Momentum of incoming lepton (in GeV)
       Ed = cl_beam_energy  

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
        print *,'CALL getSCALE-->',Ipn,xmax,xmin,tmax,tmin
        CALL getSCALE

        if(cl_radgen) then
         lst40=0  ! unpolarized case
         beam_energy=cl_beam_energy
         print *, 'call radgen_init(ctarget,beam_energy,lst40)'
        call radgen_init(ctarget,beam_energy,lst40)
c	stop
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
c
        smax=0.0
        la0=0
        dx=(xmax-xmin)/nx
        dq=(q2max-q2min)/nq
        dt=(tmax-tmin)/nt
        dfe=2*pi/nfe
        dfg=2*pi/nfg
        if(cl_pol.ge.1) la0=-1
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

       if(cl_dvcs) then
        CALL bmkxsec(xbd, Q2d, del2d, phield,phigd,dstot)
       else if(cl_pi0) then
        CALL excpi0(xbd, Q2d, del2d, phield,phigd,dstot)
       else
        print *,'ERROR! at least one reaction should be ON'
        stop
       endif

             stotint=stotint+dstot*cix*cit*ciq*cife*cifg
c	if (dstot.ne.0.) print*,stotint,dstot,cix,cit,ciq,
c     +	cife, cifg, xbd, Q2d, del2d, phield
c
           if(smax.lt.dstot) then
            smax=dstot
            if(cl_verblev.gt.1) print *,
     6'smax',heli,helpi,xbd,q2d,del2d,stotint,dstot,smax
c,etime()
           endif
c
            ENDDO
           ENDDO
        ENDDO
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
        common/radhs/radtheta
        real radtheta
ccc- Hyon-Suk

C..
C..
        integer nx,nq,nt2,nfe,nfg
        real  xmax,xmin,q2max,q2min,tmax,tmin 
        common /xmax/  xmax,xmin,q2max,q2min,tmax,tmin 
     6                 ,nx,nq,nt2,nfe,nfg
C..     generate a DVCS event
      double precision xbd,Q2d,del2d,phield,phigd,dstot,yb,w2
c
c       radiative stuff
c
        real VPGEN(4),VPRAD(4),PhRAD(4),q2tr,utr,weight
        REAL PolarTheta, Azimuthalphi

        REAL VEC(8),e_photon,prob_max,prob_curr,ggam2p
        INTEGER I,ITRY,LEN
        logical outoflimits

        LEN=8
        ITRY=0
 10     CALL RANLUX(VEC,LEN)
        ITRY=ITRY+1
        xbd=xmin+(xmax-xmin)*VEC(1)
        Q2d=q2min+(q2max-q2min)*VEC(2)
        
        if(cl_radgen) then
ccc- Hyon-Suk
         PhRAD(1)=0.
         PhRAD(2)=0.
         PhRAD(3)=0.
         PhRAD(4)=0.
ccc- Hyon-Suk         
         radq2=Q2d
         radnu=Q2d/2.0/Mp/xbd
         if( outoflimits(xbd,Q2d) ) then
          if(itry .gt.10000.and.cl_verblev.eq.1000) 
     6print *,'stuck before radgen:',itry,istatus,xbd,q2d
          goto 10 
         endif
         call getgammastar(xbd,Q2d,phield,VPGEN)
c     print *,'call RADGEN(cl_beam_energy,VPGEN,VPRAD,PhRAD,q2tr,
c     utr,weight)'
	 call RADGEN(cl_beam_energy,VPGEN,VPRAD,PhRAD,q2tr,utr,weight)
         Q2d=q2tr
         xbd=Q2d/2.0/Mp/utr
         if( outoflimits(xbd,Q2d) ) then
          if(itry .gt.10000.and.cl_verblev.eq.1000) 
     6print *,'stuck after radgen:',itry,istatus,xbd,q2d
          goto 10 
         endif
c
         ggam2p=PhRAD(1)**2+PhRAD(2)**2+PhRAD(3)**2
         PhRAD(4)=sqrt(ggam2p)
         if(ggam2p.gt.0) then
           gradt=PolarTheta(PhRAD(1),PhRAD(2),PhRAD(3))
           gradf=Azimuthalphi(PhRAD(1),PhRAD(2))
           grade=PhRAD(4)
          else 
           gradt=0
           gradf=0
           grade=0
         endif   
           radweight=weight !del2d
           gelebits=ita
ccc- Hyon-Suk
         radtheta=gradt
ccc- Hyon-Suk
        endif
        
c        if((PhRAD(4).gt.0.1).and.(radtheta.gt.0.15)) print *,'genevent : PhRAD = ',PhRAD
        
        del2d=tmin+(tmax-tmin)*VEC(3)
        phield=2*pi*VEC(4)
        phigd=2*pi*VEC(5)
c        yb=Q2d/(2D0*Mp*xbd*cl_beam_energy)
c        w2=Mp*Mp+2.0*Mp*cl_beam_energy*yb-Q2d
ccc- Hyon-Suk
        if(radtheta.gt.0.15) then
         yb=Q2d/(2D0*Mp*xbd*cl_beam_energy)
         w2=Mp*Mp+2.0*Mp*cl_beam_energy*yb-Q2d
        else
         yb=Q2d/(2D0*Mp*xbd*(cl_beam_energy-PhRAD(4)))
         w2=Mp*Mp+2.0*Mp*(cl_beam_energy-PhRAD(4))*yb-Q2d
        endif
ccc- Hyon-Suk
        if(yb.gt.cl_ymax.or.w2.lt.cl_wmin) goto 10
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
        CALL bmkxsec(xbd, Q2d, del2d, phield,phigd,dstot)
        else
        CALL excpi0(xbd, Q2d, del2d, phield,phigd,dstot)
        endif
c
c       if(itry.GE.100000)  print *,'Check the limits',itry,smax
c     & ,heli,xbd, Q2d, del2d, phield, phigd, dstot
        if(weightOK) then
          if(dstot.le.1.0E-8) goto 10      !put in the weight 
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
       gelex=del2d
       geley=Q2d
       gelez=phigd
       gprox=xbd
       gproy=phield
       gproz=0
       return
       end

C**********************************************************************
C      random number generator
C**********************************************************************
         real function random_num()
         random_num=rndm(-1)
         return
         end


C***********************************************************************
      SUBROUTINE CLEANNT
C***********************************************************************
#include "ntupgdvcs.inc"
#include "dvcsmom.inc"
       call VZERO(gnele,12)
       call VZERO(gnpi0,35)
       call resetmom()

       return
       end
c
C***********************************************************************
         SUBROUTINE smearbos(part,part_id,pener,ptheta,pphi)
C***********************************************************************
          IMPLICIT NONE
#include "bcs.inc"
         real pmom,pener,ptheta,pphi
         integer part,part_id
         SELECT CASE (part_id)  
         CASE (3) ! electron
         pmom=pener
         CASE (14) ! proton
         pmom=pener*pener-0.88035
         if(pmom.gt.0) pmom=sqrt(pmom)
         CASE (1) ! gamma
         pmom=pener
         CASE DEFAULT
         END SELECT       
c
         rw(part+5)	= pener  ! energy (pmom=p^2) 
         rw(part+6)	= pmom*sin(ptheta)*cos(pphi)		! px
         rw(part+7)	= pmom*sin(ptheta)*sin(pphi)		! py
         rw(part+8)	= pmom*cos(ptheta)		! pz
         print *,part,part_id,pener,ptheta,pphi
         return
         end
c
c
c
C***********************************************************************
      SUBROUTINE CLASBOSFILL(IEVENT)
C***********************************************************************
#include "bcs.inc"
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
          integer ierr,lunout
	  parameter (lunout=33)
          INTEGER ind,indx,mctk,mcvx,mbank,part,boswrite
          INTEGER j,jj
          LOGICAL lOK
          real pcharge(10),z_cible
C ********************************************************
C       Write to the output event file
C ********************************************************
       call bdrop(iw,'E')
       call bgarb(iw)
       ind=mbank(iw,'HEAD',0,8,1)
       if(ind.ne.0)then
          iw(ind+1)=1
          iw(ind+2)=1
          iw(ind+3)=IEVENT
          iw(ind+4)= 100                	! event time
          iw(ind+5)=-4
          iw(ind+6)=0
          iw(ind+7)=15
          iw(ind+8)=0
       endif
      ind= nbank('TGBI',0,4,1)
       if(ind.ne.0)then
         if (cl_pol.gt.0.5) then
          iw(ind+1) = 32768
         else
          iw(ind+1) = 0 
         endif
      endif
c
c     fill lund
c
c    electron
       DO j=1,3
        p(1,j)=V3k2(j)
        p(2,j)=V3p2(j)
        p(3,j)=V3gam(j)
        if(cl_pi0) p(4,j)=V3gam2(j)     ! second photon    
       ENDDO
c
c 
        p(1,5)=mele ! mass of el
        p(2,5)=Mp    !mass of proton
        p(3,5)=0.0      ! gamma
        pcharge(1)=-1
        pcharge(2)=1
        pcharge(3)=0
       k(1,2)=11  ! electron
       k(2,2)=2212 ! proton
       k(3,2)=22   ! photon
       jj=3
       if(cl_pi0) then      ! second photon from pi-0
        k(4,2)=22
        p(4,5)=0.0
        pcharge(4)=0 
        jj=4                 ! number of particles
       endif

ccc- Hyon-Suk
c       if(phrad(4).gt.0.) then
c       print *,'=== clasbosfill : PhRAD = ',PhRAD
c       print *,'==== before PhRAD ===='
c       print *,'old pxe=',p(1,1)
c       print *,'old pye=',p(1,2)
c       print *,'old pze=',p(1,3)
c       print *,'old Ee=',sqrt(p(1,5)*p(1,5)+p(1,1)*p(1,1)+p(1,2)*p(1,2)+p(1,3)*p(1,3))
c       print *,'==== after PhRAD ===='
c       print *,'new pxe=',p(1,1)
c       print *,'new pye=',p(1,2)
c       print *,'new pze=',p(1,3)
c       print *,'pxrad=',PhRAD(1)
c       print *,'pyrad=',PhRAD(2)
c       print *,'pzrad=',PhRAD(3)
c       print *,'new Ee=',sqrt(p(1,5)*p(1,5)+p(1,1)*p(1,1)+p(1,2)*p(1,2)+p(1,3)*p(1,3))
c       print *,'Erad=',sqrt(PhRAD(1)*PhRAD(1)+PhRAD(2)*PhRAD(2)+PhRAD(3)*PhRAD(3))
c       print *,'PhRAD(4)=',PhRAD(4)
c       print *,'==== checks ===='
c       print *,'check px=',p(1,1)+PhRAD(1)
c       print *,'check py=',p(1,2)+PhRAD(2)
c       print *,'check pz=',p(1,3)+PhRAD(3)
c       print *,'check E=',sqrt(p(1,5)*p(1,5)+p(1,1)*p(1,1)+p(1,2)*p(1,2)+p(1,3)*p(1,3))+sqrt(PhRAD(1)*PhRAD(1)+PhRAD(2)*PhRAD(2)+PhRAD(3)*PhRAD(3))
c       print *,'==== without PhRAD ===='
c       print *,'!!!!!! px : -p(1,1)-p(2,1)-p(3,1) = ',-p(1,1)-p(2,1)-p(3,1)
c       print *,'pxe=',p(1,1)
c       print *,'pxp=',p(2,1)
c       print *,'pxg=',p(3,1)
c       print *,'pxrad=',PhRAD(1)
c       print *,'!!!!!! py : -p(1,2)-p(2,2)-p(3,2) = ',-p(1,2)-p(2,2)-p(3,2)
c       print *,'pye=',p(1,2)
c       print *,'pyp=',p(2,2)
c       print *,'pyg=',p(3,2)
c       print *,'pyrad=',PhRAD(2)
c       print *,'!!!!!! pz : 5.776-p(1,3)-p(2,3)-p(3,3) = ',5.776-p(1,3)-p(2,3)-p(3,3)
c       print *,'pze=',p(1,3)
c       print *,'pzp=',p(2,3)
c       print *,'pzg=',p(3,3)
c       print *,'pzrad=',PhRAD(3)
c       print *,'!!!!!! E : cl_beam_energy+Mp-sqrt(p(1,5)*p(1,5)+p(1,1)*p(1,1)+p(1,2)*p(1,2)+p(1,3)*p(1,3))-sqrt(p(2,5)*p(2,5)+p(2,1)*p(2,1)+p(2,2)*p(2,2)+p(2,3)*p(2,3))-sqrt(p(3,5)*p(3,5)+p(3,1)*p(3,1)+p(3,2)*p(3,2)+p(3,3)*p(3,3)) = ',cl_beam_energy+Mp-sqrt(p(1,5)*p(1,5)+p(1,1)*p(1,1)+p(1,2)*p(1,2)+p(1,3)*p(1,3))-sqrt(p(2,5)*p(2,5)+p(2,1)*p(2,1)+p(2,2)*p(2,2)+p(2,3)*p(2,3))-sqrt(p(3,5)*p(3,5)+p(3,1)*p(3,1)+p(3,2)*p(3,2)+p(3,3)*p(3,3))
c       print *,'==== with PhRAD ===='
c       print *,'!!!!!! px : -p(1,1)-p(2,1)-p(3,1)-p(4,1)-PhRAD(1) = ',-p(1,1)-p(2,1)-p(3,1)-p(4,1)-PhRAD(1)
c       print *,'pxe=',p(1,1)
c       print *,'pxp=',p(2,1)
c       print *,'pxg1=',p(3,1)
c       print *,'pxg2=',p(4,1)
c       print *,'pxrad=',PhRAD(1)
c       print *,'!!!!!! py : -p(1,2)-p(2,2)-p(3,2)-p(4,2)-PhRAD(2) = ',-p(1,2)-p(2,2)-p(3,2)-p(4,2)-PhRAD(2)
c       print *,'pye=',p(1,2)
c       print *,'pyp=',p(2,2)
c       print *,'pyg1=',p(3,2)
c       print *,'pyg2=',p(4,2)
c       print *,'pyrad=',PhRAD(2)
c       print *,'!!!!!! pz : 5.776-p(1,3)-p(2,3)-p(3,3)-p(4,3)-PhRAD(3) = ',5.776-p(1,3)-p(2,3)-p(3,3)-p(4,3)-PhRAD(3)
c       print *,'pze=',p(1,3)
c       print *,'pzp=',p(2,3)
c       print *,'pzg1=',p(3,3)
c       print *,'pzg2=',p(4,3)
c       print *,'pzrad=',PhRAD(3)
c       print *,'!!!!!! E : cl_beam_energy+Mp-sqrt(p(1,5)*p(1,5)+p(1,1)*p(1,1)+p(1,2)*p(1,2)+p(1,3)*p(1,3))-sqrt(p(2,5)*p(2,5)+p(2,1)*p(2,1)+p(2,2)*p(2,2)+p(2,3)*p(2,3))-sqrt(p(3,5)*p(3,5)+p(3,1)*p(3,1)+p(3,2)*p(3,2)+p(3,3)*p(3,3))-sqrt(p(4,5)*p(4,5)+p(4,1)*p(4,1)+p(4,2)*p(4,2)+p(4,3)*p(4,3))-PhRAD(4) = ',cl_beam_energy+Mp-sqrt(p(1,5)*p(1,5)+p(1,1)*p(1,1)+p(1,2)*p(1,2)+p(1,3)*p(1,3))-sqrt(p(2,5)*p(2,5)+p(2,1)*p(2,1)+p(2,2)*p(2,2)+p(2,3)*p(2,3))-sqrt(p(3,5)*p(3,5)+p(3,1)*p(3,1)+p(3,2)*p(3,2)+p(3,3)*p(3,3))-sqrt(p(4,5)*p(4,5)+p(4,1)*p(4,1)+p(4,2)*p(4,2)+p(4,3)*p(4,3))-PhRAD(4)
c       print *,'cl_beam_energy=',cl_beam_energy
c       print *,'Mp=',Mp
c       print *,'p(1,5)=',p(1,5)
c       print *,'p(2,5)=',p(2,5)
c       print *,'p(3,5)=',p(3,5)
c       print *,'p(4,5)=',p(4,5)
c       print *,'Ee=',sqrt(p(1,5)*p(1,5)+p(1,1)*p(1,1)+p(1,2)*p(1,2)+p(1,3)*p(1,3))
c       print *,'Ep=',sqrt(p(2,5)*p(2,5)+p(2,1)*p(2,1)+p(2,2)*p(2,2)+p(2,3)*p(2,3))
c       print *,'Eg1=',sqrt(p(3,5)*p(3,5)+p(3,1)*p(3,1)+p(3,2)*p(3,2)+p(3,3)*p(3,3))
c       print *,'Eg2=',sqrt(p(4,5)*p(4,5)+p(4,1)*p(4,1)+p(4,2)*p(4,2)+p(4,3)*p(4,3))
c       print *,'Erad=',sqrt(PhRAD(1)*PhRAD(1)+PhRAD(2)*PhRAD(2)+PhRAD(3)*PhRAD(3))
c       print *,'PhRAD(4)=',PhRAD(4)
c       print *,'======== what is actually written in BOS ========'
c       print *,'pxe=',p(1,1)/(sqrt(p(1,1)*p(1,1)+p(1,2)*p(1,2)+p(1,3)*p(1,3)))
c       print *,'pxp=',p(2,1)/(sqrt(p(2,1)*p(2,1)+p(2,2)*p(2,2)+p(2,3)*p(2,3)))
c       print *,'pxg1=',p(3,1)/(sqrt(p(3,1)*p(3,1)+p(3,2)*p(3,2)+p(3,3)*p(3,3)))
c       print *,'pxg2=',p(4,1)/(sqrt(p(4,1)*p(4,1)+p(4,2)*p(4,2)+p(4,3)*p(4,3)))
c       print *,'pye=',p(1,2)/(sqrt(p(1,1)*p(1,1)+p(1,2)*p(1,2)+p(1,3)*p(1,3)))
c       print *,'pyp=',p(2,2)/(sqrt(p(2,1)*p(2,1)+p(2,2)*p(2,2)+p(2,3)*p(2,3)))
c       print *,'pyg1=',p(3,2)/(sqrt(p(3,1)*p(3,1)+p(3,2)*p(3,2)+p(3,3)*p(3,3)))
c       print *,'pyg2=',p(4,2)/(sqrt(p(4,1)*p(4,1)+p(4,2)*p(4,2)+p(4,3)*p(4,3)))
c       print *,'pze=',p(1,3)/(sqrt(p(1,1)*p(1,1)+p(1,2)*p(1,2)+p(1,3)*p(1,3)))
c       print *,'pzp=',p(2,3)/(sqrt(p(2,1)*p(2,1)+p(2,2)*p(2,2)+p(2,3)*p(2,3)))
c       print *,'pzg1=',p(3,3)/(sqrt(p(3,1)*p(3,1)+p(3,2)*p(3,2)+p(3,3)*p(3,3)))
c       print *,'pzg2=',p(4,3)/(sqrt(p(4,1)*p(4,1)+p(4,2)*p(4,2)+p(4,3)*p(4,3)))
c       endif
ccc- Hyon-Suk

c       part= nbank('PART',0,13,jj) 

c       print *,heli

       part = 1
       if(part.ne.0)then
         do 999 j=1,jj
          pmom=p(j,1)*p(j,1)+p(j,2)*p(j,2)+p(j,3)*p(j,3)


          z_cible = cl_zpos+(2.*random_num()-1)*cl_zwidth	! z vertex position
c       if(part.ne.0) then
c         iw(part+1)	= lund2geantid(k(j,2))	! particle ID (GEANT)
c         rw(part+2)	= 0		! x vertex position
c         rw(part+3)	= 0		! y vertex position
c         rw(part+4)	= z_cible
c         rw(part+5)	= sqrt(p(j,5)*p(j,5)+pmom)  ! energy (pmom=p^2) 
c         rw(part+6)	= p(j,1)		! px
c         rw(part+7)	= p(j,2)		! py
c         rw(part+8)	= p(j,3)		! pz
c         rw(part+9)	= pcharge(j)            ! charge
c         iw(part+10) 	= j			! Track pointer
c         rw(part+11)	= 0
c         rw(part+12) 	= gelewgt
c         iw(part+13) 	= heli+10*helpi   ! store helicities

c         call write2hbook(part)

c        part=part+13
c       endif
 999  CONTINUE

        mcvx = NBANK('MCVX',0,5,jj)
        if(mcvx.ne.0) then
           do j=1,jj
           rw(MCVX+1) = 0
           rw(MCVX+2) = 0
           rw(MCVX+3) = z_cible
           rw(MCVX+4) = 0.
           iw(MCVX+5) = 0
           mcvx = mcvx+5
        enddo
        endif

c        mctk = nbank('MCTK',0,11,jj)

        IF( phrad(4).gt.0.010 ) THEN
        mctk = nbank('MCTK',0,11,jj+1)
	else
        mctk = nbank('MCTK',0,11,jj)
	endif
		
        if(mctk.ne.0) then
           do j=1,jj
           pmom=sqrt(p(j,1)*p(j,1)+p(j,2)*p(j,2)+p(j,3)*p(j,3))
           rw(mctk+1)=p(j,1)/pmom
           rw(mctk+2)=p(j,2)/pmom    
           rw(mctk+3)=p(j,3)/pmom
           rw(mctk+4)=pmom
           rw(mctk+5)=p(j,5)
           rw(mctk+6)=pcharge(j)
           iw(mctk+7)=k(j,2)
           iw(mctk+8)=0
           iw(mctk+9)=1
           iw(mctk+10)=0
           iw(mctk+11)=0
           mctk = mctk + 11
           enddo
        endif

ccc-Mick
c 	print*,"hello",phrad(4)
        IF( phrad(4).gt.0.010 ) THEN
           theta_radgam=PolarTheta(PhRAD(1),PhRAD(2),PhRAD(3))
           phi_radgam=Azimuthalphi(PhRAD(1),PhRAD(2))
           rw(mctk+1)=sin(theta_radgam)*cos(phi_radgam)
           rw(mctk+2)=sin(theta_radgam)*sin(phi_radgam)
           rw(mctk+3)=cos(theta_radgam)
           rw(mctk+4)=phrad(4)
           rw(mctk+5)=0.0
           rw(mctk+6)=0.0
           iw(mctk+7)=22
           iw(mctk+8)=0
           iw(mctk+9)=1
           iw(mctk+10)=0
           iw(mctk+11)=0
	   MCTK = MCTK + 11
        ENDIF
ccc-Mick

          call goodevent()     ! check the event
c
          if(dvcsgenOK) then
            if(ntOK) call HFNT(22)
c            if (datfileOK) call  writefile(jj,pmy,kmy)
          endif
c
c
       endif
c 

c
      if(dvcsgenOK.and.bosOK) call fwbos(iw,33,'E',ierr)	! write banks to file
      call bdrop(iw,'E')     		! drop the bank, reclaim the space
      call bgarb(iw)	     		! garbage collection
c
      RETURN
      END     

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
c
           if(acceg1) then
        iaccept=eg1amcaccept(cl_beam_energy,gelee,gelet,gelef)
           if(iaccept.eq.0) return
           endif
c
           if(acc12) then
        iaccept=accep11(11,gelee,gelet,gelef,2250.0)
           if(iaccept.eq.0) return
           endif

       if (gnpro.gt.0) then 
          gprotg=gamel(0,1)   ! prot_gamm e
          gprotl=gamel(0,2)   ! prot_gamm e'
          gprotv=gamel(0,3)   ! prot_gamm virt_gamma angle
          gprotmis=mismas(4)
          gprott=tchan(4)
          gprofs=phistar(4)
        endif
        if (gngam.gt.0) then
          gproz=thetagg()
        endif 
        if (gnpi0.gt.0) then
          gprox=mismas(0)
        endif

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
         subroutine write2hbook(part)
c      implicit none
#include "bcs.inc"
#include "dvcs.inc"
#include "dvcsmom.inc"
#include "ntupgdvcs.inc"



       real el04(4),elf4(4),pro4(4)
       integer part,part_id,parent_id
       real cx,cy,cz,theta,phi,pmom,px,py,pz
       real Azimuthalphi,PolarTheta,anu,gw2,pien
       real pi0x,pi0y,pi0z

       part_id=iw(part+1)
       parent_id=iw(part+13)
       px=rw(part+6)
       py=rw(part+7)
       pz=rw(part+8)
       pmom=sqrt(px*px+py*py+pz*pz)

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
         if(gnele.eq.1 .or. gelee.lt.rw(part+5) ) then
          gelee=rw(part+5)
          geleh=iw(part+13)
          gelet=theta
          gelef=phi

          if(cl_smear.eq.1) call smear(part_id,gelee,gelet,gelef)  ! smear           
          gelez=rw(part+4)
          if(bossmearOK) call smearbos(part,part_id,gelee,gelet,gelef)
          anu=cl_beam_energy-gelee
          gq2=2.*cl_beam_energy*gelee*(1.0-cos(gelet))
          gw2=Mp*Mp+2.0*Mp*anu-gq2

          if(gw2.gt.0) then
            gw=sqrt(gw2)
          else
            gw=0
          endif

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
c        print *,'q2',gq2,vdotm(qiu4,qiu4,4)          
          gxb=gq2/2.0/Mp/anu
         endif 

        CASE (14) ! proton
         gnpro=gnpro+1
         if(gnpro.eq.1 .or. gproe.lt.rw(part+5) ) then
          gproe=rw(part+5)
          gprot=theta
          gprof=phi
          if(cl_smear.eq.1) call smear(part_id,gproe,gprot,gprof)  ! smear           
      if(bossmearOK) call smearbos(part,part_id,gproe,gprot,gprof)
         endif 

        CASE (7) ! pi0
         gpi0u=gpi0u+1
         if(gnpi0.eq.1 .or. gpi0e.lt.rw(part+5) ) then
          gpi0e=rw(part+5)
          gpi0t=theta
          gpi0f=phi
         endif 

        CASE (1) ! gamma
         gngam=gngam+1
         if(gngam.eq.1 .or. ggame.lt.rw(part+5) ) then  ! write high energy gamma as gamma-1
c
           if(gngam.gt.1) then   ! if high energy is the second store the first
            ggam2e=ggame
            ggam2t=ggamt
            ggam2f=ggamf
           endif
c
          ggame=rw(part+5)
          ggamt=theta
          ggamf=phi
          if(cl_smear.eq.1) call smear(part_id,ggame,ggamt,ggamf)  ! smear           
      if(bossmearOK) call smearbos(part,part_id,ggame,ggamt,ggamf)
         else
            ggam2e=rw(part+5)
            ggam2t=theta
            ggam2f=phi
          if(cl_smear.eq.1) call smear(part_id,ggam2e,ggam2t,ggam2f)  ! smear           
      if(bossmearOK) call smearbos(part,part_id,ggam2e,ggam2t,ggam2f)
         endif
         if(gngam.gt.1) then
          gnpi0=1
          gpi0e=ggam2e+ggame
          pi0x=game*sin(gamt)*cos(gamf)+gam2e*sin(gam2t)*cos(gam2f)
          pi0y=game*sin(gamt)*sin(gamf)+gam2e*sin(gam2t)*sin(gam2f)
          pi0z=game*cos(gamt)+gam2e*cos(gam2t)
          pi0f=Azimuthalphi(pi0x,pi0y)
          pi0t=PolarTheta(pi0x,pi0y,pi0z)
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
          real dts,ps,es,ts,fs
           es=e
           ts=t
           fs=f
           sin6=9.6                  ! sin6^o
           call rannor(rn1,rn2)
           call rannor(rn3,rn4)
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

C***********************************************************************
      SUBROUTINE CLASBOSEND(recname)
C***********************************************************************
          implicit none
#include "bcs.inc"
	character*8 recname
	integer ierr
        call fwbos(iw,33,'0',ierr)
c        call bosta
        call fparm('CLOSE recname')

         print *,'END CLASBOSEND',ierr

      RETURN
      END     


C***********************************************************************
      SUBROUTINE CLASBOSINIT(recname)
C***********************************************************************
      implicit none
      
      character*8 recname
      character*132 ddl_file
#include "bcs.inc"
#include "dvcs.inc"
c
c       Open BOS input & output files
c
      CALL bnames(1000)
      call bos(iw,nbcs)
c
      call revinm ('CLAS_PARMS','clasbanks.ddl',ddl_file)
      CALL txt2bos(ddl_file) !'clasbanks.ddl')

      close (33)

      call fparm ('OPEN ' //recname// ' UNIT=33 '// 
     &     'FILE="'//bosout//'"'// 
     &     ' RECL=32760 ACTION=WRITE STATUS=NEW FORM=BINARY')

      RETURN
      END     




      SUBROUTINE TXT2BOS(FILENAME)
c
c_begin_doc
c  RCS ID string
c  $Id: dvcsgen.F,v 1.8 2003/03/14 02:29:14 avakian Exp $
c
c  Documentation for subroutine TXT2BOS
c
c  Purpose: reads in an ascii file (FILENAME) and converts that info into
c  -------- BOS bank definition.
c
c  Input Parameters:  FILENAME :C*(*): ascii file contain BOS bank definitions
c  ----------------   Format is as follows:
c anything following a `*' or a `!' is considered a comment
c   *   this is a comment, the following is an example bank definition
c   *   leading blanks are ignored
C   * the first line should contain a 4 character  NAME followed by
c   * the name of the bank and ACTION word (CREATE, WRITE, DELETE, MODIFY)
c   *  TABLE BANKname fmt  ACTION1 ACTION2
c   * The NAME line is followed by a row ordered list of format statements
c   * as defined in page 10 of the BOS manual, for example
c   *           * run and event number 
c   *           * 2 group of 3 floating points
c   * 
c   * ACTION = CREATE implies MKFMT will be called
c   * ACTION = WRITE implies BLIST(IW,'E+','bname') will be called,
c   *          absence of WRITE implies BLIST(IW,'E-','bname') will be called
c   * ACTION = DELETE implies BLIST(IW,'R+','bname') will be called
c   * ACTION = MODIFY implies THIS definition should OVERRIDE definition in
c   *                 datafile
c   * ACTION = DISPLAY implies the bank will be `shipped' to the event store
c              buffer so that CED can display it
c
c  Output Parameters:
c  -----------------
c
c  Other routines:
c  ---------------
c
c  Notes:
c  ------
c
c  Author:   Arne Freyberger      Created:  Fri Oct  6 15:21:23 EDT 1995
c  -------
c
c  Major revisions:
c  ----------------
c     
c
c_end_doc
c
      IMPLICIT NONE
      SAVE
c
c_begin_inc
c  include files :
c  ---------------------
c BOS common block  uncomment the next line for BOS include file
#include "bcs.inc"
#include "bnkfmt.inc"
c                           CLAS control module
#include "clasmdl.inc"
c_end_inc
c
c_begin_var
c  input/output variables:
c  -----------------------
c
c  Local pre-defined variables:
c  ---------------------------
c  RCS information: 
      CHARACTER*132  CFILE, CREVIS, CSTATE, CDATE, CAUTHO, CRCSID
      PARAMETER (CFILE = '$rcsRCSfile$')
      PARAMETER (CREVIS = '$rcsRevision$')
      PARAMETER (CSTATE = '$rcsState$')
      PARAMETER (CDATE = '$rcsDate$')
      PARAMETER (CAUTHO = '$rcsAuthor$')
      DATA CRCSID/   
     1'$Id: dvcsgen.F,v 1.8 2003/03/14 02:29:14 avakian Exp $'   
     2/   
c  Module information:
      CHARACTER*(*)  CRNAME, CRAUTH
      CHARACTER*100  CRMESS
      PARAMETER (CRNAME='TXT2BOS')
      PARAMETER (CRAUTH='Arne Freyberger')
c
c  Local User defined variables:
c  -----------------------------
      INTEGER   INDEXA, INDEXN
      EXTERNAL  INDEXA, INDEXN
      INTEGER ICLOC,  ISCAN, lenocc
      EXTERNAL ICLOC, ISCAN, lenocc
      character*132 spaces, strip
      external spaces, strip
c
      CHARACTER*(*) FILENAME
      character*512 LINE
      CHARACTER*256 CNAME, CWORD, CTMP
      CHARACTER*400 BFORMAT
      CHARACTER*4 BNAME
      CHARACTER*20 CWORDS(40)
      character*128  ctest
      LOGICAL LCREATE, LMODIFY, LWRITE, LDELETE
      LOGICAL LDONE, LCLEAR, LDISPLAY, LINIT, LFORMAT_DONE
      INTEGER ILUN, IRET, NLINE
      INTEGER IN_BEG, IN_END, IFORMAT_BEG, IFORMAT_END, IBEG, IEND
      INTEGER IBEGA, IBEGN, IEND1, IEND2
      INTEGER WCNT, IWORDS, ICOLUMN, iretlen
c_end_var
      DATA LINIT/.FALSE./
      DATA BFORMAT/'                              '/
c
c  executable code for routine TXT2BOS:
c----6----------------------------------------------------------------72
c
      IF (.NOT. LINIT) THEN
        LINIT = .TRUE.
        IDELETE = 0
        IWRITE = 0
        IREMOVE = 0
        ICREATE = 0
      ENDIF

c open the ascii file
 
      CALL RESLUN(CRNAME,ILUN,IRET)
      CALL REOPEN(FILENAME,ILUN,'OLD',IRET)
      if (iret .lt. 0) return
c 
      ICOLUMN = 0

c      write(crmess,1002)
 1002 format('Bankname : Format',T45,'WRITE DELETE DISPLAY')
c      call recmes(crname,'i',crmess)

c  do stuff with the file, read in line by line and 
c  parse each line into words

      LCLEAR = .TRUE.
10    CONTINUE
      READ(ILUN,11,ERR=999,END=999)LINE
11    FORMAT(A512)
      NLINE = NLINE + 1
      LINE = SPACES(LINE,1)
      IF (LINE(1:1) .EQ. '*' .OR. 
     1    LINE(1:1) .EQ. '!') GOTO 10

c  convert to upper case

      CALL CLTOU(LINE(1:256))

c construct the BOS bank, first break the line into words
c separated by blanks

      IEND1 = ICLOC('*',1,LINE,1,256)
      IEND2 = ICLOC('!',1,LINE,1,256)
      IEND = 256
      IF (IEND1*IEND2 .NE. 0) THEN
       IEND = MIN(IEND1,IEND2)
      ELSE
       IEND = MAX(IEND1,IEND2)
      ENDIF
      IF (IEND .EQ. 0) GOTO 10
      IBEGA = INDEXA(LINE(1:IEND))
      IBEGN = INDEXN(LINE(1:IEND))
      IF (IBEGN*IBEGA .EQ. 0) THEN
        IBEG = MAX(IBEGA,IBEGN)
      ELSE
        IBEG = MIN(IBEGA,IBEGN)
      ENDIF
      if (ibeg .eq. 0) goto 10
      IF (IEND .EQ. 0) IEND = 256
c
c  separate the line into words then loop over them
c
      CALL REWORD(LINE,IWORDS,CWORDS)
      IF (IWORDS .LE. 1) GOTO 10
      CWORD = CWORDS(1)
      IF (CWORD(1:1) .EQ. '!' .OR. CWORD(1:1) .EQ. '*') GOTO 10
c
c  look for TABLE as the START of the bank definition
c
      IF (CWORD(1:2) .EQ. 'TA') THEN 
        IF (.NOT. LCLEAR .and. IFORMAT_END .gt. 1) THEN
           IN_END = IFORMAT_END - 1
           iretlen = 1
           call fmcre(bformat(1:IN_END+2), ctest, iretlen)
c           WRITE(CRMESS,1001)
c     1     BNAME,BFORMAT(1:min(23,IFORMAT_END-1)), 
c     1     BNAME,ctest(1:30), 
c     2                LWRITE, LDELETE, LDISPLAY
c           CALL RECMES(' ','N',CRMESS)
c           IF (LCREATE) THEN
               if (iretlen.LT.6) then
                  CALL MKFMT(IW,BNAME,BFORMAT(1:IN_END))
               else
                  CALL MKFMT(IW,BNAME,ctest(1:iretlen))
               endif
             ICREATE = ICREATE + 1
             CREATE_NAME(ICREATE) = BNAME
             CREATE_FORMAT(ICREATE) = BFORMAT(1:IN_END)
             NCREATE_COL(ICREATE) = ICOLUMN
c           ENDIF
           IF (LWRITE) THEN
              CALL BLIST(IW,'E+',BNAME)
              IWRITE = IWRITE + 1
              CBANK_WRITE(1+4*(IWRITE-1):4+4*(IWRITE-1)) = BNAME
           ELSE
              CALL BLIST(IW,'E-',BNAME)
              IREMOVE = IREMOVE + 1
              CBANK_REMOVE(1+4*(IREMOVE-1):4+4*(IREMOVE-1)) = BNAME
           ENDIF
           IF (LDELETE) THEN 
              CALL BLIST(IW,'R+',BNAME)
              IDELETE = IDELETE + 1
              CBANK_DELETE(1+4*(IDELETE-1):4+4*(IDELETE-1)) = BNAME
           ENDIF
           IF (LDISPLAY) THEN
            NBANK_DISPLAY = NBANK_DISPLAY + 1
            CBANK_DISPLAY(NBANK_DISPLAY) = BNAME
            CBANK_FORMAT(NBANK_DISPLAY) = BFORMAT(1:IN_END)
           ENDIF
           LDONE = .FALSE.
           LFORMAT_DONE = .FALSE.
           IFORMAT_END = 0
           ICOLUMN = 0
           BNAME = '    '
           BFORMAT = '                '
        ENDIF

c  Done writting out previous bank now start over and initialize logicals

        LCLEAR = .FALSE.
        LCREATE = .FALSE.
        LDELETE = .FALSE.
        LWRITE = .FALSE.
        LMODIFY = .FALSE.
        LDISPLAY = .FALSE.
        LFORMAT_DONE = .FALSE.

c BANK name is the SECOND word!!

        BNAME = CWORDS(2)
c
        IF (BNAME .NE. '    ') THEN
c
        DO 23 WCNT=3,IWORDS
         CNAME = cwords(wcnt)
         IF (CNAME(1:3) .EQ. 'B16' .OR.
     1       CNAME(1:3) .EQ. 'B08' .OR.
     2       CNAME(1:3) .EQ. 'B32') THEN 
          BFORMAT = CNAME
          IFORMAT_END = LENOCC(BFORMAT) + 1
          LFORMAT_DONE = .TRUE.
         ELSEIF (CNAME(1:6) .EQ. 'CREATE') THEN
          LCREATE = .TRUE.
         ELSEIF (CNAME(1:5) .EQ. 'WRITE') THEN
          LWRITE = .TRUE.
         ELSEIF (CNAME(1:6) .EQ. 'DELETE') THEN
          LDELETE = .TRUE.
         ELSEIF (CNAME(1:6) .EQ. 'MODIFY') THEN
          LMODIFY = .TRUE.
         ELSEIF (CNAME(1:7) .EQ. 'DISPLAY') THEN
          LDISPLAY = .TRUE.
         ENDIF
23      CONTINUE
        ENDIF

c   check for END or END TABLE presence

      ELSEIF (CWORD(1:3) .EQ. 'END' .AND. BNAME .NE. '    ') THEN
          LDONE = .TRUE.
c
c  the following lines parse the FORMAT lines (not a TAble or END line)
c
      ELSEIF (BNAME .NE. '    ') THEN
         IFORMAT_BEG = IFORMAT_END + 1
         CNAME = CWORDS(3)
         IN_END = ISCAN(CNAME,' ') - 1
         if (in_end .eq. 0) goto 10
c
c  reject relational formats
c
         IF (CNAME(1:1) .EQ. 'D' .OR.
     1       CNAME(1:1) .EQ. 'R' .OR.
     2       CNAME(1:1) .EQ. 'M') GOTO 10
         ICOLUMN = ICOLUMN + 1
         ctmp = cwords(2)
         create_elements(icreate+1,icolumn)= ctmp(1:8)
         if (.not. lformat_done) then
          IN_BEG = 1
          IFORMAT_END = IFORMAT_BEG + IN_END - 1
          BFORMAT(IFORMAT_BEG:IFORMAT_END) = CNAME(1:IN_END)
          IFORMAT_END = IFORMAT_END + 1
          BFORMAT(IFORMAT_END:IFORMAT_END) = ','
         endif
       ENDIF
c
       IF (LDONE .AND. BNAME .NE. '    ') THEN
           IN_END = IFORMAT_END - 1
           CALL MKFMT(IW,BNAME,BFORMAT(1:IN_END))
           ICREATE = ICREATE + 1
           CREATE_NAME(ICREATE) = BNAME
           CREATE_FORMAT(ICREATE) = BFORMAT(1:IN_END)
           NCREATE_COL(ICREATE) = ICOLUMN

           IF (LWRITE) THEN
              CALL BLIST(IW,'E+',BNAME)
              IWRITE = IWRITE + 1
              CBANK_WRITE(1+4*(IWRITE-1):4+4*(IWRITE-1)) = BNAME
           ELSE
              CALL BLIST(IW,'E-',BNAME)
              IREMOVE = IREMOVE + 1
              CBANK_REMOVE(1+4*(IREMOVE-1):4+4*(IREMOVE-1)) = BNAME
           ENDIF
           IF (LDELETE) THEN 
              CALL BLIST(IW,'R+',BNAME)
              IDELETE = IDELETE + 1
              CBANK_DELETE(1+4*(IDELETE-1):4+4*(IDELETE-1)) = BNAME
           ENDIF
           IF (LDISPLAY) THEN
            NBANK_DISPLAY = NBANK_DISPLAY + 1
            CBANK_DISPLAY(NBANK_DISPLAY) = BNAME
            CBANK_FORMAT(NBANK_DISPLAY) = BFORMAT(1:IN_END)
           ENDIF
           iretlen = 1
           call fmcre(bformat, ctest, iretlen)
c           WRITE(CRMESS,1001)
c     1     BNAME,BFORMAT(1:min(23,IFORMAT_END-1)), 
c     1     BNAME,ctest, 
c     2                LWRITE, LDELETE, LDISPLAY
 1001      FORMAT('|',13x,A4,' : ',A30,' : ',T60,3(L1,6X))
c           CALL RECMES(' ','N',CRMESS)
           LDONE = .FALSE.
           LFORMAT_DONE = .FALSE.
           LCLEAR = .TRUE.
           IFORMAT_END = 0
           ICOLUMN = 0
           BNAME = '    '
           BFORMAT = '            '
       ENDIF
      GOTO 10
c
c  but format last bank before exiting
c
999   CONTINUE
      if (.not. lclear) then
       iretlen = 1
       call fmcre(bformat, ctest, iretlen)
c       WRITE(CRMESS,1001)
c     1     BNAME,BFORMAT(1:min(23,IFORMAT_END-1)), 
c     1     BNAME,ctest, 
c     2                LWRITE, LDELETE, LDISPLAY
c       CALL RECMES(' ','N',CRMESS)
       if (iretlen.LT.6) then
          CALL MKFMT(IW,BNAME,BFORMAT(1:IFORMAT_END-1))
       else
          CALL MKFMT(IW,BNAME,ctest(1:iretlen))
       endif
       ICREATE = ICREATE + 1
       CREATE_NAME(ICREATE) = BNAME
       CREATE_FORMAT(ICREATE) = BFORMAT(1:IFORMAT_END-1)
       NCREATE_COL(ICREATE) = ICOLUMN

       IF (LWRITE) THEN
           CALL BLIST(IW,'E+',BNAME)
           IWRITE = IWRITE + 1
           CBANK_WRITE(1+4*(IWRITE-1):4+4*(IWRITE-1)) = BNAME
       ELSE
           CALL BLIST(IW,'E-',BNAME)
           IREMOVE = IREMOVE + 1
           CBANK_REMOVE(1+4*(IREMOVE-1):4+4*(IREMOVE-1)) = BNAME
       ENDIF
       IF (LDELETE) THEN 
           CALL BLIST(IW,'R+',BNAME)
           IDELETE = IDELETE + 1
           CBANK_DELETE(1+4*(IDELETE-1):4+4*(IDELETE-1)) = BNAME
       ENDIF
       IF (LDISPLAY) THEN
         NBANK_DISPLAY = NBANK_DISPLAY + 1
         CBANK_DISPLAY(NBANK_DISPLAY) = BNAME
         CBANK_FORMAT(NBANK_DISPLAY) = BFORMAT(1:IN_END)
       ENDIF
c       iretlen = 30
c       call fmcre(bformat, ctest, iretlen)
c       WRITE(CRMESS,1001)
c     1     BNAME,BFORMAT(1:min(23,IFORMAT_END-1)), 
c     1     BNAME,ctest, 
c     2                LWRITE, LDELETE, LDISPLAY
c       CALL RECMES(' ','I',CRMESS)
       IFORMAT_END = 0
      endif
c      write(crmess,1002)
c      call recmes(crname,'i',crmess)
      CALL RESLUN(CRNAME,-ILUN,IRET)
c
c$$$      call recmes(crname,'i',
c$$$     1  'The following banks will be DELETED from memory before reading in the next event')
c$$$      if (idelete .ne. 0) 
c$$$     1      call recmes(crname,'i',cbank_delete(1:4+4*(IDELETE-1)))
c$$$      call recmes(crname,'i',
c$$$     1      'The following banks will NOT be written out')
c$$$      if (iremove .ne. 0) 
c$$$     1      call recmes(crname,'i',cbank_remove(1:4+4*(IREMOVE-1)))
c$$$      call recmes(crname,'i',
c$$$     1      'The following banks WILL be written out')
c$$$      if (iwrite .ne. 0) 
c$$$     1      call recmes(crname,'i',cbank_write(1:4+4*(IWRITE-1)))
      CLOSE(ILUN)
      RETURN
      END
c
c$$$c------------------------------------------------------------------------------
c$$$      COMMON /LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)      
c$$$      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
c$$$      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)

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
       real vm
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
ccc- Hyon-Suk
      common/radmick/PhRAD
      real PhRAD(4)
      common/radhs/radtheta
      real radtheta
ccc- Hyon-Suk
c
c
      istatus=0                    ! kinematic range OK
      Ed=cl_beam_energy
      E=Ed
ccc- Hyon-Suk
c      if((PhRAD(4).gt.0.).and.(radtheta.le.0.15)) E=cl_beam_energy-PhRAD(4)
ccc- Hyon-Suk
      xmin1 = Q2/(2D0*Mp*E)
      xmax1 = 1D0
      nu  = Q2/(2D0*Mp*xb)
      W2  = Mp**2 + 2D0*Mp*nu - Q2
c
      if( xb.le.xmin1 .or. xb.gt.xmax1 ) istatus=1  !    x  out of range
      if( Q2.le.cl_q2min .or. Q2.gt. cl_q2max) istatus=2  !    Q2  out of range
c      yb=nu/Ed
ccc- Hyon-Suk
      yb=nu/E
ccc- Hyon-Suk
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
      integer istat
      double precision  xmin1,xmax1
      double precision dsIunp,dsBHlp, dsIlp,dsBHtp,dsItp
      double precision temp_th_g
ccc- Hyon-Suk
      common/radmick/PhRAD
      real PhRAD(4)
      common/radhs/radtheta
      real radtheta
ccc- Hyon-Suk
c
c
      istat=0                    ! kinematic range OK
      Ed=cl_be
      E=Ed
ccc- Hyon-Suk
      if((PhRAD(4).gt.0.).and.(radtheta.le.0.15)) E=cl_be-PhRAD(4)
ccc- Hyon-Suk
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
      if( xb.le.xmin1 .or. xb.gt.xmax1 ) istat=1           !    x  out of range
      if( del2.ge.del2min .or. del2.le.del2max ) istat=2   ! delta out of range
c      yb=nu/Ed
ccc- Hyon-Suk
      yb=nu/E
ccc- Hyon-Suk
      if(yb.gt.cl_ymax) istat=3                               ! y<ymax
      call dvcsycol(del2,xb,Q2,ycol)
	
      call getphoton(xb,Q2,del2,phi_e,phi_g)
      temp_th_g=acos(V3gam(3)/(sqrt(V3gam(1)*V3gam(1)+V3gam(2)*V3gam(2)+V3gam(3)*V3gam(3))))
      if(acos((v3k2(1)*v3gam(1)+v3k2(2)*v3gam(2)+v3k2(3)*v3gam(3))/(sqrt(V3gam(1)*V3gam(1)+V3gam(2)*V3gam(2)+V3gam(3)*V3gam(3))*sqrt(V3k2(1)*V3k2(1)+V3k2(2)*V3k2(2)+V3k2(3)*V3k2(3)))).lt.0.034) istat=4
      if(temp_th_g.lt.0.034) istat=4
c      if( (ycol-yb).le.cl_ycol ) istat=4                   ! y-too big

      if (istat.eq.0) then
      Phi_gb=pi - Phi_g
      Phi_s=Phi_e
      call bhdvcs(xb,  Q2,  del2,  Phi_s,Phi_gb) 
c
      dsBH =hc0BH +hc1BH*cos(Phi_gb)+hc2BH*cos(2D0*Phi_gb)
      dsDVCS=hc0dvcs+hc1dvcs*cos(Phi_gb)+hs1dvcs*sin(Phi_gb)+hs2dvcs*sin(2*Phi_gb)
      dsIunp=hc0Iunp +heli*hs1Iunp*sin(Phi_gb)+hc1Iunp*cos(Phi_gb) 
     6     +heli*hs2Iunp*sin(2*Phi_gb)+hc2Iunp*cos(2*Phi_gb)
c      
       dsigma=dsBH +dsDVCS-Ich*dsIunp 

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
       gsin=hs1Iunp/(hc0BH+hc0dvcs+hc0Iunp)
       gsin2=hc0BH+hc0dvcs+hc0Iunp 
        gproh=helpi
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
      real HTre2,UFx,PFxq2,Tfac,u,aac2s0
      integer istat
ccc- Hyon-Suk
      common/radmick/PhRAD
      real PhRAD(4)
      common/radhs/radtheta
      real radtheta
ccc- Hyon-Suk
c
c
      istat=0                    ! kinematic range OK
ccc- Hyon-Suk
      Ed=cl_be
ccc- Hyon-Suk
      E=Ed
ccc- Hyon-Suk
      if((PhRAD(4).gt.0.).and.(radtheta.le.0.15)) E=cl_be-PhRAD(4)
ccc- Hyon-Suk
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
      if( xb.le.xmin1 .or. xb.gt.xmax1 ) istat=1           !    x  out of range
      if( del2.ge.del2min .or. del2.le.del2max ) istat=2   ! delta out of range
c      yb=nu/Ed
ccc- Hyon-Suk
      yb=nu/E
ccc- Hyon-Suk
      if(W2.lt.cl_wmin) print *,'Debiles lriv!!!!!!!!!',W2
      if (istat.eq.0) then
c
       sinpi0=2*xb*(1-xb)**5
       PFxq2=(1-xb)**2/(Mp*Mp*xb/Q2+1-xb)**2
       Tfac=exp(4*(del2-del2min))
c
c      aac2s0(xb,Q2)-leading order pdf from http://spin.riken.bnl.gov/aac/ (hep-ph/0001046)
c
       skew= xb/(2.0 - xb)
       UFx=aac2s0(skew,Q2) !HTre2(xb,del2min) !(du(xi)-dd(xi))**2
       dpi0u=0.3*Tfac*UFx*PFxq2
       dpi0u=dpi0u/Q2/(Q2+Mp**2)**2*(1.0+sinpi0*heli*sin(Phi_g))
       dsigma=dpi0u
c       print *,xb,u(xb)*xb,UFx,(u(xb)*xb)**2,UFx/(u(xb)*xb)**2
c       print *,'exc->',smax,q2,del2,'->',xb,PFxq2,UFx,dsigma
c
c      L-POL
c
       if(cl_pol.eq.1) then 
        dpi0l=0 !????    ! longpol
        dsigma=dsigma+dpi0l
       endif
c
c      T-POL
c 
       if(cl_pol.eq.2) then 
         dpi0t=0 !???? ! tpol
         dsigma=dsigma+dpi0t
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
      subroutine bhdvcs(xbd, Q2d, del2d, phield, phigd)
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
      double precision nu, k1pl, k2pl, Kfac, Jfac
      common/todvcs/ xb, yb, Q2, del2, del2min,phip,phipel,P1,P2,Kfac,ds
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
      c11_BHlp=2D0*tau-xbd*del2q2m1*(1D0-xbd+xbd*del2q2)
      c12_BHlp=1D0+xbd-(3D0-2D0*xbd)*(1D0+xbd*del2q2)-xtau*(1D0-del2q4)
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
       hc0BH=c0_BH*BHfact
       gh1dvcs=c0_BH                    !temporary
       hc1BH=c1_BH*BHfact
       hc2BH=c2_BH*BHfact
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
      double precision Kfac,Intfac
      common/todvcs/ xb, yb, Q2, del2, del2min,phip,phipel,P1,P2,Kfac,ds
      common/formfac/ F1pn(2), F2pn(2)
c
      skew = xb/(2D0 - xb)
      call amptab(skew, del2, 
     &             H1_RE, H1_IM, H1T_RE, H1T_IM,
     &             E1_RE, E1_IM, E1T_RE, E1T_IM )
*
     
c	H1_RE=H1_RE*(1.-del2*10.)
c	H1_IM=H1_IM*(1.-del2*10.)
c	H1T_RE=H1T_RE*(1.-del2*10.)
c	H1T_IM=H1T_IM*(1.-del2*10.)
c       E1_RE=E1_RE*(1.-del2*10.)
c	E1_IM=E1_IM*(1.-del2*10.)
c	E1T_RE=E1T_RE*(1.-del2*10.)
c	E1T_IM=E1T_IM*(1.-del2*10.)

C         proton/neutron
      F1 = F1pn(Ipn)
      F2 = F2pn(Ipn)
*
      deldel    = 1D0 - del2min/del2
      deldel_sq = sqrt(deldel)
      del2m2    = -del2/Mp**2
      del2m4    = -del2m2/4D0
      delm2_sq  = sqrt(del2m2)
      cy2   = 2D0 - 2D0*yb + yb**2
      Intfac=ds/(xb*yb**3*P1*P2*(-del2))
*
C  DVCS
C     
      a1 = H1_RE**2 + H1_IM**2 + H1T_RE**2 + H1T_IM**2
      a2 = 2.*( H1_RE*E1_RE  +  H1_IM*E1_IM +
     &         H1T_RE*E1T_RE + H1T_IM*E1T_IM )
      a3 =  E1_RE**2 +  E1_IM**2
      a4 = E1T_RE**2 + E1T_IM**2 
      C_DVCS = ( 4D0*(1D0-xb)*a1 - a2*xb**2 - (xb**2 + (2D0-xb)**2*del2m4)*a3 
     &           - xb**2*del2m4*a4 )/(2D0 - xb)**2
      C_DVCS_eff = -xb*C_DVCS
      c0_DVCS = 2D0*cy2*C_DVCS
      c1_DVCS = 8D0*((2D0 - yb)/(2D0 - xb))*C_DVCS_eff
      T_DVCS  = (c0_DVCS + Kfac*c1_DVCS*cos(phip))/(yb**2*Q2)
      DVCSfac=ds/(yb**2*Q2)
      hc0dvcs=c0_DVCS*DVCSfac
      hc1dvcs=Kfac*c1_DVCS*DVCSfac
C
C  INTERF
C
      C_I_re = F1*H1_RE + xb/(2D0-xb)*(F1+F2)*H1T_RE - del2m4*F2*E1_RE
      C_I_im = F1*H1_IM + xb/(2D0-xb)*(F1+F2)*H1T_IM - del2m4*F2*E1_IM
      RE2    = xb/(2D0-xb)*(H1_RE + E1_RE) + H1T_RE
      C_I_re_eff = -xb*C_I_re
      C_I_im_eff = -xb*C_I_im
      b1 = (2D0 - xb)*(1D0 - yb) 
      b1= b1 + (2D0-yb)**2/(1D0-yb)*Kfac*Kfac/del2*Q2 ! old - (1D0 - xb)*(2D0 - yb)**2*deldel
      b2 = (1D0 - yb)*xb*(F1 + F2)
      c0_I = -8D0*(2D0 - yb)*( b1*C_I_re - b2*RE2 )
      c1_I = -8D0*cy2*C_I_re
      s1_I =  8D0*yb*(2D0 - yb)*C_I_im
      c2_I = -16D0*((2D0 - yb)/(2D0 - xb))*C_I_re_eff
      s2_I =  16D0*(yb/(2D0 - xb))*C_I_im_eff
c
C     moments
c
      hs2Iunp=Kfac*Kfac*s2_I*Intfac
      hs1Iunp=Kfac*s1_I*Intfac
      ghs1dvcs=Kfac*s1_I                      !temporary
      hc2Iunp=Kfac*Kfac*c2_I*Intfac
      hc1Iunp=Kfac*c1_I*Intfac
      hc0Iunp=del2/Q2*c0_I*Intfac
c
c      print *,'mysl-dvcs',hs1Iunp,hc1Iunp,hc0Iunp
C
C  LPOL
C
      C_LP_re = (F1+F2)*skew*(H1_RE+xb/2D0*E1_RE)+
     6          F1*H1T_RE-skew*(xb/2D0*F1+del2m4*F2)*E1T_RE
c
c
      C_LP_im = (F1+F2)*skew*(H1_IM+xb/2D0*E1_IM)+
     6          F1*H1T_IM-skew*(xb/2D0*F1+del2m4*F2)*E1T_IM
c
      DC_LP_re=-skew*(F1+F2)*(H1_RE+xb/2D0*E1_RE+skew
     6*(H1T_RE+xb/2D0*E1T_RE))
      DC_LP_im=-skew*(F1+F2)*(H1_IM+xb/2D0*E1_IM+skew
     6*(H1T_IM+xb/2D0*E1T_IM))
c
      yf2=(2D0-yb)**2/(1D0-yb)+2
      c0lp_I =-8*yb*(Kfac**2*yf2*C_LP_re+(1D0-yb)*(2D0-xb)*del2/Q2
     6        *(DC_LP_re+C_LP_re))
      s1lp_I = 8*Kfac*cy2*C_LP_im
      c1lp_I = -8*Kfac*yb*(2D0-yb)*C_LP_re

C
C  TPOL
C
      xb2=xb*skew
      C_TPP_re_s = (F1+F2)*(xb2*(H1_RE+xb/2D0*E1_RE)+xb*del2m4*E1_RE)
     6          -xb2*F1*(H1T_RE+xb/2D0*E1T_RE)
      C_TPP_re_b1 =del2m4*4D0*(1D0-xb)/(2D0-xb)*F2*H1T_RE
      C_TPP_re_b21 =-del2m4*xb*F1*E1T_RE
      C_TPP_re_b22 =-del2m4*xb2*F2*E1T_RE
      C_TPP_re_b =C_TPP_re_b1+C_TPP_re_b21+C_TPP_re_b22
c
      C_TPP_im_s = (F1+F2)*(xb2*(H1_IM+xb/2D0*E1_IM)+xb*del2m4*E1_IM)
     6          -xb2*F1*(H1T_IM+xb/2D0*E1T_IM)
      C_TPP_im_b =del2m4*(4D0*(1D0-xb)/(2D0-xb)*F2*H1T_IM
     6            -(xb*F1+xb2*F2)*E1T_IM)
c
c
c
      C_TPM_re_s = 1D0/(2D0-xb)*(xb*xb*F1-(1D0-xb)*4D0*del2m4*F2)*H1_RE
     6          -xb2*(F1+F2)*(H1T_RE+del2m4*E1T_RE)
      C_TPM_re_b = (del2m4*((2D0-xb)*F1+xb2*F2) +xb2*F1)*E1_RE
c
      C_TPM_im_s = 1D0/(2D0-xb)*(xb*xb*F1-(1D0-xb)*4D0*del2m4*F2)*H1_IM
     6          -xb2*(F1+F2)*(H1T_IM+del2m4*E1T_IM)
      C_TPM_im_b = (del2m4*((2D0-xb)*F1+xb2*F2) +xb2*F1)*E1_IM
c
c
      C_TPM_re=C_TPM_re_s+C_TPM_re_b
      C_TPP_re=C_TPP_re_s+C_TPP_re_b
      C_TPM_im=C_TPM_im_s+C_TPM_im_b
      C_TPP_im=C_TPP_im_s+C_TPP_im_b
c
c
c
      DC_TPP_re=-4D0*del2m4*(F2*H1T_RE-xb/(2D0-xb)*(F1+xb/2D0*F2)*E1T_RE)
      DC_TPM_re=4D0*del2m4*(F2*H1_RE-F1*E1_RE)
c
      DC_TPP_im=-4D0*del2m4*(F2*H1T_IM-xb/(2D0-xb)*(F1+xb/2D0*F2)*E1T_IM)
      DC_TPM_im=4D0*del2m4*(F2*H1_IM-F1*E1_IM)
c
c
c
      qm8=8D0*Mp*sqrt(1D0-yb)/sqrt(Q2)
      c0tpcos_I = -qm8*Kfac*yb*(((2D0-yb)**2/(1D0-yb)+2D0)*C_TPP_re+DC_TPP_re)
      c0tpsin_I = qm8*Kfac*(2D0-yb)*((2D0-yb)**2/(1D0-yb)*C_TPM_im+DC_TPM_im)
c
      c1tpcos_I = -qm8*yb*(2D0-yb)*C_TPP_re
      c1tpsin_I = qm8*cy2*C_TPM_im
c
      s1tpcos_I = qm8*cy2*C_TPP_im
      s1tpsin_I = -qm8*yb*(2D0-yb)*C_TPM_re

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
      call amptab(skew, del2, 
     &             H1_RE, H1_IM, H1T_RE, H1T_IM,
     &             E1_RE, E1_IM, E1T_RE, E1T_IM )
C
      HTre2 = H1T_RE*H1T_RE+ H1T_IM*H1T_IM
      return
      end

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
      H1_RE = (4D0*(F1u/2D0)*h1ure + F1d*h1dre)/9D0
      H1_IM = (4D0*(F1u/2D0)*h1uim + F1d*h1dim)/9D0
          else
      H1_RE = (4D0*h1ure + h1dre)/9D0
      H1_IM = (4D0*h1uim + h1dim)/9D0
          endif
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
* 
      open(unit=11,file='gpd.dat',status='old')
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
      parameter (Mv = 0.84D0, kp = 1.79D0, kn = -1.91D0)
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
c
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
      else
       pi0mom=sqrt(Egam*Egam-metta*metta)
       costgg = (2D0*Egam*(Mp + nu) + Q2d - 2D0*Mp*nu-metta*metta)
      endif
      costgg = costgg/(2D0*pi0mom*qmod)
      if(abs(costgg).gt.1) then
c       print *,'costgg-------->',costgg   ! shit happens
       costgg=1.0
       sintgg=0.0
      else
        sintgg = sqrt(1D0 - costgg**2)
      endif

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

       CALL LUEXEC   ! decay pi0 

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
          print *,'Unknown PID!!!! ',k(3,2)
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
c*

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
*
      dimension V3k1(3), V3q(3)
      dimension  V3p1(3), V3plus(3)
      double precision V3k1,V3q,V3p1,V3plus
      double precision xbd,Q2d,del2d,phield,phigd
      double precision yb,nu,Esc,sintel,costel,costVq,sintVq,qmod
      double precision Ep,Egam,costgg,sintgg,Vgx,Vgy,Vgz
      double precision costeg,teteg,coste1g,tete1g,V3dotd
      double precision cosphe,sinphe
      integer j
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
c      if(Ed.ne.cl_beam_energy) print *,'dans getphoton : Ed=',Ed
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
      if(cl_radgen) then
       if (PhRAD(4).gt.0.) then
c        print *,'getphoton : radtheta=',radtheta
c        print *,'!!!!dans getphoton: PhRAD = ',PhRAD
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
      N=2
       DO j=1,3
        p(2,j)=V3gam(j)
       enddo
        p(2,5)=0.0      ! gamma
        k(2,2)=22   ! photon

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


      character*10 cl_target 
      REAL cl_beam_energy,cl_xpos,cl_ypos,cl_zpos,cl_zwidth,cl_be
      REAL  cl_ymin,cl_ymax,cl_tmin,cl_tmax,cl_thmin,cl_thmax
      REAL  cl_ycol,cl_scale,cl_bpi0,cl_dmin,cl_dmax,cl_wmin
      REAL  cl_q2min,cl_q2max,cl_xmin,cl_xmax
      integer cl_verblev,cl_rand,cl_nprint,cl_writef
      integer cl_triggers,cl_nmax,cl_pol,cl_gpd,cl_bh,cl_mod
      integer  cl_smear,cl_proloss,cl_seed
      real  cl_sma,cl_smb,cl_smc,cl_smd
      real  cl_sme,cl_smf,cl_smg,cl_rast
      logical dvcsgenOK,bosOK,ntOK,weightOK,cl_mom,cl_dvcs,cl_pi0,cl_eta
      logical bossmearOK,dalOK,cl_radgen,datfileOK,acce16,acceg1,acc12
      logical cl_ktcor,cl_printgpd,cl_docker
      character*7 bosout
      common /OUT_NAMES/ bosout
      common /clasdvcs1/ cl_verblev,cl_nprint,cl_proloss,cl_writef,cl_seed
      common /clasdvcs2/ cl_triggers,cl_gpd,cl_bh,cl_target,cl_nmax
      common /clasdvcs3/ cl_beam_energy,cl_pol,cl_be,cl_scale,cl_bpi0
      common /clasdvcs4/ cl_xpos, cl_ypos, cl_zpos,cl_zwidth
      common /clasdvcs5/ cl_ymin,cl_ymax,cl_tmin,cl_tmax,cl_thmin,cl_thmax
      common /clasdvcs6/ cl_ycol,cl_dmin,cl_dmax,cl_wmin
      common /clasdvcs7/ cl_q2min,cl_q2max,cl_xmin,cl_xmax
      common /clasdvcs8/ cl_smear,cl_sma,cl_smb,cl_smc,cl_smd
      common /clasdvcs9/ cl_sme,cl_smf,cl_smg,cl_rast
      common /dvcsgencont/ dvcsgenOK,bosOK,ntOK,weightOK,bossmearOK,
     6cl_mom,cl_dvcs,cl_pi0,cl_eta,cl_radgen,cl_mod,
     6datfileOK,acce16,acceg1,acc12,cl_ktcor,cl_printgpd,cl_docker
C..
      double precision Mp, mele, pi,mpi0,metta
      parameter (Mp=0.93827D0)
      parameter (mele=0.000511D0)
      parameter (pi=3.1415926536D0)
      parameter (mpi0=0.135)
      parameter (metta=0.5473)
C..
       double precision Ed,dsBH, dsDVCS, dsINT
       INTEGER Ipn,Ich
       COMMON/DVDAT2/Ed,Ipn,Ich
C..
        double precision smax
        INTEGER istatus
        common /smax/smax,istatus
C..
      integer heli,helpi
      double precision V3k2(3), V3gam(3), V3gam2(3), V3p2(3),V3el(3),V3po(3)
      COMMON/gendvcs/V3k2,V3gam,V3gam2,V3p2,V3el,V3po,heli,helpi
c
C..
        integer  Ifl, Ivar, IGPD
        common /igpd/ Ifl, Ivar, IGPD 



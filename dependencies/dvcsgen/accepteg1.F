c
        integer function eg1amcacceptgamma(gamx,gamy,gamz)
        implicit none
        real PolarTheta,Azimuthalphi
        real gamx,gamy,gamz
        real game,gamt,gamf
        real rad2deg,gamtd,gamtf
        parameter (rad2deg=57.29577951)
c
c
        gamt=PolarTheta(gamx,gamy,gamz)
        gamf=Azimuthalphi(gamx,gamy)
        gamtd=gamt*rad2deg
        if(gamtd.gt.7.and.gamtd.lt.45) then
         eg1amcacceptgamma=1
        else
         eg1amcacceptgamma=0
        endif
        return
        end
c
c
c
        integer function eg1amcaccept(ebeam,gelee,gelet,gelef)
        implicit none
c
c
        integer icutei,qe,weight_e
        real ebeam,gelee,gelet,gelef
        real pel,secte,thetaeld,phiel,torcur,del_phie
        real getsec,fittheta,fitphi
c
        real p_el,th_tr1_EL,ph_tr1_EL,sec_EL
c
        sec_EL=getsec(gelef)
        p_el=gelee
        th_tr1_EL=(gelet)*180/3.1415-fittheta(gelee)
        ph_tr1_EL=(gelef)*180/3.1415-fitphi(gelee)
        if(ph_tr1_EL.gt.360) ph_tr1_EL=ph_tr1_EL-360
        if(ph_tr1_EL.lt.0) ph_tr1_EL=ph_tr1_EL+360
c
        pel=p_el
        secte=sec_EL
        thetaeld=th_tr1_EL
        if(secte.gt.1) then
          phiel=ph_tr1_EL-(secte-1)*60+30
        else
          if (phiel.gt.330) then
           phiel=ph_tr1_EL-330
          else
           phiel=ph_tr1_EL+30
          endif       
        endif
        icutei=3
        torcur=2250.0
        qe=-1
        del_phie=0
c
        call elacc_newnew(ebeam,secte,pel,thetaeld,phiel
     6,torcur,qe,icutei,del_phie,weight_e)
c
        eg1amcaccept=weight_e
        return
        end
c
c
       real function getsec(f)
        real f,s,fd
        s=1
        fd=f*180/3.1415926
        if(fd.gt.30.and.fd.le.90) then
          s=2
        else if (fd.gt.90.and.fd.le.150) then
          s=3
        else if (fd.gt.150.and.fd.le.210) then
          s=4
        else if (fd.gt.210.and.fd.le.270) then
          s=5
        else if (fd.gt.270.and.fd.le.330) then
          s=6
        endif
        getsec=s
        return
        end
c
        real function fittheta(p)
        real p,a0,a1,a2,a3
        a0=0.6877
        a1=-0.584
        a2=0.2577
        a3=-0.058
        fittheta=(a0+a1*p)/p/p+a2+a3*p
        return
        end

c
        real function fitphi(p)
        real p,a0,a1,a2
        a0=-11.73
        a1=0.5459
        a2=-0.0743
        fitphi=a0/p+a1+a2*p
        return
        end

c
        real function hfittheta(p)
        real p,a0,a1,a2,a3
        a0=0.3263
        a1=-0.07144
        a2=0.001945
        hfittheta=(a0+a1*p)/p/p+a2
        return
        end


c
        real function hfitphi(p)
        real p,a0,a1,a2
        a0=12.53
        a1=-1.249
        a2=0.2180
        hfitphi=a0/p+a1+a2*p
        return
        end

c
        integer function eg1amcaccepthad(part_id,pare,part,parf)
        implicit none
c
        real pare,part,parf
        integer icuthi,icutho,q_h,secth,part_id,weight_h
        real p_h,thetahd,phih,phi_m,torcur
     6,del_phih
        real getsec,hfittheta,hfitphi,pmom
        real hadphi,hadtheta
        real amp2,ampi2,amka2
        parameter (amp2=0.88035,ampi2=0.019488,amka2=0.24371)
        real rad2deg
        parameter (rad2deg=57.29577951)
c
        real p_PI,th_tr1_PI,ph_tr1_PI,sec_PI
c
          hadtheta=part
          hadphi  =parf
          sec_PI=getsec(parf)       
c
         SELECT CASE (part_id)  
         CASE (2212) ! pr
          q_h=1
          pmom=sqrt(pare*pare-amp2)
         CASE (211) ! pi+
           q_h=1
           pmom=sqrt(pare*pare-ampi2)
         CASE (-211) ! pi-
           q_h=-1
           pmom=sqrt(pare*pare-ampi2)
         CASE DEFAULT
         END SELECT

        p_pi=pmom
        th_tr1_PI=(hadtheta)*rad2deg-hfittheta(pmom)
        ph_tr1_PI=(hadphi)*rad2deg-hfitphi(pmom)
        if(ph_tr1_PI.gt.360) ph_tr1_PI=ph_tr1_PI-360
        if(ph_tr1_PI.lt.0) ph_tr1_PI=ph_tr1_PI+360
c
        p_h=p_pi
        secth=sec_PI
        thetahd=th_tr1_PI
        phih=ph_tr1_PI
c
        icuthi=2
        icutho=2
        torcur=2250.0
        del_phih=0
c
        call hadacc_eg1(secth,p_h,thetahd,phih,torcur,
     *   q_h,icuthi,icutho,del_phih,phi_m,weight_h)
c
        eg1amcaccepthad=weight_h
        return
        end
c
c
        subroutine hadacc_eg1(secth,p_h,thetahd,phih,torcur,
     *                      q_h,icuthi,icutho,del_phih,phi_m,weight_h)

* Version 1.0 / December 12, 1999

* secth: sector where hadron is detected
* p_h: hadron momentum, after corrections
* thetahd: hadron scattering angle in degrees at dc layer 1
* torcur: torus current (e.g. -1500.)
* icuthi/o defines tightness of fiducial cuts for in/out-bending hadrons

*    ==>      these are currently not implemented!!!
* icuthi = 1  tight cut (~100% efficiency inside fiducial region)
* icuthi = 2  loose cut (includes some of the fringes of the CC)
* icutho = 1 tight cut for outbending hadrons
* icutho = 2 loose cut for outbending hadrons
* del_phih: hadron acceptance in phi (degrees) for each event
* weight_h: =1 if hadron is accepted,  0. if not accepted


        real ebeam, p_h,thetahd,phih, phih_mod, torcur,del_phih,
     *       thcut, pnorm,phi_m
        integer secth,q_h,weight_h


        phih_mod =  mod(phih+30., 60.)

cDefine acceptance cuts for hadrons

        weight_h=0
        del_phih= 0.

       If(torcur/q_h.gt.0.)then
          thmin = 10.  + 25.*(1.-p_h*3375./abs(torcur)/5.)**24
          th2 = 52.*(p_h/3.)**(1/16.)
          th3 = 75.
          th4 = 100.
          if(thetahd.gt.thmin.and.thetahd.lt.th2)then
             del_phih = 27.*cos((thetahd-65.)*3.14159/180.)**0.4*
     *                          (p_h*3375./torcur/5.)**(1./16.)
             phi_m    = 0
          else
             del_phih = 0.
             phi_m    = 0
          endif

      if(phih_mod.lt.(30.+del_phih).and.phih_mod.gt.(30.-del_phih))then
             weight_h=1
          else
             weight_h=0
          endif

c minimum momentum cut
          if(p_h.lt.0.2)then
             weight_h= 0
             del_phih=0.
          endif
c
       Endif
c
        return
        end
c
        subroutine proloss(es)
c
c       simplified model for losses
c
        real es
         if(es.lt.1.0) then
          es=es-0.002+(es-0.93827)/30.865       ! 0 at 1
         endif
        return
        end        

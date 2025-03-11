       	  subroutine bhradgen(ebeam,xff,q2ff,tff,phiff,vv2cut,delta,egamma,thetag,
     +  phig,ichannel,ipolff,ikeygene,stot,probn,probp,probs)
!
!         lab system; OZ is along \vec q, OXZ -- scattering plane
!
!          ebeam -- beam energy
!          xff (or x), q2ff (or q2), tff (ot r) -- kinematic invariants   
!          phiff (phi) -- phi of final proton
!          egamma, thetag, phig -- variables of simulated (second) photon, if ichannel>1
!          ichannel -- channel of scattering: 1-no second photon; 2-second photon along initial electron, 3-second photon along final electron
!          ipolff (or ipol) -- target polarization parameter (0-unpolarized, 1-longitudinally polarized alog \vec q)  
!          iapprff (or iappr) --  1-BH only, 2-BMK
!          helff (or hel) --  hel=1(0) to (not) include polarized part   

!           comment added by Sangbaek Lee, 2021
!           according to I. Akushevich and A. Ilyichev, PHYS. REV. D 98, 013005 (2018),
!           vv2cut: cuts on V. V = missing mass squared in the experimental design (GeV^2)
!           when only charge particles (electron and proton) are detected
!           and/or used for reconstruction of kinematical variables
!           delta: " a minimal energy of the photon we want to generate (i.e., calorimeter resolution)"
!           ikeygene: the number of loops to calculate the rad. probability
!           stot: total cross section.

        implicit none
#include "dvcs.inc"
#include "dvcsmom.inc"
#include "ntupgdvcs.inc"
        integer*4 iacc
        integer*4 nzd,nzdphi,ipolff,iapprff,iaddcontr_p,iaddcontr_s
	real*8 ebeam,s,q2,x,t,phi,sx,xx,mp2,ml,ml2,alpha,barn,
     +	sborn,siborn,sig0,sirad,siradtest,xi,z1m,z2m,z2cur,z1,z2,lll,be,bb,tcol,
     +  tmin,stot,vv2cut,vmax,vv2calc,vv2min_s,vv2min_p
     +  ,sitottest,sitottest_rad,sitottest_non
	real*8 tsp,tpp,zspeak,zppeak,ran2
	real*8 onmz1,onmz2,dz1,dz2,w2,t1lim,t2lim,lay,lat,layz,az,cz,
     +	sz,sqdz,cpr,spr,sbr,cphbr,sphbr,xsp
        real*8 epr,epst,egcut,z1mc,z2mc,sig0s,sig0p,sum1,sum2,vacpol,
     +	egsim,u0,w0,lauw,eta,etaz(3),phis,sum1delta,sum2delta
        real*8 sirad0,sirad000,siradsin,siradsin2,siradcos,siradcos2,
     +	siradcos3,sibor0,sibor000,siborsin,siborsin2,siborcos,siborcos2,
     +  siborcos3,siradadd,llog,sum1phi,sum2phi,sum1test,sum2test
        real*8 sum1add,sum2add,sum1tde,sum2tde,random1,random2
	 real*8 fracacc_s,naccev_s,ntotev_s,fracacc_p,naccev_p,ntotev_p

        parameter(nzd=55)
        parameter(nzdphi=5)
        real*8 deltaz1,deltaz2,ikeymc,delta,fd1(0:nzd),fd2(0:nzd),fd1phi(0:nzdphi),fd2phi(0:nzdphi),
     +	zd1(0:nzd),zd2(0:nzd),zd1phi(0:nzdphi),zd2phi(0:nzdphi)
	real*4 random_num
	integer*4 ipg,in,i,nev,itkey,itpoi,iepoi,ixq2poi,nepoi,nxq2poi,
     +	ntpoi,nphipoi,
     +	ikeyphiint,isamp
	integer*4 ikeydd,ikeyfo,ipol,iappr,ntreg1,ntreg2,ntreg3,
     +	nttot,nphipoi1,nphipoi2,izd
        real*8 v2p,v2m,a1,a2s,a2p,v2s1,v2s2,v2p1,v2p2,dds,ddp,sp  
     
	common/const/alpha,barn,mp2,ml2,ml
        common/kinpoi/s,q2,t,phi,x,xx,sborn,eta(3),ipol
        common/iappr/iappr
        data alpha/0.729735d-2/,barn/0.389379d9/
        ikeyfo=ikeygene

        x=xff
        q2=q2ff
        t=tff
        phi=phiff      
        ipol=ipolff 
        iappr=iapprff
c        hel=helff
        ml=mele!sqrt(ml2)
        ml2=ml**2
        mp2=mp**2

	 s=2d0*mp*ebeam

         lll=log(q2/ml2)
         be=2d0*alpha/pi*(LLL-1d0)
         bb=be/2d0
	 epr=ebeam-q2/2d0/mp/x
	 sx=q2/x
	 xx=s-sx
	 sp=s+xx
	 w2=sx-q2+mp2
  	 lay=sx**2+4d0*mp2*q2
	 t2lim=-0.5d0*((sx-q2)*(sx+sqrt(lay))+2d0*mp2*q2)/w2
	 t1lim=-0.5d0*((sx-q2)*(sx-sqrt(lay))+2d0*mp2*q2)/w2
         tsp=-Q2*xx/(s-q2)
         tpp=-Q2*s/(xx+q2)

 	 lat=t*(t-4d0*mp2)
	 xi=sqrt(lat/lay)

         if(ipol.eq.0)then
           eta(1)=0d0
           eta(2)=0d0
           eta(3)=0d0
         endif
         if(ipol.eq.1)then
           eta(1)=0d0
           eta(2)=0d0
           eta(3)=1d0
         endif
         if(ipol.ge.2)then
           phis=(ipol-2)*pi/2d0
           eta(1)=cos(phis)
           eta(2)=sin(phis)
           eta(3)=0d0
         endif

	 sborn=siborn(s,q2,x,t,cos(phi),sin(phi))
       if (sborn.eq.0) then
         stot = 0D0
         return
       endif
	  vmax=(sqrt(lay*lat)+sx*t)/2d0/mp2-q2+t

      vv2calc=min(vv2cut,vmax)



        a1=4.*mp*cos(phi)*sqrt(q2*(s*xx-mp2*q2))
       v2p=(sx*t+sqrt(lay*lat))/2d0/mp2-q2+t
       v2m=(sx*t-sqrt(lay*lat))/2d0/mp2-q2+t
       a2s=q2*(sp*(sx+2.*t)-lay)-t*(lay+sp*sx)
       a2p=q2*(sp*(sx+2.*t)+lay)+t*(lay-sp*sx)
        dds=a1**2*v2p*v2m+a2s**2
        ddp=a1**2*v2p*v2m+a2p**2
      zspeak=1d0-4.*v2p*v2m*lay/((v2m-v2p)*sqrt(dds)+a2s*(v2p+v2m)+2*v2p*v2m*(lay+sp*sx)) 
      zppeak=1d0-4.*v2p*v2m*lay/((v2m-v2p)*sqrt(ddp)+a2p*(v2p+v2m)+2*v2p*v2m*(lay+sp*sx)) 
        v2s1=v2p*dds/(a1**2*v2p**2+a2s**2)
        v2s2=2d0/((1./v2m+1./v2p)+a2s*(1./v2p-1./v2m)/sqrt(dds))
        v2p1=v2p*ddp/(a1**2*v2p**2+a2p**2)
        v2p2=2d0/((1./v2m+1./v2p)+a2p*(1./v2p-1./v2m)/sqrt(ddp))
        z1m=1d0-2d0*vv2calc*lay/(a1*sqrt((v2p-vv2calc)*(vv2calc-v2m))+a2s+vv2calc*(sp*sx+lay))
        z2m=1d0-2d0*vv2calc*lay/(a1*sqrt((v2p-vv2calc)*(vv2calc-v2m))+a2p+vv2calc*(sp*sx+lay))

       if(vv2calc.le.v2s2.or.cos(phi).ge.0d0)then
	  iaddcontr_s=0
      z1mc=z1m
      else
	  iaddcontr_s=1
      z1mc=zspeak
      endif
	  
      if(vv2calc.le.v2p2.or.cos(phi).ge.0d0)then
	  iaddcontr_p=0
      z2mc=z2m
      else
	  iaddcontr_p=1
      z2mc=zppeak
      endif


         deltaz1=min(delta/ebeam,1d0-z1m)  
         deltaz2=min(delta/(epr+delta),1d0-z2m)

		 
         do izd=0,nzd
           zd1(izd)=z1mc+(1d0-deltaz1-z1mc)*izd/nzd
           zd2(izd)=z2mc+(1d0-deltaz2-z2mc)*izd/nzd
           fd1(izd)=0d0
           fd2(izd)=0d0
         enddo 

         do izd=0,nzdphi
           zd1phi(izd)=z1mc+(z1m-z1mc)*izd/nzdphi
           zd2phi(izd)=z2mc+(z1m-z2mc)*izd/nzdphi
           fd1phi(izd)=0d0
           fd2phi(izd)=0d0
         enddo 

          i=ikeyfo

	  nev=2**i
	  sum1=0d0
      sum2=0d0
	  sum1delta=0d0
      sum2delta=0d0
	  sum1phi=0d0
      sum2phi=0d0
            sborn=siborn(s,q2,x,t,cos(phi),sin(phi)) 
	 do in=1,nev
c 	  z1=z1mc+(1d0-deltaz1-z1mc)*urand(iy)
 	  z1=z1mc+(1d0-1e-16-z1mc)*random_num()
          z2=1.d0
	    layz=(z1*s-xx/z2)**2+4d0*mp2*z1*q2/z2
	    az=-((z1*s-xx/z2)*t+2d0*mp2*(t-z1*q2/z2))/sqrt(lat*layz)
	    cz=(sx*(z1*s-xx/z2)+2d0*(1d0/z2+z1)*mp2*q2)/sqrt(lay*layz)
	    sz=2d0*(1d0/z2-z1)*mp*sqrt(q2*(s*xx-mp2*q2))/sqrt(lay*layz)
	    sqdz=sqrt(cz**2+sz**2*cos(phi)**2-az**2)
	    cpr=(az*cz+sqdz*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    spr=(sqdz*cz-az*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    sbr=sqrt(1d0-(cpr*cz-spr*sz*cos(phi))**2)	
	    cphbr=(cz*spr*cos(phi)+sz*cpr)/sbr
            sphbr=spr*sin(phi)/sbr 
	    xsp=z1*q2/(z1*z2*s-xx)
            etaz(1)=cz*eta(1)+sz*eta(3)  
            etaz(2)=eta(2)
            etaz(3)=-sz*eta(1)+cz*eta(3)
	    xsp=z1*q2/(z1*z2*s-xx)
c	    xb1=z1/(1d0/xb-s*(1d0-z1)/Q2)
	    sig0s=siborn(z1*s,z1*q2,xsp,t,cphbr,sphbr,1,ipol,etaz)
		if(z1.le.1d0-deltaz1)then
	    sum1add=(1d0-z1mc)*(1d0+z1**2)
     +	    *(spr/sqdz*(xsp/x)**2*sig0s   )/(1d0-z1)		
ccc     +	    *(spr/sqdz*(xsp/x)**2*sig0s-sborn)/(1d0-z1)		
     	 sum1=sum1+sum1add/dble(nev)
            do izd=1,nzd
             if(z1.le.zd1(izd))fd1(izd)=fd1(izd)+sum1add/dble(nev)
            enddo
        else
	    sum1add=(1d0-z1mc)*(1d0+z1**2)
     +	    *(spr/sqdz*(xsp/x)**2*sig0s-sborn)/(1d0-z1)		
	    sum1delta=sum1delta+sum1add/dble(nev)		
		endif	

		if(z1.lt.z1m.and.iaddcontr_s.eq.1)then 
	    cpr=(az*cz-sqdz*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    spr=(-sqdz*cz-az*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    sbr=sqrt(1d0-(cpr*cz-spr*sz*cos(phi))**2)	
	    cphbr=(cz*spr*cos(phi)+sz*cpr)/sbr
        sphbr=spr*sin(phi)/sbr 
		sig0s=siborn(z1*s,z1*q2,xsp,t,cphbr,sphbr)
	    sum1add=(1d0-z1mc)*(1d0+z1**2)*(spr/sqdz*(xsp/x)**2*sig0s    )/(1d0-z1)		
	    sum1phi=sum1phi+sum1add/dble(nev) 
            do izd=1,nzdphi
             if(z1.le.zd1phi(izd))fd1phi(izd)=fd1phi(izd)+sum1add/dble(nev)
            enddo
		endif
		
		
          z1=1.d0
c	  z2=z2mc+(1d0-deltaz2-z2mc)*urand(iy)
	  z2=z2mc+(1d0-1e-16-z2mc)*random_num()
	    layz=(z1*s-xx/z2)**2+4d0*mp2*z1*q2/z2
	    az=-((z1*s-xx/z2)*t+2d0*mp2*(t-z1*q2/z2))/sqrt(lat*layz)
	    cz=(sx*(z1*s-xx/z2)+2d0*(1d0/z2+z1)*mp2*q2)/sqrt(lay*layz)
	    sz=2d0*(1d0/z2-z1)*mp*sqrt(q2*(s*xx-mp2*q2))/sqrt(lay*layz)
	    sqdz=sqrt(cz**2+sz**2*cos(phi)**2-az**2)
	    cpr=(az*cz+sqdz*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    spr=(sqdz*cz-az*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    sbr=sqrt(1d0-(cpr*cz-spr*sz*cos(phi))**2)	
	    cphbr=(cz*spr*cos(phi)+sz*cpr)/sbr
            sphbr=spr*sin(phi)/sbr 
            etaz(1)=cz*eta(1)+sz*eta(3)  
            etaz(2)=eta(2)
            etaz(3)=-sz*eta(1)+cz*eta(3)
	    xsp=z1*q2/(z1*z2*s-xx)
c	    xb2=1d0/(1d0/xb-(1d0-z2)*s/Q2)
            sig0p=siborn(s,q2/z2,xsp,t,cphbr,sphbr)
		if(z2.le.1d0-deltaz2)then
            sum2add=(1d0-z2mc)*(1d0+z2**2)*(spr/sqdz
     +	    *(xsp/x)**2*sig0p/z2      )/(1d0-z2)		
ccc     +	    *(xsp/x)**2*sig0p/z2-sborn)/(1d0-z2)		
	    sum2=sum2+sum2add/dble(nev)
            do izd=1,nzd
             if(z2.le.zd2(izd))fd2(izd)=fd2(izd)+sum2add/dble(nev)
            enddo
		else
            sum2add=(1d0-z2mc)*(1d0+z2**2)*(spr/sqdz
     +	    *(xsp/x)**2*sig0p/z2-sborn)/(1d0-z2)		
	    sum2delta=sum2delta+sum2add/dble(nev)
        endif		
		
        if(z2.lt.z2m.and.iaddcontr_p.eq.1)then 
   	      cpr=(az*cz-sqdz*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	      spr=(-sqdz*cz-az*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	      sbr=sqrt(1d0-(cpr*cz-spr*sz*cos(phi))**2)	
	      cphbr=(cz*spr*cos(phi)+sz*cpr)/sbr
          sphbr=spr*sin(phi)/sbr 
          sig0p=siborn(s,q2/z2,xsp,t,cphbr,sphbr)
          sum2add=(1d0-z2mc)*(1d0+z2**2)*(spr/sqdz*(xsp/x)**2*sig0p/z2      )/(1d0-z2)		
   	      sum2phi=sum2phi+sum2add/dble(nev)
            do izd=1,nzdphi
             if(z2.le.zd2phi(izd))fd2phi(izd)=fd2phi(izd)+sum2add/dble(nev)
            enddo         
		 endif

		
	 enddo	


		   

          LLog= log(q2/ml2)-1d0
 		  

	 
          sirad=alpha    /Pi*vacpol(-t)*sborn
     &  	+alpha/2d0/pi*LLog*(sum1+sum2+sum1delta+sum2delta+sum1phi+sum2phi+
     &   (3d0+2.d0*log(deltaz1)+2d0*log(deltaz2)-2d0*deltaz1+0.5d0*deltaz1**2-2d0*deltaz2+0.5d0*deltaz2**2)*sborn)

	 
c          sirad= 
c     &  	+alpha/2d0/pi*(log(q2/ml2)    )
c     &          *(sum1+sum1delta+sum2+sum2delta+(z1m*(z1m+2d0)/2d0+z2m*(z2m+2d0)/2d0
c     &          +2d0*log(1d0-z1m)+2d0*log(1d0-z2m))*sborn)
	 
	 
	 
c	            write(*,'(4f8.3,4g17.8)')sirad/sborn-1.,sum1,sum1delta,sum1+sum1delta,sum2,sum2delta,sum2+sum2delta
	 
          probn=((1d0+alpha/Pi*vacpol(-t))*sborn
     &	  +alpha/2d0/pi*LLog*(sum1delta+sum2delta+(3d0+2.d0*log(deltaz1)+2d0*log(deltaz2)-2d0*deltaz1+0.5d0*deltaz1**2-2d0*deltaz2+0.5d0*deltaz2**2)*sborn))/(sborn+sirad)

          probs=
     &  	alpha/2d0/pi*LLog
     &          *(sum1+sum1phi)/(sborn+sirad)

          probp=
     &  	alpha/2d0/pi*LLog
     &          *(sum2+sum2phi)/(sborn+sirad)

c          write(47,'(g8.1,4f7.3,2g12.4,5f8.2)')delta,x,q2,t,phi,
c     1	  alpha,pi,vacpol(-t),1d2*alpha/Pi*vacpol(-t)
c	 2   ,1d2*(1.d0/(1d0-0.5*alpha/Pi*vacpol(-t))**2-1d0)
	 
	 
c	      write(*,'(6g12.4)')sum1,sum2,sum1delta,sum2delta,sum1phi,sum2phi
	 
	 
C           write(66,'(g8.1,4f7.3,3f8.2,2f10.2,2f8.2)')delta,x,q2,t,phi,
C      1	  1d2*alpha/2d0/pi*LLog*3d0
C      2	 ,1d2*alpha/Pi*vacpol(-t)
C      3   ,1d2*alpha/2d0/pi*LLog*(2.d0*log(deltaz1)+2d0*log(deltaz2)-2d0*deltaz1+0.5d0*deltaz1**2-2d0*deltaz2+0.5d0*deltaz2**2)
C      4   ,1d2*alpha/2d0/pi*LLog*(sum1+sum2+sum1phi+sum2phi)/sborn
C      5   ,1d2*alpha/2d0/pi*LLog*(sum1delta+sum2delta)/sborn
C      6   ,1d2*sirad/sborn
C      6   ,1d2*siradtest/sborn
		  

c          write(51,'(5f7.3,4f8.2,3g12.4)')delta*1d3,x,q2,t,phi
C            write(*,'(5f7.3,3g12.4)')delta*1d3,x,q2,t,phi, probn, probp, probs
c	 6   ,1d2*siradtest/sborn
c     6   ,1d2*(sitottest/sborn-1)
cc	 6   ,1d2*(sitottest_rad/sborn)
cc	 6   ,1d2*(sitottest_non/sborn-1)
c	 7   ,1d2*(exp(alpha/pi*LLog*log(deltaz1*deltaz2))-1d0-alpha/pi*LLog*log(deltaz1*deltaz2))


		  
		  
		  ntotev_s=0.
		  ntotev_p=0.
		  naccev_s=0.
		  naccev_p=0.
          random1=random_num()        
          if(random1.le.probn)then 
             ichannel=1
             egamma=0d0 
             thetag=0d0
             phig=0d0 
          elseif(random1.le.probn+probs)then
             ichannel=2
             ran2=(random1-probn)/probs
			 if(ran2.le.sum1/(sum1+sum1phi))then 
			 random2=ran2*(sum1+sum1phi)   !*sum1/sum1
             do izd=1,nzd
                if(random2.lt.fd1(izd).and.random2.ge.fd1(izd-1))then
                   z1=zd1(izd-1)+(zd1(izd)-zd1(izd-1))*(random2
     +		   -fd1(izd-1))/(fd1(izd)-fd1(izd-1))
                   Egamma=(1.-z1)*Ebeam 
                endif   
             enddo
             else
			 random2=ran2*(sum1+sum1phi)-sum1   !(ran2-sum1/(sum1+sum1phi))*(sum1+sum1phi)/sum1phi*sum1phi
             do izd=1,nzdphi
                if(random2.lt.fd1phi(izd).and.random2.ge.fd1phi(izd-1))then
                   z1=zd1phi(izd-1)+(zd1phi(izd)-zd1phi(izd-1))*(random2
     +		   -fd1phi(izd-1))/(fd1phi(izd)-fd1phi(izd-1))
                   Egamma=(1.-z1)*Ebeam 
                endif   
             enddo		 
			 endif
             phig=0d0
             thetag=acos((s*sx+2*mp2*q2)/sqrt(lay)/s)
		  else
             ichannel=3
             ran2=(random1-probn-probs)/probs
			 if(ran2.le.sum2/(sum2+sum2phi))then 
			 random2=ran2*(sum2+sum2phi)   !*sum2/sum2
             do izd=1,nzd
                if(random2.lt.fd2(izd).and.random2.ge.fd2(izd-1))then
                   z2=zd2(izd-1)+(zd2(izd)-zd2(izd-1))
     +		   *(random2-fd2(izd-1))/(fd2(izd)-fd2(izd-1))
                   Egamma=(1.-z2)/z2*Epr 
                endif   
             enddo
             else
			 random2=ran2*(sum2+sum2phi)-sum2   !(ran2-sum2/(sum2+sum2phi))*(sum2+sum2phi)/sum2phi*sum2phi
             do izd=1,nzdphi
                if(random2.lt.fd2phi(izd).and.random2.ge.fd2phi(izd-1))then
                   z2=zd2phi(izd-1)+(zd2phi(izd)-zd2phi(izd-1))
     +		   *(random2-fd2phi(izd-1))/(fd2phi(izd)-fd2phi(izd-1))
                   Egamma=(1.-z2)/z2*Epr 
                endif   
             enddo		 
			 endif
             phig=0d0
             thetag=acos((xx*sx-2*mp2*q2)/sqrt(lay)/xx)
          endif   

c           write(71,*)ichannel,Egamma,thetag

C           if(ichannel.gt.1)then
C 		  call cutacc(s,xx,q2,phi,t,ichannel,egamma,iacc)
C 		  if(ichannel.eq.2)then
C 		  if(iacc.eq.1)naccev_s=naccev_s+1. 
C 		  ntotev_s=ntotev_s+1.
C 		  endif
C 		  if(ichannel.eq.3)then
C 		  if(iacc.eq.1)naccev_p=naccev_p+1. 
C 		  ntotev_p=ntotev_p+1.
C           endif

C 		  endif
C        print *, ichannel, thetag, egamma, probn, probs, probp
		  




C       fracacc_s=naccev_s/ntotev_s
C       fracacc_p=naccev_p/ntotev_p
C 	stot=(sborn+sirad)*(probn+fracacc_s*probs+fracacc_p*probp)
       stot=sborn+sirad
	
	end
      
	  
	  subroutine cutacc(s,xx,q2,phi,t,ichannel,egamma,iacc)
       implicit none
#include "dvcs.inc"
	   real*8 s,xx,q2,phi,t,alpha,barn,mp2,ml2,ml
	   real*8 missingmass2,complanarity,complanarity0,ptfin,ang_m_c,ptfin2
	   real*8 sx,aly,sqly,cos1,sin1,cos2,sin2,e1,e2,epr,ppr,costp,sintp,scalarpr
	   Integer*4 iacc,i,iacctt
	   	  common/const/alpha,barn,mp2,ml2,ml
	   real*8 k1(0:3),k2(0:3),k(0:3),p1(0:3),p2(0:3),k_3kin(0:3),k_4kin(0:3),kgene(0:3),fin(0:3)

	   sx=s-xx
	   aly=sx**2+4.*mp2*q2
	   sqly=sqrt(aly)

	   cos1=(s*sx+2.*mp2*q2)/(s*sqly)
	   sin1=sqrt(1.-cos1**2)
	   e1=s/(2.*mp)	   
       k1(0)=e1   
       k1(1)=e1*sin1
       k1(2)=0.
       k1(3)=e1*cos1
	   
	   cos2=(xx*sx-2.*mp2*q2)/(xx*sqly)
	   sin2=sqrt(1.-cos2**2)
	   e2=xx/(2.*mp)	   
       k2(0)=e2   
       k2(1)=e2*sin2
       k2(2)=0.
       k2(3)=e2*cos2
       
	   p1(0)=mp
	   p1(1)=0.
	   p1(2)=0.
	   p1(3)=0.

	   epr=mp-t/(2.*mp)
       ppr=sqrt(epr**2-mp2)
	   costp=(2.*mp2*q2-2.*mp2*t-t*sx)/(sqly*sqrt(t*(t-4.*mp2)))
	   sintp=sqrt(1.-costp**2)
	   p2(0)=epr
	   p2(1)=ppr*sintp*cos(phi)
	   p2(2)=ppr*sintp*sin(phi)
	   p2(3)=ppr*costp
	   
	   kgene(0)=egamma
	   if(ichannel.eq.2)then
       kgene(1)=egamma*sin1
       kgene(2)=0.
       kgene(3)=egamma*cos1
	   endif

	   if(ichannel.eq.3)then
       kgene(1)=egamma*sin2
       kgene(2)=0.
       kgene(3)=egamma*cos2
	   endif
	   
	   do i=0,3
	   k_3kin(i)=k1(i)+p1(i)-k2(i)-p2(i)
	   k_4kin(i)=k1(i)+p1(i)-k2(i)-p2(i)-kgene(i)
	   fin(i)=k2(i)+p2(i)+k_4kin(i)
	   enddo
       
c	   write(*,*)' test me2 ',scalarpr(k1,k1),scalarpr(k2,k2),ml2
c	   write(*,*)' test mp2 ',scalarpr(p1,p1),scalarpr(p2,p2),mp2
c	   write(*,*)' test lam ',scalarpr(k_3kin,k_3kin),scalarpr(k_4kin,k_4kin)
c	   write(*,*)' test q2  ',2.*scalarpr(k1,k2),q2
c	   write(*,*)' test s   ',2.*scalarpr(k1,p1),s
c	   write(*,*)' test xx  ',2.*scalarpr(k2,p1),xx
c	   write(*,*)' test t   ',2.*mp2-2*scalarpr(p1,p2),t

       missingmass2=2.*scalarpr(kgene,p2) 
	   complanarity=abs(pi-acos((p2(1)*k_4kin(1)+p2(2)*k_4kin(2))/(sqrt(p2(1)**2+p2(2)**2)*sqrt(k_4kin(1)**2+k_4kin(2)**2))))*180./pi
	   complanarity0=abs(pi-acos((p2(1)*k_3kin(1)+p2(2)*k_3kin(2))/(sqrt(p2(1)**2+p2(2)**2)*sqrt(k_3kin(1)**2+k_3kin(2)**2))))*180./pi
       ptfin2=(fin(1)**2+fin(2)**2+fin(3)**2-((fin(1)*k1(1)+fin(2)*k1(2)+fin(3)*k1(3))/e1)**2)
       ang_m_c=acos((k_3kin(1)*k_4kin(1)+k_3kin(2)*k_4kin(2)+k_3kin(3)*k_4kin(3))/sqrt(k_3kin(1)**2+k_3kin(2)**2+k_3kin(3)**2)/sqrt(k_4kin(1)**2+k_4kin(2)**2+k_4kin(3)**2))*180./pi

	   
	   iacc=1
	   iacctt=10000
	   if(complanarity.gt.1.47313d0)iacctt=iacctt+100 
	   if(ptfin2.gt.0.0921189d0**2)iacctt=iacctt+10 
	   if(ang_m_c.gt.1.00133d0)iacctt=iacctt+1 
	   if(iacctt.gt.10000)iacc=0
	   
c      write(*,'(i4,6f8.4,i7)')ichannel,egamma,missingmass2,complanarity,complanarity0,1000*ptfin2,ang_m_c,iacctt	   

	   
*here are the experimental cuts used in my analysis.
*1) Missing mass squared of the (e,p) system :
*-0.226248 < MM2_ep < 0.198154 GeV2.
*2) Coplanarity, angle between the (virtual photon, proton) and (virtual photon, photon) planes : -1.47313 *< delta_phi < 1.26223 degrees.
*3) Missing transverse momentum of the (e,p,gamma) system : sqrt(p_x*p_x+p_y*p_y) < 0.0921189 GeV.
*4) Angle between predicted and measured photon : theta_gamma_X < 1.00133 degrees.
*Let me know if you need any additional information.
	   
	   end

	   real*8 function scalarpr(aaa,bbb)
	   real*8 aaa(0:3),bbb(0:3)
	   scalarpr=aaa(0)*bbb(0)-aaa(1)*bbb(1)-aaa(2)*bbb(2)-aaa(3)*bbb(3)
	   end
	   
	  
      subroutine dzsub(vv,q2,onmz,ipo,dzfun)
      implicit real*8(a-h,k-m,o-z)
#include "dvcs.inc"
	common/const/alpha,barn,mp2,ml2,ml

      lll=log(q2/ml2)
      be=2d0*alpha/pi*(LLL-1d0)
      z=1d0-onmz


      cc0=(1d0+3d0/8d0*be-be**2/48d0*(lll/3d0+pi**2-47d0/8d0))
      dga=0.5d0*be   !*(1d0-z)**(be/2d0-1d0)
     . *cc0
       if(onmz.gt.0d0)then
      zbzb=onmz**(1d0-be/2d0)
       dga=dga+zbzb*(
     . -be/4d0*(1d0+z)+be**2/32d0*(-4d0*(1d0+z)*log(onmz)
     . -(1d0+3d0*z**2)/onmz*log(z)-5d0-z) )



      endif

       if(1d0-z-4d0*mp*ml/vv.gt.0d0)then
	 zer=max(1d0-z-4d0*mp*ml/vv,0d0)
	 lll1=lll+log(onmz**2)
	 ll53=lll1-5d0/3d0
	 deen=zbzb*alpha**2/pi**2/12d0*zer**(be/2d0)/onmz
     .	  *ll53**2*(1d0+z**2+be/6d0*ll53)
       else
	 deen=0.d0
       endif
      if(ipo.eq.1)then
      dees=zbzb*alpha**2/4d0/pi**2*lll**2*(2d0*(1d0-z**3)/z+.5d0*(1d0-z)
     . +(1d0+z)*log(z))
      elseif(ipo.eq.2)then
      dees=zbzb*alpha**2/4d0/pi**2*lll**2*(5d0*(1d0-z)/2d0
     . +(1d0+z)*log(z))
      endif

      dzfun=dga+deen+dees
***	 dzfun=dzfun*(1d0-z)**(1d0-be/2d0)

      end


	double precision function siborn(s,q2,x,t,cphi,sphi)
        implicit none
#include "dvcs.inc"
	   real*8 s,q2,x,t,cphi,sphi
      real*8 phiphi
      if(sphi.ge.0d0)phiphi=acos(min(1d0,max(-1d0,cphi)))
      if(sphi.lt.0d0)phiphi=2.*pi-acos(min(1d0,max(-1d0,cphi)))
      Ed = s/2./mp
      call bmkxsec(x, Q2, t,  0, pi-phiphi,siborn) ! do not use for pol. target
      siborn = siborn*2.0*pi
      Ed = cl_beam_energy
 	return
	end
      

        double precision function peak(z)
        implicit none
#include "dvcs.inc"
        real*8 z,lay,layz,az,cz,sz,sqdz,cpr,spr,sbr,cphbr,sphbr,xsp,sigp,z1,z2,lat
        real*8 siborn,sigsp,s,q2,t,phi,x,xx,sborn,sx,alpha,barn,mp2,ml2,ml,eta,etaz(3)
        integer*4 isp,ipol  ,isubs
 	common/const/alpha,barn,mp2,ml2,ml
        common/kinpoi/s,q2,t,phi,x,xx,sborn,eta(3),ipol
        common/test/isp,isubs
         if(isp.eq.1)then
          z1=z
          z2=1d0
         endif
         if(isp.eq.2)then
          z1=1d0
          z2=z
         endif
            sx=s-xx 
 	    lat=t*(t-4d0*mp2)
            lay=sx**2+4d0*mp2*q2  
	    layz=(z1*s-xx/z2)**2+4d0*mp2*z1*q2/z2
	    az=-((z1*s-xx/z2)*t+2d0*mp2*(t-z1*q2/z2))/sqrt(lat*layz)
	    cz=(sx*(z1*s-xx/z2)+2d0*(1d0/z2+z1)*mp2*q2)/sqrt(lay*layz)
	    sz=2d0*(1d0/z2-z1)*mp*sqrt(q2*(s*xx-mp2*q2))/sqrt(lay*layz)
	    if(cz**2+sz**2*cos(phi)**2-az**2.le.0d0) then
	    peak=0d0
	    return
	    endif
	    sqdz=sqrt(cz**2+sz**2*cos(phi)**2-az**2)
c	    print*,cz**2+sz**2*cos(phi)**2-az**2,zc
c	    pause
	    cpr=(az*cz+sqdz*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    spr=(sqdz*cz-az*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    sbr=sqrt(1d0-(cpr*cz-spr*sz*cos(phi))**2)	
	    cphbr=(cz*spr*cos(phi)+sz*cpr)/sbr
            sphbr=spr*sin(phi)/sbr 
            etaz(1)=cz*eta(1)+sz*eta(3)  
            etaz(2)=eta(2)
            etaz(3)=-sz*eta(1)+cz*eta(3)
	    xsp=z1*q2/(z1*z2*s-xx)
	    sigsp=siborn(z1*s,z1*q2/z2,xsp,t,cphbr,sphbr)
	    peak=(1d0+z**2)*(spr/sqdz*(xsp/x)**2*sigsp/z2-isubs*sborn)/(1d0-z)	
        end



	

****************** vacpol *************************************

      double precision function vacpol(t)
c contribution from vacuum polarization by leptons (suml) and hadrons (sumh)
      implicit real*8(a-h,l,m,o-z)
#include "dvcs.inc"
	common/const/alpha,barn,mp2,ml2,ml
c      common/cmp/pi,alpha,amp,amp2,aml,aml2,barn
c      include 'const.inc'
      dimension am2(3)
c
c    am2 : squared masses of charge leptons
c
      data am2/.26110d-6,.111637d-1,3.18301d0/

      suml=0.
      do 10 i=1,3
	 a2=2.*am2(i)
	 sqlmi=dsqrt(t*t+2.*a2*t)
	 allmi=dlog((sqlmi+t)/(sqlmi-t))/sqlmi
  10  suml=suml+2.*(t+a2)*allmi/3.-10./9.+4.*a2*(1.-a2*allmi)/3./t
      if(t.lt.1.d0)then
	aaa = -1.345d-9
	bbb = -2.302d-3
	ccc = 4.091
      elseif(t.lt.64d0)then
	aaa = -1.512d-3
	bbb =  -2.822d-3
	ccc = 1.218
      else
	aaa = -1.1344d-3
	bbb = -3.0680d-3
	ccc = 9.9992d-1
      endif
      sumh = -(aaa+bbb*log(1.+ccc*t)) *2*pi/alpha

      vacpol=suml+sumh
c        vacpol=0.
c      print *,t,vacpol,suml,sumh,-(aaa+bbb*log(1.+ccc*t)) *2*pi/alpha,aaa,bbb,ccc,pi,alpha  

      end
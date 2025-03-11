	program dvcs
        implicit none
	real*8 ebeam,s,q2,x,t,phi,sx,xx,mp,mp2,ml,ml2,pi,alpha,barn,sborn,siborn,sig0,sirad,xi,z1m,z2m,z2cur,z1,z2,lll,be,bb,tcol,tmin
	real*8 tsp,tpp,vv2calc,vv2cut,phideg
	real*8 onmz1,onmz2,dz1,dz2,w2,t1lim,t2lim,lay,lat,layz,az,cz,sz,sqdz,cpr,spr,sbr,cphbr,sphbr,xsp
        real*8 epr,epst,egcut,z1mc,z2mc,sig0s,sig0p,sum1,sum2,vacpol,egsim,u0,w0,lauw,eta,etaz(3),phis
        real*8 sirad0,sirad000,siradsin,siradsin2,siradcos,siradcos2,siradcos3,sibor0,sibor000,siborsin,siborsin2,siborcos,siborcos2,siborcos3,siradadd,siexp
        real*8 sum1add,sum2add,Raa1,Rbb1,Rcc1,Rdd1,z1speak,z2speak,zspeak,Raa2,Rbb2,Rcc2,Rdd2,z1ppeak,z2ppeak,zppeak,sxt,alt,aly,vv2min_s,vv2max_s,vv2min_p,vv2max_p,vmax,vv2min_z1m,vv2min_z2m,vv2max_z1m,vv2max_z2m,phi_photon,tplus
	real*4 urand
	integer*4 ipg,in,i,nev,iy,itkey,itpoi,iepoi,ixq2poi,nepoi,nxq2poi,ntpoi,nphipoi,ikeyphiint,isamp,iaddcontr_s,iaddcontr_p
	integer*4 ikeydd,ikeyfo,ipol,ipol1,ipol2,ntreg1,ntreg2,ntreg3,nttot,nphipoi1,nphipoi2,ipoi
        parameter(nepoi=1)
        parameter(nxq2poi=1)
        parameter(ntpoi=1)
c        parameter(ntpoi=4)
c        parameter(nphipoi=35)
        parameter(nphipoi=1)
        real*8 are(nepoi)
        real*8 arx(nxq2poi)
        real*8 arq2(nxq2poi)
        real*8 art(ntpoi)
	common/const/alpha,barn,mp,mp2,ml2,ml,pi
        common/kinpoi/s,q2,t,phi,x,xx,sborn,eta(3),ipol
	data alpha/0.729735d-2/,barn/0.389379d9/,ml2/0.261112d-6/,mp/0.938272d0/,iy/12/
        data are/5.75d0/     ! 4.8, 5.75
        data arx/0.4d0/     
        data arq2/1.8d0/     
C        data arx/0.1d0, 0.3d0, 0.5d0, 0.3d0, 0.5d0 /     
C        data arq2/1.0d0,1.d0,1.d0,2.d0,2.d0/     
c       data art/-0.2d0/     
cc        data art/-0.5d0,-1.2d0,-1.5d0,-1.8d0/     

c        itkey =154015 ! values from art (1) or between tmin/tmax (0, ntpoi>1)
c         itkey =81008 ! values from art (1) or between tmin/tmax (0, ntpoi>1)
c         itkey =40506 ! values from art (1) or between tmin/tmax (0, ntpoi>1)
        itkey =0 ! values from art (1) or between tmin/tmax (0, ntpoi>1)
        ikeydd=0 ! calculate general eqs with DD (1) or not (0)
        ikeyfo=20 ! calculate first-order eqs (1) or not (0)
        ipol1=0
        ipol2=0
        ikeyphiint=0 ! integration over phi nphipoi must be equal 1
          if(ikeyphiint.eq.1.and.nphipoi.ne.1) stop 'ikeyphiint'
c      	 egcut=0.3
      	 egcut=999999999999.

        ml=sqrt(ml2)
        mp2=mp**2
	pi=atan(1d0)*4d0
        do iepoi=1,nepoi
         ebeam=are(iepoi)
	 s=2d0*mp*ebeam
c        do ixq2poi=1,nxq2poi
c	 x=arx(ixq2poi)
c	 q2=arq2(ixq2poi)

c           open(8,file='kinematics3_0.dat')
c           open(8,file='linked.dat')
           open(8,file='xykin.dat')
           do ipoi=1,1000  

c           read(8,*,end=111)x,q2,tplus,phi_photon		   
c           phi=pi-phi_photon*pi/180d0
 
           read(8,*,end=111)x,q2,tplus,phideg
           phi=phideg*pi/180.
		   
           t=-tplus
		   

         lll=log(q2/ml2)
         be=2d0*alpha/pi*(LLL-1d0)
         bb=be/2d0
	 epr=ebeam-q2/2d0/mp/x
	 sx=q2/x
	 xx=s-sx
	 w2=sx-q2+mp2
  	 lay=sx**2+4d0*mp2*q2
	 aly=lay
	 t2lim=-0.5d0*((sx-q2)*(sx+sqrt(lay))+2d0*mp2*q2)/w2
	 t1lim=-0.5d0*((sx-q2)*(sx-sqrt(lay))+2d0*mp2*q2)/w2
         tsp=-Q2*xx/(s-q2)
         tpp=-Q2*s/(xx+q2)
		 
c		 write(*,*)tsp,tpp
		 
c        if(itkey.gt.10000)then
c            ntreg1=itkey/10000
c            ntreg2=(itkey-ntreg1*10000)/100
c            ntreg3=itkey-ntreg1*10000-ntreg2*100
c            nttot=ntreg1+ntreg2+ntreg3            
c        else
c            nttot=ntpoi           
c        endif
cc        print *,nttot,itkey,ntreg1,ntreg2,ntreg3,t2lim,tpp,tsp,t1lim
cc        do itpoi=1,nttot
cc        do itpoi=3,3
c        do itpoi=nttot,1,-1
c         if(itkey.eq.1)then
c            t=art(itpoi)
c         elseif(itkey.eq.0)then
c            epst=(t1lim-t2lim)/100.
c            if(ntpoi.ne.1)t=t2lim+epst+(t1lim-t2lim-2.*epst)*(itpoi-1)/(ntpoi-1)
c            if(ntpoi.eq.1)t=(t2lim+t1lim)/2
c         elseif(itkey.gt.10000)then
c            if(itpoi.le.ntreg1)then
c               epst=(tpp-t2lim)/100000.
c               t=t2lim+epst+(tpp-t2lim-2.*epst)*(itpoi-1)/(ntreg1-1)  
c            elseif(itpoi.le.ntreg1+ntreg2)then
c               epst=(tsp-tpp)/100000.
c               t=tpp+epst+(tsp-tpp-2.*epst)*(itpoi-ntreg1-1)/(ntreg2-1)  
c            else
c               epst=(t1lim-tsp)/100000.
c               t=tsp+epst+(t1lim-tsp-2.*epst)*(itpoi-ntreg1-ntreg2-1)/(ntreg3-1)  
c            endif   
c         endif
		 
cc		                 t=-2.774
		 
		sxt=sx+t 
 	 lat=t*(t-4d0*mp2)
	 alt=lat
	 xi=sqrt(lat/lay)
c	do ipg=1,nphipoi
cc	do ipg=5,5
c         nphipoi1=nphipoi/3
c         nphipoi2=nphipoi-2*nphipoi1
c         if(ipg.le.nphipoi1)then
c           phi=dble(ipg-1)*0.9*pi/nphipoi1             
c         endif  
c         if(ipg.gt.nphipoi1.and.ipg.le.nphipoi1+nphipoi2)then
c           phi=0.9*pi+dble(ipg-1-nphipoi1)*(1.1-0.9)*pi/(nphipoi2-1)
c         endif  
c         if(ipg.gt.nphipoi1+nphipoi2)then
c           phi=1.1*pi+dble(ipg-nphipoi1-nphipoi2)*(2.0-1.1)*pi/nphipoi1             
c         endif  
cc	 phi=1.5/180.*pi+dble(ipg-1)*(2d0*pi-2*1.5/180.*pi)/(nphipoi-1)
cc      phi=0d0
cc      phi=pi/2d0
c                       phi=dble(ipg-1)*2d0*pi/(nphipoi-1) 
cc					   					   phi=2.793
        do ipol=ipol1,ipol2
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
	 sborn=siborn(s,q2,x,t,cos(phi),sin(phi),1,ipol,eta)
c	 print*,s,xx,sx,q2,x,t,phi,sborn
	   z1m=(t*xx-2d0*mp2*t+xi*(xx*sx-2d0*mp2*q2))/(t*s-2d0*mp2*q2+xi*(s*sx+2d0*mp2*q2))	
	   z2m=(t*xx+2d0*mp2*q2+xi*(xx*sx-2d0*mp2*q2))/(t*s+2d0*mp2*t+xi*(s*sx+2d0*mp2*q2))
	 
	 
	 Raa1 = 4.*mp2*(aly*(mp2*Q2**2+t*S*(S-Q2))+alt*Q2*sin(phi)**2*(S*xx-mp2*Q2)) 
	 Rbb1 = 4.*mp2*(t*aly*(2*mp2*Q2-xx*(2*S-Q2)-t*(Q2-S))-2.*alt*Q2*sin(phi)**2*(S*xx-mp2*Q2)) 
	 Rcc1 = 4.*mp2*(t*aly*(xx**2+t*(-xx+mp2))+alt*Q2*sin(phi)**2*(S*xx-mp2*Q2)) 
	 Rdd1 = 16.*Mp2**2*aly*alt*(aly*(t*Q2-xx*Q2-t*S)**2-4.*Q2*sin(phi)**2*(S*xx-mp2*Q2)*((t+Q2)**2*mp2+t*Sxt*(Sx-Q2)))
	 Rcc2 = 4.*mp2*(aly*(mp2*Q2**2+t*xx*(Q2+xx))+alt*Q2*sin(phi)**2*(S*xx-mp2*Q2)) 
	 Rbb2 = 4.*mp2*(t*aly*(2*mp2*Q2-S*(2*xx+Q2)-t*(Q2+xx))-2.*alt*Q2*sin(phi)**2*(S*xx-mp2*Q2)) 
	 Raa2 = 4.*mp2*(t*aly*(S**2+t*(S+mp2))+alt*Q2*sin(phi)**2*(S*xx-mp2*Q2)) 
	 Rdd2 = 16.*Mp2**2*aly*alt*(aly*(t*Q2+S*Q2+t*xx)**2-4.*Q2*sin(phi)**2*(S*xx-mp2*Q2)*((t+Q2)**2*mp2+t*Sxt*(Sx-Q2)))
	 
	   z1speak=(-Rbb1+sqrt(Rdd1))/2./Raa1
	   z2speak=(-Rbb1-sqrt(Rdd1))/2./Raa1
	   z1ppeak=(-Rbb2+sqrt(Rdd2))/2./Raa2
	   z2ppeak=(-Rbb2-sqrt(Rdd2))/2./Raa2
	 
	 zspeak=max(z1speak,z2speak)
	 zppeak=max(z1ppeak,z2ppeak)

c	  write(*,*)z1speak,z2speak,raa1,rbb1,rcc1,rdd1
c	  write(*,*)z1ppeak,z2ppeak,raa2,rbb2,rcc2,rdd2
c	  stop
	 
c	 call vv2fromz(1d0,1d0,vv2min_s,vv2max_s)
c      stop
	 call vv2fromz(zspeak,1d0,vv2min_s,vv2max_s)
	 call vv2fromz(z1m,1d0,vv2min_z1m,vv2max_z1m)
	 call vv2fromz(1d0,zppeak,vv2min_p,vv2max_p)
	 call vv2fromz(1d0,z2m,vv2min_z2m,vv2max_z2m)
	 
	 vv2cut=0.3d0
c 	 vv2cut=999999999999990.3d0
	  vmax=(sqrt(aly*alt)+sx*t)/2d0/mp2-q2+t

      vv2calc=min(vv2cut,vmax)
	 
c                write(*,*)' vmax=',vmax,vv2calc
	 
c       write(*,*)' z1m,z2m - a ',z1m,z2m
	  call z1z2vmax(z1m,z2m,vv2calc)
c       write(*,*)' z1m,z2m - b ',z1m,z2m
	 
	 
	 
c		z1m=0.99d0
c		z2m=0.99d0
	 
	 
c	 write(*,'(4hvv1 ,6g13.5)')vv2min_z1m,vv2min_s,vv2max_s,vv2max_z1m,vmax
c	 write(*,'(4hvv2 ,6g13.5)')vv2min_z2m,vv2min_p,vv2max_p,vv2max_z2m,vmax
c	 write(*,'(4hzz12,6g13.5)')z1m,z1speak,z2speak,z2m,z1ppeak,z2ppeak
c	 write(*,'(4hzz12,6g13.5)')z1m,zspeak,z2m,zppeak
	 
c	 write(*,*)S,x,q2,t,phi,sin(phi)**2*(S*X-mp2*Q2)
c	 write(*,*)z1m,z1speak,z2speak,Raa1,Rbb1,Rcc1,Rdd1
	 
c	 z1m=min(z1m,zspeak)
c	 z2m=min(z2m,zppeak)
	 
c	 write(*,*)z1m,z1speak,z2speak
c	 write(*,*)z1m,z2m,zspeak,zppeak
c	 stop
	 
******************************************************
c	             lauw=4d0*w2*(s*xx*q2-mp2*q2**2-ml2*lay)*(t1lim-t)*(t-t2lim)
c	    w0=-0.5d0*(q2+t)+0.5d0*(s+xx)/lay*(sx*(q2-t)+2d0*t*q2)+sqrt(lauw)/lay*cos(phi)	
c	    u0=w0+q2+t
c
c	 		   write(*,*)' vamx_1=',w0*(1d0-z1m)
c	 		   write(*,*)' vamx_2=',u0*(1d0/z2m-1d0)
c	 		   write(*,*)' z1m=',z1m
c	 		   write(*,*)' z2m=',z2m
******************************************************
	 
        if(ikeydd.ne.0)then
         do isamp=1,10
	 do i=ikeydd,ikeydd
	  nev=2**i
	  sirad=0d0
          if(ikeyphiint.eq.1)then
	    sirad0=0d0
	    sirad000=0d0
	    siradsin=0d0
	    siradsin2=0d0
	    siradcos=0d0
	    siradcos2=0d0
	    siradcos3=0d0
	    sibor0=0d0
	    sibor000=0d0
	    siborsin=0d0
	    siborsin2=0d0
	    siborcos=0d0
	    siborcos2=0d0
	    siborcos3=0d0
          else
            sborn=siborn(s,q2,x,t,cos(phi),sin(phi),1,ipol,eta) 
          endif
	 do in=1,nev
          if(ikeyphiint.eq.1)then
            phi=urand(iy)*2d0*pi
            lauw=4d0*w2*(s*xx*q2-mp2*q2**2-ml2*lay)*(t1lim-t)*(t-t2lim)
	    w0=-0.5d0*(q2+t)+0.5d0*(s+xx)/lay*(sx*(q2-t)+2d0*t*q2)+sqrt(lauw)/lay*cos(phi)	
	    u0=w0+q2+t
            sborn=siborn(s,q2,x,t,cos(phi),sin(phi),1,ipol,eta) 
          endif
	  onmz1=(1d0-z1m)*urand(iy)**(1d0/bb)
	  onmz2=(1d0-z2m)*urand(iy)**(1d0/bb)
	  z1=1.d0-onmz1
	  z2=1.d0-onmz2
	  z2cur=(t*xx+2d0*z1*mp2*q2+xi*(xx*sx-2d0*mp2*q2))/(z1*t*s+2d0*mp2*t+z1*xi*(s*sx+2d0*mp2*q2))
	  egsim=ebeam*(1d0-z1)+epr*(1d0/z2-1d0)
c           write(*,*)z1,z2,egsim,egcut
          if(z2.lt.z2cur.or.egsim.gt.egcut)then
c           write(*,*)z1,z2,egsim,egcut
c          if(z2.lt.z2cur)then
	    sig0=0d0
	  else
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
	    sig0=siborn(z1*s,z1*q2/z2,xsp,t,cphbr,sphbr,2,ipol,etaz)
	    call dzsub(s,q2,onmz1,1,dz1)
	    call dzsub(s,q2,onmz2,2,dz2)
	  endif	
          siradadd=(xsp/x)**2*spr/sqdz*dz1*dz2*sig0/z2
	  sirad=sirad+siradadd
          if(ikeyphiint.eq.1)then
	    sirad0=sirad0+siradadd
	    sirad000=sirad000+siradadd*u0*w0
	    siradsin=siradsin+siradadd*u0*w0*sin(phi)
	    siradsin2=siradsin2+siradadd*u0*w0*sin(2d0*phi)
	    siradcos=siradcos+siradadd*u0*w0*cos(phi)
	    siradcos2=siradcos2+siradadd*u0*w0*cos(2d0*phi)
	    siradcos3=siradcos3+siradadd*u0*w0*cos(3d0*phi)
	    sibor0=sibor0+sborn
	    sibor000=sibor000+sborn*u0*w0
	    siborsin=siborsin+sborn*u0*w0*sin(phi)
	    siborsin2=siborsin2+sborn*u0*w0*sin(2d0*phi)
	    siborcos=siborcos+sborn*u0*w0*cos(phi)
	    siborcos2=siborcos2+sborn*u0*w0*cos(2d0*phi)
	    siborcos3=siborcos3+sborn*u0*w0*cos(3d0*phi)
          endif
	 enddo	
	  sirad=sirad*(1d0-z1m)**bb*(1d0-z2m)**bb/bb**2/dble(nev)
          if(ikeyphiint.eq.1)then
	    sirad0=sirad0*(1d0-z1m)**bb*(1d0-z2m)**bb/bb**2/dble(nev)
	    sirad000=sirad000*(1d0-z1m)**bb*(1d0-z2m)**bb/bb**2/dble(nev)
	    siradsin=siradsin*(1d0-z1m)**bb*(1d0-z2m)**bb/bb**2/dble(nev)
	    siradsin2=siradsin2*(1d0-z1m)**bb*(1d0-z2m)**bb/bb**2/dble(nev)
	    siradcos=siradcos*(1d0-z1m)**bb*(1d0-z2m)**bb/bb**2/dble(nev)
	    siradcos2=siradcos2*(1d0-z1m)**bb*(1d0-z2m)**bb/bb**2/dble(nev)
	    siradcos3=siradcos3*(1d0-z1m)**bb*(1d0-z2m)**bb/bb**2/dble(nev)
	    sibor0=sibor0/dble(nev)
	    sibor000=sibor000/dble(nev)
	    siborsin=siborsin/dble(nev)
	    siborsin2=siborsin2/dble(nev)
	    siborcos=siborcos/dble(nev)
	    siborcos2=siborcos2/dble(nev)
	    siborcos3=siborcos3/dble(nev)
c  	    write(*,'(i5,14g12.3)')ipol,ebeam,q2,x,t,sibor0,sibor000,siborsin,siborsin2,siborcos,siborcos2,siborcos3
c  	    write(*,'(i5,14g12.3)')ipol,ebeam,q2,x,t,sirad0,sirad000,siradsin,siradsin2,siradcos,siradcos2,siradcos3
c  	    write(*,'(i5,14g12.3)')ipol,ebeam,q2,x,t,sibor0,sibor000,siborsin,siborsin2,siborcos,siborcos2,siborcos3
            tmin=-mp2*x**2/(1.-x+x*mp2/Q2)
            tcol=-Q2*(Q2-x*s)/(x*(q2-s))
  	    write(*,'(i5,3f8.3,f16.12,4g14.5)')ipol,ebeam,q2,x,t,sirad0,sibor0,100*(sirad0/sibor0-1d0)
  	    write(43,'(2i5,14g12.3)')ipol,itpoi,t,sibor0,sibor000,siborsin,siborsin2,siborcos,siborcos2,siborcos3
  	    write(44,'(2i5,14g12.3)')ipol,itpoi,t,sirad0,sirad000,siradsin,siradsin2,siradcos,siradcos2,siradcos3
  	    write(60+ixq2poi,'(i5,4f8.3,f16.12,4g14.6)')ipol,ebeam,egcut,q2,x,t,sirad0,sibor0,100*(sirad0/sibor0-1d0)
          endif
          if(ikeyphiint.eq.0)then
c	    write(*,'(i10,2f8.3,2g11.3,2pf6.3)')nev,t,phi,sirad,sborn,sirad/sborn-1d0
	    write(*,'(i5,4g12.3,4g15.6)')ipol,ebeam,q2,x,t,sirad,sborn,100*(sirad/sborn-1d0)
 	    write(20+itpoi,'(4g15.6)')phi,sirad,sborn,100*(sirad/sborn-1d0)
          endif
	 enddo	
         enddo
        endif 

        if(ikeyfo.ne.0)then
         do isamp=1,1
c	 z1mc=max(1d0-egcut/ebeam,z1m)
c	 z2mc=max(epr/(egcut+epr),z2m)
c	 z1mc=max(1d0-egcut/ebeam,zspeak)
c	 z2mc=max(epr/(egcut+epr),zppeak)

c      z1mc=min(z1m,zspeak)
c      z2mc=min(z2m,zppeak)
	 
      if(vv2calc.lt.vv2min_s.or.cos(phi).gt.0d0)then
	  iaddcontr_s=0
      z1mc=z1m
      else
	  iaddcontr_s=1
      z1mc=zspeak
c      z1mc=z1m
      endif
	  
      if(vv2calc.lt.vv2min_p.or.cos(phi).gt.0d0)then
	  iaddcontr_p=0
      z2mc=z2m
      else
	  iaddcontr_p=1
      z2mc=zppeak
c      z2mc=z2m
      endif

	   
	 
	 
c	 write(*,*)z1mc,z1m,zspeak,vv2calc,vv2min_s,vv2max_s
c	 write(*,*)z2mc,z2m,zppeak,vv2calc,vv2min_p,vv2max_p
	 
	 do i=ikeyfo,ikeyfo
	  nev=2**i
	  sum1=0d0
          sum2=0d0
          if(ikeyphiint.eq.1)then
	    sirad0=0d0
	    sirad000=0d0
	    siradsin=0d0
	    siradsin2=0d0
	    siradcos=0d0
	    siradcos2=0d0
	    siradcos3=0d0
	    sibor0=0d0
	    sibor000=0d0
	    siborsin=0d0
	    siborsin2=0d0
	    siborcos=0d0
	    siborcos2=0d0
	    siborcos3=0d0
          else 
            sborn=siborn(s,q2,x,t,cos(phi),sin(phi),1,ipol,eta) 
c			write(*,*)' sborn= ',sborn
          endif
	 do in=1,nev
          if(ikeyphiint.eq.1)then
            phi=urand(iy)*2d0*pi
            lauw=4d0*q2*w2*(s*xx-mp2*q2-ml2*lay)*(t1lim-t)*(t-t2lim)
	    w0=-0.5d0*(q2+t)+0.5d0*(s+xx)/lay*(sx*(q2-t)+2d0*t*q2)+sqrt(lauw)/lay*cos(phi)	
	    u0=w0+q2+t
            sborn=siborn(s,q2,x,t,cos(phi),sin(phi),1,ipol,eta) 
          endif
 	  z1=z1mc+(1d0-1e-16-z1mc)*urand(iy)
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
	    sum1add=(1d0-z1mc)*(1d0+z1**2)*(spr/sqdz*(xsp/x)**2*sig0s-sborn)/(1d0-z1)		
	    sum1=sum1+sum1add/dble(nev)
c	    sig0s=siborn(z1*s,z1*q2,xsp,t,cos(phi),1,ipol)
c	  sum1=sum1+(1d0-z1m)/dble(nev)*(1d0+z1**2)*(spr/sqdz*(xsp/x)**2*sig0s-sborn)/(1d0-z1)		

        if(z1.lt.z1m.and.iaddcontr_s.eq.1)then 
	    cpr=(az*cz-sqdz*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    spr=(-sqdz*cz-az*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    sbr=sqrt(1d0-(cpr*cz-spr*sz*cos(phi))**2)	
	    cphbr=(cz*spr*cos(phi)+sz*cpr)/sbr
        sphbr=spr*sin(phi)/sbr 
		sig0s=siborn(z1*s,z1*q2,xsp,t,cphbr,sphbr,1,ipol,etaz)
	    sum1add=(1d0-z1mc)*(1d0+z1**2)*(-spr/sqdz*(xsp/x)**2*sig0s-sborn)/(1d0-z1)		
	    sum1=sum1-sum1add/dble(nev)
         
		 endif
			


          z1=1.d0
	  z2=z2mc+(1d0-1e-16-z2mc)*urand(iy)
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
            sig0p=siborn(s,q2/z2,xsp,t,cphbr,sphbr,1,ipol,etaz)
            sum2add=(1d0-z2mc)*(1d0+z2**2)*(spr/sqdz*(xsp/x)**2*sig0p/z2-sborn)/(1d0-z2)		
	    sum2=sum2+sum2add/dble(nev)

        if(z2.lt.z2m.and.iaddcontr_p.eq.1)then 
	    cpr=(az*cz-sqdz*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    spr=(-sqdz*cz-az*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    sbr=sqrt(1d0-(cpr*cz-spr*sz*cos(phi))**2)	
	    cphbr=(cz*spr*cos(phi)+sz*cpr)/sbr
            sphbr=spr*sin(phi)/sbr 

            sig0p=siborn(s,q2/z2,xsp,t,cphbr,sphbr,1,ipol,etaz)
            sum2add=(1d0-z2mc)*(1d0+z2**2)*(-spr/sqdz*(xsp/x)**2*sig0p/z2-sborn)/(1d0-z2)		
	    sum2=sum2-sum2add/dble(nev)
         
		 endif



		siradadd=alpha    /Pi*vacpol(-t)*sborn
     &  	+alpha/2d0/pi*(log(q2/ml2)-1d0)
c     &          *(sum1add+sum2add+(z1mc*(z1mc+2d0)/2d0+z2mc*(z2mc+2d0)/2d0
c     &          +2d0*log(1d0-z1mc)+2d0*log(1d0-z2mc))*sborn)
     &          *(sum1add+sum2add+(z1m*(z1m+2d0)/2d0+z2m*(z2m+2d0)/2d0
     &          +2d0*log(1d0-z1m)+2d0*log(1d0-z2m))*sborn)
          if(ikeyphiint.eq.1)then
	    sirad0=sirad0+siradadd
	    sirad000=sirad000+siradadd*u0*w0
	    siradsin=siradsin+siradadd*u0*w0*sin(phi)
	    siradsin2=siradsin2+siradadd*u0*w0*sin(2d0*phi)
	    siradcos=siradcos+siradadd*u0*w0*cos(phi)
	    siradcos2=siradcos2+siradadd*u0*w0*cos(2d0*phi)
	    siradcos3=siradcos3+siradadd*u0*w0*cos(3d0*phi)
	    sibor0=sibor0+sborn
	    sibor000=sibor000+sborn*u0*w0
	    siborsin=siborsin+sborn*u0*w0*sin(phi)
	    siborsin2=siborsin2+sborn*u0*w0*sin(2d0*phi)
	    siborcos=siborcos+sborn*u0*w0*cos(phi)
	    siborcos2=siborcos2+sborn*u0*w0*cos(2d0*phi)
	    siborcos3=siborcos3+sborn*u0*w0*cos(3d0*phi)
          endif         
	 enddo	
caku      sirad=alpha/2d0/Pi*vacpol(-t)*sborn
          if(ikeyphiint.ne.1)then
          sirad= alpha    /Pi*vacpol(-t)*sborn
c     &  	+alpha/2d0/pi*(log(q2/ml2)-1d0)
     &  	+alpha/2d0/pi*(log(q2/ml2)    )
     &          *(sum1+sum2+(z1m*(z1m+2d0)/2d0+z2m*(z2m+2d0)/2d0
     &          +2d0*log(1d0-z1m)+2d0*log(1d0-z2m))*sborn)

          siexp=  exp(alpha/pi*log(q2/ml2)*log((1d0-z1m)*(1d0-z2m)))*sborn   + alpha    /Pi*vacpol(-t)*sborn
c     &  	+alpha/2d0/pi*(log(q2/ml2)-1d0)
     &  	+alpha/2d0/pi*(log(q2/ml2)    )
     &          *(sum1+sum2+(z1m*(z1m+2d0)/2d0+z2m*(z2m+2d0)/2d0
     &          )*sborn)

	 
c           write(*,'(4f8.3,4g11.3)')x,q2,t,phi,sborn,sirad,sirad/sborn*1d2
           write(61,'(4f8.3,f14.3,6f11.3)')x,q2,t,phi,sborn,sirad/sborn*1d2
           write(*,'(3f8.3,f8.1,10g12.4)')x,q2,t,phi*180/pi,sborn,sirad/sborn*1d2,1d0+sirad/sborn,siexp/sborn

		   call stest(sum1,sum2,z1m,z2m,zspeak,zppeak,iaddcontr_s,iaddcontr_p,sborn)
          sirad= alpha    /Pi*vacpol(-t)*sborn
     &  	+alpha/2d0/pi*(log(q2/ml2)    )
     &          *(sum1+sum2+(z1m*(z1m+2d0)/2d0+z2m*(z2m+2d0)/2d0
     &          +2d0*log(1d0-z1m)+2d0*log(1d0-z2m))*sborn)
c           write(*,'(3f8.3,f8.1,10g12.4)')x,q2,t,phi*180/pi,sborn,sirad/sborn*1d2,1d0+sirad/sborn
	 

	 
	 
c	  write(*,'(i10,2f8.3,5g11.3,f6.3)')nev,t,phi,z1mc,z2mc,sum1,sum2,sirad
c	  write(*,'(i10,2f8.3,2g11.3,2pf6.3)')nev,t,phi,sirad,sborn,sirad/sborn
ccc	  write(*,'(4g13.5)')phi,sirad,sborn,100*sirad/sborn
c           write(*,'(4f8.3,4g11.3)')x,q2,t,phi,sborn,sirad,sirad/sborn*1d2
c	  write(30+itpoi,'(4g13.5)')phi,sirad,sborn,100*sirad/sborn
c	  stop
          endif
          if(ikeyphiint.eq.1)then
	    sirad0=sirad0/dble(nev)
	    sirad000=sirad000/dble(nev)
	    siradsin=siradsin/dble(nev)
	    siradsin2=siradsin2/dble(nev)
	    siradcos=siradcos/dble(nev)
	    siradcos2=siradcos2/dble(nev)
	    siradcos3=siradcos3/dble(nev)
	    sibor0=sibor0/dble(nev)
	    sibor000=sibor000/dble(nev)
	    siborsin=siborsin/dble(nev)
	    siborsin2=siborsin2/dble(nev)
	    siborcos=siborcos/dble(nev)
	    siborcos2=siborcos2/dble(nev)
	    siborcos3=siborcos3/dble(nev)
c  	    write(*,'(i5,14g12.3)')ipol,ebeam,q2,x,t,sibor0,sibor000,siborsin,siborsin2,siborcos,siborcos2,siborcos3
c  	    write(*,'(i5,14g12.3)')ipol,ebeam,q2,x,t,sirad0,sirad000,siradsin,siradsin2,siradcos,siradcos2,siradcos3
c  	    write(*,'(i5,14g12.3)')ipol,ebeam,q2,x,t,sibor0,sibor000,siborsin,siborsin2,siborcos,siborcos2,siborcos3
ccc  	    write(*,'(i5,3f8.3,f16.12,4g14.5)')ipol,ebeam,q2,x,t,sirad0,sibor0,100*(sirad0/sibor0)


  	    write(53,'(2i5,14g12.3)')ipol,itpoi,t,sibor0,sibor000,siborsin,siborsin2,siborcos,siborcos2,siborcos3
  	    write(54,'(2i5,14g12.3)')ipol,itpoi,t,sirad0,sirad000,siradsin,siradsin2,siradcos,siradcos2,siradcos3
  	    write(70+ixq2poi,'(i5,4f8.3,f16.12,4g14.6)')ipol,ebeam,egcut,q2,x,t,sirad0,sibor0,100*(sirad0/sibor0)
          endif

	enddo	
			
	enddo	
        endif 

c	enddo	
c	enddo	
	enddo	

	
	enddo	

111        continue	
	enddo	
	
	end

	
	        subroutine z1z2vmax(z1m,z2m,vv2)
c        implicit real*8 (a-h,o-z)
        implicit none
        real*8 dlamuw,u,w,vv2,z1m,z2m,alpha,barn,mp,mp2,ml2,ml,pi,s,q2,tm,phi,xs,xx,sborn,eta
		real*8 x,t,sx,sxp,amp2,aml2,aly,sxt
		integer*4 ipol
    	common/const/alpha,barn,mp,mp2,ml2,ml,pi
        common/kinpoi/s,q2,tm,phi,xs,xx,sborn,eta(3),ipol

		x=xx 
	    t=-tm
		sx=s-xx  
		sxt=sx-t
	     amp2=mp2
	     aml2=ml2
		 sxp=s+xx
         aly=sx**2+4d0*amp2*q2
		 
        dlamuw=4d0*(S*x*q2-mp2*q2**2-ml2*aly)*((sxt*(sx-q2)-sx*vv2)*t-amp2*((q2-t-vv2)**2+4d0*q2*vv2))
        w=0.5d0*(t+vv2-q2)+sxp/2d0/aly*(sx*(q2+t+vv2)-2d0*t*q2)+sqrt(max(0d0,dlamuw))/aly*cos(phi)
        u=-0.5d0*(t+vv2-q2)+sxp/2d0/aly*(sx*(q2+t+vv2)-2d0*t*q2)+sqrt(max(0d0,dlamuw))/aly*cos(phi)

		z1m=1-vv2/w
		z2m=u/(u+vv2)

c         z1m=0.99d0
c         z2m=0.99d0
		
c		write(*,*)vv2,w,dlamuw,(S*x*q2-amp2*q2**2-aml2*aly),s,x,q2,amp2
c		stop
       
c        write(*,*)' 22',w,u,vv2,4*dlamuw,0.5d0*(t+vv2-q2)
*			if(dlamuw.lt.0d0)write(62,*)vv2,vmax,dlamuw
c            lauw=4d0*q2*w2*(s*xx-mp2*q2-ml2*lay)*(t1lim-t)*(t-t2lim)
c            w0=-0.5d0*(q2+t)+0.5d0*(s+xx)/lay*(sx*(q2-t)+2d0*t*q2)+sqrt(lauw)/lay*cos(phi)	
c	    u0=w0+q2+t
        
        end
   
	
      subroutine vv2fromz(z1,z2,vv2min,vv2max)
c	        implicit real*8(a-h,k-m,o-z)
	        implicit none
			real*8 z1,z2,vv2min,vv2max,alpha,barn,mp,mp2,ml2,ml,pi,s,q2,t,phi,x,xx,sborn,eta
			real*8 sx,lay,lat,layz,az,cz,sz,sqdz,cpr1,cpr2,vv2_1,vv2_2,sqly,sqlt,dz0
			Integer*4 ipol
	common/const/alpha,barn,mp,mp2,ml2,ml,pi
        common/kinpoi/s,q2,t,phi,x,xx,sborn,eta(3),ipol
		
		 sx=s-xx
  	     lay=sx**2+4d0*mp2*q2
    	 lat=t*(t-4d0*mp2)
		 sqlt=sqrt(lat)
		 sqly=sqrt(lay)
		layz=(z1*s-xx/z2)**2+4d0*mp2*z1*q2/z2
	    az=-((z1*s-xx/z2)*t+2d0*mp2*(t-z1*q2/z2))/sqrt(lat*layz)
	    cz=(sx*(z1*s-xx/z2)+2d0*(1d0/z2+z1)*mp2*q2)/sqrt(lay*layz)
	    sz=2d0*(1d0/z2-z1)*mp*sqrt(q2*(s*xx-mp2*q2))/sqrt(lay*layz)
	    dz0=cz**2+sz**2*cos(phi)**2-az**2
		if(dz0.lt.-1d-4) stop
		sqdz=sqrt(max(0d0,dz0))
	    cpr1=(az*cz+sqdz*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    cpr2=(az*cz-sqdz*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	 
	     VV2_1 = t*Sx/(2.*Mp2)+t-Q2+sqly*sqlt*cpr1/(2.*Mp2)
	     VV2_2 = t*Sx/(2.*Mp2)+t-Q2+sqly*sqlt*cpr2/(2.*Mp2)

		 vv2min=min(vv2_1,vv2_2)
		 vv2max=max(vv2_1,vv2_2)
		 
c		 write(*,*)' dz0= ',dz0
		 
c		 write(*,*)z1,z2,sx,s-xx,(z1*s-xx/z2) ! lay,layz,sx*(z1*s-xx/z2)+2d0*(1d0/z2+z1)*mp2*q2!,az,cz**2+sz**2,cz,sz,cz**2+sz**2*cos(phi)**2-az**2,sqdz,cpr1,cpr2
		 
		 end

	
      subroutine dzsub(vv,q2,onmz,ipo,dzfun)
      implicit real*8(a-h,k-m,o-z)
	common/const/alpha,barn,mp,mp2,ml2,ml,pi

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


	double precision function siborn(s,q2,x,t,cphi,sphi,ia,ipol,eta)
        implicit none
	real*8 s,q2,x,t,cphi,sphi,vacpol,alphar,alpha,barn,mp2,ml,mp,ml2,pi,f1,f2,u0,w0,tt1,tt2,lauw,lay,xx,sx,sp,w2,tmin,tmax,factorvac
	real*8 tt3,tt4,f3,f4,phis,etak1,etak2,etap2,eta(3),sqly,sqlsxq
        real*8 tt10,tt20,tt30,tt40,tt1m,tt2m,tt3m,tt4m
	integer*4 ia,ipol
	common/const/alpha,barn,mp,mp2,ml2,ml,pi
	sx=q2/x
	xx=s-sx
	sp=s+xx
	w2=sx-q2+mp2
	lay=sx**2+4d0*mp2*q2
        sqly=sqrt(lay)
	tmin=-0.5d0*((sx-q2)*(sx+dsqrt(lay))+2d0*mp2*q2)/w2
	tmax=-0.5d0*((sx-q2)*(sx-dsqrt(lay))+2d0*mp2*q2)/w2
	call sffun(t,f1,f2,f3,f4)
        lauw=4d0*q2*w2*(s*xx-mp2*q2-ml2*lay)*(tmax-t)*(t-tmin)
	w0=-0.5d0*(q2+t)+0.5d0*sp/lay*(sx*(q2-t)+2d0*t*q2)+sqrt(lauw)/lay*cphi	
	u0=w0+q2+t
c	if(ia.eq.1)alphar=alpha	
caku	if(ia.eq.2)alphar=alpha/(1d0-alpha/2d0/pi*vacpol(q2))	
c	if(ia.eq.2)alphar=alpha/(1d0-alpha    /pi*vacpol(q2))	
	if(ia.eq.1)factorvac=1d0	
	if(ia.eq.2)factorvac=1d0/(1d0-alpha/2d0/pi*vacpol(q2))**2	
        if(ipol.eq.0)then
	 tt10=2d0*(u0**2+w0**2-2d0*q2*t)/u0/w0
	 tt20=0.5d0*mp2*tt10+t*(s**2+xx**2-q2*sx-s*w0-xx*u0)/u0/w0
         TT1m=2.d0*t*(1d0/u0**2+1d0/w0**2)
         TT2m=mp2/2d0*TT1m+(s**2+s*t)/u0**2+(xx**2-xx*t)/W0**2
          tt1=tt10+2d0*ml2*tt1m
          tt2=tt20+2d0*ml2*tt2m
	 siborn=-barn*factorvac*alpha**3*q2*(tt1*f1+tt2*f2)/4d0/pi/s**2/x**2/t/sqly	
c	 siborn=-barn*alpha**3*q2*(tt1*f1+tt2*f2)/4d0/pi/s**2/x**2/t/sqly	
        endif
        if(ipol.gt.0)then
          sqlsxq=sqrt(S*XX*Q2-mp2*q2**2)
          etak1=-( sqlsxq/sqly*eta(1)+(S*Sx+2d0*mp2*Q2)/(2d0*mp*sqly)*eta(3) )
          etak2=-( sqlsxq/sqly*eta(1)+(xx*Sx-2d0*mp2*Q2)/(2d0*mp*sqly)*eta(3) )
          etap2=-( sqrt(lauw)/(2.*sqlsxq*sqly)*(eta(1)*cphi+eta(2)*sphi)+(-t*Sx+2d0*mp2*(Q2-t))/(2d0*mp*sqly)*eta(3) )
          tt30=4d0*(2d0*xx*(u0-q2)-2d0*s*(w0+q2)+(w0+u0)*(q2-t))*etap2*mp/((t-4d0*mp**2)*u0*w0)
          tt40=-Mp2*TT30+2d0*mp/(u0*w0)*((q2-u0)*(t*etak2+etap2*xx)+(q2+w0)*(t*etak1+etap2*s))
           tt3m=4d0*((2d0*s+t)*(1d0/u0**2+Sx/s/w0**2)+1d0/w0**2*(-q2+1d0/s*(xx**2+(xx-t)**2)))*etap2*mp/(t-4d0*mp**2)
           tt4m=-Mp2*TT3m+2*mp*((sx+t)/s/w0**2*(t*etak2+etap2*xx)-(1d0/u0**2+1d0/w0**2)*(t*etak1+etap2*s))
            tt3=tt30+2d0*ml2*tt3m
            tt4=tt40+2d0*ml	2*tt4m
          siborn=-barn*factorvac*alpha**3*Q2/(4*pi*S**2*x**2*t*sqly)*(TT3*F3+TT4*F4)
        endif

c         write(*,'(9g11.3)')t,w2,cphi,w0,u0,lauw,lauw/(4*Q2*S**2*Sx**2),siborn
c         write(*,*)' JJJ ',-0.5d0*(q2+t)+0.5d0*sp/lay*(sx*(q2-t)+2d0*t*q2)
c         write(*,*)' KKK ',lauw/(4d0*Q2**2*S**2*Sx**2),tmax
c         write(*,*)' W0 ',w0,u0
c         write(*,*)' f1,f2 ',f1,f2
c         write(*,*)' si0 ',siborn,ia
c         write(*,*)siborn,-barn*factorvac*alpha**3*q2/4d0/pi/s**2/x**2/t/sqly,tt1*f1,tt2*f2,f2*t/4d0
c         stop
	return
	end
	 
        subroutine stestold(sum1,sum2,z1mc,z2mc,sborn)
        implicit none
        external peak
        real*8 sum1,sum2,rez1,rez2,z1mc,z2mc,si0,sborn
        real*8 alpha,barn,mp,mp2,ml2,ml,pi,s,q2,t,phi,x,xx,ssborn,eta,rc
     	common/const/alpha,barn,mp,mp2,ml2,ml,pi
        common/kinpoi/s,q2,t,phi,x,xx,ssborn,eta(3),ipol
        integer*4 isp,ipol  
        common/test/isp
         si0=sborn
         isp=1 
         call simpxx(z1mc,1d0-1d-7,1000,1d-4,peak,rez1)
c         call simpxx(z1mc,0.991998482d0,10000,1d-12,peak,rez1)
         isp=2 
         call simpxx(z2mc,1d0-1d-7,1000,1d-4,peak,rez2)
c         call simpxx(z2mc,0.963736969d0,1000,1d-12,peak,rez2)
c         call simpxx(0.81711592d0,0.963736969d0,1000,1d-12,peak,rez2)
         if(abs(sum1/rez1-1d0).gt.1d-3.or.abs(sum2/rez2-1d0).gt.1d-3)write(*,'(6g12.3)')sum1,rez1,sum1/rez1,sum2,rez2,sum2/rez2
         write(*,'(6g12.3)')sum1,rez1,sum1/rez1,sum2,rez2,sum2/rez2
		 
		 rc=+alpha/2d0/pi*(log(q2/ml2)    )*(rez1+rez2)/sborn
c		 rc=+alpha/2d0/pi*(log(q2/ml2)    )*(rez1     )/sborn
	         write(*,'(6g12.3)')sborn,rc*1d2,alpha/2d0/pi*log(q2/ml2)*(rez2    )/sborn*1d2
            stop	 
        end 

        subroutine stest(sum1,sum2,z1m,z2m,zspeak,zppeak,iaddcontr_s,iaddcontr_p,sborn)
        implicit none
        external peak,peaknew
        real*8 sum1,sum2,rez1,rez2,rez2n,rez1n,z1m,z2m,zspeak,zppeak,si0,sborn
        real*8 alpha,barn,mp,mp2,ml2,ml,pi,s,q2,t,phi,x,xx,ssborn,eta,rc
     	common/const/alpha,barn,mp,mp2,ml2,ml,pi
        common/kinpoi/s,q2,t,phi,x,xx,ssborn,eta(3),ipol
        integer*4 isp,ipol  ,iaddcontr_s,iaddcontr_p
        common/test/isp
         si0=sborn
         isp=1 
         call simpxx(z1m,1d0-1d-7,1000,1d-4,peak,rez1)
           rez1n=0d0
         if(iaddcontr_s.eq.1)call simpxx(zspeak+1d-7,z1m,1000,1d-4,peaknew,rez1n)
c         call simpxx(z1mc,0.991998482d0,10000,1d-12,peak,rez1)
         isp=2 
         call simpxx(z2m,1d0-1d-7,1000,1d-4,peak,rez2)
           rez2n=0d0
         if(iaddcontr_p.eq.1)call simpxx(zppeak+1d-7,z2m,1000,1d-4,peaknew,rez2n)
c         call simpxx(z2mc,0.963736969d0,1000,1d-12,peak,rez2)
c         call simpxx(0.81711592d0,0.963736969d0,1000,1d-12,peak,rez2)
		 rez1=rez1+rez1n
		 rez2=rez2+rez2n
         if(abs(sum1/rez1-1d0).gt.1d-2.or.abs(sum2/rez2-1d0).gt.1d-2)write(*,'(6g12.3)')sum1,rez1,sum1/rez1,sum2,rez2,sum2/rez2
c         write(*,'(4h s: ,6g12.3)')zspeak,z1m
c         write(*,'(4h p: ,6g12.3)')zppeak,z2m
c         write(*,'(6g12.3)')sum1,rez1,rez1n,sum1/rez1,sum2,rez2,rez2n,sum2/rez2
		 
		 rc=+alpha/2d0/pi*(log(q2/ml2)    )*(rez1+rez2)/sborn

		 sum1=rez1
		 sum2=rez2
		 
c		 rc=+alpha/2d0/pi*(log(q2/ml2)    )*(rez1     )/sborn
c	         write(*,'(6g12.3)')sborn,rc*1d2,alpha/2d0/pi*log(q2/ml2)*(rez2    )/sborn*1d2
c            stop	 
        end 

        double precision function peak(z)
        implicit none
        real*8 z,lay,layz,az,cz,sz,sqdz,cpr,spr,sbr,cphbr,sphbr,xsp,sigp,z1,z2,lat
        real*8 siborn,sigsp,s,q2,t,phi,x,xx,sborn,sx,alpha,barn,mp,mp2,ml2,ml,pi,eta,etaz(3)
        integer*4 isp,ipol  
 	common/const/alpha,barn,mp,mp2,ml2,ml,pi
        common/kinpoi/s,q2,t,phi,x,xx,sborn,eta(3),ipol
        common/test/isp
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
	    sigsp=siborn(z1*s,z1*q2/z2,xsp,t,cphbr,sphbr,1,ipol,etaz)
	    peak=(1d0+z**2)*(spr/sqdz*(xsp/x)**2*sigsp/z2-sborn)/(1d0-z)	
c	    peak=(1d0+z**2)*(spr/sqdz*(xsp/x)**2*sigsp/z2      )/(1d0-z)	
		
*	    write(*,*)z1,z2,spr,(-sqdz*cz-az*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
		
c	    peak=(1d0+z**2)*(                            -sborn)/(1d0-z)	
c       write(*,'(4g12.4)')z,z1,z2,peak
c         stop	
        end

        double precision function peaknew(z)
        implicit none
        real*8 z,lay,layz,az,cz,sz,sqdz,cpr,spr,sbr,cphbr,sphbr,cpr_,spr_,sbr_,cphbr_,sphbr_,xsp,sigp,z1,z2,lat
        real*8 siborn,sigsp,sigsp_,s,q2,t,phi,x,xx,sborn,sx,alpha,barn,mp,mp2,ml2,ml,pi,eta,etaz(3)
        integer*4 isp,ipol  
 	common/const/alpha,barn,mp,mp2,ml2,ml,pi
        common/kinpoi/s,q2,t,phi,x,xx,sborn,eta(3),ipol
        common/test/isp
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
	    sqdz=sqrt(cz**2+sz**2*cos(phi)**2-az**2)
	    cpr=(az*cz+sqdz*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    spr=(sqdz*cz-az*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    sbr=sqrt(1d0-(cpr*cz-spr*sz*cos(phi))**2)	
	    cphbr=(cz*spr*cos(phi)+sz*cpr)/sbr
            sphbr=spr*sin(phi)/sbr 
	    cpr_=(az*cz-sqdz*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    spr_=(-sqdz*cz-az*sz*cos(phi))/(cz**2+sz**2*cos(phi)**2)
	    sbr_=sqrt(1d0-(cpr_*cz-spr_*sz*cos(phi))**2)	
	    cphbr_=(cz*spr_*cos(phi)+sz*cpr_)/sbr_
            sphbr_=spr_*sin(phi)/sbr_ 
            etaz(1)=cz*eta(1)+sz*eta(3)  
            etaz(2)=eta(2)
            etaz(3)=-sz*eta(1)+cz*eta(3)
	    xsp=z1*q2/(z1*z2*s-xx)
	    sigsp=siborn(z1*s,z1*q2/z2,xsp,t,cphbr,sphbr,1,ipol,etaz)
	    sigsp_=siborn(z1*s,z1*q2/z2,xsp,t,cphbr_,sphbr_,1,ipol,etaz)
	    peaknew=(1d0+z**2)*(spr/sqdz*(xsp/x)**2*sigsp/z2+spr_/sqdz*(xsp/x)**2*sigsp_/z2)/(1d0-z)	
		
c		write(*,*)z1,z2,spr,spr_!,sbr,sbr_,sigsp,sigsp_
		
c	    peak=(1d0+z**2)*(spr/sqdz*(xsp/x)**2*sigsp/z2      )/(1d0-z)	
		
	
		
c	    peak=(1d0+z**2)*(                            -sborn)/(1d0-z)	
c        write(*,'(4g12.4)')z,z1,z2,peak
c         stop	
        end

		
		
****************** vacpol *************************************

      double precision function vacpol(t)
c contribution from vacuum polarization by leptons (suml) and hadrons (sumh)
      implicit real*8(a-h,l,m,o-z)
	common/const/alpha,barn,mp,mp2,ml2,ml,pi
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

      end

      subroutine sffun(t,f1,f2,f3,f4)
      implicit real*8(a-h,l,m,o-z)
	common/const/alpha,barn,mp,mp2,ml2,ml,pi
c       call ffpro(-t,fd,fp)
       call nuclff(t,fd,fp)
       f1=(fd+fp)**2
       f2=4d0/t*(fd**2-t*fp**2/4d0/mp2)
       f3=f1
       f4=4d0/t*(Fd+Fp)*(Fd+t/(4*mp2)*Fp)
c      write(*,*)' ff',fd,fp,f1,f2
      end

       subroutine ffpro(t,fd,fp)

      implicit real*8(a-h,l,m,o-z)
	common/const/alpha,barn,mp,mp2,ml2,ml,pi
      gep=1.2742/(1.+t/0.6394**2)-.2742/(1.+t/1.582**2)
      gmp=(1.3262/(1.+t/0.6397**2)-.3262/(1.+t/1.3137**2))*2.7921
c     gep=1./((1.+.61*t)*(1.+2.31*t)*(1.+.04*t))
c     gmp=amm*gep
	tap=t/4d0/mp2
	fd=(gep+tap*gmp)/(1d0+tap)
	fp=(gmp-gep)/(1d0+tap)
 
      end

      subroutine nuclFF( del2,fd,fp )
      implicit double precision (A-H,k,l,m,O-Z)
      parameter (Mv = 0.843D0, kp = 1.79285D0, kn = -1.91D0, mp=0.938272d0)
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
      Fd = (GE_p - delm*GM_p)/(1D0-delm)
      Fp = (GM_p - GE_p)/(1D0-delm)         
c        write(*,*)' ttf',fd,fp 
c      F1pn(2) = (GE_n - delm*GM_n)/(1D0-delm)
c      F2pn(2) = (GM_n - GE_n)/(1D0-delm)         
      return
      end

	
c$nodebug
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
C *						      *
      FUNCTION URAND(IY)
C *						      *
C *   This is a standard pseudo-random generator      *
C *   that work on IBM-370 and IBM-PC. We don't       *
C *   know does it work on SUN? 		      *
C *						      *
C * * * * * * * * * * * * * * * * * * * * * * * * * * *
      REAL*4  URAND,S
      INTEGER*4 IY
      INTEGER*4 A,C,MASK
      PARAMETER (A  = 843314861)
      PARAMETER (C  = 453816693)
      PARAMETER (S  = 4.6566128E-10)

      IY=IAND(A*IY+C,Z'7FFFFFFF')
      URAND=FLOAT(IY)*S
      END

      subroutine simps(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
c simps
c a1,b1 -the limits of integration
c h1 -an initial step of integration
c reps1,aeps1 - relative and absolute precision of integration
c funct -a name of function subprogram for calculation of integrand +
c x - an argument of the integrand
c ai - the value of integral
c aih- the value of integral with the step of integration
c aiabs- the value of integral for module of the integrand
c this subrogram calculates the definite integral with the relative or
c absolute precision by simpson+s method with the automatical choice
c of the step of integration
c if aeps1    is very small(like 1.e-17),then calculation of integral
c with reps1,and if reps1 is very small (like 1.e-10),then calculation
c of integral with aeps1
c when aeps1=reps1=0. then calculation with the constant step h1
c
      implicit real*8(a-h,o-z)
      dimension f(7),p(5)
      h=dsign(h1,b1-a1)
      s=dsign(1.d0,h)
      a=a1
      b=b1
      ai=0.d0
      aih=0.d0
      aiabs=0.d0
      p(2)=4.d0
      p(4)=4.d0
      p(3)=2.d0
      p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=dabs(reps1)
      aeps=dabs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a
      c=0.d0
      f(1)=funct(x)/3.
    4 x0=x
      if((x0+4.*h-b)*s) 5,5,6
    6 h=(b-x0)/4.
      if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1)
      di3=dabs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
   11 f(k)=funct(x)/3.
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.*f(3)+f(5))*2.*h
      di2=di2*h
      di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=dabs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=dabs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.) 17,14,14
   17 h=2.*h
      f(1)=f(5)
      f(2)=f(6)
      f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5)
      f(3)=f(6)
      f(5)=f(7)
      f(2)=10.d16
      f(4)=10.d16
      f(6)=10.d16
      f(7)=10.d16
   18 di1=di2+(di2-di1)/15.
      ai=ai+di1
      aih=aih+di2
      aiabs=aiabs+di3
      go to 22
   21 h=h/2.
      f(7)=f(5)
      f(6)=f(4)
      f(5)=f(3)
      f(3)=f(2)
      f(2)=10.d16
      f(4)=10.d16
      x=x0
      c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end

      subroutine simpxx(a,b,np,ep,func,res)
      implicit real*8 (a-h,o-z)
      external func
      step=(b-a)/np
      call simps(a,b,step,ep,1d-18,func,ra,res,r2,r3)
      end

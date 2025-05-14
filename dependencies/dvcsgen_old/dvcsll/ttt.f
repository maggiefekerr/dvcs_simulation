      
      subroutine sifo(vv2cut,iborn,imeth)
      implicit none
      real*8 alpha,barn,mp,mp2,ml2,ml,pi
      common/const/alpha,barn,mp,mp2,ml2,ml,pi
      real*8 s,q2,t,phi,x,xx,sborn,eta
      integer*4 ipol,iy
      common/kinpoi/s,q2,t,phi,x,xx,sborn,eta(3),ipol
      real*8 vv2cut,zspeak,zppeak,siborn,sum1,sum2,z1,z2,layz,az,cz,sz,sqdz,cpr,spr,sbr,cphbr,sphbr,xsp,sig0s,sig0p,sirad,sirad_t,siexp,sum1add,sum2add,z1mc,z2mc,lat,lay,sx,etaz(3),z1m,z2m,vacpol
      integer*4 iborn,imeth,iaddcontr_s,iaddcontr_p,in,nev,urand
      
      call addcon(zspeak,zppeak,iaddcontr_s,iaddcontr_p,z1m,z2m,lat,lay,sx)
       
      if(imeth.ne.0)then
	  iy=12345678
      nev=2**abs(imeth)
      sum1=0d0
      sum2=0d0
      sborn=siborn(s,q2,x,t,cos(phi),sin(phi),1,ipol,eta) 
c            write(*,*)' sborn= ',sborn
      do in=1,nev
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
c        xb1=z1/(1d0/xb-s*(1d0-z1)/Q2)
        sig0s=siborn(z1*s,z1*q2,xsp,t,cphbr,sphbr,1,ipol,etaz)
        sum1add=(1d0-z1mc)*(1d0+z1**2)*(spr/sqdz*(xsp/x)**2*sig0s-sborn)/(1d0-z1)        
        sum1=sum1+sum1add/dble(nev)
c        sig0s=siborn(z1*s,z1*q2,xsp,t,cos(phi),1,ipol)
c      sum1=sum1+(1d0-z1m)/dble(nev)*(1d0+z1**2)*(spr/sqdz*(xsp/x)**2*sig0s-sborn)/(1d0-z1)        

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
c        xb2=1d0/(1d0/xb-(1d0-z2)*s/Q2)
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
      enddo    
     
      sirad= alpha/Pi*vacpol(-t)*sborn
c     &      +alpha/2d0/pi*(log(q2/ml2)-1d0)
     &      +alpha/2d0/pi*(log(q2/ml2)    )
     &          *(sum1+sum2+(z1m*(z1m+2d0)/2d0+z2m*(z2m+2d0)/2d0
     &          +2d0*log(1d0-z1m)+2d0*log(1d0-z2m))*sborn)

          siexp=  exp(alpha/pi*log(q2/ml2)*log((1d0-z1m)*(1d0-z2m)))*sborn   + alpha    /Pi*vacpol(-t)*sborn
c     &      +alpha/2d0/pi*(log(q2/ml2)-1d0)
     &      +alpha/2d0/pi*(log(q2/ml2)    )
     &          *(sum1+sum2+(z1m*(z1m+2d0)/2d0+z2m*(z2m+2d0)/2d0
     &          )*sborn)

      endif
     
      if(imeth.ge.0)then     
           call stest(sum1,sum2,z1m,z2m,zspeak,zppeak,iaddcontr_s,iaddcontr_p,sborn)
          sirad_t= alpha    /Pi*vacpol(-t)*sborn
     &      +alpha/2d0/pi*(log(q2/ml2)    )
     &          *(sum1+sum2+(z1m*(z1m+2d0)/2d0+z2m*(z2m+2d0)/2d0
     &          +2d0*log(1d0-z1m)+2d0*log(1d0-z2m))*sborn)
      endif 

      if(imeth.lt.0)then     
           write(61,'(4f8.3,f14.3,6f11.3)')x,q2,t,phi,sborn,sirad/sborn*1d2
           write(*,'(3f8.3,f8.1,10g12.4)')x,q2,t,phi*180/pi,sborn,sirad/sborn*1d2,1d0+sirad/sborn,siexp/sborn
      endif
      if(imeth.eq.0)then     
           write(61,'(4f8.3,f14.3,6f11.3)')x,q2,t,phi,sborn,sirad_t/sborn*1d2
           write(*,'(3f8.3,f8.1,10g12.4)')x,q2,t,phi*180/pi,sborn,sirad_t/sborn*1d2,1d0+sirad_t/sborn
      endif
      if(imeth.lt.0)then     
           write(61,'(4f8.3,f14.3,6f11.3)')x,q2,t,phi,sborn,sirad/sborn*1d2
           write(*,'(3f8.3,f8.1,10g12.4)')x,q2,t,phi*180/pi,sborn,sirad/sborn*1d2,1d0+sirad/sborn,siexp/sborn,sirad_t/sirad
      endif
      end


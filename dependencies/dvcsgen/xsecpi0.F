c exc pi0 x-section
c
cinput:
cdel2=t (negative GeV^2)
cxb,Q2 x and Q^2 (GeV^2)
c Phi_g angle in the photon frame (radians)
c E energy of the electron in GeV
c heli electron helicity -1.0 or +1.0
c
      real function xsecpi0(del2,xb,Q2, Phi_g,E,heli)
      implicit none
      real E, xb, Q2, del2,Phi_g,heli
      real Mp, mele, pi,mpi0
      parameter (Mp=0.93827)
      parameter (mele=0.000511)
      parameter (pi=3.1415926536)
      parameter (mpi0=0.135)
c
c
      real nu,W2,W,qmod,E1cm,P1cm,E2cm,P2cm,del2max,del2min
      real yb,ymax,ycol,dpi0u,sinpi0
      real  xmin1,xmax1
      real UFx,PFxq2,Tfac,dubbs,ddbbs
      real tt,tQ2,A1_PHI,A2_PHI,BT_Q2

c
      xsecpi0=0.0
      xmin1 = Q2/(2.0*Mp*E)
      xmax1 = 1.0
      nu  = Q2/(2D0*Mp*xb)
      W2  = Mp**2 + 2.0*Mp*nu - Q2
      W   = sqrt(W2)
      qmod = sqrt(nu**2 + Q2)
      tt=-del2
      tQ2=Q2
c
      E1cm = Mp*(Mp + nu)/W
      P1cm = Mp*qmod/W
      E2cm = (W2 + Mp**2 -mpi0**2)/(2D0*W)
      P2cm = (W2 - Mp**2)/(2D0*W)
      del2max = 2.0*(Mp**2 - E1cm*E2cm - P1cm*P2cm)
      del2min = 2.0*(Mp**2 - E1cm*E2cm + P1cm*P2cm)
c check the kinematics
c
      if( xb.le.xmin1 .or. xb.gt.xmax1 ) return           !    x  out of range
      if( del2.ge.del2min .or. del2.le.del2max ) return   ! delta out of range
      yb=nu/E
c
       sinpi0=2*xb*(1-xb)**5
       PFxq2=84*xb*xb*xb*(1-xb) !(1-xb)**2/(Mp*Mp*xb/Q2+1-xb)**2
       Tfac=exp(BT_Q2(tQ2)*(del2)) !-del2min))  
       UFx=(2*dubbs(xb)+ddbbs(xb))**2
       dpi0u=0.3*Tfac*UFx*PFxq2
       dpi0u=dpi0u/Q2/Q2/(Q2+Mp**2)**2 ! vpk
     6(1.0+A1_PHI(tt)*cos(Phi_g)+A2_PHI(tt)*cos(2*Phi_g)
     6    +sinpi0*heli*sin(Phi_g))
       xsecpi0=dpi0u
      return
      end
c
c     support functions
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
c  from vpk
c
c
c    dN/dt = exp(-BT_Q2(Q2))
c
      REAL FUNCTION bt_Q2(Q2)
      IMPLICIT NONE
      REAL Q2,P(3),tmp

      P(1)=2.923
      P(2)=0.5791
      P(3)=1.402

      IF(Q2.LT.2.6) THEN
        tmp=P(1)-P(2)*Q2
      ELSE
        tmp=P(3)
      ENDIF
        BT_Q2=tmp
      RETURN
      END
C
C     dN/dphi= 1.+ A1*Cos(phi) + A2*Cos(2phi)
C
      REAL FUNCTION A1_PHI(T)
      IMPLICIT NONE
      REAL T, P(3)
      P(1)=-0.5011
      P(2)=1.419
      P(3)=0.4837

      IF(T.LT.0.7) THEN
         A1_PHI=P(1)+P(2)*T
      ELSE
         A1_PHI=P(3)
      ENDIF
      RETURN
      END
c----
      REAL FUNCTION A2_PHI(T)
      IMPLICIT NONE
      REAL T, P(1)

      A2_PHI=-0.3490

      RETURN
      END

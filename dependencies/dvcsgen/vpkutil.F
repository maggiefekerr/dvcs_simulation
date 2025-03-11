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

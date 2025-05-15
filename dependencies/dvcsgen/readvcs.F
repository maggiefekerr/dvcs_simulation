          PROGRAM DVCS
c
          IMPLICIT NONE
c
          INTEGER N,K,IEV,J,iflag,jl
          REAL P,V
          COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
          
          DO IEV=1,10000000
          call GETEVENT(iflag)    ! iflag=0 good event
          if(iflag.gt.0) Then
              print *,'Total events ',IEV,IEV-1,'   flag->',iflag
              goto 999  ! no more events
          else          ! process event
            do j=1,N
             write(6,102) iev,N,k(j,2),iflag,'*',(p(j,jl),jl=1,5)
            enddo 
          endif
c
          enddo
 101      FORMAT(2x,I10)
 102      FORMAT(2x,4I5,A,10F10.4)
 999      CONTINUE
          print *,'END',IEV
          close(41)
          END
c
c
c
c
          SUBROUTINE GETEVENT(iflag)
          IMPLICIT NONE
c
          INTEGER N,K,iflag,j,jl
          REAL P,V
          COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
          LOGICAL lOK,lFirst
          DATA lFirst /.TRUE./
          save lFirst
          iflag=0
c
                if (lFirst) then   ! open the file if first time
                 OPEN(unit=41,file='dvcsgen.dat',status='unknown')
                 lFirst=.FALSE.
                endif
c
              read(41,101,end=99) N
              do  j=1,N
              read(41,102,end=99) k(j,2),(p(j,jl),jl=1,5) 
              enddo
 101      FORMAT(2x,I10)
 102      FORMAT(2x,I5,10F10.4)
          RETURN
 99       CONTINUE
          iflag=1
          RETURN
          END     

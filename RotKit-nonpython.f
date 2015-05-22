      subroutine sphcar(lat,long,x1,x2,x3)
      real*8 lat,long
      x1=cos(lat)*cos(long)
      x2=cos(lat)*sin(long)
      x3=sin(lat)
      return
      end

      subroutine carsph(x1,x2,x3,a,b)
      real*8 a,b
      parameter (pi=dacos(-1.d0))
      a = asin(x3)
      if(.not.(x1.eq.0.d0))goto 23047
      if(.not.(x2.gt.0.d0))goto 23049
      b = pi*0.5
23049 continue
      if(.not.(x2.lt.0.d0))goto 23051
      b = pi*1.5
23051 continue
      goto 23048
23047 continue
      b = atan(x2/x1)
      if(.not.(x1.lt.0.d0))goto 23053
      b = b+pi
23053 continue
      if(.not.(x1.gt.0.d0.and.x2.lt.0.d0))goto 23055
      b = b+(2*pi)
23055 continue
23048 continue
      return
      end
      
      subroutine dmatrix(c1,c2,c3,w,rmatrix)
      real*8 w
      real*8 rmatrix(3,3)
      ss = sin(w*0.5)
      r = c1*ss
      s = c2*ss
      t = c3*ss
      u = cos(w*0.5)
      rmatrix(1,1) = r*r-s*s-t*t+u*u
      rmatrix(1,2) = 2.0*(r*s-t*u)
      rmatrix(1,3) = 2.0*(t*r+s*u)
      rmatrix(2,1) = 2.0*(r*s+t*u)
      rmatrix(2,2) = s*s-t*t-r*r+u*u
      rmatrix(2,3) = 2.0*(s*t-r*u)
      rmatrix(3,1) = 2.0*(t*r-s*u)
      rmatrix(3,2) = 2.0*(s*t+r*u)
      rmatrix(3,3) = t*t-r*r-s*s+u*u
      return
      end 

      subroutine rot3(a,b,w,r,s)
      real*8 w,r,s
      call sphcar(a,b,x1,x2,x3)
      t = x1
      x1 = x1*cos(w)-x2*sin(w)
      x2 = t*sin(w)+x2*cos(w)
      call carsph(x1,x2,x3,r,s)
      return
      end
      
      subroutine rot2(a,b,w,r,s)
      real*8 w,r,s
      call sphcar(a,b,x1,x2,x3)
      t = x1
      x1 = x1*cos(w)+x3*sin(w)
      x3 = x3*cos(w)-t*sin(w)
      call carsph(x1,x2,x3,r,s)
      return
      end
    
      subroutine rotp(a,b,w,c,d,f,g)
      implicit real*8 (a-h,o-z)
      real*8 c1,d1,c2,d2,c3,d3,c4,d4,f,g
      parameter (pi=dacos(-1.d0))
      call rot3(c,d,-b,c1,d1)
      call rot2(c1,d1,a-(pi*0.5),c2,d2)
      call rot3(c2,d2,w,c3,d3)
      call rot2(c3,d3,(pi*0.5)-a,c4,d4)
      call rot3(c4,d4,b,f,g)
      return
      end

      subroutine dmatmul(mat1,row1,col1,mat2,row2,col2,mulmat)
      integer row1,col1,row2,col2,i,j,k
      real*8 mat1(row1,col1),mat2(row2,col2)
      real*8 mulmat(row1,col2)
      do23059 i = 1,row1 
      do23061 j = 1,col2 
      mulmat(i,j) = 0.0d0
23061 continue
23059 continue
      if(.not.(col1 .ne. row2))goto 23063
      write(6,1000)
1000  format(/" ERROR -- Matrices not compatible for multiplication"/)
      goto 23064
23063 continue
      do23065 i = 1,row1 
      do23067 j = 1,col2 
      do23069 k = 1,col1 
      mulmat(i,j) = mulmat(i,j) + mat1(i,k) * mat2(k,j)
23069 continue
23067 continue
23065 continue
23064 continue
      return
      end

      SUBROUTINE DJACOBI(A,N,NP,D,V,NROT)
      integer N,NP
      real*8 A(NP,NP)
      real*8 H,THETA
      integer NROT
      PARAMETER (NMAX=100)
      REAL*8 D(NP),V(NP,NP),B(NMAX),Z(NMAX)
      DO 12 IP=1,N
        DO 11 IQ=1,N
          V(IP,IQ)=0.d0
11      CONTINUE
        V(IP,IP)=1.d0        
12    CONTINUE
      DO 13 IP=1,N
        B(IP)=A(IP,IP)
        D(IP)=B(IP)
        Z(IP)=0.
13    CONTINUE
      NROT=0
      DO 24 I=1,50
        SM=0.
        DO 15 IP=1,N-1
          DO 14 IQ=IP+1,N
            SM=SM+ABS(A(IP,IQ))
14        CONTINUE
15      CONTINUE
        IF(SM.EQ.0.d0)GOTO 30
        IF(I.LT.4)THEN
          TRESH=0.2d0*SM/N**2
        ELSE
          TRESH=0.d0
        ENDIF
        DO 22 IP=1,N-1
          DO 21 IQ=IP+1,N
            G=100.*ABS(A(IP,IQ))
            IF((I.GT.4).AND.(ABS(D(IP))+G.EQ.ABS(D(IP)))
     *         .AND.(ABS(D(IQ))+G.EQ.ABS(D(IQ))))THEN
              A(IP,IQ)=0.
            ELSE IF(ABS(A(IP,IQ)).GT.TRESH)THEN
              H=D(IQ)-D(IP)
              IF(ABS(H)+G.EQ.ABS(H))THEN
                T=A(IP,IQ)/H
              ELSE
                THETA=0.5d0*H/A(IP,IQ)
                T=1.d0/(ABS(THETA)+DSQRT(1.d0+THETA**2))
                IF(THETA.LT.0.d0)T=-T
              ENDIF
              C=1.d0/DSQRT(1.0d0+T**2)
              S=T*C
              TAU=S/(1.d0+C)
              H=T*A(IP,IQ)
              Z(IP)=Z(IP)-H
              Z(IQ)=Z(IQ)+H
              D(IP)=D(IP)-H
              D(IQ)=D(IQ)+H
              A(IP,IQ)=0.d0
              DO 16 J=1,IP-1
                G=A(J,IP)
                H=A(J,IQ)
                A(J,IP)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
16            CONTINUE
              DO 17 J=IP+1,IQ-1
                G=A(IP,J)
                H=A(J,IQ)
                A(IP,J)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
17            CONTINUE
              DO 18 J=IQ+1,N
                G=A(IP,J)
                H=A(IQ,J)
                A(IP,J)=G-S*(H+G*TAU)
                A(IQ,J)=H+S*(G-H*TAU)
18            CONTINUE
              DO 19 J=1,N
                G=V(J,IP)
                H=V(J,IQ)
                V(J,IP)=G-S*(H+G*TAU)
                V(J,IQ)=H+S*(G-H*TAU)
19            CONTINUE
              NROT=NROT+1
            ENDIF
21        CONTINUE
22      CONTINUE
        DO 23 IP=1,N
          B(IP)=B(IP)+Z(IP)
          D(IP)=B(IP)
          Z(IP)=0.d0
23      CONTINUE
24    CONTINUE
      PAUSE '50 iterations should never happen'
30    DO 31 IP=1,N-1
      DO 31 IQ=IP+1,N
31    A(IP,IQ)=A(IQ,IP)
      DO 33 I=1,N-1
        K=I
        P=D(I)
        DO 41 J=I+1,N
           IF (D(J).LT.P) THEN
              K=J
              P=D(J)
           ENDIF
41       CONTINUE
        IF (K.NE.I) THEN
           D(K)=D(I)
           D(I)=P
           DO 32 J=1,N
              P=V(J,I)
              V(J,I)=V(J,K)
32             V(J,K)=P
         ENDIF
33    CONTINUE
      RETURN
      END      
      
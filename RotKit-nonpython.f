      subroutine sphcar(lat,long,x1,x2,x3)
      DOUBLE PRECISION lat,long
      x1=cos(lat)*cos(long)
      x2=cos(lat)*sin(long)
      x3=sin(lat)
      return
      end

      subroutine carsph(x1,x2,x3,a,b)
      DOUBLE PRECISION a,b
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
      DOUBLE PRECISION w
      DOUBLE PRECISION rmatrix(3,3)
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
      DOUBLE PRECISION w,r,s,a,b
      call sphcar(a,b,x1,x2,x3)
      t = x1
      x1 = x1*cos(w)-x2*sin(w)
      x2 = t*sin(w)+x2*cos(w)
      call carsph(x1,x2,x3,r,s)
      return
      end
      
      subroutine rot2(a,b,w,r,s)
      DOUBLE PRECISION w,r,s,a,b
      call sphcar(a,b,x1,x2,x3)
      t = x1
      x1 = x1*cos(w)+x3*sin(w)
      x3 = x3*cos(w)-t*sin(w)
      call carsph(x1,x2,x3,r,s)
      return
      end
    
      subroutine rotp(a,b,w,c,d,f,g)
      implicit DOUBLE PRECISION (a-h,o-z)
      DOUBLE PRECISION c1,d1,c2,d2,c3,d3,c4,d4,f,g,a
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
      DOUBLE PRECISION mat1(row1,col1),mat2(row2,col2)
      DOUBLE PRECISION mulmat(row1,col2)
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
      DOUBLE PRECISION A(NP,NP)
      DOUBLE PRECISION H,THETA
      integer NROT
      PARAMETER (NMAX=100)
      DOUBLE PRECISION D(NP),V(NP,NP),B(NMAX),Z(NMAX)
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

c Subroutines from this point up to squash from Pavel Doubrovine's
c ibfrlib0.3.f minus the finb2 routine. Preamble below from that file.
c-----------------------------------------------------------------------
c
c       Subroutines used to interpolate between finite rotations
c
c       Written by: Pavel V. Doubrovine (2006-2010)
c
c       Release 0.3
c       Last modified: 22.12.2010
c
c       NOTE: This software comes with no warranty, use it on your
c       own risk. Please DO NOT REDISTRIBUTE without contacting
c       the author. PLEASE REPORT BUGS to: paveld@fys.uio.no
c
c       The interpolation method is fully described in the auxiliary
c       material of
c
c       Doubrovine, P.V., and J.A. Tarduno (2008), Linking the Late
c       Cretaceous to Paleogene Pacific plate and Atlantic bordering
c       continents using plate circuits and paleomagnetic data,
c       J. Geophys. Res., 113, B07104, doi:10.1029/2008JB005584.
c
c
c       **** BUGS FIXED in release 0.3 ***
c
c       1) Subroutine <matrixl>: The last call to <mulmm>
c       was:            call mulmm(temp1,temp3,l)
c       correct usage:  call mulmm(temp1,temp3,l,3,3)
c       2) got rid of mixed arithmetic assignments
c
c
c
c-------------------------------------------------------------------
c Calculates the rotation matrix a(3,3) from the Euler
c parameters (lat,lon,rho). All angular values in radians.
c Formulas from Cox and Hart (1986), p. 227.
c-------------------------------------------------------------------
        SUBROUTINE ep2r(lat,lon,rho,a)
        DOUBLE PRECISION lat,lon,rho,a(3,3),ev(3)
        ev(1)=cos(lat)*cos(lon)
        ev(2)=cos(lat)*sin(lon)
        ev(3)=sin(lat)
        a(1,1)=ev(1)*ev(1)*(1.-cos(rho))+cos(rho)
        a(1,2)=ev(1)*ev(2)*(1.-cos(rho))-ev(3)*sin(rho)
        a(1,3)=ev(1)*ev(3)*(1.-cos(rho))+ev(2)*sin(rho)
        a(2,1)=ev(2)*ev(1)*(1.-cos(rho))+ev(3)*sin(rho)
        a(2,2)=ev(2)*ev(2)*(1.-cos(rho))+cos(rho)
        a(2,3)=ev(2)*ev(3)*(1.-cos(rho))-ev(1)*sin(rho)
        a(3,1)=ev(3)*ev(1)*(1.-cos(rho))-ev(2)*sin(rho)
        a(3,2)=ev(3)*ev(2)*(1.-cos(rho))+ev(1)*sin(rho)
        a(3,3)=ev(3)*ev(3)*(1.-cos(rho))+cos(rho)
        return
        END
c-------------------------------------------------------------------
c Matrix transpose: b=a**t (a(n,n), b(n,n))
c-------------------------------------------------------------------
        SUBROUTINE matt(a,b,n,np)
        INTEGER n,np
        DOUBLE PRECISION a(np,np), b(np,np)
        INTEGER i,j
        do 11 i=1,n
          do 12 j=1,n
            b(i,j)=a(j,i)
12        continue
11      continue
        return
        END
c-------------------------------------------------------------------
c Multiplication a*b=c (a(n,n), b(n,n), b(n,n))
c-------------------------------------------------------------------
        SUBROUTINE mulmm(a,b,c,n,np)
        INTEGER n,np
        DOUBLE PRECISION a(np,np), b(np,np), c(np,np)
        INTEGER i,j,k
        do 11 i=1,n
          do 12 j=1,n
            c(i,j)=0.
            do 13 k=1,n
              c(i,j)=c(i,j)+a(i,k)*b(k,j)
13          continue
12        continue
11      continue
        return
        END
c-------------------------------------------------------------------
c Subroutine to convert a rotation matrix (a(3,3)) into the
c Euler pole (lat,lon) and rotation angle (rho). Returns
c lat, lon, rho, and rotation pseudovector r(3).
c All angular values are in radians.
c-------------------------------------------------------------------
        SUBROUTINE rpv(a,lat,lon,rho,r)
        INTEGER i
        DOUBLE PRECISION a(3,3), lat, lon, rho, r(3)
        DOUBLE PRECISION tra, x(3), sroot, pi
        PARAMETER (pi=3.14159265358979323846264D+0)
        tra=0.
        do 11 i=1,3
          tra=tra+a(i,i)
11      continue
        x(1)=a(3,2)-a(2,3)
        x(2)=a(1,3)-a(3,1)
        x(3)=a(2,1)-a(1,2)
        sroot=sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
c
        if(abs(tra-1.).eq.0.) then
            rho=pi*0.5
          else
            rho=atan(sroot/(tra-1.))
            if(rho.lt.0.) rho=rho+pi
        endif
c
        if(sroot.eq.0.) then
            lat=0.
          else
            lat=asin(x(3)/sroot)
        endif
c
        if(sroot.eq.0) then
          lon=0.
        else
          if(x(1).eq.0.) then
            if(x(2).ge.0.) then
              lon=0.5*pi
            else
              lon=-0.5*pi
            endif
          else
            lon=atan(x(2)/x(1))
            if(x(1).lt.0) lon=lon+pi
          endif
        endif
c
        r(1)=rho*cos(lat)*cos(lon)
        r(2)=rho*cos(lat)*sin(lon)
        r(3)=rho*sin(lat)
        return
        END
c-------------------------------------------------------------------
c Calculates rotation matrix a=Phi(h) from the Euler pole lat,lon
c and rotation angle
c-------------------------------------------------------------------
        SUBROUTINE phi(lat,lon,rho,a)
        DOUBLE PRECISION lat,lon,rho,p(3),a(3,3),b,c,d
        p(1)=cos(lat)*cos(lon)
        p(2)=cos(lat)*sin(lon)
        p(3)=sin(lat)
        c=cos(rho)
        b=sin(rho)
        d=1.-c
        a(1,1)=p(1)*p(1)*d+c
        a(1,2)=p(1)*p(2)*d-p(3)*b
        a(1,3)=p(1)*p(3)*d+p(2)*b
        a(2,1)=p(2)*p(1)*d+p(3)*b
        a(2,2)=p(2)*p(2)*d+c
        a(2,3)=p(2)*p(3)*d-p(1)*b
        a(3,1)=p(3)*p(1)*d-p(2)*b
        a(3,2)=p(3)*p(2)*d+p(1)*b
        a(3,3)=p(3)*p(3)*d+c
        return
        END
c-------------------------------------------------------------------
c Calculates matrix L=(A1**t)*(Fs1**-1)*Ds1*(Ds**-1)*Fs*A1.
c Uses MATRIXF, MATRIXD, MATRIXFINV, MATRIDINV, MULMM, MATT
c Input a1(3,3), rs(3), rs1(3); Output: l(3,3)
c-------------------------------------------------------------------
        SUBROUTINE matrixl(a1,rs,rs1,l)
        DOUBLE PRECISION a1(3,3),rs(3),rs1(3),l(3,3)
        DOUBLE PRECISION temp1(3,3),temp2(3,3),temp3(3,3)
        call matrixf(rs,temp1)
        call mulmm(temp1,a1,temp2,3,3)
        call matrixdinv(rs,temp1)
        call mulmm(temp1,temp2,temp3,3,3)
        call matrixd(rs1,temp1)
        call mulmm(temp1,temp3,temp2,3,3)
        call matrixfinv(rs1,temp1)
        call mulmm(temp1,temp2,temp3,3,3)
        call matt(a1,temp1,3,3)
        call mulmm(temp1,temp3,l,3,3)
        return
        END
c-------------------------------------------------------------------
c Calculates matrix of partial derivatives D(i,j)=dsr(i)/dr(j)
c Input: rotation vector r(3); Output: matrix d(3,3)
c-------------------------------------------------------------------
        SUBROUTINE matrixd(r,d)
        INTEGER i,j
        DOUBLE PRECISION r(3),d(3,3),normr,sf1,cf2,sf2
        normr=sqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
        if(normr.ne.0) then
          sf1=sin(normr/2.)/normr
          cf2=cos(normr/2.)/(2.*normr*normr)
          sf2=sin(normr/2.)/(normr*normr*normr)
          do 11 i=1,3
            do 12 j=1,3
              if(i.eq.j) then
                d(i,j)=sf1+(cf2-sf2)*r(i)*r(j)
              else
                d(i,j)=(cf2-sf2)*r(i)*r(j)
              endif
12          continue
11        continue
        else
          do 13 i=1,3
            do 14 j=1,3
              if(i.eq.j) then
                d(i,j)=0.5
              else
                d(i,j)=0.
              endif
14          continue
13        continue
        endif
        return
        END
c-------------------------------------------------------------------
c Calculates inverse matrix of partial derivatives
c D-1(i,j)=dr(i)/dsr(j). Uses SQUASH.
c Input: rotation vector r(3); Output: matrix dinv(3,3)
c-------------------------------------------------------------------
        SUBROUTINE matrixdinv(r,dinv)
        INTEGER i,j
        DOUBLE PRECISION r(3),dinv(3,3),sr(3),normsr,sf1,cf2,sf2
        call squash(r,sr)
        normsr=sqrt(sr(1)*sr(1)+sr(2)*sr(2)+sr(3)*sr(3))
        if(normsr.ne.0) then
          sf1=2.*asin(normsr)/normsr
          cf2=2./(normsr*normsr*sqrt(1.-normsr*normsr))
          sf2=2.*asin(normsr)/(normsr*normsr*normsr)
          do 11 i=1,3
            do 12 j=1,3
              if(i.eq.j) then
                dinv(i,j)=sf1+(cf2-sf2)*sr(i)*sr(j)
              else
                dinv(i,j)=(cf2-sf2)*sr(i)*sr(j)
              endif
12          continue
11        continue
        else
          do 13 i=1,3
            do 14 j=1,3
              if(i.eq.j) then
                dinv(i,j)=2.
              else
                dinv(i,j)=0.
              endif
14          continue
13        continue
        endif
        return
        END
c-------------------------------------------------------------------
c Calculates matrix F=sqrt(1-|sr|**2)I+M(sr), where sr(3) is the
c squashed rotation vector. Uses SQUASH.
c Input: rotation vector r(3); Output: matrix f(3,3)
c-------------------------------------------------------------------
        SUBROUTINE matrixf(r,f)
        DOUBLE PRECISION r(3),f(3,3),sr(3),sroot
        call squash(r,sr)
        sroot=sqrt(1.-(sr(1)*sr(1)+sr(2)*sr(2)+sr(3)*sr(3)))
        do 11 i=1,3
          f(i,i)=sroot
11      continue
        f(1,2)=-sr(3)
        f(1,3)= sr(2)
        f(2,1)= sr(3)
        f(2,3)=-sr(1)
        f(3,1)=-sr(2)
        f(3,2)= sr(1)
        return
        END
c-------------------------------------------------------------------
c Calculates the inverse of matrix F. Uses SQUASH.
c Input: rotation vector r(3); Output: matrix finv(3,3)
c-------------------------------------------------------------------
        SUBROUTINE matrixfinv(r,finv)
        DOUBLE PRECISION r(3),finv(3,3),sr(3),sroot
        call squash(r,sr)
        sroot=sqrt(1.-(sr(1)*sr(1)+sr(2)*sr(2)+sr(3)*sr(3)))
        finv(1,1)= sroot+sr(1)*sr(1)/sroot
        finv(1,2)= sr(3)+sr(1)*sr(2)/sroot
        finv(1,3)=-sr(2)+sr(1)*sr(3)/sroot
        finv(2,1)=-sr(3)+sr(2)*sr(1)/sroot
        finv(2,2)= sroot+sr(2)*sr(2)/sroot
        finv(2,3)= sr(1)+sr(2)*sr(3)/sroot
        finv(3,1)= sr(2)+sr(3)*sr(1)/sroot
        finv(3,2)=-sr(1)+sr(3)*sr(2)/sroot
        finv(3,3)= sroot+sr(3)*sr(3)/sroot
        return
        END
c-------------------------------------------------------------------
c Calculates "squashed" rotation pseudovector
c Input: rotation vector r(3); Output: squashed vector sr(3)
c-------------------------------------------------------------------
        SUBROUTINE squash(r,sr)
        INTEGER i
        DOUBLE PRECISION r(3),sr(3),normr,sfact
        normr=sqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
        if(normr.ne.0.) then
            sfact=sin(normr/2.)/normr
          else
            sfact=0.5
        endif
        do 11 i=1,3
          sr(i)=sfact*r(i)
11      continue
        return
        END

      
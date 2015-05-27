! -*- f90 -*-
! File rotatep.f
      subroutine rotatepts(points,npt,poles,npol,order,rotated)
      integer npt,npol,order,ind
      real*8 points(npt,2), poles(npol,11), rotated(npt*npol,6)
! points(lat,long); poles (lat,long,rot,kappa,a-f,age),rotated(rtlat,rtlon,age,amaj,amin,az)
! order=1 for all points by each pole, order=2 for each point by every pole    
!f2py intent(in) npt,npol,order
!f2py intent(in) points,poles
!f2py intent(out)rotated
      if (order .eq. 2) then
        ind=0
        do 50000 i=1,npt
        do 50001 j=1,npol
        ind=ind+1
        call ptrotate(points(i,1),points(i,2),poles(j,1),poles(j,2),poles(j,3),poles(j,4),poles(j,5),poles(j,6),poles(j,7),poles(j,8),poles(j,9),poles(j,10),rotated(ind,1),rotated(ind,2),rotated(ind,4),rotated(ind,5),rotated(ind,6))
		rotated(ind,3)=poles(j,11)
50001   continue
50000   continue
      else if (order .eq. 1) then
        ind=0
        do 60000 j=1,npol
        do 60001 i=1,npt
        ind=ind+1
        call ptrotate(points(i,1),points(i,2),poles(j,1),poles(j,2),poles(j,3),poles(j,4),poles(j,5),poles(j,6),poles(j,7),poles(j,8),poles(j,9),poles(j,10),rotated(ind,1),rotated(ind,2),rotated(ind,4),rotated(ind,5),rotated(ind,6))
		rotated(ind,3)=poles(j,11)
60001   continue
60000   continue
      endif
      return
      end
      subroutine ptrotate(rlat,rlon,polat,polon,omega,rkaphat,a,b,c,d,e,f,rtlat,rtlon,amaj,amin,az)
      real*8 a,b,c,d,e,f
      real*8 n_disp,rlat,rlon,polat,polon,omega,rkaphat,rtlat,rtlon,amaj,amin,az
      real*8 r(3,3),w(3,3),allmat(3,3),simmat(3,3),ell_ned(3,3)
      real*8 eval(2),covt(3,3),rotmat(3,3),a1t(3,3),rtmt(3,3)
      real*8 rtmttp(3,3),tmp(3,3),ell_mat(3,3),ev(3,2)
      real*8 eig_val(3),eig_vec(3,3),evned(3,2)
      character*30 ians2*1
      parameter (pi=3.14159265)
	  parameter (raddeg =  180. / pi )
	  parameter (degrad = pi / 180. )
!f2py intent(in) rlat,rlon,polat,polon,omega,rkaphat,a,b,c,d,e,f
!f2py intent(out) rtlat,rtlon,amaj,amin,az
      do23000 i = 1,3 
      do23002 j = 1,3 
      r(i,j) = 0.00d0
      allmat(i,j) = 0.00d0
      ell_ned(i,j) = 0.00d0
23002 continue
23000 continue
	  ians=2
	  ians2='n'
	  covt(1,1)=a/rkaphat
	  covt(1,2)=b/rkaphat
	  covt(1,3)=c/rkaphat
	  covt(2,1)=b/rkaphat
	  covt(2,2)=d/rkaphat
	  covt(2,3)=e/rkaphat
	  covt(3,1)=c/rkaphat
	  covt(3,2)=e/rkaphat
	  covt(3,3)=f/rkaphat
      if(.not.(polon .lt. 0.0)) goto 23004
      polon = polon + 360.0
23004 continue
      call sphcar(polat*degrad,polon*degrad,unit_wx,unit_wy,unit_wz)
!      print*, rlat,rlon,polat,polon,omega,unit_wx,unit_wy,unit_wz
      call dmatrix(unit_wx,unit_wy,unit_wz,omega*degrad,rotmat)
      wx = unit_wx * omega * degrad
      wy = unit_wy * omega * degrad
      wz = unit_wz * omega * degrad
      do23006 i = 1,3 
      do23008 j = 1,3
      rtmttp(j,i) = rotmat(i,j)
23008 continue
23006 continue
! In original code a whole rescaling thing related to text prompts in original. w prob redundant now.
      do23020 i = 1,3 
      do23022 j = 1,3
      w(i,j) = covt(i,j)
23022 continue
23020 continue
      if(.not.(rlon .lt. 0.0))goto 23031
      rlon = rlon + 360.0
23031 continue
      call rotp(polat*degrad,polon*degrad,omega*degrad,rlat*degrad,rlon*degrad,rtlat,rtlon)
      rtlon=rtlon*raddeg
      rtlat=rtlat*raddeg
      if(.not.(rtlon .lt. 0.0))goto 23033
      rtlon = rtlon + 360.0
23033 continue
      call sphcar(rlat*degrad,rlon*degrad,x,y,z)
      call sphcar(rtlat*degrad,rtlon*degrad,xr,yr,zr)
      n_disp = 111.13*(rtlat-rlat)
      e_disp = (rtlon - rlon)/degrad
      disp = 111.13*sqrt((xr-x)**2 + (yr-y)**2 + (zr-z)**2)/degrad
      dd_wrt_wx = 111.13*(wx*z*z + wx*y*y - wy*x*y - wz*x*z)/(degrad*disp)
      dd_wrt_wy = 111.13*(wy*x*x + wy*z*z - wx*x*y - wz*y*z)/(degrad*disp)
      dd_wrt_wz = 111.13*(wz*y*y + wz*x*x - wx*x*z - wy*y*z)/(degrad*disp)
      sig_xx = w(1,1) * (dd_wrt_wx)**2
      sig_yy = w(2,2) * (dd_wrt_wy)**2
      sig_zz = w(3,3) * (dd_wrt_wz)**2
      cov_xy = w(1,2) * dd_wrt_wx * dd_wrt_wy
      cov_xz = w(1,3) * dd_wrt_wx * dd_wrt_wz
      cov_yz = w(2,3) * dd_wrt_wy * dd_wrt_wz
      sig_disp = 111.13*sqrt(sig_xx + sig_yy + sig_zz + 2.0d0*(cov_xy+cov_xz+cov_yz))/degrad
      r(1,2) = z 
      r(2,1) = -z
      r(1,3) = -y 
      r(3,1) = y
      r(3,2) = -x 
      r(2,3) = x
      call dmatmul(r,3,3,w,3,3,allmat)
      r(1,2) = -z 
      r(2,1) = z
      r(1,3) = y 
      r(3,1) = -y
      r(3,2) = x 
      r(2,3) = -x
      call dmatmul(allmat,3,3,r,3,3,a1t)
      call dmatmul(rotmat,3,3,a1t,3,3,tmp)
      call dmatmul(tmp,3,3,rtmttp,3,3,ell_mat)
      simmat(1,1) = -sin(rtlat*degrad)*cos(rtlon*degrad)
      simmat(1,2) = -sin(rtlat*degrad)*sin(rtlon*degrad)
      simmat(1,3) = cos(rtlat*degrad)
      simmat(2,1) = -sin(rtlon*degrad)
      simmat(2,2) = cos(rtlon*degrad)
      simmat(2,3) = 0.0d0
      simmat(3,1) = -cos(rtlat*degrad)*cos(rtlon*degrad)
      simmat(3,2) = -cos(rtlat*degrad)*sin(rtlon*degrad)
      simmat(3,3) = -sin(rtlat*degrad)
      call dmatmul(simmat,3,3,ell_mat,3,3,tmp)
      do23035 i = 1,3 
      do23037 j = 1,3
      rtmt(i,j) = simmat(j,i)
23037 continue
23035 continue      
      call dmatmul(tmp,3,3,rtmt,3,3,ell_ned)
      call djacobi(ell_ned,3,3,eig_val,eig_vec,nrot)
      ev(1,1) = eig_vec(1,3)
      ev(2,1) = eig_vec(2,3)
      ev(1,2) = eig_vec(1,2) 
      ev(2,2) = eig_vec(2,2)
      eval(1) = eig_val(3) 
      eval(2) = eig_val(2)
      ax1 = sqrt(eval(1))/degrad
      ax2 = sqrt(eval(2))/degrad
      az = atan2(ev(2,1),ev(1,1))/degrad
      amaj = dsqrt(2.0d0) * ax1
      amin = dsqrt(2.0d0) * ax2
      if(.not.(ians .eq. 2))goto 23039
      amaj = amaj * dsqrt(3.0d0)
      amin = amin * dsqrt(3.0d0)
23039 continue
!      print*, rtlat,rtlon,amaj,amin,az
      return
      end

      subroutine make_ell(amaj,amin,az,alat,alon,elat,elon)
      real*8 amaj,amin,az,alat,alon,vlat,vlon
      real*8 xlat,xlon,plat,plon,temp1,temp2,zero,ytenin
      real*8 elat(101),elon(101)  
!f2py intent(in) amaj,amin,az,alat,alon
!f2py intent(out) elat,elon
      parameter (pi=3.14159265)
	  parameter (raddeg =  180. / pi )
	  parameter (degrad = pi / 180. )
      ang = (alon-az)*degrad
      do 23041 loop=1,100 
      x = (amaj*degrad)*cos(2.0*pi*(dble(loop-1)/100.0))
      y = (amin*degrad)*sin(2.0*pi*(dble(loop-1)/100.0))
      z = sqrt(1.-x**2-y**2)
      call carsph(x,y,z,vlat,vlon)
      temp1 = (alon+90.d0)*degrad
      temp2 = (90.d0-alat)*degrad
      zero = 0.d0
      ytenin = 0.5*pi
      call rotp(ytenin,zero,ang,vlat,vlon,xlat,xlon)
      call rotp(zero,temp1,temp2,xlat,xlon,plat,plon)
      if(.not.(plon .lt. 0.d0))goto 23043
      plon = plon+(pi*2)
23043 continue
      if(.not.(loop .eq. 1))goto 23045
      elat(101) = plat*raddeg 
      elon(101) = plon*raddeg
23045 continue
      elat(loop)=plat*raddeg
      elon(loop)=plon*raddeg
23041 continue
      return
      end
c finb2 routine from Pavel Doubrovine's ibfrlib0.3.f other routines moved to 
c RotKit-nonpython      
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
c-----------------------------------------------------------------------
c Interpolation between two consecutive finite rotations:
c
c       Input:
c       Euler parameters for rotations 1 & 2 (ep1 and ep2)
c       ep(3)=(lat,lon,rho)
c       lat - Euler pole latitude [degrees N]
c       lon - Euler pole latitude [egrees E]
c       rho - rotation angle (degrees CCW)
c
c       Covariance matrices for rotations 1 & 2 (scaled by kappa),
c       cov1(3,3) and cov2(3,3) [radian**2]
c
c       xi=(t-t1)/(t2-t1) - interpolation parameter. t1 and t2 are
c       the ages of bounding rotations; t is the age of interpolated
c       rotation, t1<t<t2.
c
c       Output:
c       Euler parameters (ep) and covariance matrix (cov) for the
c       interpolated rotation
c
c Uses subroutines: EP2R, MATT, MULMM, RPV, PHI, MATRIXL
c-----------------------------------------------------------------------
        SUBROUTINE finb2(ep1,cov1,ep2,cov2,xi,ep,cov)
        INTEGER i,j
        DOUBLE PRECISION ep1(3),cov1(3,3),ep2(3),cov2(3,3)
        DOUBLE PRECISION xi,ep(3),cov(3,3)
        DOUBLE PRECISION deg,ahat1(3,3),ahat2(3,3),s(3,3),ahat(3,3),
     &  ps(3),p(3),ps1(3),lats,lons,rhos,lat,lon,rho,a(3,3),ml(3,3),
     &  l1(3,3),l2(3,3),temp1(3,3),temp2(3,3),c1(3,3),c2(3,3)
        PARAMETER (deg=.0174532925199432958D+0)
!f2py intent(in) ep1,ep2,cov1,cov2,xi
!f2py intent(out) ep(3)
!f2py intent(out) cov(3,3)
c
c INTERPOLATE FINITE ROTATION
c
c Check if 2 rotations are the same (i.e. no motion between t1 and t2)
c -> if TRUE set ep=ep1 and go to uncertainties
c
        if(ep1(1).eq.ep2(1).and.ep1(2).eq.ep2(2).and.
     &     ep1(3).eq.ep2(3)) then
          do 11 i=1,3
            ep(i)=ep1(i)
11        continue
          goto 222
        endif
c
c Calculate rotation matrices A1, A2
        call ep2r(ep1(1)*deg,ep1(2)*deg,ep1(3)*deg,ahat1)
        call ep2r(ep2(1)*deg,ep2(2)*deg,ep2(3)*deg,ahat2)
c
c Calculate the stage rotation S=A2*A1**t
        call matt(ahat1,a,3,3)
        call mulmm(ahat2,a,s,3,3)
c
c Calculate rotation pseudovector (ps) for S (S=Phi(ps))
        call rpv(s,lats,lons,rhos,ps)
c
c Calculate the interpolated rotation A=S'*A1, S'=Phi(xi*ps)
        rhos=xi*rhos
        call phi(lats,lons,rhos,a)
        call mulmm(a,ahat1,ahat,3,3)
c
c Calculate the pole and angle for A
        call rpv(ahat,lat,lon,rho,p)
        ep(1)=lat/deg
        ep(2)=lon/deg
        ep(3)=rho/deg
c
c
222     continue
c
c
c UNCERTAINTY OF THE INTERPOLTED ROTATION
c
c
c Check if 2 poles and uncertainties are the same (i.e.,
c we used the same pole, extrapolated rotation angle for ep2,
c and chose to set cov2=cov1) -> if TRUE set cov=cov1 and return
c
        if(ep1(1).eq.ep2(1).and.ep1(2).eq.ep2(2).and.
     &     cov1(1,1).eq.cov2(1,1).and.cov1(1,2).eq.cov2(1,2).and.
     &     cov1(1,3).eq.cov2(1,3).and.cov1(2,2).eq.cov2(2,2).and.
     &     cov1(2,3).eq.cov2(2,3).and.cov1(3,3).eq.cov2(3,3)) then
          do 21 i=1,3
            do 22 j=1,3
              cov(i,j)=cov1(i,j)
22          continue
21        continue
          goto 333
        endif
c
c Calculate the rotation pseudovector for S'
        do 141 i=1,3
          ps1(i)=ps(i)*xi
141     continue
c
c Calculate matrix L and cov(h) for A
        call matrixl(ahat1,ps,ps1,ml)
c
        do 145 i=1,3
          do 146 j=1,3
            if(i.eq.j) then
              l1(i,j)=1.-xi*ml(i,j)
              l2(i,j)=xi*ml(i,j)
            else
              l1(i,j)=-xi*ml(i,j)
              l2(i,j)=xi*ml(i,j)
            endif
146       continue
145     continue
c
        call matt(l1,temp1,3,3)
        call mulmm(cov1,temp1,temp2,3,3)
        call mulmm(l1,temp2,c1,3,3)
c
        call matt(l2,temp1,3,3)
        call mulmm(cov2,temp1,temp2,3,3)
        call mulmm(l2,temp2,c2,3,3)
c
        do 147 i=1,3
          do 148 j=1,3
            cov(i,j)=c1(i,j)+c2(i,j)
148       continue
147     continue
c
333     return
        END
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
      parameter (pi=dacos(-1.d0))
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
      parameter (pi=dacos(-1.d0))
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
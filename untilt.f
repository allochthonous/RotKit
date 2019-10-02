! -*- f90 -*-
! File untilt.f
      subroutine dir2cart(dir,xyz)
      real*8 dir(2),xyz(3)
      parameter (pi=dacos(-1.d0))
      parameter (degrad =  pi / 180.)
!f2py intent(in) dir
!f2py intent(out) xyz
      xyz(1)=cos(dir(1)*degrad)*cos(dir(2)*degrad)
      xyz(2)=cos(dir(2)*degrad)*sin(dir(1)*degrad)
      xyz(3)=sin(dir(2)*degrad)
      return
      end

      subroutine cart2dir(xyz,dir)
      real*8 dir(2),xyz(3)
!f2py intent(in) xyz
!f2py intent(out) dir
      parameter (pi=dacos(-1.d0))
      parameter (raddeg =  180. / pi)
      dir(2) = asin(xyz(3))
      if(.not.(xyz(1).eq.0.d0))goto 23047
      if(.not.(xyz(2).gt.0.d0))goto 23049
      dir(1) = pi*0.5
23049 continue
      if(.not.(xyz(2).lt.0.d0))goto 23051
      dir(1) = pi*1.5
23051 continue
      goto 23048
23047 continue
      dir(1) = atan(xyz(2)/xyz(1))
      if(.not.(xyz(1).lt.0.d0))goto 23053
      dir(1) = dir(1)+pi
23053 continue
      if(.not.(xyz(1).gt.0.d0.and.xyz(2).lt.0.d0))goto 23055
      dir(1) = dir(1)+(2*pi)
23055 continue
23048 continue
      dir(1)=dir(1)*raddeg
      dir(2)=dir(2)*raddeg
      return
      end


      subroutine untilt(dir,bedstrike,beddip,rotdir)
! adapted from the python routine dotilt in PMagPy, by Lisa Tauxe
      real*8 dir(2),bedstrike,beddip
      real*8 sa,ca,cdp,sdp
      real*8 newxyz(3)
      real*8 xyz(3),rotdir(2)
!f2py intent(in) dir
!f2py intent(in) bedstrike,beddip
!f2py intent(out)rotdir
      parameter (pi=3.14159265)
	  parameter (degrad =  pi / 180.)
      call dir2cart(dir,xyz)
!      print*, xyz
! get sines and cosines of new coordinate system
      sa=-sin((bedstrike+90.)*degrad)
      ca=cos((bedstrike+90.)*degrad)
      cdp=cos(beddip*degrad)
      sdp=sin(beddip*degrad)
! do the rotation
      newxyz(1)=xyz(1)*(sa*sa+ca*ca*cdp)+xyz(2)*(ca*sa*(1.-cdp))+xyz(3)*sdp*ca
      newxyz(2)=xyz(1)*ca*sa*(1.-cdp)+xyz(2)*(ca*ca+sa*sa*cdp)-xyz(3)*sa*sdp
      newxyz(3)=-(xyz(1)*ca*sdp-xyz(2)*sdp*sa-xyz(3)*cdp)
! convert from cartesian back to direction
!      print*, newxyz
      call cart2dir(newxyz,rotdir)
      return
      end
      
      
      subroutine untilt_dirs(dirs,n,bedstrike,beddip,rotdirs)
      integer n
      real*8 dirs(n,2),rotdirs(n,2),dir(2),rotdir(2)
      real*8 bedstrike,beddip
!f2py intent(in) dirs, n
!f2py intent(in) bedstrike,beddip
!f2py intent(out)rotdirs
      do 50000 i=1,n
      do 50002 j=1,2
      dir(j)=dirs(i,j)
50002 continue       
      call untilt(dir,bedstrike,beddip,rotdir) 
      do 50001 j=1,2
      rotdirs(i,j)=rotdir(j)
50001 continue       
50000 continue
      return
      end
      
      subroutine fshdev(k,dec,inc)
      real*8 k, dec, inc
      real*8 R1, R2, L, a, fac
      parameter (pi=dacos(-1.d0))
      parameter (raddeg =  180. / pi)
!f2py intent(in) k
!f2py intent(out) dec,inc
      R1=rand()
      R2=rand()
      L=exp(-2.*k)
      a=R1*(1-L)+L
      fac=sqrt((-log(a))/(2*k))
      inc=90.-2*asin(fac)*raddeg
      dec=2*pi*R2*raddeg
      return
      end

      subroutine get_fish(tdec,tinc,k,n,dirs)
      integer n
      real*8 tdec,tinc,k
      real*8 rotmat(3,3), dirs(n,2)
      real*8 dir(2),xyz(3),result(3)
      real*8 dec,inc
!f2py intent(in) tdec,tinc,k,n
!f2py intent(out) dirs
! Is probably not the most elegant way to make a rotation matrix...
C     dir(1)=tdec
C     dir(2)=90.-tinc
C     call dir2cart(dir,xyz)
C     rotmat(1,1)=xyz(1)
C     rotmat(1,2)=xyz(2)
C     rotmat(1,3)=xyz(3)
C     dir(1)=tdec-90.
C     dir(2)=0.
C     call dir2cart(dir,xyz)	
C     rotmat(2,1)=xyz(1)
C     rotmat(2,2)=xyz(2)
C     rotmat(2,3)=xyz(3)     
C     dir(1)=tdec
C     dir(2)=tinc
C     call dir2cart(dir,xyz)
C     rotmat(3,1)=xyz(1)
C     rotmat(3,2)=xyz(2)
C     rotmat(3,3)=xyz(3)
      do 60000 i=1,n
      call fshdev(k,dec,inc)
      dir(1)=dec
      dir(2)=inc
      call untilt(dir,tdec-90,90-tinc,dir)
C     call dir2cart(dir,xyz)
C     result(1)=xyz(1)*rotmat(1,1)+xyz(2)*rotmat(2,1)+xyz(3)*rotmat(3,1)
C     result(2)=xyz(1)*rotmat(1,2)+xyz(2)*rotmat(2,2)+xyz(3)*rotmat(3,2)
C     result(3)=xyz(1)*rotmat(1,3)+xyz(2)*rotmat(2,3)+xyz(3)*rotmat(3,3)
C     call cart2dir(result,dir)
      dirs(i,1)=dir(1)
      dirs(i,2)=dir(2)
60000 continue      
      return
      end
               
           
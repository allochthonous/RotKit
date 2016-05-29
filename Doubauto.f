        PROGRAM Doubauto
c
        DOUBLE PRECISION t1,ep1(3),kappa1,cov1(3,3),t2,ep2(3),kappa2,
     &                   cov2(3,3),t,xi,ep(3),cov(3,3)
        INTEGER i,j
        CHARACTER(3) rpt
c
1       write(*,'(a)') 'Input rotation 1 [age,EPlat,EPlon,angle]:'
        read(*,*,err=1) t1,ep1(1),ep1(2),ep1(3)
2       write(*,'(a)') 'Input covariance 1 [k,c11,c12,c13,c22,c23,c33]:'
        read(*,*,err=2) kappa1,cov1(1,1),cov1(1,2),cov1(1,3),cov1(2,2),
     &                  cov1(2,3),cov1(3,3)
        cov1(2,1)=cov1(1,2)
        cov1(3,1)=cov1(1,3)
        cov1(3,2)=cov1(2,3)
        do 3 i=1,3
          do 4 j=1,3
             cov1(i,j)=cov1(i,j)/kappa1
4         continue
3       continue
c
11      write(*,'(a)') 'Input rotation 2 [age,EPlat,EPlon,angle]: '
        read(*,*,err=11) t2,ep2(1),ep2(2),ep2(3)
        if(t2.le.t1) then
          write(*,*) 't2 must be > t1'
          goto 11
        endif
12      write(*,'(a)') 'Input covariance 2 [k,c11,c12,c13,c22,c23,c33]:'
        read(*,*,err=12) kappa2,cov2(1,1),cov2(1,2),cov2(1,3),cov2(2,2),
     &                   cov2(2,3),cov2(3,3)
        cov2(2,1)=cov2(1,2)
        cov2(3,1)=cov2(1,3)
        cov2(3,2)=cov2(2,3)
        do 13 i=1,3
          do 14 j=1,3
             cov2(i,j)=cov2(i,j)/kappa2
14        continue
13      continue
c
21      write(*,'(a)') 'Input intepolation age:'
        read(*,*,err=21) t
        if(t.lt.t1.or.t.gt.t2) then
          write(*,'(a)') 't must be between t1 and t2'
          goto 21
        endif
        xi=(t-t1)/(t2-t1)
        call finb2(ep1,cov1,ep2,cov2,xi,ep,cov)
c
c        write(*,*)
c        write(*,'(a,f7.3,a)') 'Rotation 1, ',t1,' Ma:'
c        write(*,'(3f18.12)') (ep1(i),i=1,3)
c        write(*,'(a)')'Covariance 1:'
c        do 31 i=1,3
c          write(*,'(3e21.12)') (cov1(i,j),j=1,3)
c31      continue
c        write(*,*)
c        write(*,'(a,f7.3,a)') 'Rotation 2, ',t2,' Ma:'
c        write(*,'(3f18.12)') (ep2(i),i=1,3)
c        write(*,'(a)')'Covariance 2:'
c        do 32 i=1,3
c          write(*,'(3e21.12)') (cov2(i,j),j=1,3)
c32      continue
c        write(*,*)
c        write(*,'(a,f7.3,a)')'Interpolated rotation, ',t,' Ma:'
        write(*,'(f7.3, 3f9.3, a, 6ES12.3E2 )') t, (ep(i),i=1,3), '  1.00', (cov(1,i),i=1,3), cov(2,2), cov(2,3), cov(3,3)
c        write(*,'(a)')'Covariance:'
c        do 33 i=1,3
c          write(*,'(3e21.12)') (cov(i,j),j=1,3)
c33      continue
cc
c        write(*,*)
c        write(*,'(a)') 'Another interpolation? [yes/no]:'
c        read(*,'(a)') rpt
c        if(rpt.eq.'yes') goto 1
c        stop
        END 

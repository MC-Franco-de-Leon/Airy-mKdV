ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   THIS FILE CONTAINS THE VELOCITY EVALUATION FOR H-SOF LAMBDA'S
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine lambda(m,h,p2,dcut,theta,dkap,un,utan
     &,sl,sigma) 
        parameter (maxn=4800,maxnintf=1)
        implicit double precision (a-h,o-z)
        dimension theta(0:4800),un(0:4800),utan(0:4800)
        dimension wt(0:4800),dkap(0:4800),tempi(0:4800)
        dimension stheta(0:4800),sigma(0:4800),sutan(0:4800)
        dimension rlambdal(0:4800)
c        common /essen/ ag,tau
c        common /init/ x0,y0
c        common /initvel/unorm
c        common /blob/dell

c      subroutine intsolve(nintf,m,h,p2,x,y,sigma,romega,
c     &                     dcut,tol,tol1,theta)
      
       dimension thetad(0:4800)
       dimension  vel(0:4800)

c*******************************************************************
c  gmres stuff

       parameter(nmax=maxn*4*maxnintf, liwork=20, maxl=50)
       parameter(lrwork=1+nmax*(maxl+6)+maxl*(maxl+3))
       dimension rhs(nmax), soln(nmax)
       dimension iwork(liwork), rwork(lrwork)
       dimension sb(nmax),sx(nmax),ia(1),ja(1),a(1)
       common /maxiter/itermax

c********************** more common stuff *************************

      common /par/hd,ddcut,p2d,md
      common /space/thetad,sld
      common /precon/ mpre
      common /iterations/iter
      external matvec_l,msolve_l

c****************** generate common stuff ************************

       mpre=5
       tol1=1.d-10
       tol=1.d-10
       dcut=1.d-10
       mtot=0
       h1=p2/dble(m)
       sld=sl

       p2d=p2
       nintfd=nintf
       fluxd=flux
       ddcut=dcut

       hd=h
       md=m
       do j=0,m
        thetad(j)=theta(j)
        stheta(j)=dkap(j)*sl
       end do

       call righthandside(un,utan,rhs)


c total number of unknowns (since periodic), extra one for lambdaA

       n=m+1

c make first guess.

       indx=0 
       do j=0,m       
        indx=indx+1
        soln(indx)=0.d0
       end do

c  Solve linear system using GMRES.
 
        nelt=1
        itol = 0
        itermax = 10000
        do i=1,20
           iwork(i) = 0
        end do
        iwork(1) = maxl
        iwork(4) = 1

c       iunit=6  write message from mres
c       iunit=0  don't write

        iunit=0

        call dgmres(n, rhs, soln, nelt, ia, ja, a, isym, matvec_l,
     &             msolve_l, itol, tol1, itmax, iter, err, ierr, iunit,
     &             sb,sx, rwork, lrwork, iwork, liwork, rrwork, iiwork)
 
        if (ierr .ne. 0 .or. iter .gt. itermax-1) then
          write(8,*)';** error diff GMRES, ierr=',ierr,' iter=',iter
          write(6,*)' ** error diff GMRES, ierr=',ierr,' iter=',iter
          return
        end if


c  now write solutions out into sigma.

        indx=0    
        do j=0,m-1
         indx=indx+1
         sigma(j)=soln(indx)
        end do
        sigma(m)=sigma(0)
        sigma(m+1)=soln(m+1)
c        romega=soln(m+1)
cc        write(31,*)romega

CC double check the solutionof lambda's
c        call eq1(theta,sl,sigma,vel)
c        do j=0,m-1
c           write(93,*)sl,theta(j),vel(j)-rhs(j+1),vel(m+1)-rhs(m+1)
c        end do
c        stop
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine righthandside(un,utan,rhs)
        implicit double precision (a-h,o-z)

        parameter (maxn=4800,maxnintf=1)
        parameter(nmax=maxn*4*maxnintf)
        dimension rhs(nmax)
        dimension theta(0:4800),un(0:4800),utan(0:4800)
        dimension stheta(0:4800),sutan(0:4800),tempi(0:4800)
        dimension salpha2(0:4800),sx(0:4800),sy(0:4800)
        dimension x(0:4800),y(0:4800)

        common /par/h,dcut,p2,m
        common /space/theta,sl
        common /area/aa0

       call fd1x(m,theta,stheta,h,dcut)       

c T perodic 

       call fd1(m,utan,sutan,h,dcut)


c consider robust contraints

      call recon(m,h,p2,sl,theta,x,y,avx,avy,dcut)
      call fd1(m,x,sx,h,dcut)
      call fd1(m,y,sy,h,dcut)

      do j=0,m-1
        salpha2(j)=(sx(j)**2.0+sy(j)**2.0)
      end do  

      aa=0.d0
      do i=0,m-1   
         aa=aa+(sin(theta(i))*x(i)-cos(theta(i))*y(i)
     $ )/dble(2.0)*sl/dble(m)
      end do


c form the RHS of the two equations
c now add them together, pay attention to the sign before each item.

c      rll=0.d0
c      rla=-0.d0

      rll=10.d0
      rla=-10.d0
      

       indx=0
       do j=0,m-1
         indx=indx+1
         rhs(indx)=-sutan(j)/sl-stheta(j)/sl*un(j)
     $-rll/2.d0*(1.d0-sl*sl/salpha2(j))
       end do

        call fin(m,h,p2,tempi,un)
        rhs(m+1)=-sl*tempi(m)
     $-rla/2.d0*(aa-aa0*aa0/aa)
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  THIS FILE CONTAINS MATVEC AND MSOLVE FOR THE 1ST
c  KIND INTEGRAL EQUATION.
c
c  this is in special form for the eigenvalue solve.
c
c  NOW DONE DIRECTLY IN TERMS OF LOG
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  THIS ROUTINE DOES THE DIAGONAL PRECONDITIONING.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine msolve_l(n, r, z)

      parameter (maxn=4800)
      implicit double precision (a-h,o-z)
 
      dimension r(1), z(1)

c      dimension x(0:maxn,maxnintf),y(0:maxn,maxnintf)
      dimension vel(0:maxn),stheta(0:maxn)
      dimension sigma(0:maxn),risigma(0:maxn)
c      dimension smean(maxnintf),vel1(0:maxn,maxnintf)     
      dimension ssigma(0:maxn), theta(0:maxn)
 

      common /par/h,dcut,p2,m
      common /space/theta,sl
      common /precon/ mpre
      common /essen/ ag,tau,epslon0,rmiu,eps
c      common /ffdot/ f,fdot,flux,area
      common /theta_smax/rmin,rmax


c do diagonal preconditioning if mpre=0

      mpre=5

       if (mpre .eq. 0)then
           do i=1,n
             z(i)=r(i)
           end do
       else

c do FFT transform preconditioning.

         call fd1x(m,theta,stheta,h,dcut) 

         rmin=stheta(1)
         rmax=stheta(1)
         do j=0,m-1
            rmin=min(stheta(j),rmin)
            rmax=max(stheta(j),rmax)
         end do
         
         indx=0
         do j=0,m-1
          indx=indx+1
          ssigma(j)=r(indx)
         end do
         ssigma(m)=ssigma(0)

         call fast(ssigma,m)

         do j=1,m/2+1
           k=2*(j-1)
           ssigma(k)=ssigma(k)/(dble(j-1)*dble(j-1)+rmax*rmax)
           ssigma(k+1)=ssigma(k+1)/(dble(j-1)*dble(j-1)+rmax*rmax)
         end do

         z(m+1)=(r(m+1)-rmin*ssigma(0))/sl
         call fsst(ssigma,m)
         indx=0
         do j=0,m-1
           indx=indx+1
           z(indx)=ssigma(j)
         end do
       end if
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  THIS DOES MATRIX MULTIPLICATION
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine matvec_l(n, xx, yy)

      implicit double precision (a-h,o-z)
      parameter (maxn=4800,maxnintf=1)

      dimension xx(1), yy(1)

c      dimension x(0:maxn,maxnintf),y(0:maxn,maxnintf)
c      dimension h(maxnintf)
c      dimension m(maxnintf)

      dimension vel(0:maxn),sigma(0:maxn),theta(0:4800)
 
      common /par/h,dcut,p2,m
      common /space/theta,sl
c      common /essen/ ag,tau,epslon0,rmiu,eps
c      common /precon/ mpre
c      common /ffdot/ f,fdot,flux,area
   
      indx=0
c      bsigma=0.d0
      do j=0,m-1
       indx=indx+1 
       sigma(j)=xx(indx)
c       bsigma=bsigma+sigma(j)
      end do

      sigma(m)=sigma(0)
      sigma(m+1)=xx(m+1)

c      romega=xx(m+1)
c      bsigma=bsigma/dble(m)

c now compute various integrals

       call eq1(theta,sl,sigma,vel)

c now put result into right place.

      indx=0

      do j=0,m-1
       indx=indx+1
       yy(indx)=vel(j)
      end do
      yy(m+1)=vel(m+1)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c THIS CODE COMPUTES THE SOLUTION TO THE HELE-SHAW PROBLEM WITH
c SURFACE TENSION AND DENSITY STRATIFICATION
c
c THIS CODE USES EXPLICIT ADAMS-BASHFORTH FOR THE SL EQUATION AND
c HAS AN INTEGRATING FACTOR FOR THE THETA EQUATION. BOTH ARE IN THE
c EQUAL ARCLENGTH FRAME. IT CALLS:
c                                   hsforce.f
c                                   libhs.f
c
c  THIS CODE COMPUTES THE CLOSED INTERFACE PROBLEM WITH PERIODIC
c  DATA ON [0,1]
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
      subroutine Principal(IniPert,deltafixt,t,N,dt,ncut,tp,eps,
     1  mode,willmore)    
      implicit double precision(a-h,o-z),integer(i-n)
      dimension xn(0:4800),yn(0:4800)
      dimension theta(0:4800),thetan(0:4800),thetann(0:4800)
      dimension alpha(0:4800),dalpha(0:4800)
      dimension fnthetan(0:4800),ftheta(0:4800),fntheta(0:4800)
      dimension temp1(0:4800),temp2(0:4800),temp3(0:4800),temp4(0:4800)
      integer N,ncut,nstep
      common /essen/ ag,tau
      common /init/ x0,y0
      common  a0,sl0

      common /initvel/unorm
      common /symm/aniso,flux,tcomp,wc,md
      common /blob/del


      open(33,file='Deltap.m')
      open(32,file='R1p.m')
      open(31,file='rPert1p.m')
      open(30,file='r0p.m')
      open(29,file='deloradp.m')
      open(28,file='areap.m')
      open(1,file='xp.m')
      open(2,file='yp.m')
      open(10,file='curp.m')
      open(11,file='tp.m')
c      open(13,file='wnorm.m')
      open(14,file='slp.m')
      open(15,file='thetap.m')
      open(55,file='temp.m')

      print*,''
      read(5,*)t,N,dt,tau,ncut,tp,ag,eps,del,md,aniso,flux,wc

      write(1,*)'%  t,N,dt'
      write(1,*)'% ', t,N,dt
      write(1,*)'% tau,ncut,tp,ag ='
      write(1,*)'% ', tau,ncut,tp,ag

      del=0.d0
      dcut=10.d0**(-ncut)
      pi=4.d0*datan(1.d0)
      p2=2.d0*pi
      h=1.d0/dble(N)
      nstep=int(t/dt+0.00001d0)
cccc      eps=0.1d0
      tol=1.d-13
      N2=N/2
      m=N
      dto=dt
      ktp=int(tp/dt+0.00001d0)
      write(11,*)'t=['
      write(14,*)'sl=['
	write(29,*)'delorad=['
        write(30,*)'r0=['
        write(31,*)'rPert1=['
        write(32,*)'R1=['
   	write(33,*)'Delta=['
c     Initialization.
        
        ktime=0
      call initiall(N,h,p2,sl,theta,eps,dcut,tol,alpha,xn,yn,dalpha,
     1              Inipert)

	
       	   
	   a0=0.d0
      do i=0,m-1   
         a0=a0+(sin(theta(i))*xn(i)-cos(theta(i))*yn(i))
     $ /dble(2.0)*sl/dble(m)
       end do 
   	sl0 = sl

c  First step is done by forward Euler method.
     
      write(6,*)'time = ',real(dt),'   sl = ',sl 

      call forcsl(m,h,p2,del,dcut,sl,theta,dalpha,fsl)
      call forctheta(m,h,p2,del,dcut,sl,theta,dalpha,ftheta)

      sln=sl+dt*fsl


      do j=0,m-1
       thetan(j)=theta(j)+dt*ftheta(j)
      end do
c      thetan(m)=thetan(0)


c    update initial interface position using Euler.
 
        rncomp2=dcos(theta(0))
        rncomp2o=dcos(theta(0))
 
        rncomp1=-sin(theta(0))
        rncomp1o=-sin(theta(0))

        unormo=0.d0
 
        x0=x0+.5d0*dt*(2.d0*unorm*rncomp1-0*unormo*rncomp1o)
        y0=y0+.5d0*dt*(2.d0*unorm*rncomp2-0*unormo*rncomp2o)

        unormo=unorm


      call fthetaim(m,h,p2,del,dcut,sl,theta,dalpha,fntheta)

c  write results if the time is right
 
        if (mod(ktime,ktp) .eq. 0) then
          tcomp=dt*dble(ktime)
          np=ktime/ktp+1
          call writout(N,h,p2,np,tcomp,dcut,sl,theta,ktime,pi
     &              ,usum1,dt)
        endif

c  now begin time iterations

      do 100 ktime=2,nstep
      
        tcomp=dt*dble(ktime)

c first update sl by Adams-Bashforth.

        call forcsl(m,h,p2,del,dcut,sln,thetan,dalpha,fsln)
        write(6,*)'time=',real(tcomp),'sl =',sln
        slnn=sln+.5d0*dt*(3.d0*fsln-fsl)

c now update theta using an integrating factor

        call fthetaim(m,h,p2,del,dcut,sln,thetan,dalpha,fnthetan)
 
        do j=0,m-1
         temp1(j)=fnthetan(j)
         temp2(j)=thetan(j)-p2*dble(j)*h
         temp3(j)=fntheta(j)
        end do
        temp1(m)=temp1(0)
        temp2(m)=temp2(0)
        temp3(m)=temp3(0)

        call fast(temp1,m)
        call fast(temp2,m)
        call fast(temp3,m)

        do j=1,m/2+1

         k=2*(j-1)
         rk=p2*dble(j-1)




c Local model	 
         rsl=wc*(rk/sl)**4
         rslnn=wc*(rk/slnn)**4
         rsln=wc*(rk/sln)**4

         d1=.5d0*dt*(rsln+rslnn)
         d2=.5d0*dt*(rsl+rslnn)+dt*rsln
         
         d1=dexp(-d1)
         d2=dexp(-d2)

         temp4(k)=d1*temp2(k)+.5d0*dt*
     &                               (3.d0*d1*temp1(k)-d2*temp3(k))
         temp4(k+1)=d1*temp2(k+1)+.5d0*dt*
     &                           (3.d0*d1*temp1(k+1)-d2*temp3(k+1))

        end do

        call fsst(temp4,m)

        do j=0,m-1
         thetann(j)=temp4(j)+p2*dble(j)*h
        end do
        thetann(m)=thetann(0)+p2*dble(m)*h



c    update initial interface position using AB2.

        rncomp2=dcos(thetan(0))
        rncomp2o=dcos(theta(0))

        rncomp1=-sin(thetan(0))
        rncomp1o=-sin(theta(0))

        x0n=x0+.5d0*dt*(3.d0*unorm*rncomp1-unormo*rncomp1o)
        y0n=y0+.5d0*dt*(3.d0*unorm*rncomp2-unormo*rncomp2o)

c  update values

        x0=x0n
        y0=y0n

        unormo=unorm

        sl=sln
        sln=slnn
        fsl=fsln

        do j=0,m-1
          theta(j)=thetan(j)
          thetan(j)=thetann(j)
          fntheta(j)=fnthetan(j)
        end do
        theta(m)=theta(0)
        thetan(m)=thetan(0)
        fntheta(m)=fntheta(0)


c  write results if the time is right

        if (mod(ktime,ktp) .eq. 0) then
          tcomp=dt*dble(ktime)
          np=ktime/ktp+1
        call writout(N,h,p2,np,tcomp,dcut,sln,thetan,ktime,
     &        pi,usum1,dt)
ccc          call writout(N,h,p2,np,tcomp,dcut,sln,thetan)
        endif

	
 100  continue
      write(11,*)'];'
      write(14,*)'];'
	write(29,*)'];'
        write(30,*)'];'
	write(31,*)'];'
        write(32,*)'];'
	write(33,*)'];'
c      stop
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c    THIS ROUTINE PRINTS THE OUTPUT
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine writout(N,h,p2,ktp,tcomp,dcut,sl,theta,ktime,pi
     &     ,usum1,dt,deltafixt)
      implicit double precision (a-h,o-z),integer(i-n)
      dimension theta(0:4800),dkap(0:4800),x(0:4800)
      dimension rPert1(0:4800),R1(0:4800)
      dimension y(0:4800),sdkap(0:4800),ytemp(0:4800)
      dimension wtemp(0:4800),w(0:4800), epslon(0:4800), term1(0:4800)
      integer N,i
      common /essen/ ag,tau
      common /symm/aniso,md

c  construct necessary relevant quantities


      call recon(N,h,p2,sl,theta,x,y,avx,avy,dcut)
      call acurv(N,sl,theta,dkap,dcut,h)
      call fd1(N,dkap,sdkap,h,dcut)

c****************calculating perturbation delta
	do i=0,N-1
		rPert1(i)=dsqrt(x(i)**2+y(i)**2)
	end do	
	call fin(N,h,p2,R1,rPert1)
	
	do i=0,N-1
		Delta=max(abs(rPert1(i)-R1(N)),Delta)
	end do
c*************this is delta for specific time 25
        if (tcomp .eq. 25) then
		deltafixt=Delta
	endif

c*****************************
c Change here    
      do j=0,N-1
       w(j)=tau*sdkap(j)+sl*dsin(theta(j))
       
        epslon(j)=1.d0
     &-(dble(md*md)-1.d0)*rmiu*
     &dcos(dble(md)*(theta(j)-theta0))
       term1(j)=-epslon(j)*dkap(j)
      end do
      w(N)=w(0)


      do j=0,N-1
        ytemp(j)=y(j)
        wtemp(j)=w(j)
      end do
      ytemp(N)=ytemp(0)
      wtemp(N)=wtemp(0)
       
      
      call vmax(N,dkap,curmax)
      curinv=1.d0/curmax
      

 10   continue
      call fast(ytemp,N)
      call fast(wtemp,N)
	
      write(1,*)'% t=',tcomp
      write(1,*)'x(:,',ktp,')=['
      write(2,*)'% t=',tcomp
      write(2,*)'y(:,',ktp,')=['
c      write(3,*)'% t=',tcomp
c      write(3,*)'w(:,',ktp,')=['
      write(10,*)'% t=',tcomp
      write(10,*)'cur(:,',ktp,')=['

c      write(13,*)'% t=',tcomp
c      write(13,*)'wnorm(:,',ktp,')=['
      write(15,*)'% t=',tcomp
      write(15,*)'theta(:,',ktp,')=['
      
      write(55,*)'% t=', tcomp
      write(55,*)'term1(:,', ktp,')=['
     
 
      write(11,*)tcomp
      write(14,*)sl


      do 12 k=0,N,2
      wt1=dsqrt((wtemp(k)/dble(N))**2+(wtemp(k+1)/dble(N))**2)
      yt1=dsqrt((ytemp(k)/dble(N))**2+(ytemp(k+1)/dble(N))**2)
      wt2=log(wt1+1.d-20)/log(10.d0)
      yt2=log(yt1+1.d-20)/log(10.d0)
 12   continue


c************* Figure out cnetroid and output (x,y) ***********

      sum2=dble(0)
      sumxx=dble(0)
      sumyy=dble(0)
      x(N)=x(0)
      y(N)=y(0)

      do 13 i=0,N-1   
         sum2=sum2+(sin(theta(i))*x(i)-cos(theta(i))*y(i)
     $ )/dble(2.0)*sl/N
         sumxx=sumxx+(sin(theta(i))*x(i)**2-cos(theta(i))*0.d0
     $ )/dble(2.0)*sl/N
         sumyy=sumyy+(sin(theta(i))*0.d0-cos(theta(i))*y(i)**2
     $ )/dble(2.0)*sl/N         
c         write(3,*)w(i)
         write(10,*)dkap(i)
c         write(13,*)w(i)/sl
         write(15,*)theta(i)
 13   continue

c************ compute centroid *******************

      xcc=sumxx/sum2
      ycc=sumyy/sum2


c********** calculate max( \delta/R) and ouput (x,y) ********

      r0=dsqrt(sum2/pi)
      delmax1=dble(0)
      delmax2=dble(0)
 
      do 14 i=0,N-1
         dx=x(i)-xcc
         dy=y(i)-ycc
         write(1,*)x(i)-xcc
         write(2,*)y(i)-ycc
         write(55,*)term1(i)
         delmax1=min((dsqrt(dx*dx+dy*dy)-r0),delmax1)
         delmax2=max((dsqrt(dx*dx+dy*dy)-r0),delmax2)
	 
14     continue
     
       write(1,*)x(N)-xcc
       write(2,*)y(N)-ycc

      delmax=max(abs(delmax1),abs(delmax2))
      delorad=delmax/r0
	write(90,*),delmax
c********* write out everything we need ******************

        write(29,*)delorad
        write(30,*)r0
  	write(31,*)rPert1
	write(32,*)R1
	write(33,*)Delta
        write(28,*)sum2

       write(1,*)'];'
       write(2,*)'];'
       write(3,*)'];'
       write(10,*)'];'
c       write(13,*)'];'
       write(15,*)'];'
       write(55,*)'];'
       
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



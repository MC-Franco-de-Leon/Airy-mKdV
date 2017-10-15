ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Airy mKdV. 10/15/17  
c
c  THIS CONTAINS THE LIBRARY OF ROUTINES FOR THE COMPUTATION
c  OF THE TANGENT ANGLE FORMULATION OF SOME FREE-SURFACE FLOWS.
c  PARTICULARLY FOR AIRY FLOW. THIS ASSUMES THAT A CLOSED
c  CURVE IS GIVEN AS INITIAL DATA. THIS IS PERIODIC ON [0,1].
c  THIS FILE CONTAINS:
c
c      arc       - computes the equal arclength points.
c      recon     - reconstructs the x,y positions from theta,sl
c      fin       - computes the indefinite integral of a per. function
c      kfilter   - filters a periodic function
c                      smoothing is used, last 2 modes set to 0
c      kfilterx  - filters x=alpha+s
c      kfiltern  - sets n/2 mode to 0
c      kfilternx - sets n/2 mode to 0 in x=alpha+s
c      fd1       - compute deriv of a per. function
c      acurv     - compute curvature from theta,sl
c      uset      - compute norm and tan vel. from ave vel, this calls
c                     recon, velocity, rntvel. This part is for vortex sheets,
c                     and must be modified for Hele-Shaw.
c      rntvel    - compute normal and tangential vel from ave vel.
c      rtvel     - compute tan. vel. for equal arclength frame
c      velocity  - compute ave. velocity
c      initiall  - initialize the solution vector
c      stfilter  - filter x so that x_alpha=1+s_alpha
c      tfilter   - not used, filters theta
c      vmax      - compute max of a vector
c      thetsolve - compute tan angle from curvature.
c
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   THIS SUBROUTINE, GIVEN A PARAMETRIZATION BETA(J)=j*h, GIVES A
c   SET OF EQUALLY SPACED POINTS IN ARCLENGTH
c   THIS IS DONE VIA NEWTON ITERATION
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine arc(m,h,ds,alpha,tol)
       implicit double precision (a-h,o-z)
       dimension alpha(0:4800),dsi(0:4800)
       dimension ds(0:4800),f(0:4800),fp(0:4800)
       dimension a(0:4800),b(0:4800),c(0:4800),d(0:4800)

       m1=m+1
       pi=4.d0*datan(1.d0)
       p2=2.d0*pi

c  integrate ds and set-up functions for interpolation and iteration
       call fin(m,h,p2,dsi,ds)
       ave=dsi(m)
       do j=0,m
         f(j)=dble(j)*h*ave-dsi(j)
         fp(j)=-ds(j)
       end do 

       do j=0,m
         a(j)=f(j)/dble(m)
         b(j)=0.d0
         c(j)=fp(j)/dble(m)
         d(j)=0.d0
       end do

       call fft842(0,m,a,b)
       call fft842(0,m,c,d)

c  for each point, iterate.
       alpha(0)=0.d0
       do j=1,m 

          ao=alpha(j-1)+h*(ave/ds(j))

c  begin iterations

          do l=0,100
             
             aoo=p2*ao
             call synth(f,v,aoo,m1,a,b)
             call synth(fp,vp,aoo,m1,c,d)
             v=v+(dble(j)*h-ao)*ave

             dalp=v/vp
             an=ao-dalp

             err=dabs(dalp)
c             write(6,*)dalp

             if (err .le. tol) goto 1000

             ao=an

          end do

 1000     continue
          alpha(j)=an

       end do
         
       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  GIVEN AN ARCLENGTH SL, AND A TAN ANGLE THETA, THIS CODE
c  RECONSTRUCTS THE CURVE AND GIVES X,Y
c
c  NEED EQUAL ARCLENGTH TO DO
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine recon (m,h,p2,sl,theta,x,y,avx,avy,dcut)
       implicit double precision (a-h,o-z)
       dimension theta(0:4800),x(0:4800),y(0:4800)
       dimension temp(0:4800),temp1(0:4800)
       common /init/ x0,y0
c       common /centroid/ xc,yc  

       do j=0,m-1
        temp(j)=sl*dcos(theta(j))
        temp1(j)=sl*dsin(theta(j))
       end do
      
       call fin(m,h,p2,x,temp)
       call fin(m,h,p2,y,temp1)

       avx=x(m)
       avy=y(m)
 
c       avx=1.d0
c       avy=0.d0

       do j=1,m
        x(j)=x(j)-dble(j)*h*avx+x0
        y(j)=y(j)-dble(j)*h*avy+y0
       end do
       x(0)=x0
       y(0)=y0

c       call kfilter(m,dcut,x)
c       call kfilter(m,dcut,y)
       
       avy=0.d0
       avx=1.d0

       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   THIS CODE FORMS THE INDEFINITE INTEGRAL OF A PERIODIC FUNCTION
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine fin (n,h,p2,yi,y)
       implicit double precision (a-h,o-z)
       dimension yi(0:4800),y(0:4800),a(0:4800)
       dimension b(0:4800)

       do j=0,n-1
        a(j)=y(j)
       end do
       a(n)=a(0)
 

       call fast(a,n)
       b(0)=0.d0
       b(1)=0.d0

       do j=1,n/2

        k=2*j
        b(k)=a(k+1)/(p2*dble(j))
        b(k+1)=-a(k)/(p2*dble(j))

       end do

       call fsst(b,n)
       b(n)=b(0)

       do j=0,n
        yi(j)=(a(0)/dble(n))*dble(j)*h + b(j)-b(n)
       end do

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine fdi(n,h,p2,dy,y)
       implicit double precision (a-h,o-z)
       dimension dy(0:4800),y(0:4800),a(0:4800)
       dimension b(0:4800)

       do j=0,n-1
        a(j)=y(j)
       end do
       a(n)=a(0)
 
         call fast(a,n)
         b(0)=0.d0
         b(1)=0.d0
         do  j=1,n/2
         k=2*j
         b(k)=-a(k+1)*dfloat(j)*dble(p2)
         b(k+1)=a(k)*dfloat(j)*dble(p2)
         end do  
         
         call fsst(b,n)
         b(n)=b(0)

       do i=0,n-1
        dy(i)=b(i)
       end do
       dy(n)=dy(0)

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine filterfdi(n,h,p2,dy,y,dcut)
       implicit double precision (a-h,o-z)
       dimension dy(0:4800),y(0:4800),a(0:4800)
       dimension b(0:4800),fb(0:4800)
       n2=n/2
       do j=0,n-1
        a(j)=y(j)
       end do
       a(n)=a(0)
 
         call fast(a,n)
         b(0)=0.d0
         b(1)=0.d0
         do  j=1,n2
         k=2*j
         b(k)=-a(k+1)*dfloat(j)*dble(p2)
         b(k+1)=a(k)*dfloat(j)*dble(p2)
         end do  
c         b(n)=0
c         b(n+1)=0
         call fsst(b,n)
         b(n)=b(0)
         call kfilter(n,h,p2,fb,b,dcut)
       do i=0,n-1
        dy(i)=fb(i)
       end do
       dy(n)=dy(0)

       return
       end
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   fourier filter
c
c   USING FAST.f FOURIER TRANSFORM
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine kfilter(n,h,p2,b,y,dcut)
       implicit double precision (a-h,o-z)
       dimension fy(0:4800),y(0:4800),a(0:4800)
       dimension b(0:4800)
       almbd=.5
       auxrho=dfloat(1)
       n2=n/2         
       do j=0,n-1
        a(j)=y(j)
       end do
       a(n)=a(0)
 
         call fast(a,n)
         b(0)=a(0)
         b(1)=a(1)
         do  j=0,n/2
         k=2*j
        
         
         
          avar=dble(j)/n2
          if (avar .le. almbd) then
             auxrho=1
          else if ((avar .gt. almbd) .and. (avar .le. 1)) then
             aux3=(1-((avar-almbd)/(1-almbd)))**4
             aquotient=exp(1-dfloat(1)/aux3)
             auxrho=aquotient          

c filter 0
c             aux3=1-((avar-almbd)/(1-almbd))**2
c             aquotient=exp(-dfloat(1)/aux3)
c             auxrho=exp(dfloat(1))*aquotient          
cfilter ceniceros
c              ax=1+((-avar+almbd)/(1-almbd))
c              auxrho=(35-84*ax+70*ax**2-20*ax**3)*ax**4

          end if
          
c         auxrho=1 

         b(k)=a(k)*auxrho
         b(k+1)=a(k+1)*auxrho
         
         
         test=(b(k)/dble(n))**2 + (b(k+1)/dble(n))**2
         test=dsqrt(test)
         auxtol=10.d0**(-9)
         if (test .le. auxtol) then
           b(k)=0.d0
           b(k+1)=0.d0
         end if
         
         end do  
         
         
         
         call fsst(b,n)
         b(n)=b(0)


       return
       end       

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c  subroutine to compute SPECTRAL derivatives in space
c  for periodic data
c
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine fd1(n,x,sx,h,dcut)
       implicit double precision(a-h,o-z)
       dimension x(0:4800),sx(0:4800),x1(0:4800),sx1(0:4800)
      

       do i=0,n-1
        x1(i)=x(i)
       end do
       x1(n)=x(0)

       call fdiff(x1,sx1,n+1,h)

       do i=0,n-1
        sx(i)=sx1(i)
       end do
       sx(n)=sx(0)



       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  remove linear part.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine fd1x(n,x,sx,h,dcut)
       implicit double precision(a-h,o-z)
       dimension x(0:4800),sx(0:4800),x1(0:4800),sx1(0:4800)

       pi=4.d0*datan(1.d0)
       p2=2.d0*pi
      

       do i=0,n-1
        x1(i)=x(i)-p2*dble(i)*h

       end do
       x1(n)=x(0)


       call fdiff(x1,sx1,n+1,h)

       do i=0,n-1
        sx(i)=sx1(i)+p2

       end do
       sx(n)=sx(0)


       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c THIS ROUTINE COMPUTES THE CURVATURE FROM THE TANGENT
c ANGLE FORMULATION
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine acurv(m,sl,theta,dkap,dcut,h)
       implicit double precision (a-h,o-z)
       dimension theta(0:4800),dkap(0:4800)
       dimension stheta(0:4800),b(0:4800)

       call fd1x(m,theta,stheta,h,dcut)
c       call kfilter(m,dcut,stheta)

       do j=0,m-1
        dkap(j)=stheta(j)/sl
       end do
       dkap(m)=dkap(0)

       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c THIS ROUTINE COMPUTES THE CURVATURE FROM THE TANGENT
c ANGLE FORMULATION
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine Dacurv(theta,dkap,m,sl,dcut,h)
       implicit double precision (a-h,o-z)
       dimension theta(0:4800),dkap(0:4800)
       dimension stheta(0:4800)
       dimension x1(0:4800),sx1(0:4800),sx(0:4800)
       pi=4.d0*datan(1.d0)
       p2=2.d0*pi
        
       do i=0,m-1
        x1(i)=theta(i)-p2*dble(i)*h
       end do
       x1(m)=theta(0)
       
       call filterfdi(m,h,p2,sx1,x1,dcut)
c       call fdiff(x1,sx1,m+1,h) 
       
       
       do i=0,m-1
        stheta(i)=sx1(i)+p2
       end do
       stheta(m)=stheta(0)       
       do j=0,m-1
        dkap(j)=stheta(j)/sl
       end do
       dkap(m)=dkap(0)
        
       return
       end




   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  THIS CODE COMPUTES THE NORMAL VELOCITY GIVEN THE U,V FROM THE
c  VELOCITY ROUTINE
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rntvel(m,h,p2,x,y,sl,theta,u,v,un,ut,dcut)
       implicit double precision (a-h,o-z)
       dimension x(0:4800),y(0:4800)
       dimension u(0:4800),v(0:4800),un(0:4800),ut(0:4800)
       dimension sx(0:4800),sy(0:4800),theta(0:4800)

c       call kfilter(m,dcut,u)
c       call kfilter(m,dcut,v)

C       write(22,*)'sx=, sy,sl',sl
       do j=0,m-1

        sx(j)=cos(theta(j))
        sy(j)=sin(theta(j))
C         write(22,*)sx(j),sy(j)
 
       end do


       do j=0,m-1

        un(j)=(sx(j)*v(j)-sy(j)*u(j))
        ut(j)=(sx(j)*u(j)+sy(j)*v(j))

       end do
       un(m)=un(0)
       ut(m)=ut(0)

c       call kfilter(m,dcut,un)
c       call kfilter(m,dcut,ut)
 
       return
       end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  THIS ROUTINE COMPUTES THE TANGENT VELOCITY GIVEN THE NORMAL
c  VELOCITY. IT ASSUMeS THAT THE TANGENTIAL VELOCITY AT ALPHA=0
c  IS 0
c
c  IT IS TO BE USED IN CONJUNCTION WITH THE TANGENT ANGLE FORMULATION
c
c     UNCONSTRAINTED
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine rtvel_un(m,h,p2,theta,unorm,utan,dcut)
       implicit double precision (a-h,o-z)
       dimension theta(0:4800),unorm(0:4800),utan(0:4800)
       dimension stheta(0:4800),temp(0:4800)

       call fd1x(m,theta,stheta,h,dcut)
c       call kfilter(m,dcut,stheta)

       do j=0,m-1
        temp(j)=stheta(j)*unorm(j)  
        utan(j)=0.d0
       end do
       utan(m)=0.d0
       temp(m)=temp(0)

       call fin(m,h,p2,utan,temp)

c original       do j=0,m-1
c        utan(j)=utan(j)-utan(m)*dble(j)*h
c       end do
c        utan(m)=utan(0)

c new tangential velocity, constant stuff
       
       do j=0,m-1
c        utan(j)=utan(j)-utan(m)*dble(j)*h
          utan(j)=0.d0
       end do
        utan(m)=utan(0)

c       call kfilter(m,dcut,utan)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  THIS ROUTINE COMPUTES THE TANGENT VELOCITY GIVEN THE NORMAL
c  VELOCITY. IT ASSUMeS THAT THE TANGENTIAL VELOCITY AT ALPHA=0
c  IS 0
c
c  IT IS TO BE USED IN CONJUNCTION WITH THE TANGENT ANGLE FORMULATION
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine rtvel(m,h,p2,theta,un,ut,dcut,sl)       

       implicit double precision (a-h,o-z)
       dimension theta(0:4800),un(0:4800),ut(0:4800)
       dimension stheta(0:4800),temp(0:4800) 
       dimension dkap(0:4800)


       call Dacurv(theta,dkap,m,sl,dcut,h)
       do j=0,m-1
         ut(j)=.5d0*dkap(j)**2
       end do
        ut(m)=ut(0)

       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c   THIS CODE FORMS THE NORMAL AND TAN VELOCITY GIVEN THE ANGLE AND ARCLENGTH 
c   and returns the normal and tan vel. from the ave. vel. integral
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine uset(m,h,p2,dcut,sl,theta,dalpha,w,un,ut)
       implicit double precision (a-h,o-z)
       dimension theta(0:4800),un(0:4800),ut(0:4800),x(0:4800)
       dimension wt(0:4800),dkap(0:4800),sdkap(0:4800)
       dimension temp2(0:4800),temp4(0:4800),dalpha(0:4800)
       dimension stheta(0:4800),sstheta(0:4800)
       double precision unorm,utan
       common /Nvel/ unorm
        common /length/salphainv,salpha

        
        do j=0,m-1
         temp2(j)=theta(j)-p2*dble(j)*h
        end do
        temp2(m)=temp2(0)
        call fast(temp2,m)
        
        do j=1,m/2
         k=2*j
         rk=p2*dble(j)
c Local model	 
         rsl=(rk/sl)**2
         
         temp4(k)=rsl*temp2(k)
         temp4(k+1)=rsl*temp2(k+1)
        end do
        call fsst(temp4,m)       

        do j=0,m-1
        un(j)=temp4(j) 
        end do
         un(m)=un(0)
         
         unorm=un(0) 
    
        return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  THIS SUBROUTINE GIVES INITIAL DATA IN X,Y.
c  THIS ROUTINE IS FOR CLOSED CURVES.
c
c ***************************************************
c    2 PLACES NEED TO BE MODIFIED TO CHANGE INITIAL CONDITIONS
c ***************************************************
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine initiall(m,h,p2,sl,theta,dcut,tol,alpha,xn,yn,
     &                         dalpha)
       implicit double precision (a-h,o-z)
       dimension slo(0:4800),dalpha(0:4800)
       dimension sx(0:4800),sy(0:4800),alpha(0:4800),theta(0:4800)
       dimension xn(0:4800),yn(0:4800),sxn(0:4800),syn(0:4800)
       dimension x(0:4800),y(0:4800)
       
       dimension area2(0:4800),areai(0:4800)
       dimension sumxx(0:4800),sumyy(0:4800)
       dimension sumxxi(0:4800),sumyyi(0:4800)


        double precision r0,rad


       common /init/ x0,y0
       common /centroid/ xc,yc  


       pi=4.d0*datan(1.d0)

       do j=0,m

c*************************************
c THIS IS 1st PLACE FOR INITIAL COND




c Hsiao-Fan testing 3, Ellipse
       x(j)=2.d0*cos(1.d0*p2*dble(j)*h)
       y(j)=sin(1.d0*p2*dble(j)*h)


c**************************************


        end do

        call fd1(m,y,sy,h,dcut)
        call fd1(m,x,sx,h,dcut)

        do j=0,m-1
         slo(j)=dsqrt(sx(j)**2+sy(j)**2)
        end do
        slo(m)=slo(0)


        call arc(m,h,slo,alpha,tol)


        do j=0,m

c***************************
c THIS IS 2nd PLACE TO CHANGE
c       


c Hsiao-Fan testing 3, Ellipse
        xn(j)=2.d0*cos(1.d0*p2*alpha(j))
        yn(j)=sin(1.d0*p2*alpha(j))


        end do

   
     
       call fd1(m,xn,sxn,h,dcut)
       call fd1(m,yn,syn,h,dcut)

       sl=0.d0
       do j=0,m-1
         sl=sl+dsqrt(sxn(j)**2+syn(j)**2)
       end do
       sl = sl/dble(m)

c  now find theta, recover it from curvature.



        acount=0

        if (dabs(sxn(0)) .ge. tol) then
          t0=datan(syn(0)/sxn(0))
          else  if (sxn(0) .ge. 0) then
             acount=acount+1;
          end if
          if  (sxn(m-1) .ge. 0) then
             acount=acount+1;
          end if
          if  (sxn(m-2) .ge. 0) then
             acount=acount+1;
          end if
          if  (acount .ge. 2) then
             t0=p2/4.d0          
          else
             t0=-p2/4.d0          
        end if

       call thetsolve(m,h,p2,tol,alpha,sxn,syn,theta,sl,t0,dcut)
        
        
        x0=xn(0)
        y0=yn(0)
     
        
       write(6,*)'sxn(0)',sxn(0),'syn(0)',syn(0)
       write(6,*)'syn(0)/sxn(0)=',syn(0)/sxn(0),'t0=',t0
       write(6,*)'theta: ',theta(0)
       write(6,*)'x0: ',x0,'y0 :',y0

                
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c    THIS ROUTINE PROJECTS SL SO THAT
c
c    X_alpha=SL*COS(THETA)=1 + PERIODIC
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine stfilter(m,h,p2,dcut,sl,theta)
       implicit double precision (a-h,o-z)
       dimension theta(0:4800)
       dimension tempx(0:4800),tempxi(0:4800)

       do j=0,m-1
        tempx(j)=sl*cos(theta(j))
       end do
       tempx(m)=tempx(0)

       call fin(m,h,p2,tempxi,tempx)

       avx=tempxi(m)

       temp=sl/avx
       sl=temp

       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c    THIS ROUTINE FILTERS JUST **THETA**  TO ENSURE THAT THE
c    ALPHA COEFF is EQUAL TO 1:
c
c    X_alpha=SL*COS(THETA)=1 + PERIODIC
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tfilter(m,h,p2,dcut,sl,theta)
       implicit double precision (a-h,o-z)
       dimension theta(0:4800)
       dimension tempx(0:4800),tempxi(0:4800)
       dimension tempy(0:4800),tempyi(0:4800)

       do j=0,m-1
        tempx(j)=cos(theta(j))
        tempy(j)=sin(theta(j))
       end do
       tempx(m)=tempx(0)
       tempy(m)=tempy(0)

       call fin(m,h,p2,tempxi,tempx)
       call fin(m,h,p2,tempyi,tempy)

       avx=tempxi(m)
       avy=tempyi(m)

       do j=0,m-1
         tempx(j)=tempx(j)+ 1.d0/sl -avx
         tempy(j)=tempy(j)-avy
       end do
       tempx(m)=tempx(0)
       tempy(m)=tempy(0)

       do j=0,m-1
         theta(j)=datan(tempy(j)/tempx(j))
       end do
       theta(m)=theta(0)

c       call kfilter(m,dcut,sl)
c       call kfilter(m,dcut,theta)

       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   find max value of a vector
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine vmax(m,v,vmx)
       implicit double precision (a-h,o-z)
       dimension v(0:4800)

       vmx=0.d0
       do j=0,m
        vmx=dmax1(v(j),vmx)
       end do

       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    solve for tangent angle from curvature.
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine thetsolve(m,h,p2,tol,alpha,sxn,syn,theta,sl,t0,dcut)
        implicit double precision (a-h,o-z)
        dimension sxn(0:4800),syn(0:4800),theta(0:4800)
        dimension alpha(0:4800),dkap(0:4800)
        dimension sdkap(0:4800),dalpha(0:4800)
        dimension ssxn(0:4800),ssyn(0:4800)
        dimension stheta(0:4800)

        dimension temp2(0:4800),temp4(0:4800),temp1(0:4800)

        dimension un(0:4800),utt(0:4800),sun(0:4800)
        dimension xwy(0:4800),sxwy(0:4800),error(0:4800)
        dimension aux(0:4800)
        dimension un1(0:4800),un2(0:4800),aux3(0:4800)
        dimension dkaps(0:4800),un3(0:4800)
        dimension fntheta(0:4800)
        common /init/ x0,y0
        common /centroid/ xc,yc      
        common  a0,sl0        
        double precision w


        open(41,file='k1p.m')
        open(42,file='k2p.m')
        open(43,file='n1p.m')
        open(44,file='n2p.m')
        open(45,file='ns1p.m')
        open(46,file='ns2p.m')

        open(8,file='nl1p.m')
        open(9,file='nl2p.m')
        
        open(20,file='lin1p.m')
        open(21,file='lin2p.m')

        open(19,file='theta0p.m')
        
        write(41,*)'k1=['
        write(42,*)'k2=['
        write(43,*)'n1=['
        write(44,*)'n2=['
        write(45,*)'ns1=['
        write(46,*)'ns2=['
        
        write(19,*)'theta0=['

        
        write(8,*)'nl1=['
        write(9,*)'nl2=['
        
        write(20,*)'lin1=['
        write(21,*)'lin2=['


c  compute curvature from x,y



        call fd1(m,sxn,ssxn,h,dcut)
        call fd1(m,syn,ssyn,h,dcut)


        do j=0,m-1
         dkap(j)=(sxn(j)*ssyn(j)-syn(j)*ssxn(j))/(sl**3.d0)
        end do
        dkap(m)=dkap(0)
        write(6,*)dkap(j)

        write(6,*)'Initial curvature from x,y = ',dkap(0)

c  theta_alpha = sl* dkap

        do j=0,m
          stheta(j)=sl*dkap(j)-p2
        end do

c  integrate up. 



        call fin(m,h,p2,theta,stheta)

        write(6,*)'ave=',theta(m)/p2
        write(6,*)'+++++++++++++++++++++++++++*****************&&!!!***'
        write(6,*)'+++++++++++++++++++++++++++*****************&&!!!***'

c  put in constant
   
        
        do j=1,m
           theta(j)=theta(j)+t0+p2*dble(j)*h
        end do
        theta(0)=t0

c testing***************************************************************************
           call uset(m,h,p2,dcut,sl,theta,dalpha,w,un,utt)

         
         
        write(41,*)'];'
        write(42,*)'];'
        write(43,*)'];'
        write(44,*)'];'
        write(45,*)'];'
        write(46,*)'];'

        write(19,*)'];'

        write(8,*)'];'
        write(9,*)'];'
        
        write(20,*)'];'
        write(21,*)'];'
        

        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


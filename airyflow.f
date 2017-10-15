cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c This Computes the evolution of a curve whose normal vel. follows 
c                                  AIRY FLOW
c
c THIS CODE USES EXPLICIT ADAMS-BASHFORTH FOR THE SL EQUATION AND
c HAS AN INTEGRATING FACTOR FOR THE THETA EQUATION. BOTH ARE IN THE
c EQUAL ARCLENGTH FRAME. IT CALLS:
c                                   nl.f
c                                   library.f
c
c  THIS CODE COMPUTES THE CLOSED INTERFACE PROBLEM WITH PERIODIC
c  DATA ON [0,1]
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
        implicit double precision(a-h,o-z),integer(i-n)

        dimension xn(0:4800),yn(0:4800)

        dimension theta(0:4800),thetan(0:4800),thetann(0:4800)
        dimension stheta(0:4800)
        dimension ut(0:4800),un(0:4800),utt(0:4800)        
        dimension fnthetan(0:4800),ftheta(0:4800),fntheta(0:4800)
        
        dimension alpha(0:4800),dalpha(0:4800),sun(0:4800)
        
        dimension temp11(0:4800),temp12(0:4800)
        dimension temp13(0:4800),temp14(0:4800)
        dimension temp15(0:4800),temp16(0:4800)
        
        dimension area2(0:4800),areai(0:4800)
        dimension sumxx(0:4800),sumyy(0:4800)
        dimension sumxxi(0:4800),sumyyi(0:4800)
        
        dimension H1n(0:4800),H2n(0:4800),H3n(0:4800)

c for testing
        dimension auxsin(0:4800),auxsind(0:4800)
        dimension auxsini(0:4800),aauxsini(0:4800)



c for evolving the centroid

        dimension auxx(0:4800),auxy(0:4800)
        dimension auxxi(0:4800),auxyi(0:4800)


c this is for the first Euler step
        dimension sx(0:4800),ssx(0:4800)
        dimension sy(0:4800),ssy(0:4800)
        dimension dkap(0:4800)
        
        parameter(lda=33000)
        double precision temp1(lda) ,temp2(lda) 
        double precision temp3(lda) ,temp4(lda) 
        double precision temp5(lda) ,temp6(lda) 
        
        double precision temp11,temp12
        double precision temp13,temp14
        double precision temp15,temp16
        
        
        double precision unorm,utan,constant
        double precision auxa,auxb,aux1,aux2
        double precision aux1i,aux1r,aux2i,aux2r

        double precision avex,avey

        
        integer N,ncut,nstep



        common /init/ x0,y0
        common /length/salphainv,salpha
        common /centroid/ xc,yc         
        common  a0,sl0
        common  /Conserv/H10,H20,H30

        open(34,file='Deltap.m')
        open(32,file='yy0p.m')
        open(31,file='xx0p.m')
        open(30,file='r0p.m')
        open(29,file='deloradp.m')
        open(28,file='areap.m')

        open(1,file='xp.m')
        open(2,file='yp.m')

        open(3,file='theta1p.m')
        open(4,file='thetaEp.m')

        open(10,file='thetap.m')
        open(11,file='tp.m')
        open(12,file='slp.m')
        open(13,file='kp.m')
        open(14,file='ksp.m')
        open(15,file='kssp.m')
        open(16,file='Nlp.m')



        open(18,file='thetaADBp.m')
        open(22,file='theta3p.m')

        open(23,file='H1p.m')
        open(24,file='H2p.m')
        open(25,file='H3p.m')
      
        open(26,file='Spectp.m')
     
        open(70,file='kEp.m')
        open(71,file='ksEp.m')
        open(72,file='kssEp.m')

        open(73,file='kADBp.m')
        open(74,file='ksADBp.m')
        open(75,file='kssADBp.m')

        open(76,file='NlEp.m')
        
        open(90,file='xEp.m')
        open(91,file='yEp.m')
        open(92,file='xADBp.m')
        open(93,file='yADBp.m')
      

        open(96,file='dftheta1p.m')
        open(97,file='dftheta2p.m')


        print*,'t,N,dt,ncut,tp'
        read(5,*)t,N,dt,ncut,tp
      
        tau=1
      
        
        write(70,*)'kE=['
        write(71,*)'ksE=['
        write(72,*)'kssE=['
     
        write(73,*)'kADB=['
        write(74,*)'ksADB=['
        write(75,*)'kssADB=['
        
        
        write(90,*)'xE=['
        write(91,*)'yE=['
        write(92,*)'xADB=['
        write(93,*)'yADB=['
        
        write(96,*)'dftheta1=['
        write(97,*)'dftheta2=[' 

        write(1,*)'%  t,N,dt'
        write(1,*)'% ', t,N,dt
        write(1,*)'% ncut,tp ='
        write(1,*)'% ', ncut,tp

   
       dcut=10.d0**(-ncut)
       write(6,*)'dcut',dcut
       pi=4.d0*datan(1.d0)
       p2=2.d0*pi
       h=1.d0/dble(N)
       nstep=int(t/dt+0.00001d0)

       tol=1.d-14
       N2=N/2
       m=N
       dto=dt
       ktp=int(tp/dt+0.00001d0)
      
        write(11,*)'t=['
        write(12,*)'sl=['
        write(18,*)'thetaADB=['
        write(22,*)'theta3=['
        write(23,*)'H1=['
        write(24,*)'H2=['
        write(25,*)'H3=['
        write(29,*)'delorad=['
        write(30,*)'r0=['
        write(31,*)'xx0=['
        write(32,*)'yy0=['
        write(3,*)'theta1=['
        write(4,*)'thetaE=['
        write(34,*)'Delta=['
	


c ********    Initialization of the curve in (x, y) coordinates with equal arc-length ******************************************        
        ktime=0
        sallphainv=dble(0) 
        sallpha=dble(0) 

        call initiall(N,h,p2,sl,theta,dcut,tol,alpha,xn,yn,dalpha)

        call recon(N,h,p2,sl,theta,xn,yn,avx,avy,dcut)
        xn(N)=xn(0)
        yn(N)=yn(0)


        sl0 = sl
        write(6,*)'al ',sl

        salphainv=(p2/sl)
        salpha=sl/p2
        
        
        write(6,*)'SIgma inverse 0: ',salphainv



c ********initialize area, radius and centroid**************

	a0=0.d0
	r0=dble(0)
         Delta=dble(0)
         do i=0,m-1   
          a0=a0+((sin(theta(i))*xn(i)-cos(theta(i))*yn(i))
     $ /dble(2.0))*(sl/dble(m))
         end do 

        avex=dble(0)
        avey=dble(0)


        do  i=0,N-1   
          avex=avex+((sin(theta(i))*xn(i)**2 )/dble(2.0))
          avey=avey-cos(theta(i))*(yn(i)**2)/dble(2.0)
        end do
         avex=avex*(sl/dble(N))
         avey=avey*(sl/dble(N))
         r0=dsqrt(a0/pi)

         xc=avex/a0
         yc=avey/a0

c perturbation ******************

         
         delmax1=dble(0)
         aux1=xn(0)-xc
         aux2=yn(0)-yc
        
         delmax1=dsqrt(aux1*aux1+aux2*aux2)-r0
         delmax2=dble(0)
 
 
          do  i=0,m-1
           dx=xn(i)-xc
           dy=yn(i)-yc
           delmax1=min((dsqrt(dx*dx+dy*dy)-r0),delmax1)
           delmax2=max((dsqrt(dx*dx+dy*dy)-r0),delmax2)	 
	 end do
         
          delmax=max(abs(delmax1),abs(delmax2))
          Delta=delmax   
     
        
c compute relevan quantities for H quantities *************************       
         call Dacurv(theta,dkap,N,sl,dcut,h)
c          call acurv(N,sl,theta,dkap,dcut,h)
        
        call uset(N,h,p2,dcut,sl,theta,dalpha,w,un,utt)
        call fd1(N,un,sun,h,dcut)
        do j=0,N-1
            sun(j)=sun(j)/sl
	end do
	sun(N)=sun(0)
c******* Initial conserved quantities for mKdV***********************


        H10=dble(0)
        H20=dble(0)
        H30=dble(0)
        
        do j=0,N-1
         H10=H10+dkap(j)
         H20=H20+(dkap(j))**2
         H30=H30+(.5*un(j)**2-(dkap(j)**4)/8)
        end do        
        H10=H10*sl/dble(N)
        H20=H20*sl/dble(N)
        H30=H30*sl/dble(N)
                
        

c        do j=0,N-1
c         temp11(j)=dkap(j)
c         temp12(j)=(dkap(j))**2
c         temp13(j)=(.5*un(j)**2-(dkap(j)**4)/8)
c        end do        
c         temp11(N)= temp11(0)
c         temp12(N)= temp12(0)
c         temp13(N)= temp13(0) 
c        call fin(N,h,p2,H1n,temp11)
c        call fin(N,h,p2,H2n,temp12)
c        call fin(N,h,p2,H3n,temp13)      
c        H10=H1n(N)
c        H20=H2n(N)
c        H30=H3n(N)
  
  
        write(6,*)'H10= ',H10,'H20=',H20,'H30=',H30
  

        

c ********    Print/Save initialization******************************************    
        	write(6,*)'INITIAL VALUES, STEP ZERO '

          tcomp=dt*dble(ktime)
          np=ktime/ktp+1
            
          call screen(N,h,p2,np,tcomp,dcut,theta,ktime,dt,sl)
          call writout(N,h,p2,np,tcomp,dcut,theta,ktime,dt,sl)

		
        	write(6,*)'TESTTTTTTT'


c ********************************************************************         
c Here we start the second step using Euler  + integrating factor.
         ktime=1 

        	write(6,*)'theta(0)  ',theta(0)
     
         call fthetaim(m,h,p2,dcut,sl,theta,dalpha,fntheta)                 
        do j=0,m-1
         temp11(j)=fntheta(j)
         temp12(j)=theta(j)-p2*dble(j)*h
        end do
        temp11(m)=temp11(0)
        temp12(m)=temp12(0)

        call fast(temp11,m)
        call fast(temp12,m)
         

         aux1=dt*temp11(0)+temp12(0)
         aux2=dt*temp11(1)+temp12(1)

         temp13(0)=aux1
         temp13(1)=aux2
         
         do  j=1,m/2
         k=2*j
         rk=p2*dble(j)

c Local model	 
         rsl=(rk/sl)**3
         
         angle1=1.d0*dt*rsl
         
         angle1r=dcos(angle1)  
         angle1i=-dsin(angle1)  

         aux1=dt*temp11(k)+temp12(k)
         aux2=dt*temp11(k+1)+temp12(k+1)

         temp13(k)=angle1r*aux1-angle1i*aux2
         temp13(k+1)=angle1i*aux1+angle1r*aux2
        end do
              
        call fsst(temp13,m)        

        do j=0,m-1
         thetan(j)=temp13(j)+p2*dble(j)*h
        end do
        thetan(m)=thetan(0)+p2*dble(m)*h

        	write(6,*)'thetan(0)  ',thetan(0)
        
        
c    update initial interface position (x0,y0) using Euler **************
         call uset(N,h,p2,dcut,sl,theta,dalpha,w,un,utt) 
        rncomp2=-dcos(theta(0))
        rncomp1=dsin(theta(0))

        unormo=0.d0
        x0=x0+.5d0*dt*(2.d0*unorm*rncomp1)
        y0=y0+.5d0*dt*(2.d0*unorm*rncomp2)
        unormo=unorm        

c ** save initial and Euler step values for theta***
             do j=0,m-1
              aux=p2*dble(j)*h+dsin(p2*dble(j)*h+0*ktime)
              write(3,*)theta(j)
              aux=p2*dble(j)*h+dsin(p2*dble(j)*h+1*ktime)
              write(4,*)thetan(j)
             end do   
             
c save the second non linearity        
         call fthetaim(m,h,p2,dcut,sl,thetan,dalpha,fnthetan) 
         do j=0,m-1
            write(76,*)fnthetan(j)
         end do


      
c  ***********************

        	write(6,*)'EULER: STEP ONE  '
        do j=0,m-1
          auxx(j)=dabs(thetan(j)-theta(j))
        end do
        auxx(m)=auxx(0)
        aux2 = MAXVAL (auxx)
        write(6,*)'Max variation of theta0-thetEuler= ',aux2

          tcomp=dt*dble(ktime)
          np=ktime/ktp+1
        if (mod(ktime,ktp) .eq. 0) then
          call writout(N,h,p2,np,tcomp,dcut,thetan,ktime,dt,sl)
        endif        
        call screen(N,h,p2,np,tcomp,dcut,thetan,ktime,dt,sl)


        



c  now begin time iterations**********************

       do 100 ktime=2,nstep

c update theta using an integrating factor
        call fthetaim(m,h,p2,dcut,sl,thetan,dalpha,fnthetan) 
        
        if (ktime .eq. 2) then
        do j=0,m-1
            write(96,*)fnthetan(j)
        end do
        endif
        
        if (ktime .eq. 2) then
        do j=0,m-1
            write(97,*)fnthetan(j)
        end do
        endif

        do j=0,m-1
         temp11(j)=fnthetan(j)
         temp12(j)=thetan(j)-p2*dble(j)*h
         temp13(j)=fntheta(j)
        end do
        temp11(m)=temp11(0)
        temp12(m)=temp12(0)
        temp13(m)=temp13(0)

        call fast(temp11,m)
        call fast(temp12,m)
        call fast(temp13,m)
        


         auxa=3.d0*temp11(0)-temp13(0)
         auxb=3.d0*temp11(1)-temp13(1)
         
         temp14(0)=temp12(0)+.5d0*dt*auxa
         temp14(1)=temp12(1)+.5d0*dt*auxb

        do j=1,m/2

         k=2*j
         rk=p2*dble(j)

c Local model	 
         rsl=(rk/sl)**3

         angle1=-1.d0*dt*rsl
         angle2=-2.d0*dt*rsl
         
         angle1r=dcos(angle1)  
         angle1i=dsin(angle1)  
         
         angle2r=dcos(angle2)
         angle2i=dsin(angle2) 
          
         aux1=angle1r*temp12(k)-angle1i*temp12(k+1)
         aux4=angle1i*temp12(k)+angle1r*temp12(k+1)
         
         aux2=angle1r*temp11(k)-angle1i*temp11(k+1)
         aux3=angle2r*temp13(k)-angle2i*temp13(k+1)

         aux5=angle1i*temp11(k)+angle1r*temp11(k+1)
         aux6=angle2i*temp13(k)+angle2r*temp13(k+1)
         
         auxa=3.d0*(aux2)-(aux3)
         auxb=3.d0*(aux5)-(aux6)

         temp14(k)=aux1+.5d0*dt*auxa
         temp14(k+1)=aux4+.5d0*dt*auxb

        end do

        call fsst(temp14,m)

        do j=0,m-1
         thetann(j)=temp14(j)+p2*dble(j)*h
        end do
c        thetann(m)=thetann(0)
        thetann(m)=thetann(0)+p2*dble(m)*h
        
      
      
c    update initial interface position (x0,y0) using AB2.
         call uset(N,h,p2,dcut,sl,thetann,dalpha,w,un,utt)

        rncomp2=-dcos(thetan(0))
        rncomp2o=-dcos(theta(0))

        rncomp1=dsin(thetan(0))
        rncomp1o=dsin(theta(0))

        x0n=x0+.5d0*dt*(3.d0*unorm*rncomp1-1.d0*unormo*rncomp1o)
        y0n=y0+.5d0*dt*(3.d0*unorm*rncomp2-1.d0*unormo*rncomp2o)

c  update values

        x0=x0n
        y0=y0n

        unormo=unorm

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
          write(6,*)'WRITEOUT :D'
          call writout(N,h,p2,np,tcomp,dcut,thetan,ktime,dt,sl)
          write(6,*)'Screen :D'
          call screen(N,h,p2,np,tcomp,dcut,thetan,ktime,dt,sl)
        endif


c we save first ADB step to check*********	
	 if (ktime .eq. 2) then
	   write(6,*)'WE SAVE AND PRINT ADB third Step.'
            tcomp=dt*dble(ktime)
            np=ktime/ktp+1	 
            call screen(N,h,p2,np,tcomp,dcut,thetan,ktime,dt,sl)
	 
            call recon(N,h,p2,sl,thetan,xn,yn,avx,avy,dcut)
            xn(N)=xn(0)
            yn(N)=yn(0)
            call Dacurv(theta,dkap,N,sl,dcut,h)
c            call acurv(N,sl,theta,dkap,dcut,h)
            
           call uset(N,h,p2,dcut,sl,thetan,dalpha,w,un,utt)
           call fd1(N,un,sun,h,dcut)
          do j=0,N-1
            sun(j)=sun(j)/sl
	end do
	sun(N)=sun(0)
	do j=0,m-1
	   write(18,*)thetan(j)
            write(73,*)dkap(j)
            write(74,*)un(j)
            write(75,*)sun(j)
        end do

        endif
 100  continue



       write(3,*)'];'
       write(4,*)'];'
       write(11,*)'];'
       write(12,*)'];'
       write(18,*)'];'
       write(22,*)'];'
       
       write(23,*)'];'
       write(24,*)'];'
       write(25,*)'];'
      
       write(29,*)'];'
       write(30,*)'];'
       write(31,*)'];'
       write(32,*)'];'


       write(33,*)'];'
       write(34,*)'];'


       
        write(70,*)'];'
        write(71,*)'];'
        write(72,*)'];'
        write(73,*)'];'
        write(74,*)'];'
        write(75,*)'];'

        write(76,*)'];'



        write(90,*)'];'
        write(91,*)'];'
        write(92,*)'];'
        write(93,*)'];'

        write(96,*)'];'
        write(97,*)'];' 



       stop
       end
       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c    THIS ROUTINE PRINTS THE OUTPUT IN SCREEN
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine screen(N,h,p2,ktp,tcomp,dcut,theta,ktime,dt,sl)
        implicit double precision (a-h,o-z),integer(i-n)
        
        dimension xn(0:4800),yn(0:4800)
        dimension theta(0:4800),dkap(0:4800)
        dimension sun(0:4800),un (0:4800),w(0:4800)
        dimension utt(0:4800),dalpha(0:4800)
        dimension area2(0:4800),areai(0:4800)
        dimension sumxx(0:4800),sumyy(0:4800)
        dimension sumxxi(0:4800),sumyyi(0:4800)
        dimension fntheta(0:4800)
        
        dimension temp1(0:4800),temp2(0:4800),temp3(0:4800)
        dimension h1n(0:4800),h2n(0:4800),h3n(0:4800)
        double precision r0,ave,avex,avey,Delta,xc,yc

        common /length/salphainv,salpha


         integer N,i
     

        common /init/ x0,y0
c        common  a0,sl0
        common  /Conserv/H10,H20,H30        
        
        pi=4.d0*datan(1.d0)


        call recon(N,h,p2,sl,theta,xn,yn,avx,avy,dcut)
        xn(N)=xn(0)
        yn(N)=yn(0)
        

       r0=dble(3)
       Delta=dble(0)        
       ave=dble(0)
       avex=dble(0)
       avey=dble(0)
       xc=dble(0)
       yc=dble(0)
       
        do  i=0,N-1   
         ave=ave+((sin(theta(i))*xn(i)-cos(theta(i))*yn(i))
     $ /dble(2.0))*(sl/dble(N))
         avex=avex+sin(theta(i))*xn(i)**2/dble(2.0)
         avey=avey-cos(theta(i))*yn(i)**2/dble(2.0) 
 	end do
          avex=avex*(sl/dble(N))
          avey=avey*(sl/dble(N))

         xc=avex/ave
         yc=avey/ave

c************* Figure out perturbation and radius ***********
        r0=dsqrt(ave/pi)
        
         Delta=dble(0)
         delmax=dble(0)
         delmax1=dble(0)
         aux1=xn(0)-xc
         aux2=yn(0)-yc
        
         delmax1=dsqrt(aux1*aux1+aux2*aux2)-r0
         delmax2=dble(0)
 
 
          do  i=0,N-1
           dx=xn(i)-xc
           dy=yn(i)-yc
           delmax1=min((dsqrt(dx*dx+dy*dy)-r0),delmax1)
           delmax2=max((dsqrt(dx*dx+dy*dy)-r0),delmax2)	 
	 end do
         
          delmax=max(abs(delmax1),abs(delmax2))
          Delta=delmax  


           write(6,*)'Perturbation  = ',Delta
        
c compute relevan quantities *************************       
         call fthetaim(N,h,p2,dcut,sl,theta,dalpha,fntheta)                 

        call Dacurv(theta,dkap,N,sl,dcut,h)
c        call acurv(N,sl,theta,dkap,dcut,h)
 
        call uset(N,h,p2,dcut,sl,theta,dalpha,w,un,utt)
        call fd1(N,un,sun,h,dcut)
        do j=0,N-1
            sun(j)=sun(j)/sl
	end do
	sun(N)=sun(0)
c******* check conserved quantities for mKdV***********************
 
        ah1=dble(0)
        ah2=dble(0)
        ah3=dble(0)
        do j=0,N-1
         ah1=ah1+dkap(j)
         ah2=ah2+(dkap(j))**2
         ah3=ah3+(.5*un(j)**2-(dkap(j)**4)/8)
        end do        
        ah1=ah1*sl/dble(N)
        ah2=ah2*sl/dble(N)
        ah3=ah3*sl/dble(N)
        
c        ah1=dble(3)
c        ah2=dble(3)
c        ah3=dble(3)
c        do j=0,N-1
c         temp1(j)=dkap(j)
c         temp2(j)=(dkap(j))**2
c         temp3(j)=(.5*un(j)**2-(dkap(j)**4)/8)
c        end do        
c         temp1(N)= temp1(0)
c         temp2(N)= temp2(0)
c         temp3(N)= temp3(0)        
c        call fin(N,h,p2,h1n,temp1)
c        call fin(N,h,p2,h2n,temp2)
c        call fin(N,h,p2,h3n,temp3)
c        ah1=h1n(N)
c        ah2=h2n(N)
c        ah3=h3n(N)
        
        aux=(ah1-H10)/H10    
        ah1=dabs(aux)
        aux=(ah2-H20)/H20
        ah2=dabs(aux)
        aux=(ah3-H30)/H30
        ah3=dabs(aux)


        
        write(6,*)'*******new step*****++++++++++++++++++++++++++'
        write(6,*)'ktp = ',ktp,'tcomp=', tcomp              
        write(6,*)'H10= ',ah1,'H20=',ah2,'H30=',ah3
        write(6,*)'Length = ',sl,'Area = ',ave
        write(6,*)'x0 = ',xn(0),'y0 = ',yn(0),'xc = ',xc,'yc = ',yc
        write(6,*)'Perturbation = ',Delta,'r0=',r0
        
c here we print other  relevant quantities        
        
       	aux2 = MAXVAL (fntheta)
	aux1=MINVAL (fntheta)
      	write(6,*)'Maximum NL= ',aux2,'Min NL= ',aux1

        write(6,*)'thetan0 = ',theta(0)
        write(6,*)'Curvature 0 = ', dkap(0)
        write(6,*)'Normal 0 = ', un(0)
        write(6,*)'sun 0 = ', sun(0)

       	aux2 = MAXVAL (dkap)
	aux1=MINVAL (dkap)
      	write(6,*)'Maximum curvature= ',aux2,'Min curvature= ',aux1

       	aux2 = MAXVAL (un)
	aux1=MINVAL (un)
      	write(6,*)'Maximum normal= ',aux2,'Min normal= ',aux1

		
       	aux2 = MAXVAL (sun)
	aux1=MINVAL (sun)
      	write(6,*)'Maximum ks= ',aux2,'Min ks= ',aux1


        
       end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c    THIS ROUTINE SAVES THE OUTPUT  
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine writout(N,h,p2,ktp,tcomp,dcut,theta,ktime,dt,sl)
        implicit double precision (a-h,o-z),integer(i-n)
       dimension xn(0:4800),yn(0:4800)
       dimension theta(0:4800),dkap(0:4800),x(0:4800)
       dimension w(0:4800),sun(0:4800),un (0:4800)
       dimension utt(0:4800),dalpha(0:4800),ftheta(0:4800)
       dimension area2(0:4800),areai(0:4800)
       dimension sumxx(0:4800),sumyy(0:4800)
       dimension sumxxi(0:4800),sumyyi(0:4800)
       dimension temp1(0:4800),temp2(0:4800),temp3(0:4800)
       dimension h1n(0:4800),h2n(0:4800),h3n(0:4800)

       dimension y(0:4800),ytemp(0:4800),fnthetan(0:4800)
        integer N,i

c for evolving the centroid

       dimension auxx(0:4800),auxy(0:4800)
       dimension auxxi(0:4800),auxyi(0:4800)

c for recomputing x,y

       dimension xNew(0:4800),yNew(0:4800)
       dimension xtemp0(0:4800),ytemp0(0:4800)
       dimension xtemp1(0:4800),ytemp1(0:4800)
       dimension xtemp2(0:4800),ytemp2(0:4800)
       dimension xtemp3(0:4800),ytemp3(0:4800)

        double precision r0,ave,avex,avey,Delta,xc,yc

        common /init/ x0,y0
        common /length/salphainv,salpha
        common  a0,sl0
        common  /Conserv/H10,H20,H30

       pi=4.d0*datan(1.d0)

       write(1,*)'% t=',tcomp
       write(1,*)'x(:,',ktp,')=['
       write(2,*)'% t=',tcomp
       write(2,*)'y(:,',ktp,')=['

       write(10,*)'% t=',tcomp
       write(10,*)'theta(:,',ktp,')=['
      
       write(13,*)'% t=',tcomp
       write(13,*)'k(:,',ktp,')=['

       write(14,*)'% t=',tcomp
       write(14,*)'ks(:,',ktp,')=['
      
       write(15,*)'% t=',tcomp
       write(15,*)'kss(:,',ktp,')=['
      
       write(16,*)'% t=',tcomp
       write(16,*)'Nl(:,',ktp,')=['
      
       write(26,*)'% t=',tcomp
       write(26,*)'Spect(:,',ktp,')=['

c  construct necessary relevant quantities


        call recon(N,h,p2,sl,theta,x,y,avx,avy,dcut)

        x(N)=x(0)
        y(N)=y(0)
     

c ********reconstruct area radius and centroid**************

       r0=dble(3)
       Delta=dble(0)        
       ave=dble(0)
       avex=dble(0)
       avey=dble(0)
       xc=dble(0)
       yc=dble(0)
       
       do  i=0,N-1   
         ave=ave+((sin(theta(i))*x(i)-cos(theta(i))*y(i))
     $ /dble(2.0))*(sl/dble(N))
         avex=avex+sin(theta(i))*x(i)**2/dble(2.0)
         avey=avey-cos(theta(i))*y(i)**2/dble(2.0)
 	end do
	
	avex=avex*(sl/dble(N))
	avey=avey*(sl/dble(N))

c          ave=ave/dble(N)
c          avex=avex/dble(N)
c	 avey=avey/dble(N)
c          area2(N)=area2(0)
c          sumxx(N)=sumxx(0)
c          sumyy(N)=sumyy(0)
c         call fin(N,h,p2,areai,area2)
c         ave=areai(N)
c         call fin(N,h,p2,sumxxi,sumxx)
c         avex=sumxxi(N)
c         call fin(N,h,p2,sumyyi,sumyy)
c         avey=sumyyi(N)
 
         xc=avex/ave
         yc=avey/ave

        do 14 i=0,N-1
          write(1,*)x(i)-xc
          write(2,*)y(i)-yc
14     continue
        write(1,*)x(N)-xc
        write(2,*)y(N)-yc

     

c************* Figure out perturbation and output (x,y) ***********
        r0=dsqrt(ave/pi)
        

         


         Delta=dble(0)
         delmax=dble(0)
         delmax1=dble(0)
         aux1=x(0)-xc
         aux2=y(0)-yc
        
         delmax1=dsqrt(aux1*aux1+aux2*aux2)-r0
         delmax2=dble(0)
 
 
          do  i=0,N-1
           dx=x(i)-xc
           dy=y(i)-yc
           delmax1=min((dsqrt(dx*dx+dy*dy)-r0),delmax1)
           delmax2=max((dsqrt(dx*dx+dy*dy)-r0),delmax2)	 
	 end do
         
          delmax=max(abs(delmax1),abs(delmax2))
          Delta=delmax  

          delorad=delmax/r0



       call Dacurv(theta,dkap,N,sl,dcut,h)
c         call acurv(N,sl,theta,dkap,dcut,h)
     
      
       call uset(N,h,p2,dcut,sl,theta,dalpha,w,un,utt)
       call fd1(N,un,sun,h,dcut)
       call forctheta(N,h,p2,dcut,sl,theta,dalpha,ftheta)             
      
       call fthetaim(N,h,p2,dcut,sl,theta,dalpha,fnthetan) 

      
        do j=0,N-1
            sun(j)=sun(j)/sl
            write(10,*)theta(j)
            write(13,*)dkap(j)
            write(14,*)un(j)
            write(15,*)sun(j)
             aux=ftheta(j)+sun(j)
            write(16,*)fnthetan(j)
        end do
        sun(N)=sun(0)
        
           
       write(11,*)tcomp
       write(12,*)sl
       write(31,*)x0
       write(32,*)y0

c******* check conserved quantities for mKdV***********************

        
 
        ah1=dble(0)
        ah2=dble(0)
        ah3=dble(0)
        do j=0,N-1
         ah1=ah1+dkap(j)
         ah2=ah2+(dkap(j))**2
         ah3=ah3+(.5*un(j)**2-(dkap(j)**4)/8)
        end do        
        ah1=ah1*sl/dble(N)
        ah2=ah2*sl/dble(N)
        ah3=ah3*sl/dble(N)
   
c        ah1=dble(3)
c        ah2=dble(3)
c        ah3=dble(3)
c        do j=0,N-1
c         temp1(j)=dkap(j)
c         temp2(j)=(dkap(j))**2
c         temp3(j)=(.5*un(j)**2-(dkap(j)**4)/8)
c        end do        
c         temp1(N)= temp1(0)
c         temp2(N)= temp2(0)
c         temp3(N)= temp3(0)        
c        call fin(N,h,p2,h1n,temp1)
c        call fin(N,h,p2,h2n,temp2)
c        call fin(N,h,p2,h3n,temp3)
c        ah1=h1n(N)
c        ah2=h2n(N)
c        ah3=h3n(N)   
   
     
        aux=(ah1-H10)/H10    
        ah1=dabs(aux)
        aux=(ah2-H20)/H20
        ah2=dabs(aux)
        aux=(ah3-H30)/H30
        ah3=dabs(aux)
        
        write(23,*)ah1
        write(24,*)ah2
        write(25,*)ah3
        
        
        
c*******  Here we save theta spectrum (norm of coefficients in Fourier space)*****

        do j=0,N-1
         xtemp2(j)=theta(j)-p2*dble(j)*h
        end do
        xtemp2(N)=xtemp2(0)
        call fast(xtemp2,N)
        do j=1,N/2+1
         k=2*(j-1)
         aux=dsqrt(xtemp2(k)**2+xtemp2(k+1)**2)
         write(26,*)aux
        end do

c printing*****************


        write(6,*)'-********_________________*************-'
        write(6,*)
        write(6,*)'ktp = ',ktp,'tcomp=', tcomp
        write(6,*)'Length = ',sl,'Area = ',ave
        write(6,*)'x0 = ',x0,'y0 = ',y0,'xc = ',xc,'yc = ',yc
        write(6,*)'Perturbation = ',Delta,'r0=',r0
        write(6,*)'H10= ',ah1,'H20=',ah2,'H30=',ah3
        
        write(6,*)'thetan0 = ',theta(0)
        write(6,*)'check curvature and normal'
        write(6,*)'-*********************-'          
        write(6,*)'Curvature 0 = ', dkap(0)
        write(6,*)'Normal 0 = ', un(0)
        write(6,*)'sun 0 = ', sun(0)
	aux2 = MAXVAL (dkap)
        	write(6,*)'Max curvature= ',aux2
        	aux2 = MINVAL (dkap)
        	write(6,*)'Min  curvature= ',aux2
	aux2 = MAXVAL (un)
        	write(6,*)'Max normal= ',aux2
        	aux2 = MINVAL (un)
        	write(6,*)'Min  normal= ',aux2
   	aux2 = MAXVAL (sun)
        	write(6,*)'Max sun= ',aux2
        	aux2 = MINVAL (sun)
        	write(6,*)'Min  sun= ',aux2



        write(29,*)delorad
        write(30,*)r0
         	
        write(34,*)Delta
        write(28,*)sum2

        write(1,*)'];'
        write(2,*)'];'

        write(10,*)'];'
       
        write(13,*)'];'
        write(14,*)'];'
        write(15,*)'];'
        write(16,*)'];'

        write(26,*)'];'
       
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


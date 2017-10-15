ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   THIS ROUTINE COMPUTES THE FORCING ON THE RHS OF THE EQN
c   FOR SL.
c                       periodic on [0,1]
c                       and for closed curves
c
c   10/15/17
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine forcsl(m,h,p2,dcut,sl,theta,dalpha,fsl)
        implicit double precision (a-h,o-z)
        dimension theta(0:4800),dalpha(0:4800)
        dimension stheta(0:4800), un(0:4800),utt(0:4800)
        dimension temp(0:4800),tempi(0:4800),w(0:4800)
        double precision unorm,utan

        common /init/ x0,y0
        common /Nvel/ unorm

        call fd1x(m,theta,stheta,h,dcut)
c       call kfilter(m,dcut,stheta)

        call uset(m,h,p2,dcut,sl,theta,dalpha,w,un,utt)

        do j=0,m-1
          temp(j)=stheta(j)*un(j)
        end do
        temp(m)=temp(0)

        call fin(m,h,p2,tempi,temp)

        fsl=tempi(m)

        return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  THIS ROUTINE COMPUTES THE FORCING ON THE RHS OF THE EQUATION
c  FOR THETA, for the explicit method
c
c  10/15/17   periodic on [0,1]
c             and for closed curves
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine forctheta(m,h,p2,dcut,sl,theta,dalpha,ftheta)       
       
        implicit double precision(a-h,o-z)
        dimension theta(0:4800),ftheta(0:4800),dalpha(0:4800)
        dimension stheta(0:4800), un(0:4800),utt(0:4800)
        dimension sun(0:4800),ut(0:4800),w(0:4800)
        dimension dkap(0:4800),dkaps(0:4800)

        dimension temp2(0:4800),temp4(0:4800)

        double precision unorm,utan

        common /Nvel/ unorm
        common /init/ x0,y0
 
        call fd1x(m,theta,stheta,h,dcut)


        call uset(m,h,p2,dcut,sl,theta,dalpha,w,un,utt)
        call fd1(m,un,sun,h,dcut)

        call rtvel(m,h,p2,theta,un,ut,dcut,sl)       
        call acurv(m,sl,theta,dkap,dcut,h)


 
 

        do j=0,m-1
          ftheta(j)=-sun(j)/sl+dkap(j)*ut(j)
        end do
        ftheta(m)=ftheta(0) 
 
 
 

c       call kfilter(m,dcut,ftheta)
 
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   THIS ROUTINE COMPUTES THE MODIFIED FORCING OF THE
c   RHS FOR THE EQUATION FOR THETA FOR USE WITH THE
c   IMPLICIT SCHEME.
c
c     periodic on [0,1] and for closed curves.
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine fthetaim(m,h,p2,dcut,sl,theta,dalpha,fntheta) 
        implicit double precision (a-h,o-z)
        dimension theta(0:4800),ftheta(0:4800),dalpha(0:4800)
        dimension fntheta(0:4800),temp(0:4800),tempt(0:4800)
        dimension stheta(0:4800),ut(0:4800),temp2(0:4800)
        dimension un(0:4800),utt(0:4800),w(0:4800),compare(0:4800)
        dimension sun(0:4800),dkap(0:4800)


        double precision unorm,utan
        common /Nvel/ unorm    
        common /init/ x0,y0

            
        call uset(m,h,p2,dcut,sl,theta,dalpha,w,un,utt)

        call rtvel(m,h,p2,theta,un,ut,dcut,sl)       

        call Dacurv(theta,dkap,m,sl,dcut,h)


         do j=0,m-1
           fntheta(j)=dkap(j)*ut(j)
         end do
         fntheta(m)=fntheta(0)


        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

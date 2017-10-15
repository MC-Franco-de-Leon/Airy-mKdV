ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  compute the LHS of equation1
c  
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine eq1(theta,sl,sigma,vel)
        parameter (maxn=4800)
        implicit double precision(a-h,o-z)

        dimension temp(0:maxn),temp1(0:maxn),temp2(0:maxn),xi(1:4800)
        dimension temp3(0:maxn),temp4(0:maxn),temp11(0:maxn),tt(0:400)
        dimension sigma(0:maxn),vel(0:maxn),dkap(0:4800),qq(0:4800)
        dimension theta(0:maxn),stheta(0:maxn),tempi(0:maxn),a(1:4800)

      common /par/h,dcut,p2,m
c      common /space/theta,sl
c      common /precon/mpre

	call acurv(m,sl,theta,dkap,dcut,h)   
c make sure velocity is zero initially

          do j=0,m
           vel(j)=0.d0
          end do

c compute 1st term

       do j=0,m-1
        temp(j)=sigma(j)
        stheta(j)=dkap(j)*sl/p2
       end do
       temp(m)=temp(0)
       stheta(m)=stheta(0)

       call fast(temp,m)

       do j=1,m/2+1
        k=2*(j-1)
        temp1(k)=-temp(k)*dfloat(j-1)*dfloat(j-1)*p2*p2
        temp1(k+1)=-temp(k+1)*dfloat(j-1)*dfloat(j-1)*p2*p2
       end do
       temp1(m)=temp1(0)

       call fsst(temp1,m)

c compute 2nd term and 3rd term

       do j=0,m-1
          temp2(j)=stheta(j)*stheta(j)*sigma(j)
          temp3(j)=sl**stheta(j)*sigma(m+1)
       end do
       temp2(m)=temp2(0)
       temp3(m)=temp3(0)

c compute 4th term in  m+1  equation

       do j=0,m-1
        temp4(j)=stheta(j)*sigma(j)
       end do
       temp4(m)=temp4(0)

       call fin(m,h,p2,tempi,temp4)

       temp44=tempi(m)

c compute 5th term in m+1 equation   

       temp5=sigma(m+1)*sl

c now put all terms together

          do j=0,m-1
            vel(j)=-temp1(j)+temp2(j)+temp3(j)
          end do

          vel(m)=vel(0)
          vel(m+1)=temp44+temp5

c       call kfilter(m,dcut,fntheta)

       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


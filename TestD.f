
      implicit double precision(a-h,o-z),integer(i-n)
      common /essen/ ag,tau
      common /init/ x0,y0
      common  a0,sl0

      common /initvel/unorm
      common /symm/aniso,flux,tcomp,wc,md
      common /blob/del

      open(34,file='DFixtp.m')

      write(34,*)'DeltaFix=['

      Inipert=0.d0
      h=.5/10

      do i=1,10
          Inipert=Inipert+i*h
          call principal(IniPert,deltafixt,t,N,dt,ncut,tp,eps,mode,
     1                  willmore) 
          write(34,*)deltafixt
      end do

      write(34,*)'];'
      stop
      end

      

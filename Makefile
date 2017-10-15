#
#
#FFLAGS = -fast
#FFLAGS = -C -g
$LFLAGS = -g
LINK = gfortran 

hsclosed : FAST1.o airyflow.o nl.o library.o uset_un.o lambda11.o dgmres.o mat_L.o eq1.o 
         
	gfortran -O2 -o hsclosed FAST1.o airyflow.o nl.o library.o uset_un.o lambda11.o dgmres.o mat_L.o eq1.o 



FAST1.o : FAST1.f
	gfortran -O2 -c  FAST1.f

airyflow.o : airyflow.f
	gfortran -O2 -c  airyflow.f

nl.o : nl.f
	gfortran -O2 -c  nl.f


lambda11.o : lambda11.f
	gfortran -O2 -c lambda11.f

dgmres.o : dgmres.f
	gfortran -O2 -c dgmres.f

mat_L.o : mat_L.f
	gfortran -O2 -c  mat_L.f

eq1.o : eq1.f
	gfortran -O2 -c eq1.f

library.o : library.f
	gfortran -O2 -c  library.f




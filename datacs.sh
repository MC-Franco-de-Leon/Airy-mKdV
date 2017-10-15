#!/bin/bash

vars[0]="1.024,1024,.00008,13, .1024,\n
 t,   N,   dt, ncut, tp"
vars[1]="1.024,1024,.00004,13, .1024,\n
 t,   N,   dt, ncut, tp"
vars[2]="1.024,1024,.00002,13, .1024,\n
 t,   N,   dt, ncut, tp"
vars[3]="1.024,1024,.00001,13, .1024,\n
 t,   N,   dt, ncut, tp"
vars[4]="1.024,1024,.0000025,13, .1024,\n
 t,   N,   dt, ncut, tp"
vars[5]="1.024,1024,.000005,13, .1024,\n
 t,   N,   dt, ncut, tp"

#vars=('t1' 't2' 't3' 't4' 't5' 't6')

echo 'Iniciando compilación... '
echo '(Préndete un churro ke esto se lleva un rato)'
echo

count=0
for i in $( ls -d dir*/ ); do
    count=$((count+1))
    echo "Entrando a $i"
    cd $i
    echo "Dando caña en " $(pwd)
    echo -e ${vars[$count]}>TUPPUT.txt
    echo "Saliendo de $i"
    cd ..
    echo "Estamos en " $(pwd)
done
echo "Hecho."

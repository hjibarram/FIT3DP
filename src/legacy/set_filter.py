import numpy as np
import matplotlib.pyplot as plt
from numpy import float_
name='GALEX_NUV.dat'
name='GALEX_FUV.dat'
name='JOHNSON_B.dat'
name='JOHNSON_U.dat'
name='COUSINS_R.dat'
name='COUSINS_I.dat'
name='2MASS_H.dat'
name='2MASS_J.dat'
name='2MASS_Ks.dat'
name2='NUV_GALEX.txt'
name2='FUV_GALEX.txt'
name2='B_Johnson.txt'
name2='U_Johnson.txt'
name2='R_Cousins.txt'
name2='I_Cousins.txt'
name2='H_2MASS.txt'
name2='J_2MASS.txt'
name2='K_2MASS.txt'
con=0
f=open(name,'r')
f1=open(name2,'w')
for line in f:
    line=line.replace('\n','')
    if '2MASS' in name:
        fac=10000.
    else:
        fac=1.0
    data=line.split(' ')
    data=filter(None,data)
    if float_(data[1]) < 0:
        data[1]='0.0'
    data[0]=str(float_(data[0])*fac)
    data[1]=str(float_(data[1])/1.0)
    if (con % 6) == 0:
        f1.write(str('1  ')+data[0]+' '+data[1]+' \n')
    con=con+1
f.close()
f1.close()
   

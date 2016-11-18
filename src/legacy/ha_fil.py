import numpy as np
import matplotlib.pyplot as plt

vo=900.0
lo=7014.0
vel_light=299792.458
z=0.07691907746666148 #lo/6562.8-1.0 
xo=6562.8*(1+z)
wo=6400*(1+z)
w1=6700*(1+z)
wave=np.arange(wo, w1, 5.0)
bw=40
f=0.6*np.exp(-((xo-wave)/bw)**4./2.)**2.0
plt.plot(wave,f)
plt.show()
f1=open('ha_filter_sh.txt','w')
for i in range(0, len(f)):
    f1.write('1 '+str(wave[i])+'  '+str(f[i])+' \n')
f1.close() 

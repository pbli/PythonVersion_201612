import numpy as np
import matplotlib.pyplot as plt
import math

objN=40
result=np.genfromtxt("Scenario_crossing_left_to_right_50mph.csv",delimiter=',',skip_header=1)
#sizeA=result.shape  #size of array
r=result[:,2] #range
rdot=result[:,3]/3.6 #range rate
theta=result[:,4] #angle
t=result[:,6]# time

error=np.zeros(40)
for i in range(0,40):
    error[i]=np.mean(np.absolute(np.sqrt(result[:,7+i]*result[:,7+i]+result[:,i+47]*result[:,i+47])-r))#1D
print(np.argmin(error)) #object number is 32, it is 33rd object

RadarPx=result[:,7+32]
RadarPy=result[:,7+32+40]
RadarVx=result[:,7+32+80]
RadarVy=result[:,7+32+120]
RefPx=r*np.cos(-theta*math.pi/180)
RefPy=r*np.sin(-theta*math.pi/180)

plt.figure()
plt.plot(RadarPx,RadarPy,color='Red',linewidth=2.0,label='Radar')
plt.plot(RefPx,RefPy,color='Blue', linewidth=2.0,label='Reference')
plt.xlabel('x(m)')
plt.ylabel('y(m)')
plt.ylim([-1,5.5])
plt.legend()
plt.show(block=False)
plt.savefig('Trajectory.pdf')

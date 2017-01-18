import numpy as np
import matplotlib.pyplot as plt
import csv
from kalman import kalman

result=np.genfromtxt("Scenario_crossing_left_to_right_50mph.csv",delimiter=',',skip_header=1)
#reference data, polar coordinates
r=result[:,2] #m
rDot=result[:,3]/3.6 # m/s
theta=-result[:,4]*np.pi/180 #radius
t=result[:,6]# time in s
for i in range(0,len(t)-1):
    if np.abs(t[i+1]-t[i])>1:
        t[i+1]=t[i]+(t[i+2]-t[i])/2
t=t-t[0]

plt.plot(theta,label='theta')
theta=kalman(theta,0.0007,0.01,10,t)
plt.plot(theta,color='Red',label='filted theta',linewidth=2)
plt.legend()
plt.show(block=False)
plt.savefig('filtedTheta.pdf')

thetaDot=np.diff(theta)/np.diff(t)
plt.figure()
plt.plot(thetaDot,label='thetaDot')
thetaDot=kalman(thetaDot,0.00002,0.005,2,t[0:len(t)-1])
plt.plot(thetaDot,color='Red',label='filted thetaDot',linewidth=2)
plt.legend()
plt.show(block=False)
plt.savefig('filted thetaDot.pdf')

plt.figure()
plt.plot(rDot,label='rDot')
rDot=kalman(rDot,0.00001,0.1,1,t)
plt.plot(rDot,color='Red',label='filted rDot',linewidth=2)
plt.legend()
plt.show(block=False)
plt.savefig('filted rDot.pdf')

#find objectNO
error=np.zeros(40)
for i in range(0,40):
    error[i]=np.mean(np.absolute(np.sqrt(result[:,7+i]*result[:,7+i]+result[:,i+47]*result[:,i+47])-r))#1D
print(np.argmin(error)) #object number is 32, it is the 33rd object

#radar p v x y, ref p v x y
RadarPx=result[:,7+32]
RadarPy=result[:,7+32+40]
RadarVx=result[:,7+32+80]
RadarVy=result[:,7+32+120]
RefPx=r*np.cos(theta)
RefPy=r*np.sin(theta)
n=len(r)
RefVx=rDot[0:n-1]*np.cos(theta[0:n-1])-r[0:n-1]*np.sin(theta[0:n-1])*thetaDot
RefVy=rDot[0:n-1]*np.sin(theta[0:n-1])+r[0:n-1]*np.cos(theta[0:n-1])*thetaDot

plt.figure()
plt.plot(RefPx,RefPy,color='Blue', linewidth=2.0,label='Reference')
plt.plot(RadarPx,RadarPy,color='Red',linewidth=2.0,label='Radar')
plt.xlabel('x(m)')
plt.ylabel('y(m)')
plt.ylim([-1,5.5])
plt.legend()
plt.show(block=False)
plt.savefig('Trajectory.pdf')

f,axes=plt.subplots(2,2,sharex=True)
axes[0,0].plot(t,RefPx,color='Blue', linewidth=2.0,label='Ref')
axes[0,0].plot(t,RadarPx,color='Red',linewidth=2.0,label='Rad')
axes[0,0].legend(loc="upper left" )
axes[0,0].set_xlabel('t(s)')
axes[0,0].set_title('Px(m)')
axes[0,0].set_ylim([0,80])


axes[0,1].plot(t,RadarPx-RefPx, linewidth=2.0)
axes[0,1].set_title('ePx(m)')
axes[0,1].set_xlabel('t(s)')
axes[0,1].set_ylim([-0.3,0.4])


axes[1,0].plot(t,RefPy,color='Blue', linewidth=2.0,label='Ref')
axes[1,0].plot(t,RadarPy,color='Red',linewidth=2.0,label='Rad')
axes[1,0].set_xlabel('t(s)')
axes[1,0].set_title('Py(m)')
axes[1,0].legend(loc="lower left")
axes[1,0].set_ylim([-2,5])

axes[1,1].plot(t,RadarPy-RefPy, linewidth=2.0)
axes[1,1].set_xlabel('t(s)')
axes[1,1].set_title('ePY(m)')
plt.show(block=False)
plt.savefig('Position and error at X and Y.pdf')

f,axes=plt.subplots(2,2,sharex=True)
axes[0,0].plot(t[0:n-1],RefVx,color='Blue', linewidth=2.0,label='Ref')
axes[0,0].plot(t,RadarVx,color='Red',linewidth=2.0,label='Rad')
axes[0,0].set_xlabel('t(s)')
axes[0,0].set_title('Vx(m/s)')
axes[0,0].legend(loc="lower left")

axes[0,1].plot(t[0:n-1],RadarVx[0:n-1]-RefVx, linewidth=2.0)
axes[0,1].set_xlabel('t(s)')
axes[0,1].set_title('eVx(m)')

axes[1,0].plot(t[0:n-1],RefVy,color='Blue', linewidth=2.0,label='Ref')
axes[1,0].plot(t,RadarVy,color='Red',linewidth=2.0,label='Rad')
axes[1,0].set_xlabel('t(s)')
axes[1,0].set_title('Vy(m/s)')
axes[1,0].set_ylim([-2.5,2.5])
axes[1,0].legend(loc="upper left")

axes[1,1].plot(t[0:n-1],RadarVy[0:n-1]-RefVy, linewidth=2.0)
axes[1,1].set_xlabel('t(s)')
axes[1,1].set_title('eVy(m)')
plt.show(block=False)
plt.savefig('Velocity and error at X and Y.pdf')

meanPVXY=np.array([np.mean(RadarPx-RefPx),np.mean(RadarPy-RefPy),np.mean(RadarVx[0:n-1]-RefVx),np.mean(RadarVy[0:n-1]-RefVy)])
stdPVXY=np.array([np.std(RadarPx-RefPx),np.std(RadarPy-RefPy),np.std(RadarVx[0:n-1]-RefVx),np.std(RadarVy[0:n-1]-RefVy)])
itemName=['Px', 'Py', 'Vx', 'Vy']
output=zip(itemName,meanPVXY,stdPVXY)
with open('report.csv','w', newline='') as f:
    writer=csv.writer(f)
    writer.writerow(['Term','mean','std'])
    for meb in output:
        writer.writerows([meb])
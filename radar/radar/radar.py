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

plt.figure()
plt.plot(t,RefPx,color='Blue', linewidth=2.0,label='ReferencePositionX')
plt.plot(t,RadarPx,color='Red',linewidth=2.0,label='RadarPositionX')
plt.xlabel('t(s)')
plt.ylabel('Px(m)')
plt.ylim([0,80])
plt.legend()
plt.show(block=False)
plt.savefig('Position X.pdf')

plt.figure()
plt.plot(t,RadarPx-RefPx, linewidth=2.0,label='position error X')
plt.xlabel('t(s)')
plt.ylabel('ePx(m)')
plt.ylim([-0.3,0.4])
plt.legend()
plt.show(block=False)
plt.savefig('Position error X.pdf')

plt.figure()
plt.plot(t,RefPy,color='Blue', linewidth=2.0,label='ReferencePositionY')
plt.plot(t,RadarPy,color='Red',linewidth=2.0,label='RadarPositionY')
plt.xlabel('t(s)')
plt.ylabel('Py(m)')
plt.legend()
plt.ylim([-2,5])
plt.show(block=False)
plt.savefig('Position Y.pdf')

plt.figure()
plt.plot(t,RadarPy-RefPy, linewidth=2.0,label='position error Y')
plt.xlabel('t(s)')
plt.ylabel('ePY(m)')
plt.legend()
plt.show(block=False)
plt.savefig('Position error Y.pdf')

plt.figure()
plt.plot(t[0:n-1],RefVx,color='Blue', linewidth=2.0,label='ReferenceVelocityX')
plt.plot(t,RadarVx,color='Red',linewidth=2.0,label='RadarVelocityX')
plt.xlabel('t(s)')
plt.ylabel('Vx(m/s)')
plt.legend()
plt.show(block=False)
plt.savefig('Velocity X.pdf')

plt.figure()
plt.plot(t[0:n-1],RadarVx[0:n-1]-RefVx, linewidth=2.0,label='Velocity error X')
plt.xlabel('t(s)')
plt.ylabel('eVx(m)')
plt.legend()
plt.show(block=False)
plt.savefig('Velocity error X.pdf')

plt.figure()
plt.plot(t[0:n-1],RefVy,color='Blue', linewidth=2.0,label='ReferenceVelocityY')
plt.plot(t,RadarVy,color='Red',linewidth=2.0,label='RadarVelocityY')
plt.xlabel('t(s)')
plt.ylabel('Vy(m/s)')
plt.ylim([-2.5,2.5])
plt.legend()
plt.show(block=False)
plt.savefig('Velocity Y.pdf')
plt.figure()

plt.plot(t[0:n-1],RadarVy[0:n-1]-RefVy, linewidth=2.0,label='VelocityErrorY')
plt.xlabel('t(s)')
plt.ylabel('eVy(m)')
plt.legend()
plt.show(block=False)
plt.savefig('Velocity error at Y.pdf')

meanPVXY=np.array([np.mean(RadarPx-RefPx),np.mean(RadarPy-RefPy),np.mean(RadarVx[0:n-1]-RefVx),np.mean(RadarVy[0:n-1]-RefVy)])
stdPVXY=np.array([np.std(RadarPx-RefPx),np.std(RadarPy-RefPy),np.std(RadarVx[0:n-1]-RefVx),np.std(RadarVy[0:n-1]-RefVy)])
itemName=['Px', 'Py', 'Vx', 'Vy']
output=zip(itemName,meanPVXY,stdPVXY)
with open('report.csv','w', newline='') as f:
    writer=csv.writer(f)
    writer.writerow(['Term','mean','std'])
    for meb in output:
        writer.writerows([meb])
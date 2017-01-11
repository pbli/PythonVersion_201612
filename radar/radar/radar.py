import numpy as np
import matplotlib.pyplot as plt

objN=40
result=np.genfromtxt("Scenario_crossing_left_to_right_50mph.csv",delimiter=',',skip_header=1)
r=result[:,2] #range m
rDot=result[:,3]/3.6 #range rate m/s
theta=-result[:,4]*np.pi/180 #angle in radius
t=result[:,6]# time in s
for i in range(0,len(t)-1):
    if np.abs(t[i+1]-t[i])>1:
        t[i+1]=t[i]+(t[i+2]-t[i])/2

thetaDot=theta*0;
thetaDot[1:]=np.diff(theta)/np.diff(t)

# apply filters to smooth theta,thetaDot
'''
#median
for i in range(1,len(t)-1):
    thetaDot[i]=np.median(thetaDot[i-1:i+2])
    theta[i]=np.median(theta[i-1:i+2])
plt.plot(theta)
plt.show(block=False)
plt.figure()
plt.plot(thetaDot)
plt.show()
'''



'''
#kalman
T=np.mean(np.diff(t))
q=0.001 #control smoothness
X1=np.zeros([3,len(t)])
X1[1,1]=theta[1]
F=np.array([[1, T, 0.5*T*T], [0, 1, T], [0, 0, 1]])#system matrix
P=np.zeros([3,3])
H=np.array([[1,0,0]])
Q=np.array([[q,0,0],[0,q,0],[0, 0, q]])
R=0.0000001;#evaluate the measurement
for i in range(2,len(t)):
    X1[:,i]=F.dot(X1[:,i-1])
    P=F.dot(P.dot(F.T))+Q
    K=P.dot((H.T).dot(np.linalg.inv(H.dot(P.dot(H.T))+R)))
    X1[:,i]=X1[:,i]+K.dot(theta[i]-H.dot(X1[:,i]))
    P=(np.eye(3)-K.dot(H)).dot(P)
plt.plot(X1[1,:])
plt.plot(theta)
plt.show(block=False)
plt.figure()
plt.plot(thetaDot)
plt.plot(X1[2,:])
plt.show()
'''


#find objectNO
error=np.zeros(40)
for i in range(0,40):
    error[i]=np.mean(np.absolute(np.sqrt(result[:,7+i]*result[:,7+i]+result[:,i+47]*result[:,i+47])-r))#1D
print(np.argmin(error)) #object number is 32, it is the 33rd object

RadarPx=result[:,7+32]
RadarPy=result[:,7+32+40]
RadarVx=result[:,7+32+80]
RadarVy=result[:,7+32+120]
RefPx=r*np.cos(theta)
RefPy=r*np.sin(theta)
RefVx=rDot*np.cos(theta)-r*np.sin(theta)*thetaDot
RefVy=rDot*np.sin(theta)+r*np.cos(theta)*thetaDot


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
plt.plot(t,RefPx,color='Blue', linewidth=2.0,label='ReferenceVelocityX')
plt.plot(t,RadarPx,color='Red',linewidth=2.0,label='RadarVelocityX')
plt.xlabel('t(s)')
plt.ylabel('Px(m)')
plt.legend()
plt.show(block=False)
plt.savefig('Position X.pdf')

plt.figure()
plt.plot(t,RefPy,color='Blue', linewidth=2.0,label='ReferenceVelocityX')
plt.plot(t,RadarPy,color='Red',linewidth=2.0,label='RadarVelocityX')
plt.xlabel('t(s)')
plt.ylabel('Py(m)')
plt.legend()
plt.show(block=False)
plt.savefig('Position Y.pdf')

plt.figure()
plt.plot(t,RefVx,color='Blue', linewidth=2.0,label='ReferenceVelocityX')
plt.plot(t,RadarVx,color='Red',linewidth=2.0,label='RadarVelocityX')
plt.xlabel('t(s)')
plt.ylabel('Vx(m/s)')
plt.legend()
plt.show(block=False)
plt.savefig('Velocity X.pdf')

plt.figure()
plt.plot(t,RefVy,color='Blue', linewidth=2.0,label='ReferenceVelocityY')
plt.plot(t,RadarVy,color='Red',linewidth=2.0,label='RadarVelocityY')
plt.xlabel('t(s)')
plt.ylabel('Vy(m/s)')
plt.legend()
plt.show(block=False)
plt.savefig('Velocity Y.pdf')

meanPVXY=np.array([np.mean(RadarPx-RefPx),np.mean(RadarPy-RefPy),np.mean(RadarVx-RefVx),np.mean(RadarVy-RefVy)])
stdPVXY=np.array([np.std(RadarPx-RefPx),np.std(RadarPy-RefPy),np.std(RadarVx-RefVx),np.std(RadarVy-RefVy)])
OUTPUT=np.array([meanPVXY,stdPVXY])
np.savetxt('OUTPUT.csv',OUTPUT,delimiter=',',fmt='%1.5f')


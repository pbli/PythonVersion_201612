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
t=t-t[0]
#use filter to smooth theta
q=0.0007 #control smoothness
R=0.01#evaluate the measurement
X1=np.zeros([3,len(t)])
X1[0,0]=theta[0]
P=np.eye(3)*10
H=np.array([[1,0,0]])
Q=np.array([[q,0,0],[0,q,0],[0, 0, q]])
for i in range(1,len(theta)):
    T=t[i]-t[i-1]
    F=np.array([[1, T, 0.5*T*T], [0, 1, T], [0, 0, 1]])#system matrix
    X1[:,i]=F.dot(X1[:,i-1])
    P=F.dot(P.dot(F.T))+Q
    K=P.dot((H.T).dot(np.linalg.inv(H.dot(P.dot(H.T))+R)))
    X1[:,i]=X1[:,i]+K.dot(theta[i]-H.dot(X1[:,i]))
    P=(np.eye(3)-K.dot(H)).dot(P)
plt.plot(theta,label='theta')
plt.plot(X1[0,:],color='Red',label='filted')
plt.legend()
plt.show(block=False)
plt.savefig('theta.pdf')
theta=X1[0,:]

thetaDot=np.diff(theta)/np.diff(t)
#use filter to smooth theta
q=1/10000 #control smoothness
R=1/10#evaluate the measurement
X1=np.zeros([3,len(thetaDot)])
X1[0,0]=thetaDot[0]
P=np.eye(3)*1
H=np.array([[1,0,0]])
Q=np.array([[q,0,0],[0,q,0],[0, 0, q]])
for i in range(1,len(thetaDot)):
    T=t[i+1]-t[i]
    F=np.array([[1, T, 0.5*T*T], [0, 1, T], [0, 0, 1]])#system matrix
    X1[:,i]=F.dot(X1[:,i-1])
    P=F.dot(P.dot(F.T))+Q
    K=P.dot((H.T).dot(np.linalg.inv(H.dot(P.dot(H.T))+R)))
    X1[:,i]=X1[:,i]+K.dot(thetaDot[i]-H.dot(X1[:,i]))
    P=(np.eye(3)-K.dot(H)).dot(P)
plt.figure()
plt.plot(thetaDot,label='thetadot')
plt.plot(X1[0,:],color='Red',label='filted')
plt.legend()
plt.show(block=False)
plt.savefig('filted thetaDot.pdf')
thetaDot=X1[0,:]

#use filter to smooth rdot
q=1/100000 #control smoothness
R=1/10#evaluate the measurement
X1=np.zeros([3,len(rDot)])
X1[0,0]=thetaDot[0]
P=np.eye(3)*1
H=np.array([[1,0,0]])
Q=np.array([[q,0,0],[0,q,0],[0, 0, q]])
for i in range(1,len(rDot)):
    T=t[i]-t[i-1]
    F=np.array([[1, T, 0.5*T*T], [0, 1, T], [0, 0, 1]])#system matrix
    X1[:,i]=F.dot(X1[:,i-1])
    P=F.dot(P.dot(F.T))+Q
    K=P.dot((H.T).dot(np.linalg.inv(H.dot(P.dot(H.T))+R)))
    X1[:,i]=X1[:,i]+K.dot(rDot[i]-H.dot(X1[:,i]))
    P=(np.eye(3)-K.dot(H)).dot(P)
plt.figure()
plt.plot(rDot,label='rdot')
plt.plot(X1[0,:],color='Red',label='filted rdot')
plt.legend()
plt.show(block=False)
plt.savefig('filted rdot.pdf')
rDot=X1[0,:]

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
plt.legend()
plt.show(block=False)
plt.savefig('Position X.pdf')

plt.figure()
plt.plot(t,RefPy,color='Blue', linewidth=2.0,label='ReferencePositionY')
plt.plot(t,RadarPy,color='Red',linewidth=2.0,label='RadarPositionX')
plt.xlabel('t(s)')
plt.ylabel('Py(m)')
plt.legend()
plt.show(block=False)
plt.savefig('Position Y.pdf')

plt.figure()
plt.plot(t[0:n-1],RefVx,color='Blue', linewidth=2.0,label='ReferenceVelocityX')
plt.plot(t,RadarVx,color='Red',linewidth=2.0,label='RadarVelocityX')
plt.xlabel('t(s)')
plt.ylabel('Vx(m/s)')
plt.legend()
plt.show(block=False)
plt.savefig('Velocity X.pdf')

plt.figure()
plt.plot(t[0:n-1],RefVy,color='Blue', linewidth=2.0,label='ReferenceVelocityY')
plt.plot(t,RadarVy,color='Red',linewidth=2.0,label='RadarVelocityY')
plt.xlabel('t(s)')
plt.ylabel('Vy(m/s)')
plt.legend()
plt.show(block=False)
plt.savefig('Velocity Y.pdf')

meanPVXY=np.array([np.mean(RadarPx-RefPx),np.mean(RadarPy-RefPy),np.mean(RadarVx[0:n-1]-RefVx),np.mean(RadarVy[0:n-1]-RefVy)])
stdPVXY=np.array([np.std(RadarPx-RefPx),np.std(RadarPy-RefPy),np.std(RadarVx[0:n-1]-RefVx),np.std(RadarVy[0:n-1]-RefVy)])
OUTPUT=np.array([meanPVXY,stdPVXY])
np.savetxt('OUTPUT.csv',OUTPUT,delimiter=',',fmt='%1.5f')


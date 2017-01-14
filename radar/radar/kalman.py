#three order kalman filter, assume the 3rd derivative is constant
def kalman(x1,q,R,p,t):
    #x1: variables,
    #q: process noise
    #R: measurement noise
    #p: control P matrix
    #t: time series
    import numpy as np
    X1=np.zeros([3,len(t)])
    X1[0,0]=x1[0]
    P=np.eye(3)*p
    H=np.array([[1,0,0]])
    Q=np.array([[q,0,0],[0,q,0],[0, 0, q]])
    for i in range(1,len(t)):
        T=t[i]-t[i-1]
        F=np.array([[1, T, 0.5*T*T], [0, 1, T], [0, 0, 1]])#system matrix
        X1[:,i]=F@X1[:,i-1]
        P=F@P@(F.T)+Q
        K=P@(H.T)@np.linalg.inv(H@P@(H.T)+R)
        X1[:,i]=X1[:,i]+K@(x1[i]-H@X1[:,i])
        P=(np.eye(3)-K@H)@P
    return X1[0,:]
import csv
import numpy as np
import matplotlib.pyplot as plt
reader=csv.reader(open("Scenario_crossing_left_to_right_50mph.csv","r"),delimiter=',')
keys=next(reader)   #key of data
x=list(reader)      #data
objN=40
result=np.array(x).astype('float')#store data in an array as float type
#print(keys)
sizeA=result.shape  #size of array
for i in range(0,sizeA[1]):# plot all
    plt.figure(i)
    plt.plot(result[:,i])
    plt.show(block=False)
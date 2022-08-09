#!/usr/bin/python

# Acquire surface elevations based on python2 

import os 

pathname = os.path.abspath('.')
savePath = os.path.join(pathname,'gaugesVOF')
if not os.path.isdir(savePath):
    os.makedirs(savePath)

postPath = os.path.join(pathname,'sets')
if os.path.isdir(postPath):
    postPath = 'sets'
else:
    postPath = os.path.join(pathname,'postProcessing/surfaceElevation')
    if os.path.isdir(postPath):
        postPath = 'postProcessing/surfaceElevation'
    else:
        postPath = 'postProcessing/sampleDict'

# List of time dirs in order
a = os.listdir('./'+postPath)
a.sort(lambda a,b: cmp(float(a), float(b)))

# Get number of sensors
dir1 = os.path.join(pathname,postPath,a[int(len(a)/2.0)])
b = os.listdir(dir1)
nSens = 0
index = []
for i in range(len(b)):
    test1 = b[i].find('VOF') + 1
    test2 = b[i].find('alpha') + 1
    if test1 and test2:
        index.append(i)
        nSens += 1

first = True

for i in range(nSens):
    # Create files to write
    fileName = b[index[i]][0:b[index[i]].find('_')]
    fileW = open(os.path.join(savePath,fileName), 'w')
    print 'Sensor ' + '%i' % int(i+1) + ' of ' + '%i' % nSens + '.'

    # Read files time by time
    for j in range(len(a)):
        directory = os.path.join(pathname,postPath,a[j])
        try:
            fileR = open(os.path.join(directory,b[index[i]]), 'r')
        except:
            print 'WARNING - File not present: ' + os.path.join(directory,b[index[i]])
        else:
            data = fileR.read()
            fileR.close()
            data = data.split('\n')
                      
            if first: # First time step
                coord = j
                first = False
            
            x = []
            y = []
            z = []
            alpha = []
    
            # x y z alpha1 calculation
            for k in range(len(data)-1):
                line = data[k]
                line = line.split('\t') 
                
                z.append(float(line[0]))
                alpha.append(float(line[1]))
                
            # Integrate in Z
            wLevel = z[0]
            for k in range(len(z)-1):
                if (alpha[k]>0.5) and (alpha[k+1]<0.5):
                       wLevel = z[k] + (alpha[k]-0.5)*(z[k+1]-z[k])/(alpha[k]-alpha[k+1])
                if (alpha[k]==0.5):
                       wLevel = z[k]
            # Write to file
            time = a[j]
            fileW.write(time + ' ' + '%.6f' % wLevel + '\n')

    fileW.close()

print 'Done'

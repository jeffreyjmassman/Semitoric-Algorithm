#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Buggy, needs to be fixed
# UPDATE: Appears to be working!

from math import gcd
import numpy as np

p0 = np.array([0,1])
p1 = np.array([2,2])
p2 = np.array([3,3])
p3 = np.array([5,4])
p4 = np.array([6,4])
p5 = np.array([8,3])

h = np.array([0.5,0,4,0.3])

P = [p0,p1,p2,p3,p4,p5]

'''
Function that calculates the bound on packing a triangle on a vertex to the LEFT of an arbitrarily long contiguous sequence of
FAKE corners with UP cuts; however, this algorithm can most likely be generalized to all such scenarios

Parameters: P, a list of vertices of the polygon such that, if it is 1-indexed and of length m, then the 2nd element is
the corner to be packed at; the 1st element is the vertex 'before' this packing corner; elements 3 to m-1 are the marked 
vertices; and point m is the vertex after the last marked vertex in the sequence

            h, a vector of length m-3 corresponding to the 'heights' (rather, the positive distance below) of marked points
            for each marked point, in order
            
Possible errors with program: should not work when v1 - w1 = 0; this corresponds to a packed triangle with a vertical edge; 
as such , it will never 'cross over' the first cut since it is parallel to it
'''

def getConstraint(P,h):
    
    n = len(h)
    v = (P[0] - P[1]) / gcd(P[1][0] - P[0][0], P[1][1] - P[0][1]) 
    w = (P[2] - P[1]) / gcd(P[2][0] - P[1][0], P[2][1] - P[1][1]) # w and v are the primitive vectors
    v1 = v[0] 
    w1 = w[0]
    v2 = v[1]
    w2 = w[1]
    an = h[-1]
    beta = an * np.abs(v1 - w1)
    
    if beta > gcd(P[-1][0] - P[-2][0], P[-1][1] - P[-2][1]): # if we go beyond the last edge, scale it back
        an = gcd(P[-1][0] - P[-2][0], P[-1][1] - P[-2][1]) / np.abs(v1 - w1)
    alphas = np.array([an])
#     print(f'the first calculated alpha (an) is {an}')
    
    for k in range(n-2,-1,-1): 
        p = P[k+3]
        q = P[k+2]
        ak = alphas[-1] + (q[1] - p[1]) - ((w2 - v2)/(w1 - v1) - (k+1)) * (q[0] - p[0]) # k+1 due to off-by-one
        
#         print(f'previous alpha: {alphas[-1]}')
#         print(f'change in y: {q[1] - p[1]}')
#         print(f'slope: {((w2 - v2)/(w1 - v1) - (k+1))}')
#         print(f'change in x: {q[0] - p[0]}')
#         print(f' k = {k} and ak = {ak}')
        
        if ak > h[k]: # if ak is too big, we have to slide everything up
            d = ak - h[k]
#             print(f'd is {d}')
            
            if any(alphas - d < 0): # if we go outside the polygon, we have to reset...
#                 print(f'for ak = {ak}, we tried to slide up but failed; some alpha went out; have to start over')
                alphas = np.array([])
                beta = h[k] * np.abs(w1 - v1) # same formula Jason and I discovered
                
                if beta > gcd(p[0] - q[0], p[1] - q[1]):
                    ak = gcd(p[0] - q[0], p[1] - q[1]) / np.abs(w1-v1)
#                     print('after restarting, had to rescale ak')
                else:
                    ak = h[k]
                alphas = np.append(alphas, ak)
                
            else:
#                 print('Successfully slid up')
                alphas = np.append(alphas, ak) # ...otherwise, we slide everything up by d
                alphas = alphas - d
                
        elif ak < 0: # If alpha goes outside the polygon right away, we have to reset; nothing we can do
#             print(f'for ak = {ak}, we went outside the polygon; have to start over')
            alphas = np.array([])
            beta = h[k] * np.abs(w1 - v1) # same formula Jason and I discovered
                
            if beta > gcd(p[0] - q[0], p[1] - q[1]):
                    ak = gcd(p[0] - q[0], p[1] - q[1]) / np.abs(w1-v1)
                   
            else:
                ak = h[k]
            alphas = np.append(alphas, ak)
                
        else: # otherwise, we are good
#             print('so far so good...')
            alphas = np.append(alphas, ak)
#         print(alphas)
#     print(np.abs(w1 - v1))
    return gcd(P[2][0] - P[1][0], P[2][1] - P[1][1]) + alphas[-1] * np.abs(w1 - v1)


# In[3]:


np.array([1,2]) > 0


# In[6]:


# Test 1

p0 = np.array([12,0])
p1 = np.array([0,0])
p2 = np.array([8,4])
p3 = np.array([10,3])
p4 = np.array([12,0])

#p5 = np.array([8,3])

h = np.array([2,1])

P = [p0,p1,p2,p3,p4]
getConstraint(P,h)


# In[2]:


'''
Modified code found on https://stackoverflow.com/questions/9170838/surface-plots-in-matplotlib to help plot surface
for different marked point heights

DISCLAIMER: This is not a valid semitoric polygon; p0 = p4 does not satisfy the delzant condition. However, all the other
vertices are fine, and in fact that is enough; we assume p1 satisfies the Delzant condition
'''

from mpl_toolkits.mplot3d import Axes3D  
import matplotlib.pyplot as plt

# Polygon
p0 = np.array([12,0])
p1 = np.array([0,0])
p2 = np.array([8,4])
p3 = np.array([10,3])
p4 = np.array([12,0])

P = [p0,p1,p2,p3,p4]

H1 = np.arange(0, 4.0, 0.05)
H2 = np.arange(0, 3.0, 0.05)

lmda = np.zeros(shape = (len(H2),len(H1)))
                
col = 0                
for h1 in H1:
    row = 0
    for h2 in H2:
        h = np.array([h1,h2])
        lmda[row, col] = getConstraint(P,h)
        row += 1
    col += 1

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(H1, H2)

ax.plot_surface(X, Y, lmda)

ax.set_xlabel('h1')
ax.set_ylabel('h2')
ax.set_zlabel('constraint')

ax.view_init(30, 210) # change the view of the plot; first parameter rotates about h1, second about lambda

plt.show()


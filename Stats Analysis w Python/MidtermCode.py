import numpy as np
from scipy.linalg import null_space 
A = np.array([[2,0,0,0,1,0,-2,-3,0],
             [0,2,0,0,-1,0,0,0,-1],
             [0,0,2,0,-1,0,-2,0,0],
             [0,0,0,2,0,0,0,-1,-2],
             [0,0,0,0,-1,2,0,-2,0]]) #### Change this to correct A matrix #### 
ns = null_space(A) 
for i in ns:
    for j in i:
         print(round(j,0))
    

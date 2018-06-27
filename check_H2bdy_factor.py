# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 09:45:46 2018
@author: suksmono
script to check consistencies between the s-domain and q-domain
compensation factor in transforming k-body -> 2-body interactions
"""

from sympy import *
import numpy as np
# define s matrix
S3=[ [ 1, 1, 1], [ 1, 1,-1], [ 1,-1, 1], [ 1,-1,-1], \
     [-1, 1, 1], [-1, 1,-1], [-1,-1, 1], [-1,-1,-1] ]
Q3=[ [ 0, 0, 0], [ 0, 0, 1], [ 0, 1, 0], [ 0, 1, 1], \
     [ 1, 0, 0], [ 1, 0, 1], [ 1, 1, 0], [ 1, 1, 1] ]
# define function vq2s
def vq2s(x):
    return(1-2*x)
# --
delta=5;
# Compare energy levels in q-domain vs s-domain
print('\nLevels of energy in q-domain vs s-domain:\n')
for m in range(0, 7):
     v1=Q3[m];
     qn=v1[0]; qi=v1[1]; qj=v1[2];
     # sVec=[qn qi qj]; % h2*h1 -> h3
     #%% q-domain: eq.(8)
     H2q = 3*float(qn)+ float(qi and qj) - 2*float(qi and qn) \
          - 2*float(qj and qn);
     #%% s-domain: eq.(10)mod
     sn=vq2s(qn); si=vq2s(qi);sj=vq2s(qj);
     '''
     --------------------------------------------------------
       H-compensation on s-domain:  H^(sn, si, sj, delta);
       use delta = 5
     --------------------------------------------------------
     '''
     H2s3 =  (5*si)/4 + (5*sj)/4 - (5*sn)/2 + (5*si*sj)/4 \
            -(5*si*sn)/2 - (5*sj*sn)/2 + 15/4;
     print(m,': q =', v1,'->s =',S3[m],'Ground',(qn==(qi and qj)),'|', '|Eq=',H2q,'|Es=', H2s3 )

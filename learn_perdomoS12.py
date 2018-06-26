"""
Boolean reduction of k-body spin interaction to 2-body
for implementation of quantum computing on D-Wave
Created by: suksmono@{STEI-ITB,MDR Inc.}
Example: Case taken from Protein Foldings's of A. Perdomo-Ortiz paper, 
"Finding low-energy conformations of lattice protein models
by quantum annealing", Sci.Rep., 2012, 2:571
"""
from sympy import *
# define used symbols
q1, q2, q3, q4, q5, q6 = symbols('q1 q2 q3 q4 q5 q6')
s1, s2, s3, s4, s5, s6 = symbols('s1 s2 s3 s4 s5 s6')

'''
THE PROBLEM
'''
S12= 4 -3*q1 + 4*q2 -4*q1*q2 -q3 + q1*q3 -2*q2*q3 +4*q4 -2*q1*q4 - 8*q2*q4 + 5*q1*q2*q4 -2*q3*q4 + 5*q2*q3*q4 - q1*q2*q3*q4


''' 
k-body -> 2-body
function H2sub= H2body(x1,x2,y,d12)
% H_5(q1,q2,Q5;d5)=d5*(3*Q5+q1*q2-2*q1*Q5-2*q2*q5)
H2sub=d12*(3*y+x1*x2-2*x1*y-2*x2*y);
'''
def H2sub(x1,x2,y,d12):
    return(d12*(3*y+x1*x2-2*x1*y-2*x2*y))
    
d12=6; d34=4
S12_2bdy =simplify( S12.subs(
        {q1*q2:q5, q3*q4:q6})
    + H2sub(q1,q2,q5,d12)
    + H2sub(q3,q4,q6,d34)
    )

'''
RESULT
q-domain, 2-body interaction:
 6*q1*q2 + q1*q3 - 2*q1*q4 - 12*q1*q5 
 - 3*q1 - 2*q2*q3 - 8*q2*q4 - 12*q2*q5 
 + 5*q2*q6 + 4*q2 + 4*q3*q4 - 8*q3*q6 
 - q3 + 5*q4*q5 - 8*q4*q6 + 4*q4 
 - q5*q6 + 14*q5 + 10*q6 + 4
''' 
print('\n q-domain, 2-body interaction:\n',S12_2bdy)    

# define functions: q2s, s2q

def q2s(x):
#    return((1/2)*(1-x))
    return(1/2 -x/2)
    
S12s=simplify( S12_2bdy.subs(
        {q1:q2s(s1), q2:q2s(s2), q3:q2s(s3),
         q4:q2s(s4), q5:q2s(s5), q6:q2s(s6) }) )
print('\n s-domain, 2-body interaction:\n',S12s)


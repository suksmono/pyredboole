"""
------------------------------------------------------------
symbolic compution of factorization problem
Created by: Suksmono
problem adopted from:
https://qiita.com/YuichiroMinato/items/9c9fba2b5bc978ec9e73
------------------------------------------------------------
"""
from sympy import *
# define used symbols
q0, q1, q2, q3, d = symbols('q0 q1 q2 q3 d')

'''
we will calculate the Hamiltonian
H=(15 - p*q)^2, where p=(1+2q0+41q), q=(1+2q2)
'''
# it is known that the prime factors are odd numbers
p=1+2*q0+4*q1;
q=1+2*q2;

# Hamiltonian of the problem
H=(15-p*q)**2;
EQ0=expand(H);
# display result for checking
print('EQ0=',EQ0);

# substitute qi^2->qi
EQ1=simplify( \
        EQ0.subs({q0**2:q0, q1**2:q1, q2**2:q2}) \
    );

print('EQ1=\n',EQ1);


# k-body to 2-bsyms q3 d;ody transform
def H2sub(x1,x2,y,d12):
    return(d12*(3*y+x1*x2-2*x1*y-2*x2*y))
#
H2b=simplify( \
          EQ1.subs({q1*q2:q3}) \
        + H2sub(q1,q2,q3,d) \
    );
print('H2b= \n',H2b);

'''
COMPARE RESULTS
> Output of this program:
  H2b= d*(q1*q2 - 2*q1*q3 - 2*q2*q3 + 3*q3) + 16*q0*q1 - 56*q0*q2 
        + 128*q0*q3 - 52*q0 - 96*q1 - 52*q2 - 48*q3 + 196
> Hamiltonian on the web
  H = 128q_0q_3 + 16q_0q_1 - 56q_0q_2 - 52q_0 - 48q_1q_2 
     - 96q_1 - 52q_2 + 196 + \delta(3q_3+q_1q_2-2q_1q_3-2q_2q_3)
'''
"""
translate q-domain compensation factor H(qi, qj, qn;delta)
tp s-domain H(si, sj, qn; delta)
@author: suksmono@{STEI-ITB, MDR Inc.}
"""
from sympy import *
qi, qj, qn = symbols('qi qj qn')
si, sj, sn = symbols('si sj sn')

dij=symbols('dij')
# q-domain H2bdy
H2q = dij*(3*qn + qi*qj - 2*qi*qn - 2*qj*qn)

#define transform function
def q2s(x):
    return(1/2 -x/2)

# do transform
H2s =simplify( H2q.subs(                    \
        {qi:q2s(si), qj:q2s(sj), qn:q2s(sn),}) \
        )
#
print('\nq-domain H2body:\n',H2q)
print('\ns-domain H2body:\n',H2s)
# we can substitute dij with numeric value
# eg, for dij=5, we write H2s.subs({dij:5})

# CHECK consistencies, assume dij=5
# q-domain
Q3=[ [ 0, 0, 0], [ 0, 0, 1], [ 0, 1, 0], [ 0, 1, 1], \
     [ 1, 0, 0], [ 1, 0, 1], [ 1, 1, 0], [ 1, 1, 1] ]
# define function vq2s
def vq2s(x):
    return(1-2*x)
# --
#dij=5;
# Compare energy levels in q-domain vs s-domain
print('\nLevels of energy in q-domain vs s-domain:\n')
vq_val=[0,0,0]; # init vs array
vs_val=[0,0,0]; # init vs array
for m in range(0, 8):
     vq_val=Q3[m];
     qn_val=vq_val[0]; qi_val=vq_val[1]; qj_val=vq_val[2];

     H2q_val =  H2q.subs({qi:qi_val, qj:qj_val, qn:qn_val, dij:5})
     #H2q.subs({qi:qi_val, qj:qj_val, qn:qn_val, dij:5})

     # convert to s
     sn_val=vq2s(qn_val); si_val=vq2s(qi_val); sj_val=vq2s(qj_val)
     vs_val[0]=si_val; vs_val[1]=sj_val; vs_val[2]=sn_val
     '''
     --------------------------------------------------------
       H-compensation on s-domain:  H^(sn, si, sj, delta);
       use delta = 5
     --------------------------------------------------------
     '''
     H2s_val =  H2s.subs({si:si_val, sj:sj_val, sn:sn_val, dij:5})

     print(m,': q =', vq_val,'->s =',vs_val,'Ground',(qn_val==(qi_val and qj_val)), \
           '|', '|Eq=',H2q_val,'|Es=', H2s_val )

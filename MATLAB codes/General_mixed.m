clear all;
close all;
clc;
syms F a b c real
assume((0<F)&(F<1))

H = [1;0]; %|0>
V = [0;1]; %|1>
Id = eye(2,2);
Id2=kron(Id,Id);
sigmax = [0 1; 1 0];
rx2=kron(sigmax,sigmax);
xI=kron(sigmax,Id);
Ix=kron(Id,sigmax);

A=1/sqrt(2)*(kron(H, H) + kron(V, V));
B=1/sqrt(2)*(kron(H, V) - kron(V, H));
C=1/sqrt(2)*(kron(H, V) + kron(V, H));
D=1/sqrt(2)*(kron(H, H) - kron(V, V));

in=F*(A*A')+a*(D*D')+c*(B*B')+b*(C*C')
%%%%%%%%%%%%%%%%%%%%%%%
aqubit=[1,0;0,0];% Ancilla qubit
CNOT=[1,0,0,0;0,1,0,0;0,0,0,1;0,0,1,0]; %CNot gate between second qubit of the entangled pair and the ancilla qubit.
c2=kron(Id,CNOT);

CCNOT = [1, 0, 0, 0, 0, 0, 0, 0; %CNot gate between first qubit of the entangled pair and the ancilla qubit.
         0, 1, 0, 0, 0, 0, 0, 0;
         0, 0, 1, 0, 0, 0, 0, 0;
         0, 0, 0, 1, 0, 0, 0, 0;
         0, 0, 0, 0, 0, 1, 0, 0;
         0, 0, 0, 0, 1, 0, 0, 0;
         0, 0, 0, 0, 0, 0, 0, 1;
         0, 0, 0, 0, 0, 0, 1, 0];

m1=kron(Id,kron(Id,V*V')); % Measuring ancilla qubit in the state 1
m0=kron(Id,kron(Id,H*H')); % Measuring ancilla qubit in the state 0

had=[1/sqrt(2) 1/sqrt(2); 1/sqrt(2) -1/sqrt(2)]; % Hadamard gate
H2=kron(had,had);

%%%%%%%%%%%%%%%%%%%%%%% First PC(Z):
%rho_g1=m0*CCNOT*c2*kron(in,aqubit)*c2'*CCNOT'*m0'; 
rho_g1=m1*CCNOT*c2*kron(in,aqubit)*c2'*CCNOT'*m1';
rho_gf1=(PartialTrace(rho_g1,[3])); % the state of the system after measuring and tracing out the ancilla qubit
g1=simplify(trace(rho_gf1)) %probability of obtaining this result
rho_gfn1=simplify(rho_gf1/g1) %Normalized state

%%%%%%%%%%%%%%%%%%%%%% Second PC(Z):
inh2=simplify(H2*rho_gfn1*H2) % employing Hadamard gate on both qubits of the entangled pair

%rho_g2=m0*CCNOT*c2*kron(inh2,aqubit)*c2'*CCNOT'*m0'; 
rho_g2=m1*CCNOT*c2*kron(inh2,aqubit)*c2'*CCNOT'*m1';
rho_gf2=(PartialTrace(rho_g2,[3])); % the state of the system after measuring and tracing out the ancilla qubit
g2=simplify(trace(rho_gf2)) %probability of obtaining this result
rho_gfn1=simplify(rho_gf2/g2) %Normalized state
g_total=simplify(g1*g2)
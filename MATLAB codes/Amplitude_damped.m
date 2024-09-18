clear all;
close all;
clc;
syms r real
assume((0<r)&(r<1))

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

in=A*A';
%%% Amplitude damping noise:
e0 = [1 0; 0 sqrt(1-r)];
e1 = [0 sqrt(r); 0 0];

%%% One qubit of the entangled pair is affected by Pauli Z error:
%E0=kron(e0,Id);
%E1=kron(e1,Id);
%per=simplify(E0*in*E0'+E1*in*E1')

%%% Two qubits of the entangled pair is affected by Pauli Z error:
E0=kron(e0,e0);
E1=kron(e0,e1);
E2=kron(e1,e0);
E3=kron(e1,e1);
per2=simplify(E0*in*E0'+E1*in*E1'+E2*in*E2'+E3*in*E3')

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

inh=simplify(H2*per2*H2) % employing Hadamard gate on both qubits of the entangled pair

%rho_g1=m0*CCNOT*c2*kron(inh,aqubit)*c2'*CCNOT'*m0'; 
rho_g1=m1*CCNOT*c2*kron(inh,aqubit)*c2'*CCNOT'*m1';
rho_gf1=(PartialTrace(rho_g1,[3])); % the state of the system after measuring and tracing out the ancilla qubit
g1=simplify(trace(rho_gf1)) %probability of obtaining this result
rho_gfn1=simplify(rho_gf1/g1) %Normalized state

%%%%%%%%%%%%%%%%%%%%%% Second PC(Z):
inh2=simplify(H2*rho_gfn1*H2) % employing Hadamard gate on both qubits of the entangled pair

rho_g2=m0*CCNOT*c2*kron(inh2,aqubit)*c2'*CCNOT'*m0'; 
%rho_g2=m1*CCNOT*c2*kron(inh2,aqubit)*c2'*CCNOT'*m1';
rho_gf2=(PartialTrace(rho_g2,[3])); % the state of the system after measuring and tracing out the ancilla qubit
g2=simplify(trace(rho_gf2)) %probability of obtaining this result
rho_gfn1=simplify(rho_gf2/g2) %Normalized state
g_total=simplify(g1*g2)
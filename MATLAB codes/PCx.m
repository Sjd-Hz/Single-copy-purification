clear all;
close all;
clc;
syms F real
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

in=A*A';
%%% Bit flip or X error:
e0 = [sqrt(F) 0; 0 sqrt(F)];
e1 = [0 sqrt(1-F); sqrt(1-F) 0];

%%% One qubit of the entangled pair is affected by Pauli X error:
E0=kron(e0,Id);
E1=kron(e1,Id);
per=simplify(E0*in*E0'+E1*in*E1')

%%% Two qubits of the entangled pair is affected by Pauli X error:
%E0=kron(e0,e0);
%E1=kron(e0,e1);
%E2=kron(e1,e0);
%E3=kron(e1,e1);
%per2=simplify(E0*in*E0'+E1*in*E1'+E2*in*E2'+E3*in*E3')

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

%rho_g1=m0*CCNOT*c2*kron(per,aqubit)*c2'*CCNOT'*m0'; 
rho_g1=m1*CCNOT*c2*kron(per,aqubit)*c2'*CCNOT'*m1';
rho_gf1=(PartialTrace(rho_g1,[3])); % the state of the system after measuring and tracing out the ancilla qubit
g1=simplify(trace(rho_gf1)) %probability of obtaining this result
rho_gfn1=simplify(rho_gf1/g1) %Normalized state
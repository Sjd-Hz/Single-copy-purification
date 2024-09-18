clear all;
close all;
clc;

H = [1;0]; %|0>
V = [0;1]; %|1>

phi_p=1/sqrt(2)*(kron(H, H) + kron(V, V));
si_m=1/sqrt(2)*(kron(H, V) - kron(V, H));
si_p=1/sqrt(2)*(kron(H, V) + kron(V, H));
phi_m=1/sqrt(2)*(kron(H, H) - kron(V, V));

in=phi_p*phi_p'

%%%%%%%%%%%%%%%%%
%%%Puali error
pi=sqrt(0.6);
px=sqrt(0.3);
pz=sqrt(0.1);
e{1}=pi*[1 0;0 1];
e{2}=px*[0 1;1 0];
e{3}=pz*[1 0;0 -1];

E=cell(1,9);
k=1;
for j=1:3
    ii=1;
E{k}=kron(e{j},e{ii});
E{k+1}=kron(e{j},e{ii+1});
E{k+2}=kron(e{j},e{ii+2});
k=k+3;
end

rer=0;
for kk=1:9
    rer=(rer+E{kk}*in*E{kk}');
end
disp(rer)
F=trace(rer*(phi_p*phi_p'))
a=trace(rer*(si_p*si_p'))
b=trace(rer*(phi_m*phi_m'))
c=trace(rer*(si_m*si_m'))

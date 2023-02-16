function [Phi_Phi,Phi_F,Phi_R,A_e, B_e,C_e]=mpcgain(Ap,Bp,Cp,Nc,Np);
%Nc control horizon
%Np prediction horizon
%program created for assignment number 2
[m1,n1]=size(Cp);
[n1,n_in]=size(Bp);
%nl-- dimension of the state variables;
%ml-- number of outputs;
%n_in-- number of inputs.
%%%%%%%%%%%%%%%%
%Augment state equations
%%%%%%%%%%%%%%%%
A_e=eye(n1+m1,n1+m1); 
A_e(1:n1,1:n1)=Ap;
A_e(n1+1:n1+m1,1:n1)=Cp*Ap;

B_e=zeros(n1+m1,n_in);
B_e(1:n1,:)=Bp;
B_e(n1+1:n1+m1,:)=Cp*Bp;
C_e=zeros(m1,n1+m1);
C_e(:,n1+1:n1+m1)=eye(m1,m1);
%dimension of the extended state space
n=n1+m1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h(1,:)=C_e;   %scalar for single-input and single-output system
F(1,:)=C_e*A_e;   % row vector (1xn) one row all column 
for kk=2:Np
    h(kk,:)=h(kk-1,:)*A_e;   
    F(kk,:)= F(kk-1,:)*A_e;
end
v=h*B_e;
Phi=zeros(Np,Nc);  %declare the dimension of Phi
Phi(:,1)=v;        % first column
for i=2:Nc
   Phi(:,i)=[zeros(i-1,1);v(1:Np-i+1,1)];  %Toplitz matrix
end
BarRs=ones(Np,1);
Phi_Phi= Phi'*Phi;
Phi_F= Phi'*F;
Phi_R=Phi'*BarRs;


  
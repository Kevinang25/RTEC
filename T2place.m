function K=T2place(A,B,P);
%A,B are system matrix
%P is the desired characteristic polynomial
%(P=s^n+a_{n-1}s^{n-1}+\ldots+a_0)
%dimension of P is n+1
n=length(B);
Lc=ctrb(A,B);
m=rank(Lc);
if (m<n), disp('Uncontrollable'); return; end
[V1,D1,U1]=svd(Lc);
InLc=U1*inv(D1)*V1';
gamma=InLc(n,:);
Tinv=ctrb(A',gamma')';
[V2,D2,U2]=svd(Tinv);
T=U2*inv(D2)*V2';
Astf=Tinv*A*T;
for kk=1:n;
    Pi(kk)=P(n+2-kk);
end
Khat=Pi+Astf(n,:);
K=Khat*Tinv;
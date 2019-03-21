%[ccdf,pm] = es4mgh(x,a0,a,A,C,mu,gam,lam,chi,psi)
%Survivor function P(L>x) = 1-F(x) and
%tail conditional mean, E(L|L>x), of
%
%L = a0+a'*X+X'*A*X, where:
%
%   X = mu+W*gam+sqrt(W)*C*Z,
%   Z~N(0,I), and W~GIG(lam,chi,psi), i.e.,
%   X is distributed as multivariate GHyp.
%
%The actual computation is done in a mex file (es4mghint.mexw32 on a
%Windows machine). The mex file depends on libpthread-2.dll, which must
%reside in the same directory or somewhere on the Windows search path.
%pthread is distributed under the terms of the GNU LGPL, which can be
%found at http://www.gnu.org/licenses/lgpl.html
%
% (c) 2010 S.A. Broda
function [ccdf,pm]=es4mgh(xvec,a0,a,A,C,mu,gam,lam,chi,psi)
CAC=C'*A*C;CAC=(CAC+CAC')/2;
[P,L]=eig(CAC);omega=diag(L)';
muA=mu'*A;CP=C*P;gA=gam'*A;
c=a'*gam+2*muA*gam;
d=a'*CP+2*muA*CP;
e=2*gA*CP;
k=gA*gam;kk=a0+a'*mu+muA*mu;
qvec=xvec-kk;
ccdf=es4mghint(qvec,omega,k,c,d,e,lam,chi,psi,1);%use mex file
if nargout==2,
    pm=es4mghint(qvec,omega,k,c,d,e,lam,chi,psi,0);%use mex file
    pm=pm./ccdf+kk;
end
end

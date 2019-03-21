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
% (c) 2010 S.A. Broda

function [ccdf,pm] = es4mgh_vectorized(x,a0,a,A,C,mu,gam,lam,chi,psi)
nq = numel(x);ccdf = zeros(1,nq);CAC = C'*A*C;CAC = (CAC+CAC')/2;
[P,L] = eig(CAC);omega = diag(L)';muA = mu'*A;CP = C*P;gA = gam'*A;
c = a'*gam+2*muA*gam;d = a'*CP+2*muA*CP;e = 2*gA*CP;
de = d.*e;
d2 = d.*d;e2 = e.*e;k = gA*gam;kk = a0+a'*mu+muA*mu;q = x-kk;
LK2 = lklam(lam,chi,psi);
rp = 0.318309886183791;   %1/pi;
for ii = 1:nq
    ccdf(ii) = quadgk(@(s)imag(Q0Mcgf(1i*s,-1i*s*q(ii),1))./s,0,inf);
end
ccdf = .5+rp*ccdf;
if nargout>1
    pm = zeros(1,nq);
    for ii = 1:nq
        pm(ii) = quadgk(@(s)imag(Q0Mcgf(1i*s,-1i*s*q(ii),0))./s,0,inf);
    end
    pm = (.5*Q0Mcgf(0,0,0)+rp*pm)./ccdf+kk;
end
    function [M] = Q0Mcgf(s,t,docdf)
        nu = 1./(1-2*omega'*s);
        se2n = e2*nu;sd2n = d2*nu;sden = de*nu;s2 = s.*s;
        lrho = s*c+s2.*sden+.5*sum(log(nu));
        alpha1 = k*s+.5*s2.*se2n;alpha2 = .5*s2.*sd2n;
        cat=chi-2*(alpha2+t);pa=psi-2*alpha1;LK2kr=-LK2+lrho;
        M = exp(lklam(lam,cat,pa)+LK2kr);
        if docdf == 0%get M0s
            nu2 = nu.*nu;
            lrhop = c+2*s.*sden+2*s2.*(de.*omega*nu2)+omega*nu;
            alpha1p = k+s.*se2n+s2.*(e2.*omega*nu2);
            alpha2p = s.*sd2n+s2.*(d2.*omega*nu2);
            M0 = exp(lklam(lam+1,cat,pa)+LK2kr);%A&S, p. 483, Eq. 11.3.4
            dM0da1 = exp(lklam(lam+2,cat,pa)+LK2kr);
            M = M.*alpha2p+dM0da1.*alpha1p+M0.*lrhop;
        end
    end
end
function k = lklam(lam,chi,psi)
l2=0.693147180559945;%log(2)
if all(chi == 0)
    k = -lam*log(.5*psi)+gammaln(lam);
elseif all(psi == 0)
    k = lam*log(.5*chi)+gammaln(-lam);
else
    scp=sqrt(chi.*psi);
    k = l2+.5*lam*log(chi./psi)+log(besselk(lam,scp,1))-scp;
end
end
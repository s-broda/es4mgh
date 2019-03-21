!must be linked against libslatec (depends on dsort, zbesk, dgamln and dqag)
!when compiling with gfortran -fopenmp, link also against libgomp
MODULE es4mgh_mod
CONTAINS
      subroutine es4mgh(i, qvec,omega,k,c,d,e,lambda,chi,psi,n,ne,docdf)
          USE mexinterface
          IMPLICIT NONE
          INTEGER(kind=mwSize), INTENT(IN) :: n, ne
          DOUBLE PRECISION, INTENT(IN) :: qvec(:), k, c, d(:), e(:), lambda, chi, psi,  docdf
          DOUBLE PRECISION, INTENT (OUT) :: i(1:n)
          DOUBLE PRECISION, TARGET, INTENT(IN) :: omega(:)
          DOUBLE PRECISION, TARGET :: e2(1:ne), de(1:ne), d2(1:ne)
          DOUBLE PRECISION, POINTER :: p_omega_(:), p_d2_(:), p_e2_(:), p_de_(:)
          INTEGER, PARAMETER :: gk=1, limit=650, lenw=4*650
          DOUBLE PRECISION, PARAMETER :: pi=3.141592653589793D0, epsabs=1.0D-10, epsrel=1.0D-6
          DOUBLE PRECISION :: oabs(1:ne), dummy(1:ne), lk2_, themean, ub, ubnew, q_,&
            k_, c_,  lambda_, chi_, psi_, abserr, last, work(1:lenw)
          INTEGER :: is, docdf_, trans_, m,  neval, ier, iwork(1:limit)
          INTEGER(kind=mwSize) :: ne_
          !define global vars to pass extra params to integrand
          COMMON /qcommon/ q_
          !$OMP THREADPRIVATE(/qcommon/)
          COMMON /rest/ p_omega_,k_, c_, p_d2_,p_e2_,p_de_,lambda_,chi_,psi_,ne_, lk2_, trans_, docdf_
          k_=k;c_=c;lambda_=lambda;chi_=chi;psi_=psi;ne_=ne;docdf_=docdf;
          p_omega_=>omega;
          d2=d**2;p_d2_=>d2
          e2=e**2;p_e2_=>e2
          de=d*e;p_de_=>de
          lk2_=lklam(lambda,chi+(0.0D0,0.0D0),psi+(0.0D0,0.0D0))
          oabs=ABS(omega)
          CALL dsort(oabs,dummy,ne,-1)
          ub=1
          IF (docdf==1) THEN!bound integrand
              DO m=1,ne
                  IF (oabs(m)==0) EXIT
                  ubnew=(-1/DBLE(m))*(2.0D0*LOG(pi)+2*LOG(DBLE(m))+2.0D0*LOG(epsabs)&
                    -LOG(4.0D0)+SUM(LOG(2.0D0*oabs(1:m))))
                  ubnew=EXP(ubnew)
                  ubnew=SQRT(ubnew)/(1+SQRT(ubnew))!upper bound after transformation
                  IF (ubnew<ub) ub=ubnew
              END DO
          END IF
          trans_=1!transform integrand onto 0-1 to weaken singularity at left endpoint
          !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(is,abserr,neval,ier,last,iwork,work)
          DO is=1,n
              q_=qvec(is)
              CALL dqag(integrand, 0.0D0, ub, epsabs, epsrel,&
                gk, i(is), abserr, neval, ier, limit, lenw, last, iwork, work)
          END DO
          !$OMP END PARALLEL DO
          IF (docdf==1) THEN
              i=i/pi+.5D0
          ELSE
              trans_=0
              q_=0
              docdf_=2
              themean=integrand(0.0D0)
              i=i/pi+themean/2.0D0
          END IF
      END SUBROUTINE

      FUNCTION lklam(lambda, chi, psi)
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(IN) :: lambda
          COMPLEX*16, INTENT(IN) :: chi, psi
          COMPLEX*16 :: lklam, bk
          DOUBLE PRECISION, EXTERNAL :: dgamln
          INTEGER ierr
          IF (chi==0) THEN
               lklam=-lambda*LOG(psi/2.0D0)+dgamln(lambda,ierr)
          ELSEIF (psi==0) THEN
               lklam=lambda*LOG(chi/2.0D0)+dgamln(-lambda,ierr)
          ELSE
               bk = besselk(lambda,SQRT(chi*psi),1)
               lklam=LOG(2.0D0)+(lambda/2.0D0)*LOG(chi/psi)+LOG(bk)-SQRT(chi*psi)
          END IF
      END FUNCTION


      FUNCTION besselk(nu, z, scale)
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(IN) :: nu
          COMPLEX*16, INTENT(IN) :: z
          INTEGER, INTENT(IN) :: scale
          COMPLEX*16 :: besselk
          DOUBLE PRECISION :: zr, zi, cyr(1:1), cyi(1:1), fnu
          INTEGER :: nf, nz, ierr
          CALL zbesk(DBLE(z), AIMAG(z), ABS(nu), scale+1, 1, cyr, cyi, nz, ierr)
          besselk=cyr(1)+(0.0D0,1.0D0)*cyi(1)
      END FUNCTION


      FUNCTION integrand(svec) RESULT(i)
          USE mexinterface
          IMPLICIT NONE
          DOUBLE PRECISION, INTENT(IN) :: svec
          INTEGER(kind=mwSize) :: ne, ne_
          INTEGER ::  docdf_, trans_
          DOUBLE PRECISION :: q_, omega(1:ne_), k_, c_,  d2(1:ne_), e2(1:ne_), de(1:ne_),&
            lambda_, chi_, psi_, lk2_, i, ss, s2s
          COMPLEX*16 :: s, t, nu(1:ne_), se2n, sd2n, sden, sln, lrho,&
            a1, a2, lrhop, a1p, a2p, mgf, mgf0, dm0da1
          DOUBLE PRECISION, POINTER :: p_omega_(:), p_d2_(:), p_e2_(:), p_de_(:)
          COMMON /rest/ p_omega_,k_,c_,p_d2_,p_e2_,p_de_,lambda_,chi_,psi_,ne_, lk2_, trans_, docdf_
          COMMON /qcommon/ q_
          !$OMP THREADPRIVATE(/qcommon/)
          IF (trans_==1) THEN
            ss = svec / (1-svec)
            s2s =  ss**2
          ELSE
            s2s=svec
          END IF
          omega=p_omega_!dereference global pointers to local variables
          d2=p_d2_
          e2=p_e2_
          de=p_de_
          ne=ne_
          s=(0.0D0,1.0D0)*s2s
          t=-q_*s
          nu=1/(1-2.0D0*s*omega)
          se2n=SUM(e2*nu)
          sd2n=SUM(d2*nu)
          sden=SUM(de*nu)
          sln=SUM(LOG(nu))
          lrho=s*c_+s**2*sden+0.5D0*sln
          a1=k_*s+.5D0*s**2*se2n
          a2=.5D0*s**2*sd2n
          lrhop=c_+2.0D0*s*sden+2.0D0*s**2*SUM(de*omega*nu**2)+SUM(omega*nu)
          a1p=k_+s*se2n+s**2*SUM(e2*omega*nu**2)
          a2p=s*sd2n+s**2*SUM(d2*omega*nu**2)
          IF (docdf_==1) THEN !get M/s
                i=exp(lklam(lambda_,chi_-2.0D0*(a2+t),psi_-2.0D0*a1)-lk2_+lrho)/s !I is real, so here we implicitly take the real part
          ELSE
                mgf=exp(lklam(lambda_,chi_-2.0D0*(a2+t),psi_-2.0D0*a1)-lk2_+lrho)
                mgf0=exp(lklam(lambda_+1,chi_-2.0D0*(a2+t),psi_-2.0D0*a1)-lk2_+lrho) !A&S, p. 483, Eq. 11.3.4
                dm0da1=exp(lklam(lambda_+2,chi_-2.0D0*(a2+t),psi_-2.0D0*a1)-lk2_+lrho)
                IF (docdf_==0) THEN !get M0s/s
                    i=(mgf*a2p+dm0da1*a1p+mgf0*lrhop)/s !here too
                ELSE !get M0s (for the mean)
                    i=(mgf*a2p+dm0da1*a1p+mgf0*lrhop) !and here
                END IF
          END IF
          IF (trans_==1) THEN
            i =  2.0D0*ss * i / (1 - svec)**2
          END IF
      END FUNCTION
END MODULE

SUBROUTINE mexFunction(nlhs, plhs, nrhs, prhs)
    USE mexinterface
    IMPLICIT NONE
    integer(kind=mwSize) :: n, m, ne
    integer(kind=mwPointer) ::  plhs(nlhs), prhs(nrhs)
    INTEGER ::  nlhs, nrhs
    n = mxgetn(prhs(1))
    ne = mxgetn(prhs(2))
    plhs(1) = mxcreatedoublematrix(1_8,n,0)
    call resargs(%val(mxgetpr(plhs(1))),%val(mxgetpr(prhs(1))),%val(mxgetpr(prhs(2))), &
    %val(mxgetpr(prhs(3))),%val(mxgetpr(prhs(4)))  , %val(mxgetpr(prhs(5))), &
    %val(mxgetpr(prhs(6))),%val(mxgetpr(prhs(7))),%val(mxgetpr(prhs(8))),  &
    %val(mxgetpr(prhs(9))),n,ne,%val(mxgetpr(prhs(10))))
END SUBROUTINE
SUBROUTINE resargs(i, qvec,omega,k,c,d,e,lambda,chi,psi,n,ne,docdf)
    USE es4mgh_mod
    USE mexinterface
    IMPLICIT NONE
    DOUBLE PRECISION, intent(in) ::  qvec(n), omega(ne), k, c, d(ne), e(ne), lambda, chi, psi, docdf
    DOUBLE PRECISION, intent(out) :: i(n)
    INTEGER(kind=mwSize) :: n, ne
    call es4mgh(i, qvec,omega,k,c,d,e,lambda,chi,psi,n,ne,docdf)
END SUBROUTINE

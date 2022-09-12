
MODULE alib
	USE nrtype
	USE nrutil
	USE nr
	IMPLICIT NONE
	real(dp), PARAMETER :: SQRTPI=1.7724539_dp

	INTERFACE logit
	MODULE PROCEDURE logit_scalar,logit_vec,logit_scalar_dble,logit_vec_dble
	END INTERFACE
	
	
	INTERFACE min2pls
	MODULE PROCEDURE min2pls_scalar,min2pls_vec,min2pls_scalar_dble,min2pls_vec_dble
	END INTERFACE

	INTERFACE multinom
	MODULE PROCEDURE multinom_sngl, multinom_dble
	END INTERFACE


	INTERFACE condmom
	MODULE PROCEDURE condmom1_sngl, condmom2_sngl,condmom1_dble, condmom2_dble
	END INTERFACE

	INTERFACE one
	MODULE PROCEDURE one0, one1, one2, one3
	END INTERFACE

CONTAINS
	!Pdf of normal with mean mu and stdev sigma
	FUNCTION normpdf(x,mu,sigma)
	real(dp), DIMENSION(:), INTENT(IN) :: x
	real(dp), INTENT(IN) :: mu,sigma
	real(dp), DIMENSION(size(x)) :: normpdf
	normpdf=1.0_dp/ sigma*sqrt(2.0_dp)*SQRTPI * exp(-(x-mu)**2/2.0_dp*sigma**2 ) 
	END FUNCTION normpdf

	!Pdf of log-normal with mean mu and stdev sigma
	FUNCTION normpdfl(x,mu,sigma) 
	real(dp), DIMENSION(:), INTENT(IN) :: x
	real(dp), INTENT(IN) :: mu,sigma
	real(dp), DIMENSION(size(x)) :: normpdfl
	normpdfl=1.0_dp/ x*sigma*sqrt(2.0_dp)*SQRTPI * exp(-(log(x)-mu)**2/2.0_dp*sigma**2 ) 
	END FUNCTION normpdfl

	!For calculating conditional expectation of a normal random variable using nr's gauleg. Used by gauleg. 
	FUNCTION funcn(x,mu,sigma) 
	real(dp), DIMENSION(:), INTENT(IN) :: x
	real(dp), INTENT(IN) :: mu,sigma
	real(dp), DIMENSION(size(x)) :: funcn
	funcn=x*normpdf(x,mu,sigma)
	END FUNCTION funcn

	!Uses GAUSSBIN from SLgenlib to get the prob weights on the variables that are 
	!joint normal with mean 0 and sig(1:2) and ro and using a given grid
	SUBROUTINE anormjnt2(npoint,sig,ro,grid,wgtcond)
	INTEGER(I4B), INTENT(IN) :: NPOINT
	real(dp), INTENT(IN) :: sig(2),ro
	real(dp), DIMENSION(NPOINT), INTENT(IN) :: grid
	real(dp), DIMENSION(NPOINT,2,2) :: wgtcond
	INTEGER(I4B) :: i
	real(dp), PARAMETER :: mu(2)=0.0_dp
	real(dp) :: mq,sq,varcov(2,2),wgt(NPOINT,2)
	varcov(1,1)=sig(1)**2.0
	varcov(2,2)=sig(2)**2.0 
	varcov(1,2)=ro*sig(1)*sig(2)
	varcov(2,1)=varcov(1,2)
	do i=1,NPOINT
		mq=grid(i)*varcov(1,2)/varcov(1,1)
		sq=varcov(1,1)*(1.0_dp-ro**2)
		sq=sq**0.5_dp
		!call gaussbin(mq,sq,grid,wgtcond(:,i,1))

		mq=grid(i)*varcov(1,2)/varcov(2,2)
		sq=varcov(2,2)*(1.0_dp-ro**2)
		sq=sq**0.5_dp
		!call gaussbin(mq,sq,grid,wgtcond(:,i,2))
	end do 
	write(*,'(/1x,t3,a,t10,a,t22,a/)') '#','GRID(I)','WGT1COND2(:,I)','WGT2COND1(:,I)'
	do i=1,NPOINT
		write(*,'(1x,i2,3f12.6)') i,grid(i),wgtcond(:,i,1),wgtcond(:,i,2)
	end do
	!The below is to check the results of the above 
	! gaussbin- takes in paramters of continuous normal distribution (mu, sigma) and vector
	! of allowed values (values) and calculates fraction of outcomes closest to
	! each element in vector, i.e. a discrete approximation to that normal distribution
	!call gaussbin(mu(1),sig(1),grid,wgt(:,1))
	!call gaussbin(mu(2),sig(2),grid,wgt(:,2))		
	write(*,'(/1x,t3,a,t10,a,t22,a/)') '#','WGT(i,1)','SUM'
		do i=1,NPOINT
			write(*,'(I4,2F12.4)') i,wgt(i,1),sum(wgt(:,2)*wgtcond(:,i,1))
		end do 
	END SUBROUTINE anormjnt2

	! transforms number on real line into [0,1] using logit transformation
	FUNCTION logit_scalar(xx)
	real, INTENT(IN) :: xx
	real :: logit_scalar
	logit_scalar=EXP(xx)/(1.0_dp+EXP(xx))
	END FUNCTION logit_scalar

	FUNCTION logit_vec(xx)
	real, DIMENSION(:), INTENT(IN) :: xx
	real, DIMENSION(size(xx)) :: logit_vec
	logit_vec=EXP(xx)/(1.0_dp+EXP(xx))
	END FUNCTION logit_vec

	FUNCTION logit_scalar_dble(xx)
	real(dp), INTENT(IN) :: xx
	real(dp) :: logit_scalar_dble
	logit_scalar_dble=EXP(xx)/(1.0_dp+EXP(xx))
	END FUNCTION logit_scalar_dble

	FUNCTION logit_vec_dble(xx)
	real(dp), DIMENSION(:), INTENT(IN) :: xx
	real(dp), DIMENSION(size(xx)) :: logit_vec_dble
	logit_vec_dble=EXP(xx)/(1.0_dp+EXP(xx))
	END FUNCTION logit_vec_dble
    
	
	! inverse of logit transformation
	FUNCTION logitinv(xx)
	REAL(DP), INTENT(IN) :: xx  
	REAL(DP) :: logitinv
	REAL(DP) :: xmin,xmax
	xmin=0.000001_dp
	xmax=0.999999_dp     
	IF (xx<xmin) THEN
		logitinv=LOG(xmin/(1.0_dp-xmin))
	ELSEIF (xx>xmax) THEN
		logitinv=LOG(xmax/(1.0_dp-xmax))
	ELSE
		logitinv=LOG(xx/(1.0_dp-xx))
	ENDIF        
	END FUNCTION


	FUNCTION min2pls_scalar(xx)
	! transforms any real number into something at the [-1.0,1.0] interval
	real, INTENT(IN) :: xx
	real :: min2pls_scalar
	min2pls_scalar=2.0_dp*(1.0_dp/(1.0_dp+exp(-xx)))-1.0_dp 
	END FUNCTION
	FUNCTION min2pls_vec(xx)
	! transforms any real number into something at the [-1.0,1.0] interval
	real, DIMENSION(:), INTENT(IN) :: xx
	real, DIMENSION(size(xx)) :: min2pls_vec
	min2pls_vec=2.0_dp*(1.0_dp/(1.0_dp+exp(-xx)))-1.0_dp 
	END FUNCTION


	FUNCTION min2pls_scalar_dble(xx)
	! transforms any real number into something at the [-1.0,1.0] interval
	real(dp), INTENT(IN) :: xx
	real(dp) :: min2pls_scalar_dble
	min2pls_scalar_dble=2.0_dp*(1.0_dp/(1.0_dp+exp(-xx)))-1.0_dp 
	END FUNCTION 
	FUNCTION min2pls_vec_dble(xx)
	! transforms any real number into something at the [-1.0,1.0] interval
	real(dp), DIMENSION(:), INTENT(IN) :: xx
	real(dp), DIMENSION(size(xx)) :: min2pls_vec_dble
	min2pls_vec_dble=2.0_dp*(1.0_dp/(1.0_dp+exp(-xx)))-1.0_dp 
	END FUNCTION

	FUNCTION min2plsinv(yy)	
	! transforms any real number into something at the [-1.0,1.0] interval
	real(dp), INTENT(IN) :: yy
	real(dp) :: min2plsinv
	min2plsinv=logitinv(  (yy+1)/2.0_dp ) 
	END FUNCTION
    
    
	SUBROUTINE cholesky (a,n,p)
	implicit none
	INTEGER(I4B),intent(in)::n
	real(dp), DIMENSION(n,n),intent(in)::a
	real(dp), DIMENSION(n,n),intent(out)::p
	real(dp), DIMENSION(n,n)::aa
	real(dp) :: sum
	INTEGER(I4B) i,j,k
	sum = 0.0_dp
	aa(1:n,1:n)=a(1:n,1:n)
	do 13 i = 1,n
		do 12 j = i,n
			sum=aa(i,j)
			do 11 k = i-1,1,-1
				sum = sum - aa(i,k)*aa(j,k)
11			continue 
			if (i.eq.j) then
				if (sum.le.0.0_dp) stop 'choldc failed'
				p(i,i) = dsqrt(sum)
			else
				aa(j,i) = sum/p(i,i)
				p(j,i) = aa(j,i)
			end if
12		continue
13	continue
	return
	end SUBROUTINE cholesky

	SUBROUTINE getcholesky(sig,ro,cd)	
	real(dp), INTENT(IN) :: sig(2),ro  !,rho(2,2)	!rho is the correlation matrix and cd is the cholesky decomposition which is declared in global
	real(dp), INTENT(OUT) :: cd(2,2)
	INTEGER(I4B) :: j,k 
	real(dp) :: rho(2,2) 
	!**********************************************************************
	!*Call Cholesky routine that does the following:                      *
	!*Construct the covariance matrix btw male and female wage draws      *
	!*Take Cholesky Decomposition of the covariance matrix                *
	!**********************************************************************
	cd=0.0_dp 
	!Construct the correlation matrix. And then using the corr matrix, construct the cov matrix
	!sigma(1)=sig(1)   
	!sigma(2)=sig(2)	
	RHO(1,1)=1.0_dp		
	RHO(2,2)=1.0_dp		
	RHO(2,1)=ro		
	RHO(1,2)=RHO(2,1)
	!if (myid.eq.0) then
	!	print*, "Here are the sigmas: ", sig_w1,sig_w2,ro
	!	print*, "Here is the correlation matrix: "
	!	print*, RHO(1,1),RHO(1,2)
	!	print*, RHO(2,1),RHO(2,2)
	!end if
	!Now turn correlation matrix into varcov matrix
	do j=1,2
		do k=1,2
			RHO(j,k)=RHO(j,k)*sig(j)*sig(k)
		end do 
	end do 
	!if (myid.eq.0) then
	!	print*, "Here is the covariance matrix: "
	!	print*, RHO(1,1),RHO(1,2)
	!	print*, RHO(2,1),RHO(2,2)
	!end if
	!now get the cholesky decomposition
	call cholesky (RHO,2,CD)
	!!if (myid.eq.0) then
	!!	print*, "Here is the cholesky decomposition: ", sig_w1,sig_w2,ro
	!!	print*, CDM(1,1),CDM(1,2)
	!!	print*, CDM(2,1),CDM(2,2)
	!!end if
	END SUBROUTINE getcholesky

	! convert from n-dimensional index to linear one
	FUNCTION ndim2lin(dims,nindex)
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dims !like outpit of size(A)- vector
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: nindex !n-dimensional index- row vector
	INTEGER(I4B) :: ndim2lin
	INTEGER(I4B), DIMENSION(SIZE(dims)-1) :: cumprod ! running cumulative product of dims
	INTEGER(I4B) :: ii,ndims 
	ndims=SIZE(dims)
	cumprod(1)=dims(1)
	DO ii=2,ndims-1
		cumprod(ii)=dims(ii)*cumprod(ii-1)
	ENDDO
	ndim2lin=nindex(1)+SUM(cumprod*(nindex(2:ndims)-1))
	END FUNCTION ndim2lin
    
    ! convert from linear index to n-dimensional one
	FUNCTION lin2ndim(dims,linindex)
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: dims ! like outpit of size(A)
	INTEGER(I4B), INTENT(IN) :: linindex ! linear index
	INTEGER(I4B) :: lin2ndim(SIZE(dims)) 
	INTEGER(I4B) :: nn, prod, linindexcopy,ndims        
	ndims=SIZE(dims)
	linindexcopy=linindex
	DO nn=NDims,2,-1
		!WRITE(*,*) linindexcopy
		prod=PRODUCT(dims(1:nn-1))
		lin2ndim(nn)=CEILING(REAL(linindexcopy)/prod)
		linindexcopy=MOD(linindexcopy,prod)
		IF (linindexcopy==0) THEN
			linindexcopy=prod
		ENDIF
	ENDDO
	lin2ndim(1)=linindexcopy        
	END FUNCTION lin2ndim

	RECURSIVE FUNCTION recsum(n)  result(rec)
	!-----Recsum------------------------------------------------------
	!  FUNCTION to calculate sums recursively
	!---------------------------------------------------------------------
	INTEGER(I4B) :: rec
	INTEGER(I4B), intent(in) :: n
	if (n == 0) then
		rec = 0
	else
		rec = n + recsum(n-1)
	END if 
	END FUNCTION recsum

	! draw from a multinomal distribution.
	! Inputs are the vector of probabilities for each outcome and a U[0,1] random draw
	! Output is the index of the outcome
	FUNCTION multinom_sngl(probvec,randdraw)
		real(sp), DIMENSION(:), INTENT(IN) :: probvec ! probabilities of each state
		real(sp), INTENT(IN) :: randdraw ! U[0,1]
		real(sp) :: multinom_sngl
		INTEGER(I4B) :: nn
		real(sp) :: cutoff
		real(sp), DIMENSION(SIZE(probvec)) :: probvecnorm
		probvecnorm=probvec/SUM(probvec) ! normalize if probabilities don't sum to one 
		! cutoff is upper bound of interval that if randdraw falls in that interval, get current outcome
		cutoff=probvecnorm(1)
		nn=1
		DO WHILE (randdraw>cutoff)
			cutoff=cutoff+probvecnorm(nn+1)
			nn=nn+1
		ENDDO
		multinom_sngl=nn 
	END FUNCTION
	FUNCTION multinom_dble(probvec,randdraw)
		real(dp), DIMENSION(:), INTENT(IN) :: probvec ! probabilities of each state
		real(dp), INTENT(IN) :: randdraw ! U[0,1]
		real(dp) :: multinom_dble
		INTEGER(I4B) :: nn
		real(dp) :: cutoff
		real(dp), DIMENSION(SIZE(probvec)) :: probvecnorm
		probvecnorm=probvec/SUM(probvec) ! normalize if probabilities don't sum to one 
		! cutoff is upper bound of interval that if randdraw falls in that interval, get current outcome
		cutoff=probvecnorm(1)
		nn=1
		DO WHILE (randdraw>cutoff)
			cutoff=cutoff+probvecnorm(nn+1)
			nn=nn+1
		ENDDO
		multinom_dble=nn 
	END FUNCTION


	! condmom are subroutines that calculate conditional moments from arrays of data
	! condmom1 is for moments defined over a 1-dim array of observations
	! condmom2 is for moments defined over a 2-dim array of observations
	SUBROUTINE condmom1_sngl(im,cond,xx,moment,countmom,varmom)
	INTEGER(I4B), INTENT(IN) :: im
	LOGICAL, DIMENSION(:), INTENT(IN) :: cond ! condition that must be true for observation to contribute to the moment
	real(sp), DIMENSION(:), INTENT(IN) :: xx ! moment to be calculated
	real(sp), DIMENSION(:), INTENT(OUT) :: moment ! value of conditional moment
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: countmom ! number of observations contributing to each moment
	real(sp), DIMENSION(:), INTENT(OUT) :: varmom ! now CONDITIONAL variance of each moemnt	
	real(sp),PARAMETER :: d1=1.0_sp	
	countmom(im)=SUM(one(cond))
	moment(im)=SUM(xx,cond)/MAX(d1*countmom(im),0.1_sp)
	varmom(im)=SUM(xx**2,cond)/MAX(d1*countmom(im),0.1_sp) - (SUM(xx,cond)/MAX(d1*countmom(im),0.1_sp))**2
	END SUBROUTINE
	SUBROUTINE condmom2_sngl(im,cond,xx,moment,countmom,varmom)
	INTEGER(I4B), INTENT(IN) :: im
	LOGICAL, DIMENSION(:,:), INTENT(IN) :: cond ! condition that must be true for observation to contribute to the moment
	real(sp), DIMENSION(:,:), INTENT(IN) :: xx ! moment to be calculated
	real(sp), DIMENSION(:), INTENT(OUT) :: moment ! value of conditional moment
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: countmom ! number of observations contributing to each moment
	real(sp), DIMENSION(:), INTENT(OUT) :: varmom ! now CONDITIONAL variance of each moemnt	
	real(sp), PARAMETER :: d1=1.0_sp
	countmom(im)=SUM(one(cond))
	moment(im)=SUM(xx,cond)/MAX(d1*countmom(im),0.1_sp)
	varmom(im)=SUM(xx**2,cond)/MAX(d1*countmom(im),0.1_sp) - (SUM(xx,cond)/MAX(d1*countmom(im),0.1_sp))**2
	!varmom(im)=SUM(xx**2,cond)/SIZE(xx) - (SUM(xx,cond)/SIZE(xx))**2
	!ahu 081212: adding that whenever there is noone in that cell, the moment should be -1, not 0 as it is now 
	!because otherwise it is never clear, whether it is really 0 or that is just because there is noone in that cell 
	!ahu 082112 taking this out again because it causes the objective function to be too non-smooth
	!ahu 082112 if (countmom(im)==0) then
	!ahu 082112 	moment(im)=-1
	!ahu 082112 end if 
	!ahu 082112 taking this out again because it causes the objective function to be too non-smooth
	!ahu 101012: putting it back in temporarily
	!ahu 111412 if (countmom(im)==0) then
	!ahu 111412  	moment(im)=-1
	!ahu 111412 end if 
	END SUBROUTINE

	SUBROUTINE condmom1_dble(im,cond,xx,moment,countmom,varmom)
	INTEGER(I4B), INTENT(IN) :: im
	LOGICAL, DIMENSION(:), INTENT(IN) :: cond ! condition that must be true for observation to contribute to the moment
	real(dp), DIMENSION(:), INTENT(IN) :: xx ! moment to be calculated
	real(dp), DIMENSION(:), INTENT(OUT) :: moment ! value of conditional moment
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: countmom ! number of observations contributing to each moment
	real(dp), DIMENSION(:), INTENT(OUT) :: varmom ! now CONDITIONAL variance of each moemnt	
	real(dp),PARAMETER :: d1=1.0_dp	
	countmom(im)=SUM(one(cond))
	moment(im)=SUM(xx,cond)/MAX(d1*countmom(im),0.1_dp)
	varmom(im)=SUM(xx**2,cond)/MAX(d1*countmom(im),0.1_dp) - (SUM(xx,cond)/MAX(d1*countmom(im),0.1_dp))**2
	END SUBROUTINE
	SUBROUTINE condmom2_dble(im,cond,xx,moment,countmom,varmom)
	INTEGER(I4B), INTENT(IN) :: im
	LOGICAL, DIMENSION(:,:), INTENT(IN) :: cond ! condition that must be true for observation to contribute to the moment
	real(dp), DIMENSION(:,:), INTENT(IN) :: xx ! moment to be calculated
	real(dp), DIMENSION(:), INTENT(OUT) :: moment ! value of conditional moment
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: countmom ! number of observations contributing to each moment
	real(dp), DIMENSION(:), INTENT(OUT) :: varmom ! now CONDITIONAL variance of each moemnt	
	real(dp), PARAMETER :: d1=1.0_dp	
	countmom(im)=SUM(one(cond))
	moment(im)=SUM(xx,cond)/MAX(d1*countmom(im),0.1_dp)
	varmom(im)=SUM(xx**2,cond)/MAX(d1*countmom(im),0.1_dp) - (SUM(xx,cond)/MAX(d1*countmom(im),0.1_dp))**2
	END SUBROUTINE



	pure function one0(x)
	logical, intent(in):: x
	integer(i4b) one0
	one0=0
	if (x) one0=1
	end function one0

	pure function one1(x)
	logical, intent(in):: x(:)
	integer(I4B) one1(size(x))
	one1=0
	where (x) one1=1
	end function one1

	pure function one2(x)
	logical, intent(in):: x(:,:)
	integer(I4B) one2(size(x,1),size(x,2))
	one2=0
	where (x) one2=1
	end function one2

	pure function one3(x)
	logical, intent(in):: x(:,:,:)
	integer(I4B) one3(size(x,1),size(x,2),size(x,3))
	one3=0
	where (x) one3=1
	end function one3

	SUBROUTINE precisiontests
	! look at p.457 in the fortran book
	REAL(SP), PARAMETER :: dum=-10000000_sp,dum2=-1000000000_sp,dum5=0.00000000005_sp
	REAL(DP), PARAMETER :: dum5dp=0.00000000005_dp
	REAL(SP), PARAMETER :: dum6sp=0.12345972985_sp
	REAL(DP), PARAMETER :: dum6dp=0.12345972985_dp
	REAL(SP), PARAMETER :: dum7sp=12345972985.0_sp
	REAL(DP), PARAMETER :: dum7dp=12345972985.0_dp
	real(sp) :: test8sp=0.000000000050000000_sp !6.666666666666666 !, test2=0_dp, test3=0
	real(dp) :: test8dp=0.000000000050000000_dp !6.666666666666666 !, test2=0_dp, test3=0
	real(sp) :: test1sp=0.005_sp  !5_sp !0.00000000005_sp !  0 000 000_sp !6.666666666666666 !, test2=0_dp, test3=0
	real(dp) :: test1dp=0.005_dp  !5_sp !0.00000000005_sp !  0 000 000_sp !6.666666666666666 !, test2=0_dp, test3=0
	real(dp) :: test2=6.66666666666666 !, test2=0_dp, test3=0
	real(dp) :: test3=6.666666666666666_dp !, test2=0_dp, test3=0
	real(sp) :: test4=6.666666_sp
	REAL(SP) :: test5=6.60_sp,test6=6.60 !6.6666666_sp
	REAL(DP) :: test7=6.60_dp,test8=6.60,tempd !6.6666666_sp
	REAL(SP) :: temp,val(2)

	write(*,'(F20.15)') 2.0_dp+1.0_sp
	write(*,*) "Here is dum "
	write(*,'(2F20.2)') dum,dum2
	write(*,'(2F20.15)') dum5,dum5dp
	write(*,'(2F20.15,2I4)') dum6sp,dum6dp,precision(dum6sp),precision(dum6dp)
	write(*,'(2F20.2)') dum7sp,dum7dp
	write(*,*) kind(0.0),kind(0.0_sp),kind(0._sp)
	write(*,*) 2_sp*test1sp-0.01_sp
	write(*,*) "Here is test1sp and test1dp (0.005) "
	write(*,'(F10.5)') test8sp
	write(*,'(2F10.5)') test1sp,test1dp  !,kind(test1),precision(test1),kind(0.0),kind(0.0d0)  !,test2,test3
	write(*,*) test2,test3,kind(test2),kind(test3),precision(test2),precision(test3)
	write(*,*) test4,kind(test4),precision(test4)
	write(*,*) "test 5 and test 6 : "
	write(*,'(3F20.15)') test5,test6,test5-test6
	write(*,'(3F12.6)') test5,test6,test5-test6
	write(*,'(3F12.7)') test5,test6,test5-test6
	write(*,*) "test 7 and test 8 : "
	write(*,'(2F12.6,F20.15)') test7,test8,test7-test8
	write(*,'(3F20.15)') test7,test8,test7-test8
	write(*,'(3F12.6)') test7,test8,test7-test8
	write(*,'(3F12.7)') test7,test8,test7-test8
	write(*,*) test7,test8,test7-test8
	!cdfcheck=2.0_sp
	!print*, "Here sp ", cdfcheck,	sqrt(cdfcheck)
	!print*, "Here sp ", 2.0_sp,	sqrt(2.0_sp)
	!cdfcheck=2d0
	!print*, "Here d0 ", cdfcheck,	sqrt(cdfcheck)
	!print*, "Here d0 ", 2d0,	sqrt(2d0)
	!print*, "Here sp*d0 ", 1.0_sp*2d0,	sqrt(1.0_sp*2d0)
	!print*, "Here epsilon ", epsilon(1d0),cdfcheck
	!call precisiontests
	!write(*,*) 0.	
	!write(*,*) 0.0
	!write(*,*) 0._sp
	!write(*,*) 0.0_sp
	!write(*,*) PI,PI_D

	VAL = 12345783123.8976
	!TG = 12345678.8976_DP
	!WRITE(*,*) val(1),precision(val(1)),precision(1.),precision(ttt),kind(val(1))
	!WRITE(*,'(es14.8)') val(1)
	!WRITE(*,*) TG
	val(2) = 123456789.		! since val is declared as SP it will lose the last digit and will round off the 8 before it to 9 so that you'll only get 12345679 for this number (since the mantissa is only 8 digits)
	WRITE(*,*) val(2)    
	WRITE(*,'(es20.4)') val(2)
	WRITE(*,'(es20.7)') val(2)	! this is the best one, because it is exactly using 1 main and 7 decimal points and this is good since 8 is exactly the size of the mantissa. 
					! so then it doesn't truncate (see the one above this) or fill in the rest with weird unnecessary numbers (see the e ones below this)
	WRITE(*,'(es20.8)') val(2)
	WRITE(*,'(es20.9)') val(2)
	WRITE(*,'(es20.10)') val(2)
	WRITE(*,'(es15.7)') val(2)	! this is the best one too because it saves space 
	WRITE(*,'(e15.7)')  val(2)	
	WRITE(*,'(e15.8)')  val(2)	! this is the best one also 
	!WRITE(*,'(e15.9)')  val(2)	! this doesn't work 
	WRITE(*,*) "STOP "
	WRITE(*,*) val(2)		! for sp, mantissa can only be 8 digits so that the computer reads val as 12345679 (where 9 is for roudded off 8). so this writes 1.2345679E+08
	WRITE(*,'(F11.2)') val(2)	! the width of the field is not large enough for 123456789 (which has 9 digits). so this writes **************
	WRITE(*,'(F12.2)') val(2)	! width of the field is fine, but now it will fill up that 9th digit with SOME number (in this case 2) so that it will write 123456792.00
	WRITE(*,'(F13.2)') val(2)	! width of the field is fine, but now it will fill up that 9th digit with SOME number (in this case 2) so that it will write 123456792.00 (but with more space in the beginning than the above one) 
	WRITE(*,'(F20.2)') val(2)	! width of the field is fine, but now it will fill up that 9th digit with SOME number (in this case 2) so that it will write 123456792.00 (but with more space in the beginning than the above one) 				
	WRITE(*,*) "STOP "
	val(2)=12345678912345.
	WRITE(*,'(F12.2)') val(2)	! same explanations above (i.e. the comp will just add a different number 2. only now, the number does the width of field. 
	WRITE(*,'(F13.2)') val(2)   
	WRITE(*,'(F20.2)') val(2)   

	!write(*,'(3E15.8)') val
	!write(*,'(3E16.9)') val
	!write(*,'(3E20.12)') val
	!write(*,'(3F10.2)') val
	!write(*,'(3F15.2)') val
	!write(*,'(4EN10.1)') val,val(3)+100000.
	!write(*,'(4EN12.3)') val,val(3)+100000.

	write(*,*) kind(1d0),precision(1d0)
	PRINT*, TINY(TEMP),TINY(TEMPD) 
	PRINT*, EPSILON(TEMP),EPSILON(TEMPD)
	PRINT*, (TEMP >= TINY(TEMP))
	PRINT*, (TEMP >= 0.)
	PRINT*, (TEMP >= 0D0)

	write(*,*) kind(pi),precision(pi)
	write(*,*) kind(1d0),precision(1d0)
	print*, pi 
	!print*, rp(6:7),temp,eps,temp-eps
	temp=7365.370_sp
	print*, temp
	write(*,'(G10.3)') temp
	write(*,'(G14.7)') temp
	temp=7.124370_sp
	print*, temp
	write(*,'(G10.3)') temp
	write(*,'(G14.7)') temp
	stop
	temp=7365.370_sp
	write(*,'(F10.3)') temp
	write(*,'(F13.5)') temp
	write(*,'(F18.9)') temp
	write(*,'(F20.10)') temp
	write(*,'(F30.20)') temp
	write(*,'(F30.20)') 7365.370_sp
	stop
	temp=7365.370000_sp
	print*, temp
	write(*,'(F13.5)') temp
	write(*,'(F18.9)') temp
	write(*,'(F20.10)') temp
	write(*,'(F30.20)') temp
	write(*,'(F30.20)') 7365.370_sp
	stop
	temp=7365.3701171875_sp
	write(*,'(F13.5)') temp
	write(*,'(F18.9)') temp
	write(*,'(F20.10)') temp
	write(*,'(F30.20)') temp
	write(*,'(F30.20)') 7365.3701171875_sp
	temp=7365.3701171875123456789_sp
	write(*,'(F13.5)') temp
	write(*,'(F18.9)') temp
	write(*,'(F20.10)') temp
	write(*,'(F30.20)') temp
	tempd=6.66666666666666
	write(*,*) tempd
	tempd=6.6666_dp
	write(*,*) tempd
	write(*,'(F30.20)') tempd
	write(*,'(F30.15)') tempd
	write(*,'(F20.14)') tempd
	temp=6.6666_sp
	write(*,*) temp
	write(*,'(F30.20)') temp
	write(*,'(F30.15)') temp
	write(*,'(F20.14)') temp
	write(*,'(F20.7)') temp
	write(*,'(F20.6)') temp
	write(*,'(F20.5)') temp
	temp=6.666666_sp
	write(*,*) temp
	write(*,'(F30.20)') temp
	write(*,'(F30.15)') temp
	write(*,'(F20.14)') temp
	write(*,'(F20.7)') temp
	write(*,'(F20.6)') temp
	write(*,'(F20.5)') temp
	stop
	tempd=6.6666666666666_dp
	write(*,*) tempd
	write(*,'(F30.20)') tempd
	tempd=6.66666666666666_dp
	write(*,*) tempd
	write(*,'(F30.20)') tempd
	write(*,'(F40.30)') tempd
	tempd=6.66660000000_dp
	write(*,*) tempd
	write(*,'(F20.14)') tempd
	write(*,'(F30.20)') tempd

	stop

	tempd=6.66666666666666_dp
	write(*,*) tempd
	tempd=6.666666666666666_dp
	write(*,*) tempd
	tempd=6.66666666666666666_dp
	write(*,*) tempd
	tempd=6.666666666666666_dp
	write(*,*) tempd
	tempd=6.666666666666666666_dp
	write(*,*) tempd
	write(*,*) " ********************** TEMPD   ******************************** "
	tempd=7365.370_dp
	print*, tempd
	write(*,'(F13.5)') tempd
	write(*,'(F18.9)') tempd
	write(*,'(F20.10)') tempd
	write(*,'(F30.20)') tempd


	!VAL = 12345783123.8976
	!print*, val(2) 

	!TG = 12345678.8976_DP
	!WRITE(*,*) val(1),precision(val(1)),precision(1.),precision(ttt),kind(val(1))
	!WRITE(*,'(es14.8)') val(1)
	!WRITE(*,*) TG
	val(2) = 123456789.		! since val is declared as SP it will lose the last digit and will round off the 8 before it to 9 so that you'll only get 12345679 for this number (since the mantissa is only 8 digits)
	WRITE(*,*) val(2)    
	WRITE(*,'("es20.4",es20.4)') val(2)
	WRITE(*,'("es20.7",es20.7)') val(2)	! this is the best one, because it is exactly using 1 main and 7 decimal points and this is good since 8 is exactly the size of the mantissa. 
					! so then it doesn't truncate (see the one above this) or fill in the rest with weird unnecessary numbers (see the e ones below this)
	WRITE(*,'("es20.8",es20.8)') val(2)
	WRITE(*,'("es20.9",es20.9)') val(2)
	WRITE(*,'("es20.10",es20.10)') val(2)
	WRITE(*,'("es15.7",es15.7)') val(2)	! this is the best one too because it saves space 
	WRITE(*,'("e20.4",e15.7)')  val(2)	
	WRITE(*,'("e15.8",e15.8)')  val(2)	! this is the best one also 
	!WRITE(*,'(e15.9)')  val(2)	! this doesn't work 
	!WRITE(*,*) "STOP "
	WRITE(*,*) val(2)		! for sp, mantissa can only be 8 digits so that the computer reads val as 12345679 (where 9 is for roudded off 8). so this writes 1.2345679E+08
	WRITE(*,'("F11.2",F11.2)') val(2)	! the width of the field is not large enough for 123456789 (which has 9 digits). so this writes **************
	WRITE(*,'("F12.2",F12.2)') val(2)	! width of the field is fine, but now it will fill up that 9th digit with SOME number (in this case 2) so that it will write 123456792.00
	WRITE(*,'("F13.2",F13.2)') val(2)	! width of the field is fine, but now it will fill up that 9th digit with SOME number (in this case 2) so that it will write 123456792.00 (but with more space in the beginning than the above one) 
	WRITE(*,'("F20.2",F20.2)') val(2)	! width of the field is fine, but now it will fill up that 9th digit with SOME number (in this case 2) so that it will write 123456792.00 (but with more space in the beginning than the above one) 				
	WRITE(*,*) "STOP "
	val(2)=12345678912345.
	WRITE(*,'(F12.2)') val(2)	! same explanations above (i.e. the comp will just add a different number 2. only now, the number does the width of field. 
	WRITE(*,'(F13.2)') val(2)   
	WRITE(*,'(F20.2)') val(2)   
	END SUBROUTINE precisiontests


SUBROUTINE compwage(ai,aj,gamma,w) 
	real(dp), INTENT(IN) :: ai,aj,gamma 
	real(dp), INTENT(OUT) :: w
	w = (ai - aj)** (-gamma) - gamma * ai * (ai - aj)** (-gamma-1)
	!w = 1./exp(ai - aj) - ai * (1./exp(ai-aj) )
	!w = ( exp(ai - aj) ) ** (-gamma)  *  (1.-gamma*ai) 
END SUBROUTINE 

subroutine rand2num(mu,sig,rand,num) 
	real(dp), INTENT(IN) :: rand(:),mu,sig
	real(dp), INTENT(OUT) :: num(:)
	num=sqrt(2.0_Dp) * SIG * rand(:) + MU
end subroutine rand2num

END MODULE alib

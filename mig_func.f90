MODULE FUNC
	USE params
	IMPLICIT NONE
CONTAINS
	FUNCTION fnwge(dg,dtyp,dl,dw,de,dr)						
	INTEGER(I4B), INTENT(IN) :: dg,dtyp,dl,de,dr						! gender,typ,location,education,experience
	REAL(RP), INTENT(IN) :: dw								! wage draw
	REAL(RP) :: fnwge
	if (dg==1) then 
		fnwge=exp(ALF10(dl)+ALF11*abs(de==2) + ALF12*dr + ALF13*(dr**2) + dw ) 
	else if (dg==2) then 
		fnwge=exp(ALF20(dl)+ALF21*abs(de==2) + ALF22*dr + ALF23*(dr**2) + dw ) 
	end if 
	END FUNCTION fnwge

	FUNCTION fnprof(dw0,de,dsex)
	INTEGER(I4B), INTENT(IN) :: dw0,de,dsex
	REAL(RP), DIMENSION(3) :: fnprof
	fnprof=0.0_rp
	if ( dw0 <= NP ) then 
		if ( de==1 .and. dsex==1 ) then 
			fnprof(1:2)=exp(PSIO(1:2)) !exp( PSIO(1) + PSIO(2) * abs(dsex==1) + PSIO(3) * abs(de==2) )	! offer 	 
		else if ( de==2 .and. dsex==1 ) then 
			fnprof(1:2)=exp(PSIO(3:4)) 
		else if ( de==1 .and. dsex==2 ) then 
			fnprof(1:2)=exp(PSIO(5:6)) 
		else if ( de==2 .and. dsex==2 ) then 
			fnprof(1:2)=exp(PSIO(7:8)) 
		end if 
		fnprof(3)=exp( 0.0_rp )		! nothing happens												
	else if (dw0 == NP1) then 
		if ( de==1 .and. dsex==1 ) then 
			fnprof(1)=exp(PSIO(9)) !exp( PSIO(1) + PSIO(2) * abs(dsex==1) + PSIO(3) * abs(de==2) )	! offer 	 
		else if ( de==2 .and. dsex==1 ) then 
			fnprof(1)=exp(PSIO(10) )
		else if ( de==1 .and. dsex==2 ) then 
			fnprof(1)=exp(PSIO(11)) 
		else if ( de==2 .and. dsex==2 ) then 
			fnprof(1)=exp(PSIO(12)) 
		end if 
		fnprof(2)=0.0_rp		! 0 since you can't get laid off if you don't have a job! 
		fnprof(3)=exp(0.0_rp)		! nothing happens												
	else  
		print*, "in fnprof: dw0 > NP1 which doesnt' make sense as that's a state variable " 
		stop
	end if 
	fnprof(1:3)=fnprof/sum(fnprof)
	!fnprof=0.0_rp
	!fnprof(1)=1.0_rp						
	if (skriv) then 
		if ( abs(  sum(fnprof) - 1.0_rp  ) > eps) then ; print*, "error in getfnprof : offer does not add up " , sum(fnprof) ; stop ; end if 
	end if 
	END FUNCTION fnprof

	FUNCTION fnprloc(orig)
	INTEGER(I4B), INTENT(IN) :: orig		! Origin location
	REAL(RP), DIMENSION(NL) :: fnprloc
	REAL(RP) :: sum_sans_orig 
	INTEGER(I4B) :: j
	fnprloc=0.0_rp
	fnprloc(orig)=logit(PSIL(1))
	sum_sans_orig=0.0_rp
	do j=1,NL	
		if (j /= orig) then 
			fnprloc(j)= exp( PSIL(2) * DISTANCE(j,orig) + PSIL(3) * POPSIZE(j) ) 
			sum_sans_orig = sum_sans_orig + fnprloc(j) 
		end if 
	end do 
	do j=1,NL	
		if (j /= orig) then 
			fnprloc(j)=(1.0_rp-fnprloc(orig) ) * fnprloc(j)/sum_sans_orig
		end if 
	end do 
	if ( abs(  sum(fnprloc) - 1.0_rp  ) > eps) then ; print*, "error in getfnprloc : offer does not add up " , sum(fnprloc) ; stop ; end if 
	END FUNCTION fnprloc

	FUNCTION fnprhc(dr,dw)
	INTEGER(I4B), INTENT(IN) :: dr,dw		! experience and employment: w<=NP work, w==NP1 not work,  w=NP2 nothing/can't be a state variable here so if you get this, there's something wrong
	REAL(RP), DIMENSION(NEXP) :: fnprhc
	INTEGER(I4B) :: j
	if (skriv) then 
		if ( dw > NP1 ) then ; print*, "in fnprof: dw0 > NP1 which doesnt' make sense as that's a state variable " ; stop ; end if 
	end if 
	fnprhc=0.0_rp
	do j=1,NEXP	
		if ( dw <= NP ) then 
			if (j==dr) then 
				fnprhc(j)=exp(0.0_rp)
			else if ( j-dr == -1 ) then  
				fnprhc(j)=exp(  PSIH(1)   )
			else if ( j-dr == +1 ) then  
				fnprhc(j)=exp(  PSIH(2)   )
			else 
				fnprhc(j) = 0.0_rp
			end if 
		else if ( dw == NP1 ) then 
			if (j==dr) then 
				fnprhc(j)=exp(0.0_rp)
			else if ( j-dr == -1 ) then  
				fnprhc(j)=exp(  PSIH(3)   )
			else if ( j-dr == +1 ) then  
				fnprhc(j)=exp(  PSIH(4)   )
			else 
				fnprhc(j) = 0.0_rp
			end if 
		end if 
	end do 	
	fnprhc(:)=fnprhc(:)/sum(fnprhc)
	if (skriv) then 
		if ( abs(sum(fnprhc(:))-1.0_rp) > eps ) then ; print*, " error in fnprhc: prhc does not add up " , dw , sum(fnprhc(:)) ; stop ; end if 
	end if 
	END FUNCTION fnprhc

	FUNCTION fnprkid(kid0)
	INTEGER(I4B), INTENT(IN) :: kid0
	REAL(RP), DIMENSION(0:MAXKID) :: fnprkid
	INTEGER(I4B) :: j
	fnprkid=0.0_rp
	do j=kid0,MAXKID
		fnprkid(j)=exp(  PSIK * (j-kid0) )
		if ( abs(j-kid0) > 1 ) then 
			fnprkid(j)=0.0_rp	!can only move one step up or down
		end if 
	end do 						
	fnprkid(0:MAXKID)=fnprkid(0:MAXKID)/sum(fnprkid(0:MAXKID))
	if (skriv) then 
		if ( abs(sum(fnprkid(0:MAXKID))-1.0_dp)>eps ) then ; print*, "error in fnprkid: prkid does not add up " , kid0 , sum(fnprkid(0:MAXKID)) ; stop ; end if 
	end if 
	END FUNCTION fnprkid

	FUNCTION fnmove(kid) 
	INTEGER(I4B), INTENT(IN) :: kid
	REAL(RP) :: fnmove
	fnmove = CST + KCST * abs(kid>0)
	!fnmove = fnmove / DIV
	END FUNCTION fnmove

	SUBROUTINE q2wloc(dq,dw,dl)
	! extract indeces w,l from q 
	INTEGER(I4B), INTENT(IN) :: dq		
	INTEGER(I4B), INTENT(OUT) :: dw,dl
	INTEGER(I4B), DIMENSION(2) :: indeces	
		indeces=lin2ndim( (/ NP2 , NL /) , dq )
		dw=indeces(1)
		dl=indeces(2)
		if (skriv) then  
			if ( dq > NQS ) then ; print*, "q2wl: q > NQS", dq, NQS,indeces ; stop ; end if  
			if ( dw > NP2 ) then ; print*, "q2wl: w > NP2" ; stop ; end if  
			if ( dl > NL  ) then ; print*, "q2wl: l > NL" ; stop ; end if  
		end if 
	END SUBROUTINE
	SUBROUTINE wloc2q(dq,dw,dl)
	! construct combined q from w,l
	INTEGER(I4B), INTENT(OUT) :: dq		
	INTEGER(I4B), INTENT(IN) :: dw,dl		
    	dq = ndim2lin( (/ NP2 , NL /),(/ dw,dl /) )
		if (skriv) then 		
			if ( dq > NQS ) then ; print*, "wl2q: q > NQS" ; stop ; end if  
			if ( dw > NP2 ) then ; print*, "wl2q: w > NP2" ; stop ; end if  
			if ( dl > NL  ) then ; print*, "wl2q: l > NL" ; stop ; end if  
		end if 
	END SUBROUTINE

	SUBROUTINE x2edexp(dx,de,dr)
	! extract indeces educ,experience from x
	INTEGER(I4B), INTENT(IN) :: dx		
	INTEGER(I4B), INTENT(OUT) :: de,dr
	INTEGER(I4B), DIMENSION(2) :: indeces	
		indeces=lin2ndim( (/ NEDUC, NEXP /) , dx )
		de=indeces(1)
		dr=indeces(2)
	end SUBROUTINE
	SUBROUTINE edexp2x(dx,de,dr)
	!construct combined x from educ,experience
	INTEGER(I4B), INTENT(OUT) :: dx		
	INTEGER(I4B), INTENT(IN) :: de,dr
		dx=ndim2lin( (/ NEDUC, NEXP /),(/ de,dr /) )
	END SUBROUTINE

	SUBROUTINE index2cotyphome(index,co,typ,home)
	! extract indeces cohort,type,educ,homeloc from combined index
	INTEGER(I4B), INTENT(in) :: index		
	INTEGER(I4B), INTENT(out) :: co,typ,home	
	INTEGER(I4B), DIMENSION(3) :: indeces	
	if (groups) then 
		!indeces=lin2ndim((/NCO,NTYP,NL/),index)
		co=1 !indeces(1)
		typ=1 !indeces(2)
		home=1 !indeces(3)
	else 
		indeces=lin2ndim((/NCO,NTYP,NL/),index)
		co=indeces(1)
		typ=indeces(2)
		home=indeces(3)
	end if 
	END SUBROUTINE index2cotyphome
	SUBROUTINE cotyphome2index(index,co,typ,home)
	!construct combined index from co,typ,home
	INTEGER(I4B), INTENT(out) :: index		! combined index
	INTEGER(I4B), INTENT(in) :: co,typ,home 
	if (groups) then 
		index=1 !ndim2lin((/NCO,NTYP,NL/),(/co,typ,home/))
	else 
		index=ndim2lin((/NCO,NTYP,NL/),(/co,typ,home/))
	end if 
	END SUBROUTINE cotyphome2index
END MODULE FUNC


!ahu 083012 REAL(SP), dimension(size(Vmsepcho,1),size(Vmsepcho,2)) :: Vmsepchot,Vfsepchot 
!if (sts%kid/=0) then ; print*, "Error in getdecrule_s. There are kids here ",sts%kid ; END if  this is ok since sometimes this is called by getdivchoice
!if (h==1.and.offer==indlayoff) then ; print*, "Error! hprev is 1 and offer is indlayoff! " ; stop ; END if 
!fortran tip: time : you don't need this following assignment to penalty if statement since you assign penalty values to the relevant elements of these arrays in the subsequent do loop in any case
!if (sex==M) then ; Vmsepchot(:,:)=penalty !Vmsepchot(:,:)=penalty this takes the same time as Vmsepchot=penalty ....
!no need for these and they take time. maxchos=initi ; maxvs=init

!ahu 083012 REAL(SP), dimension(size(Vmcho,1),size(Vmcho,2),size(Vmcho,3)) :: Vmchot  !,nashprod 
!ahu 083012 REAL(SP), dimension(size(Vfcho,1),size(Vfcho,2),size(Vfcho,3)) :: Vfchot 

!no need for these and it takes time. they get assigned in any case. maxcho=initi ; maxv=init ; mufin=initi 
!fortran tip: time: you don't need this following assignment to penalty since you assign penalty values to the relevant elements in the subsequent do loop in any case 
!Vmchot=penalty ; Vfchot=penalty ; nashprod=penalty !Vmchot(:,:,:)=penalty ; Vfchot(:,:,:)=penalty this takes the same time as Vmchot=penalty
!when the dimensions of Vmchot was (loccho,hmcho,hfcho,ieps,mui,i3,theta)
!Vmchot=penalty 
!takes much longer time than 
!Vmchot(:,:,:,,ieps,mui,i3,theta)=penalty 
!but when the dimensions of Vmchot is (loccho,hmcho,hfcho) 
!Vmchot=penalty 
!takes the same time as 
!Vmchot(:,:,:)=penalty 
!musum=penalty !!!ahu 091912



!SUBROUTINE getshares(vce,vs,PI,SHARE)	
!	REAL(SP), INTENT(in) :: vce(2,nce),vs(2)
!	REAL(SP), INTENT(out) :: PI(nce),SHARE(2,nce)
!	REAL(SP) :: vdif(2),sum
!	INTEGER(I4B) :: cho,qcho

!	do cho=1,nce
!		qcho=chovec(cho,iq,iq0)						
!		if (qcho>0) then 
!			vdif(:)=vce(:,cho)-vs(:)
!			pi(cho)=totinc(qcho,ix)
!			if (  pi(cho)-abs(vdif(2)-vdif(1)) >  eps   ) then		!Interior Solution
!				sum= pi(cho) + vdif(2) - vdif(1)
!				SHARE(1,cho)=sum*0.5d0/pi(cho) 
!				SHARE(2,cho)=1d0-SHARE(1,cho)
!			else if (pi(cho)-(vdif(1)-vdif(2)) < eps ) then		!Corner Solution where she gets to consume all income
!				SHARE(:,cho)=0d0
!			else if (pi(cho)-(vdif(2)-vdif(1)) < eps ) then		!Corner Solution where he gets to consume all income
!				SHARE(:,cho)=1d0
!			END if 
!		else 
!			SHARE(:,cho)=penalty
!			pi(cho)=penalty
!		END if 
!	END do 
!END SUBROUTINE getshares
!SUBROUTINE getdecmar(sex,qs,xs,ia,index,ie,CHECKNB,MAXV)
!	INTEGER(I4B), INTENT(in) :: sex,qs(2),xs(2),ia,index,ie
!	type(choice), INTENT(out) :: MAXC
!	REAL(SP), INTENT(out) :: VAL(4)
!	INTEGER(I4B) :: iq,ix,ie,iqm,ixm,iqf,ixf
	
!	do i=1,2
!		vs(i)=vsep(i,qs(i),xs(i),ia,index)
!		vc=vmar(i,iq,ix,ia,index)+mariegrid(ie)
!	END do 
!	vdif=vc-vs
!	pi=totinc(iq,ix)
!	call getnashprod(vdif,pi,CHECKNB,PROD)
!	if (checknb>0) then 
!		MAXC%ic=0
!		MAXC%iq=iq
!		MAXC%ix=ix
!	else	
!		MAXC=choice0
!	END if 
!END SUBROUTINE getdecmar

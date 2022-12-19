! information about the value functions
! get vmsep/vfsep=values of being single for males
! this is the value of remaining single, either because you started the period single and you reject the partner you meet or because you leave a relationship at
! the start of the period. not the value of starting period single, as it does not include the value of potentially entering into a relationship		  			
!get each spouse's consumption transfers for each j alternative using foc to the nb problem to get:
!share(2,cho): spouses' respective consumption shares of total income (pi) for each alternative cho
!these share solutions are conditional on the nb problem being well-defined. otherwise they don't mean anything
!i check for whether the nb problem is well-defined in the getnb subroutine later because i wait for marie
! checknb: check if the nb problem is well-defined meaning that there exists some allocation in the utility possibility set that is mutually beneficial
	!testing(1,1)=1 ; testing(1,2)=2 ; testing(1,3)=3 ; testing(2,1)=4 ; testing(2,2)=5 ; testing(2,3)=6
	!print*, "maxloc(testing,1),maxloc(testing,2)",maxloc(testing,1),maxloc(testing,2),"maxval",maxval(testing(2,:),1)   !,maxval(testing,2)
	!print*, nqs*nqs,nq*nq,count(pps(:,:,1)),count(ppc(:,:))
module sol
	use params
	use share
	use myaz
	implicit none
	real(sp) :: begintime,time(5)   
    real(dp) :: vsingtest(2) !ag090822 agsept2022
contains		
	subroutine solve
	real(dp) :: vmax(2),val(2),vsum(2),valso(2),probmuq(nx,nq,nq),probmux(nx,nx,nq)
    real(dp), dimension(nxs,nqs) :: vm0_s,vf0_s,vm_s,vf_s,prob_s,vcontm_s,vcontf_s
	real(dp), dimension(nx,nq) :: vm_c,vf_c,prob,vcontm_c,vcontf_c
	real(dp), dimension(nepsmove,nxs,nqs,nqs) :: vmr,vfr
	integer(i4b) :: ia,q0,q,x0,x,z,iepsmove,index,trueindex,dd(12),ed(2),qm,xm,qf,xf,callfrom
    !print*, 'Here is iter',iter
	begintime=secnds(0.0)
    yaz=.false.
    if (groups) then 
        call cotyphome2index(myindex,myco,mytyp,myhome)
    else 
        myindex=0
    end if 
    do x0=1,nx
    do q0=1,nq
        do q=1,nq
        probmuq(x0,q,q0)= ppcq(q,q0,x0)
        end do 
        do x=1,nx
        probmux(x0,x,q0)=ppcx(x,q0,x0) !ppcq(q,q0,x0)*ppcx(x,q0,x0)
        end do 
    end do 
end do 


    call get_util_w
    if (onlysingles) then ; utilc=pen ; end if 
	insol=.true.
	!prob=0.0_dp
	decm0_s = ipen	; decf0_s = ipen  ; decm_s = ipen	; decf_s = ipen ; dec_mar=ipen ! dum is a dummy for yaz_checknb
    vm    = pen ; vf   = pen ; vm0_c = pen ; vf0_c = pen
    do trueindex=1,ninp
        if (groups) then 
            index=1
        else 
            index=trueindex
        end if
        ind: if ( (.not.groups) .or.  (myindex==trueindex) ) then !(myco==init(nn)%co.and.mytyp==typsim.and. myhome==init(nn)%hme)  ) then	  
		time(1)=secnds(0.0)
		emaxm_s=0.0_dp	; emaxf_s=0.0_dp ; emaxm_c=0.0_dp	; emaxf_c = 0.0_dp
		do ia=mxa,mna,-1
			!print*, "mysay,ia,trueindex ",iter,mysay,ia,trueindex !ahu 030622
			vmr  = pen ; vfr = pen	
			time(1)=secnds(0.0)
			if (ia==mxa) then 
				vm0_s = utils(1,:,:,trueindex) + wsnet(1,:,:,trueindex)		! v0: value function without movecost 
				vf0_s = utils(2,:,:,trueindex) + wsnet(2,:,:,trueindex)		! v0: value function without movecost
				vm0_c(:,:,ia,index) = utilc(1,:,:,trueindex)	! v0: value function without movecost, umar, consumption
				vf0_c(:,:,ia,index) = utilc(2,:,:,trueindex)	! v0: value function without movecost, umar, consumption
                if (terminalval) then 
                    vcontm_s=vm0_s
                    vcontf_s=vf0_s
                    vcontm_c=vm0_c(:,:,ia,index)
                    vcontf_c=vf0_c(:,:,ia,index)
                    vm0_s=vm0_s+delta*vcontm_s+(delta**2)*vcontm_s+(delta**3)*vcontm_s+(delta**4)*vcontm_s+(delta**5)*vcontm_s+(delta**6)*vcontm_s+(delta**7)*vcontm_s+(delta**8)*vcontm_s+(delta**9)*vcontm_s+(delta**10)*vcontm_s
                    vf0_s=vf0_s+delta*vcontf_s+(delta**2)*vcontf_s+(delta**3)*vcontf_s+(delta**4)*vcontf_s+(delta**5)*vcontf_s+(delta**6)*vcontf_s+(delta**7)*vcontf_s+(delta**8)*vcontf_s+(delta**9)*vcontf_s+(delta**10)*vcontf_s
                    vm0_c(:,:,ia,index)=vm0_c(:,:,ia,index)+delta*vcontm_c+(delta**2)*vcontm_c+(delta**3)*vcontm_c+(delta**4)*vcontm_c+(delta**5)*vcontm_c+(delta**6)*vcontm_c+(delta**7)*vcontm_c+(delta**8)*vcontm_c+(delta**9)*vcontm_c+(delta**10)*vcontm_c
                    vf0_c(:,:,ia,index)=vf0_c(:,:,ia,index)+delta*vcontf_c+(delta**2)*vcontf_c+(delta**3)*vcontf_c+(delta**4)*vcontf_c+(delta**5)*vcontf_c+(delta**6)*vcontf_c+(delta**7)*vcontf_c+(delta**8)*vcontf_c+(delta**9)*vcontf_c+(delta**10)*vcontf_c                                
                end if 
            else 
				vm0_s = utils(1,:,:,trueindex) + wsnet(1,:,:,trueindex)		+ delta * emaxm_s(:,:,ia+1)	! v0: value function without movecost 
				vf0_s = utils(2,:,:,trueindex) + wsnet(2,:,:,trueindex)		+ delta * emaxf_s(:,:,ia+1)	! v0: value function without movecost 
				vm0_c(:,:,ia,index) = utilc(1,:,:,trueindex)	+ delta * emaxm_c(:,:,ia+1)	! v0: value function without movecost, umar, consumption
				vf0_c(:,:,ia,index) = utilc(2,:,:,trueindex)	+ delta * emaxf_c(:,:,ia+1)	! v0: value function without movecost, umar, consumption
			end if 		
            
			
			!if (skriv.and.trueindex==1.and.(ia<=19.or.ia==45)) then ; yaz=.true. ; else ; yaz=.false. ; end if !ahu 0327 trueindex==2
            dd = -1 ; dd(1:2) = (/ia,trueindex/) 
			dd(7)=1 ; call getdec_s(dd,vm0_s,decm0_s(:,:,:,:,ia,index),vm(:,:,:,:,ia,index)) 
			yaz=.false.
			dd(7)=2 ; call getdec_s(dd,vf0_s,decf0_s(:,:,:,:,ia,index),vf(:,:,:,:,ia,index))	

			if (onlysingles) then 	
				decm_s(:,:,:,:,ia,index)=decm0_s(:,:,:,:,ia,index)
				decf_s(:,:,:,:,ia,index)=decf0_s(:,:,:,:,ia,index)
				vmr=vm(:,:,:,:,ia,index)
				vfr=vf(:,:,:,:,ia,index) 
			else 
                !for loop opt 
                do q=1,nq
                    do x=1,nx
                        vm0ctemp(q,x)=vm0_c(x,q,ia,index) 
                        vf0ctemp(q,x)=vf0_c(x,q,ia,index) 
                        wmctemp(q,x,trueindex)=wcnet(1,x,q,trueindex)
                        wfctemp(q,x,trueindex)=wcnet(2,x,q,trueindex)
                    end do 
                end do 

				do q0=1,nq	
					vm_c=pen ; vf_c=pen
					prob=0.0_dp
					do q=1,nq
                    do x=1,nx	
						
                            !ahumarch1122 
                            if (skriv.and.(ia==50.or.ia==48.or.ia==18).and.q0==4.and.q==92.and.x==19.and.trueindex==1) then ; yaz=.true. ; else ; yaz=.false. ; end if !ahu030822
                            !ahu030822 Make sure you don't have iepsmove in the if statements here

							if ( ppc(q,q0) ) then
                                valso=pen
								val=0.0_dp
								umar: do z=nz,1,-1		
                                    moveshocks: do iepsmove=1,nepsmove
                                        !callfrom is for telling getdec and therefore yaz_getdec where we are when yaz_getdec is called from within getdec
                                        !(if skriv and yaz true) when callfrom=40, then within getdec I call yaz_getdec to write in 200 the altspecific value functions and bargaining stuff
                                        !in getdec when callfrom is 40 or 80, then dd(8) altj and dd(9) altq are filled in within the choice loop
                                        callfrom=40 !ag090122 agsept2022
                                        dd = (/ia,trueindex,q,x,z,q0,callfrom,-1,-1,-1,iepsmove,-1 /) 	! (/ ia,index,q,x,z,q0,callfrom,altj/jmax,altq/qmax,relmax,iepsmove,gender /)  	                                        
                                        call getdec(dd,vmax,valso)	! dd is inout and dd(8:10) tells you jmax,qmax,relmax 
                                        if (yaz) then ; call yaz_decision(dd,vmax) ; end if	!callfrom tells yaz_decision where we are 									
		                                val=val+ppso(iepsmove)*mgt(z)*vmax
                                    end do moveshocks
                                end do umar
								vm_c(x,q)=val(1) ; vf_c(x,q)=val(2)
							end if !pc0(q0)=.true. i.e. w(1) or w(2) are not np2 and pc
						end do !x
					end do !q
					!if (icheck_eqvmvf) then
					!	call check_eqvmvf(dd,vm_c,vf_c) 
					!end if 


                        emaxm_c(:,q0,ia)=0.0_dp
						emaxf_c(:,q0,ia)=0.0_dp
                        if ( maxval(qq2w(:,q0)) <= np1 ) then !state variable part of the q space i.e. w <= np1
						    !prob=matmul( reshape(ppcq(:,q0,x0),(/nq,1/)) , reshape(ppcx(:,q0,x0),(/1,nx/)) )
						    !emaxm_c(x0,q0,ia)=sum(prob*vm_c)
						    !emaxf_c(x0,q0,ia)=sum(prob*vf_c)                    
                            do q=1,nq
                            do x=1,nx
                            do x0=1,nx
                                    emaxm_c(x0,q0,ia)=emaxm_c(x0,q0,ia)+probmuq(x0,q,q0)*probmux(x0,x,q0)*vm_c(x,q)
						            emaxf_c(x0,q0,ia)=emaxf_c(x0,q0,ia)+probmuq(x0,q,q0)*probmux(x0,x,q0)*vf_c(x,q)
                            end do 
                            end do 
                            end do	 				
                        end if !state spce check 
                end do !q0
				!if (Trueindex==1) print*, "I even came here and nothing happened"

                !timeline: when single, you first decide on loc and work. and then when the day comes to an end, you go to the marriage market in whatever location you chose for that period.
				vm_s=0.0_dp ; vf_s=0.0_dp 
				do qm=1,nqs	
                do xm=1,nxs					
                    ed(1)=x2e(xm)   !ahu summer18 051418: adding ed dimension to nonlabinc
						!no need for this ahu040518 vsingle=0.0_dp
						do qf=1,nqs
                        do xf=1,nxs			
                            ed(2)=x2e(xf)   !ahu summer18 051418: adding ed dimension to nonlabinc
							
                                !i = dd(8) ; n=dd(9) 
                                valso(1) = vm0_s(xm,qm)
                                !i = dd(10) ; n=dd(11) 
                                valso(2) = vf0_s(xf,qf)
                                vsingtest(1:2)=valso(1:2) !ag090822
                    
                                if (skriv.and.(ia==18).and.(qm==18).and.xm==1) then ; yaz=.true. ; else ; yaz=.false. ; end if     !ahu 0327 trueindex==2								
                                !no need for this ahu040518 sex: do g=1,2
                                !no need for this ahu040518 if (g==1) then ; qm=qr ; qf=qp ; xm=xr ; xf=xp ; end if ! determine the match and the joint q from the combination
                                !no need for this ahu040518 if (g==2) then ; qm=qp ; qf=qr ; xm=xp ; xf=xr ; end if ! determine the match and the joint q from the combination
                                if (ppmeetq(qm,qf)>0.0_dp .and. ppmeetx(xm,xf)>0.0_dp) then 
                                    q0=-1	; q = q2qq(qm,qf)  ; x = x2xx(xm,xf) 
                                
                                    if (q2w(qm)>np1.or.q2w(qf)>np1.or.q==0.or.x==0) stop   !if (ps0(qm) .and. ps0(qf) .and. q > 0 .and. x >0) then ! all w should be <= np1	and q,x should exist
                                    marshock: do z=1,nz  
                                    callfrom=50
                                    dd = (/ ia,trueindex,q,x,z,q0,callfrom,-1,q,-1,-1,-1 /) ! (/ ia,index,q,x,z,q0,callfrom,altj,altq,rel,iepsmove,gender /)  	                                        
                                    !callfrom is for telling getdec and therefore yaz_getdec where we are when yaz_getdec is called from within getdec
                                    !(if skriv and yaz true) when callfrom=50, then within getdec I call yaz_getdec to write in 200 the altspecific value functions and bargaining stuff
                                    !in getdec when callfrom is 40 or 80, then dd(8) altj and dd(9) altq are filled in within the choice loop
                                    !in getdec when callfrom is 50, then dd(8) altj is -1 and dd(9) altq is just q (same as dd(3) ) 
                                                        ! g and q0 is just -1 here (/ ia,index,q,x,z,q0,gender,j,q,rel,iepsmove /) 
                                                        ! when in the marriage market, the index that corresponds 
                                                        ! to altq (i.e. i in yaz for example or dim 9 in dd) is just set to q since there is no altq in marriage market
                                    call getdec(dd,vmax,valso) 
                                    dec_mar(z,x,q,ia,index)=dd(10)
                                    val=vmax(1:2)         
                                    if (yaz) then ; call yaz_decision(dd,vmax) ; end if	!callfrom tells yaz_decision where we are 				

                                    !if (  icheck_eqvmvf.and.qr==qp.and.xr==xp.and.  (abs(vec(1)-vec(2)) > eps .or. abs(vsum(1)-vsum(2)) > eps) ) then 
                                    !	print*, "vm not eq to vf!", ia,qr,xr,vec(1:4),vsum
                                    !	stop 
                                    !end if 
                                    !no need for this ahu040518 vsingle(g,qp,xp,z)=val(g)	
                                    vm_s(xm,qm)=vm_s(xm,qm) + mgt(z) * ppmeetq(qm,qf) * ppmeetx(xm,xf) * val(1)
                                    vf_s(xf,qf)=vf_s(xf,qf) + mgt(z) * ppmeetq(qm,qf) * ppmeetx(xm,xf) * val(2)
										end do marshock 
									end if 
							        !no need for this ahu040518 end do sex 
							        !if (icheck_eqvmvf) then 
							        !	do z=1,nz
							        !		if ( abs(vsingle(1,qp,xp,z) - vsingle(2,qp,xp,z))>eps ) then 
							        !			write(*,'("vsingle not equal!!!")')
							        !			stop
							        !		end if 
							        !	end do 
							        !end if 
							end do 
						end do 							
						!moving this after the loop ahu040518 vm_s(xr,qr) = pmeet*vm_s(xr,qr) + (1.0_dp-pmeet)*vm0_s(xr,qr)
						!moving this after the loop ahu040518 vf_s(xr,qr) = pmeet*vf_s(xr,qr) + (1.0_dp-pmeet)*vf0_s(xr,qr)
						!if (  icheck_eqvmvf.and.abs(vm_s(xr,qr)-vf_s(xr,qr)) > eps ) then 
						!	print*, "vm not equal to vf!", ia,qr,xr,vm_s(xr,qr),vf_s(xr,qr)
						!	stop 
						!end if 
					end do 
				end do
				vm_s(:,:) = pmeet*vm_s(:,:) + (1.0_dp-pmeet)*vm0_s(:,:) !ahu 040518
				vf_s(:,:) = pmeet*vf_s(:,:) + (1.0_dp-pmeet)*vf0_s(:,:) !ahu 040518
                !vm_s=vm0_s
                !vf_s=vf0_s
                !dec_mar=0
				yaz=.false.
				if (skriv.and.(ia==mxa.or.ia==29).and.trueindex==1) yaz=.false.
				dd = -1 ; dd(1:2) = (/ia,trueindex/) 							! (/ ia,index,q,x,z,q0,g,j,altq,rel,iepsmove /)
				dd(7)=1 ; call getdec_s(dd,vm_s,decm_s(:,:,:,:,ia,index),vmr(:,:,:,:)     ) 
				dd(7)=2 ; call getdec_s(dd,vf_s,decf_s(:,:,:,:,ia,index),vfr(:,:,:,:)     ) 			
			end if 	! only singles or no 

			do q0=1,nqs
                if (q2w(q0)<=np1) then !state variable part of the q space i.e. w <= np1
                do x0=1,nxs			
						do q=1,nqs
                            do x=1,nxs 
                                do iepsmove=1,nepsmove
                                    !prob_s=matmul( reshape(ppsq(:,q0,x0,1),(/nqs,1/)) , reshape(ppsx(:,q0,x0),(/1,nxs/)) )
						            !emaxm_s(q0,x0,ia)	= sum( prob_s * vmr(:,:,q0) )		
                                    !prob_s=matmul( reshape(ppsq(:,q0,x0,2),(/nqs,1/)) , reshape(ppsx(:,q0,x0),(/1,nxs/)) )
						            !emaxf_s(q0,x0,ia)	= sum( probf_s * vfr(:,:,:,q0) )
                                    emaxm_s(x0,q0,ia)	= emaxm_s(x0,q0,ia) + ppso(iepsmove) * ppsq(q,q0,x0,1) * ppsx(x,q0,x0) * vmr(iepsmove,x,q,q0) 
                                    emaxf_s(x0,q0,ia)	= emaxf_s(x0,q0,ia) + ppso(iepsmove) * ppsq(q,q0,x0,2) * ppsx(x,q0,x0) * vfr(iepsmove,x,q,q0) 	
                                end do 
                            end do 
                        end do
				end do 
                    end if !state variable part of the q space i.e. w <= np1
			end do 

			yaz=.false.			
			!if (skriv) call yaz1(index,ia) !ahu 0327. this was trueindex and corrected it to index.  
		end do ! age 
        end if ind
	end do !index
	!if ( (.not.chkstep).and.(.not.optimize) ) print*, "Finished solve and mysay is: ", mysay
	if (onlysingles) then
		dec_mar=0
	end if 
	insol=.false.
	end subroutine solve



    !ag090122 agsept2022 Changed getdec_c to getdec in order to include mar market decisinos in getdec as well
    !                     so that everything is all in one same place and if decisino protocol changes are made they have to be made to only one place
	!                     OLDER VERSIONS OF GETDEC_C CAN BE FOUND IN MIGSOLGETDEC.F90 FILE
    subroutine getdec(dd,vmax,vsing) 
        integer(i4b), intent(inout) :: dd(:)
        real(dp), intent(out) :: vmax(2) 
        real(dp), intent(in) :: vsing(2) 
        real(dp) :: vecj(5,nc),vec(5),surplusj(nc),surplus,transfers(2),vdif(2),mcost(2),asum,yazvec(5)
        integer(i4b) :: ia,index,q,x,z,q0,g,jmax,qmax,relmax,iepsmove,i,i0,n,trueindex,j,de(1),ed(2),locch,loc0,callfrom
        logical :: haveenough(nc),haveenoughtomakeindiff,haveenoughtomakeindiff_alt,intsol(3),pc(2),pc_alt(2)
        !dd = (/ia,trueindex,q,x,z,q0,gender/callfrom,-1,-1,-1,        iepsmove /) 	 
        !dd=  (/ia,trueindex,q,x,z,q0,gender/callfrom,jmax,qmax,relmax,iepsmove /)  	
        callfrom=dd(7)
        ia=dd(1) 
        trueindex=dd(2) 
        if (groups) then 
            index=1
        else 
            index=trueindex
        end if
        q=dd(3) 
        x=dd(4) ; ed(:)=xx2e(:,x)     
        z=dd(5) 
        q0=dd(6)  
        vec=pen				! initialize
        jmax=-1				! initialize
        qmax=-1				! initialize
        relmax=-1            ! initialize
        vmax=pen			! initialize
        transfers=pen 
        if (Callfrom==40.or.callfrom==80) then !calling from sol/couples/getdec (40) OR simul/couples/getdec (80) 
            iepsmove=dd(11)   
            vec=pen
            i = qq2q(1,q) ; n=xx2x(1,x) ; i0 = qq2q(1,q0) 
            vec(1) = vm(iepsmove,n,i,i0,ia,index) + divpenalty
            mcost(1)=movecost(n,i0,trueindex)  
            i = qq2q(2,q) ; n=xx2x(2,x) ; i0 = qq2q(2,q0) 	
            vec(2) = vf(iepsmove,n,i,i0,ia,index) + divpenalty	
            mcost(2)=movecost(n,i0,trueindex)
            loc0=qq2l(1,q0)     !ahumarch2022 ahu032022
            !ahu030822 adding the below in order to check the mar rates decreasing with mumar situation to compare chk2's better
            !if (yazmax) then !ahu030822
                !write(201,'("********************************")')
                !write(201,'(6x,tr5,"age",tr6,"q0",tr7,"q",tr7,"z",tr7,"j",tr4,"altq",tr1,"iepsmve",tr1,"trueind",tr1,"sex",TR3,"X",tr3,"def")' ) !ahu030822
                    !write(201,'(6x,8i8,i4,I4,L6)')	ia,q0,q,z,j,i,iepsmove,trueindex,-1,X,defj(j)    ! altq is just the q that altrnative j corresponds to       !ahu030822
            !end if !ahu030822
            vecj=pen ; surplusj=pen ; haveenough=.FALSE. ; yazvec=pen
            choice: do j=1,nc	
                i = ch(j,q,q0)	!alternative q
                if (i>0 ) then		
                    locch=qq2l(1,i)     !ahumarch2122 ahu032122
                    !if (qq2l(1,i).ne.qq2l(2,i)) then ; print*, 'something wrong in dec_c' ; stop; end if 
                    !if (qq2l(1,q0).ne.qq2l(2,q0)) then ; print*, 'something wrong in dec_c' ; stop; end if 
                    !ahumarch2122 ahu032122 vec(3) = vm0_c9726(i,x,ia,index)  + mg(z,trueindex) + one( qq2l(1,i) /= qq2l(1,q0)) * mcost(1) + one( qq2l(1,i) /= qq2l(1,q0)) * moveshock_m(iepsmove)  !fnmove(kid)     
                    !ahumarch2122 ahu032122 vec(4) = vf0_c(i,x,ia,index)  + mg(z,trueindex) + one( qq2l(2,i) /= qq2l(2,q0)) * mcost(2) + one( qq2l(2,i) /= qq2l(2,q0)) * moveshock_f(iepsmove)  !fnmove(kid)     
                    if (callfrom==40) then !calling from sol
                    vecj(3,j) = vm0ctemp(i,x)  + mg(z,trueindex) + one( locch /= loc0 ) * (mcost(1) + moveshock_m(iepsmove) ) !ahumarch2022 ahu032022     
                    vecj(4,j) = vf0ctemp(i,x)  + mg(z,trueindex) + one( locch /= loc0 ) * (mcost(2) + moveshock_f(iepsmove) ) !ahumarch2022 ahu032022
                    else if (callfrom==80) then !calling from simulation
                    vecj(3,j) = vm0_c(x,i,ia,index)  + mg(z,trueindex) + one( locch /= loc0 ) * (mcost(1) + moveshock_m(iepsmove) ) !ahumarch2022 ahu032022     
                    vecj(4,j) = vf0_c(x,i,ia,index)  + mg(z,trueindex) + one( locch /= loc0 ) * (mcost(2) + moveshock_f(iepsmove) ) !ahumarch2022 ahu032022
                    end if 
                    !vec(5) = wc(1,x,i,trueindex) + wc(2,x,i,trueindex) + one( qq2w(1,q0)<=np )* ubc(1,x,i,trueindex) + one( qq2w(2,q0)<=np )* ubc(2,i,x,trueindex) + nonlabinc + nonlabinc 	 !ahu summer18 050318: added the ubc
                    vecj(5,j) = wmctemp(i,x,trueindex) + wfctemp(i,x,trueindex) + nonlabinc(ed(1)) + nonlabinc(ed(2)) 
                    yazvec(1:5)=vecj(1:5,j) !just for yaz
                    vecj(1,j)=vecj(3,j)-vec(1)   !ahumarch1522 ahu031522 adding cornersol
                    vecj(2,j)=vecj(4,j)-vec(2)   !ahumarch1522 ahu031522 adding cornersol
                    !ahumarch2122 ahu032122 replacing this with vdif for saving time surplusj(j) = vec(5) + vec(3) - vec(1) + vec(4) - vec(2)  
                    surplusj(j) = vecj(5,j) + sum(vecj(1:2,j))  !ahumarch2122 ahu032122 replacing with vdif for saving time    
                    pc(1:2)	= ( vecj(1:2,j) + eps >= 0.0_dp )
                    asum = sum(  one(.not. pc)  *   abs( vecj(1:2,j) )   ) 
                    haveenough(j)=(  vecj(5,j) + eps - asum  >= 0.0_dp  )
                    if (yaz) then 
                        dd(8)=j ; dd(9)=i !to tell yaz_getdec which alternative is being evaluated (where j is the choice and i is the corresponding q)
                        call yaz_getdec(dd,yazvec(1:5),surplusj(j),pc(1:2),asum,haveenough(j))
                    end if
                end if      
            end do choice    
            de=maxloc(surplusj,MASK=haveenough) 
            if (de(1)>0) then 
                jmax=de(1) ; qmax=ch(jmax,q,q0) ; relmax=1
                surplus=surplusj(jmax)                
                vdif(1:2)=vecj(1:2,jmax)    
                vec(3:5)=vecj(3:5,jmax)
                if (yaz) then ; write(200,*)  ; write(200,*) ; write(400,*) ; write(400,*) ; end if
                if (yaz.and.callfrom==40) write(200,*) "in sol/getdec couples de(1)>0: jmax,qmax,relmax,surplusJ(jmax),vecj(1:5,jmax) NOW GO TO 2013"
                if (yaz.and.callfrom==40) write(200,'(I4,I8,I4,6F9.2)') jmax,qmax,relmax,surplusj(jmax),vecj(1:5,jmax)
                if (yaz.and.callfrom==80) write(400,*) "in sim/getdec couples de(1)>0: jmax,qmax,relmax,surplusJ(jmax),vecj(1:5,jmax) NOW GO TO 2013"
                if (yaz.and.callfrom==80) write(400,'(I4,I8,I4,6F9.2)') jmax,qmax,relmax,surplusj(jmax),vecj(1:5,jmax)
                GO TO 2013
            else 
                jmax=-1 ; qmax=-1  ; relmax=0 ; vmax(1:2)=vec(1:2)
                if (yaz) then ; write(200,*)  ; write(200,*) ; write(400,*) ; write(400,*) ; end if
                if (yaz.and.callfrom==40) write(200,*) "in sol/getdec couples de(1)<=0: jmax,qmax,relmax,de(1) NOW GO TO 2017"
                if (yaz.and.callfrom==40) write(200,*) jmax,qmax,relmax,de(1)
                if (yaz.and.callfrom==80) write(400,*) "in sim/getdec couples de(1)<=0: jmax,qmax,relmax,de(1) NOW GO TO 2017"
                if (yaz.and.callfrom==80) write(400,*) jmax,qmax,relmax,de(1)
                GO TO 2017
            end if 
        end if 

        if (callfrom==50) then !calling from sol/marmkt (this getdec instance is only called in sol (not sim because I save all decisions to dec array)
            iepsmove=-1 
            vec(1:2)=vsingtest(1:2) !ag090822 agsept2022 vsing(1:2) 
            vec(3) = vm0ctemp(q,x) + mg(z,trueindex) 
            vec(4) = vf0ctemp(q,x) + mg(z,trueindex) 
            vec(5) = wmctemp(q,x,trueindex) +  wfctemp(q,x,trueindex) + nonlabinc(ed(1)) + nonlabinc(ed(2))                                                    
            vdif(1)=vec(3)-vec(1)   !ahumarch1522 ahu031522 adding cornersol
            vdif(2)=vec(4)-vec(2)   !ahumarch1522 ahu031522 adding cornersol
            !ahumarch2122 ahu032122 replacing this with vdif for saving time surplus = vec(5) + vec(3) - vec(1) + vec(4) - vec(2)  
            surplus = vec(5) + sum(vdif(1:2))  !ahumarch2122 ahu032122 replacing with vdif for saving time
            pc(1:2)	= ( vdif + eps >= 0.0_dp )	!pc(1:2)    = ( vec(3:4) - vec(1:2) >= 0.0_dp )						        
            !ahu032122 ahumarch2122 commenting out to save time pc_alt(1:2)=( vdif >= 0.0_dp )	
            asum = sum(  one(.not. pc)  *   abs( vdif )   ) 
            haveenoughtomakeindiff=(  vec(5) + eps - asum  >= 0.0_dp  )
            de(1)=one( surplus+eps > 0.0_dp .and. haveenoughtomakeindiff )             
            if (yaz) then 
                !callfrom=50 is sol/marmkt/getdec and before calling getdec, dd(9) is just set to q i.e. equal to dd(3) and dd dd(8)=-1 ; dd(9)=-1 !to tell yaz_getdec which alternative is being evaluated (where j is the choice and i is the corresponding q)
                !but when calling from singles mar market there is no j or i so set equal to -1 
                call yaz_getdec(dd,vec(1:5),surplus,pc(1:2),asum,haveenoughtomakeindiff) !called from sol/getdec/marmkt 
            end if 
            if (de(1)>0) then 
                jmax=-1 ; qmax=-1 ; relmax=1
                if (yaz) write(200,*) "in sol/getdec marmkt de(1)>0: jmax,qmax,relmax,surplus,vec(1:5) NOW GO TO 2013"
                if (yaz) write(200,'(I4,I8,I4,6F9.2)') jmax,qmax,relmax,surplus,vec(1:5) 
                GO TO 2013 !have to calculate vmax and transfers before going to 2017 
            else 
                jmax=-1 ; qmax=-1  ; relmax=0 ; vmax(1:2)=vec(1:2)
                if (yaz) write(200,*) "in sol/getdec marmkt de(1)<=0: jmax,qmax,relmax,surplus,vec(1:5) NOW GO TO 2017"
                if (yaz) write(200,'(I4,I8,I4,6F9.2)') jmax,qmax,relmax,surplus,vec(1:5)                 
                GO TO 2017 ! go directly to 2017
            end if 
        end if 

        2013 if ( vec(5) + eps >= abs(vdif(1)-vdif(2)) )  then 
            transfers(1) = alf * surplus -  vdif(1)                                                           
            transfers(2) = (1.0_dp-alf) * (vec(5) + vdif(1) ) - alf*vdif(2)  
            !vmax(1:2)=vec(3:4)+transfers(1:2)
            vmax(1:2)=vec(1:2)+alf*surplus 
            GO TO 2017
        else if ( vec(5) <= vdif(1)-vdif(2) + eps  )  then 
            transfers(1)=0.0_dp                         !ahumarch1522 adding cornersol
            transfers(2)=vec(5)                         !ahumarch1522 adding cornersol
            vmax(1:2)=vec(3:4)+transfers(1:2)
            GO TO 2017
        else if ( vec(5) <= vdif(2)-vdif(1) + eps  )   then 
            transfers(1)=vec(5)                         !ahumarch1522 adding cornersol
            transfers(2)=0.0_dp                         !ahumarch1522 adding cornersol
            vmax(1:2)=vec(3:4)+transfers(1:2)
            GO TO 2017
        !else if (relmax==0) then 
        !    vmax(1:2)=vec(1:2)
        !    GO TO 2017
        else    
            print*, "why here"
            stop
        end if 
        2017 dd(8)=jmax
        dd(9)=qmax
        dd(10)=relmax 
        if (yaz.and.callfrom==40) write(200,*) "in sol/getdec couples end: jmax,qmax,relmax,surplus),vec(1:5),transfers"
        if (yaz.and.callfrom==40) write(200,'(I4,I8,I4,8F9.2)') jmax,qmax,relmax,surplus,vec(1:5),transfers 
        if (yaz.and.callfrom==50) write(200,*) "in sol/getdec mar mkt end: jmax,qmax,relmax,surplus),vec(1:5),transfers"
        if (yaz.and.callfrom==50) write(200,'(I4,I8,I4,8F9.2)') jmax,qmax,relmax,surplus,vec(1:5),transfers 
        if (yaz.and.callfrom==80) write(400,*) "in sim/getdec couples end: jmax,qmax,relmax,surplus),vec(1:5),transfers"
        if (yaz.and.callfrom==80) write(400,'(I4,I8,I4,8F9.2)') jmax,qmax,relmax,surplus,vec(1:5),transfers  
    end subroutine getdec

        !!!!!if (relmax==1) then
            !ahu032122 ahumarch2122 commenting out to save time haveenoughtomakeindiff_alt=(  vec(5) - asum  >= 0.0_dp  )
            !ii is the max of surplusj. 
                !check if that ii option is welldef and interior/corner
                !welldef ----> you have found your best ii and set jmax equal to that
                !not welldef ----> set surplusj(ii)=pen and update nn and go on to the next best ii 
                !    
            !In orer for NB to be well defined, we need: 1) surplus>=0 2) there exists some transfer such that transfer1>=0 and transfer2>=0
            !Condition 2 follows from the requirement that the utility transfer has to only come from current wsum and not the V's. In other words, 
            !w has to be such that it can cover any utility transfers that are needed o have the PC's hold. 
            !See notes in black notebook about why Condition 2 translates into that asum condition. 
            !Once we establish that NB is well defined, then we can go ahead and see if the optimum is interior or corner. 
            !intsol(1)=( vec(5) + eps >= abs(vdif(1)-vdif(2)) )  !interior
            !intsol(2)=( vec(5) <= vdif(1)-vdif(2) + eps  )      !corner where c1=0 and c2=wsum !ahumarch2122 ahu032122 moving eps to the other side
            !intsol(3)=( vec(5) <= vdif(2)-vdif(1) + eps  )      !corner where c1=wsum and c2=0 !ahumarch2122 ahu032122 moving eps to the other side
        !!!!!    if ( intsol(1) ) then     !IF welldef, then check whether the solution is interior or corner  
        !!!!!        transfers(1) = alf * surplus -  vdif(1)                     !ahumarch2122 ahu032122 replacing with vdif to save time                                           
        !!!!!        transfers(2) = (1.0_dp-alf) * (vec(5) + vdif(1) ) - alf*vdif(2) !ahumarch2122 ahu032122 replacing with vdif to save time   
                !vmarioj(1:2,j)=vec(1:2)+alf*surplusj(j)
                !vmario(1:2)=vec(3:4)+transfers(1:2)
                !if (vec(5)-abs(vdif(1)-vdif(2)) < 0.0_dp ) then !ahumarch1522 if the transfers statement is correct, then this should be correct too!
                !    print*, "There is something wrong!", vec(5),abs(vdif(1)-vdif(2))  !ahumarch1522
                !    stop                                                    !ahumarch1522
                !end if                                                      !ahumarch1522 
                !if ( minval(transfers) + eps2 < 0.0_dp .or. maxval(transfers) > vec(5)+eps2 ) then ! if (minval(transfers)>0.0_dp.and.maxval(transfers)<vec(5) ) then 
                !    print*, "Int cond fine but transfers not fine!", transfers(1:2),vec(5)  !ahumarch1522
                !    stop                                                                    !ahumarch1522                            
                !end if 
                !if ( minval(transfers) < 0.0_dp .or. maxval(transfers) > vec(5) ) then ! if (minval(transfers)>0.0_dp.and.maxval(transfers)<vec(5) ) then 
                !    print*, "Int cond fine but transfers not fine 2!", transfers(1:2),vec(5)  !ahumarch1522
                !    stop                                                                    !ahumarch1522                            
                !end if 
        !!!!!    else if ( intsol(2) ) then                      !ahumarch1522 adding cornersol
        !!!!!        transfers(1)=0.0_dp                         !ahumarch1522 adding cornersol
        !!!!!        transfers(2)=vec(5)                         !ahumarch1522 adding cornersol
                !vmario(1:2)=vec(3:4)+transfers(1:2)
                !vec(1:2)=vec(3:4)+transfers(1:2)
                !if (vec(5) >= vdif(1)-vdif(2) ) then        !ahumarch1522 adding cornersol
                !    print*, "There is somethingwrong in corner1!", vec(5),vdif(1)-vdif(2)  !ahumarch1522 adding cornersol
                !    stop
                !end if                                       !ahumarch1522 adding cornersol
        !!!!!    else if ( intsol(3) ) then                      !ahumarch1522 adding cornersol
        !!!!!        transfers(1)=vec(5)                         !ahumarch1522 adding cornersol
        !!!!!        transfers(2)=0.0_dp                         !ahumarch1522 adding cornersol
                !vmario(1:2)=vec(3:4)+transfers(1:2)
                !if (vec(5) >= vdif(2)-vdif(1) ) then        !ahumarch1522 adding cornersol
                !    print*, "There is somethingwrong in corner2!", vec(5),vdif(1)-vdif(2)  !ahumarch1522 adding cornersol
                !    stop
                !end if                                       !ahumarch1522 adding cornersol
        !!!!!    else 
        !!!!!        print*, "I should not end up here because we already established that we haveenoughtomakeindiff"
        !!!!!        stop
        !!!!!    end if
        !!!!!    vmax(1:2)=vec(3:4)+transfers(1:2)
        !!!!!    if (yazmax) then !ahu030822
        !!!!!        write(201,'(2x,tr1,"age",tr2,"q0",tr3,"q",tr2,"qt",tr2,"ie",TR3,"X",tr1,"rel")' ) !ahu030822
                !write(201,'(2x,7i4)')	ia,q0,q,qmax,iepsmove,X,relmax    ! altq is just the q that altrnative j corresponds to       !ahu030822
        !!!!!    end if !ahu030822
        !!!!!else  if (relmax==0) then !Second condition for welldef NB not well defined either because don't have enough to make indiff
        !!!!!    vmax(1:2)=vec(1:2)   
        !!!!!else 
        !!!!!    print*, "I should not end up here because we already established that we haveenoughtomakeindiff"
        !!!!!    stop   
        !!!!!end if 
    


        subroutine checkdecmar(vec,vdif,transfers,intsol)
            real(dp), intent(in) :: vec(5),vdif(2),transfers(2)
            logical :: intsol(3)
            !ahumarch1522 Now check if the interior condition is really indeed satisfied but get rid of this later
              if (vec(5)-abs(vdif(1)-vdif(2)) < 0.0_dp ) then !ahumarch1522 if the transfers statement is correct, then this should be correct too!
              !   print*, "There is something wrong!", ia,q0,q,x,vec(5),abs(vdif(1)-vdif(2))  !ahumarch1522
                  print*, "intsol1",intsol(1)
                  print*, "vec(5)",vec(5)
                  print*, "Vdif(1)",Vdif(1)
                  print*, "Vdif(2)",Vdif(2)
                  print*, "Vdif1-Vdif2",Vdif(1)-Vdif(2)
                  print*, "vec5 >= abs(Vdif1-Vdif2)?",vec(5),abs(Vdif(1)-Vdif(2))
                  print*, "minval(transfers) >= 0?",minval(transfers),0.0_dp                         
                  print*, "maxval(transfers) <= vec5?",maxval(transfers),vec(5)                     
                  stop                                                    !ahumarch1522
              end if                                                      !ahumarch1522 
              if ( minval(transfers) + eps2 < 0.0_dp .or. maxval(transfers) > vec(5)+eps2 ) then ! if (minval(transfers)>0.0_dp.and.maxval(transfers)<vec(5) ) then 
            !      print*, "Int cond fine but transfers not fine!", ia,q0,q,x,z,transfers(1:2),vec(5)  !ahumarch1522
                  print*, "intsol1",intsol(1)
                  print*, "vec(5)",vec(5)
                  print*, "Vdif(1)",Vdif(1)
                  print*, "Vdif(2)",Vdif(2)
                  print*, "Vdif1-Vdif2",Vdif(1)-Vdif(2)
                  print*, "vec5 >= abs(Vdif1-Vdif2)?",vec(5),abs(Vdif(1)-Vdif(2))
                   print*, "minval(transfers) >= 0?",minval(transfers),0.0_dp                         
                   print*, "maxval(transfers) <= vec5?",maxval(transfers),vec(5)                     
                   stop                                                                    !ahumarch1522                            
               end if 
               if ( minval(transfers) < 0.0_dp .or. maxval(transfers) > vec(5) ) then ! if (minval(transfers)>0.0_dp.and.maxval(transfers)<vec(5) ) then 
            !       print*, "Int cond fine but transfers not fine 2!", ia,q0,q,x,z,transfers(1:2),vec(5)  !ahumarch1522
                   print*, "intsol1",intsol(1)
                   print*, "vec(5)",vec(5)
                   print*, "Vdif(1)",Vdif(1)
                   print*, "Vdif(2)",Vdif(2)
                   print*, "Vdif1-Vdif2",Vdif(1)-Vdif(2)
                   print*, "vec5 >= abs(Vdif1-Vdif2)?",vec(5),abs(Vdif(1)-Vdif(2))
                   print*, "minval(transfers) >= 0?",minval(transfers),0.0_dp                         
                   print*, "maxval(transfers) <= vec5?",maxval(transfers),vec(5)   
                   !!transfers(2) = (1.0_dp-alf) * (vec(5) + vec(3) - vec(1) ) - alf*( vec(4)-vec(2) )   
                   !print*, "1-alf",1.0_dp-alf
                   !print*, "(vec5+vdif(1))",  (vec(5) + vec(3) - vec(1) )
                   !print*, "(1-alf)(vec5+vdif(1))",  (1.0_dp-alf) * (vec(5) + vec(3) - vec(1) )
                   !print*, "alf(vdif(2))", alf*( vec(4)-vec(2) )   
                   !print*, "transfers2", (1.0_dp-alf) * (vec(5) + vec(3) - vec(1) ) - alf*( vec(4)-vec(2) )   
                   !print*, "transfers2 alt",  (vec(5) + vec(3) - vec(1) ) - alf*( vec(5) + vec(3) - vec(1) + vec(4)-vec(2)  )   
                   !print*, "test1",  ( vec(5) + vec(3) - vec(1) + vec(4)-vec(2)  )   
                   !print*, "test2",  alf*( vec(5) + vec(3) - vec(1) + vec(4)-vec(2)  )   
                    !write(*,'(4F30.20)'),  (vec(5) + vec(3) - vec(1) ),vec(5),vec(3),vec(1)
                    ! testing(1)=vec(5)+vdif(1)
                    !write(*,'(F30.20)'), testing(1)
                    !testing(2)=alf*surplus
                    !write(*,'(F30.20)'), testing(2)
                    ! testing(3)=testing(1)-testing(2)
                    !write(*,'(F30.20)'), testing(3)                             
                !intsol1 T
                ! vec(5)   5327.9026198192287
                ! Vdif(1)   8749.1068802403315
                ! Vdif(2)   3421.2042604211028
                ! Vdif1-Vdif2   5327.9026198192287
                ! vec5 >= abs(Vdif1-Vdif2)?   5327.9026198192287        5327.9026198192287
                ! minval(transfers) >= 0?   0.0000000000000000        0.0000000000000000
                ! maxval(transfers) <= vec5?   5327.9026198192296        5327.9026198192287
                ! 1-alf  0.50000000000000000
                ! (vec5+vdif(1))   14077.009500059561
                ! (1-alf)(vec5+vdif(1))   7038.5047500297806
                ! alf(vdif(2))   1710.6021302105514
                ! transfers2   5327.9026198192296
                ! transfers2 alt   5327.9026198192296
                ! test1   17498.213760480663
                ! test2   8749.1068802403315
                !"out.txt" 133L, 10998B                        
                   stop                                                                    !ahumarch1522                            
            end if 

            !**********************************************************************************
            !THE BELOW ARE IF STATEMENTS FOR CHECKING AT THE BEGINNING OF SOL ROUTINE (AT THE THE BEGINNING OF COUPLES DO LOOP)
            !if (skriv.and.ia==48.and.q0==4.and.q==92.and.x==19.and.trueindex==1) then ; yaz=.true. ; else ; yaz=.false. ; end if !ahu030822
            !if (skriv.and.ia==18.and.q0==4.and.q==92.and.x==19.and.trueindex==1) then ; yaz=.true. ; else ; yaz=.false. ; end if !ahu030822              
            !if (skriv.and.(ia==48).and.(q0==4).and.q<=92.and.x==19.and.iepsmove==10.and.trueindex==1) then ; yaz=.true. ; else ; yaz=.false. ; end if  !ahu030622 !ahu030822
            !if (skriv.and.(ia==48).and.(q0==4).and.(q>=190.and.q<=201).and.x==1.and.trueindex==1) then ; yaz=.true. ; else ; yaz=.false. ; end if  !ahu030622 !ahu030822
                
            !if (skriv.and.(ia<=18.or.ia==45.or.ia==49.or.ia==50).and.(q0<=15).and.(q>=190.and.q<=201).and.x==1.and.trueindex==1) then ; yaz=.true. ; else ; yaz=.false. ; end if  !ahu030622 !ahu030822
                !if (skriv.and.(ia==18.or.ia==49.or.ia==50).and.x==1.and.trueindex==1.and.q0<=50) then ; yazmax=.true. ; else ; yazmax=.false. ; end if
            
            !if (skriv.and.ia==48.and.q0==4.and.q==92.and.iepsmove==10.and.x==19.and.trueindex==1) then ; yaz=.true. ; else ; yaz=.false. ; end if !ahu030822            
            !if (skriv.and.(ia==50).and.(q0<=15).and.(q>=190.and.q<=201).and.x==1.and.trueindex==1) then ; yaz=.true. ; else ; yaz=.false. ; end if  !ahumarch1122

            !ahumarch1122 if ( skriv.and.q0==10.and. q2w( qq2q(1,q) )==np2 .and. q2w( qq2q(2,q) )==np2 .and. q2l( qq2q(1,q) )==1 .and. q2l( qq2q(2,q) )==1 ) then  !ahu030622 
                !ahumarch1122 print*, 'Is this it', q,q2w( qq2q(1,q) ),q2w( qq2q(2,q) ), q2l( qq2q(1,q) ) !ahu030622 
            !ahumarch1122 end if  !ahu030622 
            !if (skriv.and.ia==48.and.trueindex==1.and.x==19.and.q0==4.and.(q==92) ) then ; yazmax=.false. ; else ; yazmax=.false. ; end if
            !if (skriv.and.(ia==18.or.ia==50).and.(q>33.and.q<=37).and.(q0>=32.and.q0<=35).and.x==2) then ; yaz=.true. ; else ; yaz=.false. ; end if  
            !if (skriv.and.(q0>=32.and.q0<=35).and.x==2.and.(q>=33.and.q<=37).and.(ia==mxa.or.ia==47)) then ; yaz=.false. ; else ; yaz=.false. ; end if  
            !if (skriv.and.(ia>=40).and.(q0==18).and.q==24.and.x==1.and.trueindex==1) then ; yaz=.true. ; else ; yaz=.false. ; end if  !ahu 0327 trueindex==2    
           !**********************************************************************************


        end subroutine 

	subroutine getdec_c(dd,vmax) 
	integer(i4b), intent(inout) :: dd(:)
	real(dp), intent(out) :: vmax(2) 
	real(dp) :: mcost(2),vsum(2),vec(5),surplusj(nc),transfers(2),vcheck(2),val(2),nashprod(nc),vdif(2),asum,vmario(2,nc)
	real(dp) :: vdifj(2,nc),vecj(5,nc),vmarioj(2,nc)
    integer(i4b) :: ia,index,q,x,z,q0,g,jmax,qmax,relmax,iepsmove,i,i0,n,trueindex,j,de(1),ed(2),locch,loc0,nn !ahu032022 ahumarch2022
	logical :: welldef,defj(nc),haveenoughtomakeindiff,haveenoughtomakeindiff_alt,intsol(3),pc(2),pc_alt(2)
    !dd= (/ ia,index,q,x,z,q0,g,jmax,qmax,relmax,iepsmove /)  		
    !dd = (/ia,trueindex,q,x,z,q0,-1,-1,-1,-1,iepsmove /) 	! (/ ia,index,q,x,z,q0,gender,jmax,qmax,relmax,iepsmove /)  	
	ia=dd(1) 
	trueindex=dd(2) 
    if (groups) then 
        index=1
    else 
        index=trueindex
    end if
    q=dd(3) 
	x=dd(4) ; ed(:)=xx2e(:,x)     
	z=dd(5) 
	q0=dd(6)  
	iepsmove=dd(11)   
    vec=pen				! initialize
    !jmax=0				! initialize
    !qmax=0				! initialize
    !relmax=0            ! initialize
    !vmax=pen			! initialize
    !surplusmax=pen      ! initialize  
    surplusj=pen
    nashprod=pen        !ahu summer18 042318
    transfers=pen
    vmario=pen          !ahumarch1522 addin,g cornersol
    defj=.FALSE.
    i = qq2q(1,q) ; n=xx2x(1,x) ; i0 = qq2q(1,q0) 
    vec(1) = vm(iepsmove,i,n,i0,ia,index) + divpenalty
    mcost(1)=movecost(i0,n,trueindex)    
    i = qq2q(2,q) ; n=xx2x(2,x) ; i0 = qq2q(2,q0) 	
    vec(2) = vf(iepsmove,i,n,i0,ia,index) + divpenalty	
    mcost(2)=movecost(i0,n,trueindex)    
    loc0=qq2l(1,q0)     !ahumarch2022 ahu032022
	!ahu030822 adding the below in order to check the mar rates decreasing with mumar situation to compare chk2's better
	!if (yazmax) then !ahu030822
		!write(201,'("********************************")')
		!write(201,'(6x,tr5,"age",tr6,"q0",tr7,"q",tr7,"z",tr7,"j",tr4,"altq",tr1,"iepsmve",tr1,"trueind",tr1,"sex",TR3,"X",tr3,"def")' ) !ahu030822
			!write(201,'(6x,8i8,i4,I4,L6)')	ia,q0,q,z,j,i,iepsmove,trueindex,-1,X,defj(j)    ! altq is just the q that altrnative j corresponds to       !ahu030822
	!end if !ahu030822
    vecj=pen ; vdifj=pen ; surplusj=pen ; vmarioj=pen
    choice: do j=1,nc	
        i = ch(j,q,q0)	!alternative q
        if (i>0 ) then		
            locch=qq2l(1,i)     !ahumarch2122 ahu032122
            dd(8:9)=(/j,i/)     !ahumarch2122 ahu022122 moved this up
            !if (qq2l(1,i).ne.qq2l(2,i)) then ; print*, 'something wrong in dec_c' ; stop; end if 
            !if (qq2l(1,q0).ne.qq2l(2,q0)) then ; print*, 'something wrong in dec_c' ; stop; end if 
            !ahumarch2122 ahu032122 vec(3) = vm0_c(i,x,ia,index)  + mg(z,trueindex) + one( qq2l(1,i) /= qq2l(1,q0)) * mcost(1) + one( qq2l(1,i) /= qq2l(1,q0)) * moveshock_m(iepsmove)  !fnmove(kid)     
            !ahumarch2122 ahu032122 vec(4) = vf0_c(i,x,ia,index)  + mg(z,trueindex) + one( qq2l(2,i) /= qq2l(2,q0)) * mcost(2) + one( qq2l(2,i) /= qq2l(2,q0)) * moveshock_f(iepsmove)  !fnmove(kid)     
            vecj(3,j) = vm0_c(i,x,ia,index)  + mg(z,trueindex) + one( locch /= loc0 ) * (mcost(1) + moveshock_m(iepsmove)/moveshockdiv ) !ahumarch2022 ahu032022     
            vecj(4,j) = vf0_c(i,x,ia,index)  + mg(z,trueindex) + one( locch /= loc0 ) * (mcost(2) + moveshock_f(iepsmove)/moveshockdiv ) !ahumarch2022 ahu032022
            !vec(5) = wc(1,i,x,trueindex) + wc(2,i,x,trueindex) + one( qq2w(1,q0)<=np )* ubc(1,x,i,trueindex) + one( qq2w(2,q0)<=np )* ubc(2,x,i,trueindex) + nonlabinc + nonlabinc 	 !ahu summer18 050318: added the ubc
            vecj(5,j) = wc(1,i,x,trueindex) + wc(2,i,x,trueindex) + nonlabinc(ed(1)) + nonlabinc(ed(2)) 
            vdifj(1,j)=vecj(3,j)-vec(1)   !ahumarch1522 ahu031522 adding cornersol
            vdifj(2,j)=vecj(4,j)-vec(2)   !ahumarch1522 ahu031522 adding cornersol
            !ahumarch2122 ahu032122 replacing this with vdif for saving time surplusj(j) = vec(5) + vec(3) - vec(1) + vec(4) - vec(2)  
            surplusj(j) = vecj(5,j) + sum(vdifj(1:2,j))  !ahumarch2122 ahu032122 replacing with vdif for saving time    
        end if   
    end do choice    


    !ii is the max of surplusj. 
    !check if that ii option is welldef and interior/corner
    !welldef ----> you have found your best ii and set jmax equal to that
    !not welldef ----> set surplusj(ii)=pen and update nn and go on to the next best ii 
    !
    jmax=-1 ; nn=0 !initiate jmax and the below do loop goes until this jmax is changed to a positive integer value
    if (maxval(surplusj)>0.0_dp) then 
    do while (jmax<0.and.nn<=nc+1.and.maxval(surplusj)>0.0_dp)
        !check if that max is well defined or not. if it is  then move on. if not then lookat the rest of the surplusj's. 
        !if (trueindex==1)   print*, "HEre it isat first",nn,jmax,surplusj(1:nc)
        de=maxloc(surplusj)
        j=de(1)
        !if (trueindex==1.and.ia==49.and.q0==4.and.x==19.and.q>=40.and.q<=100)   print*, "HEre it isat first",nn,jmax,j,surplusj(1:nc),maxval(surplusj),(maxval(surplusj)>0.0_dp)
        if (surplusj(j)+eps > 0.0_dp ) then !First condition for welldef !ahumarch2822 ahu032822
            i = ch(j,q,q0)	!alternative q
            !locch=qq2l(1,i)     
            !dd(8:9)=(/j,i/)     
            pc(1:2)	= ( vdifj(:,j) + eps >= 0.0_dp )
            asum = sum(  one(.not. pc)  *   abs( vdifj(:,j) )   ) 
            haveenoughtomakeindiff=(  vecj(5,j) + eps - asum  >= 0.0_dp  )
            intsol(1)=( vecj(5,j) + eps >= abs(vdifj(1,j)-vdifj(2,j)) )  
            intsol(2)=( vecj(5,j) <= vdifj(1,j)-vdifj(2,j) + eps  )      !corner where c1=0 and c2=wsum !ahumarch2122 ahu032122 moving eps to the other side
            intsol(3)=( vecj(5,j) <= vdifj(2,j)-vdifj(1,j) + eps  )      !corner where c1=wsum and c2=0 !ahumarch2122 ahu032122 moving eps to the other side
            if (haveenoughtomakeindiff ) then !Second condition for welldef (haveenoughtomakeindiff) ahumarch2822 ahu032822
                defj(j)=.TRUE.    ; jmax=j      

                if ( intsol(1) ) then     !IF welldef, then check whether the solution is interior or corner  
                    transfers(1) = alf * surplusj(j) -  vdifj(1,j)                     !ahumarch2122 ahu032122 replacing with vdif to save time                                           
                    transfers(2) = (1.0_dp-alf) * (vecj(5,j) + vdifj(1,j) ) - alf*vdifj(2,j) !ahumarch2122 ahu032122 replacing with vdif to save time   
                    vmarioj(1:2,j)=vecj(3:4,j)+transfers(1:2)
                    !vmarioj(1:2,j)=vec(1:2)+alf*surplusj(j)
                    !if (vec(5)-abs(vdif(1)-vdif(2)) < 0.0_dp ) then !ahumarch1522 if the transfers statement is correct, then this should be correct too!
                    !    print*, "There is something wrong!", vec(5),abs(vdif(1)-vdif(2))  !ahumarch1522
                    !    stop                                                    !ahumarch1522
                    !end if                                                      !ahumarch1522 
                    !if ( minval(transfers) + eps2 < 0.0_dp .or. maxval(transfers) > vec(5)+eps2 ) then ! if (minval(transfers)>0.0_dp.and.maxval(transfers)<vec(5) ) then 
                    !    print*, "Int cond fine but transfers not fine!", transfers(1:2),vec(5)  !ahumarch1522
                    !    stop                                                                    !ahumarch1522                            
                    !end if 
                    !if ( minval(transfers) < 0.0_dp .or. maxval(transfers) > vec(5) ) then ! if (minval(transfers)>0.0_dp.and.maxval(transfers)<vec(5) ) then 
                    !    print*, "Int cond fine but transfers not fine 2!", transfers(1:2),vec(5)  !ahumarch1522
                    !    stop                                                                    !ahumarch1522                            
                    !end if 
                else if ( intsol(2) ) then                      !ahumarch1522 adding cornersol
                    transfers(1)=0.0_dp                         !ahumarch1522 adding cornersol
                    transfers(2)=vecj(5,j)                         !ahumarch1522 adding cornersol
                    vmarioj(1:2,j)=vecj(3:4,j)+transfers(1:2)
                    !if (vec(5) >= vdif(1)-vdif(2) ) then        !ahumarch1522 adding cornersol
                    !    print*, "There is somethingwrong in corner1!", vec(5),vdif(1)-vdif(2)  !ahumarch1522 adding cornersol
                    !    stop
                    !end if                                       !ahumarch1522 adding cornersol
                else if ( intsol(3) ) then                      !ahumarch1522 adding cornersol
                    transfers(1)=vecj(5,j)                         !ahumarch1522 adding cornersol
                    transfers(2)=0.0_dp                         !ahumarch1522 adding cornersol
                    vmarioj(1:2,j)=vecj(3:4,j)+transfers(1:2)
                    !if (vec(5) >= vdif(2)-vdif(1) ) then        !ahumarch1522 adding cornersol
                    !    print*, "There is somethingwrong in corner2!", vec(5),vdif(1)-vdif(2)  !ahumarch1522 adding cornersol
                    !    stop
                    !end if                                       !ahumarch1522 adding cornersol
                else 
                    print*, "I should not end up here because we already established that we haveenoughtomakeindiff"
                    stop
                end if
            else    !Second condition for welldef NB not well defined either because don't have enough to make indiff
                defj(j)=.FALSE. ; jmax=-1 ;    surplusj(j)=pen  
                !surplusj(j)=pen
                !nashprod(j)=pen  !ahumarch1122 uncomenting this out
                !transfers=pen
            end if !Second condition for welldef (haveenoughtomakeindiff) ahumarch2822 ahu032822
        else !First condition for welldef does not hold (surplus is not positive) ahumarch2822 ahu 032822
            defj(j)=.FALSE. ; jmax=-1   ;    surplusj(j)=pen   
        end if  !First condition for well def ( is surplus positive? surplus>0?) ahumarch2822 ahu032822 
        nn=nn+1
        !if (trueindex==1.and.ia==49.and.q0==4.and.x==19.and.q>=40.and.q<=100)  print*, "HEre it is",nn,jmax,j,defj(j),vmarioj(1:2,j)
        !if done looking for maximum then you will be done 
    end do !if done looking for max, then exit ! If not, then next check to see if the next best if weldef and interior or corner
    end if !maxval(surplusj)>0
    !if (trueindex==1) print*, "I barely got out but then not really!",nn,jmax
    if (jmax>0) then 
        qmax=ch(jmax,q,q0)
        relmax=1
        !vmax(1:2)=vec(1:2)+0.5_dp*surplusj(jmax)
        vmax(1:2)=vmarioj(1:2,jmax)
    else
        jmax=0 ; qmax=0
        relmax=0
        vmax(1:2)=vec(1:2)
    end if
    dd(8)=jmax
    dd(9)=qmax
    dd(10)=relmax 
	!if (yazmax) then !ahu030822
	!	write(201,'(2x,tr1,"age",tr2,"q0",tr3,"q",tr2,"qt",tr2,"ie",TR3,"X",tr1,"rel")' ) !ahu030822
	!	write(201,'(2x,7i4)')	ia,q0,q,qmax,iepsmove,X,relmax    ! altq is just the q that altrnative j corresponds to       !ahu030822
	!end if !ahu030822
    end subroutine getdec_c


    ! Look at dropbox/familymig/v130814/Source2 and Source11 and temp_010413 for the old versions of this
	subroutine checknb( dd, vec , def ,vsum)   
	integer, intent(in) :: dd(:) 
	real(dp), intent(in) :: vec(5)        
	logical, intent(out) :: def
	real(dp), intent(out) :: vsum(2)
	real(dp) :: surplus,transfers(2)
	vsum = pen ; transfers=pen 
	def  = .FALSE.
    surplus = vec(5) + vec(3) - vec(1) + vec(4) - vec(2)  
	if ( surplus + eps >= 0.0_dp )  then
		transfers(1) = alf * surplus - ( vec(3)-vec(1) )                                                    ! wage
		transfers(2) = (1.-alf) * (vec(5) + vec(3) - vec(1) ) - alf*( vec(4)-vec(2) )                       ! firm payoff
        !if (transfers(1) < 0.0_dp) cornersol=(/ 0._dp,vec(5) /) 
		!if (transfers(2) < 0.0_dp) cornersol=(/ vec(5),0._dp /) 
        !if not corner if statement. but not sure if this is right. look at the other tries below. 
        !because the below checks for corner solution o rinterior solution (checking transfers>0 or <w afterwards is the same as 
        !by the way checking the condition for checking FOC at corners 0 and w. see notes p.3 . 
        !but before this you need to figure otu whether there are 
        !any feasible utility allocations through doing that checking: 
        !asum = sum(  abs(.not. pc)  *   abs( vdif )   )
	    !def  = ( w + eps - asum  >= 0.0_dp )
    
		if ( minval(transfers) + eps >= 0.0_dp .and. maxval(transfers) < vec(5) ) then ! if (minval(transfers)>0.0_dp.and.maxval(transfers)<vec(5) ) then 
            def=.TRUE.
		    vsum(1) = transfers(1) + vec(3)      ! vsum is the value function after all transfers are made. i.e. vsum = wage + beta*EV_worker         ! vec(3) is just beta*EV i.e. what the worker has before transfers are made
		    vsum(2) = transfers(2) + vec(4)      ! vsum = firmpayoff + beta*EV_firm     ! " 
            if ( vsum(1) + eps < vec(1) .or. vsum(2) + eps < vec(2) ) then
                print*, "Transfers positive and less than wages but someone is getting a bad deal "
		        write(*,'(2x,2(tr6,"vbar"),2(tr8,"vc"),tr8,"wc",tr6,"surp")' )    
		        write(*,'(2x,6f10.2)') vec,surplus
		        write(*,'(2x,2(tr5,"trans"),2(tr6,"vsum") )' )    
		        write(*,'(2x,4f10.2)') transfers,vsum
                write(*,*) vsum(1)-vec(1),vsum(2)-vec(2)
            end if 
        end if
	end if
	end subroutine checknb

	!if ( pc(1) .and. pc(2) ) then	
	!	def  = .true.				
	!else 
	!	asum = sum(  abs(.not. pc)  *   abs( vdif )   )
	!	def  = ( w + eps - asum  >= 0.0_dp )
	!end if 
	!if (def) then 
	!	vsum	= vec(1:2) + 0.5_dp * ( w + eps + sum(vdif) )
	!	if ( vsum(1) < vec(1) .or. vsum(2) < vec(2) ) then 
	!		print*, "error in checknb " ,vsum(1)-vec(1),vsum(2)-vec(2),( w - asum  >= 0.0_dp ),w,asum,pc,w + sum(vdif) 
	!		print*, "vdif(2),wf_s,wf_s+vdif(2),wf_s+vdif(2) ", vdif(2)  !,wagetemp(2),wagetemp(2) + vdif(2),wagetemp(2) + vec(4) - vec(2) 
	!		stop 
	!	end if
	!end if 
	!end subroutine checknb

	!vsum = penaltyval 
	!transfers = penaltyval
	!def  = .FALSE.
	!criter=.FALSE.
	!surplus = vec(5) + vec(3) - vec(1) + vec(4) - vec(2)  
	! CASE 1: NO CONSTRAINT ON TRANSFERS. I.E. WAGES AND FIRM SHARES CAN BE NEGATIVE. 
	!         When there are no limits on what can be transferred, the only criterion for NB to be well-defined is simply that the match surplus is positive.  
	!if (nonneg_constraint==0) then                      ! no nonneg constraints for wages nor firm payoff  
	!	def = ( surplus + eps >= 0.0 )                  ! When there are no limits on what can be transferred btw parties, the only criterion for NB to be well defined is that surplus is positive. 
		
	! CASE 2: NONNEGATIVITY CONSTRAINTS ON TRANSFERS: 
	!           So that transfers can only come from current period resources. Cannot borrow from the continuation values to pay the other party today. 
	!           This limits the set of of feasible allocations and therefore the number of times NB is welldefined is less in this case. Because now any divsion of total surplus can only be implemented by using current period resources, which might not be enough sometimes. 
	!else if (nonneg_constraint==2) then                 ! see notes for how the nonneg constraints mean that the match surplus (outputnet) needs to satisfy the following three criteria
	!	criter(1) = ( vec(5) + eps >= vec(1) - vec(3) ) 
	!	criter(2) = ( vec(5) + eps >= vec(2) - vec(4) ) 
	!	criter(3) = ( vec(5) + eps >= vec(1) - vec(3) + vec(2) - vec(4) ) 
	!	def = ( criter(1) .and. criter(2) .and. criter(3) ) 
	!end if 
	

	! Check interior opt conditions 
	! If NB is well defined, then take FOC to get optimal transfers (wage and firm payoffs)
	!if (def) then 
	!	transfers(1) = alpha * surplus - ( vec(3)-vec(1) )                                                      ! wage
	!	transfers(2) = (1.-alpha) * (vec(5) + vec(3) - vec(1) ) - alpha*( vec(4)-vec(2) )                       ! firm payoff
	!	if (transfers(1) < 0.) transfers=(/ 0.,vec(5) /) 
	!	if (transfers(2) < 0.) transfers=(/ vec(5),0. /) 
	!	vsum(1) = transfers(1) + vec(3)      ! vsum is the value function after all transfers are made. i.e. vsum = wage + beta*EV_worker         ! vec(3) is just beta*EV i.e. what the worker has before transfers are made
	!	vsum(2) = transfers(2) + vec(4)      ! vsum = firmpayoff + beta*EV_firm     ! " 
	!	if (writeval) then ! if want to check whether things look right: check whether vsum's and the below temp's are equal. they should be since they are the same way of calculating value functions after the transfers.        
	!   end if 
	!end if 
	

	subroutine getdec_s_old(dd,vs,qs,vs1)
	integer(i4b), intent(in) :: dd(:)
	real(dp), dimension(nqs,nxs), intent(in) :: vs
	integer(i4b), dimension(nepsmove,nqs,nxs,nqs), intent(out) :: qs
	real(dp), dimension(nepsmove,nqs,nxs,nqs), intent(out) :: vs1
	real(dp), dimension(np2,nl) :: dumv0,dumv1,dumv2,dumv
    real(dp) :: v0,v_cond_u,vbar,totmovecost(nl),moveshock(nepsmove)
	integer(i4b) :: q0,q,qbar,qstar,q_cond_u,l,l_cond_u,w,x,ss(size(dd,1)),iepsmove,iepskid,sex,index,trueindex

    trueindex=dd(2)
    sex=dd(7)
	qs = -1 ; vs1 = 0.0_dp
    if (groups) then 
        index=1
    else 
        index=trueindex
    end if
    if (sex==1) then 
        moveshock=moveshock_m
    else if (sex==2) then 
        moveshock=moveshock_f
    else 
        print*, 'decs something wrong'
        stop
    end if 
	do q0=1,nqs
		if ( q2w(q0)<=np1 ) then 
			do x=1,nxs
                dumv0(:,:) = reshape( vs(:,x) , (/ np2, nl /) )	                ! turn q index into its components (w,l)
                !dumv1(:,:) = reshape( ubs(sex,:,x,trueindex) , (/ np2, nl /) )	! turn q index into its components (w,l)    !ahu summer18 050318
                !dumv2 = dumv0 + dumv1 * one(q2w(q0) <= np)                      !!ahu summer18 050318: !if they were working the previous period, then they get to collect unemp benefit but not otherwise
                dumv2=dumv0
                
                moveshocks: do iepsmove=1,nepsmove
				    !Before drawing the rest of the shocks (ofloc and wagedraw), determine some max options conditional on unemployment and given status quo (no need to know the ofloc and wagedraw for this)
                    totmovecost(:)=movecost(q0,x,trueindex) ; totmovecost(q2l(q0))=moveshock(iepsmove) !0.0_dp
                    do w=1,np1
					    dumv(w,:) = dumv2(w,:) + totmovecost(:) !ahu summer18 050318 turned dumv0 into dumv2
				    end do 
                    
				    l_cond_u=maxloc(dumv(np1,:),1)	! max location option conditional on unemployment 
				    q_cond_u=wl2q(np1,l_cond_u )    ! turn that into q
				    v_cond_u=dumv(np1,l_cond_u)		! max val conditional on unemployment
				    v0=dumv(q2w(q0),q2l(q0))		! value at status quo
				    if (v0>v_cond_u) then			! choice btw status quo (vs(q0,x)) and max unemployment option (vs(i,x) calculted above:
					    vbar=v0 ; qbar=q0
				    else 
					    vbar=v_cond_u	  ; qbar=q_cond_u
				    end if 					
				    if ( q2w(qbar)==np2) then ; print*, "in getdec_s and q2w(j) is np2 " ; stop ; end if 
				    if (qbar==0) then ; print*, "in getdec_s: j is 0 ", dd(1),dd(2) ; stop ; end if 
				    ofloc: do l=1,nl											
					    !r = locate( real(dumv(1:np,l)) , real(vs(j,x))   )	! res wage for each location l, above which you accept the offer from location l
					    wagedraw: do w=1,np2
						    q = wl2q(w,l)
						    if (w <= np) then 
							    qstar = qbar*one( dumv(w,l)<=vbar) + q*one(dumv(w,l)>vbar)
						    else if (w == np1) then 
							    qstar  = q_cond_u	! if you get laid off, the only things to choose from are unemployment option and the max thing to do is qu in that case which is calculated above 
						    else if (w == np2) then 
							    qstar  = qbar		! if nothing happens then choice is q0u which is the maximum btw status quo and the best unemployment option 
						    end if 		
						    if (q2w(qstar)==np2) then ; print*, "in getdec_s and q2w(qstar) is np2 " ; stop ; end if 
						    qs(iepsmove,q,x,q0)  =  qstar				
						    vs1(iepsmove,q,x,q0) =	dumv(q2w(qstar),q2l(qstar))
                            !if (     abs( vs(qstar,x) - dumv(q2w(qstar),q2l(qstar))+ totmovecost(q2l(qstar))    )   > eps ) then ; print*, "in getdec_s and vs/=dumv ",vs(qstar,x),dumv(q2w(qstar),q2l(qstar))-totmovecost(q2l(qstar)) ; stop ; end if  		
						    if (yaz.and.q2w(q0)==1.and.l==8.and.x==1) then
							    ss = dd  !(/ ia,index,q,x,z,q0,gender,j,altq,iepsmove /)  													
							    ss(3:6)= (/q,x,-1,q0/)
                                ss(11)=iepsmove
							    !call yaz_decs(ss,qstar,v0,v_cond_u,vbar,dumv)   !look here
						    end if 
                        end do wagedraw
				    end do ofloc
                end do moveshocks
			end do  !x 				
		end if		! w0<=np1 i.e. ps0(q0) is true
	end do			! q0
	end subroutine getdec_s_old


	subroutine getdec_s(dd,vs,qs,vs1)
	integer(i4b), intent(in) :: dd(:)
	real(dp), dimension(nxs,nqs), intent(in) :: vs
	integer(i4b), dimension(nepsmove,nxs,nqs,nqs), intent(out) :: qs
	real(dp), dimension(nepsmove,nxs,nqs,nqs), intent(out) :: vs1
	integer(i4b) :: q0,q,l,l0,w,x,i,j,iepsmove,sex,index,trueindex,jstar(1),age,ss(size(dd,1)),ia
    real(dp) :: mcost,vcho(ncs),moveshock(nepsmove),bshock(nepsmove)

    trueindex=dd(2)
    sex=dd(7) ; ia=dd(1)
	qs = -1 ; vs1 = 0.0_dp
    if (groups) then 
        index=1
    else 
        index=trueindex
    end if
    if (sex==1) then 
        moveshock=moveshock_m
        bshock=bshock_m
    else if (sex==2) then 
        moveshock=moveshock_f
        bshock=bshock_f
    else 
        print*, 'decs something wrong'
        stop
    end if 
	do q0=1,nqs
        l0=q2l(q0)
		if ( q2w(q0)<=np1 ) then 
            ofloc: do l=1,nl											
            wagedraw: do w=1,np2
            q = wl2q(w,l)
            do x=1,nxs
                mcost=movecost(x,q0,trueindex)        
                moveshocks: do iepsmove=1,nepsmove

                vcho=pen
                choice: do j=1,ncs	
                    i = chs(j,q,q0)	!alternative q
                    if (i>0 ) then		
                        if (q2w(i)<=np) then
                            vcho(j) = vs(x,i) + one(q2l(i)/=l0) * mcost + one( q2l(i)/=l0) * moveshock(iepsmove)  !fnmove(kid)  
                        else if (q2w(i)==np1) then 
                            vcho(j) = vs(x,i) + one(q2l(i)/=l0) * mcost + one( q2l(i)/=l0) * moveshock(iepsmove) + bshock(iepsmove)   
                        end if                                         
                    end if 
                end do choice
                
                            				
                vs1(iepsmove,x,q,q0) =	maxval(vcho)
                jstar(1)=maxloc(vcho,1)
                qs(iepsmove,x,q,q0)  =  chs(jstar(1),q,q0)
                if ( vs1(iepsmove,x,q,q0) < pen+1.0_dp) then 
                    print*, 'There is something wrong in decs new new new',vs1(iepsmove,x,q,q0),dd(1),moveshock(iepsmove),mcost
                    print*, 'vcho',vcho
                    stop
                end if 
                
                if (yaz .and.q2w(q0)==1.and.l==8.and.x==1) then
                    ss = dd  !(/ ia,index,q,x,z,q0,gender,j,altq,iepsmove /)  													
                    ss(3:6)= (/q,x,-1,q0/)
                    ss(11)=iepsmove
                    !!!!!!call yaz_decs(ss,vcho)   !look here
                end if 

            end do moveshocks		
            end do !x
            end do wagedraw
            end do ofloc             
		    end if		! w0<=np1 i.e. ps0(q0) is true
	end do			! q0
	end subroutine getdec_s

    
	subroutine get_util_w
	integer(i4b) :: q,x,w(2),l(2),kid(2),ed(2),expe(2),trueindex,trueco,truetyp,truehome,g,j,k
    real(dp) :: epsw(2) !RHO(2,2),CD(2,2),
    real(sp) :: pwages(1:numbin),swages(1:numbin),staterate,fedrate,wsgross,wcgross(2)
    integer(i4b) :: bracket,bracketprev,bracketnext,sbrack,pbrack

    utils=pen
    utilc=pen
    ws=pen
    wc=pen   
    
    !if (mysay==1) then 
    !    open(unit=68857,file='checktax1.txt')
    !    open(unit=68858,file='checktax2.txt')
    !end if 

    do trueindex=1,ninp    
    call index2cotyphome(trueindex,trueco,truetyp,truehome)
		qs: do q=1,nqs 
            xs: do x=1,nxs    
        !call x2edexpkid(x,ed,exp,kid)    
            movecost(x,q,trueindex)=fnmove( q2w(q),x2kid(x),trueindex)
            do g=1,2
                w(g) = q2w(q)						! wage 
			    l(g) = q2l(q)						! location
                pwages(1:numbin)=tax(1:numbin,numbin,l(g))%pwages
			    if ( w(g) <= np ) then	
                    epsw(g)=wg(w(g),g) !sig_wge(g)*wg(w(g),g)
				    ws(g,x,q,trueindex)	= fnwge(g,truetyp, l(g),epsw(g), x2e(x), x2r(x)) 
                    wsnet(g,x,q,trueindex)=0.5_dp*ws(g,x,q,trueindex)
                    !if (mysay==1) then 
                        !if (bracket<numbin) then ; bracketnext=bracket+1 ; else ; bracketnext=numbin ; end if 
                        !if (bracket>1) then ; bracketprev=bracket-1 ; else ; bracketprev=1 ; end if 
                        !write(68857,'(tr3,"wsgross",tr1,"bracket",tr4,"prev",tr4,"here",tr4,"next")')
                        !write(68857,'(f10.1,i4,3f10.1)') wsgross,bracket,pwages(bracketprev),pwages(bracket),pwages(bracketnext)
                        !write(68857,*) 
                        !write(68858,'(tr2,"truind",tr7,"q",tr7,"x",tr1,"sex",tr1,"typ",tr3,"l",tr2,"ed",tr1,"exp",tr3,"wsgross",tr5,"wsnet",tr2,"statetax",tr4,"fedtax")')
                        !write(68858,'(3i8,5i4,4f10.2)') trueindex,q,x,g,truetyp,l(g),x2e(x),x2r(x),wsgross,wsnet(g,x,q,trueindex),staterate,fedrate
                        !write(68858,*) 
                    !end if 
                else if ( w(g) == np1 ) then 
                    ws(g,x,q,trueindex)	 = 0.0_dp
                    wsnet(g,x,q,trueindex)=0.0_dp
                end if 
                
                if (taxset==0) then
                    wsnet(g,x,q,trueindex)=ws(g,x,q,trueindex)
                end if 

            
                if ( w(g) <= np ) then						
				    utils(g,x,q,trueindex)	= uhomet(truetyp) * one(l(g)==truehome) + uhome(g) * one(l(g)==truehome) + uloc(l(g)) + nonlabinc(x2e(x))
			    else if ( w(g) == np1 ) then 
				    utils(g,x,q,trueindex)	= uhomet(truetyp) * one(l(g)==truehome) + uhome(g) * one(l(g)==truehome)  + uloc(l(g)) + alphaed(g,x2e(x))  + alphakid(g) * one(x2kid(x)>1)  + nonlabinc(x2e(x))
			    end if
                
                !ahu october2022: 
                !alphakid(g): kid(g)=1 is no kid, and kid(g)=2 and above is yes kid.   this alphakid used to have two dimensions for no reason. 
                !nkid is never above 2 though (i set it that way in data too for numkids above 2 is just always 2 in data). 
                !alphaed(g,ed(g)):   ed(g)=1 is noed and ed(g)=2 is yes ed (where  neduc is 2)    
            end do !gender
		end do xs
	end do qs
    !close(68857)
    !close(68858)

    qc: do q=1,nq
	xc: do x=1,nx
        ed(:)=xx2e(:,x)    
        expe(:)=xx2r(:,x)    
        kid(:)=xx2kid(:,x)           
        w(:) = qq2w(:,q)						! wage 
        l(:) = qq2l(:,q)						! location
        if (l(1).ne.l(2)) then ; print*, 'lm not equal to lf' ; stop ; end if 

        !******************************
        !ahu summer18 050318 
        !if (w(1)<=np) then 
        !    ubc(1,q,x,trueindex)	= 0.0_dp           
        !else if (w(1)==np1) then 
        !    epsw(1)=0.0_dp
        !    ubc(1,q,x,trueindex)	= replacement_rate*fnwge(1,truetyp, l(1),epsw(1), ed(1), expe(1) ) 
        !end if 
        !if (w(2)<=np) then 
        !    ubc(2,q,x,trueindex)	= 0.0_dp           
        !else if (w(2)==np1) then 
        !    epsw(2)=0.0_dp
        !    ubc(2,q,x,trueindex)	= replacement_rate*fnwge(2,truetyp, l(2),epsw(2), ed(2), expe(2) ) 
        !end if 
        !ahu summer18 050318 
        !******************************

        pwages(1:numbin)=tax(1:numbin,numbin,l(1))%pwages
        swages(1:numbin)=tax(numbin,1:numbin,l(2))%swages
        if ( w(1) <= np .and. w(2) <= np ) then		
            epsw(1)=wg(w(1),1) !CD(1,1)*wg(w(1),1)
            epsw(2)=wg(w(2),2) !CD(2,1)*wg(w(1),1) + CD(2,2)*wg(w(2),2)
            wc(1,x,q,trueindex)	= fnwge(1,truetyp, l(1),epsw(1), ed(1), expe(1) ) 
            wc(2,x,q,trueindex)	= fnwge(2,truetyp, l(2),epsw(2), ed(2), expe(2) )  
            wcnet(1:2,x,q,trueindex)	= 0.5_dp*wc(1:2,x,q,trueindex)   
        else if ( w(1) <= np .and. w(2) == np1 ) then		
            epsw(1)=wg(w(1),1) !sig_wge(1)*wg(w(1),1)
            wc(1,x,q,trueindex)	= fnwge(1,truetyp, l(1),epsw(1), ed(1), expe(1) ) 
            wc(2,x,q,trueindex)	= 0.0_dp
            wcnet(1:2,x,q,trueindex)	= 0.5_dp*wc(1:2,x,q,trueindex)   
        else if ( w(1) == np1 .and. w(2) <= np ) then		
            epsw(2)=wg(w(2),2) !sig_wge(2)*wg(w(2),2)
            wc(1,x,q,trueindex)	= 0.0_dp
            wc(2,x,q,trueindex)	= fnwge(2,truetyp, l(2),epsw(2), ed(2), expe(2) )   
            wcnet(1:2,x,q,trueindex)	= 0.5_dp*wc(1:2,x,q,trueindex)   
        else if ( w(1) == np1 .and. w(2) == np1 ) then		
            wc(1,x,q,trueindex)	= 0.0_dp
            wc(2,x,q,trueindex)	= 0.0_dp           
            wcnet(1:2,x,q,trueindex)	= 0.0_dp   
        end if 

        if (taxset==0) then
            wcnet(1:2,x,q,trueindex)	= wc(1:2,x,q,trueindex)
        end if 

            
        do g=1,2
            if ( w(g) <= np ) then						
                utilc(g,x,q,trueindex)	= uhomet(truetyp) * one(l(g)==truehome) +  uhome(g) * one(l(g)==truehome)   + uloc(l(g))
            else if ( w(g) == np1 ) then 
                utilc(g,x,q,trueindex)	= uhomet(truetyp) * one(l(g)==truehome) + uhome(g) * one(l(g)==truehome)   + uloc(l(g)) + alphaed(g,ed(g) ) + alphakid(g) * one(kid(g)>1) 
                !ahu october2022: 
                !alphakid(g): kid(g)=1 is no kid, and kid(g)=2 and above is yes kid.   this alphakid used to have two dimensions for no reason. 
                !nkid is never above 2 though (i set it that way in data too for numkids above 2 is just always 2 in data). 
                !alphaed(g,ed(g)):   ed(g)=1 is noed and ed(g)=2 is yes ed (where  neduc is 2)
            end if
        end do   
        
    end do xc
	end do qc
    end do !trueindex
    
    
    !Male and female wage shocks for each period and from each location
!    do it=1,nt
!	    do iloc=1,nloc
!		    do rn=1,nn
!			    up1(it,rn,iloc)=CD(1,1)*rv1(it,rn,iloc)
!			    up2(it,rn,iloc)=CD(2,1)*rv1(it,rn,iloc)+CD(2,2)*rv2(it,rn,iloc)
!		    end do 
!	    end do
!    end do 

	end subroutine get_util_w
    
    
end module sol






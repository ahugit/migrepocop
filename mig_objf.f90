
module objf
	use params
	use sol, only: solve
	use mom
	implicit none 
	!include 'mpif.h'
	real(dp) :: q_save(numit),qcont_save(nmom,numit)
	real(dp), dimension(npars,numit) :: par_save,realpar_save
	real(dp), dimension(nmom,numit) :: momdat_save,momsim_save,vardat_save
	integer(i4b), dimension(nmom,numit) :: cntdat_save,cntsim_save
	character(len=namelen), dimension(nmom) :: name
	character(len=120), dimension(nmom) :: header
	integer(i4b), dimension(nmom) :: headerloc
contains
	subroutine objfunc(parvec,objval)
	! puts everything together: takes the parameter vector, computes msm objective function at that value
	real(8), dimension(npars), intent(in) :: parvec	
	real(8), intent(out) :: objval					
	real(8), dimension(npars) :: realparvec
	type(initcond), dimension(ndata) :: init
	type(statevar), dimension(:,:), allocatable :: dat
	real(dp), dimension(nmom), save :: momdat,vardat    !deljan03
	integer(i4b), dimension(nmom), save :: cntdat       !deljan03
	real(dp), dimension(nmom) :: momsim,varsim,mymom,myvar,msm_wgt,momwgt,qcont
	integer(i4b), dimension(nmom) :: cntsim,mycnt
	integer(i4b), parameter :: datcountbar=10 !ahu jan19 010219 changing from 20 to 10 in order for emp at age 18 to count (thereis only 12 in that cell) because it identifies offer rates
	real(dp) :: time(10),mutemp1,mutemp2,mutemp3
    real(4) :: timing(6)
	integer(i4b) :: i,t,mycommrank,mpierr   
    !print*, "iter,mysay,iwritegen",iter,mysay,iwritegen !ag090522 agsept2022
    !initiate
	objval=-99.0_dp
    realparvec=-99.0_dp
    momsim=-99.0_dp
    varsim=-99.0_dp
    mymom=-99.0_Dp
    myvar=-99.0_dp
    msm_wgt=-99.0_dp 
    momwgt=-99.0_dp 
    qcont=-99.0_dp    
    cntsim=9999
    mycnt=9999
    
	time(1)=secnds(0.0)
    
    if (groups) then 
        nindex=nin
    else 
        nindex=ninp
    end if 
    allocate( decm0_s(nepsmove,nxs,nqs,nqs,mna:mxa,nindex),decf0_s(nepsmove,nxs,nqs,nqs,mna:mxa,nindex))
    allocate(decm_s(nepsmove,nxs,nqs,nqs,mna:mxa,nindex),decf_s(nepsmove,nxs,nqs,nqs,mna:mxa,nindex))
	allocate(vm(nepsmove,nxs,nqs,nqs,mna:mxa,nindex),vf(nepsmove,nxs,nqs,nqs,mna:mxa,nindex))
    allocate(dec_mar(nz,nx,nq,mna:mxa,nindex))
	allocate(vm0_c(nx,nq,mna:mxa,nindex),vf0_c(nx,nq,mna:mxa,nindex))
    allocate(vm0ctemp(nq,nx),vf0ctemp(nq,nx))

	call getpars(parvec,realparvec) 
    parsforcheck=parvec !jusr for checking where it says emax is negative. can get rid of later.
    if (iter>1) skriv=.false.
    !if (mysay==0.and.iter==3) then 
    !    skriv=.true.
    !    print*, 'here mysay,iter ',mysay, iter
    !else 
    !    skriv=.false.
    !end if
    !if (iter==1.and.iwritegen==1) then  !ahu 121118
    !    open(unit=98799, file='parameters.txt',status='replace')
    !end if 
	!if (iwritegen==1) then ; write(98799,'("iter ",i6)') iter ; end if  !ahu 121118
	!if (iwritegen==1) then !ahu 121118
    !    do i=1,npars
    !        write(98799,*) parvec(i)
    !    end do 
    !    write(98799,*)
    !end if
    !if (iter==1.and.iwritegen==1) then  !ahu 121118
    !    open(unit=987991, file='paros.txt',status='replace')
    !    do i=1,npars
    !        write(987991,'(1a15,f14.4)') parname(i),realparvec(i)
	!    end do 
    !    close(987991)
    !end if
    
    
	if (iter==1) then 
		best=huge(1.0_dp)
		momdat_save=0.0_dp
		vardat_save=0.0_dp
		cntdat_save=0
		momsim_save=0.0_dp
		cntsim_save=0
		qcont_save=0.0_dp
		q_save=0.0_dp
		par_save=0.0_dp
		realpar_save=0.0_dp
	end if 		
	if (iter<=numit) then
        !initiate
        momdat=-99.0_dp ; vardat=-99.0_dp ; cntdat=9999
		call getones                !ones
		call getdistpop	            !distance and popsize	
		call getq2q					! trans btw couple and single state space
		call getx2x					! trans btw couple and single state space
		call getch_single			! choice set cond on state variable and shock for singles
		call getch_couple			! choice set cond on state variable and shock for singles
		allocate(dat(mnad:mxa,ndata))
		call read_actualdata(init,dat)
		!  calculate moments from data, store in datamoments
		!  momentname,momentheaders,headerlocs,weights will be filled in when calculate moments from simulations
		!  so that they're not all floating around as globals. loaded into temporary variables here
		call get_mom(dat,ndata,momdat,cntdat,vardat,name,header,headerloc,momwgt)
        !print*, name(70),momdat(70),cntdat(70)
        do i=1,nmom
            if (calcvar(i)==1) then 
                if ( cntdat(i) > 0) then
                    mutemp1=momdat(i)       !/cntdat(i)
                    mutemp2=momdat(i+1)     !/cntdat(i+1)
                    momdat(i)= mutemp1   
                    momdat(i+1)= mutemp2  - mutemp1**2
                end if 
            else if (calcorr(i)==1) then 
                if ( cntdat(i) > 0) then
                    mutemp1=momdat(i)   !/cntdat(i)
                    mutemp2=momdat(i+1) !/cntdat(i+1)
                    mutemp3=momdat(i+2) !/cntdat(i+2)
                    momdat(i)= mutemp1
                    momdat(i+1)= mutemp2                    
                    momdat(i+2)= mutemp3 - mutemp1 * mutemp2   
                end if 
            end if
        end do 
		deallocate(dat)
	end if 
	call getgauss
	call getppsq 
	call getppcq
	call getppsx 
	call getppcx
	call getppmeet    
    
    
	if (skriv) call yaz0	
    timing(1)=secnds(0.0)
    call solve		
    timing(2)=secnds(timing(1))
    timing(3)=secnds(0.0)
    allocate(dat(mnad:mxa,nsim)) !ag 110416: this is allocated mna-1 rather than mna because simulation is changed to have sim(ia-1,r) at the beginning of the sim loop rather than sim(ia,r) in order to not have 0 emp at age 18
	call simulate(init,dat)
    timing(4)=secnds(timing(3))
    timing(5)=secnds(0.0)    
    if (groups) then 
		call get_mom(dat,nsim,mymom,mycnt,myvar,name,header,headerloc,momwgt)
        !deljan03
        !if ( mysay==0.and.iter==3) then 
        !    do i=1,nmom
        !        write(*,'("i,mymom,momdat,mycnt,q ",I4,2F10.2,I8,F10.2)') i,mymom(i),momdat(i),mycnt(i),q
        !    end do
        !end if        
        !deljan03
	else
		call get_mom(dat,nsim,momsim,cntsim,varsim,name,header,headerloc,momwgt)
	end if 
	deallocate(dat)
    timing(6)=secnds(timing(5))
	if (groups) then 
		call mpi_comm_rank(comm,mycommrank,mpierr)
		call mpi_allreduce(mycnt,cntsim,nmom,mpi_integer,mpi_sum,comm,mpierr)
		call mpi_allreduce(mymom*mycnt,momsim,nmom,mpi_real8,mpi_sum,comm,mpierr)
        !deljan03
        !if ( mysay==0.and.(iter==3.or.iter==2)) then         
        !    do i=1,nmom
        !        write(*,'("mysay,iter,i,momsim,momdat,cntsim,q ",2I4,I8,2F10.2,I8,F10.2)') mysay,iter,i,momsim(i),momdat(i),cntsim(i),q
        !    end do
        !end if 
        !deljan03        
        !print*, 'mygroup,mysay,mycnt,cntsim',mygroup,mysay,mycnt(112),cntsim(112)
	end if 

    !allocate(Iamtrying(mna:mxa))
    !Iamtrying=2
    !Iamtrying(mxa)=1
    !print*, 'Here it is',Iamtrying(1),Iamtrying(mxa)
    !deallocate(Iamtrying)
    
	do i=1,nmom
		!if (vardatamom(im)>0d0.or.countdatamom(im)>=datacountbar) !this "or" does not make any sense! !cohabitation correction 
		if (vardat(i)>0.0_dp .and. cntdat(i)>=datcountbar ) then
			msm_wgt(i)=vardat(i)**(-1) !1.0_dp ahu 041219    !AHU JAN19 012919
		else 
			msm_wgt(i)=0.0_dp
		end if 
		!AG090122 AGSEPT2022
		!IF WANT TO COMPARE OBJVAL BETWEEN RUNS GROUPS=TRUE AND GROUPS=FALSE THEN 
		!NEED TO RUN GROUPS=TRUE WITH THE BELOW AND THEN GROUPS=FALSE WITH THE BELOW 
		!INSTEAD OF THE LATTER IF STATEMENT (RIGHT AFTER IT)
		!BECAUSE OTHERWISE WHEN GROUPS=FALSE, IT DOES NOT CALCULATE THE MOMENTS 
		!WITH CALCVAR=1 AND CALCORR=1 THE WAY THEY ARE CALCULATED WHEN GROUPS=TRUE 
		!(I.E. THE WAY THEY ARE CALCULATED WITHIN THE IF STATEMENT WITH MUTEMP'S)
		!BUT IF COMPARING IS NOT THE PURPOSE THEN WHEN GROUPS=TRUE, 
		!THE CALCVAR AND CALCORR MOMENTS NEED TO BE CALCULATED WIHT THOSE MUTEMP'S 
		!I AM NTO JUST GETTING RID OF THESE MOMENTS BECAUSE I NEED THEM (OR THINK SO)
		!if (groups) then 
        !        if ( cntsim(i) > 0) then
		!		    momsim(i)=momsim(i)/cntsim(i)		
        !        end if 
		!end if 
		if (groups) then 
            if (calcvar(i)==0 .and. calcorr(i)==0 ) then
                if ( cntsim(i) > 0) then
				    momsim(i)=momsim(i)/cntsim(i)		!else ; simom(i)=0.0_dp  
                end if 
            else if (calcvar(i)==1) then 
                if ( cntsim(i) > 0) then
                    mutemp1=momsim(i)/cntsim(i)
                    mutemp2=momsim(i+1)/cntsim(i+1)
                    momsim(i)= mutemp1   
                    momsim(i+1)= mutemp2  - mutemp1**2
					!ag090122 agsept2022 : 
					!IF THIS IS NOT COMMENTED OUT:
					!GROUPS TRUE AND FALSE They give the same moments so they should give same objval. 
					!But they don't give same objval mainly because of this calculation above. 
					!When groups false, the momsims dont need to be dividded by cntsims and all the momsims 
					!are already calculated the way they should be by getmom (when groups false).
					!So when groups false, we do not get in these calculations within the if (groups) statement above. 
					!In this if statement, some momsims are calculated with these mutemp's (when calcvar=1 and calcorr=1) (when groups true) 
					!so there is no way to this mutemp calculation when groups true. (for calcvar=1 and calcorr=1 cases)
					!so then instead of trying to figure out how to do the calcvar and calcorr calcuations for groups true
					!i just comment them out and compare objval to groups false. 
                end if 
            else if (calcorr(i)==1) then 
                if ( cntsim(i) > 0) then
                    mutemp1=momsim(i)/cntsim(i)
                    mutemp2=momsim(i+1)/cntsim(i+1)
                    mutemp3=momsim(i+2)/cntsim(i+2)
                    momsim(i)= mutemp1
                    momsim(i+1)= mutemp2                    
                    momsim(i+2)= mutemp3 - mutemp1 * mutemp2   
                end if 
            end if
		end if 
	end do
	qcont=momwgt*msm_wgt*(momdat-momsim)**2
	objval=sum(qcont) 
    !if (iter==1) print*, 'my name is ',mysay,' and iwritegen is ',iwritegen
	if (iwritegen==1) then ; write(*,'("iter,obj: ",i6,f20.2,3f14.2)') iter,objval,timing(2),timing(4),timing(6) ; end if  
	!ahu 0317 write(*,'("iter,obj: ",3i6,f20.2,3f14.2)') mygroup,mysay,iter,q,timing(2),timing(4),timing(6)  
    
	! save the moments and objective function values from the first iteration, for comparison to the later ones: 
	t=min(numit,iter)       ! ahu 030517: will always have numit at most 2 (i.e.1 or 2). if you want more comparisons (i.e. numit>2) you have to change this. 
                            ! if numit is 1 and this is coded as t=min(2,iter), then we get out of bounds run time error for momdat_save(.,t) etc.
                            ! so be careful. but setting t=min(numit,iter) should solve that problem and also never having t=iter case as in the if statement below.
    !if (optimize .or. chkstep) then 
	!	t=min(numit,iter)
	!else 
	!	t=iter
	!end if 
	!if (t>numit) then ; print*, "error: iter>numit!!!" ; stop ; end if 
	momdat_save(:,t)=momdat
	vardat_save(:,t)=vardat
	cntdat_save(:,t)=cntdat
	momsim_save(:,t)=momsim
	cntsim_save(:,t)=cntsim
	qcont_save(:,t)=qcont
	q_save(t)=objval
	par_save(:,t)=parvec 
	realpar_save(:,t)=realparvec
	if (iwritegen==1) then
		i=maxloc(abs(par_save(:,t)-par_save(:,1)),1)
	    write(63,'(i4,". ",1a20,2f20.2,4f14.4)') i,parname(i),q_save(1),q_save(t),par_save(i,1),par_save(i,t),realpar_save(i,1),realpar_save(i,t)
    end if 
	if ((.not.optimize).or.(optimize.and.objval<best)) then	
		best=objval
		if (iwritegen==1) call writemoments(objval) 
	end if 
    
    deallocate(decm0_s,decf0_s)
    deallocate(decm_s,decf_s)
    deallocate(vm,vf)
    deallocate(dec_mar)
    deallocate(vm0_c,vf0_c)
	deallocate(vm0ctemp,vf0ctemp)

    
	iter=iter+1	
	end subroutine objfunc

	! subroutine to be called by parallel simplex optimization whenever it's time to check if have new best point
	! note that called by master. master will have current vale of 'best' from the last time this was called.
	subroutine writebest(parvector)
		real(dp), dimension(npars), intent(in) ::parvector ! transformed vector of parameters
		integer :: i
		!if (qwrite<best) then
			write(61,*) 'Found a better point'
			do i=1,npars ; write(61,*) parvector(i) ; end do
			open(unit=66,file='bestpar.txt',status='replace')
			do i=1,npars ; write(66,*) parvector(i) ; end do
			close(66)
		!end if 
	end subroutine

	subroutine writemoments(objval)
	real(8), intent(in) :: objval
	integer(i4b) :: i,t,ihead,j,k,trueindex
	open(unit=60, file=momentfile,status='replace')
    !open(unit=61 change this 61 to another number since bestval is also 61 maybe among other things, file=momentonlyfile,status='replace')
	do i=1,npars
		!write(60,'(1a15,4f12.4)') parname(i),realpar_save(i,1:4)
        write(60,'(1a15,3f10.2)') parname(i),realpar_save(i,1:numit) !ahu030622
	end do 
    write(60,*)
    write(60,'(tr2,"np",tr1,"np1",tr1,"np2",tr2,"nl",tr1,"neduc",tr2,"nexp ",tr2,"nkid",tr5,"nqs",tr6,"nq",tr6,"nx",tr5,"nxs",tr2,"nepsmv")') !ahumarch1122
	write(60,'(4i4,3(2x,i4),4i8,i4)') np,np1,np2,nl,neduc,nexp,nkid,nqs,nq,nx,nxs,nepsmove
	write(60,*) 
	write(60,'(tr2,"nz",tr2,"nh",tr1,"ncs",tr2,"nc",tr5,"ndata",tr3,"nsimech",tr6,"nsim",tr2,"ndataobs",tr6,"nmom")') 
	write(60,'(4i4,5i10)') nz,nh,ncs,nc,ndata,nsimeach,nsim,ndataobs,nmom
	write(60,*) 
	write(60,'("wage grid w:")')    !     tr6,"m(1)",tr6,"m(2)",tr6,"m(3)",tr6,"m(4)",tr6,"m(5)",tr4,"h(1)",tr4,"h(2)")') 
	write(60,*) wg(:,1)   !       ,mg(:),hgrid(:)
	write(60,*) 
    write(60,*) wg(:,2)
	write(60,*) 
	write(60,'("wage grid weights wgt:")')    
	write(60,*) wgt(:)    
    write(60,*)
    write(60,'("moveshock grid:")')    
	write(60,*) moveshock_m(:)          
	write(60,*) moveshock_f(:)          
    write(60,*)
    write(60,'("bshock grid:")')  !ahumarch1122  
	write(60,*) bshock_m(:)       !ahumarch1122   
	write(60,*) bshock_f(:)       !ahumarch1122   
    write(60,*)
    write(60,'("ppso:")')    
	write(60,*) ppso(:)    
    write(60,*)
    write(60,'("mar grid:")')    
	!write(60,*) mg(:,1)          
	!write(60,*) mg(:,2)          
	!write(60,*) mg(:,3)          
	!write(60,*) mg(:,4)  
    !do trueindex=1,ninp
    !    call index2cotyphome(trueindex,i,j,k)
	!    write(60,*) "trueindex,co,typ,home",trueindex,i,j,k
    !    write(60,*) mg(:,trueindex) 
    !end do 
   ! write(60,*)
    write(60,*)
    write(60,'("mar grid weights wgt:")')    
	write(60,*) mgt(:)    
    write(60,*)
    write(60,'("hgrid:")')    
	write(60,*) hgrid(:)          
    write(60,*)
    write(60,'("wgts:  ",1("move",tr4,"hour",tr4,"wage",tr4),tr1,"rel",tr5,"kid")') 
	write(60,'(3x,f8.2,5f8.2)') wmove,whour,wwage,wrel,wkid
	write(60,*) 
	write(60,'(tr2,"groups",tr3,"nhome",tr2,"nhomep",tr2,"onlysingles",tr2,"optimize",tr3,"chkstep",tr5,"skriv",tr1,"numit")') 
	write(60,'(2x,L6,2(3x,I4),7x,L6,3(4x,L6),I6)') groups,nhome,nhomep,onlysingles,optimize,chkstep,skriv,numit
	write(60,'(tr2,"nonneg",tr2,tr10,"eps2",tr11,"eps",tr5,"nonlabinc")') 
	write(60,*) nonneg,eps2,eps,nonlabinc
    write(60,*)
    write(60,'("objective function:")') 
	write(60,*) q_save(:)
	write(60,*) 
    write(60, '(35x,tr7,"sim",tr7,"dat",tr7,"obj",tr7,"dif",5x,tr2,"countsim",tr2,"countdat",tr4,"vardat" )' ) 
    ihead=1   
    do i=1,nmom
		if (headerloc(ihead)==i) then
			write(60,*)
            write(60,*)
			write(60,'(a120)') header(ihead)
			ihead=ihead+1
		end if
		do t=1,numit
            if (condmomcompare) then
                if (t<=2) then
    			    !write(60,'(1i5,". ",1a23,4f10.4,5x,2i10,f10.4)')	i,name(i),momsim_save(i,t),momdat_save(i,t),qcont_save(i,t),qcont_save(i,t)-qcont_save(i,1),cntsim_save(i,t),cntdat_save(i,t),vardat_save(i,t)
                else 
                    !write(60,'(1i5,". ",1a23,4f10.4,5x,2i10,f10.4)')	i,name(i),momsim_save(i,t),momdat_save(i,t),qcont_save(i,t),qcont_save(i,t)-qcont_save(i,3),cntsim_save(i,t),cntdat_save(i,t),vardat_save(i,t)
                end if 
            else 
                write(60,'(2i5,". ",1a23,2F14.4,2f14.4,5x,2i10,f14.4)')	i,t,name(i),momsim_save(i,t),momdat_save(i,t),qcont_save(i,t),qcont_save(i,t)-qcont_save(i,1),cntsim_save(i,t),cntdat_save(i,t),vardat_save(i,t)
                !write(61,'(8i5,2f20.4)')	i,t,mominfo(0:5,i),momsim_save(i,t),momdat_save(i,t)
            end if 
        end do
		!write(60,*)
	end do
	close(60)
    !close(61)
	end subroutine writemoments
end module objf

!msm_weights(im)=(countdatamom(im)**2)*vardatamom(im)**(-1) ! weighting vector (main diagonal of diaginal weighting matrix)
!making a slight modification to the above. multiplying by (n1/n)**2 rather than n1**2 (where n1 is the countdatamom of first moment here). 
!otherwise the weights blow up, as countdatamom is usually in the thousands.   
!ahu 073112 msm_weights(im)=(real(countdatamom(im))/real(ndata))*vardatamom(im)**(-1) ! weighting vector (main diagonal of diaginal weighting matrix)
!ahu 081912 msm_weights(im)=(1d0/real(ndata))*vardatamom(im)**(-1) !ahu 073112
!ahu 061813 msm_wgt(i)=1.0_dp !(real(countdatmom(i))/real(ndata))*vardatmom(i)**(-1) 
!ahu 061713 msm_wgt(i)= (real(countdatmom(i))/real(ndataobs))  *  ( vardatmom(i)**(-1) ) !ahu 061413. now using ndataobs instead of ndaata. need to check this stuff. still! 

!	if (iwritegen==1 .and. optimize .and. (obj<best.or.chkobj==2)) then
!		write(61,*) 'found a better point', obj
!		do i=1,npars ; write(61,'(f15.10)') parvec(i) ; end do
!		open(unit=66,file='parbest.txt',status='replace')
!		do i=1,npars ; write(66,'(f15.10)') parvec(i) ; end do
!		close(66)
!		best=obj	
!	end if 

		!if (datmom(i)<0.01_dp) then 
		!		datmom(i)=1000.0_dp *datmom(i)
		!		simom(i)=1000.0_dp *simom(i)
		!else if (datmom(i)>=0.01_dp.and.datmom(i)<=0.1_dp) then 
		!		datmom(i)=100.0_dp *datmom(i)
		!		simom(i)=100.0_dp *simom(i)
		!else if (datmom(i)>0.1_dp.and.datmom(i)<=1.0_dp) then 
		!		datmom(i)=10.0_dp *datmom(i)
		!		simom(i)=10.0_dp *simom(i)
		!else if (datmom(i)>4.0_dp) then 
		!		datmom(i)=0.1_dp *datmom(i)
		!		simom(i)=0.1_dp *simom(i)
		!end if 

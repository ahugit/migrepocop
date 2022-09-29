!headloc(ihead)=im; headstr(ihead)='labor market hours by gender/rel/ia';ihead=ihead+1
!ahu 061211: have to control for ia here because the two brs have different ia compositions
!ahu 061211: br 2 has no hours/kids/cohmar simultaneously in the biannual years so if you condition on all that you will just get something until they are ia 28 or something (depending on what the br grouping is)
!ahu 061211: and so if we don't control for ia, it looks as if br 2 females who are cohabiting have decreased their hours of work. but this is just a composition effect.
!ahu 061211: excluding ia 20 because, something looks weird. br 2 works too few hours at ia 20 (for females,coh,nokid). so then when i include them, it looks as if br 2 coh females with no kids work less in the later br. 

module mom
	use params
	use share
	use myaz 
	use sol, only: getdec !ag090122 agsept2022
	implicit none
contains

FUNCTION random(iseed)
! When first call, iseed must be a large positive integer.
! iseed will be changed when exit and be used for next calling.
! The range of the generated random number is between 1 and -1
!
implicit none
integer, intent(inout) :: iseed
real :: random
!
iseed = mod(8121*iseed+28411, 134456) ! 0 =< iseed < 134456
random = real(iseed)/134456. ! 0 < random < 1
!
end FUNCTION random

	! read psid data, save into global structure psiddata, calculate moments
	! notes on data
	!   record total number of observations in ndataobs
	!	need ages of observations to be strictly increasing
	!	ids should go from 1 to ndata
	!	data entry should look like: idnum, age sex rel kids edm edf incm incf hhm hhf ddm ddf rellen with <0 being blank
	!	what can be missing? anything.either entire observation is missing (rel<0)
	!		or observation is there but no record of housework (if year \=81) (housework hours <0)
	!		if not working, wage will be missing
	! ahu 071712: you can see in familymig.do that everyone starts in the data at age 16 (i drop the other people, but that's not many people)
	! this is because i wanted their location at age 16, because that's my definition of home location
	! keep this in mind, if this ever changes 
	! only want their information after they have completed school (whether high school or college)	!ahu 071712 
	! in the simulation this is automatically taken care of by the fact that the simulation starts at age age0sim(edsim(r))
	subroutine read_actualdata(init,dat)
	type(initcond), dimension(ndata), intent(out) :: init
	type(statevar), dimension(mnad:mxa,ndata), intent(out) :: dat	!data set. first entry is ia index, second observation number
	integer(i4b) :: kk,id,age,cohort,sexr,rel,kid,edr,edsp,hhr,hhsp,rellen,loc,homeloc,minage,endage,nomiss,ierr,checkminage(ndata)
	real(dp) :: wr_perhour,wsp_perhour
	dat=ones 	            ! initialize
    init=ones_init      ! initialize
    checkminage=1
	open(unit=77,file=datafilename)
	do kk=1,ndataobs
		read(77,*) id, age, cohort, sexr, rel, kid, edr, edsp, wr_perhour, wsp_perhour, hhr, hhsp, rellen,loc,minage,endage,homeloc
		if (rel/=0.and.rel/=1.and.rel/=-1) then ; print*, "data has other relationship states! ",rel ; stop ; end if !ahu 102112: count cohabiting people as married
		!if (age==mna) hme=loc      !ahu 022517: changing home back to state grew up. 
                                    !because otherwise I have to drop those who I don't observe starting from age 18 and 
                                    !then get a different cohort composition than the original sample and the avg wage per loc's end up being much smaller 
                                    !than the original numbers. 
                                    !see more details on the reason in the explanations for ahu 022517 in familymig_2.do
		dat(age,id)%id=id           !initcond		
        dat(age,id)%co=cohort       !initcond

        dat(age,id)%sexr=sexr       !initcond
		dat(age,id)%hme=homeloc     !initcond
		dat(age,id)%endage=endage   !initcond     
		dat(age,id)%edr=edr         !initcond
		dat(age,id)%expr=-99        ! (there is no exp variable in the actual data so this is always -99)
		if (kid==0) then
            dat(age,id)%kidr=1 
        else if (kid==1) then
            dat(age,id)%kidr=2 
        else if (kid==-1) then 
            dat(age,id)%kidr=-99
        else 
            print*, 'something wrong with kid in data'
            stop
        end if 
        !kid  !min(kid,maxkid)      !initcond (it should be just 0 in the beginning but might not be in the actual data so just read it from the data here as initcond)
		!ahu 030217 if (age==mna) init(id)=dat(age,id)%initcond 
        if (age==minage) then
            checkminage(id)=0
            init(id)=dat(age,id)%initcond
        end if 
        
        !ahu jan 19 010219
        !initial conditions are the state variables of age 17,which are the state variables that agents use inorder to make decisions at the beginning of age 18
        !in order to have moving rates in the data at age 17 (because we do have sim rates at age 17 in simulation since we record initconditions at age 17 and then their decisions become the age 18 variables)
        !if (age==mna) then 
        !    dat(age-1,id)%initcond=init(id)
        !    dat(age-1,id)%rel=0
        !    dat(age-1,id)%l=loc
        !    dat(age-1,id)%hr=0
        !end if 

		dat(age,id)%l=loc
		call get_dathrwge(hhr,wr_perhour,dat(age,id)%hhr,dat(age,id)%wr,dat(age,id)%logwr)        
        
		dat(age,id)%rel=rel
		if (rel==1) then 
			dat(age,id)%rellen=rellen
		else if (rel==0) then 
			dat(age,id)%rellen=-99 !99
		end if 
		
		if (rel==1) then 
			call get_dathrwge(hhsp,wsp_perhour,dat(age,id)%hhsp,dat(age,id)%wsp,dat(age,id)%logwsp)
		else if (rel==0) then 
			dat(age,id)%hhsp=-99 !99
			dat(age,id)%wsp=-99.0_dp !99.0_dp
            dat(age,id)%logwsp=-99.0_dp
		end if 

		dat(age,id)%edsp=edsp    !for some reason, don't read edsp from the data yet. just set it to -99 and do read it later. ahu 021617: now I read it!  
		dat(age,id)%expsp=-99   !no experience variable in the data 
		dat(age,id)%kidsp=-99    !kidsp is just a simulation concept
				
		call getmiss(dat(age,id),nomiss)
		dat(age,id)%nomiss=nomiss
		
		!some checks:
		if (edr.ne.1.and.edr.ne.2) then
			print*, 'There is something wrong with edr in read_actual data'
			stop 
		end if 
		if (dat(age,id)%co<0.or.dat(age,id)%sexr<0) then ; print*, "cohort and sex are negative!", age,dat(age,id)%co,cohort,dat(age,id)%sexr,sexr ; stop ; end if 

        !why is this age<agestart-1 and not age<agestart
        !because for ed, agestart is 22 in sim the initial conditions are assigned to age 21 and then 
        !the move rates at age 21 for ed people in sim turn out to be very large because they move immediately from their starting 
        !home location to the best uloc location. In the data, this moment is not there because of the setting of dat=ones 
        !for ages age<agestart. In order to let the estimation to its thing and be able to compare that large sim moving rate 
        !at age 21 to the data, I am setting the dat to ones for ages 0-20 rather than 0-21. doing this is also more consistent because 
        !now (i.e. when age<agestart case) the simulation has age 21 for the ed people but not the data. 
        !I am not doing this for noed though since for noed people there is no age read from the data that is less than mna
        !check this
        if (age<agestart(edr)-1) then 
			dat(age,id)=ones
		end if
        if (age<mna) then 
            print*, 'There is something wrong with age in read_actual data'
            stop
        end if 
    enddo 
    if (sum(checkminage)>0) then !this would mean that we don't get age=minage for some people which is not right
        print*, 'something is wrong with checkminage',checkminage
        stop
    end if
    if (minval(init(:)%id)<0.or.minval(init(:)%co)<0.or.minval(init(:)%sexr)<0.or.minval(init(:)%hme)<0.or.minval(init(:)%endage)<0.or.minval(init(:)%edr)<0 ) then 
        print*, 'something is wrong with init everything' !because this would mean things did not get set at minage. either minage was not around (which can't be right) or something else. 
        stop
    end if
	close(77)
	end subroutine read_actualdata
	
	subroutine simulate(init,sim)
	type(initcond), dimension(ndata), intent(in) :: init
	type(statevar), dimension(mnad:mxa,nsim), intent(out) :: sim
	type(shock), allocatable, dimension(:,:) :: epsim
	integer(i4b) :: q0,x0,q,x,dec(3),draw(2)
	integer(i4b) :: relnext,qnext,xnext
	integer(i4b) :: i0,n0,i,n,z
	integer(i4b) :: rel0,index,trueindex,ia,r,g,endage
	integer(i4b) :: nn,mm,typsim,hh(2),l(2),nomiss
	integer(i4b) :: p,k,imax				! for random_seed
	real(dp) :: rand,vmax(2),logw(2),inc(2),wage(2),dum(2) !dum is just a dummy variable for yaz_checknb
	logical  :: welldef,newrel0,meet,checkgroups,defj(nc)
	integer(i4b) :: qmatch,xmatch,dd(11),iepskid,iepsmove,eddummy,expdummy,ww(2),ed(2),expe(2),kid(2)
	real(dp) :: vec(5),mcost(2),surplus,surplusmax,surplusj(nc),val(2),vcheck(2),transfers(2)
	integer(i4b) :: l0,j,jmax,qmax,relmax,de(1)
	!integer, allocatable :: newseed(:)
	!integer ::     seedo
	integer(i4b), dimension(12) :: seed_init !98765
    integer(i4b) :: callfrom !ag090122 agsept2022
    real(dp) :: valso(2)

     !seedo=98765
	!seedo=seed_init    
	allocate(epsim(mna:mxa,nsim))
	
	!if (iter==1) then
	    call random_seed (size=p)
	    !p=12
	    call random_seed (put = (/(k,k=1,p)/))
       !seed_init=3
       !call random_seed (put = seed_init )
	    if (skriv) then
	        print*, 'Here is p',p
	        print*, 'Here is k',(/(k,k=1,p)/)
	    end if
	!else     
	!    allocate(newseed(p))
	!    newseed=10
	!    call random_seed( put=newseed(1:p)  )
	!    deallocate(newseed)
	!end if 

    !call random_seed()
    !call random_seed(size=seed_size) 
    !print*, "seed size for random seed is:", seed_size
    !allocate(seedonzo(seed_size))
    !call random_seed(get=seedonzo)
    !print*, "system generated random seed is:"
    !print*, seedonzo
    !seedonzo=300000
    !call random_seed(put=seedonzo)
    !print*, "my generated random seed is:"
    !print*, seedonzo
    !call random_number(rhino)
    !print*, "rhino:", rhino 
    !deallocate(seedonzo)
    
    !idumo=-1
    !rhino=ran(idumo) 
    !print*,'here is rhino today', iam,rhino
    !rhino=ran(idumo) 
    !print*,'here is rhino tomorrow', iam,rhino
   
		
	sim=ones 	! initialize to smallest valid value
	r=0
	do nn=1,ndata			! person being simulated !ahu 031911 big change	
		do mm=1,nsimeach	! number of simulations for this person
			r=r+1
			do ia=mna,mxa
				call random_number(rand)
				epsim(ia,r)%q=rand !random(seedo)
                call random_number(rand)
				epsim(ia,r)%x=rand !random(seedo) 
				call random_number(rand)
				epsim(ia,r)%marie=rand !random(seedo) 
				call random_number(rand)
				epsim(ia,r)%meet=rand !random(seedo) 
				call random_number(rand)
				epsim(ia,r)%meetq=rand !random(seedo)
				call random_number(rand)
				epsim(ia,r)%meetx=rand !random(seedo) 
				call random_number(rand)
				epsim(ia,r)%iepsmove=rand !random(seedo) 
				!call minus1(sim(ia,r))		!initialize to smallest valid value
                call random_number(rand)
                epsim(ia,r)%typ=rand !random(seedo)
			end do 
		end do 
	end do 
	!print*, 'epsim',epsim(18,5)%q,epsim(18,5)%x
	!call random_number(rand) 
    !print*, 'rand',rand
	!if (skriv) print*, 'epsim',epsim(18,5)%q,epsim(18,5)%x
	!if (skriv) call random_number(rand) 
    !if (skriv) print*, 'rand',rand
    !if (skriv) print*, 'random(seedo)',random(seedo),random(seedo)
	!if (skriv) print*, 'random(seedo)',random(seedo),random(seedo)

    !if (skriv) print*, 'random(seedo)',random(seedo),random(seedo)
	!if (skriv) print*, 'random(seedo)',random(seedo),random(seedo)
	
    checkgroups=.true.
	r=0								! overall count of total number of simulations
    typsim=-1
	iddat: do nn=1,ndata			! person being simulated !ahu 031911 big change	
		idsim: do mm=1,nsimeach		! number of simulations for this person
			r=r+1
					!sim(mna,r)%initcond=init(nn)   !co, sex, hme, endage   !ahu0115del

			if (skriv.and.r==1347) then ; yaz=.TRUE.  ; else ; yaz=.FALSE. ; end if 

        
		    IF (init(nn)%edr==1) THEN ! low ed, pick type based on ptypehs(co)
			    typsim=multinom(ptypehs(:),epsim(MNA,r)%typ) 
		    ELSE IF (init(nn)%edr==2) THEN  ! high ed type, pick type based on ptypecol(co)  
                typsim=multinom(ptypecol(:),epsim(MNA,r)%typ) 
            ELSE 
                print*, 'something is wrong with type',init(nn)%edr,typsim,ptypehs,ptypecol,epsim(MNA,r)%typ
                stop
		    ENDIF
            if (ntypp==1.and.typsim.ne.1) then
                print*, 'something wrong with ntypp and type',ntypp,typsim
                stop
            end if 
            !typ=1
            	                    
            if (init(nn)%edr.ne.1.and.init(nn)%edr.ne.2) then ; print*, 'There is something wrong with init edr' ; stop ; end if  
			call cotyphome2index(trueindex,init(nn)%co,typsim,init(nn)%hme)
            !index is set to 1 if groups, because each processor has its own value function and so that they don't need to be declared with dimension ninp
            if (groups) then 
                index=1
            else 
                index=trueindex
            end if
			!if (skriv.and.trueindex==1) then ; yaz=.TRUE.  ; else ; yaz=.FALSE. ; end if !ahu 0327 trueindex 2
			!if (skriv) then ; yaz=.TRUE.  ; else ; yaz=.FALSE. ; end if !ahu 040918 del and remove later. the above yaz was the original yaz
            if (typsim<0) then 
                print*, 'typ did not get assigned'
                stop
            end if 
			ind: if ( (.not.groups) .or.  (myindex==trueindex) ) then !(myco==init(nn)%co.and.mytyp==typsim.and. myhome==init(nn)%hme)  ) then	  
                !if (checkgroups) then 
                !    write(*,'("here I am",7I4)') mysay,init(nn)%co,typsim,init(nn)%hme,myco,mytyp,myhome
                !    checkgroups=.false.
                !end if
                !write(*,'("Here I am", 6I4)') index,trueindex,myindex,init(nn)%co,typsim,init(nn)%hme

				if (init(nn)%edr.ne.1.and.init(nn)%edr.ne.2) then ; print*, 'There is something wrong with init edr' ; stop ; end if  
				rel0=0 
				newrel0=.false.
                q0 = wl2q(np1,init(nn)%hme)
				x0 = erkid2x( init(nn)%edr, 1, 1) !ed,exp which is 1 in the beginning and kidr which is 1 in the beginning (for kid, 1 means no kid) 
				g=init(nn)%sexr
				endage=min(init(nn)%endage,mxa)				
                age: do ia=agestart(init(nn)%edr),endage  !mna,endage
        			!if (skriv.and.trueindex==1.and.ia==18) then ; yaz=.TRUE.  ; else ; yaz=.FALSE. ; end if !ahu 0327 trueindex 2 !ahu 040918 del and remove later
                    !if (skriv.and.ia<=20) then ; yaz=.TRUE.  ; else ; yaz=.FALSE. ; end if !ahu 0327 trueindex 2 !ahu 040918 del and remove later
					call qx2hrwge(g,rel0,q0,x0,trueindex,hh,l,wage,logw)
					if (rel0==0) call x2edexpkid(x0,ed(g),expe(g),kid(g))
					if (rel0==1) then
						ed(:)=xx2e(:,x0)
						expe(:)=xx2r(:,x0)
						kid(:)=xx2kid(:,x0)
					end if
                    if (ed(g).ne.init(nn)%edr) then ; print*, 'something wrong with ed and init%ed' ; stop ; end if
					if (q0<=0) print*, "q0 is no good ",r,ia-1,init(nn)%hme,wl2q(np1,init(nn)%hme)                    
					!if (yaz) then ; write(400,'(4/,2I8)') r,ia ; write(400,'("State Variables:")') ; call yaz_sim(g,rel0,q0,x0) ; end if 
                    if (yaz) then ; write(400,*) epsim(ia,r)%q,epsim(ia,r)%x,logw(g) ; end if !ahu 012019 del
					!!ag 110416: changed to have sim(ia-1,r) at the beginning of the sim loop rather than sim(ia,r) in order to not have 0 emp at age 18
					sim(ia-1,r)%initcond=init(nn)   !co, sex, hme, endage
                    sim(ia-1,r)%typ=typsim   !co, sex, hme, endage
                    if (sim(ia-1,r)%edr.ne.ed(g)) then ; print*, 'something wrong with ed' ; stop ; end if
					sim(ia-1,r)%expr=expe(g)
					sim(ia-1,r)%kidr=kid(g)                    
					sim(ia-1,r)%hhr = hh(g)
                    sim(ia-1,r)%wr = wage(g)
                    sim(ia-1,r)%logwr = logw(g)
                    sim(ia-1,r)%l = l(g)
					sim(ia-1,r)%rel=rel0
                    if (rel0>0) then                      ! r is male/sp is female or r is female/sp is male
						i = one(g==1) * 2 + one(g==2) * 1 ! r is male/sp is female or r is female/sp is male
						if (newrel0 .or. ia==agestart(init(nn)%edr) ) then	!agestart(simdat(mna,nn)%edr)   ) then		!if (newrelsimsimdat(ia,r) .or.(ia==minage)) then correct this in cohabitation also. cohabitation correction. 
							sim(ia-1,r)%rellen=1
						else
							sim(ia-1,r)%rellen=sim(ia-2,r)%rellen+1
						endif                        
						sim(ia-1,r)%edsp=ed(i)
						sim(ia-1,r)%expsp=expe(i)
						sim(ia-1,r)%kidsp=kid(i)                                            
                        sim(ia-1,r)%hhsp   = hh(i)
                        sim(ia-1,r)%wsp = wage(i)
                        sim(ia-1,r)%logwsp = logw(i)
						sim(ia-1,r)%lsp = l(i)  !just for checking purposes
                        if (l(g).ne.l(i)) then ; print*, 'lg and li not equal' ; stop ; end if 
					else if (rel0==0) then 
						sim(ia-1,r)%rellen=-99		
						sim(ia-1,r)%edsp=-99
						sim(ia-1,r)%expsp=-99
						sim(ia-1,r)%kidsp=-99                                            
                        sim(ia-1,r)%hhsp   = -99	
                        sim(ia-1,r)%wsp = -99.0_dp
						sim(ia-1,r)%logwsp = -99.0_dp
						sim(ia-1,r)%lsp = -99  !just for checking purposes
					endif
					sim(ia-1,r)%nn = nn
					sim(ia-1,r)%mm = mm
					sim(ia-1,r)%r = r
					call getmiss(sim(ia-1,r),nomiss)
					sim(ia-1,r)%nomiss=nomiss
					!if (yaz) call yaz_simpath(ia-1,nn,mm,r,sim(ia-1,r))
					!if (ia==47 .and. sim(ia-1,r)%hhr ==1) then 
                    !    print*, 'I found it!',ia-1,sim(ia-1,r)%co,sim(ia-1,r)%sexr,sim(ia-1,r)%rel
                    !    stop
                    !end if 
                        
                    
					if (rel0==0) then
						q = multinom( ppsq(:,q0,x0,g) , epsim(ia,r)%q ) 
						x = multinom( ppsx(:,q0,x0)   , epsim(ia,r)%x ) 
						meet=( epsim(ia,r)%meet<=pmeet )
						z	= multinom( mgt			, epsim(ia,r)%marie)
						iepsmove = multinom( ppso(:)   , epsim(ia,r)%iepsmove ) 
						if (g==1) dec(1) = decm_s(iepsmove,x,q,q0,ia,index) 
						if (g==2) dec(1) = decf_s(iepsmove,x,q,q0,ia,index) 
						dec(2)=x 
						relnext=0 ; qnext=dec(1) ; xnext=dec(2)    ! next period's variabeles are these unless marriage market happens:
						qmatch	= multinom( ppmeetq(:, dec(1) )	, epsim(ia,r)%meetq) 
						xmatch	= multinom( ppmeetx(:, dec(2) )	, epsim(ia,r)%meetx) 
						if (yaz) then 
                        !!!    write(400,'("Trueindex:",I4)') trueindex
						    write(400,'("Draws:",3F10.2)')  epsim(ia,r)%q, epsim(ia,r)%x, epsim(ia,r)%marie
						    write(400,'("Draws:",3I10)')  q,x,z
						    !write(400,'("Draws:")') ; call yaz_sim(g,rel0,q,x)
						    !write(400,'("Single Dec Before Mar Mkt:")') ; call yaz_sim(g,rel0,qnext,xnext)
						!!!	write(400,'("Match:")') ;  call yaz_simmatch(meet,qmatch,xmatch,z)
                        end if                         
						if (meet) then 	
							if (g==1) then ; q = q2qq(dec(1),qmatch) ; x = x2xx(dec(2),xmatch) ; end if 
							if (g==2) then ; q = q2qq(qmatch,dec(1)) ; x = x2xx(xmatch,dec(2)) ; end if 
							relnext=dec_mar(z,x,q,ia,index)
							if (relnext==1) then 
								qnext=q ; xnext=x
							else if (relnext==0) then
								qnext=dec(1) ; xnext=dec(2)
							end if 
						    if (yaz) then 
                                write(400,'("HERE IS DECMAR:",6I8)') z,q,x,ia,index,dec_mar(z,x,q,ia,index)
							    !write(400,'("Decision At Mar Mkt:")') ;  call yaz_simdecmar(relnext)
						    end if      !!!                   
                        end if !meet      
					else if (rel0==1) then 
						whereami=4		! for telling yaz where we are: Couples/Sim
						q	= multinom( ppcq(:,q0,x0)	, epsim(ia,r)%q) 
						x	= multinom( ppcx(:,q0,x0)	, epsim(ia,r)%x) 	
						z	= multinom( mgt		, epsim(ia,r)%marie)
						iepsmove = multinom( ppso(:)   , epsim(ia,r)%iepsmove ) 
                        if (yaz) then ; write(400,*) epsim(ia,r)%q,epsim(ia,r)%x,epsim(ia,r)%marie,q,x,z ; end if !ahu 012019 del
						!!!if (yaz) then 
                        !!!    write(400,'("Trueindex:",I4)') trueindex
                        !!!    write(400,'("Draws:")') ; call yaz_sim(g,rel0,q,x)
						!!!	write(400,'("z: ",I4)') z
						!!!end if                         
						!AG090122 AGSEPT2022 dd=(/ ia,trueindex,q,x,z,q0,g,-1,-1,-1,iepsmove /)			    ! (/ ia,index,q,x,z,q0,g,jmax,qmax,relmax /)  							                        
                        !ag090122 agsept2022 call getdec_c(dd,vmax)					            ! jmax=dd(8) ; qmax=dd(9) ; relmax=dd(10)   
                        valso=pen !ag090122 agsept2022
                        callfrom=80 !ag090122 agsept2022
                        dd = (/ia,trueindex,q,x,z,q0,callfrom,-1,-1,-1,iepsmove /) 	! (/ ia,index,q,x,z,q0,gender/callfrom,jmax,qmax,relmax,iepsmove /)  	                                        
                        call getdec(dd,vmax,valso)
                        relnext=dd(10)
						!!!!if (yaz) then ; call yaz_decision(dd,vmax) ; end if	    ! write down decision     
						if (relnext==1) then 
							qnext=dd(9) ; xnext=dd(4) 
						else if (relnext==0) then 
							i = qq2q(g,q) ; n=xx2x(g,x) ; i0 = qq2q(g,q0)	! translate couple q and x and q0 into singles
							if (g==1) dec(1) = decm0_s(iepsmove,n,i,i0,ia,index) 
							if (g==2) dec(1) = decf0_s(iepsmove,n,i,i0,ia,index) 
							dec(2)=n
							qnext=dec(1) ; xnext=dec(2)                   
						end if
					end if ! rel
					!!!if (yaz) then 
					!!!	write(400,'("Next: ")') ; call yaz_sim(g,relnext,qnext,xnext) 
					!!!end if 
					newrel0 = (rel0==0.and.relnext==1) 
					rel0 = relnext 
					q0 = qnext
					x0 = xnext
					relnext=-99 ; qnext=-99 ; xnext=-99 ; q=-99 ; x=-99 ; z=-99 ; qmatch=-99 ; xmatch=-99 ; dec=-99 ; meet=.FALSE.
					if (init(nn)%edr==2.and.ia<agestart(init(nn)%edr)) sim(ia,r)=ones
				end do	age
			end if ind
		end do idsim
	end do iddat
	deallocate(epsim)		
	end subroutine simulate

	subroutine qx2hrwge(g,rel,q,x,trueindex,hh,l,wage,logw)
	integer(i4b), intent(in) :: g,rel,q,x,trueindex
	integer(i4b), dimension(2) :: ed,expe
	integer(i4b), intent(out) :: hh(2),l(2)
	real(dp), intent(out) :: wage(2),logw(2) 	
	integer(i4b) :: i,w(2)
	hh=-99
    wage=-99.0_dp
	logw=-99.0_dp
	l=-99
    if (rel == 0) then			
		w(g)=q2w(q) 
		l(g)=q2l(q) 
		ed(g)=x2e(x) 
		expe(g)=x2r(x) 
		if ( w(g)<=np ) then 
			hh(g)=1
            wage(g) = ws(g,x,q,trueindex)   !fnwge( g,typ,l(g) ,wg(w(g),g) ,ed(g) ,expe(g))
			logw(g) = log( wage(g) )        !log ( fnwge( g,typ,l(g) ,wg(w(g),g) ,ed(g) ,expe(g))  )			
		else
			hh(g)=0
			wage(g) = -99.0_dp
            logw(g) = -99.0_dp
		endif
		if ( w(g)<=0.or.w(g)>np1 )  then ; print*, "error in read_dat: w<=0 or w>np1 !!" ; stop ; end if 
	else if (rel > 0) then 		
		w(:)=qq2w(:,q)
		l(:)=qq2l(:,q)
		ed(:)=xx2e(:,x)
		expe(:)=xx2r(:,x)
		do i=1,2
			if ( w(i)<=np ) then 
				hh(i)=1
				wage(i) = wc(i,x,q,trueindex)   !fnwge(i,typ,l(i) ,wg(w(i),i) ,ed(i) ,expe(i)) 		
                logw(i) = log( wage(i) )        !log ( fnwge(i,typ,l(i) ,wg(w(i),i) ,ed(i) ,expe(i)) )			
			else
				hh(i)=0
                wage(i) = -99.0_dp
				logw(i) = -99.0_dp
			endif
			if ( w(i)<=0.or.w(i)>np1 )  then ; print*, "error in read_dat: w<=0 or w>np1 !!" ; stop ; end if 
		end do 
	end if 
	end subroutine qx2hrwge
				
	subroutine getmiss(dat,nomiss)
	type(statevar), intent(in) :: dat ! data set. first entry is ia index, second observation number
	integer(i4b),intent(out) :: nomiss 
	nomiss=1
	if (dat%co<0) nomiss=0
	if (dat%sexr<0) nomiss=0
	if (dat%hme<0) nomiss=0		
	if (dat%endage<0) nomiss=0		
	if (dat%edr<0) nomiss=0
	!if (dat%expr<0) nomiss=0   !expr is alwyays missing in the data since there is no experience variable
	if (dat%rel==1.and.dat%kidr<0) nomiss=0
	if (dat%l<0) nomiss=0
	if (dat%hhr<0) nomiss=0
	if (dat%logwr<0.0_dp) nomiss=0
	if (dat%rel<0.0_dp) nomiss=0
	if (dat%rel==1.and.dat%rellen<0) nomiss=0
	if (dat%rel==1.and.dat%hhsp<0) nomiss=0    
	if (dat%rel==1.and.dat%logwsp<0.0_dp) nomiss=0
	!if (dat%edsp<0) nomiss=0   !edsp is always missing in the data because for now it is not read from the data
	!if (dat%expsp<0) nomiss=0  !expsp is alwyays missing in the data since there is no experience variable
	if (dat%kidsp<0) nomiss=0   !kidsp is always missing in the data since kidsp is just a simulation construct
	
	end subroutine getmiss 

	subroutine get_dathrwge(hr,inc_perhour,dhr,dw,dlogw)
	integer(i4b), intent(in) :: hr
	real(dp), intent(in) :: inc_perhour
	integer(i4b), intent(out) :: dhr
	real(dp), intent(out) :: dw,dlogw	
	!if (hr>0.and.wge>0) then 
        !note that the employment is decided according to whether hours is more than h_fulltime which is 1000 but annual wages is calculated using wge*h_wmult where h_wmult is 2000. 
        !This is how you did it in the original version. I do the same thing so that the wage numbers are consistent with the previous version. See page 15 last paragraph in the original text.
	!	dhr=one(hr>=h_fulltime)
	!	dwge=wge*h_wmult !ahu 021817: 1) got rid of the minmax thing 2) changed h_fulltime to h_wmult in the way annual wages is calculated      !min(max(minw,wge),maxw)*h_fulltime
	!	dwge=log(dwge)
	!	dwge=one(hr>=h_fulltime)*dwge
	!else if (hr==0.and.wge==0) then 
	!	dhr=hr 
	!	dwge=wge
	!else if (hr==-1.and.wge==-1) then 
	!	dhr=hr 
	!	dwge=wge
	!else 
	!	print*, "error in get_dathrwge ", hr,wge
	!	print*, "by construction, it should be that it's both 0, both -1 or both positive "
	!	stop
	!end if	
    
	if (hr>=0) then 
        !note that the employment is decided according to whether hours is more than h_fulltime which is 1000 but annual wages is calculated using wge*h_wmult where h_wmult is 2000. 
        !This is how you did it in the original version. I do the same thing so that the wage numbers are consistent with the previous version. See page 15 last paragraph in the original text.
		dhr=one(hr>=h_fulltime)
		if (inc_perhour>=minw.and.inc_perhour<=maxw) then
            dw=inc_perhour*h_wmult    
            dlogw=log(inc_perhour*h_wmult) !ahu 021817: 1) got rid of the minmax thing 2) changed h_fulltime to h_wmult in the way annual wages is calculated      !min(max(minw,wge),maxw)*h_fulltime            
        else 
            dw=-99.0_dp
            dlogw=-99.0_dp
        end if
    else if (hr<0) then 
        dhr=-99
        dw=-99.0_dp
        dlogw=-99.0_dp
	else 
		print*, "error in get_dathrwge ", hr,inc_perhour
		print*, "by construction, it should be that it's both 0, both -1 or both positive "
		stop
	end if	
    
    end subroutine get_dathrwge	
	
	subroutine get_mom(dat,ndat,mom,cnt,var,name,headstr,headloc,weights)
	integer(i4b), intent(in) :: ndat ! number of observations in dat array    
	type(statevar), dimension(MNAD:MXA,ndat), intent(in) :: dat ! data set. first entry is ia index, second observation number
	real(dp), dimension(nmom), intent(out) :: mom	 ! conditional mom
	integer(i4b), dimension(nmom), intent(out) :: cnt	 ! number of observations contributing to each moment
	real(dp), dimension(nmom), intent(out) :: var		 ! unconditional variance of each conditional moemnt
	character(len=namelen), dimension(nmom), intent(out) :: name ! names of mom
	character(len=500),dimension(nmom), intent(out) :: headstr ! headers for different sections of mom
	integer(i4b), dimension(nmom), intent(out) :: headloc	 ! location of headers for different sections of mom
	real(dp), dimension(nmom), intent(out) :: weights		 ! some real assoicated to each moment, maybe used as weights
	integer(i4b) :: ia,co,im,ihead,g,i,j,jj,ii,ddd
    integer(i4b) :: minsex,maxsex,maxrelo
	logical,dimension(MNAD:MXA,ndat) :: coho,cosex,cosexrel,corel
	integer(i4b), dimension(MNAD:MXA,ndat) :: norelchg,move,emph,empw,edh,edw,dur,dee,deue
	real(dp), dimension(MNAD:MXA,ndat) :: logwh,logww
    integer(i4b), dimension(ndat) :: nummove,cohogen,sexgen    !nummove_ma,nummove_si !integer(i4b), dimension(MNAD:MXA,ndat) :: nummov,nummove_mar,nummove_sin
    real(dp) :: wmovebyrel,decilegrid(ndecile)
    INTEGER(I4B),DIMENSION(MNAD:MXA,ndat) :: kidtrans,homemove,moverank
    INTEGER(I4B),DIMENSION(ndat) :: movesum
    REAL(dp),DIMENSION(MNAD:MXA,ndat) :: mean4h,deltawage4,deltawage
    INTEGER(I4B), dimension(MNAD:MXA,ndat) :: iacat,empr
    INTEGER(I4B), dimension(5,MNAD:MXA,ndat) :: etr
    LOGICAL,DIMENSION(MNAD:MXA,ndat) :: obs4h



    !initiate
	mom=-99.0_dp
	cnt=9999
	var=-99.0_dp
	headloc=-99 
	weights=0.0_dp 
    maxsex=-1
    maxrelo=-1
    coho=.FALSE.
    cosex=.FALSE.    
    cosexrel=.FALSE.
    corel=.FALSE.
    norelchg=-99
    move=-99
    emph=-99 
    empw=-99 
    edh=-99 
    edw=-99 
    dur=-99
    dee=-99
    deue=-99
    logwh=-99.0_dp
    logww=-99.0_dp
    nummove=-99 
    cohogen=-99 
    sexgen=-99 
    wmovebyrel=-99.0_dp
    calcvar=0   !declared globally
    calcorr=0   !declared globally
    mominfo=-1  !declared globally
	im=1
	ihead=1
    !decilegrid=0.0_dp
    !dur=0
    kidtrans=-99
    homemove=-99 

    headloc(ihead)=im
	!if (skriv) call yaz_getmom(dat,ndat) 
	
    do ddd=1,ndecile
        decilegrid(ddd)=8.6_dp+0.25_dp*(ddd-1)
	end do
    !print*, decilegrid
    
	WHERE ((dat(MNAD:MXAD,:)%rel>=0) .AND. dat(MNA:MXA,:)%rel>=0  )
		norelchg(MNAD:MXAD,:)=one(  dat(MNAD:MXAD,:)%rel==dat(MNA:MXA,:)%rel   ) 
	ENDWHERE
	WHERE ((dat(MNAD:MXAD,:)%l>0) .AND. (dat(MNA:MXA,:)%l>0))
		move(MNAD:MXAD,:)=one(dat(MNAD:MXAD,:)%l/=dat(MNA:MXA,:)%l)
	ENDWHERE
	WHERE ((dat(MNAD:MXA,:)%rel==1) .AND.  (dat(MNAD:MXA,:)%sexr==1)  )
		emph(MNAD:MXA,:)=dat(MNAD:MXA,:)%hhr 
		empw(MNAD:MXA,:)=dat(MNAD:MXA,:)%hhsp
		edh(MNAD:MXA,:)=dat(MNAD:MXA,:)%edr 
		edw(MNAD:MXA,:)=dat(MNAD:MXA,:)%edsp
		logwh(MNAD:MXA,:)=dat(MNAD:MXA,:)%logwr 
		logww(MNAD:MXA,:)=dat(MNAD:MXA,:)%logwsp
    ENDWHERE
	WHERE ((dat(MNAD:MXA,:)%rel==1) .AND.  (dat(MNAD:MXA,:)%sexr==2)  )
		emph(MNAD:MXA,:)=dat(MNAD:MXA,:)%hhsp
		empw(MNAD:MXA,:)=dat(MNAD:MXA,:)%hhr
		edh(MNAD:MXA,:)=dat(MNAD:MXA,:)%edsp 
		edw(MNAD:MXA,:)=dat(MNAD:MXA,:)%edr
		logwh(MNAD:MXA,:)=dat(MNAD:MXA,:)%logwsp 
		logww(MNAD:MXA,:)=dat(MNAD:MXA,:)%logwr
    ENDWHERE


    !dee for wdif conditionoing on ee transitions 
    WHERE (dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA+1:MXA,:)%logwr>=0  ) 
        dee(MNA:MXAD,:)=1
    ENDWHERE
    
    !deue for wdif conditioning on eue transitions
    WHERE ( dat(MNA:MXAD-1,:)%hhr==1 .AND. dat(MNA+1:MXAD,:)%hhr==0  .AND. dat(MNA+2:MXA,:)%hhr==1 .AND. dat(MNA:MXAD-1,:)%logwr>=0   .AND. dat(MNA+2:MXA,:)%logwr>=0  ) 
        deue(MNA:MXAD-1,:)=1
    ENDWHERE
 
	WHERE ((dat(MNAD:MXA,:)%hhr==0)   )
		dur(MNAD:MXA,:)=1
    ENDWHERE

    !kidtrans=-1.
    WHERE (  (dat(MNA:MXAD,:)%rel>0) .AND. (dat(MNA:MXAD,:)%kidr==0) .AND. (dat(MNA+1:MXA,:)%kidr>=0)  )
    kidtrans(MNA:MXAD,:)=one(dat(MNA+1:MXA,:)%kidr>0)
    ENDWHERE
    WHERE (move(MNA:MXAD,:)==1)
    homemove(MNA:MXAD,:)=one(dat(MNA+1:MXA,:)%hme==dat(MNA+1:MXA,:)%l)
    ENDWHERE

    movesum(:)=sum(move(MNA:MXAD,:),1,move(MNA:MXAD,:)>=0)
    moverank=-1
    do ia=MNA,MXAD
    moverank(ia,:)=sum(move(MNA:ia,:),1,move(MNA:ia,:)>=0)   
    end do 

    !ddd=1
    !do i=5,10
    !    print*, 'Here it is ', d1*one( (dat(25,i)%logwr>decilegrid(ddd) .and. dat(25,i)%logwr<decilegrid(ddd+1) ) )
    !end do     
    
    do i=1,ndat
        do ia=MNAD,MXAD
	        if (        (dat(ia,i)%hhr==0)    .AND. (dat(ia+1,i)%hhr==0)   ) then
		        dur(ia+1,i)=dur(ia,i)+1
            end if 
        end do 
    end do 
    
    do i=1,ndat
        !Just generating cohort and sex of someone because don't have initcond in this routine
        cohogen(i)=   maxval(dat(MNAD:MXA,i)%co)
        sexgen(i)=   maxval(dat(MNAD:MXA,i)%sexr)	                
    
        !Total number of moves overall
        !nummove(i)=sum(move(MNA:MXAI,i),  MASK=( move(MNA:MXAI,i)>0  .and. norelchg(MNA:MXAI,i)==1  )   )
        nummove(i)=sum(move(MNA:MXAD,i),  MASK=( move(MNA:MXAD,i)>=0  .and. norelchg(MNA:MXAD,i)==1  )   )
    end do 


    ddd=1
	if (skriv.and.ndat==nsim .and. iter==1) then
	open(unit=94632, file='dursim1.txt',status='replace')
    do j=1,ndat
        do ia=mnad,mxa
            if (dat(ia,j)%hme==1) then 
                write(94632,'(6i6,F10.2,3I6)' ) j,dat(ia,j)%nn,dat(ia,j)%mm,ia,dat(ia,j)%L,dat(ia,j)%hhr,dat(ia,j)%logwr,dat(ia,j)%EXPR,dat(ia,j)%HME,dat(ia,j)%sexr   !, d1*one( (dat(ia,j)%logwr>decilegrid(ddd) .and. dat(ia,j)%logwr<decilegrid(ddd+1) ) ), decilegrid(ddd),decilegrid(ddd+1)   !,dat(ia,j)%nn,dat(ia,j)%mm,dat(ia,j)%r
	        end if 
        end do 
	end do
    close(94632)
    end if 

    ddd=1
	if (skriv.and.ndat==nsim .and. iter==2) then
	open(unit=94632, file='dursim2.txt',status='replace')
    do j=1,ndat
        do ia=mnad,mxa
            if (dat(ia,j)%hme==1) then 
                write(94632,'(6i6,F10.2,3I6)' ) j,dat(ia,j)%nn,dat(ia,j)%mm,ia,dat(ia,j)%L,dat(ia,j)%hhr,dat(ia,j)%logwr,dat(ia,j)%EXPR,dat(ia,j)%HME,dat(ia,j)%sexr   !, d1*one( (dat(ia,j)%logwr>decilegrid(ddd) .and. dat(ia,j)%logwr<decilegrid(ddd+1) ) ), decilegrid(ddd),decilegrid(ddd+1)   !,dat(ia,j)%nn,dat(ia,j)%mm,dat(ia,j)%r
	        end if 
        end do 
	end do
    close(94632)
    end if 

!		if (dat(ia,j)%expr>-99.and.dat(ia,j)%expr<0) then 
	!			print*, 'dat exp is negative',ia,j,dat(ia,j)%expr,dat(ia,j)%hme,dat(ia,j)%nn,dat(ia,j)%mm,dat(ia,j)%r
	!			stop
	!		end if
    
	!if ( iter==1) then
	!do ia=mna,mxa
	!	do j=1,ndat
            !The below check does give you an error. Despite the fact that everyone starts from age 18 and there is no missing cohort variable in the data (all is 1920-49).
            !Why?
            !This is because, for those people who are ED, the agestart is 22. So in read_actualdata, after reading these people's initcond's from their age 18 (and for all of them I have all those initconds), 
            !I set all the variables from 18 to 21 to -99, including cohort (after reading cohort info into init). Hence for such people, you
            !will have a cohort variable value of -99 at age 18-21. But starting from 22, you will have the cohort value. 
            !For some, there are some ages for which there is no observation in the PSID. For example, for idnum who is ED, the person is observed from age 18-20
            !but then there is nothing until 24. So for this person, the cohort values (as well as other variable values) would be all filled in with -99 until 24. 
            !SO instead of conditioning the nummove moments on dat(mna,:)%co=co, I condition them on coho(mna,:).			
            !(I generate coho with		coho(mna:mxai,:)= (  maxval(dat(mna:mxai,:)%co)==co  )			
            !This solves the problem. 
            !if (dat(ia,j)%co<1.or.dat(ia,j)%co>nco ) then
                !print*, 'something wrong with dat%co A',ia,j,dat(ia,j)%co,ndat
                !stop
            !end if
            !if ( maxval(dat(:,j)%co) /= dat(ia,j)%co  ) then
                !print*, 'something wrong with dat%co B',ia,j,dat(ia,j)%co,ndat
                !stop
            !end if
    !    end do 
	!end do
	!end if 
    
    !if (onlysingles) then ! if doing only singles, condition all mom on singles 
	!	nomiss(mna:mxai-1,:) = (  nomiss(mna:mxai-1,:) .and. dat(mna:mxai-1,:)%rel==0 .and. dat(mna+1:mxai,:)%rel==0 ) 
	!end if 

    if (onlysingles) then 
        maxrelo=0
    else 
        maxrelo=1
    end if 
    if (onlymales.and.(.not.onlyfem)) then 
        minsex=1 ; maxsex=1
    else if ( (.not.onlymales).and.onlyfem) then 
        minsex=2 ; maxsex=2
    else if ( (.not.onlymales).and.(.not.onlyfem)) then 
        minsex=1 ; maxsex=2
    else 
        print*, 'Something is wrong in the neighborhood'
        stop
    end if 
    
    !if (onlymales) then 
    !    maxsex=1
    !else 
    !    maxsex=2
    !end if 
    
    !ag092922 sept2022: note that I got rid of all the if stateements following weights(im) i.e. if(onlysingles and j==1) weights(im)=0.0_dp
    !got rid of them because I don't need them if I have maxrel0 since maxrel0 is set to 1 if onlysingles in the above if statement before the co do loop starts. 
    !this if (onlysingles and j==1) was leading to bugs i.e. different processors would have different momwgt (all else same tho at least)
    !because I wasn't assigning j in the loops for just co or cosex (j is for rel) and when I moved around the moments do loops I forgot that j there 
    !so it was being assigned something different for each processor and they were each gibing different momwgt because of that. 
    !but just get rid of onlysingles if statements after weights(im) since you don't need that anyway (due tot he presence of maxrel)
	cohort: do co=1,nco
		coho(MNAD:MXA,:)= (  dat(MNAD:MXA,:)%co==co  )		
        corel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co  .and. dat(MNAD:MXAD,:)%rel==1 .and. norelchg(MNAD:MXAD,:)==1 )
		!cohogen(:)=  maxval(dat(mna:mxai,:)%co)  		            

    
        call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 ), d1*one(dat(MNA:MXA,:)%rel==1),mom,cnt,var)
        write(name(im),'("mar ",tr10)')	
        weights(im)=wrel !; if (onlysingles) weights(im)=0.0_dp   !forget about that if statement usually. just putting there to get im to go up by 1 for the next loop
        im=im+1

        headloc(ihead)=im; headstr(ihead)='Kids';ihead=ihead+1
        CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%rel==1).and.(norelchg(MNA:MXAD,:)==1).AND.(kidtrans(MNA:MXAD,:)>=0).AND.(iacat(MNA:MXAD,:)==1)),d1*kidtrans(MNA:MXAD,:),mom,cnt,var)
        WRITE(name(im),'("kidtrans-married ")') 
        weights(im)=0.0_dp !wkid
        im=im+1
        CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%rel==1).and.(norelchg(MNA:MXAD,:)==1).AND.(dat(MNA:MXAD,:)%kidr>=0).AND.(iacat(MNA:MXAD,:)==1)),d1*dat(MNA:MXAD,:)%kidr,mom,cnt,var)
        WRITE(name(im),'("kids-married     ")') 
        weights(im)= 0.0_dp !wkid 
        im=im+1

       !loop by only sex
        do g=minsex,maxsex
            cosex(MNAD:MXA,:)= (dat(MNAD:MXA,:)%co==co .and. dat(MNAD:MXA,:)%sexr==g  )

            headloc(ihead)=im
            if (g==1) headstr(ihead)='all men move'
            if (g==2) headstr(ihead)='all fem move'
            ihead=ihead+1

            !ahu jan19 010219: so I am still writing the age 17 (mad) mar rate. age 17 mar rate does not exit in the data since we don't record age 17 rel there. 
            !so the weight on this moment is 0 but then it's wmovebyrel
            ia=MNAD
            CALL condmom(im,( cosex(ia,:) .AND. move(ia,:)>=0 ),d1*move(ia,:),mom,cnt,var)
            WRITE(name(im),'("move by age ",I4)') ia
            weights(im)=wmove
            im=im+1
            do ia=MNA,MXAD,5
                CALL condmom(im,( cosex(ia,:) .AND. move(ia,:)>=0 ),d1*move(ia,:),mom,cnt,var)
                WRITE(name(im),'("move by age",tr3,I4)') ia
                weights(im)=wmove
                im=im+1
            end do     

            call condmom(im,(  cohogen(:)==co .and. sexgen(:)==g ) ,   d1* one( nummove(:)==0 ) ,mom,cnt,var)	
            write(name(im),'("nummove=0 ",tr5)')  
            weights(im)=wmove
            im=im+1
            call condmom(im,(  cohogen(:)==co  .and. sexgen(:)==g ),   d1* one( nummove(:)==1 ) ,mom,cnt,var)	
            write(name(im),'("nummove=1 ",tr5)')  
            weights(im)=wmove
            im=im+1
            call condmom(im,(  cohogen(:)==co  .and. sexgen(:)==g ),   d1* one( nummove(:)>=2 ) ,mom,cnt,var)	
            write(name(im),'("nummove>=2 ",tr5)')  
            weights(im)=wmove
            im=im+1


            headloc(ihead)=im
            if (g==1) headstr(ihead)='all men w and wvar at 18 and 18:19 .. for ident'
            if (g==2) headstr(ihead)='all fem w and wvar at 18 and 18:19 .. for ident'
            ihead=ihead+1

            CALL condmom(im,(   cosex(18,:) .AND.  dat(18,:)%hhr==1 .AND. dat(18,:)%edr==1 .AND. dat(18,:)%logwr>=0) ,d1*dat(18,:)%logwr,mom,cnt,var)
            WRITE(name(im),'("wnned|18 small samp size!")') 
            weights(im)=0.0_dp !; if (onlysingles.and.j==1) weights(im)=0.0_dp   !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess    
            calcvar(im)=1
            im=im+1
            CALL condmom(im,(   cosex(18,:) .AND.  dat(18,:)%hhr==1 .AND. dat(18,:)%edr==1 .AND. dat(18,:)%logwr>=0) ,d1*  (dat(18,:)%logwr**2),mom,cnt,var)
            WRITE(name(im),'("wvarned|18 small samp size!")') 
            weights(im)=0.0_dp !; if (onlysingles.and.j==1) weights(im)=0.0_dp    !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess         
            calcvar(im)=5
            im=im+1
            CALL condmom(im,(   cosex(18:19,:) .AND.  dat(18:19,:)%hhr==1 .AND. dat(18:19,:)%edr==1 .AND. dat(18:19,:)%logwr>=0) ,d1*dat(18:19,:)%logwr,mom,cnt,var)
            WRITE(name(im),'("wnned|1819 ")') 
            weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp   !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess    
            calcvar(im)=1
            im=im+1
            CALL condmom(im,(   cosex(18:19,:) .AND.  dat(18:19,:)%hhr==1 .AND. dat(18:19,:)%edr==1 .AND. dat(18:19,:)%logwr>=0) ,d1*  (dat(18:19,:)%logwr**2),mom,cnt,var)
            WRITE(name(im),'("wvarned|1819 ")') 
            weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp    !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess         
            calcvar(im)=5
            im=im+1            
            do i=1,nl,8
                CALL condmom(im,(   cosex(18,:) .AND.  dat(18,:)%hhr==1 .AND. dat(18,:)%edr==1 .AND. dat(18,:)%logwr>=0 .AND. dat(18,:)%l==i) ,d1*dat(18,:)%logwr,mom,cnt,var)
                WRITE(name(im),'("wned|18 by loc small!",i4)') i
                weights(im)=0.0_dp !; if (onlysingles.and.j==1) weights(im)=0.0_dp    !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess         
                calcvar(im)=1
                im=im+1
                CALL condmom(im,(   cosex(18,:) .AND.  dat(18,:)%hhr==1 .AND. dat(18,:)%edr==1 .AND. dat(18,:)%logwr>=0  .AND. dat(18,:)%l==i) ,d1*  (dat(18,:)%logwr**2),mom,cnt,var)
                WRITE(name(im),'("wvarned|18 by loc small!",i4)') i
                weights(im)=0.0_dp !; if (onlysingles.and.j==1) weights(im)=0.0_dp  !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess           
                calcvar(im)=5
                im=im+1
                CALL condmom(im,(   cosex(18:19,:) .AND.  dat(18:19,:)%hhr==1 .AND. dat(18:19,:)%edr==1 .AND. dat(18:19,:)%logwr>=0 .AND. dat(18:19,:)%l==i) ,d1*dat(18:19,:)%logwr,mom,cnt,var)
                WRITE(name(im),'("wned|18:19 by loc",i4)') i
                weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp    !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess         
                calcvar(im)=1
                im=im+1
                CALL condmom(im,(   cosex(18:19,:) .AND.  dat(18:19,:)%hhr==1 .AND. dat(18:19,:)%edr==1 .AND. dat(18:19,:)%logwr>=0  .AND. dat(18:19,:)%l==i) ,d1*  (dat(18:19,:)%logwr**2),mom,cnt,var)
                WRITE(name(im),'("wvarned|18:19 by loc",i4)') i
                weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp  !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess           
                calcvar(im)=5
                im=im+1
                CALL condmom(im,(   cosex(20:25,:) .AND.  dat(20:25,:)%hhr==1 .AND. dat(20:25,:)%edr==1 .AND. dat(20:25,:)%logwr>=0 .AND. dat(20:25,:)%l==i) ,d1*dat(20:25,:)%logwr,mom,cnt,var)
                WRITE(name(im),'("wned|20:25 by loc",i4)') i
                weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp    !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess         
                calcvar(im)=1
                im=im+1
                CALL condmom(im,(   cosex(20:25,:) .AND.  dat(20:25,:)%hhr==1 .AND. dat(20:25,:)%edr==1 .AND. dat(20:25,:)%logwr>=0  .AND. dat(20:25,:)%l==i) ,d1*  (dat(20:25,:)%logwr**2),mom,cnt,var)
                WRITE(name(im),'("wvarned|20:25 by loc",i4)') i
                weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp  !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess           
                calcvar(im)=5
                im=im+1
            end do 
            call condmom(im,( cosex(MNA:MXAD,:) .AND. dee(MNA:MXAD,:)==1  .AND. move(MNA:MXAD,:)==0 ),   d1*( dat(MNA+1:MXA,:)%logwr-dat(MNA:MXAD,:)%logwr ),mom,cnt,var)		
            write(name(im),'("wdif | stay ",tr2)')  
            weights(im)=wdifww  
            !calcvar(im)=1
            im=im+1 
            call condmom(im,( cosex(MNA:MXAD,:) .AND. dee(MNA:MXAD,:)==1  .AND. move(MNA:MXAD,:)==1 ),   d1*( dat(MNA+1:MXA,:)%logwr-dat(MNA:MXAD,:)%logwr ),mom,cnt,var)		
            write(name(im),'("wdif | move ",tr2)')  
            weights(im)=wdifww 
            !calcvar(im)=1
            im=im+1             
            call condmom(im,( cosex(MNA:MXAD-1,:) .AND. deue(MNA:MXAD-1,:)==1 .AND. dat(MNA+1:MXAD,:)%l==dat(MNA:MXAD-1,:)%l .AND. dat(MNA+2:MXA,:)%l==dat(MNA:MXAD-1,:)%l),   d1*( dat(MNA+2:MXA,:)%logwr-dat(MNA:MXAD-1,:)%logwr ),mom,cnt,var)		
            write(name(im),'("wdif | eue,s ",tr2)')  
            weights(im)=wdifww
            !calcvar(im)=1
            im=im+1 
            call condmom(im,( cosex(MNA:MXAD-1,:) .AND. deue(MNA:MXAD-1,:)==1   .AND. dat(MNA+1:MXAD,:)%l==dat(MNA:MXAD-1,:)%l .AND. dat(MNA+2:MXA,:)%l/=dat(MNA:MXAD-1,:)%l .AND. dat(MNA+2:MXA,:)%l>0 ),   d1*( dat(MNA+2:MXA,:)%logwr-dat(MNA:MXAD-1,:)%logwr ),mom,cnt,var)		
            write(name(im),'("wdif | eue,m ",tr2)')  
            weights(im)=wdifww 
            !calcvar(im)=1
            im=im+1 


            do i=1,nl
                CALL condmom(im,(   cosex(18:19,:) .AND.  dat(18:19,:)%hhr==1 .AND. dat(18:19,:)%edr==1 .AND. dat(18:19,:)%logwr>=0 .AND. dat(18:19,:)%l==i) ,d1*dat(18:19,:)%logwr,mom,cnt,var)
                WRITE(name(im),'("wned|18:19 by loc",i4)') i
                weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp    !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess         
                !calcvar(im)=1
                im=im+1
                !CALL condmom(im,(   cosex(18:19,:) .AND.  dat(18:19,:)%hhr==1 .AND. dat(18:19,:)%edr==1 .AND. dat(18:19,:)%logwr>=0  .AND. dat(18:19,:)%l==i) ,d1*  (dat(18:19,:)%logwr**2),mom,cnt,var)
                !WRITE(name(im),'("wvarned|18:19 by loc",i4)') i
                !weights(im)=wwage !; if (onlysingles.and.j==1) weights(im)=0.0_dp  !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess           
                !calcvar(im)=5
                !im=im+1
            end do 

                            
            headloc(ihead)=im
            if (g==1) headstr(ihead)='all men wdecile for ident'
            if (g==2) headstr(ihead)='all fem wdecile for ident'
            ihead=ihead+1

            do ddd=1,ndecile-1
                CALL condmom(im,(   cosex(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%edr==1 .AND. dat(mna:mxa,:)%logwr>=0 ) ,d1*one( (dat(mna:mxa,:)%logwr>decilegrid(ddd) .and. dat(mna:mxa,:)%logwr<decilegrid(ddd+1) ) ),mom,cnt,var)
                WRITE(name(im),'("wdecile|ned",tr1,i4)') ddd
                weights(im)=wwaged !; if (onlysingles.and.j==1) weights(im)=0.0_dp !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess 
                im=im+1
            end do
            do ddd=1,ndecile-1
                CALL condmom(im,(   cosex(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%edr==2 .AND. dat(mna:mxa,:)%logwr>=0 ) ,d1*one( (dat(mna:mxa,:)%logwr>decilegrid(ddd) .and. dat(mna:mxa,:)%logwr<decilegrid(ddd+1) ) ),mom,cnt,var)
                WRITE(name(im),'("wdecile| ed",tr1,i4)') ddd
                weights(im)=wwaged !; if (onlysingles.and.j==1) weights(im)=0.0_dp  !ag092922 sept2022 after I moved around these moms, there was still j left here and that made different procesors have different objval because momwgts were different since j was just assigned a different value by each processors I guess       
                im=im+1
            end do

            

            !ia=18
            !do i=1,nl
            !    CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==1 .AND. dat(ia,:)%logwr>=0 .AND. dat(ia,:)%l==i) ,d1*dat(ia,:)%logwr,mom,cnt,var)
            !    WRITE(name(im),'("w|noed l 18",tr1,i4)') i
            !    weights(im)=0.0_dp ; if (onlysingles.and.j==1) weights(im)=0.0_dp            
            !    im=im+1
            !end do 
            
            !ia=18
            !do i=1,nl
            !    CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==2 .AND. dat(ia,:)%logwr>=0 .AND. dat(ia,:)%l==i) ,d1*dat(ia,:)%logwr,mom,cnt,var)
            !    WRITE(name(im),'("w|ed l  18",tr1,i4)') i
            !    weights(im)=0.0_dp ; if (onlysingles.and.j==1) weights(im)=0.0_dp            
            !    im=im+1
            !end do 
        end do !g for sex
        !******************************************

        !***************************************
        !START OF LOOP BY SEX AND REL
        do g=minsex,maxsex   
        do j=maxrelo,0,-1
            if ( onlysingles ) then !  (.not.onlysingles).or.(onlysingles.and.j==0) ) then 
                cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g )			                    
            else 
                cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
            end if 
            
            if (j==1) wmovebyrel=wmovemar
            if (j==0) wmovebyrel=wmovesin


            headloc(ihead)=im
            if (g==1.and.j==0) headstr(ihead)='single men move rates'
            if (g==1.and.j==1) headstr(ihead)='married men move rates'
            if (g==2.and.j==0) headstr(ihead)='single fem move rates'
            if (g==2.and.j==1) headstr(ihead)='married fem move rates'
            ihead=ihead+1
            !ahu jan19 010219: so I am still writing the age 17 (mad) mar rate. age 17 mar rate does not exit in the data since we don't record age 17 rel there. 
            !so the weight on this moment is 0 but then it's wmovebyrel
            ia=MNAD
            CALL condmom(im,( cosexrel(ia,:) .AND. move(ia,:)>=0 ),d1*move(ia,:),mom,cnt,var)
            WRITE(name(im),'("move by age ",I4)') ia
            weights(im)=wmovebyrel
            im=im+1
            do ia=MNA,MXAD,5
                CALL condmom(im,( cosexrel(ia,:) .AND. move(ia,:)>=0 ),d1*move(ia,:),mom,cnt,var)
                WRITE(name(im),'("move by age",tr3,I4)') ia
                weights(im)=wmovebyrel 
                im=im+1
            end do     
            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0  .AND. move(MNA:MXAD,:)>=0 ),   d1* move(MNA:MXAD,:) ,mom,cnt,var)	
            write(name(im),'("move | u ",tr5)')  
            weights(im)=wmovebyrel 
            im=im+1
            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1  .AND. move(MNA:MXAD,:)>=0 ),   d1* move(MNA:MXAD,:) ,mom,cnt,var)	
            write(name(im),'("move | e ",tr5)')  
            weights(im)=wmovebyrel
            im=im+1  
            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. move(MNA:MXAD,:)==1 ),   d1*one( dat(MNA+1:MXA,:)%l==dat(MNA+1:MXA,:)%hme  ),mom,cnt,var)		
            write(name(im),'("%hme-mvs",tr3)')  !added this ahu 121718
            weights(im)=wmovebyrel
            im=im+1 


            headloc(ihead)=im
            if (g==1.and.j==0) headstr(ihead)='single men wages'
            if (g==1.and.j==1) headstr(ihead)='married men wages'
            if (g==2.and.j==0) headstr(ihead)='single fem wages'
            if (g==2.and.j==1) headstr(ihead)='married fem wages'
            ihead=ihead+1

            CALL condmom(im,(   cosexrel(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%edr==1 .AND. dat(mna:mxa,:)%logwr>=0 ) ,d1*dat(mna:mxa,:)%logwr,mom,cnt,var)
            WRITE(name(im),'("w|noed",tr1)') 
            weights(im)=wwage 
            calcvar(im)=1
            im=im+1
            CALL condmom(im,(   cosexrel(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%edr==1 .AND. dat(mna:mxa,:)%logwr>=0 ) ,d1*  (dat(mna:mxa,:)%logwr**2),mom,cnt,var)
            WRITE(name(im),'("wvar|noed",tr1)') 
            weights(im)=wwage           
            calcvar(im)=5
            im=im+1
            !CALL condmom(im,(   cosexrel(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%edr==1 .AND. dat(mna:mxa,:)%logwr>=0 ) ,d1*  (dat(mna:mxa,:)%logwr-mom(im-2))**2,mom,cnt,var)
            !WRITE(name(im),'("wrng|noed",tr1)') 
            !weights(im)=0.0_dp
            !im=im+1            
            CALL condmom(im,(   cosexrel(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%edr==2 .AND. dat(mna:mxa,:)%logwr>=0 ) ,d1*dat(mna:mxa,:)%logwr,mom,cnt,var)
            WRITE(name(im),'("w|  ed",tr1)') 
            weights(im)=wwage 
            calcvar(im)=1
            im=im+1
            CALL condmom(im,(   cosexrel(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%edr==2 .AND. dat(mna:mxa,:)%logwr>=0 ) ,d1*  (dat(mna:mxa,:)%logwr**2),mom,cnt,var)
            WRITE(name(im),'("wvar|ed",tr1)') 
            weights(im)=wwage 
            calcvar(im)=5
            im=im+1

            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dee(MNA:MXAD,:)==1  .AND. move(MNA:MXAD,:)==1 .and. dat(MNA+1:MXA,:)%l/=dat(MNA+1:MXA,:)%hme ),   d1*( dat(MNA+1:MXA,:)%logwr-dat(MNA:MXAD,:)%logwr ),mom,cnt,var)		
            write(name(im),'("wdif | hmemve=0 ",tr2)')  
            weights(im)=wdifww
            im=im+1 
            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dee(MNA:MXAD,:)==1  .AND. move(MNA:MXAD,:)==1 .and. dat(MNA+1:MXA,:)%l==dat(MNA+1:MXA,:)%hme ),   d1*( dat(MNA+1:MXA,:)%logwr-dat(MNA:MXAD,:)%logwr ),mom,cnt,var)		
            write(name(im),'("wdif | hmemve=1 ",tr2)')  
            weights(im)=wdifww
            im=im+1 

            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dee(MNA:MXAD,:)==1  .AND. move(MNA:MXAD,:)==0 ),   d1*( dat(MNA+1:MXA,:)%logwr-dat(MNA:MXAD,:)%logwr ),mom,cnt,var)		
            write(name(im),'("wdif | stay ",tr2)')  
            weights(im)=wdifww  
            calcvar(im)=1
            im=im+1 
            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dee(MNA:MXAD,:)==1 .AND. move(MNA:MXAD,:)==0 ),   d1*( dat(MNA+1:MXA,:)%logwr-dat(MNA:MXAD,:)%logwr )**2,mom,cnt,var)		
            write(name(im),'("wdif2 | stay ",tr2)')  
            weights(im)=wdifww 
            calcvar(im)=5
            im=im+1 
            !call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dee(MNA:MXAD,:)==1 .AND. move(MNA:MXAD,:)==0 ),   d1*( dat(MNA+1:MXA,:)%logwr-dat(MNA:MXAD,:)%logwr - mom(im-2) )**2,mom,cnt,var)		
            !write(name(im),'("wdif2p | stay ",tr2)')  
            !weights(im)=0.0_dp 
            !im=im+1 
            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dee(MNA:MXAD,:)==1  .AND. move(MNA:MXAD,:)==1 ),   d1*( dat(MNA+1:MXA,:)%logwr-dat(MNA:MXAD,:)%logwr ),mom,cnt,var)		
            write(name(im),'("wdif | move ",tr2)')  
            weights(im)=wdifww 
            calcvar(im)=1
            im=im+1 
            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dee(MNA:MXAD,:)==1 .AND. move(MNA:MXAD,:)==1 ),   d1*( dat(MNA+1:MXA,:)%logwr-dat(MNA:MXAD,:)%logwr )**2,mom,cnt,var)		
            write(name(im),'("wdif2 | move ",tr2)')  
            weights(im)=wdifww  
            calcvar(im)=5
            im=im+1 
            !call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dee(MNA:MXAD,:)==1 .AND. move(MNA:MXAD,:)==1 ),   d1*( dat(MNA+1:MXA,:)%logwr-dat(MNA:MXAD,:)%logwr - mom(im-2) )**2,mom,cnt,var)		
            !write(name(im),'("wdif2p | move ",tr2)')  
            !weights(im)=100.0_dp 
            !im=im+1 

            
            call condmom(im,( cosexrel(MNA:MXAD-1,:) .AND. deue(MNA:MXAD-1,:)==1 .AND. dat(MNA+1:MXAD,:)%l==dat(MNA:MXAD-1,:)%l .AND. dat(MNA+2:MXA,:)%l==dat(MNA:MXAD-1,:)%l),   d1*( dat(MNA+2:MXA,:)%logwr-dat(MNA:MXAD-1,:)%logwr ),mom,cnt,var)		
            write(name(im),'("wdif | eue,s ",tr2)')  
            weights(im)=wdifww
            calcvar(im)=1
            im=im+1 
            call condmom(im,( cosexrel(MNA:MXAD-1,:) .AND. deue(MNA:MXAD-1,:)==1  .AND. dat(MNA+1:MXAD,:)%l==dat(MNA:MXAD-1,:)%l .AND. dat(MNA+2:MXA,:)%l==dat(MNA:MXAD-1,:)%l),   d1*( dat(MNA+2:MXA,:)%logwr-dat(MNA:MXAD-1,:)%logwr )**2,mom,cnt,var)		
            write(name(im),'("wdif2 | eue,s ",tr2)')  
            weights(im)=wdifww 
            calcvar(im)=5
            im=im+1 
            !call condmom(im,( cosexrel(MNA:MXAD-1,:) .AND. deue(MNA:MXAD-1,:)==1  .AND. dat(MNA+1:MXAD,:)%l==dat(MNA:MXAD-1,:)%l .AND. dat(MNA+2:MXA,:)%l==dat(MNA:MXAD-1,:)%l),   d1*( dat(MNA+2:MXA,:)%logwr-dat(MNA:MXAD-1,:)%logwr - mom(im-2))**2,mom,cnt,var)		
            !write(name(im),'("wdif2p | eue,s ",tr2)')  
            !weights(im)=100.0_dp  ; if (onlysingles.and.j==1) weights(im)=0.0_dp
            !im=im+1 


            call condmom(im,( cosexrel(MNA:MXAD-1,:) .AND. deue(MNA:MXAD-1,:)==1   .AND. dat(MNA+1:MXAD,:)%l==dat(MNA:MXAD-1,:)%l .AND. dat(MNA+2:MXA,:)%l/=dat(MNA:MXAD-1,:)%l .AND. dat(MNA+2:MXA,:)%l>0 ),   d1*( dat(MNA+2:MXA,:)%logwr-dat(MNA:MXAD-1,:)%logwr ),mom,cnt,var)		
            write(name(im),'("wdif | eue,m ",tr2)')  
            weights(im)=wdifww 
            calcvar(im)=1
            im=im+1 
            call condmom(im,( cosexrel(MNA:MXAD-1,:) .AND. deue(MNA:MXAD-1,:)==1   .AND. dat(MNA+1:MXAD,:)%l==dat(MNA:MXAD-1,:)%l .AND. dat(MNA+2:MXA,:)%l/=dat(MNA:MXAD-1,:)%l .AND. dat(MNA+2:MXA,:)%l>0 ),   d1*( dat(MNA+2:MXA,:)%logwr-dat(MNA:MXAD-1,:)%logwr )**2,mom,cnt,var)		
            write(name(im),'("wdif2 | eue,m ",tr2)')  
            weights(im)=wdifww
            calcvar(im)=5
            im=im+1 
            !call condmom(im,( cosexrel(MNA:MXAD-1,:) .AND. deue(MNA:MXAD-1,:)==1   .AND. dat(MNA+1:MXAD,:)%l==dat(MNA:MXAD-1,:)%l .AND. dat(MNA+2:MXA,:)%l/=dat(MNA:MXAD-1,:)%l .AND. dat(MNA+2:MXA,:)%l>0 ),   d1*( dat(MNA+2:MXA,:)%logwr-dat(MNA:MXAD-1,:)%logwr -mom(im-2))**2,mom,cnt,var)		
            !write(name(im),'("wdif2p | eue,m ",tr2)')  
            !weights(im)=100.0_dp 
            !im=im+1 

            do ia=mna,mxad,10
                call condmom(im,( cosexrel(ia,:) .AND. dee(ia,:)==1  .AND. move(ia,:)==0 ),   d1*( dat(ia+1,:)%logwr-dat(ia,:)%logwr ),mom,cnt,var)		
                write(name(im),'("wdif | stay ia ",tr2,i6)')  ia
                weights(im)=wdifww
                im=im+1 
          end do   
          do ia=mna,mxad,10
                call condmom(im,( cosexrel(ia,:) .AND. dee(ia,:)==1  .AND. move(ia,:)==1 ),   d1*( dat(ia+1,:)%logwr-dat(ia,:)%logwr ),mom,cnt,var)		
                write(name(im),'("wdif | move ia ",tr2,i6)')  ia
                weights(im)=wdifww
                im=im+1 
          end do   


          headloc(ihead)=im
          if (g==1.and.j==0) headstr(ihead)='single men emp transitions'
          if (g==1.and.j==1) headstr(ihead)='married men emp transitionss'
          if (g==2.and.j==0) headstr(ihead)='single fem emp transitions'
          if (g==2.and.j==1) headstr(ihead)='married fem emp transitions'
          ihead=ihead+1

            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0 ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
            write(name(im),'("e | u stay",tr3)')  
            weights(im)=wtrans 
            im=im+1 

            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
            write(name(im),'("e | u move",tr3)')  
            weights(im)=wtrans 
            im=im+1 

            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0 ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
            write(name(im),'("e | e stay",tr3)')  
            weights(im)=wtrans 
            im=im+1 
            
            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1 ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
            write(name(im),'("e | e move",tr3)')  
            weights(im)=wtrans 
            im=im+1 

            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0  ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. move(MNA:MXAD,:)==0),mom,cnt,var)		
            write(name(im),'("e-stay | e ",tr3)')  
            weights(im)=wtrans 
            im=im+1 
            
            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0  ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. move(MNA:MXAD,:)==0),mom,cnt,var)		
            write(name(im),'("u-stay | e ",tr3)')  
            weights(im)=wtrans  
            im=im+1 
            
            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0  ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. move(MNA:MXAD,:)==1),mom,cnt,var)		
            write(name(im),'("e-move | e ",tr3)')  
            weights(im)=wtrans  
            im=im+1 
            
            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0  ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. move(MNA:MXAD,:)==1),mom,cnt,var)		
            write(name(im),'("u-move | e ",tr3)')  
            weights(im)=wtrans  
            im=im+1 

                call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0  ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. move(MNA:MXAD,:)==0),mom,cnt,var)		
            write(name(im),'("e-stay | u ",tr3)')  
            weights(im)=wtrans
            im=im+1 
            
            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0  ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. move(MNA:MXAD,:)==0),mom,cnt,var)		
            write(name(im),'("u-stay | u ",tr3)')  
            weights(im)=wtrans  
            im=im+1 
            
            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0  ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 .AND. move(MNA:MXAD,:)==1),mom,cnt,var)		
            write(name(im),'("e-move | u ",tr3)')  
            weights(im)=wtrans  
            im=im+1 
            
            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0  ),   d1*one( dat(MNA+1:MXA,:)%hhr==0 .AND. move(MNA:MXAD,:)==1),mom,cnt,var)		
            write(name(im),'("u-move | u ",tr3)')  
            weights(im)=wtrans  
            im=im+1 

 
            
            !do i=1,nl
            !    CALL condmom(im,(   cosexrel(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%edr==1 .AND. dat(mna:mxa,:)%logwr>=0 .AND. dat(mna:mxa,:)%l==i) ,d1*dat(mna:mxa,:)%logwr,mom,cnt,var)
            !    WRITE(name(im),'("w|noed l",tr1,i4)') i
            !    weights(im)=wwage    
            !    calcvar(im)=1
            !    im=im+1
            !    CALL condmom(im,(   cosexrel(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%edr==1 .AND. dat(mna:mxa,:)%logwr>=0  .AND. dat(mna:mxa,:)%l==i) ,d1*  (dat(mna:mxa,:)%logwr**2),mom,cnt,var)
            !    WRITE(name(im),'("wvar|noed l",tr1,i4)') i
            !    weights(im)=wwage         
            !    calcvar(im)=5
            !    im=im+1
            !end do 
            

            !do i=1,nl
            !    CALL condmom(im,(   cosexrel(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%edr==2 .AND. dat(mna:mxa,:)%logwr>=0 .AND. dat(mna:mxa,:)%l==i) ,d1*dat(mna:mxa,:)%logwr,mom,cnt,var)
            !    WRITE(name(im),'("w|ed l",tr1,i4)') i
            !    weights(im)=wwage           
            !    calcvar(im)=1
            !    im=im+1
            !    CALL condmom(im,(   cosexrel(mna:mxa,:) .AND.  dat(mna:mxa,:)%hhr==1 .AND. dat(mna:mxa,:)%edr==2 .AND. dat(mna:mxa,:)%logwr>=0  .AND. dat(mna:mxa,:)%l==i) ,d1*  (dat(mna:mxa,:)%logwr**2),mom,cnt,var)
            !    WRITE(name(im),'("wvar|ed l",tr1,i4)') i
            !    weights(im)=wwage     
            !    calcvar(im)=5
            !    im=im+1
            !end do 
            !do i=1,nl
            !    CALL condmom(im,(   cosexrel(18:21,:) .AND.  dat(18:21,:)%hhr==1 .AND. dat(18:21,:)%edr==1 .AND. dat(18:21,:)%logwr>=0  .AND. dat(18:21,:)%l==i) ,d1*  (dat(18:21,:)%logwr),mom,cnt,var)
            !    WRITE(name(im),'("w|1821 l",tr1,i4)') i
            !    weights(im)=wwage        
            !    im=im+1
            !    !CALL condmom(im,(   cosexrel(18:21,:) .AND.  dat(18:21,:)%hhr==1 .AND. dat(18:21,:)%edr==1 .AND. dat(18:21,:)%logwr>=0  .AND. dat(18:21,:)%l==i) ,d1*  (dat(18:21,:)%logwr**2-mom(im-1)**2),mom,cnt,var)
            !    !WRITE(name(im),'("w2|1821 l",tr1,i4)') i
            !    !weights(im)=wwage             
            !    !im=im+1
            !    !CALL condmom(im,(   cosexrel(18:21,:) .AND.  dat(18:21,:)%hhr==1 .AND. dat(18:21,:)%edr==1 .AND. dat(18:21,:)%logwr>=0  .AND. dat(18:21,:)%l==i) ,d1*  (dat(18:21,:)%logwr**3-mom(im-1)**3),mom,cnt,var)
            !    !WRITE(name(im),'("w3|1821 l",tr1,i4)') i
            !    !weights(im)=wwage     
            !    !im=im+1
            !end do 
            
            call condmom(im,( cosexrel(mna:mxa,:) .AND. dat(mna:mxa,:)%hhr>=0 ), d1*one( dat(mna:mxa,:)%hhr==1 ),mom,cnt,var)		
            write(name(im),'("e ",tr13)')			
            weights(im)=whour 
            im=im+1 
            
            call condmom(im,( cosexrel(mna:mxa,:) .AND. dat(mna:mxa,:)%hhr>=0 .AND. dat(mna:mxa,:)%edr==1 ),   d1*one( dat(mna:mxa,:)%hhr==1 ),mom,cnt,var)	
            write(name(im),'("e | noed",tr6)')  
            weights(im)=whour  
            im=im+1 
            
            call condmom(im,( cosexrel(mna:mxa,:) .AND. dat(mna:mxa,:)%hhr>=0 .AND. dat(mna:mxa,:)%edr==2 ),   d1*one( dat(mna:mxa,:)%hhr==1 ),mom,cnt,var)		
            write(name(im),'("e |   ed",tr6)')  
            weights(im)=whour 
            im=im+1 

            !Note about conditioning on age:
            !The max endage in data is 47. 
            !The way sim is done, for those whose endage is 47, the last age where variables get recorded is ia-1=46. 
            !This is because dat(ia-1,.) is recorded for each ia. So for the last age 47, the variables for 46 gets written
            !But then there is nothing after age 46, despite the fact that we do have people whose endage is 46 (namely 47).      
            !If one of the conditioning statements is norelchg, then there is nothing after age 45. 
            !Since any age 46 no relchg would need a age 47 rel. 
            !This is why emp and wage by age cond on cosexrel only goes until 45. 
            do ia=MNAd,25 
                CALL condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr>=0 ) ,d1*one(dat(ia,:)%hhr==1),mom,cnt,var)
                WRITE(name(im),'("emp by age",tr4,I2)') ia
                weights(im)=whour 
                if (ia<mna) weights(im)=0.0_dp
                im=im+1
            end do       
            do ia=agestart(NOCOLLEGE),mxad,10
                CALL condmom(im,( cosexrel(ia,:) .AND.    dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==1 .AND. dat(ia,:)%logwr>=0 ) ,d1*dat(ia,:)%logwr,mom,cnt,var)
                WRITE(name(im),'("w|noed by age",tr1,I2)') ia
                weights(im)=wwage 
                im=im+1
                !CALL condmom(im,( cosexrel(ia,:) .AND.    dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==1 .AND. dat(ia,:)%logwr>=0 ) ,d1*  (dat(ia,:)%logwr-mom(im-1))**2   ,mom,cnt,var)
                !WRITE(name(im),'("wvar|noed by age",tr1,I2)') ia
                !weights(im)=wwage
                !im=im+1
            end do            
            
            do ia=agestart(COLLEGE),mxad,10
                CALL condmom(im,( cosexrel(ia,:) .AND.   dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==2 .AND. dat(ia,:)%logwr>=0 ) ,d1*dat(ia,:)%logwr,mom,cnt,var)
                WRITE(name(im),'("w|  ed by age",tr1,I2)') ia
                weights(im)=wwage 
                im=im+1
                !CALL condmom(im,( cosexrel(ia,:) .AND.    dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==2 .AND. dat(ia,:)%logwr>=0 ) ,d1*  (dat(ia,:)%logwr-mom(im-1))**2   ,mom,cnt,var)
                !WRITE(name(im),'("wvar|ed by age",tr1,I2)') ia
                !weights(im)=wwage 
                !im=im+1
            end do            
        
            !do i=1,nl
            !    call condmom(im,( cosexrel(MNA:MXA,:) .AND. dat(MNA:MXA,:)%hhr==1 .AND. dat(MNA:MXA,:)%l==i .AND. dat(MNA:MXA,:)%logwr>=0 ),   d1*dat(MNA:MXA,:)%logwr ,mom,cnt,var)		
            !    write(name(im),'("w|loc",tr9,i4)') i			
            !    weights(im)=wwage 
            !    im=im+1 
            !end do 
        
            if (j==1) then !only do kid by age for married. for singles, it looks weird in the data (.30, .8, .16,.19) 
                do ia=MNA,MXAD,8
                    CALL condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%kidr>=1 ),d1*one( dat(ia,:)%kidr==2 ),mom,cnt,var)
                    WRITE(name(im),'("kid  by age",tr3,I4)') ia
                    weights(im)=whour 
                    im=im+1
                end do            
            end if 
            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%kidr==1  .AND. move(MNA:MXAD,:)>=0 ),   d1* move(MNA:MXAD,:) ,mom,cnt,var)	
            write(name(im),'("move | nokid ",tr1)')  
            weights(im)=wmovebyrel 
            im=im+1
            
            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%kidr==2  .AND. move(MNA:MXAD,:)>=0 ),   d1* move(MNA:MXAD,:) ,mom,cnt,var)	
            write(name(im),'("move | kid ",tr3)')  
            weights(im)=wmovebyrel
            im=im+1
            
            call condmom(im,( cosexrel(MNA:MXA,:) .AND. dat(MNA:MXA,:)%kidr==1 .AND. dat(MNA:MXA,:)%hhr>=0 ),   d1*one( dat(MNA:MXA,:)%hhr==1 ),mom,cnt,var)		
            write(name(im),'("e | nokid",tr5)') 
            weights(im)=whour 
            im=im+1 
            
            call condmom(im,( cosexrel(MNA:MXA,:) .AND. dat(MNA:MXA,:)%kidr==2 .AND. dat(MNA:MXA,:)%hhr>=0 ),   d1*one( dat(MNA:MXA,:)%hhr==1 ),mom,cnt,var)		
            write(name(im),'("e |   kid",tr5)') 
            weights(im)=whour 
            im=im+1 


            end do !j rel
        end do !g sex
        !END OF MOMENTS BY SEX AND REL 
        !*****************************************************************************************

        !*****************************************************************************************
        !START EVERYONE MOMENTS
        headloc(ihead)=im
        headstr(ihead)='getmar and getdiv rates'
        ihead=ihead+1
        do ia=mna,mxad,30 !ahu030622  changed from maxai-1 to mxad
            call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel==0 .AND. dat(ia+1,:)%rel>=0   .AND. dat(ia,:)%edr==1 ), d1*one(dat(ia+1,:)%rel==1),mom,cnt,var)
            write(name(im),'("getmarbyia,ned",tr1,i4)') ia
            weights(im)=wrel
            im=im+1
            call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel==0 .AND. dat(ia+1,:)%rel>=0   .AND. dat(ia,:)%edr==2 ), d1*one(dat(ia+1,:)%rel==1),mom,cnt,var)
            write(name(im),'("getmarbyia, ed",tr1,i4)') ia
            weights(im)=wrel
            im=im+1
        end do      
        do ia=mna,mxad,30 !ahu030622 changed from maxai-1 to mxad
            call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel==1 .AND. dat(ia+1,:)%rel>=0 ), d1*one(dat(ia+1,:)%rel==0),mom,cnt,var)
            write(name(im),'("get div by ia",tr1,i4)') ia
            weights(im)=wrel
            im=im+1
        end do            

        do i=1,nl,8            
            call condmom(im,( coho(MNA:MXAD,:) .AND. dat(MNA+1:MXA,:)%l==i .AND. move(MNA:MXAD,:)==1 ),   d1*one( dat(MNA+1:MXA,:)%l==dat(MNA+1:MXA,:)%hme  ),mom,cnt,var)		
            write(name(im),'("%hme-mvs to",tr3,i4)') i
            weights(im)=wmove 
            im=im+1 
        end do  
        
        headloc(ihead)=im; headstr(ihead)='after the first move what proportion is return to home? (all)';ihead=ihead+1     
        CALL condmom(im,( coho(MNA:MXAD,:).AND.  (moverank(MNA:MXAD,:)>1).AND.(move(MNA:MXAD,:)==1)  ),d1*one(homemove(MNA:MXAD,:)==0),mom,cnt,var)
        WRITE(name(im),'("non-home ")') 
        weights(im)=0.0_dp !wmove0
        im=im+1 
        CALL condmom(im,(coho(MNA:MXAD,:).AND.  (moverank(MNA:MXAD,:)>1).AND.(move(MNA:MXAD,:)==1)  ),d1*one(homemove(MNA:MXAD,:)==1),mom,cnt,var)
        WRITE(name(im),'("home     ")') 
        weights(im)=0.0_dp !wmove0
        im=im+1 
        
        do j=1,NL
            CALL condmom(im,( coho(MNA:MXAD,:) .AND.(move(MNA:MXAD,:)==1).AND.(dat(MNA:MXAD,:)%l>=0).AND.(norelchg(MNA:MXAD,:)==1) ),d1*one(dat(MNA:MXAD,:)%l==j),mom,cnt,var)
            WRITE(name(im),'("prop-of-moves-from        ",I4)') j
            weights(im)=0.0_dp !wmove0
            im=im+1 
        end do  
        do j=1,NL
            CALL condmom(im,( coho(MNA:MXAD,:).AND.(move(MNA:MXAD,:)==1).AND.(dat(MNA+1:MXA,:)%l>=0).AND.(norelchg(MNA:MXAD,:)==1) ),d1*one(dat(MNA+1:MXA,:)%l==j),mom,cnt,var)
            WRITE(name(im),'("prop-of-moves-to          ",I4)') j
            weights(im)=0.0_dp !wmove0
            im=im+1 
        end do  

            !headloc(ihead)=im; headstr(ihead)=' ';ihead=ihead+1		
        do i=1,nl
            call condmom(im,( coho(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%l==i  .AND. move(MNA:MXAD,:)>=0 ),   d1* move(MNA:MXAD,:) ,mom,cnt,var)		
            write(name(im),'("mvfr | loc",tr4,i4)') i
            weights(im)=wmove 
            im=im+1 
        end do

        do i=1,nl
            call condmom(im,( coho(MNA:MXAD,:) .AND. dat(MNA+1:MXA,:)%l==i .AND. move(MNA:MXAD,:)>=0 ),   d1* move(MNA:MXAD,:) ,mom,cnt,var)		
            write(name(im),'("mvto | loc",tr4,i4)') i
            weights(im)=wmove
            im=im+1 
        end do 

        do g=minsex,maxsex
            cosex(MNAD:MXA,:)= (dat(MNAD:MXA,:)%co==co .and. dat(MNAD:MXA,:)%sexr==g  )
            do i=1,5,2
                call condmom(im,( cosex(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0  .AND. dur(MNA:MXAD,:)==i),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                write(name(im),'("e|u dur,sex",2i4)') i,g
                weights(im)=whour
                im=im+1 
            end do 
            do i=1,5,2
                call condmom(im,( cosex(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%logwr>=0  .AND. dur(MNA:MXAD,:)==i ),   d1*dat(MNA+1:MXA,:)%logwr ,mom,cnt,var)		
                write(name(im),'("w|u dur,sex",2i4)') i,g
                weights(im)=wwage
                im=im+1 
            end do 
        end do !g sex

        do g=minsex,maxsex
            cosex(MNAD:MXA,:)= (dat(MNAD:MXA,:)%co==co .and. dat(MNAD:MXA,:)%sexr==g  )
            ia=MNA
            do i=1,nl
                call condmom(im,( cosex(ia,:) .AND. dat(ia,:)%hhr==1 .AND. dat(ia,:)%l==i .AND. dat(ia,:)%logwr>=0 ),   d1*dat(ia,:)%logwr ,mom,cnt,var)		
                write(name(im),'("w|loc a=18 dur,sex",2i4)') i,g	
                weights(im)=0.0_dp
                im=im+1 
            end do      
            ia=MNA
            do i=1,nl,8
                call condmom(im,( cosex(ia:ia+3,:) .AND. dat(ia:ia+3,:)%hhr==1 .AND. dat(ia:ia+3,:)%l==i .AND. dat(ia:ia+3,:)%logwr>=0 ),   d1*dat(ia:ia+3,:)%logwr ,mom,cnt,var)		
                write(name(im),'("w|loc a=18:21 dur,sex",2i4)') i,g
                weights(im)=wwage
                im=im+1 
            end do 
        end do !g sex
        !END OF EVERYONE MOMENTS
        !********************************************************

        !*****************************************************************************************
        !MOMENTS BY TYPE
        if (typemoments) then
        
            headloc(ihead)=im
            headstr(ihead)='mar by typ'
            ihead=ihead+1
            call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 ), d1*one(dat(MNA:MXA,:)%rel==1),mom,cnt,var)
            write(name(im),'("mar ",tr10)')	
            weights(im)=wrel !; if (onlysingles) weights(im)=0.0_dp   !forget about that if statement usually. just putting there to get im to go up by 1 for the next loop
            im=im+1
            do i=1,ntypp
                call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%typ==i ), d1*one(dat(MNA:MXA,:)%rel==1),mom,cnt,var)
                write(name(im),'("mar by typ ",i4)') i
                weights(im)=0.0_dp
                im=im+1
            end do         
            do i=1,ntypp
                call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%rel==1 ), d1*one(dat(MNA:MXA,:)%typ==i),mom,cnt,var)
                write(name(im),'("typ | mar ",i4)') i
                weights(im)=0.0_dp
                im=im+1
            end do 
    

    
            headloc(ihead)=im
            headstr(ihead)='move by typ,sex'
            ihead=ihead+1
            do g=minsex,maxsex
                cosex(MNAD:MXA,:)= (dat(MNAD:MXA,:)%co==co .and. dat(MNAD:MXA,:)%sexr==g  )    
                call condmom(im,( cosex(MNA:MXAD,:) .AND.  move(MNA:MXAD,:)>=0  ),   d1* move(MNA:MXAD,:) ,mom,cnt,var)	
                write(name(im),'("mv sex ",I4)') g
                weights(im)=wmove
                im=im+1
                do i=1,ntypp
                    call condmom(im,( cosex(MNA:MXAD,:) .AND.  move(MNA:MXAD,:)>=0  .AND. dat(MNA:MXAD,:)%typ==i ),   d1* move(MNA:MXAD,:) ,mom,cnt,var)	
                    write(name(im),'("mv typ,sex ",2I4)') i,g
                    weights(im)=0.0_dp
                    im=im+1
                end do !type
            end do !g sex


            headloc(ihead)=im
            headstr(ihead)='move by typ,rel,sex'
            ihead=ihead+1
            do g=minsex,maxsex
                do j=maxrelo,0,-1
                    if ( onlysingles ) then !  (.not.onlysingles).or.(onlysingles.and.j==0) ) then 
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g )			                    
                    else 
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
                    end if 
                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND.  move(MNA:MXAD,:)>=0 ),   d1* move(MNA:MXAD,:) ,mom,cnt,var)	
                    write(name(im),'("mv rel,sex ",3I4)') j,g
                    weights(im)=wmovebyrel
                    im=im+1
        
                    do i=1,ntypp
                    call condmom(im,( cosexrel(MNA:MXAD,:) .AND.  move(MNA:MXAD,:)>=0  .AND. dat(MNA:MXAD,:)%typ==i ),   d1* move(MNA:MXAD,:) ,mom,cnt,var)	
                    write(name(im),'("mv typ,rel,sex ",3I4)') i,j,g
                    weights(im)=0.0_dp 
                    im=im+1
                    end do !type
                end do !j rel
            end do !g sex 


            headloc(ihead)=im
            headstr(ihead)='wage by type,rel,sex'
            ihead=ihead+1
            do g=minsex,maxsex
                do j=maxrelo,0,-1
                    if ( onlysingles ) then !  (.not.onlysingles).or.(onlysingles.and.j==0) ) then 
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g )			                    
                    else 
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
                    end if 
    
                    CALL condmom(im,(   cosexrel(MNA:MXAD,:) .AND.  dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%edr==1 .AND. dat(MNA:MXAD,:)%logwr>=0  ) ,d1*dat(MNA:MXAD,:)%logwr,mom,cnt,var)
                    WRITE(name(im),'("w |ned ",2I4)') j,g
                    weights(im)=0.0_dp
                    im=im+1                        
                    do i=1,ntypp
                        CALL condmom(im,(   cosexrel(MNA:MXAD,:) .AND.  dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%edr==1 .AND. dat(MNA:MXAD,:)%logwr>=0  .AND. dat(MNA:MXAD,:)%typ==i ) ,d1*dat(MNA:MXAD,:)%logwr,mom,cnt,var)
                        WRITE(name(im),'("w|ned typ,rel,sex",3I4)') i,j,g
                        weights(im)=0.0_dp
                        im=im+1
                    end do

                    CALL condmom(im,(   cosexrel(MNA:MXAD,:) .AND.  dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%edr==2 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%typ==i ) ,d1*dat(MNA:MXAD,:)%logwr,mom,cnt,var)
                    WRITE(name(im),'("w|ed rel,sex ",2I4)') j,g
                    weights(im)=0.0_dp
                    im=im+1
                    do i=1,ntypp
                        CALL condmom(im,(   cosexrel(MNA:MXAD,:) .AND.  dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA:MXAD,:)%edr==2 .AND. dat(MNA:MXAD,:)%logwr>=0 .AND. dat(MNA:MXAD,:)%typ==i ) ,d1*dat(MNA:MXAD,:)%logwr,mom,cnt,var)
                        WRITE(name(im),'("w|ed typ,rel,sex",3I4)') i,j,g
                        weights(im)=0.0_dp
                        im=im+1
                    end do !typ
                
                end do !j rel
            end do !g sex




            headloc(ihead)=im
            headstr(ihead)='e|u,dur by typ,rel,sex '
            ihead=ihead+1
            do g=minsex,maxsex
                do j=maxrelo,0,-1
                    if ( onlysingles ) then !  (.not.onlysingles).or.(onlysingles.and.j==0) ) then 
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g )			                    
                    else 
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
                    end if 
                    do jj=1,5,4
                        call condmom(im,( cosex(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0  .AND. dur(MNA:MXAD,:)==jj),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                        write(name(im),'("e|u dur,rel,sex",3i4)') jj,j,g
                        weights(im)=whour
                        im=im+1 
                    end do 
                    do i=1,ntypp
                        do jj=1,5,4
                            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0  .AND. dur(MNA:MXAD,:)==jj  .AND. dat(MNA:MXAD,:)%typ==i),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
                            write(name(im),'("e|u dur,typ,rel,sex ",4i4)') jj,i,j,g
                            weights(im)=0.0_dp 
                            im=im+1 
                        end do 
                    end do 

                end do !rel
            end do !sex
            

            headloc(ihead)=im
            headstr(ihead)='w|u,dur by typ,rel,sex '
            ihead=ihead+1
            do g=minsex,maxsex
                do j=maxrelo,0,-1
                    if ( onlysingles ) then !  (.not.onlysingles).or.(onlysingles.and.j==0) ) then 
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g )			                    
                    else 
                        cosexrel(MNAD:MXAD,:)= (dat(MNAD:MXAD,:)%co==co .and. dat(MNAD:MXAD,:)%sexr==g  .and. dat(MNAD:MXAD,:)%rel==j .and. norelchg(MNAD:MXAD,:)==1 )			
                    end if 

                    do jj=1,5,4
                        call condmom(im,( cosex(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%logwr>=0  .AND. dur(MNA:MXAD,:)==jj ),   d1*dat(MNA+1:MXA,:)%logwr ,mom,cnt,var)		
                        write(name(im),'("w|u dur rel,sex",3I4)') jj,j,g
                        weights(im)=wwage
                        im=im+1 
                    end do 
                    do i=1,ntypp
                        do jj=1,5,4
                            call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%logwr>=0  .AND. dur(MNA:MXAD,:)==jj   .AND. dat(MNA:MXAD,:)%typ==i ),   d1*dat(MNA+1:MXA,:)%logwr ,mom,cnt,var)		
                            write(name(im),'("w|u dur typ,rel,sex ",4i4)') jj,i,j,g
                            weights(im)=0.0_dp
                            im=im+1 
                        end do
                    end do             
                end do !rel
            end do !sex
            

        end if !typemoments
        
        
    
    
    
    ENDDO cohort
    
    
    
    
    end subroutine get_mom
    end module mom
	

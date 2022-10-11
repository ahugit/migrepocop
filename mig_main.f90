!what are all these compile options?
!#mpif90 -O3 -fPIC -unroll -ip -axavx -xsse4.2 -openmp -vec-report -par-report -openmp-report -o train.x program.f90
!mpif90 -O3 -fPIC -unroll -ip -axavx -xsse4.2 -openmp -vec -par -openmp -o train.x program.f90
!#mpif90 train.x program.f90
	
program main 
	use params !, only: myid,iter,parname,npars,stepmin
	use objf
	use pnelder_mead  
    
	!use alib, only: logitinv
	implicit none	
	!include 'mpif.h'
	real(8), dimension(npars) :: pars1,pars2,pars3
	integer, parameter :: nj=10,nj2=2*nj
	integer :: i,j,k,pp,ierror,k1,k2
	! declarations for simulated annealing: at temperature t, a new point that is worse by an amount will be accepted with probability exp(-q/t)
	real(8) :: pars(npars),incr,realpars(npars),realpars1(npars),val
	real(sp) :: q1val(nj2),step(nj2),lb,ub
    integer :: maxfn,iprint,nloop,iquad,ifault
    real(8) :: qval,stopcr,simp,var(npars)
	!ahu 121818: the below are no longer declared as parameters. their value is set below in the program. 
    real(8) :: tstart !=15. 			! starting temperature (set to zero to turn off simulated annealing)
	real(8) :: tstep !=0.8			! fraction temperature reduced at each step
	integer :: tfreq !=30			! number of function calls between temperature reductions
	integer :: saseed !=1			! seed for random numbers
	!ahu 121818: the above are no longer declared as parameters. their value is set below in the program. 
    
    ! for the new sa routine that does grouping
	integer, dimension(0:ninp-1) :: myhomes != (/ 1,2,3,4,5,6,7,8,9 /) !,1,2,3,4,5,6,7,8,9 /)              ! types !mytyps = (/1,1,2,2/)
	integer, dimension(0:ninp-1) :: mytyps  != 1 !(/ 1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2 /)              ! types 
	integer, dimension(0:ninp-1) :: mycos != 1
	!integer, dimension(0:ninp-1) :: myindexs != (/ 1,2,3,4,5,6,7,8,9 /)  !,10,11,12,13,14,15,16,17,18 /)    ! types    !mytyps = (/1,1,2,2/)
	! keep track of time 
	real(4):: begintime,runtime
	real(8) :: temp
	integer :: numgroup,numworld,iam,mygrank0 !mygrank,ngroups,mygroup,nprocs_world,iam
	integer :: mpi_group_world,mygsize
	integer :: mpierr,mpistat(mpi_status_size)
    real(4), allocatable, dimension(:,:) :: mytime 
	call mpi_init(mpierr)
	call mpi_comm_rank(mpi_comm_world,iam,mpierr)
	call mpi_comm_size(mpi_comm_world,numworld,mpierr)
	print*, "numworld,Iam ", numworld,iam
    mysay=iam
	conditional_moments=.true.		
    if (iam==0) print*, "Here is numworld", numworld

    allocate(mytime(numworld,2))
   
	if (groups) then 	
		if (mod(numworld,ninp)>0) then ; print*, "numworld needs to be a multiple of ninp " ; end if 

		numgroup=int(numworld/ninp)		!number of groups  36/18=2 
		mygroup=mod(iam,numgroup) ! which group each process belongs to: myid mygroup: 0 0 / 1 1 / 2 0 / 3 1 / 4 0 / 5 1 / ..... / 34 0 / 35 1 
		mygrank=int(iam/numgroup) !id within group: myid mygroup: 0 0 / 1 0 / 2 1 / 3 1 / 4 2 / 5 2 / 6 3 / 7 3 / 8 4 / 9 4 / ..../ 34 17 / 35 17 
        !ahu 0317: BEWARE OF DOING THINGS LIKE BELOW. IF OYU SET THE MYGROUP THE WAY IT IS DONE BELOW THEN THE FUNCTION_DISTRIBUTE 
        !WILL NOT WORK PROPERLY. WHY? BEACAUSE MPI_BCAST DOES NOT DO WHAT IT SHOULD DO IN THIS CASE. 
        !BECAUSE THEN FOR I=1,NPROCS WHERE NPROCS IS THE NUMBER OF GROUPS, FOR I=1, THE ROOT PROCESS IS 0 AND FOR I=2, THE ROOT PROCESS IN BCAST IS 1. 
        !IF YOU DO BELLOW, WHEN BCAST SAYS SEND VALUE FROM ROOT PROCESS 1 TO EVERYONE, IT DOES NOT TAKE ROOT PROCESS 1 TO BE THE PROCEESSSES IN GROUP 2 
        !BUT IT JUST TAKES IT TO BE THE PROCESSOR WITH IAM=1. SO THEN WHEN YOU'RE GOING THROUGH I=1,NPROCS IN THE BCAST PART OF FUNCTION_DISTRIBUTE
        !YOU ARE NOT ACTUALLY GOING THROUGH GROUPS. 
        !BUT IF YOU DO THE ABOVE SETTING OF MYGROUP, THEN THE TWO ARE EQUIVALENT SINCE THE IAM'S ALSO CORRESPOND TO GROUP NUMBERS. (JUST DUE TO THE WAY MYGROUP IS SET)
        !numgroup=int(numworld/ninp)		!number of groups  36/18=2
        !mygroup=int(iam/ninp) ! which group each process belongs to: myid mygroup: 0 0 / 1 0 / 2 0 / 3 0 / 4 0 /.../17 0 / 18 1 /...../35 1 
		!mygrank=mod(iam,ninp) !id within group: myid mygroup: 0 0 / 1 1 / 2 2 / 3 3 / 4 4 / 5 5 / 6 6 / 7 7 / 8 8 ..../ 17 17 / 18 0 / 19 1 / ..../35 17 
       
        call mpi_comm_split(mpi_comm_world,mygroup,mygrank,comm,mpierr)
        !ahu 0317 call mpi_comm_split(mpi_comm_world,mygroup,iam,comm,mpierr) !ahu 0317
        !print*, "after split ", iam
        !call MPI_Comm_rank(comm, mygrank,mpierr) !ahu 0317
        !print*, "after rank ", iam
        !call MPI_Comm_size(comm, mygsize,mpierr) !ahu 0317
        !print*, "after size ", iam
        !if (mygrank0.ne.mygrank) then 
        !    print*, 'The two ways of calculating mygrank are not equal'
        !    stop
        !end if 
       ! write(*,'("iam,numgroup,mygroup,mygrank,mygsize:")') !ahu 0317
       ! write(*,'(5i4)') iam,numgroup,mygroup,mygrank,mygsize !ahu 0317
        
        pp=0
        do i=1,ncop
            do j=1,ntypp
                do k=1,nhomep
                    mycos(pp)=i
                    mytyps(pp)=j
                    myhomes(pp)=k
                    !myindexs(pp)=pp
                    pp=pp+1
                end do 
            end do 
        end do 
        if ( (pp-1).ne.(ncop*ntypp*nhomep-1)) then
            print*, 'something wrong with pp',pp,ncop*ntypp*nhomep
            stop
        end if
        if (iam==0) then 
            write(*,'("mycos,mytyps,myhomes,myindexs")') 
            if (ntypp==1) then 
                write(*,'(9I4)') mycos
                write(*,'(9I4)') mytyps
                write(*,'(9I4)') myhomes
                !write(*,'(9I4)') myindexs
            else if (ntypp==2) then
                write(*,'(18I4)') mycos
                write(*,'(18I4)') mytyps
                write(*,'(18I4)') myhomes
                !write(*,'(18I4)') myindexs
            end if
        end if
		
		!myindex = myindexs(mygrank)
		myco   = mycos(mygrank)
		mytyp  = mytyps(mygrank)
		myhome = myhomes(mygrank)
	else 
		numgroup=0 ; mygroup=0 ; mygrank=0 ; myco=0 ; mytyp=0 ; myhome=0
    end if
	begintime=secnds(0.0)
    iter=1
	iwritegen=0
	if (iam==0) iwritegen=1
	if (iwritegen==1.and.chkobj) then						! ahu april13: to check how much objval changes with each stepsize and make sure that the objval changes by similar amounts when changing each parameter. otherwise the simplex does weird things. and also to check how much the pminim routine changes objval with each iteration during estimation but the latter is not as important. 
		open(unit=63, file='chkobj.txt',status='replace')	
	end if 
	if (iam==0.and.chkstep) then 
		open(unit=64, file='chkstep.txt',status='replace')
		open(unit=65, file='step.txt',status='replace')
	end if 
	if (iam==0.and.comparepars) then 
		open(unit=64, file='chkstep.txt',status='replace')
	end if     
	!ahu030622 skriv=.false. ; if (iam==24.and.iter==1) skriv=.false.  !ahu 0327 iam=9
	skriv=.false. ; if (iam==0.and.iter==1) skriv=.false.  !ahu030622
    if (.not.groups) skriv=.false. !ahu 0327
    iwritemom=0 ; if (iam==0) iwritemom=1 
	if (skriv) then 	
		open(unit=12, file='chkdat.txt',status='replace') !written in call getmom/yaz_getmom	
		open(unit=40,file='chkq.txt',status='replace') 
		open(unit=50,file='chkpars.txt',status='replace')    !written in objf/yaz0
		open(unit=100,file='chksolgetdec_s.txt',status='replace')
        !written in sol/couples/getdec/yaz_getdec (callfrom==40), sol/marmkt/getdec/yaz_getdec (callfrom=50)
		open(unit=200,file='chksolgetdec_c.txt',status='replace') 
		open(unit=201,file='chk2b.txt',status='replace')		
        !400 is written in simulate, simulate/getdec/yaz_getdec (callfrom=80), yaz_decs,yaz_sim,yaz_simmatch,yaz_simdecmar 
        open(unit=400, file='chksim.txt',status='replace')   
		open(unit=500, file='chksimpath.txt',status='replace')	 !written by simulate/yaz_simpath
	end if 
	!call cohab_stderr(parvector,stepsize,stderrs)
nonneg=.TRUE.
!ahu030622 nonlabinc=(/ 300.0_dp,1100.0_dp /) 
nonlabinc=0.0_dp !ahu030622



    if (icheck_eqvmvf) then
		pars(3)=pars(2)
		pars(8)=pars(7)
		pars(10)=pars(9)
		pars(34:45)=pars(22:33)
		pars(47)=pars(46)
		pars(59:60)=pars(57:58)
		pars(53:56)=pars(49:52)
	end if 	
	if (icheck_eqvcvs) then 
		pars(1)=-huge(0.0)       ! sig_mar
		pars(6)=-huge(0.0)       ! divpenalty
		pars(npars)=0.0_dp          ! mu_mar
	end if 
    
      
    !open(unit=2,file='bp041019_3.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    !nonlabinc=0.0_dp
    !The below are in order to set type parameters to values they were before we changed the getpars setup 
    !(where ptype's are no longer 0 and mumar's are no longer all being set equal to mumar(1))
    !pars(75)=0.0015_dp !mumar1
    !pars(76)=0.0_dp !ptypehs2 
    !pars(77)=0.0_dp !ptypecol2
    !pars(81)=0.0015_dp !mumar2
    !pars(82)=0.0_dp !ptypehs3 
    !pars(83)=0.0_dp !ptypecol3
    !pars(87)=0.0015_dp !mumar3
    !pars(88)=0.0_dp !ptypehs4 
    !pars(89)=0.0_dp !ptypecol4
    !pars(93)=0.0015_dp !mumar4
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp     


    !open(unit=2,file='bp090922.txt',status='old',action='read') ; read(2,*) pars1	; close(2)
    nonlabinc=0.0_dp
    pars=pars1
    !changing only type parameters in pars
    !pars(70)=0 !ptype set to 0 in getpars
    !pars(71)=0 !ptype set to 0 in getpars
    pars(72)=8.9_dp !alf1t
    pars(73)=8.9_dp !alf2t
    !pars(74)=0 !cst1 this is set in getpars as it is a normalization
    pars(75)=0.1_dp !mumar1 !should this be set to 0? 

    pars(76)=0.1_dp !ptypehs2 
    pars(77)=0.1_dp !ptypecol2
    pars(78)=9.1_dp !alf1t
    pars(79)=8.9_dp !alf2t
    pars(80)=-100.0_dp !cst2 
    pars(81)=0.1_dp !mumar2    
    
    pars(82)=0.1_dp !ptypehs3 
    pars(83)=0.1_dp !ptypecol3
    pars(84)=9.8_dp !alf1t
    pars(85)=9.1_dp !alf2t
    pars(86)=-500.0_dp !cst3
    pars(87)=0.1_dp !mumar3
    
    pars(88)=-0.6_dp !ptypehs4 
    pars(89)=0.1_dp !ptypecol4
    pars(90)=9.4_dp !alf1t
    pars(91)=9.2_dp !alf2t
    pars(92)=8000.0_dp !cst4 
    pars(93)=0.1_dp !mumar4
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp        




    open(unit=2,file='bp093022.txt',status='old',action='read') ; read(2,*) pars1	; close(2)
    nonlabinc=0.0_dp
    pars=pars1
    pars(75)=0.0005_dp !mumar1
    pars(81)=0.0005_dp !mumar2
    pars(87)=0.0005_dp !mumar3
    pars(93)=0.01_dp    !mumar4
    pars(68:69)=10000.0_dp   
    !pars(74)= cst1 norm to 0 in getpars
    pars(80)=-30000.0_dp
    pars(86)=-30000.0_dp
    pars(92)=-30000.0_dp 
    pars(21)=-2.5_dp !pmeet to bring down getmar rates
    pars(10)=1.5_dp !u of m to bring up e u cond on move transitions (it's currentlt -1.2 which is 0.2 sth)
    pars(12)=0.1_dp !u of fem to bring up e u cond on move transitions (it's currently -2.7 which is 0.06)
    pars(52)=-1.8_dp !increasing alf12 because of nexp
    pars(64)=-2.2_dp !increasing alf22 because of nexp
    !pars(9)=0.5_dp !u cur m (to bring up eu move and ee move so that they don't say I'll find a job after I movey)
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp   


    !set all junk parameters to zero so that I get zeros for these parameters the next pars file after opt
    pars(17:19)=0.0_dp !psih2-4
    pars(53)=0.0_dp !alf13 
    pars(65)=0.0_dp !alf23
 
!************************increase type 3 cst to decrease getdiv (w/0 messing getmar) AND sing mv

!************************increase type 3 cst to decrease getdiv (w/0 messing getmar) AND sing mv
    onthejobsearch=.TRUE.
    pars(66)=pars1(66)+1.5_dp
    pars(67)=pars1(67)+1.0_dp
    pars(68)=50000.0_dp   
    pars(69)=20000.0_dp   
    pars(74)=-10000.0_dp !cst1 no longer norm 
    pars(75)=0.0002_dp !mumar1
    pars(80)=-10000.0_dp
    pars(81)=0.0002_dp !mumar2
    pars(86)=-10000.0_dp
    pars(87)=0.001_dp !mumar3
    pars(92)=-50000.0_dp 
    pars(93)=0.01_dp    !mumar4
    pars(27)=-3.0_dp !alphaed(m,noed)
    pars(28)=-2.0_dp !alphaed(f,noed)
    pars(29)=-5.3_dp !alphaed(m,ed)
    pars(30)=-4.0_dp !alphaed(f,ed)
    pars(86)=-50000.0_dp 
    pars(87)=0.001_dp !mumar3
    !pars(21)=-1.5_dp !pmeet
    pars(72:73)=8.8
    pars(78:79)=8.8
    pars(84:85)=8.8
    pars(90:91)=8.8
    pars(87)=0.01
    pars(93)=0.01
    pars(74)=-80000.0_dp
    pars(80)=-80000.0_dp
    pars(86)=-30000.0_dp
    pars(92)=-30000.0_dp
!call getpars(pars,realpars)
!call objfunc(pars,qval) ; realpars=realpartemp   

!pars(87)=0.01
!pars(93)=0.01
!pars(74)=-80000.0_dp
!pars(80)=-80000.0_dp
!pars(86)=-60000.0_dp
!pars(92)=-60000.0_dp
!call getpars(pars,realpars)
!call objfunc(pars,qval) ; realpars=realpartemp   

!pars(87)=0.01
!pars(93)=0.01
!pars(74)=-80000.0_dp
!pars(80)=-80000.0_dp
!pars(86)=-10000.0_dp
!pars(92)=-10000.0_dp
!call getpars(pars,realpars)
!call objfunc(pars,qval) ; realpars=realpartemp   


!pars(87)=0.001
!pars(93)=0.001
!pars(74)=-80000.0_dp
!pars(80)=-80000.0_dp
!pars(86)=-60000.0_dp
!pars(92)=-60000.0_dp
!call getpars(pars,realpars)
!call objfunc(pars,qval) ; realpars=realpartemp   


!pars(87)=0.1
!pars(93)=0.1
!pars(74)=-80000.0_dp
!pars(80)=-80000.0_dp
!pars(86)=-10000.0_dp
!pars(92)=-10000.0_dp
!call getpars(pars,realpars)
!call objfunc(pars,qval) ; realpars=realpartemp   


!pars(87)=0.1
!pars(93)=0.1
!pars(74)=-80000.0_dp
!pars(80)=-80000.0_dp
!pars(86)=-60000.0_dp
!pars(92)=-60000.0_dp
!call getpars(pars,realpars)
!call objfunc(pars,qval) ; realpars=realpartemp   


pars(16)=pars1(16)-0.5_dp
pars(52)=pars1(52)-1.0_dp
pars(64)=pars1(64)-1.0_dp
pars(72:73)=8.8+0.4_dp
pars(78:79)=8.8+0.4_dp
pars(84:85)=8.8+0.4_dp
pars(90:91)=8.8+0.4_dp
pars(87)=0.1
pars(93)=0.1
pars(74)=-80000.0_dp
pars(80)=-80000.0_dp
pars(86)=-10000.0_dp
pars(92)=-10000.0_dp
pars(26)=pars1(26)+1.0_dp
pars(28)=-4.0_dp
pars(27)=-5.3
pars(16)=pars1(16)-1.2_dp
pars(52)=pars1(52)-2.5_dp
pars(68)=30000.0_dp   !decrease sigom
pars(69)=10000.0_dp   !decrease sigof
pars(1)=3.8417_dp     !higher empcurm because of chkobj101022
pars(7)=2.4740_dp     !higher emp of f because of chkobj101022
pars(33)=1358.7599_dp !higher uloc1 because of chkobj101022
pars(35)=3354.4870_dp !higher uloc3 because of chkobj101022
pars(36)=4763.7679_dp !higher uloc4 because of chkobj101022
pars(38)=5072.9306_dp !higher uloc6 because of chkobj101022
pars(39)=5328.9859_dp !higher uloc7 because of chkobj101022
pars(41)=4298.3866_dp !higher uloc9 because of chkobj101022


if (iwritegen==1) then
    print*, pars(52),logit(pars(52)),3.0_dp*logit(pars(52))
    pars(52)=-2.0_dp
    print*, pars(52),logit(pars(52)),3.0_dp*logit(pars(52))
    pars(52)=-1.5_dp
    print*, pars(52),logit(pars(52)),3.0_dp*logit(pars(52))
    pars(52)=-1.0_dp
    print*, pars(52),logit(pars(52)),3.0_dp*logit(pars(52))
    pars(52)=-0.5_dp
    print*, pars(52),logit(pars(52)),3.0_dp*logit(pars(52))
    pars(52)=0.0_dp
    print*, pars(52),logit(pars(52)),3.0_dp*logit(pars(52))
end if

onthejobsearch=.TRUE.
pars(88:89)=10.0_dp !make everyone type4
pars(92)=0.0_dp
do j=1,15
    pars(68)=1000.0_dp*j
    call getpars(pars,realpars)
    call objfunc(pars,qval) ; realpars=realpartemp   
end do 

pars(92)=-5000.0_dp
do j=1,15
    pars(68)=1000.0_dp*j
    call getpars(pars,realpars)
    call objfunc(pars,qval) ; realpars=realpartemp   
end do 

!pars(92)=-3000.0_dp
!do j=1,5
!!    pars(68)=1000.0_dp*j
!    call getpars(pars,realpars)
!    call objfunc(pars,qval) ; realpars=realpartemp   
!end do

!*************************    
    !mytime(iam+1,1)=secnds(0.0)
    !mytime(iam+1,2)=secnds(mytime(iam+1,1))        
    !if (iam==0) print*, 'Here is qval: ', qval
    !if (iam==0) print*, 'iam,qval,mytime ', iam,qval,mytime(iam+1,2)
    deallocate(mytime)

	if (optimize) then 
		! call simplex: the call to minim (the simplex routine) is set up to work with starting vector parvector
		!if (iam==0) then ; iprint=1 ; else ; iprint=-1 ; end if 
		!open(unit=64,file='step11_10.txt',status='old',action='read') ; read(64,*) stepmin	; close(64) 			
		!open(unit=64,file='step111416.txt',status='old',action='read') ; read(64,*) stepmin	; close(64) 			
        !open(unit=64,file='stepp.txt',status='old',action='read') ; read(64,*) stepmin	; close(64) 
        !stepmin=0.5_dp

        !optimization parameters:
        stepmin=stepos !ahu 121118
        maxfn=8000
        !if (iam==0) then  
        iprint=20
        !    open(unit=6538, file='lavas.txt',status='replace')		
        !else 
        !    iprint=-1
        !end if 
        stopcr=20.0_dp
        nloop=npars+2
        iquad=1 !ag092522 agsept2022: This was 0 before but I don't think that's right so I'm trying this one. actually I don't think it matters (it only matters at end after convgence)
        simp=0.0_dp
        !sim annealing parameters: 
        tstart=0.3_dp*qval  !T0: starting temp (can set to zero to turn off sim annealing)
        tstep=0.5_dp        !fraction temp reduced at each tfreq
        tfreq=COUNT(stepmin /= zero) !number of function calls between temp reductions
        saseed=1            !seed for random numbers

        !call minim(p,    step,    nop, func,  maxfn,  iprint, stopcr, nloop, iquad,  simp, var, functn, ifault)
        !call minim(pars, stepmin, npars, qval, maxfn, iprint, stopcr, nloop, iquad,  simp, var, objfunc, writebest, ifault)    
    
    	if (groups) then 
			call pminim(pars, stepmin, npars, qval, maxfn, iprint, stopcr, nloop, iquad, simp, var, objfunc, writebest,ifault,mygroup,numgroup,tstart, tstep	, tfreq, saseed)
		else 
			call pminim(pars, stepmin, npars, qval, maxfn, iprint, stopcr, nloop, iquad, simp, var, objfunc, writebest,ifault,iam,numworld,tstart, tstep	, tfreq, saseed)
		end if 

   
		if (iam==0) then ; print*, "out of minim now and here is ifault ", ifault ; end if 	
        if (iam==0) then  
            close(6538)
        end if 
    else 
		!conditional_moments=.false.		
		conditional_moments=.true.	       
        stepmin=stepos

        !call random_seed (size=p)
	    !p=12
	    !call random_seed (put = (/(k,k=1,p)/))
	    !if (skriv) then
	    !    print*, 'Here is p',p
	    !    print*, 'Here is k',(/(k,k=1,p)/)
     !			call random_number(rand)
	!		print*, 'rand: ', rand
    ! 			call random_number(rand)
!			print*, 'rand: ', rand
 !       end if
	!else     
	!    allocate(newseed(p))
	!    newseed=10
	!    call random_seed( put=newseed(1:p)  )
	!    deallocate(newseed)
	!end if 
!		call objfunc(pars,q) ; realpars=realpartemp
 
        
        !open(unit=2,file='bestpar110516.txt',status='old',action='read') ; read(2,*) pars2	; close(2) 		
        !pars(4)=-13.7901
		!call objfunc(pars,q) ; realpars=realpartemp
        if (comparepars) then
            if (iam==0) then 
                write(64,*)
                write(64,*)
                write(64,'(tr3,"i",tr3,"j",tr2,"iter",2x,tr20,tr13,"q",tr9,"q1(j)",tr6,"pars",tr5,"pars1",5x,tr6,"realpars",5x,tr5,"realpars1",5x,tr3,"step(j)",5x,tr2,"abs(q1(j)-q)")')
            end if
            do i=1,npars
                !if (abs(stepmin(i)) >0.0_dp) 
                    pars1=pars
                    pars1(i)=pars2(i)
                    call objfunc(pars1,val) ; realpars1=realpartemp
                    j=1 ; q1val(j)=val !j is just a nuissance index here which is kept just to be consistent with the chkstep output
                    if (iam==0) then 
                    	write(64,'(2i4,i6,". ",1a20,2f14.2,2F10.2,2(5x,f14.2),(5x,f10.2),(5x,f14.2) )') i,j,iter-1,parname(i),qval,q1val(j),pars(i),pars1(i),realpars(i),realpars1(i),step(j),abs(q1val(j)-qval)
                    end if 				
                !end if 
            end do 
        end if 
        
		if (chkstep) then !don't forget to set stepmin equal to stepos (i.e. assign it some value! with getpars!)
            lb=0.4_sp ; ub=1.25_sp ; incr=0.2_dp
            do i=1,npars
                !if (  parname(i)=='cst(1)'.or.parname(i)=='cst(2)'.or.parname(i)=='ecst' ) then
                !    stepmin(i)=0.2_dp
                !else 
                !    stepmin(i)=0.0_dp
                !end if 
                if ( abs(stepmin(i)) > 0.0_dp ) then	
                !if (  parname(i)=='uloc'.or.  & 
                 ! & parname(i)=='alf13'.or.parname(i)=='alf22'.or.i==47.or.i==48 .or. parname(i)=='mu_o') then
                    !if (i==4) then !i==40.or.i==41.or.i==42.or.i==50.or.i==53.or.i==54.or.i==56.or.i==57.or.i==59.or.i==60.or.i==62.or.i==63.or.i>=65) then
					if (iam==0) then 
						write(64,*)
						write(64,*)
						write(64,'(tr3,"i",tr3,"j",tr2,"iter",2x,tr20,tr13,"q",tr9,"q1(j)",tr6,"pars",tr5,"pars1",5x,tr6,"realpars",5x,tr5,"realpars1",5x,tr3,"step(j)",5x,tr2,"abs(q1(j)-q)")')
					end if 
					step=pen 
					q1val=pen 
					k=0
					incr=0.2_dp
					jloop: do j=1,nj2
						if (parname(i)=='cst(1)'.or.parname(i)=='cst(2)') then
                            step(j)=(j-nj)*0.2_dp
						else if (parname(i)=='ecst') then
                            step(j)=(j-nj)*0.2_dp
						else if (parname(i)=='kcst') then
                            step(j)=(j-nj)*0.25_dp
						else if (parname(i)=='divpenalty') then
                            step(j)=(j-nj)*0.25_dp
						else if (parname(i)=='alphaed(m,1)') then
                            step(j)=(j-nj)*1.0_dp
						else if (parname(i)=='alphaed(m,2)') then
                            step(j)=(j-nj)*1.0_dp
						else if (parname(i)=='alphaed(f,2)') then
                            step(j)=(j-nj)*0.3_dp
						else if (parname(i)=='psio') then
                            step(j)=(j-nj)*0.3_dp
						else if (parname(i)=='psil') then
                            step(j)=(j-nj)*0.3_dp
						else if (parname(i)=='psih') then
                            step(j)=(j-nj)*0.3_dp
                        else if (parname(i)=='uloc') then 
							step(j)=(j-nj)*0.06_dp
                        else if (parname(i)=='alf11'.or.parname(i)=='alf21') then 
							step(j)=(j-nj)*0.1_dp
                        else if (parname(i)=='alf20'.and.i==47) then
                            step(j)=(j-nj)*0.05_dp
						else if (parname(i)=='alf20'.and.i==48) then
                            step(j)=(j-nj)*0.05_dp
						else if (parname(i)=='alf13') then
                            step(j)=(j-nj)*0.2_dp*abs(pars(i))
						else if (parname(i)=='alf22') then
                            step(j)=(j-nj)*0.2_dp*abs(pars(i))
                        !else if (i>=50.and.i<=67) then
						!	step(j)=(j-nj)*0.25_dp   
                        else 
                            step(j)=(j-nj)*incr*abs(pars(i))
						end if 						
						pars1=pars						
						pars1(i)=pars(i)+step(j)
						call objfunc(pars1,val) ; realpars1=realpartemp
						q1val(j)=val
						!if ( abs(q1(j)-q) < incr) incr=0.1_dp
						if (iam==0) then 
							write(64,'(2i4,i6,". ",1a20,2f14.2,2F10.2,2(5x,f14.6),(5x,f10.2),(5x,f14.2) )') i,j,iter-1,parname(i),qval,q1val(j),pars(i),pars1(i),realpars(i),realpars1(i),step(j),abs(q1val(j)-qval)
						end if 												
					end do jloop
					call sort2(q1val,step)
					j=0
					do while (lb+j*0.1_sp<=0.95_sp)
						k=locate(q1val,(lb+j*0.1_sp)*real(qval)) 
						if (iam==0) write(64,'("first try ",i6,3f14.4,i6)') j,q1val(1),q1val(nj2),(lb+j*0.1_sp)*real(qval),k
						if (k>0.and.k<nj2) exit 
						j=j+1
					end do
					j=0
					if (k==0) then 
						do while (ub-j*0.1_sp>=1.05_sp)
							k=locate(q1val,(ub-j*0.1_sp)*real(qval)) 
							if (iam==0) write(64,'("second try ",i6,3f14.4,i6)') j,q1val(1),q1val(nj2),(ub-j*0.1_sp)*real(qval),k
							if (k>0.and.k<nj2) exit 
							j=j+1
						end do 
					end if 					
					if (iam==0) write(64,'("and now it is done ",i6)') k
					if (k>0) then ; stepmin(i)=step(k) ; else ; stepmin(i)=-99.0 ; end if 
					if (iam==0) then 
						!write(64,*) "here is the best bumps for the initial simplex"
						write(64,*)
						write(64,'("lb,ub,q1min,q1max",4x,4f14.4)') lb*qval,ub*qval,q1val(1),q1val(nj2)
						write(64,*)
						if (k>0) then 
							write(64,'(2i6,4f14.4,4(5x,f14.4) )') i,k,qval,q1val(k),pars(i),pars1(i),realpars(i),realpars1(i),step(k),abs(q1val(k)-qval)
						end if 
						write(64,*) 
					end if 												

				end if 	
				if (iam==0) write(65,'(f14.4)') stepmin(i)
			end do 
		end if 
	end if 
	runtime=secnds(begintime)
	if (iam==0) then ; write(*,'("run time: ", f14.2)') runtime ; end if 
	!ag090122 agsept2022 if (groups) then 
	!ag090122 agsept2022 	call mpi_finalize(ierror)   
	!ag090122 agsept2022 end if 
    !ag090122 agsept2022 commenting out the aove groups if statment because 
    !I get mpi error when I run with groups FALSE saying some process couldnt finish because mpi finalze wasn't called etc.
    !since I invoke mpi regardless of groups, I have to call mpi finalize regardless of groups
    call mpi_finalize(ierror)   
end program main 

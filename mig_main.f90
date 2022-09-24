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
	real(8) :: pars(npars),qval,incr,realpars(npars),realpars1(npars),val
	real(sp) :: q1val(nj2),step(nj2),lb,ub
	integer, parameter :: maxfcn=8000	!ahu 121918		! max number of function calls
	integer :: iwriteminim,ind(1)
	real(8), parameter :: stopcr=0.01_dp			! convergence criterion
	integer, parameter :: nloop=npars,iquad=0			
	real(8), parameter :: simp=0.0_dp
	real(8) :: var(npars)
	integer :: ifault,track(npars)
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

	call mpi_init(mpierr)
	call mpi_comm_rank(mpi_comm_world,iam,mpierr)
	call mpi_comm_size(mpi_comm_world,numworld,mpierr)
	print*, "numworld,Iam ", numworld,iam
    mysay=iam
	conditional_moments=.true.		

   
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
        write(*,'("iam,numgroup,mygroup,mygrank,mygsize:")') !ahu 0317
        write(*,'(5i4)') iam,numgroup,mygroup,mygrank,mygsize !ahu 0317
        
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
	if (iam==0.and.optimize) then 
		open(unit=61, file='bestval.txt',status='replace')
	end if 
	if (iam==0.and.iwritegen==1) then						! ahu april13: to check how much objval changes with each stepsize and make sure that the objval changes by similar amounts when changing each parameter. otherwise the simplex does weird things. and also to check how much the pminim routine changes objval with each iteration during estimation but the latter is not as important. 
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
		open(unit=12, file='chkdat.txt',status='replace')		
		open(unit=40,file='chkq.txt',status='replace') 
		open(unit=50,file='chk0.txt',status='replace') 
		open(unit=100,file='chk1.txt',status='replace')
		open(unit=200,file='chk2.txt',status='replace')
		open(unit=201,file='chk2b.txt',status='replace')		
		!open(unit=300,file='chk3.txt',status='replace')		
        open(unit=400, file='chk4.txt',status='replace')		
		open(unit=500, file='chksimpath.txt',status='replace')		
        open(unit=88881, file='chkwage.txt',status='replace') 
        open(unit=56789, file='chktemp.txt',status='replace')
	end if 
	!call cohab_stderr(parvector,stepsize,stderrs)

!041618: STARTING NEW ESTIMATION NOW AFTER A 1 YEAR HIATUS! THE BELOW ARE AFTER TRIALS OF EPS2 PROBLEMS AND TRYING TO FIGURE OUT MOVE FOR MAR RATES ETC.
    !    nonlabinc=0.0_dp   !MAKE IT SO THAT THIS IS SET AT PARAMS
    !    nonneg=.TRUE.      !MAKE IT SO THAT THIS IS SET AT PARAMS
    !    eps2=eps           !MAKE IT SO THAT THIS IS SET AT PARAMS
    
!open(unit=2,file='bestpar041818.txt',status='old',action='read') ; read(2,*) pars1	; close(2)
!        pars1(12)=-0.1_dp     !pmeet
!open(unit=2,file='bestpar090618.txt',status='old',action='read') ; read(2,*) pars	; close(2)
!open(unit=2,file='bestpar110618.txt',status='old',action='read') ; read(2,*) pars1	; close(2)
!open(unit=2,file='bestpar112618.txt',status='old',action='read') ; read(2,*) pars1	; close(2)
!open(unit=2,file='bestpar121018.txt',status='old',action='read') ; read(2,*) pars	; close(2)
!open(unit=2,file='parnew_121118.txt',status='old',action='read') ; read(2,*) pars1	; close(2)
!open(unit=2,file='bestpar121518.txt',status='old',action='read') ; read(2,*) pars	; close(2)
!open(unit=2,file='par121718.txt',status='old',action='read') ; read(2,*) pars	; close(2)
!open(unit=2,file='bestpar121718s.txt',status='old',action='read') ; read(2,*) pars	; close(2)
!pars(22:23)=-2.3_dp
!pars(29)=-3.0_dp !alphaed(m ed)
nonneg=.TRUE.
!ahu030622 nonlabinc=(/ 300.0_dp,1100.0_dp /) 
nonlabinc=0.0_dp !ahu030622


!open(unit=2,file='bestpar121818_ntyp1_sigo.txt',status='old',action='read') ; read(2,*) pars	; close(2)
!open(unit=2,file='bestpar121918_1.txt',status='old',action='read') ; read(2,*) pars1	; close(2)
!open(unit=2,file='bestpar121918_2.txt',status='old',action='read') ; read(2,*) pars	; close(2)
!open(unit=2,file='bestpar122118_ntypp1.txt',status='old',action='read') ; read(2,*) pars1	; close(2)
!open(unit=2,file='bestpar122118_ntypp4.txt',status='old',action='read') ; read(2,*) pars2	; close(2)
!open(unit=2,file='parnew122118.txt',status='old',action='read') ; read(2,*) pars	; close(2)
!open(unit=2,file='bestpar122318_xtramom.txt',status='old',action='read') ; read(2,*) pars	; close(2)
!pars(20)=-1.0_dp !pkid because not iterating on pkid when onlysingles because the kid moments by age for singles looks weird in the data (see in get_mom the note)

!open(unit=2,file='bestpar122318_wgt0.txt',status='old',action='read') ; read(2,*) pars	; close(2)
!pars(20)=-1.0_dp !pkid because not iterating on pkid when onlysingles because the kid moments by age for singles looks weird in the data (see in get_mom the note)
!pars(75)=-0.05_dp
!pars(81)=-0.01_dp
!pars(87)=0.01_dp
!pars(93)=0.03_dp


!open(unit=2,file='bestpar122818.txt',status='old',action='read') ; read(2,*) pars	; close(2)
pars(21)=-1.2_dp    !pmeet
pars(75)=pars(93)
pars(93)=pars(81)
!pmeet
pars(21)=-1.6_dp
!psih(1)
pars(16)=-1.0_dp 
pars(17)=1.0_dp
pars(18)=1.0_dp
pars(19)=-1.0_dp
!mumar
pars(75)=0.00001_dp
pars(81)=0.0005_dp
pars(87)=0.5_dp
pars(93)=1.0_dp
!ptype
pars(70:71)=0.0_dp 
pars(76:77)=0.0_dp
pars(82:83)=0.0_dp
pars(88:89)=0.0_dp

pars(29)=pars(27)

 !   pars(90)=-1.0_dp

!open(unit=2,file='bestpar010119.txt',status='old',action='read') ; read(2,*) pars2	; close(2)


!open(unit=2,file='bestpar122318_wgt1.txt',status='old',action='read') ; read(2,*) pars2	; close(2)
!pars2(20)=-1.0_dp !pkid because not iterating on pkid when onlysingles because the kid moments by age for singles looks weird in the data (see in get_mom the note)


!open(unit=2,file='bestpar011119.txt',status='old',action='read') ; read(2,*) pars2	; close(2)


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
    
    
	!call getsteps(pars,parname,stepmin)
    !stepmin=0.0 !ahu 0317
    !stepmin(2:15)=0.5 !ahu 0317
    !stepmin(22:45)=0.0_dp !ahu 041818 not iterating on the wage function parameters except for the type ones
    !open(unit=1111,file='step110218.txt',status='old',action='read') ; read(1111,*) stepmin	; close(64) 	
    !pars(13)=3.0_dp !psil1
    !pars(14)=0.0_dp !psil2
    !pars(15)=0.0_dp !psil3
    !pars(74)=-1.0_dp
    
    !pars(16)=-3.0_dp  !psi1
    !pars(17)=0.5_dp   !psi2
    !pars(18)=-3.0_dp !psi3
    !pars(19)=0.5_dp !psi4
    !pars(52)=0.2_dp
    !pars(53)=-4.0_dp
    pars(74)=-1.2   !cst
    pars(75)=1.0_dp !mumar

!**********
pars1=pars
pars1(13)=2.5_dp
pars1(14)=0.0_dp !psil2
pars1(15)=0.0_dp !psil3
pars1(1)=0.5_dp
pars1(2)=-2.0_dp
pars1(3)=0.5_dp
pars1(4)=-2.0_dp
pars1(5)=0.5_dp
pars1(6)=-1.0_dp
pars1(7)=0.5_dp
pars1(8)=-1.0_dp
pars1(42)=pars1(42)+0.2_dp
pars1(43:44)=pars1(43:44)+0.3_dp
pars1(45:47)=pars1(45:47)+0.2_dp
pars1(48)=pars1(48)+0.2_dp
pars1(49)=pars1(49)+0.1_dp
pars1(50)=pars1(50)+0.1_dp
pars1(74)=-1.6   !cst(1)
pars1(80)=pars1(74)
pars1(86)=pars1(74)
pars1(92)=pars1(74)
pars1(75)=0.5_dp !mumar
pars1(81)=1.0_dp     
pars1(87)=1.0_dp  
pars1(93)=1.0_dp  
pars1(72:73)=0.0_dp !alf1t alf2t
pars1(78:79)=-2.4_dp     
pars1(84:85)=-1.8_dp    
pars1(90:91)=-1.2_dp    
pars=pars1
!************
pars(70:71)=0.0_dp 
pars(76:77)=0.0_dp
pars(82:83)=0.0_dp
pars(88:89)=0.0_dp
pars(72:73)=0.0_dp !alf1t alf2t
pars(78)=-0.2_dp     
pars(79)=-2.0_dp     
pars(84)=-2.0_dp    
pars(85)=-0.2_dp
pars(90)=-1.2_dp    
pars(91)=-1.2_dp
pars(75)=0.3_dp !mumar
pars(81)=0.3_dp !mumar
pars(87)=0.3_dp !mumar
pars(93)=0.3_dp !mumar
!************


pars(72:73)=0.0_dp !alf1t alf2t
pars(78)=3.0_dp     
pars(79)=-2.0_dp     
pars(84)=-2.0_dp    
pars(85)=3.0_dp
pars(90)=4.0_dp    
pars(91)=4.0_dp

    pars(66)=-4.5_dp
pars(78)=3.0_dp     
pars(79)=-2.0_dp     
pars(84)=-2.0_dp    
pars(85)=3.0_dp
pars(90)=4.0_dp    
pars(91)=4.0_dp
    
    
pars=pars2
pars(42:50)=pars(42:50)-0.4_dp
pars(78)=2.0_dp     
pars(79)=-2.0_dp     
pars(84)=-2.0_dp    
pars(85)=2.0_dp
pars(90)=2.0_dp    
pars(91)=2.0_dp
pars(52)=0.0_dp

pars(70:71)=0.0_dp 
pars(76:77)=0.0_dp
pars(82:83)=0.0_dp
pars(88:89)=0.0_dp

    pars(66)=-3.0_dp

!open(unit=2,file='bestpar011219.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    !call objfunc(pars,qval) ; realpars=realpartemp
    pars(54:62)=pars(54:62)-0.5_dp
    pars(64)=-1.0_dp
    pars(52)=-1.0_dp

!open(unit=2,file='bestpar011219_2.txt',status='old',action='read') ; read(2,*) pars1	; close(2)
    pars1(64)=-1.0_dp
    pars1(52)=-1.0_dp
    pars1(16)=-3.0_dp  !psi1
    pars1(17)=0.5_dp   !psi2
    pars1(18)=0.5_dp !psi3
    pars1(19)=-3.0_dp !psi4
    pars1(66)=-2.0_dp
    pars=pars1


!open(unit=2,file='bestpar011319.txt',status='old',action='read') ; read(2,*) pars	; close(2)
pars(52)=-1.0_dp
pars(1)=0.5_dp
pars(2)=-2.0_dp
pars(3)=0.5_dp
pars(4)=-2.0_dp
pars(5)=0.5_dp
pars(6)=-1.0_dp
pars(7)=0.5_dp
pars(8)=-1.0_dp
pars(9:10)=1.5_dp
!for the second estimation I'm sending now: 
!pars(78:79)=pars(78:79)-2.0_dp
!pars(84:85)=pars(84:85)-2.0_dp
!pars(90:91)=pars(90:91)-2.0_dp
!open(unit=2,file='bestpar011419.txt',status='old',action='read') ; read(2,*) pars	; close(2)
pars(16)=-1.0_dp
pars(52)=-2.0_dp
!fore the first estimation: in migit
!pars(78:79)=-2.0_dp
!pars(84:85)=-1.0_dp
!pars(90:91)=0.0_dp
!for the second estimation I'm sending now: in migit2
pars(78:79)=0.0_dp
pars(84:85)=1.0_dp
pars(90:91)=2.0_dp

!open(unit=2,file='bp011519_4.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    pars(18)=0.0_dp             !ahu jan19 011519 getting rid of probdown
    pars(19)=0.0_dp 			!ahu jan19 011519 getting rid of probdown
    !pars1=pars
    !pars1(52)=-2.0_dp
    !pars1(2)=0.1_dp
    !pars1(17)=-5.0_dp
    !call objfunc(pars1,qval) ; realpars=realpartemp
!open(unit=2,file='bp011719.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    pars(16)=-2.0_dp
    pars(17)=-3.0_dp
    pars(66)=-0.5_dp
    pars(74)=-3.0_dp   !cst(1)
    pars(2)=-1.0_Dp
    pars(4)=-1.0_dp 
!open(unit=2,file='bp011819_2.txt',status='old',action='read') ; read(2,*) pars1	; close(2) 
!open(unit=2,file='bp011919_2.txt',status='old',action='read') ; read(2,*) pars1	; close(2)
    !pars1(16)=-20.0_dp
    !pars1(17)=-20.0_dp
    !PARS1(66)=-20.0_dp
    !pars1(13)=0.0_dp
!open(unit=2,file='bp012119_1.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    pars(1)=1.0_dp
    pars(3)=1.0_dp
    pars(17)=-20.0_dp
    pars(66)=0.0_dp !sigma needs to be high in order to match the wvar at age 18 for L moments
    pars(16)=-3.0_dp
    pars(52)=-3.00_dp
!24.07 25.57 26.99
!open(unit=2,file='bp012519.txt',status='old',action='read') ; read(2,*) pars1	; close(2)
!open(unit=2,file='bp012619.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    !mult1=20000.0_dp    
    pars(27)=-5.0_dp
    pars(51)=-100.0_dp
    !pars(33:41)=pars1(33:41)/3.0_dp
    !VERY IMPORTANT DISCOVERIES FOR RETARDED PEOPLE: 
    !pars(72)=0.5_dp
    !pars(42:50)=pars(42:50)-0.5_dp
    pars(22:23)=3000.0_dp  !when uhome is too low then all the moves happen at beginning and then nothing afterwards. changing mvecost does not help. so it has to be around 5000 currently 
    !if on average wage differences across locations is about 2K-5K then it makes sense that this is how much uhome you would need to prevent mass moves at age 17
    !because uhome gets summed up and discounted the same way as wage differences so the avg wage difference should be the same as uhome or similar in the ballpark
    pars(24)=0.0_dp    !ecst NO MORE ITERATING ON THIS SINCE NOT IDENT PROBABLY
    pars(25)=0.0_dp    !kcst NO MORE ITERATING ON THIS SINCE NOT IDENT PROBABLY    
    pars(10)=-1.0_dp
    pars(1)=-1.0_dp
    pars(2)=-3.0_dp
    pars(3)=3.0_dp
    pars(4)=-100.0_DP   !NO MORE OFLOC LAYOFF NONSENSE
    pars(8)=-100.0_DP !NO MORE OFLOC LAYOFF NONSENSE
    pars(9)=1.0_dp
    pars(33:41)=0.0_dp
    pars(42)=9.2_dp !when this is too high then that affects moving things a lot especially in beginning.the age 18:19 loc 1 wage has few people anyway
    pars(74)=0.0_dp !-200.0_dp + I*20.0_DP
    pars(66)=-3.0_dp
    pars(69)=6000.0_dp !+ 3000.0_dp*(I-1)   !sigo !0.0_dp !8000.0_dp !1000.0_dp + 700.0_dp*I   !sigo
    pars(24)=0.0_dp !0.0_dp - 3000.0_dp*(j-1)
    pars(13)=2.5_dp
    pars(69)=6000.0_dp
    !do i=1,3
    !do j=1,3
    !pars(13)=2.5_dp+0.5_dp*(j-1)
    !pars(69)=6000.0_dp+2000.0_dp*(i-1)
    !pars(74)=0.0_dp-2000.0_dp*(j-1)
    
    
    pars(33)=-2000.0_dp !changed from -9K since you don't need it that negative large when pars42 is no longer that high
    pars(35)=-1000.0_dp
    pars(36)=-300.0_dp
    pars(37)=2000.0_dp
    pars(38)=2000.0_dp
    pars(39)=-1000.0_dp
    pars(40)=2000.0_dp
    pars(41)=-2000.0_dp
    pars(51)=-0.8_dp    
    pars(16)=0.1_dp
    pars(52)=-2.8_dp
    
    !open(unit=2,file='bp031219.txt',status='old',action='read') ; read(2,*) pars1	; close(2)
    !open(unit=2,file='bp031319.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    pars(42)=0.0_dp
    pars(43)=-0.1_dp
    pars(44)=0.0_dp
    pars(45:49)=-0.1_dp
    pars(50)=0.1_dp
    pars(72:73)=8.75_dp !9.00_dp !8.75_dp
    !pars(76:77)=-1.0_dp
    pars(78:79)=9.9_dp !9.40_dp !9.9_dp
    open(unit=2,file='bp040219_2.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    pars(22)=5000.0_dp
    pars(69)=7000.0_dp
    pars(13)=2.5_dp
    pars(1)=pars(1)+0.2_dp
    pars(2)=-2.0_dp    
    !NOW FEMALES 
    pars(5:8)=pars(1:4)
    pars(11:12)=pars(9:10)
    pars(23)=pars(22)
    pars(54:65)=pars(42:53)
    pars(67)=pars(66)
    pars(73)=pars(72)
    pars(79)=pars(78)
    !pars(28)=pars(27)
    !pars(32)=pars(31)
    !open(unit=2,file='bp040419_1.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    pars(28)=pars(28)-1.0_dp
    pars(32)=pars(32)-1.0_dp
    !pars(23)=7000.0_dp 
    !open(unit=2,file='bp040419_2.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    pars(11)=-1.0_dp    !in order to bring up eu stay
    pars(69)=10800.0_dp !in order to bring down eu move
    pars(6)=-1.3_dp     !in order to bring down ee stay

    !COMPARE MALES BEST PARAMS TO FEMALES BEST PARAMS
    open(unit=2,file='bp040219_2.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    pars(22)=5000.0_dp
    pars(69)=7000.0_dp
    pars(13)=2.5_dp
    pars(1)=pars(1)+0.2_dp
    pars(2)=-2.0_dp    
    !NOW FEMALES 
    pars(5:8)=pars(1:4)
    pars(11:12)=pars(9:10)
    pars(23)=pars(22)
    pars(54:65)=pars(42:53)
    pars(67)=pars(66)
    pars(73)=pars(72)
    pars(79)=pars(78)
    !AG 0CT 9, 2019
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp  
    pars(28)=pars(27)
    pars(32)=pars(31)   !ok I DIDN'T have alphakid's equal to each othe rin that estimation where I first added females. but that was intentional. 
    pars(68)=pars(69) !sigo trying to see whether fem will look same as male
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp  

    !open(unit=2,file='bp040519_1.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    !CHANGES
    pars(15)=0.66_dp !this is ro. It was par(68) before. 
    pars(22)=16476.0_dp
    pars(68)=10475.0_dp
    !open(unit=2,file='bp040519_2.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    !pars(14)=-1000.0_dp
    !open(unit=2,file='bp040519_3.txt',status='old',action='read') ; read(2,*) pars1	; close(2)
    !pars1(14)=-1000.0_dp  
    !open(unit=2,file='bp040819_2.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    pars(2)=-0.11_dp       
    pars(6)=-1.55_dp     
    pars(14)=0.0_dp
    pars(22)=6298.69_dp        
    pars(23)=10220.54_dp       
    pars(68)=1977.5_dp       
    pars(69)=5914.89_dp  
    !open(unit=2,file='bp041019_1.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    pars(13)=pars(13)+0.4_dp
    pars(24)=-500.0_dp
    pars(25)=-1500.0_dp
    pars(66:67)=pars(66:67)+1.0_dp
    pars(7)=0.5_dp
    
    !pars(84:85)=pars(72:73)
    !pars(90:91)=pars(78:79)
    !open(unit=2,file='bp041019_2.txt',status='old',action='read') ; read(2,*) pars1	; close(2)

    !open(unit=2,file='bp041019_3.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    pars(33:41)=pars(33:41)*2.0_dp
    pars(22:23)=pars(22:23)*1.5_dp
    pars(7)=-2.0_dp
    pars(66:67)=pars(66:67)-1.0_dp

    !open(unit=2,file='bp041119_1.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    !pars(74)=-1000.0_dp
    !pars(80)=-1000.0_dp
    !pars(86)=-1000.0_dp
    !pars(92)=-1000.0_dp
    
    !open(unit=2,file='bp041019_3.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    pars(7)=-2.0_dp
    pars(74)=-1200.0_dp
    pars(80)=-1200.0_dp
    pars(86)=-1200.0_dp
    pars(92)=-1200.0_dp
    pars(13)=-2.08_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp  
    
    !ahu 101119
     open(unit=2,file='bp040219_2.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp  
    pars(22)=5000.0_dp
    pars(69)=7000.0_dp
    pars(13)=2.5_dp
    pars(1)=pars(1)+0.2_dp
    pars(2)=-2.0_dp    
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp  
   
    
     open(unit=2,file='bp041019_3.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp  
    !pars(75)=-2.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp  
    !pars(75)=-1.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !1_mumar3
    !pars(75)=0.02_dp !ahu030622
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp   
    !pars(75)=0.29_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=2.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      

    !030822_1q010q201 
   	!multmar=100,000
    !pars(75)=0.002_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.5_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=25.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      

    !030822_2mardecrease (now that we got decreasing mar rates. will see chk2)
   	!multmar=100,000
    !pars(75)=0.002_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.5_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=25.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      

    !030822_3
    !pars(75)=0.5_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=25.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    
    !031122_1
    !pars(75)=0.00005_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.0001_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.0005_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.001_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.005_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.01_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.05_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.1_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.5_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=1.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=3.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      


	!031122_2
    !pars(75)=0.005_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.006_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.007_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.008_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.009_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.01_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.012_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.015_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
   	!pars(75)=0.02_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.025_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
   	!pars(75)=0.03_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.035_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.04_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.045_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.05_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      


	!031122_3
    !pars(75)=0.04_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.041_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.042_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.043_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.044_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.45_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.046_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.047_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
   	!pars(75)=0.048_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.049_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
   	!pars(75)=0.05_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
	!mu_mar            9998.67  10248.56  10498.46  10748.34  10998.23 110639.23  11497.97  11747.84  11997.70  12247.55  12497.40
	!In the above run, where I was trying to figure out where the div rate starting to increase (instead of decrease) occurs between mumar values changing, I all of a sudden see the below horrific view where there is indeed an increse at sme point between mumar 99998 and mumar 12497 but it is not monotonic and there is a big huge jump when mumar changes from 10998 to 110639

   !nonlabinc=0, onlysingles=false,mumar changing, skriv and yaz in getdec_c see ahu030622 or ahu 030622
 
    !031522 MARCH1522 corner solution stuff added 
    
    
    !m031722_1
    !pars(75)=0.035_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.04_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.045_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.05_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      


    !m031722_2  
    !pars(75)=0.00005_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.0001_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.001_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.01_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      

    !m031722_3
    !pars(75)=0.00005_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.0001_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.001_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.005_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp       
    !pars(75)=0.01_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.015_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.02_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      


    !m031722_4
 !   pars(75)=0.001_dp
 !   call getpars(pars,realpars)
 !   call objfunc(pars,qval) ; realpars=realpartemp      
 !   pars(75)=0.0015_dp
 !   call getpars(pars,realpars)
 !   call objfunc(pars,qval) ; realpars=realpartemp       
 !   pars(75)=0.002_dp
 !   call getpars(pars,realpars)
 !   call objfunc(pars,qval) ; realpars=realpartemp      
 !   pars(75)=0.0025_dp
 !   call getpars(pars,realpars)
 !   call objfunc(pars,qval) ; realpars=realpartemp      
 !   pars(75)=0.003_dp
 !   call getpars(pars,realpars)
 !   call objfunc(pars,qval) ; realpars=realpartemp      
 !   pars(75)=0.0035_dp
 !   call getpars(pars,realpars)
 !   call objfunc(pars,qval) ; realpars=realpartemp      
 !   pars(75)=0.004_dp
 !   call getpars(pars,realpars)
 !   call objfunc(pars,qval) ; realpars=realpartemp      
 !   pars(75)=0.0045_dp
 !   call getpars(pars,realpars)
 !   call objfunc(pars,qval) ; realpars=realpartemp      
 !   pars(75)=0.005_dp
 !   call getpars(pars,realpars)
 !   call objfunc(pars,qval) ; realpars=realpartemp      

   
    !!m031722_5
    !pars(75)=0.001_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.00105_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp       
    !pars(75)=0.0011_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.00115_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.0012_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.00125_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.0013_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.00135_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.0014_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.00145_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !pars(75)=0.0015_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    
    
    !march 23 2022
    !M032322_1 
    !pars(75)=0.0015_dp
    !nonlabinc=0.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp          
    !nonlabinc=0.5_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !nonlabinc=5.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp  
    !nonlabinc=10.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !nonlabinc=100.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !nonlabinc=500.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !nonlabinc=1000.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !nonlabinc=2000.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !nonlabinc=3000.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !nonlabinc=4000.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !nonlabinc=5000.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      
    !nonlabinc=6000.0_dp
    !call getpars(pars,realpars)
    !call objfunc(pars,qval) ; realpars=realpartemp      


    !AUGUST 28, 2022
     !pars(75)=0.0015_dp !THIS IS INHERITED FROM BEFORE. JUST RUNNING TO SEE HOW MIGSOL082822 IS DOING
     !call getpars(pars,realpars)
     !call objfunc(pars,qval) ; realpars=realpartemp  
   
    !090422 SEPTEMBER 4 2022
     pars(75)=0.0015_dp
     nonlabinc=0.0_dp
     !call getpars(pars,realpars)
     !call objfunc(pars,qval) ; realpars=realpartemp     
     !pars(75)=0.0015_dp
     !nonlabinc=2000.0_dp
     !call getpars(pars,realpars)
     !call objfunc(pars,qval) ; realpars=realpartemp      
     !pars(75)=0.0015_dp
     !nonlabinc=5000.0_dp
     !call getpars(pars,realpars)
     !call objfunc(pars,qval) ; realpars=realpartemp      
     !pars(75)=0.0015_dp
     !nonlabinc=10000.0_dp
     !call getpars(pars,realpars)
     !call objfunc(pars,qval) ; realpars=realpartemp      

     !pars(75)=0.02_dp
     !nonlabinc=0.0_dp
     !call getpars(pars,realpars)
     !call objfunc(pars,qval) ; realpars=realpartemp          
     !pars(75)=0.02_dp
     !nonlabinc=2000.0_dp
     !call getpars(pars,realpars)
     !call objfunc(pars,qval) ; realpars=realpartemp      
     !pars(75)=0.02_dp
     !nonlabinc=5000.0_dp
     !call getpars(pars,realpars)
     !call objfunc(pars,qval) ; realpars=realpartemp      
     !pars(75)=0.02_dp
     !nonlabinc=10000.0_dp
     !call getpars(pars,realpars)
     !call objfunc(pars,qval) ; realpars=realpartemp     

     !open(unit=2,file='bp090122final.txt',status='old',action='read') ; read(2,*) pars	; close(2)
     !call getpars(pars,realpars)
     !call objfunc(pars,qval) ; realpars=realpartemp     
      
    open(unit=2,file='bp041019_3.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    nonlabinc=0.0_dp
    !The below are in order to set type parameters to values they were before we changed the getpars setup 
    !(where ptype's are no longer 0 and mumar's are no longer all being set equal to mumar(1))
    pars(75)=0.0015_dp !mumar1
    pars(76)=0.0_dp !ptypehs2 
    pars(77)=0.0_dp !ptypecol2
    pars(81)=0.0015_dp !mumar2
    pars(82)=0.0_dp !ptypehs3 
    pars(83)=0.0_dp !ptypecol3
    pars(87)=0.0015_dp !mumar3
    pars(88)=0.0_dp !ptypehs4 
    pars(89)=0.0_dp !ptypecol4
    pars(93)=0.0015_dp !mumar4
    call getpars(pars,realpars)
    call objfunc(pars,qval) ; realpars=realpartemp     

 
    open(unit=2,file='bp041019_3.txt',status='old',action='read') ; read(2,*) pars	; close(2)
    nonlabinc=0.0_dp
    !changing only type parameters in pars
    pars(75)=0.0015_dp !mumar1 !should this be set to 0? 
    pars(76)=0.0_dp !ptypehs2 
    pars(77)=0.0_dp !ptypecol2
    pars(80)=1000.0_dp !cst2 
    pars(81)=0.0015_dp !mumar2    
    pars(82)=0.0_dp !ptypehs3 
    pars(83)=0.0_dp !ptypecol3
    pars(86)=2000.0_dp !cst3
    pars(87)=0.0015_dp !mumar3
    pars(88)=0.0_dp !ptypehs4 
    pars(89)=0.0_dp !ptypecol4
    pars(92)=3000.0_dp !cst4 
    pars(93)=0.0015_dp !mumar4
    call getpars(pars,realpars)
    call objfunc(pars,qval) ; realpars=realpartemp     


    !now run with new pars but again change the type parameters
    !open(unit=2,file='bp090922therm.txt',status='old',action='read') ; read(2,*) pars1	; close(2)
    open(unit=2,file='bp090922.txt',status='old',action='read') ; read(2,*) pars1	; close(2)
    nonlabinc=0.0_dp
    pars=pars1
    !pars(70)=ptypehs1  !set to 0 in getpars and kept fixed
    !pars(71)=ptypecol1 !set to 0 in getpars and kept fixed
    !pars(72)=alf1t1 
    !pars(73)=alf1t1
    !pars(74)=cst1  !is 0 (set to 0 in getpars) and should be kept fixed
    pars(75)=0.0015_dp !mumar1 !should this be set to 0? 
    
    pars(76)=0.0_dp !ptypehs2 
    pars(77)=0.0_dp !ptypecol2
    !pars(78)=alf1t2 
    !pars(79)=alf1t2
    pars(80)=1000.0_dp !cst2 
    pars(81)=0.0015_dp !mumar2
    
    pars(82)=0.0_dp !ptypehs3 
    pars(83)=0.0_dp !ptypecol3
    !pars(84)=alf1t 3
    !pars(85)=alf1t 3
    pars(86)=2000.0_dp !cst3
    pars(87)=0.0015_dp !mumar3

    pars(88)=0.0_dp !ptypehs4 
    pars(89)=0.0_dp !ptypecol4
    !pars(90)=alf1t 4
    !pars(91)=alf1t 4
    pars(92)=3000.0_dp !cst4 
    pars(93)=0.0015_dp !mumar4
    call getpars(pars,realpars)
    call objfunc(pars,qval) ; realpars=realpartemp     
    

    stepmin=stepos !ahu 121118
    tstart=qval !ag091122 agsept2022 changing thermsimp value !0.0_dp !*qval !10.0_dp !
    tstep=0.8_dp  !ag091122 agsept2022 changing tstep value !0.0_dp !0.7_dp !0.8_dp
    tfreq=2*COUNT(stepmin /= zero) !ag091122 agsept2022 changing tfreq value
    saseed=1
    !call objfunc(pars,qval) ; realpars=realpartemp  

    
  
    if (iam==0) print*, 'Here is qval,tstart,tstep,tfreq', qval,tstart,tstep,tfreq
    
    
    !do i=1,3
    !do j=1,3
    !pars(7)=pars(7)-1.0_dp
    !pars(11)=pars(11)+0.5_dp
    !pars(13)=2.5_dp-0.4_dp*i
    !pars(69)=4800.0_dp+3000.0_dp*j
    !pars(6)=pars(6)+0.5_dp
    !pars1(22)=pars1(22)+3500.0_dp
    !call objfunc(pars1,qval) ; realpars=realpartemp 
    !end do 
    !end do 
    
    !do i=1,3
    !pars1(68)=pars1(68)+2500.0_dp
    !if (iam==0) print*, 'heres',i,pars1(22),pars1(68)
    !call objfunc(pars1,qval) ; realpars=realpartemp 
    !end do 
    

    
    !pars(69)=0.0_dp 
    !call objfunc(pars,qval) ; realpars=realpartemp    !this happens for just pars33's without the other changes too. so loc 1 is special. 

    !do i=1,3
    !do j=1,3
    !pars(22)=5000.0_dp+3000.0_dp*(i-1)
    !pars(69)=7000.0_dp+6000.0_dp*(j-1)
    !pars(13)=2.5_dp-0.3_dp*(j-1)
    !call objfunc(pars,qval) ; realpars=realpartemp    !this happens for just pars33's without the other changes too. so loc 1 is special. 
    !end do 
    !end do 
    

    !pars(2)=-2.0_dp    
    !call objfunc(pars,qval) ; realpars=realpartemp    !this happens for just pars33's without the other changes too. so loc 1 is special. 
    !pars(2)=-1.5_dp    
    !call objfunc(pars,qval) ; realpars=realpartemp    !this happens for just pars33's without the other changes too. so loc 1 is special. 

    !pars(1)=pars(1)+0.2_dp
    !pars(2)=-2.0_dp    
    !call objfunc(pars,qval) ; realpars=realpartemp    !this happens for just pars33's without the other changes too. so loc 1 is special. 
    !pars(1)=pars(1)+0.2_dp
    !pars(2)=-1.5_dp    
    !call objfunc(pars,qval) ; realpars=realpartemp    !this happens for just pars33's without the other changes too. so loc 1 is special. 


    !do j=1,2
    !pars(69)=6000.0_dp-2000.0_dp*j
    !call objfunc(pars,qval) ; realpars=realpartemp    
    !end do

    
    !pars(22)=5000.0_dp
    !pars(13)=3.5_dp
    !call objfunc(pars,qval) ; realpars=realpartemp    !this happens for just pars33's without the other changes too. so loc 1 is special. 
    
    

    !pars=pars1
    !pars(33)=pars1(33)/2.0_dp
    !pars(69)=3000.0_dp
    !call objfunc(pars,qval) ; realpars=realpartemp    this happens for just pars33's without the other changes too. so loc 1 is special. 
    !pars(69)=6000.0_dp
    !call objfunc(pars,qval) ; realpars=realpartemp    
    !pars=pars1
    !pars(33)=pars1(33)
    !pars(69)=3000.0_dp
    !call objfunc(pars,qval) ; realpars=realpartemp    
    !pars(69)=6000.0_dp
    !call objfunc(pars,qval) ; realpars=realpartemp    
    !pars=pars1
    !pars(33)=pars1(33)*2.0_dp
    !pars(69)=3000.0_dp
    !call objfunc(pars,qval) ; realpars=realpartemp    
    !pars(69)=6000.0_dp
    !call objfunc(pars,qval) ; realpars=realpartemp    

    
    
    !pars(33:41)=pars1(33:41)/10.0_dp
    !pars(69)=0.0_dp
    !call objfunc(pars,qval) ; realpars=realpartemp    
    !pars1(69)=3000.0_dp
    !call objfunc(pars1,qval) ; realpars=realpartemp    
    !pars1(69)=6000.0_dp
    !call objfunc(pars1,qval) ; realpars=realpartemp    

    
    
    !pars1(22)=1500.0_dp
    !pars1(69)=3000.0_dp
    !call objfunc(pars1,qval) ; realpars=realpartemp   
    !pars1(41)=-1000.0_dp
    !pars1(38)=3000.0_dp
    !pars1(33:37)=1000.0_dp
    !pars1(38:41)=-1000.0_dp
    !pars1(33:41)=0.0_dp
    !pars1(22)=5000.0_dp
    !pars1(69)=6000.0_dp
    !call objfunc(pars1,qval) ; realpars=realpartemp   
    !pars1(33:37)=10000.0_dp
    !pars1(38:41)=-10000.0_dp
    !pars1(69)=10000.0_dp
    !call objfunc(pars1,qval) ; realpars=realpartemp   
    !pars1(69)=15000.0_dp
    !call objfunc(pars1,qval) ; realpars=realpartemp   
    !pars1(69)=30000.0_dp
    !call objfunc(pars1,qval) ; realpars=realpartemp   

    !pars1(33:41)=0.0_dp
    !pars1(22)=-10.0_dp
    !pars1(69)=0.5_dp
    !call objfunc(pars1,qval) ; realpars=realpartemp   
    !pars1(33:41)=0.0_dp
    !pars1(22)=-3.0_dp
    !pars1(69)=0.5_dp
    !call objfunc(pars1,qval) ; realpars=realpartemp   
    !pars1(33:41)=0.0_dp
    !pars1(22)=-2.0_dp
    !pars1(69)=0.5_dp
    !call objfunc(pars1,qval) ; realpars=realpartemp   
    !pars1(33:41)=0.0_dp
    !pars1(22)=-1.0_dp
    !pars1(69)=0.5_dp
    !call objfunc(pars1,qval) ; realpars=realpartemp   

    

    
	if (optimize) then 
		! call simplex: the call to minim (the simplex routine) is set up to work with starting vector parvector
		!if (iam==0) then ; iwriteminim=1 ; else ; iwriteminim=-1 ; end if 
		!open(unit=64,file='step11_10.txt',status='old',action='read') ; read(64,*) stepmin	; close(64) 			
		!open(unit=64,file='step111416.txt',status='old',action='read') ; read(64,*) stepmin	; close(64) 			
        !open(unit=64,file='stepp.txt',status='old',action='read') ; read(64,*) stepmin	; close(64) 
        !stepmin=0.5_dp

        
        iwriteminim=30
		if (groups) then 
			call pminim(pars, stepmin, npars, qval, maxfcn, iwriteminim, stopcr, nloop, iquad, simp, var, objfunc, writebest,ifault,mygroup,numgroup,tstart, tstep	, tfreq, saseed)
		else 
			call pminim(pars, stepmin, npars, qval, maxfcn, iwriteminim, stopcr, nloop, iquad, simp, var, objfunc, writebest,ifault,iam,numworld,tstart, tstep	, tfreq, saseed)
		end if 
		if (iam==0) then ; print*, "out of minim now and here is ifault ", ifault ; end if 	
	else 
		!conditional_moments=.false.		
		conditional_moments=.true.	       

        !pars(6)=-0.1_dp
        !pars(68)=0.9_dp
        !pars(69)=1.3_dp
        !pars(78)=-100.0_dp
        !nonneg=.FALSE.
        !nonlabinc=(/ 300.0_dp,1100.0_dp /) 
        
        !pars(7)=-4.0_dp      !alphaed(m,noed)
        !pars(9)=-4.0_dp      !alphaed(m,ed)
        !call objfunc(pars,q) ; realpars=realpartemp        
        !pars(2)=0.8_dp
        !call objfunc(pars,q) ; realpars=realpartemp

        !sig_o=4000.0_dp
        !call objfunc(pars,qval) ; realpars=realpartemp
        !sig_o=7000.0_dp
        !call objfunc(pars,qval) ; realpars=realpartemp
        

        !pars(24)=1.5_dp
        !call objfunc(pars,q) ; realpars=realpartemp
        !pars1=pars
        !pars1(1)=1.0_dp
        !pars1(2)=-1.0_dp
        !call objfunc(pars1,q) ; realpars=realpartemp

        !nonneg=.FALSE.
        !nonlabinc=(/ 0.0_dp,0.0_dp /) 

        !mom11418_2_pmeet14:
        !pars(12)=-1.8_dp     !pmeet
        !pars(4)=0._dp     !cst(1)
        !pars(5)=0._dp     !kcst
        !pars(75)=0._dp    !cst(ntypp) 
        !pars(76)=0._dp    !ecst
        !nonneg=.TRUE.
        !nonlabinc=(/ 300.0_dp,1100.0_dp /) 
        !pars(68:69)=0.001_dp
        !call objfunc(pars,q) ; realpars=realpartemp
        !pars(68:69)=0.01_dp
        !call objfunc(pars,q) ; realpars=realpartemp
        !pars(68:69)=0.05_dp
        !call objfunc(pars,q) ; realpars=realpartemp
        !pars(68:69)=0.1_dp
        !call objfunc(pars,q) ; realpars=realpartemp
        !pars(68:69)=0.15_dp
        !call objfunc(pars,q) ; realpars=realpartemp
        !pars(68:69)=0.2_dp
        !call objfunc(pars,q) ; realpars=realpartemp

        !mom11418_3_pmeet14_alphaedhigh:
        !pars(4)=0._dp     !cst(1)
        !pars(5)=0._dp     !kcst
        !pars(75)=0._dp    !cst(ntypp) 
        !pars(76)=0._dp    !ecst
        !pars(68:69)=0.001_dp
        !call objfunc(pars,q) ; realpars=realpartemp
        !pars(68:69)=0.01_dp
        !call objfunc(pars,q) ; realpars=realpartemp
        !pars(68:69)=0.05_dp
        !call objfunc(pars,q) ; realpars=realpartemp
        !pars(68:69)=0.1_dp
        !call objfunc(pars,q) ; realpars=realpartemp
        !pars(68:69)=0.15_dp
        !call objfunc(pars,q) ; realpars=realpartemp

        !if (iam==0) then 
		!    open(unit=2225, file='chksimplex.txt',status='replace')
        !end if 
        !do i=1,npars
        !    pars1=pars
        !    pars1(i)=pars(i)+stepmin(i)
        !    call objfunc(pars1,val) ; realpars1=realpartemp
        !    if (iam==0) then 
        !        write(2225,*)
        !        write(2225,*)
        !        write(2225,'(tr3,"i",tr3,"j",tr2,"iter",2x,tr20,tr13,"q",tr9,"q1(j)",tr6,"pars",tr5,"pars1",5x,tr6,"realpars",5x,tr5,"realpars1",5x,tr3,"step(j)",5x,tr2,"abs(q1(j)-q)")')
        !        write(2225,'(2i4,i6,". ",1a20,2f14.2,2F10.2,2(5x,f14.2),(5x,f10.2),(5x,f14.2) )') i,0,iter-1,parname(i),q,val,pars(i),pars1(i),realpars(i),realpars1(i),stepmin(i),abs(val-q)			
        !    end if
        !end do
        !close(2225)
        
        

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
        
		if (chkstep) then
            lb=0.4_sp ; ub=1.25_sp ; incr=0.1_dp
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
!		if (groups) then
!			do i=0,ninp-1
!				ranks(i)=i
!			enddo
!			call mpi_comm_group(mpi_comm_world,mpi_group_world,mpierr)
!			call mpi_group_incl(mpi_group_world,ninp,ranks,group,mpierr)
!			call mpi_comm_create(mpi_comm_world,group,comm,mpierr)
 !       end if 

		!splitkey=mod(myid,nprocs/ninp)
		!call mpi_comm_split(mpi_comm_world,splitkey,myid,comm,mpierr)

 	!call mpi_bcast(stepmin,npars,mpi_real8,0,mpi_comm_world,mpierr)
	!if (iam==15) then 
	!	print*, "here is stepmin(6:8) "  , mygroup,iam,stepmin(6:8)
	!end if 
	
	!print*, "here everyone is done now ", mygroup,iam 
	!print*, 1.0_dp/(1.0_dp+exp(pars(64))+exp(pars(65))),exp(pars(64))/(1.0_dp+exp(pars(64))+exp(pars(65))),exp(pars(65))/(1.0_dp+exp(pars(64))+exp(pars(65)))
	!print*, 1.0_dp/(1.0_dp+exp(pars(66))+exp(pars(67))),exp(pars(66))/(1.0_dp+exp(pars(66))+exp(pars(67))),exp(pars(67))/(1.0_dp+exp(pars(66))+exp(pars(67)))	

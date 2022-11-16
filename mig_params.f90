	module params
	!insert space and keep tabs option in visual studio
    ! testing
	use nrtype 
	use nr 
	use alib, only: one,logit,logitinv,min2pls,min2plsinv,lin2ndim,ndim2lin,sqrtpi,multinom,condmom
	implicit none 
    include 'mpif.h'
    integer :: mysay,mygroup
	integer :: iter,comm,iwritegen 
	!real(dp) :: one=1.0_dp
	!integer(i4b), parameter :: rp = kind(1.0d0)			! kind(1.0) !!!
    integer(i4b) :: policytax !set in main to determine whether running policy experiment and which one
    real(dp) :: moveshockdiv !assigned in main 
    real(dp), parameter :: replacement_rate=0.4_dp          !ahu summer18 050318: added replacement rate
    integer(i4b), parameter :: nl=9,ndecile=10
    !ahu030622	logical, parameter :: groups=.true.,onlysingles=.true.,onlymales=.false.,onlyfem=.false.,optimize=.true.,chkstep=.false.,condmomcompare=.false.,comparepars=.false.,extramoments=.true.
    integer(i4b), parameter :: numit=2
    logical, parameter :: groups=.true.,onlysingles=.false.,onlymales=.false.,onlyfem=.false.
    logical, parameter :: optimize=.false.,chkstep=.false.,chkobj=.true.,condmomcompare=.false.,comparepars=.false.
    logical, parameter :: typemoments=.true.
    logical :: nonneg,terminalval
    logical :: onthejobsearch=.TRUE. !set in main
    real(dp), dimension(2) :: nonlabinc !=(/ 0.0_dp,0.0_dp /) !(/ 300.0_dp,1100.0_dp /) !ahu summer18 051418: changing it back to parameter and changing dimension to 2 (not educ and educ) !ahu summer18 042318 changing this so it is set at main again
	real(dp), parameter :: eps = 1.0d-6,zero=0.0_dp,epstest=2.0_dp					! could do tiny(.) but that gives a number that is way too small and therefore not appropriate for those places where we check the inequalities or equalities	
	real(dp), parameter :: eps2= 1.0d-6
	integer(i4b), parameter :: nhome=1,nhomep=nl
	logical :: conditional_moments		! can set this in main
	logical :: skriv,yaz,insol,yazmax   
	character(len=1), parameter :: runid='r'		! string inserted into output filenames to identify which run !ahu 062413 set this in main instead 
	integer(i4b), parameter :: nco=1,ncop=1
	integer(i4b), parameter :: ntyp=1,ntypp=4   ! types !ahu030622 changed ntypp to 1 (was 4)
	integer(i4b), parameter :: nin  = nco * ntyp * nhome
	integer(i4b), parameter :: ninp = ncop * ntypp * nhomep
    integer(i4b) :: nindex !determined in objf according to groups, it's either nin or ninp
	integer(i4b) :: iwritemom,myhome,mytyp,myco,myindex,mygrank
	logical, parameter :: indsimplexwrite=.false.		! in optimiation with parallel simplex, re-solve model and write moments to file whenever check for,find new best point 
	logical, parameter :: write_welfare=.false.			! write things needed for welfare calculations
	logical, parameter :: grid_usage=.false.			! keep track of how much of each grid is being used in simulations
	logical, parameter :: icheck_eqvmvf=.false.,icheck_eqvcvs=.false.,icheck_probs=.false.
	integer(i4b), parameter :: npars    = 93
    character(len=15), dimension(npars) :: parname ! names of each parameter   !ahu 121118 now declkare as global it here instead of getsteps
    character(len=25), dimension(nl) :: locname
    real(dp), dimension(npars) :: stepmin,realpartemp,parsforcheck,stepos !ahu 121118
	!character(len=15), dimension(npars) :: parname 
	integer(i4b), parameter :: mna=18,mxa=50   !,mxai=50		!ahu 070713 40	!if you ever change mina from 16 that might be problematic, because of psiddata%home and simdata%home definitions. look in read_data and read simdata for this
    integer(i4b), parameter :: MNAD=MNA-1,MXAD=MXA-1            !ahu jan19 010219
    integer(i4b), parameter :: nh=2,nexp=2,nsimeach=10,neduc=2,nkid=2 !kid is 1 if no kid,2 if yes kid !ahu 0327 changed nsimeach from 10 to 5
	integer(i4b), parameter :: np=3,nz=1 !ag090122 agsept2022 changed nz frmo 1 to 5 !ahu 121818 changed from 3 to 6 !ahu 0327 changed np from 5 to 2
	integer(i4b), parameter :: nqs = (np+2) *  nl
	integer(i4b), parameter :: nq  = (np+2) * (np+2) * nl
	integer(i4b), parameter :: np1=np+1, np2=np+2   !w=np1 is getting laid off in the shock space q and w=np2 is nothing happening. In the state space q, np1 is unemployment and np2 is not in the state space q (I try to check that it is not there, at various points in code)
	integer(i4b), parameter :: nxs = neduc * nexp * nkid
	integer(i4b), parameter :: nx  = nxs * nxs
	integer(i4b), parameter :: ncs = 3 !ahu october2022 HUGE MAJOR CHANGE HUGE MAJOR CHANGE nl+2
	integer(i4b), parameter :: nc  = 9   !nl+8 !ahu october2022 HUGE MAJOR CHANGE HUGE MAJOR CHANGE
    integer(i4b), parameter :: nepsmove=3, nepskid=2 !ag090122 agsept2022 changed nepsmove frmo 2 to 5 !ahumarch1022 changed nepsmove to 2 from 13
	integer(i4b) :: numperdat,numperobsdat,numpersim !previously ndata,ndataobs,nsim ALL set in main now 
	!integer(i4b), parameter :: ndataobs = 84507  set in main now 
	!integer(i4b), parameter :: nsim     = ndata*nsimeach  set in main now
	integer(i4b), parameter :: nmom     = 2000 !ahu summer18 050418: changed from 4200 to 498
    integer(i4b) :: calcvar(nmom),calcorr(nmom)
	integer(i4b), parameter :: maxrellength=10
	integer(i4b), parameter :: namelen=90					!if you change this, don't forget to also change a100 in writemoments	
	integer(i4b), parameter :: ma=1,fe=2
    INTEGER(I4B), PARAMETER :: NOCOLLEGE=1,COLLEGE=2
	integer(i4b), parameter, dimension(2) :: agestart=(/ mna,22 /)		!changed this from 18,22 !chanage this back ahu 070312 (/18,22/) !starting age for simulations for each education level
	real(dp), parameter :: mult1=10000.0_dp !ahu jan19 012519
    real(dp), parameter :: multmar=50000.0_dp,multsigo=300000.0_dp,multdiv=5000.0_dp,multcst=30000.0_dp !ahu jan19 012019  !ahu030622 VERY IMPORTANT CHANGE MULTMAR
	real(dp), parameter :: maxhgrid=8.0_dp 
	real(dp), parameter :: tottime=16.0_dp
	real(dp), parameter :: hhours_conv=250.0_dp					! multuiply hours per day by this to get hours per year for hours worked
	real(dp), parameter :: maxh=hhours_conv*maxhgrid				! hhours_conv*hmgrid(ntwork) ! truncation point of labor market hours in data !ahu 071712
	real(dp), parameter, dimension(nh) :: hgrid=(/ 0.0_dp,maxhgrid /) 
	real(dp), parameter, dimension(nh) :: hgrid_ann=hgrid*hhours_conv
	real(dp), parameter :: d1=1.0_dp						! note that now doing calculations as single real(sp) not double real (i modified module1 to allow this)
	real(dp), parameter :: h_parttime=1000.0_dp					! number of hours to be considred part time !ahu 071712
	real(dp), parameter :: h_fulltime=1000.0_dp					! ahu 071712 changed from 1000 too 2000 !ahu 062812 changed this to 1000. 2000. ! number of hours to be considred full time 
        !ahu 021817: note that the employment is decided according to whether hours is more than h_fulltime which is 1000 but annual wages is calculated using wge*h_wmult where h_wmult is 2000. 
        !This is how you did it in the original version. I do the same thing so that the wage numbers are consistent with the previous version. See page 15 last paragraph in the original text.	
    real(dp), parameter :: h_wmult=2000.0_dp                    !what you multiply hourly wages with in order to turn them into annual wages
    real(dp), parameter :: hbar=h_parttime
	real(dp), parameter :: minw=1.0_dp					! lower truncation point of male log wage
	real(dp), parameter :: maxw=150.0_dp                ! upper truncation point of male log wage
	real(dp), parameter :: pen=-99999999.0_dp
	integer(i4b), parameter :: ipen=-99999	
    real(dp), parameter :: wtrans=10.0_dp,wwaged=1.0_dp,wdifww=1.0_dp,wrel=1.0_dp,wmove=1.0_dp,whour=1.0_dp,wwvar=10.0_dp
    real(dp), parameter :: wwage=1.0_dp,wkid=1.0_dp,wmovemar=1.0_dp,wmovesin=1.0_dp,wwagebymove=1.0_dp		!ahu 121918 changed wmove to 10 from 1 and changed wmovemar from 10 to 100		! weights for moments for married couples. set in objfunc.
!    real(dp), parameter :: wtrans=100.0_dp,wwaged=10.0_dp,wdifww=100.0_dp,wrel=1.0_dp,wmove=10.0_dp,whour=1.0_dp,wwvar=100.0_dp
!    real(dp), parameter :: wwage=1.0_dp,wkid=1.0_dp,wmovemar=1.0_dp,wmovesin=1.0_dp,wwagebymove=1.0_dp		!ahu 121918 changed wmove to 10 from 1 and changed wmovemar from 10 to 100		! weights for moments for married couples. set in objfunc.
    character(len=23) :: datafilename  != 'familymigpsid.txt' ! data filename set in main now
	!character(len=23), parameter :: initcondfile= 'familymiginit111113.txt' ! filename of initial conditions
	character(len=23) :: momentfile='mom.txt'
    character(len=23) :: momentonlyfile='momonly.txt'
    integer(i4b) :: mominfo(0:5,nmom)
	!model parameters 
	!real(sp):: umdist_condf(np,np),ufdist_condm(np,np)
	real(dp), parameter :: delta=0.96_dp,alf=0.5_dp   !ahumarch1022 delta=0.96_dp changing to 0 to figure out the mumardecrease problem
	real(dp), parameter :: mu_wge(2)=0.0_dp
	real(dp) :: sig_wge(2),mu_mar(ntypp),sig_mar,ro,mu_o , sigo_m,sigo_f
	real(dp) :: uhome(2),uhomet(ntypp),alphaed(2,neduc),alphakid(nkid) !ahu october2022: changing alphakid so that it doesn't have that obsolete dimension anymore
	real(dp) :: cst(ntypp),kcst,divpenalty,uloc(nl),sig_uloc
	real(dp) :: alf10(nl),alf11,alf12,alf13,alf1t(ntypp)            ! types
	real(dp) :: alf20(nl),alf21,alf22,alf23,alf2t(ntypp)            ! types
	real(dp) :: ptype,pmeet,omega(2),ptypehs(ntypp),ptypecol(ntypp) ! types
	real(dp) :: pkid,psio(12),psil(2),psih !getting rid of psih2,3,4 so that now psih is just a scalar
	real(dp) :: popsize(nl)
	integer(i4b) :: distance(nl,nl)
	real(dp) :: wg(np,2),wgt(np),mg(nz,ninp),mgt(nz),best
	type :: initcond
        integer(i4b) :: id
		integer(i4b) :: co				
		integer(i4b) :: sexr			
		integer(i4b) :: hme
		integer(i4b) :: endage		
		integer(i4b) :: edr			
    end type
	type, extends(initcond) :: statevar
        integer(i4b) :: expr
        integer(i4b) :: kidr      
        integer(i4b) :: hhr		! annual hours worked by r (turned into discrete in read_data)
		real(dp) :: logwr,wr	! wr is annual income and logwr is log of annual income. hourly wage (wr_perhour or wsp_perhour) is read from the data (see familymig_2.do to see how it's turned into hourly) and that hourly wage is turned into annual by multiplying it by h_wmult	  
        integer(i4b) :: l		! location		
        integer(i4b) :: rel		! relationship status. -1: not observed, 0: single, 1: married, 2: cohabiting
		integer(i4b) :: rellen	! length of current relationship starting from 1 in first period
        integer(i4b) :: edsp
        integer(i4b) :: expsp
        integer(i4b) :: kidsp
        integer(i4b) :: hhsp	! annual hours worked by spouse (turned into discrete in read_data)
		real(dp) :: logwsp,wsp  ! see explanation for logwr,wr.  
        integer(i4b) :: lsp
        integer(i4b) :: nomiss
        integer(i4b) :: nn,mm,r,typ
	end type
	type :: shock
		real(dp) :: meet 
		real(dp) :: marie
		real(dp) :: meetq 
		real(dp) :: meetx 
		real(dp) :: q
		real(dp) :: x
        real(dp) :: iepsmove
        real(dp) :: typ
	end type	
	type(statevar) :: ones
    type(initcond) :: ones_init
    type :: taxo
        real(dp) :: pwages
        real(dp) :: swages 
        real(dp) :: statesin
        real(dp) :: fedsin
        real(dp) :: statemar
        real(dp) :: fedmar
    end type 
    integer(i4b), parameter :: numbin=31
    integer(i4b), parameter :: numtaxes=nl*numbin*numbin 
    integer(i4b) :: taxset !set in main
    type(taxo) :: tax(numbin,numbin,nl) !pwages,swages.myreg
    real(dp) :: pbracket(numbin),sbracket(numbin)
	type(initcond), dimension(:), allocatable :: init !array for initial conditions (size just num of persons in actual data)
contains

	! get parameters from transformed values. input is free
	! to be any real(sp) resulting parameters are appropriately constrained
	subroutine getpars(par,realpar)
	real(dp), dimension(npars), intent(in)  :: par ! vector of parameters
	real(dp), dimension(npars), intent(out) :: realpar ! vector of parameters
	integer(i4b) :: g,i,j,ed,indust1(ntypp),indust2(ntypp)
    integer(i4b) :: dw0,de,dsex
    real(dp) :: temprob(nl),junk

    stepos=0.0_dp
    indust1=0
    indust2=0
    temprob=0.0_dp 

	realpar=pen 
    parname=''
    j=1
    !ahu jan19 012819: not iterating on ed offers anymore. replacing them with curloc and ofloc offers instead 
    !note that realpar's for psio parameters are reassigned at the end of this file just for visual purpoes, to write those in writemoments.
    !but the actual values that are used are assigned to psio right here and those are the ones that are used in fnprof.
	realpar(j)=par(j)               ; parname(j)='emp,cur,m' ; stepos(j)=-1.0_dp  ; if (onlyfem) stepos(j)=0.0_dp !psio is for the offer function
	psio(1)=realpar(j)	            ; j=j+1
	realpar(j)=par(j)               ; parname(j)='emp,cur,m' ; stepos(j)=-1.0_dp  ; if (onlyfem) stepos(j)=0.0_dp !psio is for the offer function
	psio(2)=realpar(j)	            ; j=j+1
	realpar(j)=par(j)               ; parname(j)='emp,of,m' ; stepos(j)=-1.0_dp  ; if (onlyfem) stepos(j)=0.0_dp!psio is for the offer function
	psio(3)=realpar(j)	            ; j=j+1
	realpar(j)=0.0_dp               ; parname(j)='emp,of,m' ; stepos(j)=0.0_dp  ; if (onlyfem) stepos(j)=0.0_dp !psio is for the offer function !NO MORE OFLOC LAYOFF NONSENSE
	psio(4)=realpar(j)	            ; j=j+1
	realpar(j)=par(j)               ; parname(j)='emp,cur,f' ; stepos(j)=1.0_dp ; if (onlymales) stepos(j)=0.0_dp !psio is for the offer function
	psio(5)=realpar(j)	            ; j=j+1
	realpar(j)=par(j)               ; parname(j)='emp,cur,f' ; stepos(j)=1.0_dp ; if (onlymales) stepos(j)=0.0_dp  !psio is for the offer function
	psio(6)=realpar(j)	            ; j=j+1
	realpar(j)=par(j)               ; parname(j)='emp,of,f' ; stepos(j)=1.0_dp ; if (onlymales) stepos(j)=0.0_dp  !psio is for the offer function
	psio(7)=realpar(j)	            ; j=j+1
	realpar(j)=0.0_dp               ; parname(j)='emp,of,f' ; stepos(j)=0.0_dp ; if (onlymales) stepos(j)=0.0_dp  !psio is for the offer function !NO MORE OFLOC LAYOFF NONSENSE
	psio(8)=realpar(j)	            ; j=j+1
	realpar(j)=par(j)               ; parname(j)='u,cur,m' ; stepos(j)=-1.0_dp  ; if (onlyfem) stepos(j)=0.0_dp !psio is for the offer function
	psio(9)=realpar(j)	            ; j=j+1
	realpar(j)=par(j)               ; parname(j)='u,of,m' ; stepos(j)=-1.0_dp  ; if (onlyfem) stepos(j)=0.0_dp !psio is for the offer function
	psio(10)=realpar(j)	            ; j=j+1
	realpar(j)=par(j)               ; parname(j)='u,cur,f' ; stepos(j)=1.0_dp  ; if (onlymales) stepos(j)=0.0_dp !psio is for the offer function
	psio(11)=realpar(j)	            ; j=j+1
	realpar(j)=par(j)               ; parname(j)='u,of,f' ; stepos(j)=1.0_dp  ; if (onlymales) stepos(j)=0.0_dp !psio is for the offer function
	psio(12)=realpar(j)	            ; j=j+1
    !print*, 'Here is psio12',j-1
    !note that realpar's for psio parameters are reassigned at the end of this file just for visual purpoes, to write those in writemoments.
    !but the actual values that are used are assigned to psio right here and those are the ones that are used in fnprof.

	realpar(j)=par(j)               ; parname(j)='psil(1)' ; stepos(j)=0.0_dp 
	psil(1)=realpar(j)	            ; j=j+1
	realpar(j)= par(j)             ; parname(j)='uhomet 1'	; stepos(j)=1.0_dp*par(j) 
	uhomet(1)=realpar(j)            ; j=j+1
	realpar(j)= par(j)             ; parname(j)='uhomet 2'	; stepos(j)=1.0_dp*par(j) !2.0_dp*(1.0_dp/(1.0_dp+exp(-par(j))))-1.0_dp 
	uhomet(2)=realpar(j)            ; j=j+1
    !if (iwritegen==1) print*, "Here is ro", ro, par(j-1)
    !note that realpar's for psih parameters are reassigned at the end of this file just for visual purpoes, to write those in writemoments.
    !but the actual values that are used are assigned to psih right here and those are the ones that are used in fnprhc.
	realpar(j)=par(j)               ; parname(j)='p(ex=2|ex=1),e' ; stepos(j)=1.0_dp !this is psih, the only one that governs fnprhc.  
	psih=realpar(j)	            ; j=j+1
	realpar(j)=0.0_dp               ; parname(j)='p(ex=1|ex=1),e' ; stepos(j)=0.0_dp !this is just for visuals. not a parameter anymore. 
	junk=realpar(j)	            ; j=j+1
	realpar(j)=0.0_dp               ; parname(j)='p(ex=2|ex=2),e' ; stepos(j)=0.0_dp !this is just for visuals. not a parameter anymore. 
	junk=realpar(j)	            ; j=j+1
	realpar(j)=0.0_dp               ; parname(j)='p(ex=1|ex=2),e' ; stepos(j)=0.0_dp !this is just for visuals. not a parameter anymore. 
	junk=realpar(j)	            ; j=j+1
    !note that realpar's for psih parameters are reassigned at the end of this file just for visual purpoes, to write those in writemoments.
    !but the actual values that are used are assigned to psih right here and those are the ones that are used in fnprhc.

    
    realpar(j)=logit(par(j))        ; parname(j)='pkid' ; stepos(j)=0.0_dp  ; if (onlysingles) stepos(j)=0.0_dp !20
	pkid=realpar(j)                 ; j=j+1
    realpar(j)=logit(par(j))	    ; parname(j)='pmeet' ; stepos(j)=0.0_dp ; if (onlysingles) stepos(j)=0.0_dp !21
	pmeet=realpar(j)                ; j=j+1


    !realpar(j) = mult1c * logit(par(j))              ; parname(j)='uhome(1)' ; stepos(j)=0.2_dp	!mult3*logit(par(2:3)) !22:23
	!uhome(1)=realpar(j)                             ; j=j+1
    !realpar(j) = mult1c * logit(par(j))              ; parname(j)='uhome(2)' ; stepos(j)=0.2_dp	 ; if (onlymales) stepos(j)=0.0_dp !mult3*logit(par(2:3)) !22:23
	!uhome(2)=realpar(j)                             ; j=j+1

    realpar(j) = par(j)             ; parname(j)='uhome(1)' ; stepos(j)=0.5_dp*PAR(J)	 ; if (onlyfem) stepos(j)=0.0_dp !mult3*logit(par(2:3)) !22:23
	uhome(1)=realpar(j)                             ; j=j+1
    realpar(j) = par(j)              ; parname(j)='uhome(2)' ; stepos(j)=0.5_dp*PAR(J)	 ; if (onlymales) stepos(j)=0.0_dp !mult3*logit(par(2:3)) !22:23
	uhome(2)=realpar(j)                             ; j=j+1
    realpar(j)=par(j)            ; parname(j)='uhomet 3'	; stepos(j)=0.5_dp*par(j) !24 !-1.0_dp*mult1c * logit(par(j)) !ahu 112718 changing to only minus from: mult1 * min2pls(par(j))     ! types
    uhomet(3)=realpar(j)                                     ; j=j+1               ! types
	realpar(j) = par(j)          ; parname(j)='kcst'	; stepos(j)=2.0_dp*(-1000.0_dp) !*par(j) !25 !-1.0_dp*mult1c * logit(par(j)) !ahu 112718 changing to only minus from: mult1 * min2pls(par(5)) !mult2*logit(par(4:6))	
	kcst=realpar(j)                                     ; j=j+1
	realpar(j) = -1.0_dp*multdiv * logit(par(j))          ; parname(j)='divpenalty'	; stepos(j)=2.0_dp ; if (onlysingles) stepos(j)=0.0_dp !26 !ahu 112718 changing to only minus from: mult1 * min2pls(par(6))                         !ahu summer18 050418: changed from 1000 to 10,000 (mult to mult1)
	divpenalty=realpar(j)                               ; j=j+1
    !print*, 'Here is divpenalty',j-1,divpenalty 

    realpar(j:j+1) = mult1 * logit(par(j:j+1))          ; parname(j)='alphaed(m,ned)' ; parname(j+1)='alphaed(f,ned)'    !27:28   !ahu jan19 012719 changing it yet again back to logit because there is not that much of different in objval between alpha=0 and alpha=-49000    !ahu jan19 012019 changing it back to min2pls  ! noed !ahu 112718 changing to only plus from: mult1*min2pls(par(7:8))   !mult1 * logit(par(7))	
	stepos(j)=2.0_dp            ; if (onlyfem) stepos(j)=1.0_dp 
    stepos(j+1)=2.0_dp          ; if (onlymales) stepos(j+1)=1.0_dp 
    alphaed(:,1)=realpar(j:j+1)                         ; j=j+2 !alphaed(m:f,noed)  [educ=1 noed, educ=2 ed]  mult1 * min2pls(par(j:j+1))
    
    realpar(j:j+1) = mult1 * logit(par(j:j+1))          ; parname(j)='alphaed(m,ed)' ; parname(j+1)='alphaed(f,ed)'    !27:28   !ahu jan19 012719 changing it yet again back to logit because there is not that much of different in objval between alpha=0 and alpha=-49000    !ahu jan19 012019 changing it back to min2pls  ! noed !ahu 112718 changing to only plus from: mult1*min2pls(par(7:8))   !mult1 * logit(par(7))	
	stepos(j)=2.0_dp            ; if (onlyfem) stepos(j)=1.0_dp 
    stepos(j+1)=2.0_dp          ; if (onlymales) stepos(j+1)=1.0_dp 
    alphaed(:,2)=realpar(j:j+1)                         ; j=j+2 !alphaed(m:f,ed)  [educ=1 noed, educ=2 ed]  mult1 * min2pls(par(j:j+1))
    
    realpar(j:j+1)=mult1 * logit(par(j:j+1))            ; parname(j)='alphakid(m)' ; parname(j+1)='alphakid(f)'          !31:32           !ahu 112718 changing to only plus from: mult1 * min2pls(par(j:j+1))	 !mult1 * logit(par(9:10))	
    stepos(j)=2.0_dp            ; if (onlyfem) stepos(j)=0.0_dp ; 	stepos(j+1)=1.0_dp   ; if (onlymales) stepos(j:j+1)=0.0_dp 
    alphakid(:)=realpar(j:j+1)                        ; j=j+2         
    !print*, 'Here is uloc',j
	
    !uloc: 33-41
    !do ed=1,2
        do i=1,nl
            if (i==2) then
			    realpar(j) = 0.0_dp  ; stepos(j)=0.0_dp
			    uloc(i)=0.0_dp
		    else 
			    realpar(j) = par(j) ; stepos(j)=4.0_dp*PAR(J)    !mult1 * min2pls( par(j) )
			    uloc(i)=realpar(j)
		    end if 
            parname(j)='uloc' 
            j=j+1
        end do
	!end do
    !print*, 'Here is alf10',j
	!wage 42: 65
    do i=1,nl
        if (i==5) then
            realpar(j)=0.0_dp !1.5_dp*min2pls(par(j))+8.5_dp 
            alf10(i)=realpar(j)
        else 
            realpar(j)=par(j) !1.5_dp*min2pls(par(j))+8.5_dp 
            alf10(i)=realpar(j)
        end if     
        parname(j)='alf10' ; stepos(j)=0.3_dp ; if (i==3) stepos(j)=0.0_dp  ; if (onlyfem) stepos(j)=0.0_dp !dont iterate on alpha10(1) for a second
        j=j+1
    end do 
    !print*, 'Here is alf11 etc',j
    realpar(j)=logit(par(j))                        ; parname(j)='alf11' ; stepos(j)=0.5_dp  ; if (onlyfem) stepos(j)=0.0_dp
	alf11=realpar(j)                                ; j=j+1
    !print*, 'Here is alf12',j	
    realpar(j)=logit(par(j))                 ; parname(j)='alf12' ; stepos(j)=0.5_dp  ; if (onlyfem) stepos(j)=0.0_dp
    alf12=realpar(j)                                ; j=j+1
    !print*, 'Here is alf13',j	
    realpar(j)=0.0_dp                               ; parname(j)='alf13' ; stepos(j)=0.0_dp   ; if (onlyfem) stepos(j)=0.0_dp  !-1.0_dp*logit(par(j)) 
    alf13=realpar(j)	                            ; j=j+1
    !print*, 'Here is alf20',j
    do i=1,nl
        if (i==5) then
            realpar(j)=0.0_dp !1.5_dp*min2pls(par(j))+8.5_dp 
            alf20(i)=realpar(j)
        else 
            realpar(j)=par(j) !1.5_dp*min2pls(par(j))+8.5_dp 
            alf20(i)=realpar(j)
        end if     
        parname(j)='alf20' ; stepos(j)=0.3_dp  ; if (i==3) stepos(j)=0.0_dp ;  if (onlymales) stepos(j)=0.0_dp 
		j=j+1
	end do 
    !print*, 'Here is alf21 etc',j
	realpar(j)=logit(par(j))                        ; parname(j)='alf21' ; stepos(j)=0.3_dp  ; if (onlymales) stepos(j)=0.0_dp 
	alf21=realpar(j)                                ; j=j+1
    !print*, 'Here is alf22',j	    
    realpar(j)=logit(par(j))                 ; parname(j)='alf22' ; stepos(j)=0.3_dp  ; if (onlymales) stepos(j)=0.0_dp 
	alf22=realpar(j)                                ; j=j+1
    !print*, 'Here is alf23',j	
    realpar(j)=par(j)                              ; parname(j)='uhomet 4' ; stepos(j)=par(j)  ; if (onlymales) stepos(j)=0.0_dp !-1.0_dp*logit(par(j)) 
	uhomet(4)=realpar(j)	                            ; j=j+1
	
    realpar(j:j+1)=logit(par(j:j+1))                ; parname(j:j+1)='sig_wge'	; stepos(j)=1.0_dp ; stepos(j+1)=1.0_dp	  ; if (onlyfem) stepos(j)=0.0_dp  ; if (onlymales) stepos(j+1)=0.0_dp !66:67
	sig_wge(1:2)=realpar(j:j+1)                     ; j=j+2
    !sigom and sigof: 68:69
    realpar(j)=multsigo * logit(par(j))                               ; parname(j)='sigo_m'	; stepos(j)=0.5_dp ; if (nepsmove==1) stepos(j)=0.0_dp ; if (onlyfem) stepos(j)=0.0_dp
    !print*, "Here it is sigom", j,par(j),realpar(j)
    sigo_m=realpar(j)                                ; j=j+1
    realpar(j)=multsigo * logit(par(j))                               ; parname(j)='sigo_f'	; stepos(j)=0.5_dp ; if (nepsmove==1) stepos(j)=0.0_dp ; if (onlymales) stepos(j)=0.0_dp
    !print*, "Here it is sigof", j,par(j),realpar(j)
    sigo_f=realpar(j)                                ; j=j+1

    
    do i=1,ntypp !The below are parameters 70 to 93. So 70-75 is type1, 76-81 is type2, 82-87 is type3, 88-93 is type4.
        if (i==1) then 
            realpar(j)=0.0_dp                           ; parname(j)='ptypehs' ; stepos(j)=1.0_dp
            ptypehs(i)=exp(realpar(j))                  ; indust1(i)=j ; j=j+1
            realpar(j)=0.0_dp                           ; parname(j)='ptypecol'  ; stepos(j)=1.0_dp
            ptypecol(i)=exp(realpar(j))                 ; indust2(i)=j ; j=j+1
            realpar(j)=par(j)                           ; parname(j)='alf1t'     ; stepos(j)=0.3_dp  ; if (onlyfem) stepos(j)=0.0_dp
            alf1t(i)=realpar(j)                         ; j=j+1
            realpar(j)=par(j)                          ; parname(j)='alf2t'     ; stepos(j)=0.3_dp  ; if (onlymales) stepos(j)=0.0_dp
	        alf2t(i)=realpar(j)                         ; j=j+1
            !realpar(j)= -1.0_dp*mult1c * logit(par(j))   ; parname(j)='cst'       ; stepos(j)=0.5_dp
            !cst(i)=realpar(j)                           ; j=j+1 
            realpar(j)=par(j)                          ; parname(j)='cst'       ; stepos(j)=2.0_dp*(-5000.0_dp) !not iterating on this anymore. see notes. under cost vs. sigo. they are just not sep ident I think. 
            cst(i)=realpar(j)                           ; j=j+1 
            !ahu082822 august2022 print*, 'mumar(1)',j,par(j),multmar, min2pls(par(j)),multmar*min2pls(par(j))
            realpar(j)=multmar * logit(par(j))          ; parname(j)='mu_mar'     ; stepos(j)=1.0_dp    ; if (onlysingles) stepos(j)=0.0_dp 	    
            mu_mar(i)=realpar(j)                        ; j=j+1      
        else
            realpar(j)=par(j)                            ; parname(j)='ptypehs' ; stepos(j)=1.0_dp
            ptypehs(i)=exp(realpar(j))                  ; indust1(i)=j ; j=j+1
            realpar(j)=par(j)                            ; parname(j)='ptypecol'  ; stepos(j)=1.0_dp
            ptypecol(i)=exp(realpar(j))                 ; indust2(i)=j ; j=j+1
            realpar(j)=par(j)                           ; parname(j)='alf1t'     ; stepos(j)=0.3_dp  ; if (onlyfem) stepos(j)=0.0_dp
            alf1t(i)=realpar(j)                         ; j=j+1
            realpar(j)=par(j)                           ; parname(j)='alf2t'     ; stepos(j)=0.3_dp  ; if (onlymales) stepos(j)=0.0_dp
	        alf2t(i)=realpar(j)                         ; j=j+1
            !if (i<4) then
            !    realpar(j)= par(j)                         ; parname(j)='cst'       ; stepos(j)=1.0_dp*par(j)
            !    cst(i)=realpar(j)                           ; j=j+1 
            !else if (i==4) then
            realpar(j)= par(j)                         ; parname(j)='cst'       ; stepos(j)=2.0_dp*(-5000.0_dp)
            cst(i)=realpar(j)                           ; j=j+1 
            !end if
            realpar(j)=multmar * logit(par(j))                         ; parname(j)='mu_mar'     ; stepos(j)=1.0_dp   ; if (onlysingles) stepos(j)=0.0_dp 	    
            mu_mar(i)=realpar(j)                        ; j=j+1      
        end if 
    	!print*, 'Here is cost',cst(i),par(j-1),realpar(j-1)
        !ahu 122818 changed mult1 to multmar 
    end do 
    ptypehs(:)=ptypehs(:)/sum(ptypehs)
    ptypecol(:)=ptypecol(:)/sum(ptypecol)
    do i=1,ntypp
        realpar(indust1(i))=ptypehs(i)
        realpar(indust2(i))=ptypecol(i)
    end do 

    !******************** ahu october2022 **********************************************
    ! writing of realpar for psio parameters in getpars by calling fnprof 
    ! writing of realpar for psih parameters in getpars by calling fnprof 
    ! note that the below is just for visual purposes. realpar is not used in the solution or simulation. the actual parameters psio are used. 
    ! and here I am not reassigning values of psio. just reassignign values of realpar so I can write it in writemoments. 
    ! see writemoments file for more detailed writing of the fnprof parameters. 
    temprob(1:3)=fnprof(np,5,1) !emp curloc m
    realpar(1:2)=temprob(1:2) 
    temprob(1:3)=fnprof(np,10,1) !emp ofloc m
    realpar(3:4)=temprob(1:2) 
    temprob(1:3)=fnprof(np,5,2) !emp curloc f 
    realpar(5:6)=temprob(1:2) 
    temprob(1:3)=fnprof(np,10,2) !emp ofloc f
    realpar(7:8)=temprob(1:2) 
    temprob(1:3)=fnprof(np1,5,1) !unemp curloc m
    realpar(9)=temprob(1) 
    temprob(1:3)=fnprof(np1,10,1) !unemp ofloc m
    realpar(10)=temprob(1) 
    temprob(1:3)=fnprof(np1,5,2) !unemp curloc f 
    realpar(11)=temprob(1) 
    temprob(1:3)=fnprof(np1,10,2) !unemp ofloc f
    realpar(12)=temprob(1) 
    if (nexp>2) then ; print*, "it is only ok to write this way if nexp is 2! so beware" ; stop ; end if
    temprob(1:nexp)=fnprhc(1,np) !when experience is 1 and when working
    realpar(16)=temprob(2) !this is prob of moving to experience=2 when your experience is 1. !
    realpar(17)=temprob(1) !this is prob of moving to experience=1 when your experience is 1. 
    temprob(1:nexp)=fnprhc(nexp,np) !when experience is 2 and when working
    realpar(18)=temprob(2) !this is prob of moving to experience=2 when your experience is 2.  
    realpar(19)=temprob(1) !this is prob of moving to experience=1 when your experience is 2.
    !note that we are not writing the fnprhc(.,np1) because that is just prob of staying where you are is 1. Check this in the writing of fnprhc in writemoments. 
    !******************** ahu october2022 **********************************************
    temprob(1:nl)=fnprloc(1) !if origin location is loc1, what is the probability of drawing location 1
    realpar(13)=temprob(1)
    !******************** ahu october2022 **********************************************
    !Assign Cendiv names 
    locname(1)='New England        '
    locname(2)='Mid Atlantic       '
    locname(3)='East North Central '
    locname(4)='West North Central '
    locname(5)='South Atlantic     '
    locname(6)='East South Central '
    locname(7)='West South Central '
    locname(8)='Mountain           '
    locname(9)='Pacific            '
    !******************** ahu october2022 **********************************************

    !alphaed(:,2)=alphaed(:,1)
    mu_o=0.0_dp
    sig_mar=0.0_dp
    ro=0.0_dp !no longer a parameter 

    !***********************
    !ahu 041118 del and remove later:
    !alphaed(2,:)=alphaed(1,:)
    !psio(5:8)=psio(1:4)
    !psio(11:12)=psio(9:10)
    !alphakid(2,:)=alphakid(1,:)
    
    !alf20=alf10
    !alf21=alf11
    !alf22=alf12
    !alf23=alf13
    !uhome(2)=uhome(1)
    !***********************
    
    !psio(1:4)=psio(5:8)
    !psio(9:10)=psio(11:12)
    
    !pkid=0.0_dp
    !alphakid=0.0_dp
    !kcst=0.0_dp 
    !ro=0.0_dp
    !sig_mar=0.0_dp
    !scst=0.0_dp
    
    !alf11=0.0_dp
    !alf21=0.0_dp
    !psio(3:4)=psio(1:2)
    !psio(7:8)=psio(5:6)
    !psio(10)=psio(9)
    !psio(12)=psio(11)
    !uloc(:,2)=uloc(:,1)
    !alphaed(:,2)=alphaed(:,1)
    !alf1t(2)=0.0_dp
    !alf2t(2)=0.0_dp
    !ptypecol=ptypehs
    
    
    !cst=0.0_dp
    !kcst=0.0_dp
    
    !ro=0.0_dp !0.98_dp
    !alpha=0.0_dp
    !kcst=0.0_dp
    !pkid=0.0_dp
    !alf1t(2)=alf1t(1)
    !alf2t(2)=alf2t(1)
    !cst(1)=cst(2)
    !uhome=0.0_dp

    !uhome=0.0_dp
    !cst=-150000.0_dp
    !kcst=-150000.0_dp
    !divpenalty=0.0_dp
    !pkid=0.0_dp
    !kcst=0.0_dp
    !alpha(:,2)=alpha(:,1)
    
    !sig_o=sig_mar    
    !alpha(1,:)=alpha(2,:) 
    	!if ((.not.optimize).and.(.not.chkstep) ) print*, "sig_mar,mu_mar,npars;,; ", sig_mar,mu_mar,j
    
	!print*, logitinv(alf11),logitinv(alf12),logitinv(alf13)
	!print*, logitinv(alf21),logitinv(alf22),logitinv(alf23)

	!if (j/=npars) then ; print*, "something wrong in getpar! ",j,npars ; stop ; end if
    

    
        realpartemp=realpar
	end subroutine getpars

	subroutine getsteps(par,step)
	real(dp), dimension(:), intent(in) :: par 
	!character(len=15), dimension(:), intent(out) :: name ! names of each parameter
	real(dp), dimension(:), intent(out) :: step 
	integer(i4b) :: i,j
    step=0.0_dp
	end subroutine getsteps

	subroutine getdistpop
	integer(i4b) :: i,j
	distance=0
	popsize=0.0_dp 
    !ahu 030717: redid this adjacence thing. see map and appendix of the draft. 
	distance(1,2)=1 
	distance(2,1)=1 
	distance(2,5)=1 
	distance(2,3)=1 
	distance(3,2)=1 
	distance(3,4)=1 
	distance(3,6)=1 
	distance(3,5)=1 
	distance(4,3)=1 
	distance(4,8)=1 
	distance(4,7)=1 
	distance(4,6)=1 
	distance(5,8)=1 
	distance(5,3)=1 
	distance(5,6)=1 
	distance(6,5)=1 
	distance(6,7)=1 
	distance(6,4)=1 
	distance(6,3)=1 
	distance(7,6)=1 
	distance(7,8)=1 
	distance(7,4)=1 
	distance(8,7)=1 
	distance(8,9)=1 
	distance(8,4)=1 
	distance(9,8)=1 
    do i=1,nl 
		do j=1,nl 
			if (j==i) then 
				if (distance(j,i)==1) then 		! ahu 061513: for some locaitons, you did have that their distance(i,i) was 1 so corrected this
                    print*, 'distance(i,i) should not be 1 because that is not adjacence'
                    stop
                end if 
			end if 
		end do 
	end do 

	popsize(1)=0.9112_dp	!new england
	popsize(2)=2.695_dp		!middle atlantic 
	popsize(3)=2.9598_dp	!east north central
	popsize(4)=1.2468_dp	!west north central
	popsize(5)=2.5736_dp	!south atlantic
	popsize(6)=1.0089_dp	!east south central
	popsize(7)=1.6121_dp	!west south central
	popsize(8)=0.7661_dp	!mountain
	popsize(9)=2.2757_dp	!pacific
	end subroutine getdistpop

	subroutine getones
		ones%co=-99
		ones%sexr=-99
		ones%hme=-99
		ones%endage=-99
		ones%edr=-99
		ones%expr=-99
        ones%kidr=-99
		ones%hhr=-99
        ones%logwr=-99.0_dp
        ones%wr=-99.0_dp
        ones%wsp=-99.0_dp
        ones%l=-99
        ones%rel=-99
        ones%rellen=-99
        ones%edsp=-99
        ones%expsp=-99
        ones%kidsp=-99
        ones%hhsp=-99
        ones%logwsp=-99.0_dp
		ones%lsp=-99        
		ones%nomiss=-99
        ones%nn=-99
        ones%mm=-99
        ones%r=-99
        ones%typ=-99
        
        ones_init%id=-99
        ones_init%co=-99
        ones_init%sexr=-99
        ones_init%hme=-99
        ones_init%endage=-99
        ones_init%edr=-99
        end subroutine getones	

        function fnwge(dg,dtyp,dl,dw,de,dr)		!de is educ here but in fnprof it's no longer educ				
            integer(i4b), intent(in) :: dg,dtyp,dl,de,dr						! gender,typ,location,education,experience
            real(dp), intent(in) :: dw								! wage draw
            real(dp) :: fnwge
            if (dg==1) then 
                fnwge=exp(alf1t(dtyp)+alf10(dl)+alf11*one(de==2) + alf12*(dr-1) + alf13*((dr-1)**2) + dw ) 
            else if (dg==2) then 
                fnwge=exp(alf2t(dtyp)+alf20(dl)+alf21*one(de==2) + alf22*(dr-1) + alf23*((dr-1)**2) + dw ) 
            end if 
            end function fnwge
        
            function fnprof(dw0,de,dsex) !ahu october2022: note that de was ed before but now it's wtr it's curloc or ofloc (takes on values 5 or 10)
            integer(i4b), intent(in) :: dw0,de,dsex
            real(dp), dimension(3) :: fnprof
            if (onthejobsearch) then 
                fnprof=0.0_dp
                if ( dw0 <= np ) then  
                    if ( de==5 .and. dsex==1 ) then 
                        fnprof(1:2)=exp(psio(1:2)) !exp( psio(1) + psio(2) * abs(dsex==1) + psio(3) * abs(de==2) )	! offer 	 
                    else if ( de==10 .and. dsex==1 ) then 
                        fnprof(1)=exp(psio(3)) 
                        fnprof(2)=0.0_dp !psio4 nuissance parameter
                    else if ( de==5 .and. dsex==2 ) then 
                        fnprof(1:2)=exp(psio(5:6)) 
                    else if ( de==10 .and. dsex==2 ) then 
                        fnprof(1)=exp(psio(7)) 
                        fnprof(2)=0.0_dp !psio8 nuissance parameter
                    end if 
                    fnprof(3)=exp( 0.0_dp )		! nothing happens												
                else if (dw0 == np1) then 
                    if ( de==5 .and. dsex==1 ) then 
                        fnprof(1)=exp(psio(9)) !exp( psio(1) + psio(2) * abs(dsex==1) + psio(3) * abs(de==2) )	! offer 	 
                    else if ( de==10 .and. dsex==1 ) then 
                        fnprof(1)=exp(psio(10) )
                    else if ( de==5 .and. dsex==2 ) then 
                        fnprof(1)=exp(psio(11)) 
                    else if ( de==10 .and. dsex==2 ) then 
                        fnprof(1)=exp(psio(12)) 
                    end if 
                    fnprof(2)=0.0_dp		! 0 since you can't get laid off if you don't have a job! 
                    fnprof(3)=exp(0.0_dp)		! nothing happens												
                else  
                    print*, "in fnprof: dw0 > np1 which doesnt' make sense as that's a state variable " , dw0,de,dsex
                    stop
                end if 
                fnprof(1:3)=fnprof/sum(fnprof)
            else 
                fnprof=0.0_dp
                if ( dw0 <= np ) then  
                    if ( de==5 .and. dsex==1 ) then 
                        fnprof(1:2)=0.0_dp !exp(psio(1:2)) !exp( psio(1) + psio(2) * abs(dsex==1) + psio(3) * abs(de==2) )	! offer 	 
                    else if ( de==10 .and. dsex==1 ) then 
                        fnprof(1)=exp(psio(3)) 
                        fnprof(2)=0.0_dp !psio4 nuissance parameter
                    else if ( de==5 .and. dsex==2 ) then 
                        fnprof(1:2)=0.0_dp !exp(psio(5:6)) 
                    else if ( de==10 .and. dsex==2 ) then 
                        fnprof(1)=exp(psio(7)) 
                        fnprof(2)=0.0_dp !psio8 nuissance parameter
                    end if 
                    fnprof(3)=exp( 0.0_dp )		! nothing happens												
                else if (dw0 == np1) then 
                    if ( de==5 .and. dsex==1 ) then 
                        fnprof(1)=exp(psio(9)) !exp( psio(1) + psio(2) * abs(dsex==1) + psio(3) * abs(de==2) )	! offer 	 
                    else if ( de==10 .and. dsex==1 ) then 
                        fnprof(1)=exp(psio(10) )
                    else if ( de==5 .and. dsex==2 ) then 
                        fnprof(1)=exp(psio(11)) 
                    else if ( de==10 .and. dsex==2 ) then 
                        fnprof(1)=exp(psio(12)) 
                    end if 
                    fnprof(2)=0.0_dp		! 0 since you can't get laid off if you don't have a job! 
                    fnprof(3)=exp(0.0_dp)		! nothing happens												
                else  
                    print*, "in fnprof: dw0 > np1 which doesnt' make sense as that's a state variable " , dw0,de,dsex
                    stop
                end if 
                fnprof(1:3)=fnprof/sum(fnprof)
            end if 
            if (skriv) then 
                if ( abs(  sum(fnprof) - 1.0_dp  ) > eps) then ; print*, "error in getfnprof : offer does not add up " , sum(fnprof) ; stop ; end if 
            end if 
            end function fnprof
        
        
            
            
            function fnprloc(orig)
            integer(i4b), intent(in) :: orig		! origin location
            real(dp), dimension(nl) :: fnprloc
            integer(i4b) :: j
            real(dp), dimension(nl) :: temp
            fnprloc=0.0_dp
            temp=0.0_dp
            !ahu 030717 fnprloc(orig)=logit(psil(1)) !ahu 030717: putting psil(1) inside the exp instead because otherwise very high prloc from origin and low from others and 
                                                     !            we get very low moving rates especially for married people. 
            do j=1,nl	
                !ahu 030717 if (j /= orig) then 
                !ahu 030717 	fnprloc(j)= exp( psil(2) * distance(j,orig) + psil(3) * popsize(j) ) 
                !ahu 030717 	sum_sans_orig = sum_sans_orig + fnprloc(j) 
                !ahu 030717 end if 
                temp(j)= exp( psil(1) * one(j==orig) )   ! + psil(2) * distance(j,orig) )    !+ psil(3) * popsize(j) ) 
            end do 
            do j=1,nl	
                !ahu 030717 if (j /= orig) then 
                !ahu 030717 	fnprloc(j)=(1.0_dp-fnprloc(orig) ) * fnprloc(j)/sum_sans_orig
                !ahu 030717 end if 
                fnprloc(j)=temp(j)/sum(temp(:))
            end do 
            if ( abs(  sum(fnprloc) - 1.0_dp  ) > eps) then ; print*, "error in getfnprloc : offer does not add up " , sum(fnprloc) ; stop ; end if 
            end function fnprloc
            
            function fnprhc(dr,dw)
                integer(i4b), intent(in) :: dr,dw		! experience and employment: w<=np work, w==np1 not work,  w=np2 nothing/can't be a state variable here so if you get this, there's something wrong
                real(dp), dimension(nexp) :: fnprhc
                integer(i4b) :: j
                if (skriv) then 
                    if ( dw > np1 ) then ; print*, "in fnprof: dw0 > np1 which doesnt' make sense as that's a state variable " ; stop ; end if 
                end if 
                !ahu october2022: note that when j=nexp, there will be no j such that j-dr=+1, so fnprhc(nexp) will be 1/0+1+1 = 1 and all other fnprhc(j)'s are 0. 
                fnprhc=0.0_dp
                do j=1,nexp	
                    if ( dw <= np ) then 
                        if ( j-dr == +1 ) then  
                            fnprhc(j)=exp(psih)      !ahu jan19 011719 changing to logit
                        else if (j==dr) then 
                            fnprhc(j)=exp(0.0_dp)  
                        else if ( j-dr == -1 ) then  !ahu jan19 011519 getting rid of probdown
                            fnprhc(j)=0.0_dp
                        else 
                            fnprhc(j)=0.0_dp
                        end if 
                    else if ( dw == np1 ) then !ahu october2022 no exp increase or decrease if unemp. so j such that j=dr is 1/0+1+0 =1 and all other fnprhc(j)'s are 0. 
                        if ( j-dr == +1 ) then  
                            fnprhc(j)= 0.0_dp    !ahu jan19 011719 changing to logit
                        else if (j==dr) then 
                            fnprhc(j)=exp(0.0_dp)      !exp(0.0_dp)   !ahu jan19 011719 changing to logit
                        else if ( j-dr == -1 ) then  !ahu jan19 011519 getting rid of probdown
                            fnprhc(j)= 0.0_dp
                        else 
                            fnprhc(j) = 0.0_dp
                        end if 
                    end if 
                end do 	
                fnprhc(:)=fnprhc(:)/sum(fnprhc)
                !print*, dr,fnprhc(:)
                
                if (skriv) then 
                    if ( abs(sum(fnprhc(:))-1.0_dp) > eps ) then ; print*, " error in fnprhc: prhc does not add up " , dw , sum(fnprhc(:)) ; stop ; end if 
                end if 
                end function fnprhc

 
                function fnmove(empo,kid,trueindex) 
                    integer(i4b), intent(in) :: empo,kid,trueindex
                    real(dp) :: fnmove
                    integer(i4b) :: c,t,h
                    call index2cotyphome(trueindex,c,t,h)			
                    fnmove = cst(t)  + kcst * one(kid>1) !+ kcst * one(empo==np1)  !kid 1 is no kid, kid 2 is yes kid
                    !ahu october2022: no idea why this was kcst * one(empo==np1) 
                    !fnmove = fnmove / div
                end function fnmove

            !function fnprkid(kid0)
            !integer(i4b), intent(in) :: kid0
            !real(dp), dimension(0:maxkid) :: fnprkid
            !integer(i4b) :: j
            !fnprkid=0.0_dp
            !do j=kid0,maxkid
            !	fnprkid(j)=exp(  pkid * (j-kid0) )
            !	if ( abs(j-kid0) > 1 ) then 
            !		fnprkid(j)=0.0_dp	!can only move one step up or down
            !	end if 
            !end do 						
            !fnprkid(0:maxkid)=fnprkid(0:maxkid)/sum(fnprkid(0:maxkid))
            !if (skriv) then 
            !	if ( abs(sum(fnprkid(0:maxkid))-1.0_dp)>eps ) then ; print*, "error in fnprkid: prkid does not add up " , kid0 , sum(fnprkid(0:maxkid)) ; stop ; end if 
            !end if 
            !end function fnprkid


                subroutine q2wloc(dq,dw,dl)
                    ! extract indeces w,l from q 
                    integer(i4b), intent(in) :: dq		
                    integer(i4b), intent(out) :: dw,dl
                    integer(i4b), dimension(2) :: indeces	
                        indeces=lin2ndim( (/ np2 , nl /) , dq )
                        dw=indeces(1)
                        dl=indeces(2)
                        if (skriv) then  
                            if ( dq > nqs ) then ; print*, "q2wl: q > nqs", dq, nqs,indeces ; stop ; end if  
                            if ( dw > np2 ) then ; print*, "q2wl: w > np2" ; stop ; end if  
                            if ( dl > nl  ) then ; print*, "q2wl: l > nl" ; stop ; end if  
                        end if 
                    end subroutine
                    subroutine wloc2q(dq,dw,dl)
                    ! construct combined q from w,l
                    integer(i4b), intent(out) :: dq		
                    integer(i4b), intent(in) :: dw,dl		
                        dq = ndim2lin( (/ np2 , nl /),(/ dw,dl /) )
                        if (skriv) then 		
                            if ( dq > nqs ) then ; print*, "wl2q: q > nqs" ; stop ; end if  
                            if ( dw > np2 ) then ; print*, "wl2q: w > np2" ; stop ; end if  
                            if ( dl > nl  ) then ; print*, "wl2q: l > nl" ; stop ; end if  
                        end if 
                    end subroutine
                
                    subroutine x2edexpkid(dx,de,dr, dkid)
                    ! extract indeces educ,experience from x
                    integer(i4b), intent(in) :: dx		
                    integer(i4b), intent(out) :: de,dr,dkid
                    integer(i4b), dimension(3) :: indeces	
                        indeces=lin2ndim( (/ neduc, nexp, nkid /) , dx )
                        de=indeces(1)
                        dr=indeces(2)
                        dkid=indeces(3)
                    end subroutine
                    subroutine edexpkid2x(dx,de,dr,dkid)
                    !construct combined x from educ,experience
                    integer(i4b), intent(out) :: dx		
                    integer(i4b), intent(in) :: de,dr,dkid
                        dx=ndim2lin( (/ neduc, nexp, nkid /),(/ de,dr,dkid /) )
                    end subroutine
                
                    subroutine index2cotyphome(index,co,typ,home)
                    ! extract indeces cohort,type,educ,homeloc from combined index
                    integer(i4b), intent(in) :: index		
                    integer(i4b), intent(out) :: co,typ,home	
                    integer(i4b), dimension(3) :: indeces	
                    !if (groups) then 
                        !indeces=lin2ndim((/nco,ntyp,nl/),index)
                        !print*, 'this should not be called if groups!'
                        !stop
                        !co=myco !indeces(1)
                        !typ=mytyp !indeces(2)
                        !home=myhome !indeces(3)
                    !else 
                        indeces=lin2ndim((/ncop,ntypp,nhomep/),index)
                        co=indeces(1)
                        typ=indeces(2)
                        home=indeces(3)
                    !end if 
                    end subroutine index2cotyphome
                    subroutine cotyphome2index(index,co,typ,home)
                    !construct combined index from co,typ,home
                    integer(i4b), intent(out) :: index		! combined index
                    integer(i4b), intent(in) :: co,typ,home 
                    !if (groups) then 
                    !	index=1 !ndim2lin((/nco,ntyp,nl/),(/co,typ,home/))
                    !else 
                        index=ndim2lin((/ncop,ntypp,nhomep/),(/co,typ,home/))
                    !end if 
                    end subroutine cotyphome2index
                 
end module params

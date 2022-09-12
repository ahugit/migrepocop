	if (parcheck) then 
		if (myid==0) open(unit=10,file='par0.txt',status='replace') 
		if (myid==1) open(unit=11,file='par1.txt',status='replace')
		if (myid==2) open(unit=12,file='par2.txt',status='replace')
		if (myid==3) open(unit=13,file='par3.txt',status='replace')
		if (myid==4) open(unit=14,file='par4.txt',status='replace')
		if (myid==5) open(unit=15,file='par5.txt',status='replace')
		if (myid==6) open(unit=16,file='par6.txt',status='replace')
		if (myid==7) open(unit=17,file='par7.txt',status='replace')
		if (myid==8) open(unit=18,file='par8.txt',status='replace')
		if (myid==9) open(unit=19,file='par9.txt',status='replace') 
		if (myid==10) open(unit=20,file='par10.txt',status='replace')
		if (myid==11) open(unit=21,file='par11.txt',status='replace')
		if (myid==12) open(unit=22,file='par12.txt',status='replace')
		if (myid==13) open(unit=23,file='par13.txt',status='replace')
		if (myid==14) open(unit=24,file='par14.txt',status='replace')
		if (myid==15) open(unit=25,file='par15.txt',status='replace')
		if (myid==16) open(unit=26,file='par16.txt',status='replace')
		if (myid==17) open(unit=27,file='par17.txt',status='replace')
	end if 


	if (parcheck) then 
		if (myid==0) then 
			write(10,'("ITER ", I4)') iter	
			do i=1,npars ; write(10,'(F14.5)') parvec(i) ; end do
		else if (myid==1) then 
			write(11,'("ITER ", I4)') iter	
			do i=1,npars ; write(11,'(F14.5)') parvec(i) ; end do
		else if (myid==2) then  
			write(12,'("ITER ", I4)') iter	
			do i=1,npars ; write(12,'(F14.5)') parvec(i) ; end do
		else if (myid==3) then  
			write(13,'("ITER ", I4)') iter	
			do i=1,npars ; write(13,'(F14.5)') parvec(i) ; end do
		else if (myid==4) then  
			write(14,'("ITER ", I4)') iter	
			do i=1,npars ; write(14,'(F14.5)') parvec(i) ; end do
		else if (myid==5) then  
			write(15,'("ITER ", I4)') iter	
			do i=1,npars ; write(15,'(F14.5)') parvec(i) ; end do
		else if (myid==6) then  
			write(16,'("ITER ", I4)') iter	
			do i=1,npars ; write(16,'(F14.5)') parvec(i) ; end do
		else if (myid==7) then  
			write(17,'("ITER ", I4)') iter	
			do i=1,npars ; write(17,'(F14.5)') parvec(i) ; end do
		else if (myid==8) then 
			write(18,'("ITER ", I4)') iter	
			do i=1,npars ; write(18,'(F14.5)') parvec(i) ; end do
		else if (myid==9) then 
			write(19,'("ITER ", I4)') iter	
			do i=1,npars ; write(19,'(F14.5)') parvec(i) ; end do
		else if (myid==10) then  
			write(20,'("ITER ", I4)') iter	
			do i=1,npars ; write(20,'(F14.5)') parvec(i) ; end do
		else if (myid==11) then  
			write(21,'("ITER ", I4)') iter	
			do i=1,npars ; write(21,'(F14.5)') parvec(i) ; end do
		else if (myid==12) then  
			write(22,'("ITER ", I4)') iter	
			do i=1,npars ; write(22,'(F14.5)') parvec(i) ; end do
		else if (myid==13) then  
			write(23,'("ITER ", I4)') iter	
			do i=1,npars ; write(23,'(F14.5)') parvec(i) ; end do
		else if (myid==14) then  
			write(24,'("ITER ", I4)') iter	
			do i=1,npars ; write(24,'(F14.5)') parvec(i) ; end do
		else if (myid==15) then  
			write(25,'("ITER ", I4)') iter	
			do i=1,npars ; write(25,'(F14.5)') parvec(i) ; end do
		else if (myid==16) then 
			write(26,'("ITER ", I4)') iter	
			do i=1,npars ; write(26,'(F14.5)') parvec(i) ; end do
		else if (myid==17) then 
			write(27,'("ITER ", I4)') iter	
			do i=1,npars ; write(27,'(F14.5)') parvec(i) ; end do
		end if 
	end if ! parcheck 


	! get vector of moments from data or simulations
	SUBROUTINE get_moments(dat,ndat,moments,countmom,varmom,momentsname,headstr,headloc,weights)
	INTEGER(I4B), INTENT(IN) :: NDAT ! number of observations in dat array    
	TYPE(DATT), DIMENSION(MNA:MXAI,NDAT), INTENT(IN) :: dat ! data set. first entry is ia index, second observation number
	real(rp), DIMENSION(NMOM), INTENT(OUT) :: moments	 ! conditional moments
	INTEGER(I4B), DIMENSION(NMOM), INTENT(OUT) :: countmom	 ! number of observations contributing to each moment
	real(rp), DIMENSION(NMOM), INTENT(OUT) :: varmom		 ! UNCONDITIONAL variance of each conditional moemnt
	CHARACTER(LEN=namelen), DIMENSION(NMOM), INTENT(OUT) :: momentsname ! names of moments
	CHARACTER(LEN=120),DIMENSION(NMOM), INTENT(OUT) :: headstr ! headers for different sections of moments
	INTEGER(I4B), DIMENSION(NMOM), INTENT(OUT) :: headloc	 ! location of headers for different sections of moments
	real(rp), DIMENSION(NMOM), INTENT(OUT) :: weights	 ! some real assoicated to each moment, maybe used as weights
	INTEGER(I4B),DIMENSION(MNA:MXAI,ndat) :: kidtrans,move,homemove,moverank,NoWorkLast3,NoPartLast3,disolve,norelchg
	INTEGER(I4B),DIMENSION(ndat) :: movesum
	LOGICAL,DIMENSION(MNA:MXAI,ndat) :: ObsLast3,obs4h
	real(rp),DIMENSION(MNA:MXAI,ndat) :: mean4h,mean4hsq,frac4h0,deltawage4,deltawage
	INTEGER(I4B), dimension(MNA:MXAI,ndat) :: iacat,move_single,move_mar,empr
	INTEGER(I4B), dimension(MNA:MXAI,ndat) :: movesum_single,movesum_mar
	INTEGER(I4B), dimension(5,MNA:MXAI,ndat) :: etr
	INTEGER(I4B) :: ia,j,co,im,ihead,g,i
	real(rp), PARAMETER :: miss=-100.0_rp
	weights=0.0_rp 
	moments=miss
	countmom=miss
	varmom=miss
	headloc=miss 
	kidtrans=miss ; move=miss ; movesum=miss ; homemove=miss ; NoWorkLast3=miss ; NoPartLast3=miss
	disolve=miss ; norelchg=miss
	ObsLast3=miss ; obs4h=miss
	mean4h=miss ; mean4hsq=miss ; frac4h0=miss ; deltawage4=miss ; deltawage=miss
	move_single=miss ; move_mar=miss ; movesum_single=miss ; movesum_mar=miss ; etr=miss
	iacat=0
	iacat(18:40,:)=1
	!NOBS=SIZE(dat,2) ! number of observations.
	!deltawage4=-1d0
	! Averia hours worked past four years, averia of square of hours worked last four years, fraction of
	!   last four years out of labor force. Interact these with current wages and/or wage growth
	empr=abs( dat%hhr>=H_FULLTIME )
	DO ia=MNA+4,MXAI
		! are there observations of hours in the past 4 years 
		obs4h(ia,:)=(SUM(ABS(1.*dat(ia-4:ia-1,:)%hhr>=0),1)>0)
		WHERE (obs4h(ia,:))
			mean4h(ia,:)=SUM(empr(ia-4:ia-1,:),1,empr(ia-4:ia-1,:)>=0)/SUM(ABS(1.*empr(ia-4:ia-1,:)>=0),1)
			!mean4hsq(ia,:)=SUM((dat(ia-4:ia-1,:)%hhr)**2,1,dat(ia-4:ia-1,:)%hhr>=0)/SUM(ABS(1.*dat(ia-4:ia-1,:)%hhr>=0),1)
			!frac4h0(ia,:)=SUM(ABS(1.*dat(ia-4:ia-1,:)%hhr==0),1)/SUM(ABS(1.*dat(ia-4:ia-1,:)%hhr>=0),1)
			deltawage4(ia,:)=dat(ia,:)%logwr-dat(ia-4,:)%logwr
		ENDWHERE
	ENDDO  

	!kidtrans=-1.
	WHERE (  (dat(MNA:MXAI-1,:)%rel>0) .AND. (dat(MNA:MXAI-1,:)%kid==0) .AND. (dat(MNA+1:MXAI,:)%kid>=0)  )
		kidtrans(MNA:MXAI-1,:)=ABS(dat(MNA+1:MXAI,:)%kid>0)
	ENDWHERE
	!disolve=-1.
	WHERE ((dat(MNA:MXAI-1,:)%rel>0) .AND. (dat(MNA+1:MXAI,:)%rel>=0))
		disolve(MNA:MXAI-1,:)=ABS(  dat(MNA+1:MXAI,:)%rel==0   .OR.  dat(MNA+1:MXAI,:)%rellen==1  )
	ENDWHERE
	WHERE ((dat(MNA:MXAI-1,:)%rel>=0) .AND. dat(MNA+1:MXAI,:)%rel>=0  )
		norelchg(MNA:MXAI-1,:)=ABS(  dat(MNA:MXAI-1,:)%rel==dat(MNA+1:MXAI,:)%rel   ) 
	ENDWHERE
	!move=-1.
	WHERE ((dat(MNA:MXAI-1,:)%l>0) .AND. (dat(MNA+1:MXAI,:)%l>0))
		move(MNA:MXAI-1,:)=ABS(dat(MNA:MXAI-1,:)%l/=dat(MNA+1:MXAI,:)%l)
	ENDWHERE
	WHERE (move(MNA:MXAI-1,:)==1)
		homemove(MNA:MXAI-1,:)=ABS(dat(MNA+1:MXAI,:)%home==dat(MNA+1:MXAI,:)%l)
	ENDWHERE
	
	movesum(:)=sum(move(MNA:MXAI-1,:),1,move>=0)
	moverank=-1
	do ia=MNA,MXAI-1
		moverank(ia,:)=sum(move(MNA:ia,:),1,move>=0)   
	end do 

	WHERE (move(MNA:MXAI-1,:)>=0.and.dat(MNA:MXAI-1,:)%rel>=0)
		move_single(MNA:MXAI-1,:)=ABS(   move(MNA:MXAI-1,:)==1  .and.  dat(MNA:MXAI-1,:)%rel==0   )
		move_mar(MNA:MXAI-1,:)=ABS(   move(MNA:MXAI-1,:)==1   .and.   dat(MNA:MXAI-1,:)%rel==1  )
       ENDWHERE
	do ia=MNA,MXAI-1
		movesum_single(ia,:)=sum(move_single(MNA:ia,:),1,move_single>=0)
		movesum_mar(ia,:)=sum(move_mar(MNA:ia,:),1,move_mar>=0)
		where (movesum_single(ia,:)>=1)
			movesum_single(ia,:)=1
		elsewhere (movesum_single(ia,:)==0)
			movesum_single(ia,:)=0
		endwhere 
		where (movesum_mar(ia,:)>=1)
			movesum_mar(ia,:)=1
		elsewhere (movesum_mar(ia,:)==0)
			movesum_mar(ia,:)=0
		endwhere 
        end do

	!deltawage=-1.	
	WHERE ((dat(MNA:MXAI-1,:)%logwr>0) .AND. (dat(MNA+1:MXAI,:)%logwr>0))
		deltawage(MNA:MXAI-1,:)=dat(MNA+1:MXAI,:)%logwr-dat(MNA:MXAI-1,:)%logwr
	ENDWHERE
	!etr=-1.
	!ahu 062512 given the new definition of hours (was discrete and now it's continuous, changing the definition of etr) 
	!  <FULLTIME TO >=FULLTIME
	WHERE ( (dat(MNA:MXAI-1,:)%hhr<H_FULLTIME).AND.(dat(MNA:MXAI-1,:)%hhr>=0) .AND. (dat(MNA+1:MXAI,:)%hhr>=0) ) !ahu 062512 check this -1 
		etr(1,MNA:MXAI-1,:)=ABS(dat(MNA+1:MXAI,:)%hhr>=H_FULLTIME)
	ENDWHERE
	WHERE ( (dat(MNA:MXAI-1,:)%hhr>=H_FULLTIME).AND.(dat(MNA:MXAI-1,:)%hhr>=0) .AND. (dat(MNA+1:MXAI,:)%hhr>=0) ) !ahu 062512 check this -1 
		etr(2,MNA:MXAI-1,:)=ABS(dat(MNA+1:MXAI,:)%hhr>=H_FULLTIME)
	ENDWHERE
	!  <PARTTIME TO >=PARTTIME
	WHERE ( (dat(MNA:MXAI-1,:)%hhr<H_PARTTIME).AND.(dat(MNA:MXAI-1,:)%hhr>=0) .AND. (dat(MNA+1:MXAI,:)%hhr>=0) ) !ahu 062512 check this -1 
		etr(3,MNA:MXAI-1,:)=ABS(dat(MNA+1:MXAI,:)%hhr>=H_PARTTIME)
	ENDWHERE
	WHERE ( (dat(MNA:MXAI-1,:)%hhr>=H_PARTTIME).AND.(dat(MNA:MXAI-1,:)%hhr>=0) .AND. (dat(MNA+1:MXAI,:)%hhr>=0) ) !ahu 062512 check this -1 
		etr(4,MNA:MXAI-1,:)=ABS(dat(MNA+1:MXAI,:)%hhr>=H_PARTTIME)
	ENDWHERE
	!  <FULLTIME TO <FULLTIME
	WHERE ( (dat(MNA:MXAI-1,:)%hhr<H_FULLTIME).AND.(dat(MNA:MXAI-1,:)%hhr>=0) .AND. (dat(MNA+1:MXAI,:)%hhr>=0) ) !ahu 062512 check this -1 
		etr(5,MNA:MXAI-1,:)=ABS(dat(MNA+1:MXAI,:)%hhr<H_FULLTIME)
	ENDWHERE

	
	if (yaz) then
	if (ndat==nsim) then
		DO j=1,nsim
			! write(12,*) ' id  ia sex rel  ed edp kid   hh  wage    dd incsp ln je cm sc sm kt ds O3'
			DO ia=MNA,MXAI
				write(12,'(7I12,6F25.5,7I12,2F25.5,9I12)') j,ia,dat(ia,j)%co,dat(ia,j)%sexr,&
				& dat(ia,j)%rel,dat(ia,j)%kid,dat(ia,j)%edr, &
				& dat(ia,j)%logwr,dat(ia,j)%logwsp,dat(ia,j)%hhr, & 
				& dat(ia,j)%hhsp,dat(ia,j)%ddr,dat(ia,j)%ddsp,&
				& dat(ia,j)%rellen,-1,kidtrans(ia,j),& 
				& disolve(ia,j),ObsLast3(ia,j),dat(ia,j)%l,& 
				& move(ia,j),deltawage(ia,j),deltawage4(ia,j),etr(1:5,ia,j),moverank(ia,j),homemove(ia,j),& 
				& norelchg(ia,j)
			ENDDO
		ENDDO 
	end if 
	end if 

	! set weights equal to one but can adjust weight on any particular moment by setting weights(im)=.. below
	! SL052311: NOTE THAT DEFAULT WEIGHTS ARE ZERO
	weights=0.
	im=1
	ihead=1	
	do co=1,NCO
		headloc(ihead)=im
		IF (co==1) THEN
			headstr(ihead)='br 1';ihead=ihead+1
		ELSE
			headstr(ihead)='br 2';ihead=ihead+1
		ENDIF	
		CALL condmom(im,((dat%co==co).AND.(dat%rel>=0)),abs(d1*(dat%rel==1)),moments,countmom,varmom)
		WRITE(momentsname(im),'("frac-married/all")') 
		weights(im)=wrel
		im=im+1
		CALL condmom(im,((dat%co==co).AND.(dat%rel>=0).AND.((dat%edr==1))),abs(d1*(dat%rel==1)),moments,countmom,varmom)
		WRITE(momentsname(im),'("frac-married/hs")') 
		weights(im)=wrel
		im=im+1
		CALL condmom(im,((dat%co==co).AND.(dat%rel>=0).AND.((dat%edr==2))),abs(d1*(dat%rel==1)),moments,countmom,varmom)
		WRITE(momentsname(im),'("frac-married/col")') 
		weights(im)=wrel
		im=im+1

		headloc(ihead)=im; headstr(ihead)='Marriia rates by ia ';ihead=ihead+1	
		DO ia=19,MXAI,4
			CALL condmom(im,((dat(ia,:)%co==co).AND.(dat(ia,:)%rel>=0)),abs(d1*(dat(ia,:)%rel==1)),moments,countmom,varmom)
			WRITE(momentsname(im),'("frac-married ",I4)') ia 
			weights(im)=wrel
			im=im+1
		ENDDO 

		headloc(ihead)=im; headstr(ihead)='Disolution Rates ';ihead=ihead+1	
		CALL condmom(im,((dat%co==co).AND.(disolve>=0)),d1*disolve,moments,countmom,varmom)
		WRITE(momentsname(im),'("disolve         ")') 
		weights(im)=wrel
		im=im+1
		CALL condmom(im,((dat%co==co).AND.(disolve>=0).and.(dat%kid==0)),d1*disolve,moments,countmom,varmom)
		WRITE(momentsname(im),'("disolve | nokid ")') 
		weights(im)=wrel
		im=im+1
		CALL condmom(im,((dat%co==co).AND.(disolve>=0).and.(dat%kid>0)),d1*disolve,moments,countmom,varmom)
		WRITE(momentsname(im),'("disolve | kid   ")') 
		weights(im)=wrel
		im=im+1
		CALL condmom(im,((dat%co==co).AND.(disolve>=0).and.(dat%hhr<H_FULLTIME)),d1*disolve,moments,countmom,varmom)
		WRITE(momentsname(im),'("disolve | u     ")') 
		weights(im)=wrel
		im=im+1
		CALL condmom(im,((dat%co==co).AND.(disolve>=0).and.(dat%hhr>=H_FULLTIME)),d1*disolve,moments,countmom,varmom)
		WRITE(momentsname(im),'("disolve | e     ")') 
		weights(im)=wrel
		im=im+1
		DO j=1,3
			CALL condmom(im,((dat%co==co).AND.(dat%rellen>3*(j-1)).AND.(dat%rellen<=3*j).AND.(disolve>=0)),d1*disolve,moments,countmom,varmom)
			WRITE(momentsname(im),'("disolve/dur ",2I4)') 3*(j-1)+1,3*j
			weights(im)=wrel
			im=im+1
		ENDDO

		!headloc(ihead)=im; headstr(ihead)='Move by total number of moves (all ages)';ihead=ihead+1
		!movesum(:)=sum(move(MNA:MXAI-1,:),1,move>=0)  
		!do j=0,4
		!	CALL condmom(im,((dat%co==co)),d1*move(MNA:MXAI-1,:),moments,countmom,varmom)
			!CALL condmom(im,((movesum>=0)),d1*abs(movesum==j),moments,countmom,varmom)
		!	WRITE(momentsname(im),'("movesum ",I4)') j
		!	weights(im)=wmove
		!	im=im+1 
		!end do 
		!do j=0,4
		!	CALL condmom(im,((movesum>=0)),d1*abs(movesum==j),moments,countmom,varmom)
		!	WRITE(momentsname(im),'("movesum ",I4)') j
		!	weights(im)=wmove
		!	im=im+1 
		!end do 


		headloc(ihead)=im; headstr(ihead)='After the first move what proportion is return to home? (all ages)';ihead=ihead+1		
		CALL condmom(im,((moverank>1).AND.(move==1)  ),d1*abs(homemove==0),moments,countmom,varmom)
		WRITE(momentsname(im),'("non-home ")') 
		weights(im)=wmove0
		im=im+1 
		CALL condmom(im,((moverank>1).AND.(move==1)  ),d1*abs(homemove==1),moments,countmom,varmom)
		WRITE(momentsname(im),'("home     ")') 
		weights(im)=wmove0
		im=im+1 

		headloc(ihead)=im; headstr(ihead)='Move Rates ';ihead=ihead+1
		CALL condmom(im,((dat%co==co).AND.(dat%rel==0).and.(norelchg==1).AND.(move>=0)),d1*move,moments,countmom,varmom)
		WRITE(momentsname(im),'("move/single")') 
		weights(im)=wmove0 	    		
		im=im+1 
		CALL condmom(im,((dat%co==co).AND.(dat%rel==1).and.(norelchg==1).AND.(move>=0)),d1*move,moments,countmom,varmom)
		WRITE(momentsname(im),'("move/married")') 
		weights(im)=wmove1 	    		
		im=im+1 
		do ia=19,MXAI,4
			CALL condmom(im,((dat(ia,:)%co==co).AND.(dat(ia,:)%rel==0).and.(norelchg(ia,:)==1).AND.(move(ia,:)>=0)),d1*move(ia,:),moments,countmom,varmom)
			WRITE(momentsname(im),'("move/single ",I4)') ia
			weights(im)=wmove0 	    		
			im=im+1 
		end do 

		do ia=19,MXAI,4
			CALL condmom(im,((dat(ia,:)%co==co).AND.(dat(ia,:)%rel==1).and.(norelchg(ia,:)==1).AND.(move(ia,:)>=0)),d1*move(ia,:),moments,countmom,varmom)
			WRITE(momentsname(im),'("move/married ",I4)') ia
			weights(im)=wmove1 	    		
			im=im+1 ! married
		end do 

		headloc(ihead)=im; headstr(ihead)='Move by kid (all ages) ';ihead=ihead+1
		CALL condmom(im,((dat%co==co).AND.(dat%rel==1).and.(norelchg==1).and.(dat%kid==0).AND.(move>=0)),d1*move,moments,countmom,varmom)
		WRITE(momentsname(im),'("move/nokid")') 
		weights(im)=wmove1  		
		im=im+1 
		CALL condmom(im,((dat%co==co).AND.(dat%rel==1).and.(norelchg==1).and.(dat%kid>0).AND.(move>=0)),d1*move,moments,countmom,varmom)
		WRITE(momentsname(im),'("move/kid")') 
		weights(im)=wmove1
		im=im+1 

		headloc(ihead)=im; headstr(ihead)='Employment by gender/rel (all ages) ';ihead=ihead+1
		do g=1,2
			do j=0,1
				CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==j).AND.(dat%hhr>=0).AND.(iacat==1)),d1*abs(dat%hhr>=H_FULLTIME),moments,countmom,varmom)
				WRITE(momentsname(im),'("emp | gender/rel ",2I4)') g,j
				if (j==0) weights(im)=whour0 
				if (j==1) weights(im)=whour1 
				im=im+1
			end do 
		end do 
		headloc(ihead)=im; headstr(ihead)='Employment by gender/kid (married and all age)';ihead=ihead+1
		do g=1,2
			do j=0,1
				CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==1).AND.(dat%kid==j).AND.(dat%hhr>=0).AND.(iacat==1)),d1*abs(dat%hhr>=H_FULLTIME),moments,countmom,varmom)
				WRITE(momentsname(im),'("emp | gender/kid ",2I4)') g,j
				weights(im)=whour1 
				im=im+1 
			end do 
		end do 
		headloc(ihead)=im; headstr(ihead)='Move by whether they are working full time and gender (all ages) ';ihead=ihead+1
		do g=1,2
			CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).and.(norelchg==1).AND.(dat%hhr>=0).AND.(dat%hhr<H_FULLTIME).AND.(move>=0)),d1*move,moments,countmom,varmom)
			WRITE(momentsname(im),'("move|u ",I4)') g
			weights(im)=wmove0
			im=im+1
			CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).and.(norelchg==1).AND.(dat%hhr>=0).AND.(dat%hhr>=H_FULLTIME).AND.(move>=0)),d1*move,moments,countmom,varmom)
			WRITE(momentsname(im),'("move|e ",I4)') g
			weights(im)=wmove0
			im=im+1
		end do 
		headloc(ihead)=im; headstr(ihead)='Move by whether they are working full time and gender and rel (all ages) ';ihead=ihead+1
		do g=1,2
			do j=0,1
				CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==j).and.(norelchg==1).AND.(dat%hhr>=0).AND.(dat%hhr<H_FULLTIME).AND.(move>=0)),d1*move,moments,countmom,varmom)
				WRITE(momentsname(im),'("move|u ",2I4)') g,j
				if (j==0) weights(im)=wmove0 
				if (j==1) weights(im)=wmove1 
				im=im+1
				CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==j).and.(norelchg==1).AND.(dat%hhr>=0).AND.(dat%hhr>=H_FULLTIME).AND.(move>=0)),d1*move,moments,countmom,varmom)
				WRITE(momentsname(im),'("move|e ",2I4)') g,j
				if (j==0) weights(im)=wmove0 
				if (j==1) weights(im)=wmove1 
				im=im+1
			end do 
		end do 

		headloc(ihead)=im; headstr(ihead)='U2EMP by gender/rel - FULL TIME (all ages)';ihead=ihead+1
		do g=1,2
			do j=0,1
				CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==j).and.(norelchg==1).AND.(etr(1,:,:)>=0).AND.(iacat==1)),d1*etr(1,:,:),moments,countmom,varmom)
				WRITE(momentsname(im),'("u2emp          ",2I4)') g,j
				if (j==0) weights(im)=whour0 
				if (j==1) weights(im)=whour1 
				im=im+1 
				CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==j).and.(norelchg==1).AND.(etr(1,:,:)>=0).AND.(move==0).AND.(iacat==1)),d1*etr(1,:,:),moments,countmom,varmom)
				WRITE(momentsname(im),'("u2emp | stay   ",2I4)') g,j
				if (j==0) weights(im)=whour0 
				if (j==1) weights(im)=whour1 
				im=im+1 
			end do 
		end do 

		headloc(ihead)=im; headstr(ihead)='EMP2EMP by gender/rel - FULL TIME (all ages)';ihead=ihead+1
		do g=1,2
			do j=0,1
				CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==j).and.(norelchg==1).AND.(etr(2,:,:)>=0).AND.(iacat==1)),d1*etr(2,:,:),moments,countmom,varmom)
				WRITE(momentsname(im),'("emp2emp        ",2I4)') g,j
				if (j==0) weights(im)=whour0 
				if (j==1) weights(im)=whour1 
				im=im+1 
				CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==j).and.(norelchg==1).AND.(etr(2,:,:)>=0).AND.(move==0).AND.(iacat==1)),d1*etr(2,:,:),moments,countmom,varmom)
				WRITE(momentsname(im),'("emp2emp | stay ",2I4)') g,j
				if (j==0) weights(im)=whour0 
				if (j==1) weights(im)=whour1 
				im=im+1 
				!CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==j).and.(norelchg==1).AND.(etr(2,:,:)>=0).AND.(move==1).AND.(iacat==1)),d1*etr(2,:,:),moments,countmom,varmom)
				!WRITE(momentsname(im),'("emp2emp | move ",2I4)') g,j
				!weights(im)=whour
				!im=im+1 		
			end do 
		end do 
		headloc(ihead)=im; headstr(ihead)='Mean wage by gender/loc - FULL TIME (all ages)';ihead=ihead+1
		do g=1,2
			do j=1,NL
				CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%l==j).AND.(dat%hhr>=H_FULLTIME).AND.(dat%logwr>=0).AND.(iacat==1)),d1*dat%logwr,moments,countmom,varmom)
				WRITE(momentsname(im),'("logwage/gender/loc ",2I4)') g,j
				weights(im)=wwage0
				im=im+1
			end do 
		end do 

		headloc(ihead)=im; headstr(ihead)='Mean wage by gender/loc/ed - FULL TIME (all ages)';ihead=ihead+1
		do g=1,2
			do j=1,NL
				CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%edr==1).AND.(dat%l==j).AND.(dat%hhr>=H_FULLTIME).AND.(dat%logwr>=0).AND.(iacat==1)),d1*dat%logwr,moments,countmom,varmom)
				WRITE(momentsname(im),'("logwage/gender/loc | hs ",2I4)') g,j
				weights(im)=wwage0
				im=im+1
				CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%edr==2).AND.(dat%l==j).AND.(dat%hhr>=H_FULLTIME).AND.(dat%logwr>=0).AND.(iacat==1)),d1*dat%logwr,moments,countmom,varmom)
				WRITE(momentsname(im),'("logwage/gender/loc | col  ",2I4)') g,j
				weights(im)=wwage0
				im=im+1
			end do 
		end do 


		headloc(ihead)=im; headstr(ihead)='Mean wage by gender/ed - FULL TIME (all ages)';ihead=ihead+1
		do g=1,2
			do j=1,2	
				CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%edr==j).AND.(dat%hhr>=H_FULLTIME).AND.(dat%logwr>=0).AND.(iacat==1)),d1*dat%logwr,moments,countmom,varmom)
				WRITE(momentsname(im),'("logw/gender/ed ",2I4)') g,j
				weights(im)=wwage0
				im=im+1
			end do 
		end do 

		headloc(ihead)=im; headstr(ihead)='Mean wage by gender/ed/ia - FULL TIME (all ages)';ihead=ihead+1
		do g=1,2
			do j=1,2	
				do ia=MNA+1,MXAI-4,4
					CALL condmom(im,((dat(ia,:)%co==co).AND.(dat(ia,:)%sexr==g).AND.(dat(ia,:)%edr==j).AND.(dat(ia,:)%hhr>=H_FULLTIME).AND.(dat(ia,:)%logwr>=0)),d1*dat(ia,:)%logwr,moments,countmom,varmom)
					WRITE(momentsname(im),'("logw/gender/ed/ia",3I4)') g,j,ia
					weights(im)=wwage0
					im=im+1
				end do 
			end do 
		end do 

		headloc(ihead)=im; headstr(ihead)='wage-growth for stay and move by gender/rel (all ages)';ihead=ihead+1
		do g=1,2
			do j=0,1	
				CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==j).AND.(move==0).AND.(dat%logwr>=0).AND.(deltawage>miss).AND.(iacat==1)),d1*deltawage,moments,countmom,varmom)
				WRITE(momentsname(im),'("wage-growth/stay/gender/rel ",2I4)') g,j
				if (j==0) weights(im)=wwage0 
				if (j==1) weights(im)=wwage1 
				im=im+1
				CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==j).AND.(move==1).AND.(dat%logwr>=0).AND.(deltawage>miss).AND.(iacat==1)),d1*deltawage,moments,countmom,varmom)
				WRITE(momentsname(im),'("wage-growth/move/gender/rel ",2I4)') g,j
				if (j==0) weights(im)=wwage0 
				if (j==1) weights(im)=wwage1 
				im=im+1
			end do 
		end do 
			
		headloc(ihead)=im; headstr(ihead)='education by location';ihead=ihead+1
		do j=1,NL
			CALL condmom(im,((dat%co==co).AND.(dat%l==j).AND.(dat%edr>=0)),abs(d1*(dat%edr==2)),moments,countmom,varmom)
			WRITE(momentsname(im),'("prop-col ",I4)') j
			weights(im)=wwage0
			im=im+1
		end do 

		headloc(ihead)=im; headstr(ihead)='mean log wage by gender/ia - FULL TIME';ihead=ihead+1
		do g=1,2
		DO ia=MNA+1,MXAI,4
			CALL condmom(im,((dat(ia,:)%co==co).AND.(dat(ia,:)%sexr==g).AND.(dat(ia,:)%hhr>=H_FULLTIME).AND.(dat(ia,:)%logwr>=0)),d1*dat(ia,:)%logwr,moments,countmom,varmom)
			WRITE(momentsname(im),'("mean-log-wage/gender/ia ",2I4)') g,ia
			weights(im)=wwage0
			im=im+1
		ENDDO
		end do 

		headloc(ihead)=im; headstr(ihead)='mean log wage and log wage sq by gender and by hours';ihead=ihead+1
		do g=1,2
			CALL condmom(im,((dat%co==co).AND.(dat%sexr==g)	.AND.(dat%hhr>=H_FULLTIME).AND.(dat%logwr>=0).AND.(iacat==1)),d1*(dat%logwr),moments,countmom,varmom)
			WRITE(momentsname(im),'("mean-log-wage/full-time/gender ",I4)') g
			weights(im)=wwage0
			im=im+1
			CALL condmom(im,((dat%co==co).AND.(dat%sexr==g)	.AND.(dat%hhr>=H_FULLTIME).AND.(dat%logwr>=0).AND.(iacat==1)),d1*(dat%logwr)**2,moments,countmom,varmom)
			WRITE(momentsname(im),'("mean-log-wgsq/full-time/gender ",I4)') g
			weights(im)=wwage0
			im=im+1
		end do 

		headloc(ihead)=im; headstr(ihead)='wages by gender/rel - FULL TIME';ihead=ihead+1
		do g=1,2
			do j=0,1
				CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==j) .AND.(dat%hhr>=H_FULLTIME).AND.(dat%logwr>=0).AND.(iacat==1)),d1*dat%logwr,moments,countmom,varmom)
				WRITE(momentsname(im),'("logwage/gender/rel ",2I4)') g,j
				if (j==0) weights(im)=wwage0 
				if (j==1) weights(im)=wwage1 
				im=im+1
			end do 
		end do 

		headloc(ihead)=im; headstr(ihead)='log wage spXr - FULL TIME ';ihead=ihead+1
		CALL condmom(im,((dat%co==co).AND.(dat%rel==1).AND.(dat%hhr>=H_FULLTIME).AND.(dat%hhsp>=H_FULLTIME).AND.(dat%logwr>=0).AND.(dat%logwsp>=0).AND.(iacat==1)),d1*dat%logwr*dat%logwsp,moments,countmom,varmom)
		WRITE(momentsname(im),'("log-wage-spXr")') 
		weights(im)=wwage1
		im=im+1

		! wages-growth interacted with hours worked and hours squared by gender
		headloc(ihead)=im; headstr(ihead)='wages-growth, interacted with hours worked and hours squared by gender';ihead=ihead+1
		CALL condmom(im,((dat(22:MXAI,:)%co==co).AND.(dat(22:MXAI,:)%sexr==1).AND.(dat(22:MXAI,:)%logwr>=0).AND.(dat(18:MXAI-4,:)%logwr>=0).AND.(deltawage4(22:MXAI,:)>=0)),d1*deltawage4(22:MXAI,:),moments,countmom,varmom)
		WRITE(momentsname(im),'("wage-growth/male             ")') 
		CALL condmom(im+1,((dat(22:MXAI,:)%co==co).AND.(dat(22:MXAI,:)%sexr==1).AND.(dat(22:MXAI,:)%logwr>=0).AND.(dat(18:MXAI-4,:)%logwr>=0).AND.(deltawage4(22:MXAI,:)>=0).AND.(mean4h(22:MXAI,:)>=0)),d1*deltawage4(22:MXAI,:)*mean4h(22:MXAI,:),moments,countmom,varmom)
		WRITE(momentsname(im+1),'("wage-growth-hours/male     ")') 
		CALL condmom(im+2,((dat(22:MXAI,:)%co==co).AND.(dat(22:MXAI,:)%sexr==2).AND.(dat(22:MXAI,:)%logwr>=0).AND.(dat(18:MXAI-4,:)%logwr>=0).AND.(deltawage4(22:MXAI,:)>=0)),d1*deltawage4(22:MXAI,:),moments,countmom,varmom)
		WRITE(momentsname(im+2),'("wage-growth/female         ")') 
		CALL condmom(im+3,((dat(22:MXAI,:)%co==co).AND.(dat(22:MXAI,:)%sexr==2).AND.(dat(22:MXAI,:)%logwr>=0).AND.(dat(18:MXAI-4,:)%logwr>=0).AND.(deltawage4(22:MXAI,:)>=0).AND.(mean4h(22:MXAI,:)>=0)),d1*deltawage4(22:MXAI,:)*mean4h(22:MXAI,:),moments,countmom,varmom)
		WRITE(momentsname(im+3),'("wage-growth-hours/female   ")') 
		weights(im:im+3)=wwage0
		im=im+4


		headloc(ihead)=im; headstr(ihead)='Prop by loc and Prop of moves by loc (all age) ';ihead=ihead+1
		!ahu 082012 do j=1,NL
		!ahu 082012 CALL condmom(im,((dat(MNA:MXAI-1,:)%co==co).AND.(dat(MNA:MXAI-1,:)%l==j).AND.(move(MNA:MXAI-1,:)>=0)),d1*move(MNA:MXAI-1,:),moments,countmom,varmom)
		!ahu 082012 WRITE(momentsname(im),'("move-rates-from-loc",3I4)') 0,0,j
		!ahu 082012 weights(im)=wmove
		!ahu 082012 im=im+1 
		!ahu 082012 end do 
		!ahu 082012 do j=1,NL
		!ahu 082012 	CALL condmom(im,((dat(MNA:MXAI-1,:)%co==co).AND.(dat(MNA+1:MXAI,:)%l==j).AND.(move(MNA:MXAI-1,:)>=0)),d1*move(MNA:MXAI-1,:),moments,countmom,varmom)
		!ahu 082012 	WRITE(momentsname(im),'("move-rates-to-loc",3I4)') 0,0,j
		!ahu 082012 	weights(im)=wmove
		!ahu 082012 	im=im+1 
		!ahu 082012 end do 
		do j=1,NL
			CALL condmom(im,((dat%co==co)),d1*abs(dat%l==j),moments,countmom,varmom)
			WRITE(momentsname(im),'("prop-loc                  ",I4)') j
			weights(im)=wmove0
			im=im+1 
		end do 	
		do j=1,NL
			CALL condmom(im,((dat(MNA:MXAI-1,:)%co==co).AND.(move(MNA:MXAI-1,:)==1).AND.(dat(MNA:MXAI-1,:)%l>=0)),d1*abs(dat(MNA:MXAI-1,:)%l==j),moments,countmom,varmom)
			WRITE(momentsname(im),'("prop-of-moves-from        ",I4)') j
			weights(im)=wmove0
			im=im+1 
		end do 	
		do j=1,NL
			CALL condmom(im,((dat(MNA:MXAI-1,:)%co==co).AND.(move(MNA:MXAI-1,:)==1).AND.(dat(MNA+1:MXAI,:)%l>=0)),d1*abs(dat(MNA+1:MXAI,:)%l==j),moments,countmom,varmom)
			WRITE(momentsname(im),'("prop-of-moves-to          ",I4)') j
			weights(im)=wmove0
			im=im+1 
		end do 	
		do j=1,NL
			CALL condmom(im,((dat(MNA:MXAI-1,:)%co==co).AND.(dat(MNA+1:MXAI,:)%l==j).AND.(homemove(MNA:MXAI-1,:)>=0)),d1*abs(homemove(MNA:MXAI-1,:)==0),moments,countmom,varmom)
			WRITE(momentsname(im),'("prop-of-non-home-moves-to ",I4)') j
			weights(im)=wmove0
			im=im+1 
			CALL condmom(im,((dat(MNA:MXAI-1,:)%co==co).AND.(dat(MNA+1:MXAI,:)%l==j).AND.(homemove(MNA:MXAI-1,:)>=0)),d1*abs(homemove(MNA:MXAI-1,:)==1),moments,countmom,varmom)
			WRITE(momentsname(im),'("prop-of-home-moves-to     ",I4)') j
			weights(im)=wmove0
			im=im+1 
		end do 	
		do i=1,NL
			do j=1,NL
				CALL condmom(im,((dat(MNA:MXAI-1,:)%co==co).AND.(move(MNA:MXAI-1,:)==1).AND.(dat(MNA:MXAI-1,:)%l==i)),d1*abs(dat(MNA+1:MXAI,:)%l==j),moments,countmom,varmom)
				WRITE(momentsname(im),'("matloc                    ",2I4)') j,i
				weights(im)=wmove0
				im=im+1 
			end do 	
		end do 
		do i=1,NL
		CALL condmom(im,((dat%co==co).AND.(move>=0).AND.(dat%l==i)),d1*abs(move==1),moments,countmom,varmom)
			WRITE(momentsname(im),'("move | loc                ",I4)') i
			weights(im)=wmove0
			im=im+1 
		end do 
		!headloc(ihead)=im; headstr(ihead)='Labor market hours by gender/rel/ia';ihead=ihead+1
		!ahu 061211: have to control for ia here because the two brs have different ia compositions
		!ahu 061211: br 2 has no hours/kids/cohmar simultaneously in the biannual years so if you condition on all that you will just get something until they are ia 28 or something (depending on what the br grouping is)
		!ahu 061211: and so if we don't control for ia, it looks as if br 2 females who are cohabiting have decreased their hours of work. but this is just a composition effect.
		!ahu 061211: excluding ia 20 because, something looks weird. br 2 works too few hours at ia 20 (for females,coh,nokid). so then when I include them, it looks as if br 2 coh females with no kids work less in the later br. 
		!do g=1,2
		!	do j=0,1
		!		CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==j).AND.(dat%hhr>=0).AND.(iacat==1)),dat%hhr,moments,countmom,varmom)
		!		WRITE(momentsname(im),'("hrs/gender/rel ",2I4)') g,j
		!		weights(im)=whour
		!		im=im+1 
		!	end do 
		!end do 
		
		headloc(ihead)=im; headstr(ihead)='Kids';ihead=ihead+1
		CALL condmom(im,((dat%co==co).AND.(dat%rel==1).and.(norelchg==1).AND.(kidtrans>=0).AND.(iacat==1)),d1*kidtrans,moments,countmom,varmom)
		WRITE(momentsname(im),'("kidtrans-married ")') 
		weights(im)=wkid
		im=im+1
		CALL condmom(im,((dat%co==co).AND.(dat%rel==1).and.(norelchg==1).AND.(dat%kid>=0).AND.(iacat==1)),d1*dat%kid,moments,countmom,varmom)
		WRITE(momentsname(im),'("kids-married     ")') 
		weights(im)= wkid 
		im=im+1
	ENDDO ! end br loop
	!print*, "Here is NMOM ", im-1
	
	!IF YOU ADD A MOMENT DON'T FORGET TO ADD IM=IM+1 TO THE LAST MOMENT ABOVE
	!!if (myid==0) then ; print*, 'Here is im in get_moments',im-1 ; end if 
	END SUBROUTINE get_moments


	SUBROUTINE get_extramoments(dat,ndat,moments,countmom,varmom,momentsname,headstr,headloc,weights)
	! Do other moments as well for display purposes- need to be included in total count of NMOM
	INTEGER(I4B), INTENT(IN) :: NDAT					! number of observations in dat array    
	TYPE(DATT), DIMENSION(MNA:MXAI,NDAT), INTENT(IN) :: dat		! data set. first entry is ia index, second observation number
	real(rp), DIMENSION(NMOM), INTENT(OUT) :: moments			! conditional moments
	INTEGER(I4B), DIMENSION(NMOM), INTENT(OUT) :: countmom		! number of observations contributing to each moment
	real(rp), DIMENSION(NMOM), INTENT(OUT) :: varmom			! UNCONDITIONAL variance of each conditional moemnt
	CHARACTER(LEN=namelen), DIMENSION(NMOM), INTENT(OUT) :: momentsname ! names of moments
	CHARACTER(LEN=120),DIMENSION(NMOM), INTENT(OUT) :: headstr		! headers for different sections of moments
	INTEGER(I4B), DIMENSION(NMOM), INTENT(OUT) :: headloc		! location of headers for different sections of moments
	real(rp), DIMENSION(NMOM), INTENT(OUT) :: weights			! some real assoicated to each moment, maybe used as weights
	INTEGER(I4B),DIMENSION(MNA:MXAI,ndat) :: move,norelchg
	INTEGER(I4B), dimension(MNA:MXAI,ndat) :: movesum_single,movesum_mar,iacat
	INTEGER(I4B), dimension(5,MNA:MXAI,ndat) :: etr
	INTEGER(I4B) :: g,j,ia,im,co,ihead

	weights=0.0_rp
	im=1
	ihead=1	
	cohort: DO co=1,NCO 
		headloc(ihead)=im
		IF (co==1) THEN
			headstr(ihead)='Cohort 1: Other Moments';ihead=ihead+1
		ELSE
			headstr(ihead)='Cohort 2: Other Moments';ihead=ihead+1
		ENDIF
		headloc(ihead)=im; headstr(ihead)='mean log wage by gender/movesum_rel/ia (FULL TIME)';ihead=ihead+1		
		do g=1,2
			do j=0,1
				do ia=MNA+1,MXAI-4,4
					CALL condmom(im,((dat(ia,:)%co==co).AND.(dat(ia,:)%sexr==g).AND.(dat(ia,:)%rel==0).AND.(movesum_single(ia,:)==j).AND.(dat(ia,:)%hhr>=H_FULLTIME).AND.(dat(ia,:)%logwr>=0)),d1*dat(ia,:)%logwr,moments,countmom,varmom)
					WRITE(momentsname(im),'("mean-log-wage-by-gender/movesing/ia ",3I4)') g,j,ia
					weights(im)=0.0_rp
					im=im+1
				end do 
			end do 
		end do 
		do g=1,2
			do j=0,1
				do ia=MNA+1,MXAI-4,4
					CALL condmom(im,((dat(ia,:)%co==co).AND.(dat(ia,:)%sexr==g).AND.(dat(ia,:)%rel==1).AND.(movesum_mar(ia,:)==j).AND.(dat(ia,:)%hhr>=H_FULLTIME).AND.(dat(ia,:)%logwr>=0)),d1*dat(ia,:)%logwr,moments,countmom,varmom)
					WRITE(momentsname(im),'("mean-log-wage-by-gender/movemar/ia  ",3I4)') g,j,ia
					weights(im)=0.0_rp
					im=im+1
				end do 
			end do 
		end do 
		headloc(ihead)=im; headstr(ihead)='U2EMP by gender/home - SINGLES - FULL TIME (all ages)';ihead=ihead+1
		do g=1,2
			CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==0).and.(norelchg==1).AND.(move==0).AND.(etr(1,:,:)>=0).AND.(iacat==1)),d1*etr(1,:,:),moments,countmom,varmom)
			WRITE(momentsname(im),'("u2emp | stay          ",I2)') g
			weights(im)=0.0_rp
			im=im+1 
			CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==0).and.(norelchg==1).AND.(move==1).AND.(dat%home==dat%l).AND.(etr(5,:,:)>=0).AND.(iacat==1)),d1*etr(5,:,:),moments,countmom,varmom)
			WRITE(momentsname(im),'("u2u | move,athme1    ",I2)') g
			weights(im)=0.0_rp
			im=im+1 
			CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==0).and.(norelchg==1).AND.(move==1).AND.(dat%home==dat%l).AND.(etr(1,:,:)>=0).AND.(iacat==1)),d1*etr(1,:,:),moments,countmom,varmom)
			WRITE(momentsname(im),'("u2emp | move,athme1  ",I2)') g
			weights(im)=0.0_rp
			im=im+1 
			CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==0).and.(norelchg==1).AND.(move==1).AND.(dat%home/=dat%l).AND.(etr(1,:,:)>=0).AND.(iacat==1)),d1*etr(1,:,:),moments,countmom,varmom)
			WRITE(momentsname(im),'("u2emp | move,athme0  ",I2)') g
			weights(im)=0.0_rp
			im=im+1 
		end do 


		headloc(ihead)=im; headstr(ihead)='EMP2EMP by gender/home - SINGLES - FULL TIME (all ages)';ihead=ihead+1
		do g=1,2
			CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==0).and.(norelchg==1).AND.(move==0).AND.(etr(2,:,:)>=0).AND.(iacat==1)),d1*etr(2,:,:),moments,countmom,varmom)
			WRITE(momentsname(im),'("emp2emp | stay        ",I2)') g
			weights(im)=0.0_rp
			im=im+1 
			CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==0).and.(norelchg==1).AND.(move==1).AND.(dat%home==dat%l).AND.(etr(2,:,:)>=0).AND.(iacat==1)),d1*etr(2,:,:),moments,countmom,varmom)
			WRITE(momentsname(im),'("emp2emp | move,athme1",I2)') g
			weights(im)=0.0_rp
			im=im+1 
			CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==0).and.(norelchg==1).AND.(move==1).AND.(dat%home/=dat%l).AND.(etr(2,:,:)>=0).AND.(iacat==1)),d1*etr(2,:,:),moments,countmom,varmom)
			WRITE(momentsname(im),'("emp2emp | move,athme0",I2)') g
			weights(im)=0.0_rp
			im=im+1 
		end do 

		headloc(ihead)=im; headstr(ihead)='move rates by gender,at home or not at home and employment status (less than full time or fulltime) ';ihead=ihead+1
		do g=1,2
			CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%hhr>=0).AND.(dat%hhr<H_FULLTIME).AND.(dat%l/=dat%home).AND.(move>=0)),d1*move,moments,countmom,varmom)
			WRITE(momentsname(im),'("move | u,athme0 ",I2)') g
			weights(im)=0.0_rp
			im=im+1
			CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%hhr>=0).AND.(dat%hhr<H_FULLTIME).AND.(dat%l==dat%home).AND.(move>=0)),d1*move,moments,countmom,varmom)
			WRITE(momentsname(im),'("move | u,athme1 ",I2)') g
			weights(im)=0.0_rp
			im=im+1
			CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%hhr>=0).AND.(dat%hhr>=H_FULLTIME).AND.(dat%l/=dat%home).AND.(move>=0)),d1*move,moments,countmom,varmom)
			WRITE(momentsname(im),'("move | e,athme0 ",I2)') g
			weights(im)=0.0_rp
			im=im+1
			CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%hhr>=0).AND.(dat%hhr>=H_FULLTIME).AND.(dat%l==dat%home).AND.(move>=0)),d1*move,moments,countmom,varmom)
			WRITE(momentsname(im),'("move | u,athme1 ",I2)') g
			weights(im)=0.0_rp
			im=im+1
		end do 
			
		headloc(ihead)=im; headstr(ihead)='Employment by gender and spouse employment status ';ihead=ihead+1
		do g=1,2
			CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==1).AND.(dat%hhsp>=0).AND.(dat%hhsp<H_FULLTIME).AND.(dat%hhr>=0).AND.(iacat==1)),d1*abs(dat%hhr>=H_FULLTIME),moments,countmom,varmom)
			WRITE(momentsname(im),'("emp | sp<fulltime ",I4)') g
			weights(im)=0.0_rp
			im=im+1
			CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==1).AND.(dat%hhsp>=H_FULLTIME).AND.(dat%hhr>=0).AND.(iacat==1)),d1*abs(dat%hhr>=H_FULLTIME),moments,countmom,varmom)
			WRITE(momentsname(im),'("emp | sp->=fulltime ",I4)') g
			weights(im)=0.0_rp
			im=im+1
		end do 
		!headloc(ihead)=im; headstr(ihead)='Mar rates by education, ias 25-35';ihead=ihead+1
		CALL condmom(im,((dat(25:35,:)%co==co).AND.(dat(25:35,:)%rel>=0).AND.(dat(25:35,:)%edr==1)),abs(d1*(dat(25:35,:)%rel==1)),moments,countmom,varmom)
		WRITE(momentsname(im),'("frac married, hs 25:35")')
		weights(im)=0.0_rp	    		
		im=im+1
		CALL condmom(im,((dat(25:35,:)%co==co).AND.(dat(25:35,:)%rel>=0).AND.(dat(25:35,:)%edr==2)),abs(d1*(dat(25:35,:)%rel==1)),moments,countmom,varmom)
		WRITE(momentsname(im),'("frac married, col , 25:35")')
		weights(im)=0.0_rp
		im=im+1
		headloc(ihead)=im
		IF (co==1) THEN
			headstr(ihead)='br 1: Marriage Rates by Sex and ia';ihead=ihead+1
		ELSE
			headstr(ihead)='br 2: Marriage Rates by Sex and ia';ihead=ihead+1
		ENDIF
		do ia=25,35,2
			CALL condmom(im,((dat(ia,:)%co==co).AND.(dat(ia,:)%rel>=0).AND.((dat(ia,:)%edr==1))),abs(d1*(dat(ia,:)%rel>0)),moments,countmom,varmom)
			WRITE(momentsname(im),'("frac get together, hs  ",1I2)') ia
			weights(im)=0.0_rp 
			im=im+1
			CALL condmom(im,((dat(ia,:)%co==co).AND.(dat(ia,:)%rel>=0).AND.((dat(ia,:)%edr==2))),abs(d1*(dat(ia,:)%rel>0)),moments,countmom,varmom)
			WRITE(momentsname(im),'("frac get together, col  ",1I2)') ia
			weights(im)=0.0_rp 
			im=im+1
		end do  
		headloc(ihead)=im; headstr(ihead)='Move to rates by gender and relstat ';ihead=ihead+1
		do j=1,NL
			CALL condmom(im,((dat(MNA:MXAI-1,:)%co==co).AND.(dat(MNA:MXAI-1,:)%sexr==1).AND.(dat(MNA:MXAI-1,:)%rel==0).AND.(move(MNA:MXAI-1,:)==1).AND.(dat(MNA+1:MXAI,:)%l>=0)),d1*abs(dat(MNA+1:MXAI,:)%l==j),moments,countmom,varmom)
			WRITE(momentsname(im),'("male,single,prop of moves to each loc ",I2)') j
			weights(im)=0.0_rp
			im=im+1 
		end do 	
		do j=1,NL
			CALL condmom(im,((dat(MNA:MXAI-1,:)%co==co).AND.(dat(MNA:MXAI-1,:)%sexr==1).AND.(dat(MNA:MXAI-1,:)%rel==1).AND.(move(MNA:MXAI-1,:)==1).AND.(dat(MNA+1:MXAI,:)%l>=0)),d1*abs(dat(MNA+1:MXAI,:)%l==j),moments,countmom,varmom)
			WRITE(momentsname(im),'("male,   mar,prop of moves to each loc ",I2)') j
			weights(im)=0.0_rp
			im=im+1 
		end do 	
		do j=1,NL
			CALL condmom(im,((dat(MNA:MXAI-1,:)%co==co).AND.(dat(MNA:MXAI-1,:)%sexr==2).AND.(dat(MNA:MXAI-1,:)%rel==0).AND.(move(MNA:MXAI-1,:)==1).AND.(dat(MNA+1:MXAI,:)%l>=0)),d1*abs(dat(MNA+1:MXAI,:)%l==j),moments,countmom,varmom)
			WRITE(momentsname(im),'("fem, single,prop of moves to each loc ",I2)') j
			weights(im)=0.0_rp
			im=im+1 
		end do 	
		do j=1,NL
			CALL condmom(im,((dat(MNA:MXAI-1,:)%co==co).AND.(dat(MNA:MXAI-1,:)%sexr==2).AND.(dat(MNA:MXAI-1,:)%rel==1).AND.(move(MNA:MXAI-1,:)==1).AND.(dat(MNA+1:MXAI,:)%l>=0)),d1*abs(dat(MNA+1:MXAI,:)%l==j),moments,countmom,varmom)
			WRITE(momentsname(im),'("fem,    mar,prop of moves to each loc ",I2)') j
			weights(im)=0.0_rp
			im=im+1 
		end do 	
		headloc(ihead)=im; headstr(ihead)='move from rates extras ';ihead=ihead+1
		do j=1,NL
			CALL condmom(im,((dat(MNA:MXAI-1,:)%co==co).AND.(move(MNA:MXAI-1,:)==1).AND.(dat(MNA:MXAI-1,:)%l==j).AND.(dat(MNA:MXAI-1,:)%home>=0)),d1*abs(dat(MNA:MXAI-1,:)%home/=j),moments,countmom,varmom)
			WRITE(momentsname(im),'("prop of non-home moves from loc     ",I2)') j
			weights(im)=0.0_rp
			im=im+1 
			CALL condmom(im,((dat(MNA:MXAI-1,:)%co==co).AND.(move(MNA:MXAI-1,:)==1).AND.(dat(MNA:MXAI-1,:)%l==j).AND.(dat(MNA:MXAI-1,:)%home>=0)),d1*abs(dat(MNA:MXAI-1,:)%home==j),moments,countmom,varmom)
			WRITE(momentsname(im),'("prop of home moves from loc         ",I2)') j
			weights(im)=0.0_rp
			im=im+1 
		end do 	
		! wages for working women interacted with fraction of last 4 years out of labor force, by ed
		!headloc(ihead)=im; headstr(ihead)='wages for working women interacted with fraction of last 4 years out of labor force, by ed';ihead=ihead+1
		!DO ia=22,34,4
		!	CALL condmom(im,((dat(ia,:)%co==co).AND.(dat(ia,:)%sexr==2).AND.(dat(ia,:)%logwr>0).AND.(dat(ia,:)%edr==1).AND.(obs4h(ia,:))),d1*dat(ia,:)%logwr*frac4h0(ia,:),moments,countmom,varmom)
		!	WRITE(momentsname(im),'("logw-years out, fem hs,",1I2)') ia
		!	weights(im)=0.
		!	im=im+1
		!ENDDO
		!DO ia=26,34,4
		!	CALL condmom(im,((dat(ia,:)%co==co).AND.(dat(ia,:)%sexr==2).AND.(dat(ia,:)%logwr>0).AND.(dat(ia,:)%edr==2).AND.(obs4h(ia,:))),d1*dat(ia,:)%logwr*frac4h0(ia,:),moments,countmom,varmom)
		!	WRITE(momentsname(im),'("logw-years out, fem col,",1I2)') ia
		!	weights(im)=0.
		!	im=im+1
		!ENDDO
	ENDDO cohort
	END SUBROUTINE get_extramoments


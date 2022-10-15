            !ahu 122818 commenting out the loc moment
            !do i=1,nl
            !    call condmom(im,( cosexrel(MNA:MXA,:) .AND. dat(MNA:MXA,:)%hhr>=0 .AND. dat(MNA:MXA,:)%l>0 ),   d1*one( dat(MNA:MXA,:)%l==i ),mom,cnt,var)		
            !    write(name(im),'("loc",tr11,i2)') i			
            !    weights(im)=wmovebyrel ; if (onlysingles.and.j==1) weights(im)=0.0_dp
            !    im=im+1 
            !end do 
        
            !do i=1,nl
            !    call condmom(im,( cosexrel(MNA:MXA,:) .AND. dat(MNA:MXA,:)%hhr>=0 .AND. dat(MNA:MXA,:)%l==i ),   d1*one( dat(MNA:MXA,:)%hhr==1 ),mom,cnt,var)		
            !    write(name(im),'("e|loc",tr9,i4)') i			
            !    weights(im)=whour  ; if (onlysingles.and.j==1) weights(im)=0.0_dp
            !    im=im+1
            !end do 

!*****************************
do i=1,ndat
    !Just generating cohort and sex of someone because don't have initcond in this routine
    cohogen(i)=   maxval(dat(MNAD:MXA,i)%co)
    sexgen(i)=   maxval(dat(MNAD:MXA,i)%sexr)	                

    !Total number of moves overall
    !nummove(i)=sum(move(MNA:MXAI,i),  MASK=( move(MNA:MXAI,i)>0  .and. norelchg(MNA:MXAI,i)==1  )   )
    nummove(i)=sum(move(MNA:MXAD,i),  MASK=( move(MNA:MXAD,i)>=0  .and. norelchg(MNA:MXAD,i)==1  )   )

    !do ia=MNA,MXAI
    !    nummov(ia,i)=sum(move(MNA:ia,i),  MASK=( move(MNA:ia,i)>0  )   )
    !end do        
    !nummove_mar: Total number of moves that occurred when married, up to ia
    !nummove_sin: Total number of moves that occurred when single, up to ia
    !nummove_ma(i)=sum(move(MNA:MXAI,i),  MASK=( move(MNA:MXAI,i)>0.and.dat(MNA:MXAI,i)%rel==1)  )
    !nummove_si(i)=sum(move(MNA:MXAI,i),  MASK=( move(MNA:MXAI,i)>0.and.dat(MNA:MXAI,i)%rel==0)  )
    !**************************************
    !ahu jan19 010119: commenting out the below and not doing the nummove by rel moments 
    !this is because for example in data nummove-0 is 0.84 while both nummove_mar and _sin are 0.92 something. 
    !I think the norelchg requirement is the reason. With that requirement it is not clear what these moments mean and 
    !they might be hindering the matching of the move by age moments. For example, there are times when nummove=0 is understated by a whole lot in sim
    !whereas movey age is overstated 
    !do ia=MNA,MXAI
    !    nummove_mar(ia,i)=sum(move(MNA:ia,i),  MASK=( move(MNA:ia,i)>0.and.dat(MNA:ia,i)%rel==1 .and. norelchg(MNA:ia,i)==1 )  )
    !    nummove_sin(ia,i)=sum(move(MNA:ia,i),  MASK=( move(MNA:ia,i)>0.and.dat(MNA:ia,i)%rel==0 .and. norelchg(MNA:ia,i)==1 )  )
    !end do
    !if (iter==1) then
    !    if ( nummove_mar(MXAI,i) /= nummove_mar(MXAI-1,i) .OR. minval(nummove_mar(:,i))<0 ) then
    !        print*, 'something wrong with nummove_mar'
    !        stop
    !    end if 
    !    if ( nummove_sin(MXAI,i) /= nummove_sin(MXAI-1,i) .OR. minval(nummove_sin(:,i))<0 ) then
    !        print*, 'something wrong with nummove_sin'
    !        stop
    !    end if 
    !end if 
    !**************************************

end do 

!*********************************
if (extramoments) then 
    mean4h=-99.0_dp 
    deltawage4=-99.0_dp
    deltawage=-99.0_dp
    obs4h=.FALSE.
    iacat=-99
    empr=-99
    moverank=-99
    iacat=0
    iacat(18:40,:)=1
    !NOBS=SIZE(dat,2) ! number of observations.
    !deltawage4=-1d0
    ! Averia hours worked past four years, averia of square of hours worked last four years, fraction of
    !   last four years out of labor force. Interact these with current wages and/or wage growth
    empr(MNAD:MXA,:)=one( dat(MNAD:MXA,:)%hhr>=H_FULLTIME )
    DO ia=MNA+4,MXA
        ! are there observations of hours in the past 4 years 
        obs4h(ia,:)=(SUM(one(1.*dat(ia-4:ia-1,:)%hhr>=0),1)>0)
        WHERE (obs4h(ia,:))
            mean4h(ia,:)=SUM(empr(ia-4:ia-1,:),1,empr(ia-4:ia-1,:)>=0)/SUM(one(1.*empr(ia-4:ia-1,:)>=0),1)
            !mean4hsq(ia,:)=SUM((dat(ia-4:ia-1,:)%hhr)**2,1,dat(ia-4:ia-1,:)%hhr>=0)/SUM(one(1.*dat(ia-4:ia-1,:)%hhr>=0),1)
            !frac4h0(ia,:)=SUM(one(1.*dat(ia-4:ia-1,:)%hhr==0),1)/SUM(one(1.*dat(ia-4:ia-1,:)%hhr>=0),1)
            deltawage4(ia,:)=dat(ia,:)%logwr-dat(ia-4,:)%logwr
        ENDWHERE
    ENDDO  
    
    !WHERE (move(MNA:MXAD,:)>=0.and.dat(MNA:MXAD,:)%rel>=0)
    !    move_single(MNA:MXAD,:)=one(   move(MNA:MXAD,:)==1  .and.  dat(MNA:MXAD,:)%rel==0   )
    !    move_mar(MNA:MXAD,:)=one(   move(MNA:MXAD,:)==1   .and.   dat(MNA:MXAD,:)%rel==1  )
    !ENDWHERE
    !do ia=MNA,MXAD
    !    movesum_single(ia,:)=sum(move_single(MNA:ia,:),1,move_single>=0)
    !    movesum_mar(ia,:)=sum(move_mar(MNA:ia,:),1,move_mar>=0)
    !    where (movesum_single(ia,:)>=1)
    !        movesum_single(ia,:)=1
    !    elsewhere (movesum_single(ia,:)==0)
    !        movesum_single(ia,:)=0
    !    endwhere 
    !    where (movesum_mar(ia,:)>=1)
    !        movesum_mar(ia,:)=1
    !    elsewhere (movesum_mar(ia,:)==0)
    !        movesum_mar(ia,:)=0
    !    endwhere 
    !    end do
    
    !deltawage=-1.	
    WHERE ((dat(MNA:MXAD,:)%logwr>0) .AND. (dat(MNA+1:MXA,:)%logwr>0))
        deltawage(MNA:MXAD,:)=dat(MNA+1:MXA,:)%logwr-dat(MNA:MXAD,:)%logwr
    ENDWHERE
    
   ! headloc(ihead)=im ; headstr(ihead)='everyone misc' ; ihead=ihead+1
    headloc(ihead)=im; headstr(ihead)='employment by gender/rel (all ages) ';ihead=ihead+1
    do g=1,2
        do j=0,maxrelo
            CALL condmom(im,((dat(MNA:MXA,:)%co==co).AND.(dat(MNA:MXA,:)%sexr==g).AND.(dat(MNA:MXA,:)%rel==j).AND.(dat(MNA:MXA,:)%hhr>=0).AND.(iacat(MNA:MXA,:)==1)),d1*one(dat(MNA:MXA,:)%hhr>=H_FULLTIME),mom,cnt,var)
            WRITE(name(im),'("emp | gender/rel ",2I4)') g,j
            if (j==0) weights(im)=0.0_dp !whour0 
            if (j==1) weights(im)=0.0_dp !whour1 
            im=im+1
        end do 
    end do 
    headloc(ihead)=im; headstr(ihead)='employment by gender/kid (married and all age)';ihead=ihead+1
    do g=1,2
        do j=0,1
            CALL condmom(im,((dat(MNA:MXA,:)%co==co).AND.(dat(MNA:MXA,:)%sexr==g).AND.(dat(MNA:MXA,:)%rel==1).AND.(dat(MNA:MXA,:)%kidr==j).AND.(dat(MNA:MXA,:)%hhr>=0).AND.(iacat(MNA:MXA,:)==1)),d1*one(dat(MNA:MXA,:)%hhr>=H_FULLTIME),mom,cnt,var)
            WRITE(name(im),'("emp | gender/kid ",2I4)') g,j
            weights(im)=0.0_dp !whour1 
            im=im+1 
        end do 
    end do 
    
    
    !headloc(ihead)=im ; headstr(ihead)='everyone misc' ; ihead=ihead+1
    headloc(ihead)=im; headstr(ihead)='wage-growth for stay and move by gender/rel (all ages)';ihead=ihead+1
    do g=1,2
        do j=0,maxrelo 
            CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%sexr==g).AND.(dat(MNA:MXAD,:)%rel==j).AND.(move(MNA:MXAD,:)==0).AND.(dat(MNA:MXAD,:)%logwr>=0).AND.(deltawage(MNA:MXAD,:)>-100.0_dp).AND.(iacat(MNA:MXAD,:)==1)),d1*deltawage(MNA:MXAD,:),mom,cnt,var)
            WRITE(name(im),'("wage-growth/stay/gender/rel ",2I4)') g,j
            if (j==0) weights(im)=0.0_dp !wwage0 
            if (j==1) weights(im)=0.0_dp !wwage1 
            im=im+1
            CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%sexr==g).AND.(dat(MNA:MXAD,:)%rel==j).AND.(move(MNA:MXAD,:)==1).AND.(dat(MNA:MXAD,:)%logwr>=0).AND.(deltawage(MNA:MXAD,:)>-100.0_dp).AND.(iacat(MNA:MXAD,:)==1)),d1*deltawage(MNA:MXAD,:),mom,cnt,var)
            WRITE(name(im),'("wage-growth/move/gender/rel ",2I4)') g,j
            if (j==0) weights(im)=0.0_dp !wwage0 
            if (j==1) weights(im)=0.0_dp !wwage1 
            im=im+1
        end do 
    end do 
        

    
    !headloc(ihead)=im ; headstr(ihead)='everyone misc' ; ihead=ihead+1
    headloc(ihead)=im; headstr(ihead)='education by location';ihead=ihead+1
    do j=1,NL
        CALL condmom(im,((dat(MNA:MXA,:)%co==co).AND.(dat(MNA:MXA,:)%l==j).AND.(dat(MNA:MXA,:)%edr>=0)),d1*one((dat(MNA:MXA,:)%edr==2)),mom,cnt,var)
        WRITE(name(im),'("prop-col ",I4)') j
        weights(im)=0.0_dp !wwage0
        im=im+1
    end do 

    
    !headloc(ihead)=im ; headstr(ihead)='everyone misc' ; ihead=ihead+1
    headloc(ihead)=im; headstr(ihead)='Prop by loc and Prop of moves by loc (all age) ';ihead=ihead+1
    do j=1,NL
        CALL condmom(im,((dat(MNA:MXA,:)%co==co)),d1*one(dat(MNA:MXA,:)%l==j),mom,cnt,var)
        WRITE(name(im),'("prop-loc                  ",I4)') j
        weights(im)=0.0_dp !wmove0
        im=im+1 
    end do 

    headloc(ihead)=im
    IF (co==1) THEN
        headstr(ihead)='marriage rates by sex and ia';ihead=ihead+1
    ELSE
        headstr(ihead)='marriage rates by sex and ia';ihead=ihead+1
    ENDIF
    do ia=25,35,5
        CALL condmom(im,(     (dat(ia,:)%co==co).AND.(dat(ia,:)%rel>=0).AND.(dat(ia,:)%edr==1)   ),d1*one(dat(ia,:)%rel>0),mom,cnt,var)
        WRITE(name(im),'("frac get together, hs  ",1I2)') ia
        weights(im)=0.0_dp !0.0_rp 
        im=im+1
        CALL condmom(im,((dat(ia,:)%co==co).AND.(dat(ia,:)%rel>=0).AND.(dat(ia,:)%edr==2)),d1*one(dat(ia,:)%rel>0),mom,cnt,var)
        WRITE(name(im),'("frac get together, col  ",1I2)') ia
        weights(im)=0.0_dp !0.0_rp 
        im=im+1
    end do  
    headloc(ihead)=im; headstr(ihead)='move from rates extras ';ihead=ihead+1
    do j=1,NL
        CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(move(MNA:MXAD,:)==1).AND.(dat(MNA:MXAD,:)%l==j).AND.(dat(MNA:MXAD,:)%hme>=0)),d1*one(dat(MNA:MXAD,:)%hme/=j),mom,cnt,var)
        WRITE(name(im),'("prop of non-home moves from loc     ",I2)') j
        weights(im)=0.0_dp !0.0_rp
        im=im+1 
        CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(move(MNA:MXAD,:)==1).AND.(dat(MNA:MXAD,:)%l==j).AND.(dat(MNA:MXAD,:)%hme>=0)),d1*one(dat(MNA:MXAD,:)%hme==j),mom,cnt,var)
        WRITE(name(im),'("prop of home moves from loc         ",I2)') j
        weights(im)=0.0_dp !0.0_rp
        im=im+1 
    end do  
    
    
end if !extramoments

!*********************
        !do ia=mna,mxai
        !    call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel>=0  .AND. dat(ia,:)%edr==1 ), d1*one(dat(ia,:)%rel==1),mom,cnt,var)
        !    write(name(im),'("mar by ia,ned",tr1,i4)') ia
        !    weights(im)=0.0_dp
        !    im=im+1
        !    call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel>=0  .AND. dat(ia,:)%edr==2 ), d1*one(dat(ia,:)%rel==1),mom,cnt,var)
        !    write(name(im),'("mar by ia, ed",tr1,i4)') ia
        !    weights(im)=0.0_dp
        !    im=im+1
        !end do            
                         
        !do ia=mna,mxai
        !    call condmom(im,( coho(ia,:) .AND. dat(ia,:)%rel>=0 ), d1*one(dat(ia,:)%rel==1),mom,cnt,var)
        !    write(name(im),'("mar by ia",tr5,i4)') ia
        !    weights(im)=wrel
        !    im=im+1
        !end do

!*************************
        !etr=-1.
            !ahu 062512 given the new definition of hours (was discrete and now it's continuous, changing the definition of etr) 
            !  <FULLTIME TO >=FULLTIME
            WHERE ( (dat(MNA:MXAD,:)%hhr<H_FULLTIME).AND.(dat(MNA:MXAD,:)%hhr>=0) .AND. (dat(MNA+1:MXA,:)%hhr>=0) ) !ahu 062512 check this -1 
            etr(1,MNA:MXAD,:)=one(dat(MNA+1:MXA,:)%hhr>=H_FULLTIME)
        ENDWHERE
        WHERE ( (dat(MNA:MXAD,:)%hhr>=H_FULLTIME).AND.(dat(MNA:MXAD,:)%hhr>=0) .AND. (dat(MNA+1:MXA,:)%hhr>=0) ) !ahu 062512 check this -1 
            etr(2,MNA:MXAD,:)=one(dat(MNA+1:MXA,:)%hhr>=H_FULLTIME)
        ENDWHERE
        !  <PARTTIME TO >=PARTTIME
        WHERE ( (dat(MNA:MXAD,:)%hhr<H_PARTTIME).AND.(dat(MNA:MXAD,:)%hhr>=0) .AND. (dat(MNA+1:MXA,:)%hhr>=0) ) !ahu 062512 check this -1 
            etr(3,MNA:MXAD,:)=one(dat(MNA+1:MXA,:)%hhr>=H_PARTTIME)
        ENDWHERE
        WHERE ( (dat(MNA:MXAD,:)%hhr>=H_PARTTIME).AND.(dat(MNA:MXAD,:)%hhr>=0) .AND. (dat(MNA+1:MXA,:)%hhr>=0) ) !ahu 062512 check this -1 
            etr(4,MNA:MXAD,:)=one(dat(MNA+1:MXA,:)%hhr>=H_PARTTIME)
        ENDWHERE
        !  <FULLTIME TO <FULLTIME
        WHERE ( (dat(MNA:MXAD,:)%hhr<H_FULLTIME).AND.(dat(MNA:MXAD,:)%hhr>=0) .AND. (dat(MNA+1:MXA,:)%hhr>=0) ) !ahu 062512 check this -1 
            etr(5,MNA:MXAD,:)=one(dat(MNA+1:MXA,:)%hhr<H_FULLTIME)
        ENDWHERE

!*************************************
        !headloc(ihead)=im ; headstr(ihead)='everyone misc' ; ihead=ihead+1
        headloc(ihead)=im; headstr(ihead)='move rates (all age, both sex';ihead=ihead+1
        CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%rel==0).and.(norelchg(MNA:MXAD,:)==1).AND.(move(MNA:MXAD,:)>=0)),d1*move(MNA:MXAD,:),mom,cnt,var)
        WRITE(name(im),'("move/single")') 
        weights(im)=0.0_dp !wmove0              
        im=im+1 
        CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%rel==1).and.(norelchg(MNA:MXAD,:)==1).AND.(move(MNA:MXAD,:)>=0)),d1*move(MNA:MXAD,:),mom,cnt,var)
        WRITE(name(im),'("move/married")') 
        weights(im)=0.0_dp !wmove1              
        im=im+1
!*************************************



do i=1,5,4
    call condmom(im,( cosex(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0  .AND. dur(MNA:MXAD,:)==i .AND. dat(MNA:MXAD,:)%edr==1 ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
    write(name(im),'("e|u by dur ned",tr5,i2)') i
    weights(im)=whour
    im=im+1 
end do 
do i=1,5,4
    call condmom(im,( cosex(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0  .AND. dur(MNA:MXAD,:)==i .AND. dat(MNA:MXAD,:)%edr==2 ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
    write(name(im),'("e|u by dur  ed",tr5,i2)') i
    weights(im)=whour
    im=im+1 
end do 
do i=1,5,4
    call condmom(im,( cosex(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%logwr>=0  .AND. dur(MNA:MXAD,:)==i  .AND. dat(MNA:MXAD,:)%edr==1),   d1*dat(MNA+1:MXA,:)%logwr ,mom,cnt,var)		
    write(name(im),'("w|u by dur ned",tr5,i2)') i
    weights(im)=wwage
    im=im+1 
end do 
do i=1,5,4
    call condmom(im,( cosex(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%logwr>=0  .AND. dur(MNA:MXAD,:)==i  .AND. dat(MNA:MXAD,:)%edr==2),   d1*dat(MNA+1:MXA,:)%logwr ,mom,cnt,var)		
    write(name(im),'("w|u by dur  ed",tr5,i2)') i
    weights(im)=wwage
    im=im+1 
end do 




!************************************
        !do i=1,ntypp
        !    call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%typ==i .and. dat(MNA:MXA,:)%edr==1 ), d1*one(dat(MNA:MXA,:)%rel==1),mom,cnt,var)
        !    write(name(im),'("mar|ned by typ ",i4)') i
        !    weights(im)=0.0_dp
        !    im=im+1
        !end do 

        !  do i=1,ntypp
        !    call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%typ==i .and. dat(MNA:MXA,:)%edr==2 ), d1*one(dat(MNA:MXA,:)%rel==1),mom,cnt,var)
        !    write(name(im),'("mar| ed by typ ",i4)') i
        !    weights(im)=0.0_dp
        !    im=im+1
        !end do 

        
        !do i=1,ntypp
        !    call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%rel==1  .and. dat(MNA:MXA,:)%edr==1 ), d1*one(dat(MNA:MXA,:)%typ==i),mom,cnt,var)
        !    write(name(im),'("typ | mar,noed ",tr1,i4)') i
        !    weights(im)=0.0_dp
        !    im=im+1
        !end do 

        !do i=1,ntypp
        !    call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 .AND. dat(MNA:MXA,:)%rel==1  .and. dat(MNA:MXA,:)%edr==2 ), d1*one(dat(MNA:MXA,:)%typ==i),mom,cnt,var)
        !    write(name(im),'("typ | mar,  ed ",tr1,i4)') i
        !    weights(im)=0.0_dp
        !    im=im+1
        !end do     
        

            !do ia=mnaD,18,1
            !    do i=1,ntypp
        !		    CALL condmom(im,( cosexrel(ia,:) .AND. move(ia,:)>=0 .AND. dat(ia,:)%typ==i ),d1*move(ia,:),mom,cnt,var)
        !		    WRITE(name(im),'("mv by a,tp",I8,I4)') ia,i
        !		    weights(im)=0.0_dp
        !		    im=im+1
            !       end do
            !   end do            
    !		do ia=19,25,2
    !            do i=1,ntypp
    !				    CALL condmom(im,( cosexrel(ia,:) .AND. move(ia,:)>=0 .AND. dat(ia,:)%typ==i ),d1*move(ia,:),mom,cnt,var)
    !			    WRITE(name(im),'("mv by a,tp",I8,I4)') ia,i
    !			    weights(im)=0.0_dp
    !			    im=im+1
        !           end do
    !        end do            




                
        
                
				!do ia=agestart(NOCOLLEGE),45,6 !mxai-1
                !    if ( j==0) then 
    		!			CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==1 .AND. dat(ia,:)%wr>=0 .AND. nummove_sin(ia,:)==0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
	        !        else if (j==1) then
    	!				CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==1 .AND. dat(ia,:)%wr>=0 .AND. nummove_mar(ia,:)==0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
         !           end if                    
          !          WRITE(name(im),'("inedmvrel=0 ",I2)') ia ; mominfo(0,im)=11 ; mominfo(1,im)=g ; mominfo(2,im)=j ; mominfo(3,im)=1  ; mominfo(4,im)=0 ; mominfo(5,im)=ia
			!		weights(im)=wwagebymove ; if (onlysingles.and.j==1) weights(im)=0.0_dp
			!		im=im+1
		!		end do            
            
		!		do ia=agestart(NOCOLLEGE),45 !mxai-1
         !           if ( j==0) then 
    		!			CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==1 .AND. dat(ia,:)%wr>=0 .AND. nummove_sin(ia,:)>0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
            !       else if (j==1) then
    	!				CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==1 .AND. dat(ia,:)%wr>=0 .AND. nummove_mar(ia,:)>0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
         !           end if                    
          !          WRITE(name(im),'("inedmvrel>0 ",I2)') ia ; mominfo(0,im)=11  ; mominfo(1,im)=g ; mominfo(2,im)=j ; mominfo(3,im)=1  ; mominfo(4,im)=1 ; mominfo(5,im)=ia
			!		weights(im)=wwagebymove ; if (onlysingles.and.j==1) weights(im)=0.0_dp
			!		im=im+1
		!		end do            

		!		do ia=agestart(COLLEGE),45 !mxai-1
         !           if ( j==0) then 
    		!			CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==2 .AND. dat(ia,:)%wr>=0 .AND. nummove_sin(ia,:)==0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
	         !       else if (j==1) then
    			!		CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==2 .AND. dat(ia,:)%wr>=0 .AND. nummove_mar(ia,:)==0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
                !    end if                    
                !    WRITE(name(im),'("i edmvrel=0 ",I2)') ia ; mominfo(0,im)=11  ; mominfo(1,im)=g ; mominfo(2,im)=j ; mominfo(3,im)=2 ; mominfo(4,im)=0 ; mominfo(5,im)=ia
			!		weights(im)=wwagebymove ; if (onlysingles.and.j==1) weights(im)=0.0_dp
		!			im=im+1
		!		end do            
				
		!		do ia=agestart(COLLEGE),45 !mxai-1
         !           if ( j==0) then 
    		!			CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==2 .AND. dat(ia,:)%wr>=0 .AND. nummove_sin(ia,:)>0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
	         !       else if (j==1) then
    			!		CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==2 .AND. dat(ia,:)%wr>=0 .AND. nummove_mar(ia,:)>0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
                !    end if                    
                !    WRITE(name(im),'("i edmvrel>0 ",I2)') ia ; mominfo(0,im)=11  ; mominfo(1,im)=g ; mominfo(2,im)=j ; mominfo(3,im)=2  ; mominfo(4,im)=1 ; mominfo(5,im)=ia
			!		weights(im)=wwagebymove ; if (onlysingles.and.j==1) weights(im)=0.0_dp
		!			im=im+1
		!		end do            
                
				!do ia=agestart(NOCOLLEGE),45 !mxai-1
                !    CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==1 .AND. dat(ia,:)%wr>=0 .AND. nummove(:)==0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
                !    WRITE(name(im),'("inedmvrel=0 ",I2)') ia ; mominfo(0,im)=8 ; mominfo(1,im)=g ; mominfo(2,im)=j ; mominfo(3,im)=1  ; mominfo(4,im)=0 ; mominfo(5,im)=ia
				!	weights(im)=wwage
				!	im=im+1
				!end do            
            
				!do ia=agestart(NOCOLLEGE),45 !mxai-1
                !    CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==1 .AND. dat(ia,:)%wr>=0 .AND. nummove(:)>0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
                !    WRITE(name(im),'("inedmvrel>0 ",I2)') ia ; mominfo(0,im)=8 ; mominfo(1,im)=g ; mominfo(2,im)=j ; mominfo(3,im)=1  ; mominfo(4,im)=1 ; mominfo(5,im)=ia
				!	weights(im)=wwage
				!	im=im+1
				!end do            

				!do ia=agestart(COLLEGE),45 !mxai-1
                !    CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==2 .AND. dat(ia,:)%wr>=0 .AND. nummove(:)==0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
                !    WRITE(name(im),'("i edmvrel=0 ",I2)') ia ; mominfo(0,im)=8 ; mominfo(1,im)=g ; mominfo(2,im)=j ; mominfo(3,im)=2 ; mominfo(4,im)=0 ; mominfo(5,im)=ia
				!	weights(im)=wwage
				!	im=im+1
				!end do            
				
				!do ia=agestart(COLLEGE),45 !mxai-1
                !    CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==2 .AND. dat(ia,:)%wr>=0 .AND. nummove(:)>0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
                !    WRITE(name(im),'("i edmvrel>0 ",I2)') ia ; mominfo(0,im)=8 ; mominfo(1,im)=g ; mominfo(2,im)=j ; mominfo(3,im)=2  ; mominfo(4,im)=1 ; mominfo(5,im)=ia
				!	weights(im)=wwage
				!	im=im+1
				!end do            
                

                
                
                
                
				!do ia=agestart(NOCOLLEGE),45 !mxai-1
                !    CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==1 .AND. dat(ia,:)%wr>=0 .AND. nummov(ia,:)==0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
                !    WRITE(name(im),'("inedmvrel=0 ",I2)') ia ; mominfo(0,im)=9 ; mominfo(1,im)=g ; mominfo(2,im)=j ; mominfo(3,im)=1  ; mominfo(4,im)=0 ; mominfo(5,im)=ia
				!	weights(im)=wwage
				!	im=im+1
				!end do            
            
				!do ia=agestart(NOCOLLEGE),45 !mxai-1
                !    CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==1 .AND. dat(ia,:)%wr>=0 .AND. nummov(ia,:)>0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
                !    WRITE(name(im),'("inedmvrel>0 ",I2)') ia ; mominfo(0,im)=9 ; mominfo(1,im)=g ; mominfo(2,im)=j ; mominfo(3,im)=1  ; mominfo(4,im)=1 ; mominfo(5,im)=ia
				!	weights(im)=wwage
				!	im=im+1
				!end do            

				!do ia=agestart(COLLEGE),45 !mxai-1
                !    CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==2 .AND. dat(ia,:)%wr>=0 .AND. nummov(ia,:)==0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
                !    WRITE(name(im),'("i edmvrel=0 ",I2)') ia ; mominfo(0,im)=9 ; mominfo(1,im)=g ; mominfo(2,im)=j ; mominfo(3,im)=2 ; mominfo(4,im)=0 ; mominfo(5,im)=ia
				!	weights(im)=wwage
				!	im=im+1
				!end do            
				
				!do ia=agestart(COLLEGE),45 !mxai-1
                !    CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==2 .AND. dat(ia,:)%wr>=0 .AND. nummov(ia,:)>0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
                !    WRITE(name(im),'("i edmvrel>0 ",I2)') ia ; mominfo(0,im)=9 ; mominfo(1,im)=g ; mominfo(2,im)=j ; mominfo(3,im)=2  ; mominfo(4,im)=1 ; mominfo(5,im)=ia
				!	weights(im)=wwage
				!	im=im+1
				!end do            
                
                

                            
                
                
                
				!do ia=agestart(NOCOLLEGE),45 !mxai-1
                !    if ( j==0) then 
    			!		CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==1 .AND. dat(ia,:)%wr>=0 .AND. nummove_si(:)==0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
	            !    else if (j==1) then
    			!		CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==1 .AND. dat(ia,:)%wr>=0 .AND. nummove_ma(:)==0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
                !    end if                    
                !    WRITE(name(im),'("inedmvrel=0 ",I2)') ia ; mominfo(0,im)=10 ; mominfo(1,im)=g ; mominfo(2,im)=j ; mominfo(3,im)=1  ; mominfo(4,im)=0 ; mominfo(5,im)=ia
				!	weights(im)=wwage
				!	im=im+1
				!end do            
            
				!do ia=agestart(NOCOLLEGE),45 !mxai-1
                !    if ( j==0) then 
    			!		CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==1 .AND. dat(ia,:)%wr>=0 .AND. nummove_si(:)>0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
                !   else if (j==1) then
    			!		CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==1 .AND. dat(ia,:)%wr>=0 .AND. nummove_ma(:)>0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
                !    end if                    
                !    WRITE(name(im),'("inedmvrel>0 ",I2)') ia ; mominfo(0,im)=10  ; mominfo(1,im)=g ; mominfo(2,im)=j ; mominfo(3,im)=1  ; mominfo(4,im)=1 ; mominfo(5,im)=ia
				!	weights(im)=wwage
				!	im=im+1
				!end do            

				!do ia=agestart(COLLEGE),45 !mxai-1
                !    if ( j==0) then 
    			!		CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==2 .AND. dat(ia,:)%wr>=0 .AND. nummove_si(:)==0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
	            !    else if (j==1) then
    			!		CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==2 .AND. dat(ia,:)%wr>=0 .AND. nummove_ma(:)==0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
                !    end if                    
                !    WRITE(name(im),'("i edmvrel=0 ",I2)') ia ; mominfo(0,im)=10  ; mominfo(1,im)=g ; mominfo(2,im)=j ; mominfo(3,im)=2 ; mominfo(4,im)=0 ; mominfo(5,im)=ia
				!	weights(im)=wwage
				!	im=im+1
				!end do            
				
				!do ia=agestart(COLLEGE),45 !mxai-1
                !    if ( j==0) then 
    			!		CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==2 .AND. dat(ia,:)%wr>=0 .AND. nummove_si(:)>0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
	            !    else if (j==1) then
    			!		CALL condmom(im,(   cosexrel(ia,:) .AND.  dat(ia,:)%hhr==1 .AND. dat(ia,:)%edr==2 .AND. dat(ia,:)%wr>=0 .AND. nummove_ma(:)>0 ) ,d1*dat(ia,:)%wr,mom,cnt,var)
                !    end if                    
                !    WRITE(name(im),'("i edmvrel>0 ",I2)') ia ; mominfo(0,im)=10  ; mominfo(1,im)=g ; mominfo(2,im)=j ; mominfo(3,im)=2  ; mominfo(4,im)=1 ; mominfo(5,im)=ia
				!	weights(im)=wwage
				!	im=im+1
				!end do                          
 

        !headloc(ihead)=im; headstr(ihead)='what proportion of people does each location have?';ihead=ihead+1
        !headloc(ihead)=im; headstr(ihead)='proportion of college graduates in each location (all ages)';ihead=ihead+1
        !headloc(ihead)=im; headstr(ihead)='among all moves, what is the % that is from location j';ihead=ihead+1
        !headloc(ihead)=im; headstr(ihead)='among all moves, what is the % that is	to location j';ihead=ihead+1
        !headloc(ihead)=im; headstr(ihead)='among all moves to loc j, what is the % that is and is not a back-to-home move?';ihead=ihead+1
        !headloc(ihead)=im; headstr(ihead)='after the first move, what proportion is return to home? (all ages)';ihead=ihead+1		
        !headloc(ihead)=im; headstr(ihead)='among all people in loc j, what is the % who moved?';ihead=ihead+1		

        !if (iwritegen==1) then ; print*, 'here is im-1 in get_mom',im-1 ; end if 
        !if ( (.not.chkstep).and.(.not.optimize) ) print*, 'here is im-1 in get_mom',im-1 
	
    
        !**************************************

  do j=1,NL
    CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA+1:MXA,:)%l==j).AND.(homemove(MNA:MXAD,:)>=0).AND.(norelchg(MNA:MXAD,:)==1) ),d1*one(homemove(MNA:MXAD,:)==0),mom,cnt,var)
    WRITE(name(im),'("%nonhme-mvs-to ",I4)') j
    weights(im)=0.0_dp !wmove0
    im=im+1 
    CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA+1:MXA,:)%l==j).AND.(homemove(MNA:MXAD,:)>=0)).AND.(norelchg(MNA:MXAD,:)==1) ,d1*one(homemove(MNA:MXAD,:)==1),mom,cnt,var)
    WRITE(name(im),'("%hme-mvs-to     ",I4)') j
    weights(im)=0.0_dp !wmove0
    im=im+1 
end do  
        !do ia=MNA,19 
        !    call condmom(im,( dat(ia,:)%rel==0 .AND. dat(ia+1,:)%rel>=0   .AND. dat(ia,:)%sexr==1 ), d1*one(dat(ia+1,:)%rel==1),mom,cnt,var)
        !    write(name(im),'("getmarbyia,m  ",tr1,i4)') ia
        !    weights(im)=wrel
        !    im=im+1
        !    call condmom(im,( dat(ia,:)%rel==0 .AND. dat(ia+1,:)%rel>=0   .AND. dat(ia,:)%sexr==2 ), d1*one(dat(ia+1,:)%rel==1),mom,cnt,var)
        !    write(name(im),'("getmarbyia,f  ",tr1,i4)') ia
        !    weights(im)=wrel
        !    im=im+1
        !end do            

        !call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0  .AND. dat(MNA:MXA,:)%edr==1 ), d1*one(dat(MNA:MXA,:)%rel==1),mom,cnt,var)
        !write(name(im),'("mar ned",tr7)')	
        !weights(im)=0.0_dp
        !im=im+1

        !call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0  .AND. dat(MNA:MXA,:)%edr==2 ), d1*one(dat(MNA:MXA,:)%rel==1),mom,cnt,var)
        !write(name(im),'("mar  ed",tr7)')	
        !weights(im)=0.0_dp
        !im=im+1
    
    
        !call condmom(im,( coho(MNA:MXA,:) .AND. dat(MNA:MXA,:)%rel>=0 ), d1*one(dat(MNA:MXA,:)%rel==0),mom,cnt,var)  
        !write(name(im),'("sin ",tr10)')              
        !weights(im)=wrel
        !im=im+1

            

        !do ii=1,2
        !    do jj=1,2
        !    
        !        call condmom(im,(  corel(MNA:MXAD,:) .AND. edh(MNA:MXAD,:)==ii  .AND. edw(MNA:MXAD,:)==jj  .AND. emph(MNA:MXAD,:)==1  .AND. empw(MNA:MXAD,:)==1  .AND. logwh(MNA:MXAD,:)>=0 .AND. logww(MNA:MXAD,:)>=0 ),   d1*logwh(MNA:MXAD,:)  ,mom,cnt,var)	
        !        write(name(im),'("corr1 ",tr4,2i4)')  ii,jj
        !        weights(im)=0.0_dp
        !        calcorr(im)=1
        !        im=im+1
        !        call condmom(im,(  corel(MNA:MXAD,:) .AND. edh(MNA:MXAD,:)==ii  .AND. edw(MNA:MXAD,:)==jj  .AND. emph(MNA:MXAD,:)==1  .AND. empw(MNA:MXAD,:)==1  .AND. logwh(MNA:MXAD,:)>=0 .AND. logww(MNA:MXAD,:)>=0 ),   d1*logww(MNA:MXAD,:)  ,mom,cnt,var)	
        !        write(name(im),'("corr2 ",tr4,2i4)')  ii,jj
        !        weights(im)=0.0_dp
        !        calcorr(im)=5
        !        im=im+1
        !        call condmom(im,(  corel(MNA:MXAD,:) .AND. edh(MNA:MXAD,:)==ii  .AND. edw(MNA:MXAD,:)==jj  .AND. emph(MNA:MXAD,:)==1  .AND. empw(MNA:MXAD,:)==1  .AND. logwh(MNA:MXAD,:)>=0 .AND. logww(MNA:MXAD,:)>=0 ),   d1*logwh(MNA:MXAD,:)*logww(MNA:MXAD,:)  ,mom,cnt,var)	
        !        write(name(im),'("corr3 ",tr4,2i4)')  ii,jj
        !        weights(im)=wwage
        !        calcorr(im)=5
        !        im=im+1
        !        call condmom(im,(  corel(MNA:MXAD,:) .AND. edh(MNA:MXAD,:)==ii  .AND. edw(MNA:MXAD,:)==jj  .AND. emph(MNA:MXAD,:)==1  .AND. empw(MNA:MXAD,:)==1  .AND. logwh(MNA:MXAD,:)>=0 .AND. logww(MNA:MXAD,:)>=0 ),   d1*(logwh(MNA:MXAD,:)-mom(im-3))*(logww(MNA:MXAD,:)-mom(im-2))  ,mom,cnt,var)	
        !        write(name(im),'("wrngcorr ",tr4,2i4)')  ii,jj
        !        weights(im)=0.0_dp
        !        im=im+1

            !   end do 
        !end do 
        
        !headloc(ihead)=im
        !headstr(ihead)='MARRIED MOVE'
        !ihead=ihead+1
        
        !call condmom(im,( corel(MNA:MXAD,:) .AND. emph(MNA:MXAD,:)==1  .AND. empw(MNA:MXAD,:)==1 .AND. move(MNA:MXAD,:)>=0 ),   d1* move(MNA:MXAD,:) ,mom,cnt,var)	
        !write(name(im),'("move | ee ",tr5)')  
        !weights(im)=wmovemar   !ahu 122818 changed from wmove to wmovemar2 
        !im=im+1
        !!call condmom(im,( corel(mna:mxai-1,:) .AND. dat(MNA:MXAI-1,:)%hhr==1 .AND. dat(MNA:MXAI-1,:)%hhsp==1  .AND. dat(MNA:MXAI-1,:)%sexr==1 .AND. move(mna:mxai-1,:)>=0 ),   d1* move(mna:mxai-1,:) ,mom,cnt,var)	

        !call condmom(im,( corel(MNA:MXAD,:) .AND. emph(MNA:MXAD,:)==1 .AND. empw(MNA:MXAD,:)==0 .AND. move(MNA:MXAD,:)>=0 ),   d1* move(MNA:MXAD,:) ,mom,cnt,var)	
        !write(name(im),'("move | eu ",tr5)')  
        ! weights(im)=wmovemar   !ahu 122818 changed from wmove to wmovemar2
        !im=im+1
        !!call condmom(im,( corel(mna:mxai-1,:) .AND. dat(MNA:MXAI-1,:)%hhr==1 .AND. dat(MNA:MXAI-1,:)%hhsp==0  .AND. dat(MNA:MXAI-1,:)%sexr==1 .AND. move(mna:mxai-1,:)>=0 ),   d1* move(mna:mxai-1,:) ,mom,cnt,var)	

        !call condmom(im,( corel(MNA:MXAD,:) .AND. emph(MNA:MXAD,:)==0 .AND. empw(MNA:MXAD,:)==1 .AND. move(MNA:MXAD,:)>=0 ),   d1* move(MNA:MXAD,:) ,mom,cnt,var)	
        !write(name(im),'("move | ue ",tr5)')  
        ! weights(im)=wmovemar   !ahu 122818 changed from wmove to wmovemar2
        !im=im+1
        !!call condmom(im,( corel(mna:mxai-1,:) .AND. dat(MNA:MXAI-1,:)%hhr==0 .AND. dat(MNA:MXAI-1,:)%hhsp==1  .AND. dat(MNA:MXAI-1,:)%sexr==1 .AND. move(mna:mxai-1,:)>=0 ),   d1* move(mna:mxai-1,:) ,mom,cnt,var)	
    
        !call condmom(im,( corel(MNA:MXAD,:) .AND. emph(MNA:MXAD,:)==0 .AND. empw(MNA:MXAD,:)==0 .AND. move(MNA:MXAD,:)>=0 ),   d1* move(MNA:MXAD,:) ,mom,cnt,var)	
        !write(name(im),'("move | uu ",tr5)')  
        ! weights(im)=wmovemar   !ahu 122818 changed from wmove to wmovemar2
        !im=im+1
        !!call condmom(im,( corel(mna:mxai-1,:) .AND. dat(MNA:MXAI-1,:)%hhr==0 .AND. dat(MNA:MXAI-1,:)%hhsp==0  .AND. dat(MNA:MXAI-1,:)%sexr==1 .AND. move(mna:mxai-1,:)>=0 ),   d1* move(mna:mxai-1,:) ,mom,cnt,var)	
    

        !call condmom(im,(  corel(MNA:MXAD,:) .AND. edh(MNA:MXAD,:)==2  .AND. edw(MNA:MXAD,:)==2 .AND. move(MNA:MXAD,:)>=0 ),   d1* move(MNA:MXAD,:) ,mom,cnt,var)	
        !write(name(im),'("move|eded ",tr4)')  
        ! weights(im)=wmovemar   !ahu 122818 changed from wmove to wmovemar2
        !im=im+1
        !!call condmom(im,( corel(mna:mxai-1,:) .AND. dat(MNA:MXAI-1,:)%edr==2 .AND. dat(MNA:MXAI-1,:)%edsp==2  .AND. dat(MNA:MXAI-1,:)%sexr==1 .AND. move(mna:mxai-1,:)>=0 ),   d1* move(mna:mxai-1,:) ,mom,cnt,var)	

        !call condmom(im,(  corel(MNA:MXAD,:) .AND. edh(MNA:MXAD,:)==2  .AND. edw(MNA:MXAD,:)==1 .AND. move(MNA:MXAD,:)>=0 ),   d1* move(MNA:MXAD,:) ,mom,cnt,var)	
        !write(name(im),'("move|ednoed ",tr2)')  
        ! weights(im)=wmovemar  !ahu 122818 changed from wmove to wmovemar2
        !im=im+1
        !!call condmom(im,( corel(mna:mxai-1,:) .AND. dat(MNA:MXAI-1,:)%edr==2 .AND. dat(MNA:MXAI-1,:)%edsp==1  .AND. dat(MNA:MXAI-1,:)%sexr==1 .AND. move(mna:mxai-1,:)>=0 ),   d1* move(mna:mxai-1,:) ,mom,cnt,var)	
    
        !call condmom(im,(  corel(MNA:MXAD,:) .AND. edh(MNA:MXAD,:)==1  .AND. edw(MNA:MXAD,:)==2 .AND. move(MNA:MXAD,:)>=0 ),   d1* move(MNA:MXAD,:) ,mom,cnt,var)	
        !write(name(im),'("move|noeded ",tr2)')  
        ! weights(im)=wmovemar   !ahu 122818 changed from wmove to wmovemar2
        !im=im+1
        !!call condmom(im,( corel(mna:mxai-1,:) .AND. dat(MNA:MXAI-1,:)%edr==1 .AND. dat(MNA:MXAI-1,:)%edsp==2  .AND. dat(MNA:MXAI-1,:)%sexr==1 .AND. move(mna:mxai-1,:)>=0 ),   d1* move(mna:mxai-1,:) ,mom,cnt,var)	          
    
        !call condmom(im,(  corel(MNA:MXAD,:) .AND. edh(MNA:MXAD,:)==1  .AND. edw(MNA:MXAD,:)==1 .AND. move(MNA:MXAD,:)>=0 ),   d1* move(MNA:MXAD,:) ,mom,cnt,var)	
        !write(name(im),'("move|noednoed ")')  
        ! weights(im)=wmovemar   !ahu 122818 changed from wmove to wmovemar2
        !im=im+1
        !!call condmom(im,( corel(mna:mxai-1,:) .AND. dat(MNA:MXAI-1,:)%edr==1 .AND. dat(MNA:MXAI-1,:)%edsp==1  .AND. dat(MNA:MXAI-1,:)%sexr==1 .AND. move(mna:mxai-1,:)>=0 ),   d1* move(mna:mxai-1,:) ,mom,cnt,var)	            
        
            
        !ahu summer18 041118: the below getmar by gender and ia was added just to figure out the eps2 problem (why marriage rates were so sensitive to eps2 and 
        !why they were so different by gender). I figured it out now (i.e. the gender disparity in marriage rates go away when I set fem's psio equal to male's
        !but when I was looking at these getmar rates I also noticed the following: 
        !if I am recording age17 rel as 0 and when they make their decision at age 18 I record the rel for age18 as potentiall ymarried
        !but that is not consistent with data is it?
        !ahu jan19 010219: so I am still writing the age 17 (mad) mar rate. age 17 mar rate does not exit in the data since we don't record age 17 rel there. 
        !so the weight on this moment is 0 but then it's wrel for age=mna and age=19
        !do ia=mnad,mnad
        !    call condmom(im,( dat(ia,:)%rel==0 .AND. dat(ia+1,:)%rel>=0   .AND. dat(ia,:)%sexr==1 ), d1*one(dat(ia+1,:)%rel==1),mom,cnt,var)
        !    write(name(im),'("getmarbyia,m  ",tr1,i4)') ia
        !    weights(im)=0.0_dp
        !    im=im+1
        !    call condmom(im,( dat(ia,:)%rel==0 .AND. dat(ia+1,:)%rel>=0   .AND. dat(ia,:)%sexr==2 ), d1*one(dat(ia+1,:)%rel==1),mom,cnt,var)
        !    write(name(im),'("getmarbyia,f  ",tr1,i4)') ia
        !    weights(im)=0.0_dp
        !    im=im+1
        !end do            
    
        !Note about conditioning on age:
        !The max endage in data is 47. 
        !The way sim is done, for those whose endage is 47, the last age where variables get recorded is ia-1=46. 
        !This is because dat(ia-1,.) is recorded for each ia. So for the last age 47, the variables for 46 gets written
        !But then there is nothing after age 46, despite the fact that we do have people whose endage is 46 (namely 47). 

        !call condmom(im,( coho(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%rel==1  .AND. dat(MNA+1:MXA,:)%rel>=0 .AND. move(MNA:MXAD,:)==0 ), d1*one(dat(MNA+1:MXA,:)%rel==0),mom,cnt,var)
        !write(name(im),'("getdiv | nomv",tr1,i4)') 
        !weights(im)=wrel
        !im=im+1
        !call condmom(im,( coho(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%rel==1  .AND. dat(MNA+1:MXA,:)%rel>=0 .AND. move(MNA:MXAD,:)==1 ), d1*one(dat(MNA+1:MXA,:)%rel==0),mom,cnt,var)
        !write(name(im),'("getdiv |   mv",tr1,i4)') 
        !weights(im)=wrel
        !im=im+1

                                        

        !************************
        !ahu jan19 010119: commenting out the below and not doing the nummove by rel moments 
        !this is because for example in data nummove-0 is 0.84 while both nummove_mar and _sin are 0.92 something. 
        !I think the norelchg requirement is the reason. With that requirement it is not clear what these moments mean and 
        !they might be hindering the matching of the move by age moments. For example, there are times when nummove=0 is understated by a whole lot in sim
        !whereas movey age is overstated 
        !call condmom(im,(  cohogen(:)==co ) ,   d1* one( nummove_mar(MXAI,:)==0 ) ,mom,cnt,var)	
        !write(name(im),'("nummove_mar=0 ",tr1)')  
        !weights(im)=wmove
        !im=im+1
        !call condmom(im,(  cohogen(:)==co ),   d1* one( nummove_mar(MXAI,:)==1 ) ,mom,cnt,var)	
        !write(name(im),'("nummove_mar=1 ",tr1)')  
        !weights(im)=wmove
        !im=im+1
        
        !call condmom(im,(  cohogen(:)==co ),   d1* one( nummove_mar(MXAI,:)==2 ) ,mom,cnt,var)	
        !write(name(im),'("nummove_mar=2 ",tr1)')  
        !weights(im)=wmove
        !im=im+1
        
        !call condmom(im,(  cohogen(:)==co ),   d1* one( nummove_mar(MXAI,:)>=3 ) ,mom,cnt,var)	
        !write(name(im),'("nummove_mar>=3 ")')  
        !weights(im)=wmove
        !im=im+1

        !call condmom(im,(  cohogen(:)==co ) ,   d1* one( nummove_sin(MXAI,:)==0 ) ,mom,cnt,var)	
        !write(name(im),'("nummove_sin=0 ",tr1)')  
        !weights(im)=wmove
        !im=im+1

        !call condmom(im,(  cohogen(:)==co ),   d1* one( nummove_sin(MXAI,:)==1 ) ,mom,cnt,var)	
        !write(name(im),'("nummove_sin=1 ",tr1)')  
        !weights(im)=wmove
        !im=im+1
        
        !call condmom(im,(  cohogen(:)==co ),   d1* one( nummove_sin(MXAI,:)==2 ) ,mom,cnt,var)	
        !write(name(im),'("nummove_sin=2 ",tr1)')  
        !weights(im)=wmove
        !im=im+1
        
        !call condmom(im,(  cohogen(:)==co ),   d1* one( nummove_sin(MXAI,:)>=3 ) ,mom,cnt,var)	
        !write(name(im),'("nummove_sin>=3 ")')  
        !weights(im)=wmove
        !im=im+1
        !********************

 

             
            
            
            !do jj=1,2
            !do ia=mnad,MXAD,5
                !   call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr==0 .AND. dat(ia+1,:)%hhr>=0 .AND. move(ia,:)==1  .AND. dat(ia,:)%edr==jj ),   d1*one( dat(ia+1,:)%hhr==1 ),mom,cnt,var)		
                !   write(name(im),'("eumv",tr3,2i4)')  jj,ia
                !   weights(im)=0.0_dp  ; if (onlysingles.and.j==1) weights(im)=0.0_dp
                !   im=im+1 
            !end do 
            !end do 

        !do jj=1,2
        !    do ia=mnad,MXAD,5
        !        call condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr==1 .AND. dat(ia+1,:)%hhr>=0 .AND. move(ia,:)==1  .AND. dat(ia,:)%edr==jj ),   d1*one( dat(ia+1,:)%hhr==1 ),mom,cnt,var)		
    !		    write(name(im),'("eemv",tr3,2i4)')  jj,ia
    !		    weights(im)=0.0_dp  ; if (onlysingles.and.j==1) weights(im)=0.0_dp
        !           im=im+1                 
        !       end do 
        !       end do 


            
        !      do jj=1,2
        !          call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0  .AND. dat(MNA:MXAD,:)%edr==jj ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
    !		    write(name(im),'("e | u stay edu",tr3,i4)')  jj
!			    weights(im)=0.0_dp  ; if (onlysingles.and.j==1) weights(im)=0.0_dp
!                im=im+1 

    !               call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==0 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1  .AND. dat(MNA:MXAD,:)%edr==jj ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
    !		    write(name(im),'("e | u move edu",tr3,i4)')  jj
    !		    weights(im)=0.0_dp  ; if (onlysingles.and.j==1) weights(im)=0.0_dp
    !            im=im+1 

        !           call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==0  .AND. dat(MNA:MXAD,:)%edr==jj ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
    !		    write(name(im),'("e | e stay edu",tr3,i4)')  jj
    !		    weights(im)=0.0_dp  ; if (onlysingles.and.j==1) weights(im)=0.0_dp
        !           im=im+1 
            
        !          call condmom(im,( cosexrel(MNA:MXAD,:) .AND. dat(MNA:MXAD,:)%hhr==1 .AND. dat(MNA+1:MXA,:)%hhr>=0 .AND. move(MNA:MXAD,:)==1  .AND. dat(MNA:MXAD,:)%edr==jj ),   d1*one( dat(MNA+1:MXA,:)%hhr==1 ),mom,cnt,var)		
        !	    write(name(im),'("e | e move edu",tr3,i4)')  jj
        !	    weights(im)=0.0_dp  ; if (onlysingles.and.j==1) weights(im)=0.0_dp
            !       im=im+1                 
            !   end do 
            
            

            !do ia=MNAD,25 
        !		CALL condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr>=0 .and. dat(ia,:)%edr==1) ,d1*one(dat(ia,:)%hhr==1),mom,cnt,var)
        !		WRITE(name(im),'("emp a,ned",tr4,I2)') ia
        !		weights(im)=whour ; if (onlysingles.and.j==1) weights(im)=0.0_dp
            !       if (ia<mna) weights(im)=0.0_dp
        !		im=im+1
        !	end do            
            !   do ia=MNAD,25 
            !	CALL condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr>=0 .and. dat(ia,:)%edr==2) ,d1*one(dat(ia,:)%hhr==1),mom,cnt,var)
        !		WRITE(name(im),'("emp a,ed",tr4,I2)') ia
        !		weights(im)=whour ; if (onlysingles.and.j==1) weights(im)=0.0_dp
            !       if (ia<agestart(college)) weights(im)=0.0_dp
            !      im=im+1
        !	end do            
            !do ia=mna,45 !mxai-1
            !	CALL condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr>=0  .AND. dat(ia,:)%edr==1 ) ,d1*one(dat(ia,:)%hhr==1),mom,cnt,var)
            !	WRITE(name(im),'("e|ned by age",tr2,I2)') ia
            !	weights(im)=whour ; if (onlysingles.and.j==1) weights(im)=0.0_dp
            !	im=im+1
            !end do            

            !do ia=mna,45 !mxai-1
            !	CALL condmom(im,( cosexrel(ia,:) .AND. dat(ia,:)%hhr>=0  .AND. dat(ia,:)%edr==2 ) ,d1*one(dat(ia,:)%hhr==1),mom,cnt,var)
            !	WRITE(name(im),'("e| ed by age",tr2,I2)') ia
            !	weights(im)=whour ; if (onlysingles.and.j==1) weights(im)=0.0_dp
            !	im=im+1
            !end do


            !do i=1,nl
            !    call condmom(im,( cosexrel(mna:mxai,:) .AND. dat(mna:mxai,:)%hhr>=0 .AND. dat(mna:mxai,:)%l==i  .AND. dat(mna:mxai,:)%edr==1 ),   d1*one( dat(mna:mxai,:)%hhr==1 ),mom,cnt,var)		
            !    write(name(im),'("e|loc noed",tr4,i4)') i			
            !    weights(im)=whour  ; if (onlysingles.and.j==1) weights(im)=0.0_dp
            !    im=im+1 
            !end do 
        
            !do i=1,nl
            !    call condmom(im,( cosexrel(mna:mxai,:) .AND. dat(mna:mxai,:)%hhr==1 .AND. dat(mna:mxai,:)%l==i  .AND. dat(mna:mxai,:)%edr==1 .AND. dat(mna:mxai,:)%logwr>=0 ),   d1*dat(mna:mxai,:)%logwr ,mom,cnt,var)		
            !    write(name(im),'("w|loc noed",tr4,i4)') i			
            !    weights(im)=wwage ; if (onlysingles.and.j==1) weights(im)=0.0_dp
            !    im=im+1 
            !end do 

            !do i=1,nl
            !    call condmom(im,( cosexrel(mna:mxai,:) .AND. dat(mna:mxai,:)%hhr>=0 .AND. dat(mna:mxai,:)%l==i  .AND. dat(mna:mxai,:)%edr==2 ),   d1*one( dat(mna:mxai,:)%hhr==1 ),mom,cnt,var)		
            !    write(name(im),'("e|loc   ed",tr4,i4)') i			
            !    weights(im)=whour  ; if (onlysingles.and.j==1) weights(im)=0.0_dp
            !    im=im+1 
            !end do 
        
            !do i=1,nl
            !    call condmom(im,( cosexrel(mna:mxai,:) .AND. dat(mna:mxai,:)%hhr==1 .AND. dat(mna:mxai,:)%l==i  .AND. dat(mna:mxai,:)%edr==2 .AND. dat(mna:mxai,:)%logwr>=0 ),   d1*dat(mna:mxai,:)%logwr ,mom,cnt,var)		
            !    write(name(im),'("w|loc   ed",tr4,i4)') i			
            !    weights(im)=wwage ; if (onlysingles.and.j==1) weights(im)=0.0_dp
            !    im=im+1 
            !end do 
    

        

            
            
            
        

        !ahu summer18: before the below moments were conditioned on cosexrel instead of just cosex. 
        !I remove the rel conditioning because sometimes these seem to have noone in the cells and they lead to jumpiness 
        !See notes for 042118. 
        !ia=MNAd
        !call condmom(im,( cosex(ia,:) .AND. dat(ia,:)%hhr==0 .AND. dat(ia+1,:)%hhr>=0 ),   d1*one( dat(ia+1,:)%hhr==1 ),mom,cnt,var)		
        !write(name(im),'("e|u by ia",tr5,i2)') ia
        !weights(im)=0.0_dp
        !im=im+1 
        !call condmom(im,( cosex(ia,:) .AND. dat(ia,:)%hhr==0 .AND. dat(ia+1,:)%hhr==1 .AND. dat(ia+1,:)%logwr>=0  ),   d1*dat(ia+1,:)%logwr ,mom,cnt,var)		
        !write(name(im),'("w|u by ia",tr5,i2)') ia
        !weights(im)=0.0_dp 
        !im=im+1 

            !do ia=MNA,MXAD,8
            !    call condmom(im,( cosex(ia,:) .AND. dat(ia,:)%hhr==0 .AND. dat(ia+1,:)%hhr>=0 ),   d1*one( dat(ia+1,:)%hhr==1 ),mom,cnt,var)		
            !    write(name(im),'("e|u by ia",tr5,i2)') ia
            !    weights(im)=0.0_dp  
            !    im=im+1 

            !    call condmom(im,( cosex(ia,:) .AND. dat(ia,:)%hhr==0 .AND. dat(ia+1,:)%hhr==1 .AND. dat(ia+1,:)%logwr>=0  ),   d1*dat(ia+1,:)%logwr ,mom,cnt,var)		
            !    write(name(im),'("w|u by ia",tr5,i2)') ia
            !    weights(im)=0.0_dp
            !    im=im+1 
            !end do 


        !********************
        !Everyone MISC



            
            !if (yaz) then
            !if (ndat==nsim) then
            !    DO j=1,nsim
            !        ! write(12,*) ' id  ia sex rel  ed edp kid   hh  wage    dd incsp ln je cm sc sm kt ds O3'
            !        DO ia=MNA,MXAD
            !            write(12,'(7I12,6F25.5,7I12,2F25.5,9I12)') j,ia,dat(ia,j)%co,dat(ia,j)%sexr,&
            !            & dat(ia,j)%rel,dat(ia,j)%kidr,dat(ia,j)%edr, &
            !            & dat(ia,j)%logwr,dat(ia,j)%logwsp,dat(ia,j)%hhr, & 
            !            & dat(ia,j)%hhsp,dat(ia,j)%ddr,dat(ia,j)%ddsp,&
            !            & dat(ia,j)%rellen,-1,kidtrans(ia,j),& 
            !            & disolve(ia,j),ObsLast3(ia,j),dat(ia,j)%l,& 
            !            & move(ia,j),deltawage(ia,j),deltawage4(ia,j),etr(1:5,ia,j),moverank(ia,j),homemove(ia,j),& 
            !            & norelchg(ia,j)
            !        ENDDO
            !    ENDDO 
            !end if 
            !end if 
            
            
            
            
            !headloc(ihead)=im; headstr(ihead)='Move by total number of moves (all ages)';ihead=ihead+1
            !movesum(:)=sum(move(MNA:MXAI-1,:),1,move>=0)  
            !do j=0,4
            !	CALL condmom(im,((dat%co==co)),d1*move(MNA:MXAI-1,:),mom,cnt,var)
                !CALL condmom(im,((movesum>=0)),d1*one(movesum==j),mom,cnt,var)
            !	WRITE(name(im),'("movesum ",I4)') j
            !	weights(im)=wmove
            !	im=im+1 
            !end do 
            !do j=0,4
            !	CALL condmom(im,((movesum>=0)),d1*one(movesum==j),mom,cnt,var)
            !	WRITE(name(im),'("movesum ",I4)') j
            !	weights(im)=wmove
            !	im=im+1 
            !end do 

            !headloc(ihead)=im; headstr(ihead)='Move by whether they are working full time and gender (all ages) ';ihead=ihead+1
            !do g=1,2
            !    CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).and.(norelchg==1).AND.(dat%hhr>=0).AND.(dat%hhr<H_FULLTIME).AND.(move>=0)),d1*move,mom,cnt,var)
            !    WRITE(name(im),'("move|u ",I4)') g
            !    weights(im)=wmove0
            !    im=im+1
            !    CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).and.(norelchg==1).AND.(dat%hhr>=0).AND.(dat%hhr>=H_FULLTIME).AND.(move>=0)),d1*move,mom,cnt,var)
            !    WRITE(name(im),'("move|e ",I4)') g
            !    weights(im)=wmove0
            !    im=im+1
            !end do 
            !headloc(ihead)=im; headstr(ihead)='Move by whether they are working full time and gender and rel (all ages) ';ihead=ihead+1
            !do g=1,2
            !    do j=0,1
            !        CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==j).and.(norelchg==1).AND.(dat%hhr>=0).AND.(dat%hhr<H_FULLTIME).AND.(move>=0)),d1*move,mom,cnt,var)
            !        WRITE(name(im),'("move|u ",2I4)') g,j
            !        if (j==0) weights(im)=wmove0 
            !        if (j==1) weights(im)=wmove1 
            !        im=im+1
            !        CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==j).and.(norelchg==1).AND.(dat%hhr>=0).AND.(dat%hhr>=H_FULLTIME).AND.(move>=0)),d1*move,mom,cnt,var)
            !        WRITE(name(im),'("move|e ",2I4)') g,j
            !        if (j==0) weights(im)=wmove0 
            !        if (j==1) weights(im)=wmove1 
            !        im=im+1
            !    end do 
            !end do 
        



                !headloc(ihead)=im; headstr(ihead)='Mean wage by gender/loc - FULL TIME (all ages)';ihead=ihead+1
                !do g=1,2
                !    do j=1,NL
                !        CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%l==j).AND.(dat%hhr>=H_FULLTIME).AND.(dat%logwr>=0).AND.(iacat==1)),d1*dat%logwr,mom,cnt,var)
                !        WRITE(name(im),'("logwage/gender/loc ",2I4)') g,j
                !        weights(im)=wwage0
                !        im=im+1
                !    end do 
                !end do 
            
                !headloc(ihead)=im; headstr(ihead)='Mean wage by gender/loc/ed - FULL TIME (all ages)';ihead=ihead+1
                !do g=1,2
                !    do j=1,NL
                !        CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%edr==1).AND.(dat%l==j).AND.(dat%hhr>=H_FULLTIME).AND.(dat%logwr>=0).AND.(iacat==1)),d1*dat%logwr,mom,cnt,var)
                !        WRITE(name(im),'("logwage/gender/loc | hs ",2I4)') g,j
                !        weights(im)=wwage0
                !        im=im+1
                !        CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%edr==2).AND.(dat%l==j).AND.(dat%hhr>=H_FULLTIME).AND.(dat%logwr>=0).AND.(iacat==1)),d1*dat%logwr,mom,cnt,var)
                !        WRITE(name(im),'("logwage/gender/loc | col  ",2I4)') g,j
                !        weights(im)=wwage0
                !        im=im+1
                !    end do 
                !end do 

                ! wages-growth interacted with hours worked and hours squared by gender
                !headloc(ihead)=im; headstr(ihead)='wages-growth, interacted with hours worked and hours squared by gender';ihead=ihead+1
                !CALL condmom(im,((dat(MNA+4:MXA,:)%co==co).AND.(dat(MNA+4:MXA,:)%sexr==1).AND.(dat(MNA+4:MXA,:)%logwr>=0).AND.(dat(MNA:MXA-4,:)%logwr>=0).AND.(deltawage4(MNA+4:MXA,:)>=0)),d1*deltawage4(MNA+4:MXA,:),mom,cnt,var)
                !WRITE(name(im),'("wage-growth/male             ")') 
                !CALL condmom(im+1,((dat(MNA+4:MXA,:)%co==co).AND.(dat(MNA+4:MXA,:)%sexr==1).AND.(dat(MNA+4:MXA,:)%logwr>=0).AND.(dat(MNA:MXA-4,:)%logwr>=0).AND.(deltawage4(MNA+4:MXA,:)>=0).AND.(mean4h(MNA+4:MXA,:)>=0)),d1*deltawage4(MNA+4:MXA,:)*mean4h(MNA+4:MXA,:),mom,cnt,var)
                !WRITE(momentsname(im+1),'("wage-growth-hours/male     ")') 
                !CALL condmom(im+2,((dat(MNA+4:MXA,:)%co==co).AND.(dat(MNA+4:MXA,:)%sexr==2).AND.(dat(MNA+4:MXA,:)%logwr>=0).AND.(dat(MNA:MXA-4,:)%logwr>=0).AND.(deltawage4(MNA+4:MXA,:)>=0)),d1*deltawage4(MNA+4:MXA,:),mom,cnt,var)
                !WRITE(momentsname(im+2),'("wage-growth/female         ")') 
                !CALL condmom(im+3,((dat(MNA+4:MXA,:)%co==co).AND.(dat(MNA+4:MXA,:)%sexr==2).AND.(dat(MNA+4:MXA,:)%logwr>=0).AND.(dat(MNA:MXA-4,:)%logwr>=0).AND.(deltawage4(MNA+4:MXA,:)>=0).AND.(mean4h(MNA+4:MXA,:)>=0)),d1*deltawage4(MNA+4:MXA,:)*mean4h(MNA+4:MXA,:),mom,cnt,var)
                !WRITE(momentsname(im+3),'("wage-growth-hours/female   ")') 
                !weights(im:im+3)=0.0_dp !wwage0
                !im=im+4
            
            
                !ahu 082012 do j=1,NL
                !ahu 082012 CALL condmom(im,((dat(MNA:MXAI-1,:)%co==co).AND.(dat(MNA:MXAI-1,:)%l==j).AND.(move(MNA:MXAI-1,:)>=0)),d1*move(MNA:MXAI-1,:),mom,cnt,var)
                !ahu 082012 WRITE(name(im),'("move-rates-from-loc",3I4)') 0,0,j
                !ahu 082012 weights(im)=wmove
                !ahu 082012 im=im+1 
                !ahu 082012 end do 
                !ahu 082012 do j=1,NL
                !ahu 082012     CALL condmom(im,((dat(MNA:MXAI-1,:)%co==co).AND.(dat(MNA+1:MXAI,:)%l==j).AND.(move(MNA:MXAI-1,:)>=0)),d1*move(MNA:MXAI-1,:),mom,cnt,var)
                !ahu 082012     WRITE(name(im),'("move-rates-to-loc",3I4)') 0,0,j
                !ahu 082012     weights(im)=wmove
                !ahu 082012     im=im+1 
                !ahu 082012 end do 



            !do i=1,NL
            !    do j=1,NL
            !        CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(move(MNA:MXAD,:)==1).AND.(dat(MNA:MXAD,:)%l==i)),d1*one(dat(MNA+1:MXA,:)%l==j),mom,cnt,var)
            !        WRITE(name(im),'("matloc                    ",2I4)') j,i
            !        weights(im)=wmove0
            !        im=im+1 
            !    end do     
            !end do 
            !do i=1,NL
            !CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(move(MNA:MXAD,:)>=0).AND.(dat(MNA:MXAD,:)%l==i)),d1*one(move(MNA:MXAD,:)==1),mom,cnt,var)
            !    WRITE(name(im),'("move | loc                ",I4)') i
            !    weights(im)=wmove0
            !    im=im+1 
            !end do 
        !headloc(ihead)=im; headstr(ihead)='Labor market hours by gender/rel/ia';ihead=ihead+1
        !ahu 061211: have to control for ia here because the two brs have different ia compositions
        !ahu 061211: br 2 has no hours/kids/cohmar simultaneously in the biannual years so if you condition on all that you will just get something until they are ia 28 or something (depending on what the br grouping is)
        !ahu 061211: and so if we don't control for ia, it looks as if br 2 females who are cohabiting have decreased their hours of work. but this is just a composition effect.
        !ahu 061211: excluding ia 20 because, something looks weird. br 2 works too few hours at ia 20 (for females,coh,nokid). so then when I include them, it looks as if br 2 coh females with no kids work less in the later br. 
        !do g=1,2
        !   do j=0,1
        !       CALL condmom(im,((dat%co==co).AND.(dat%sexr==g).AND.(dat%rel==j).AND.(dat%hhr>=0).AND.(iacat==1)),dat%hhr,mom,cnt,var)
        !       WRITE(name(im),'("hrs/gender/rel ",2I4)') g,j
        !       weights(im)=whour
        !       im=im+1 
        !   end do 
        !end do 
            

            !headloc(ihead)=im; headstr(ihead)='mean log wage by gender/movesum_rel/ia (FULL TIME)';ihead=ihead+1       
            !do g=1,2
            !    do j=0,1
            !        do ia=MNA+1,MXAI-4,4
            !            CALL condmom(im,((dat(ia,:)%co==co).AND.(dat(ia,:)%sexr==g).AND.(dat(ia,:)%rel==0).AND.(movesum_single(ia,:)==j).AND.(dat(ia,:)%hhr>=H_FULLTIME).AND.(dat(ia,:)%logwr>=0)),d1*dat(ia,:)%logwr,mom,cnt,var)
            !            WRITE(name(im),'("mean-log-wage-by-gender/movesing/ia ",3I4)') g,j,ia
            !            weights(im)=0.0_rp
            !            im=im+1
            !        end do 
            !    end do 
            !end do 
            !do g=1,2
            !    do j=0,1
            !        do ia=MNA+1,MXAI-4,4
            !            CALL condmom(im,((dat(ia,:)%co==co).AND.(dat(ia,:)%sexr==g).AND.(dat(ia,:)%rel==1).AND.(movesum_mar(ia,:)==j).AND.(dat(ia,:)%hhr>=H_FULLTIME).AND.(dat(ia,:)%logwr>=0)),d1*dat(ia,:)%logwr,mom,cnt,var)
            !            WRITE(name(im),'("mean-log-wage-by-gender/movemar/ia  ",3I4)') g,j,ia
            !            weights(im)=0.0_rp
            !            im=im+1
            !        end do 
            !    end do 
            !end do 
            !headloc(ihead)=im; headstr(ihead)='Move to rates by gender and relstat ';ihead=ihead+1
            !do j=1,NL
            !    CALL condmom(im,((dat(MNA:MXAI-1,:)%co==co).AND.(dat(MNA:MXAI-1,:)%sexr==1).AND.(dat(MNA:MXAI-1,:)%rel==0).AND.(move(MNA:MXAI-1,:)==1).AND.(dat(MNA+1:MXAI,:)%l>=0)),d1*one(dat(MNA+1:MXAI,:)%l==j),mom,cnt,var)
            !    WRITE(name(im),'("male,single,prop of moves to each loc ",I2)') j
            !    weights(im)=0.0_rp
            !    im=im+1 
            !end do     
            !do j=1,NL
            !    CALL condmom(im,((dat(MNA:MXAI-1,:)%co==co).AND.(dat(MNA:MXAI-1,:)%sexr==1).AND.(dat(MNA:MXAI-1,:)%rel==1).AND.(move(MNA:MXAI-1,:)==1).AND.(dat(MNA+1:MXAI,:)%l>=0)),d1*one(dat(MNA+1:MXAI,:)%l==j),mom,cnt,var)
            !    WRITE(name(im),'("male,   mar,prop of moves to each loc ",I2)') j
            !    weights(im)=0.0_rp
            !    im=im+1 
            !end do     
            !do j=1,NL
            !    CALL condmom(im,((dat(MNA:MXAI-1,:)%co==co).AND.(dat(MNA:MXAI-1,:)%sexr==2).AND.(dat(MNA:MXAI-1,:)%rel==0).AND.(move(MNA:MXAI-1,:)==1).AND.(dat(MNA+1:MXAI,:)%l>=0)),d1*one(dat(MNA+1:MXAI,:)%l==j),mom,cnt,var)
            !    WRITE(name(im),'("fem, single,prop of moves to each loc ",I2)') j
            !    weights(im)=0.0_rp
            !    im=im+1 
            !end do     
            !do j=1,NL
            !    CALL condmom(im,((dat(MNA:MXAI-1,:)%co==co).AND.(dat(MNA:MXAI-1,:)%sexr==2).AND.(dat(MNA:MXAI-1,:)%rel==1).AND.(move(MNA:MXAI-1,:)==1).AND.(dat(MNA+1:MXAI,:)%l>=0)),d1*one(dat(MNA+1:MXAI,:)%l==j),mom,cnt,var)
            !    WRITE(name(im),'("fem,    mar,prop of moves to each loc ",I2)') j
            !    weights(im)=0.0_rp
            !    im=im+1 
            !end do     


        ! wages for working women interacted with fraction of last 4 years out of labor force, by ed
            !headloc(ihead)=im; headstr(ihead)='wages for working women interacted with fraction of last 4 years out of labor force, by ed';ihead=ihead+1
            !DO ia=22,34,4
            !   CALL condmom(im,((dat(ia,:)%co==co).AND.(dat(ia,:)%sexr==2).AND.(dat(ia,:)%logwr>0).AND.(dat(ia,:)%edr==1).AND.(obs4h(ia,:))),d1*dat(ia,:)%logwr*frac4h0(ia,:),mom,cnt,var)
            !   WRITE(name(im),'("logw-years out, fem hs,",1I2)') ia
            !   weights(im)=0.
            !   im=im+1
            !ENDDO
            !DO ia=26,34,4
            !   CALL condmom(im,((dat(ia,:)%co==co).AND.(dat(ia,:)%sexr==2).AND.(dat(ia,:)%logwr>0).AND.(dat(ia,:)%edr==2).AND.(obs4h(ia,:))),d1*dat(ia,:)%logwr*frac4h0(ia,:),mom,cnt,var)
            !   WRITE(name(im),'("logw-years out, fem col,",1I2)') ia
            !   weights(im)=0.
            !   im=im+1
            !ENDDO

            !do ia=19,MXAD,4
            !    CALL condmom(im,((dat(ia,:)%co==co).AND.(dat(ia,:)%rel==0).and.(norelchg(ia,:)==1).AND.(move(ia,:)>=0)),d1*move(ia,:),mom,cnt,var)
            !    WRITE(name(im),'("move/single ",I4)') ia
            !    weights(im)=0.0_dp !wmove0              
            !    im=im+1 
            !end do 
        
            !do ia=19,MXAD,4
            !    CALL condmom(im,((dat(ia,:)%co==co).AND.(dat(ia,:)%rel==1).and.(norelchg(ia,:)==1).AND.(move(ia,:)>=0)),d1*move(ia,:),mom,cnt,var)
            !    WRITE(name(im),'("move/married ",I4)') ia
            !    weights(im)=0.0_dp !wmove1              
            !    im=im+1 ! married
            !end do 
        
            !headloc(ihead)=im; headstr(ihead)='Mean wage by gender/ed/ia - FULL TIME (all ages)';ihead=ihead+1
            !do g=1,2
            !    do j=1,2    
            !        do ia=MNAD,24
            !            CALL condmom(im,((dat(ia,:)%co==co).AND.(dat(ia,:)%sexr==g).AND.(dat(ia,:)%edr==j).AND.(dat(ia,:)%hhr>=H_FULLTIME).AND.(dat(ia,:)%logwr>=0)),d1*dat(ia,:)%logwr,mom,cnt,var)
            !            WRITE(name(im),'("logw/gender/ed/ia",3I4)') g,j,ia
            !            weights(im)=0.0_dp !wwage0
            !            im=im+1
            !        end do 
            !        do ia=24,MXA,4
            !            CALL condmom(im,((dat(ia,:)%co==co).AND.(dat(ia,:)%sexr==g).AND.(dat(ia,:)%edr==j).AND.(dat(ia,:)%hhr>=H_FULLTIME).AND.(dat(ia,:)%logwr>=0)),d1*dat(ia,:)%logwr,mom,cnt,var)
            !            WRITE(name(im),'("logw/gender/ed/ia",3I4)') g,j,ia
            !            weights(im)=0.0_dp !wwage0
            !            im=im+1
            !        end do 
            !    end do 
            !end do 



            !headloc(ihead)=im; headstr(ihead)='mean log wage by gender/ia - FULL TIME';ihead=ihead+1
            !do g=1,2
            !DO ia=MNA,MXA,4
            !    CALL condmom(im,((dat(ia,:)%co==co).AND.(dat(ia,:)%sexr==g).AND.(dat(ia,:)%hhr>=H_FULLTIME).AND.(dat(ia,:)%logwr>=0)),d1*dat(ia,:)%logwr,mom,cnt,var)
            !    WRITE(name(im),'("mean-log-wage/gender/ia ",2I4)') g,ia
            !    weights(im)=0.0_dp !wwage0
            !    im=im+1
            !ENDDO
            !end do 



            !THESE DON'T WORK. PROBABLY NEED TO CHANGE FULLTIME VARIABLE IN GLOBAL
            
            headloc(ihead)=im ; headstr(ihead)='everyone misc' ; ihead=ihead+1
            headloc(ihead)=im; headstr(ihead)='U2EMP by gender/rel - FULL TIME (all ages)';ihead=ihead+1
            do g=1,2
                do j=0,1
                    CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%sexr==g).AND.(dat(MNA:MXAD,:)%rel==j).and.(norelchg(MNA:MXAD,:)==1).AND.(etr(1,MNA:MXAD,:)>=0).AND.(iacat(MNA:MXAD,:)==1)),d1*etr(1,MNA:MXAD,:),mom,cnt,var)
                    WRITE(name(im),'("u2emp          ",2I4)') g,j
                    if (j==0) weights(im)=0.0_dp !whour0 
                    if (j==1) weights(im)=0.0_dp !whour1 
                    im=im+1 
                    CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%sexr==g).AND.(dat(MNA:MXAD,:)%rel==j).and.(norelchg(MNA:MXAD,:)==1).AND.(etr(1,MNA:MXAD,:)>=0).AND.(move(MNA:MXAD,:)==0).AND.(iacat(MNA:MXAD,:)==1)),d1*etr(1,MNA:MXAD,:),mom,cnt,var)
                    WRITE(name(im),'("u2emp | stay   ",2I4)') g,j
                    if (j==0) weights(im)=0.0_dp !whour0 
                    if (j==1) weights(im)=0.0_dp !whour1 
                    im=im+1 
                    CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%sexr==g).AND.(dat(MNA:MXAD,:)%rel==j).and.(norelchg(MNA:MXAD,:)==1).AND.(etr(1,MNA:MXAD,:)>=0).AND.(move(MNA:MXAD,:)==1).AND.(iacat(MNA:MXAD,:)==1)),d1*etr(1,MNA:MXAD,:),mom,cnt,var)
                    WRITE(name(im),'("u2emp | move   ",2I4)') g,j
                    if (j==0) weights(im)=0.0_dp !whour0 
                    if (j==1) weights(im)=0.0_dp !whour1 
                    im=im+1 
                    CALL condmom(im,((dat(MNAD:MXA,:)%co==co).AND.(dat(MNAD:MXA,:)%sexr==g).AND.(dat(MNAD:MXA,:)%rel==j).and.(norelchg(MNAD:MXA,:)==1).AND.(etr(1,MNAD:MXA,:)>=0).AND.(iacat(MNAD:MXA,:)==1)),d1*etr(1,MNAD:MXA,:),mom,cnt,var)
                    WRITE(name(im),'("u2emp          ",2I4)') g,j
                    if (j==0) weights(im)=0.0_dp !whour0 
                    if (j==1) weights(im)=0.0_dp !whour1 
                    im=im+1 
                    CALL condmom(im,((dat(MNAD:MXA,:)%co==co).AND.(dat(MNAD:MXA,:)%sexr==g).AND.(dat(MNAD:MXA,:)%rel==j).and.(norelchg(MNAD:MXA,:)==1).AND.(etr(1,MNAD:MXA,:)>=0).AND.(move(MNAD:MXA,:)==0).AND.(iacat(MNAD:MXA,:)==1)),d1*etr(1,MNAD:MXA,:),mom,cnt,var)
                    WRITE(name(im),'("u2emp | stay   ",2I4)') g,j
                    if (j==0) weights(im)=0.0_dp !whour0 
                    if (j==1) weights(im)=0.0_dp !whour1 
                    im=im+1 
                    CALL condmom(im,((dat(MNAD:MXA,:)%co==co).AND.(dat(MNAD:MXA,:)%sexr==g).AND.(dat(MNAD:MXA,:)%rel==j).and.(norelchg(MNAD:MXA,:)==1).AND.(etr(1,MNAD:MXA,:)>=0).AND.(move(MNAD:MXA,:)==1).AND.(iacat(MNAD:MXA,:)==1)),d1*etr(1,MNAD:MXA,:),mom,cnt,var)
                    WRITE(name(im),'("u2emp | move   ",2I4)') g,j
                    if (j==0) weights(im)=0.0_dp !whour0 
                    if (j==1) weights(im)=0.0_dp !whour1 
                    im=im+1 
                end do !rel j
            end do  !sex g



            !THESE DON'T WORK. PROBABLY NEED TO CHANGE FULLTIME VARIABLE IN GLOBAL
            
            headloc(ihead)=im ; headstr(ihead)='everyone misc' ; ihead=ihead+1
            headloc(ihead)=im; headstr(ihead)='EMP2EMP by gender/rel - FULL TIME (all ages)';ihead=ihead+1
            do g=1,2
                do j=0,1
                    CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%sexr==g).AND.(dat(MNA:MXAD,:)%rel==j).and.(norelchg(MNA:MXAD,:)==1).AND.(etr(2,MNA:MXAD,:)>=0).AND.(iacat(MNA:MXAD,:)==1)),d1*etr(2,MNA:MXAD,:),mom,cnt,var)
                    WRITE(name(im),'("emp2emp        ",2I4)') g,j
                    if (j==0) weights(im)=0.0_dp !whour0 
                    if (j==1) weights(im)=0.0_dp !whour1 
                    im=im+1 
                    CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%sexr==g).AND.(dat(MNA:MXAD,:)%rel==j).and.(norelchg(MNA:MXAD,:)==1).AND.(etr(2,MNA:MXAD,:)>=0).AND.(move(MNA:MXAD,:)==0).AND.(iacat(MNA:MXAD,:)==1)),d1*etr(2,MNA:MXAD,:),mom,cnt,var)
                    WRITE(name(im),'("emp2emp | stay ",2I4)') g,j
                    if (j==0) weights(im)=0.0_dp !whour0 
                    if (j==1) weights(im)=0.0_dp !whour1 
                    im=im+1 
                    CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%sexr==g).AND.(dat(MNA:MXAD,:)%rel==j).and.(norelchg(MNA:MXAD,:)==1).AND.(etr(2,MNA:MXAD,:)>=0).AND.(move(MNA:MXAD,:)==1).AND.(iacat(MNA:MXAD,:)==1)),d1*etr(2,MNA:MXAD,:),mom,cnt,var)
                    WRITE(name(im),'("emp2emp | move ",2I4)') g,j
                    if (j==0) weights(im)=0.0_dp !whour0 
                    if (j==1) weights(im)=0.0_dp !whour1 
                    im=im+1 
        
                    CALL condmom(im,((dat(MNAD:MXA,:)%co==co).AND.(dat(MNAD:MXA,:)%sexr==g).AND.(dat(MNAD:MXA,:)%rel==j).and.(norelchg(MNAD:MXA,:)==1).AND.(etr(2,MNAD:MXA,:)>=0).AND.(iacat(MNAD:MXA,:)==1)),d1*etr(2,MNAD:MXA,:),mom,cnt,var)
                    WRITE(name(im),'("emp2emp        ",2I4)') g,j
                    if (j==0) weights(im)=0.0_dp !whour0 
                    if (j==1) weights(im)=0.0_dp !whour1 
                    im=im+1 
                    CALL condmom(im,((dat(MNAD:MXA,:)%co==co).AND.(dat(MNAD:MXA,:)%sexr==g).AND.(dat(MNAD:MXA,:)%rel==j).and.(norelchg(MNAD:MXA,:)==1).AND.(etr(2,MNAD:MXA,:)>=0).AND.(move(MNAD:MXA,:)==0).AND.(iacat(MNAD:MXA,:)==1)),d1*etr(2,MNAD:MXA,:),mom,cnt,var)
                    WRITE(name(im),'("emp2emp | stay ",2I4)') g,j
                    if (j==0) weights(im)=0.0_dp !whour0 
                    if (j==1) weights(im)=0.0_dp !whour1 
                    im=im+1 
                    CALL condmom(im,((dat(MNAD:MXA,:)%co==co).AND.(dat(MNAD:MXA,:)%sexr==g).AND.(dat(MNAD:MXA,:)%rel==j).and.(norelchg(MNAD:MXA,:)==1).AND.(etr(2,MNAD:MXA,:)>=0).AND.(move(MNAD:MXA,:)==1).AND.(iacat(MNAD:MXA,:)==1)),d1*etr(2,MNAD:MXA,:),mom,cnt,var)
                    WRITE(name(im),'("emp2emp | move ",2I4)') g,j
                    if (j==0) weights(im)=0.0_dp !whour0 
                    if (j==1) weights(im)=0.0_dp !whour1 
                    im=im+1         
                end do 
            end do 
            


            !THESE DON'T WORK. PROBABLY NEED TO CHANGE FULLTIME VARIABLE IN GLOBAL
            
            headloc(ihead)=im ; headstr(ihead)='everyone misc' ; ihead=ihead+1
            headloc(ihead)=im; headstr(ihead)='U2EMP by gender/home - SINGLES - FULL TIME (all ages)';ihead=ihead+1
            do g=1,2
                CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%sexr==g).AND.(dat(MNA:MXAD,:)%rel==0).and.(norelchg(MNA:MXAD,:)==1).AND.(move(MNA:MXAD,:)==0).AND.(etr(1,MNA:MXAD,:)>=0).AND.(iacat(MNA:MXAD,:)==1)),d1*etr(1,MNA:MXAD,:),mom,cnt,var)
                WRITE(name(im),'("u2emp | stay          ",I2)') g
                weights(im)=0.0_dp !0.0_rp
                im=im+1 
                CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%sexr==g).AND.(dat(MNA:MXAD,:)%rel==0).and.(norelchg(MNA:MXAD,:)==1).AND.(move(MNA:MXAD,:)==1).AND.(dat(MNA:MXAD,:)%hme==dat(MNA:MXAD,:)%l).AND.(etr(5,MNA:MXAD,:)>=0).AND.(iacat(MNA:MXAD,:)==1)),d1*etr(5,MNA:MXAD,:),mom,cnt,var)
                WRITE(name(im),'("u2u | move(MNA:MXAD,:),athme1    ",I2)') g
                weights(im)=0.0_dp !0.0_rp
                im=im+1 
                CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%sexr==g).AND.(dat(MNA:MXAD,:)%rel==0).and.(norelchg(MNA:MXAD,:)==1).AND.(move(MNA:MXAD,:)==1).AND.(dat(MNA:MXAD,:)%hme==dat(MNA:MXAD,:)%l).AND.(etr(1,MNA:MXAD,:)>=0).AND.(iacat(MNA:MXAD,:)==1)),d1*etr(1,MNA:MXAD,:),mom,cnt,var)
                WRITE(name(im),'("u2emp | move(MNA:MXAD,:),athme1  ",I2)') g
                weights(im)=0.0_dp !0.0_rp
                im=im+1 
                CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%sexr==g).AND.(dat(MNA:MXAD,:)%rel==0).and.(norelchg(MNA:MXAD,:)==1).AND.(move(MNA:MXAD,:)==1).AND.(dat(MNA:MXAD,:)%hme/=dat(MNA:MXAD,:)%l).AND.(etr(1,MNA:MXAD,:)>=0).AND.(iacat(MNA:MXAD,:)==1)),d1*etr(1,MNA:MXAD,:),mom,cnt,var)
                WRITE(name(im),'("u2emp | move(MNA:MXAD,:),athme0  ",I2)') g
                weights(im)=0.0_dp !0.0_rp
                im=im+1 
            end do 
        
        


            !THESE DON'T WORK. PROBABLY NEED TO CHANGE FULLTIME VARIABLE IN GLOBAL
            
            headloc(ihead)=im ; headstr(ihead)='everyone misc' ; ihead=ihead+1
            headloc(ihead)=im; headstr(ihead)='EMP2EMP by gender/home - SINGLES - FULL TIME (all ages)';ihead=ihead+1
            do g=1,2
                CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%sexr==g).AND.(dat(MNA:MXAD,:)%rel==0).and.(norelchg(MNA:MXAD,:)==1).AND.(move(MNA:MXAD,:)==0).AND.(etr(2,MNA:MXAD,:)>=0).AND.(iacat(MNA:MXAD,:)==1)),d1*etr(2,MNA:MXAD,:),mom,cnt,var)
                WRITE(name(im),'("emp2emp | stay        ",I2)') g
                weights(im)=0.0_dp !0.0_rp
                im=im+1 
                CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%sexr==g).AND.(dat(MNA:MXAD,:)%rel==0).and.(norelchg(MNA:MXAD,:)==1).AND.(move(MNA:MXAD,:)==1).AND.(dat(MNA:MXAD,:)%hme==dat(MNA:MXAD,:)%l).AND.(etr(2,MNA:MXAD,:)>=0).AND.(iacat(MNA:MXAD,:)==1)),d1*etr(2,MNA:MXAD,:),mom,cnt,var)
                WRITE(name(im),'("emp2emp | move(MNA:MXAD,:),athme1",I2)') g
                weights(im)=0.0_dp !0.0_rp
                im=im+1 
                CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%sexr==g).AND.(dat(MNA:MXAD,:)%rel==0).and.(norelchg(MNA:MXAD,:)==1).AND.(move(MNA:MXAD,:)==1).AND.(dat(MNA:MXAD,:)%hme/=dat(MNA:MXAD,:)%l).AND.(etr(2,MNA:MXAD,:)>=0).AND.(iacat(MNA:MXAD,:)==1)),d1*etr(2,MNA:MXAD,:),mom,cnt,var)
                WRITE(name(im),'("emp2emp | move(MNA:MXAD,:),athme0",I2)') g
                weights(im)=0.0_dp !0.0_rp
                im=im+1 
            end do 





            !THESE DON'T WORK. PROBABLY NEED TO CHANGE FULLTIME VARIABLE IN GLOBAL
            
            headloc(ihead)=im ; headstr(ihead)='everyone misc' ; ihead=ihead+1
            headloc(ihead)=im; headstr(ihead)='move(MNA:MXAD,:) rates by gender,at home or not at home and employment status (less than full time or fulltime) ';ihead=ihead+1
            do g=1,2
                CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%sexr==g).AND.(dat(MNA:MXAD,:)%hhr>=0).AND.(dat(MNA:MXAD,:)%hhr<H_FULLTIME).AND.(dat(MNA:MXAD,:)%l/=dat(MNA:MXAD,:)%hme).AND.(move(MNA:MXAD,:)>=0)),d1*move(MNA:MXAD,:),mom,cnt,var)
                WRITE(name(im),'("move(MNA:MXAD,:) | u,athme0 ",I2)') g
                weights(im)=0.0_dp !0.0_rp
                im=im+1
                CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%sexr==g).AND.(dat(MNA:MXAD,:)%hhr>=0).AND.(dat(MNA:MXAD,:)%hhr<H_FULLTIME).AND.(dat(MNA:MXAD,:)%l==dat(MNA:MXAD,:)%hme).AND.(move(MNA:MXAD,:)>=0)),d1*move(MNA:MXAD,:),mom,cnt,var)
                WRITE(name(im),'("move(MNA:MXAD,:) | u,athme1 ",I2)') g
                weights(im)=0.0_dp !0.0_rp
                im=im+1
                CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%sexr==g).AND.(dat(MNA:MXAD,:)%hhr>=0).AND.(dat(MNA:MXAD,:)%hhr>=H_FULLTIME).AND.(dat(MNA:MXAD,:)%l/=dat(MNA:MXAD,:)%hme).AND.(move(MNA:MXAD,:)>=0)),d1*move(MNA:MXAD,:),mom,cnt,var)
                WRITE(name(im),'("move(MNA:MXAD,:) | e,athme0 ",I2)') g
                weights(im)=0.0_dp !0.0_rp
                im=im+1
                CALL condmom(im,((dat(MNA:MXAD,:)%co==co).AND.(dat(MNA:MXAD,:)%sexr==g).AND.(dat(MNA:MXAD,:)%hhr>=0).AND.(dat(MNA:MXAD,:)%hhr>=H_FULLTIME).AND.(dat(MNA:MXAD,:)%l==dat(MNA:MXAD,:)%hme).AND.(move(MNA:MXAD,:)>=0)),d1*move(MNA:MXAD,:),mom,cnt,var)
                WRITE(name(im),'("move(MNA:MXAD,:) | u,athme1 ",I2)') g
                weights(im)=0.0_dp !0.0_rp
                im=im+1
            end do 
        
        

            !THESE DON'T WORK. PROBABLY NEED TO CHANGE FULLTIME VARIABLE IN GLOBAL
            headloc(ihead)=im ; headstr(ihead)='everyone misc' ; ihead=ihead+1
            headloc(ihead)=im; headstr(ihead)='Employment by gender and spouse employment status ';ihead=ihead+1
            do g=1,2
                CALL condmom(im,((dat(MNA:MXA,:)%co==co).AND.(dat(MNA:MXA,:)%sexr==g).AND.(dat(MNA:MXA,:)%rel==1).AND.(dat(MNA:MXA,:)%hhsp>=0).AND.(dat(MNA:MXA,:)%hhsp<H_FULLTIME).AND.(dat(MNA:MXA,:)%hhr>=0).AND.(iacat(MNA:MXA,:)==1)),d1*one(dat(MNA:MXA,:)%hhr>=H_FULLTIME),mom,cnt,var)
                WRITE(name(im),'("emp | sp<fulltime ",I4)') g
                weights(im)=0.0_dp !0.0_rp
                im=im+1
                CALL condmom(im,((dat(MNA:MXA,:)%co==co).AND.(dat(MNA:MXA,:)%sexr==g).AND.(dat(MNA:MXA,:)%rel==1).AND.(dat(MNA:MXA,:)%hhsp>=H_FULLTIME).AND.(dat(MNA:MXA,:)%hhr>=0).AND.(iacat(MNA:MXA,:)==1)),d1*one(dat(MNA:MXA,:)%hhr>=H_FULLTIME),mom,cnt,var)
                WRITE(name(im),'("emp | sp->=fulltime ",I4)') g
                weights(im)=0.0_dp !0.0_rp
                im=im+1
            end do 
            !headloc(ihead)=im; headstr(ihead)='Mar rates by education, ias 25-35';ihead=ihead+1
            CALL condmom(im,((dat(25:35,:)%co==co).AND.(dat(25:35,:)%rel>=0).AND.(dat(25:35,:)%edr==1)),d1*one(dat(25:35,:)%rel==1),mom,cnt,var)
            WRITE(name(im),'("frac married, hs 25:35")')
            weights(im)=0.0_dp !0.0_rp              
            im=im+1
            CALL condmom(im,((dat(25:35,:)%co==co).AND.(dat(25:35,:)%rel>=0).AND.(dat(25:35,:)%edr==2)),d1*one(dat(25:35,:)%rel==1),mom,cnt,var)
            WRITE(name(im),'("frac married, col , 25:35")')
            weights(im)=0.0_dp !0.0_rp
            im=im+1

            headloc(ihead)=im ; headstr(ihead)='everyone misc' ; ihead=ihead+1
            headloc(ihead)=im; headstr(ihead)='wages by gender/rel - FULL TIME';ihead=ihead+1
            do g=1,2
                do j=0,1
                    CALL condmom(im,((dat(MNA:MXA,:)%co==co).AND.(dat(MNA:MXA,:)%sexr==g).AND.(dat(MNA:MXA,:)%rel==j) .AND.(dat(MNA:MXA,:)%hhr>=H_FULLTIME).AND.(dat(MNA:MXA,:)%logwr>=0).AND.(iacat(MNA:MXA,:)==1)),d1*dat(MNA:MXA,:)%logwr,mom,cnt,var)
                    WRITE(name(im),'("logwage/gender/rel ",2I4)') g,j
                    if (j==0) weights(im)=0.0_dp !wwage0 
                    if (j==1) weights(im)=0.0_dp !wwage1 
                    im=im+1
                end do 
            end do 
        
            headloc(ihead)=im; headstr(ihead)='log wage spXr - FULL TIME ';ihead=ihead+1
            CALL condmom(im,((dat(MNA:MXA,:)%co==co).AND.(dat(MNA:MXA,:)%rel==1).AND.(dat(MNA:MXA,:)%hhr>=H_FULLTIME).AND.(dat(MNA:MXA,:)%hhsp>=H_FULLTIME).AND.(dat(MNA:MXA,:)%logwr>=0).AND.(dat(MNA:MXA,:)%logwsp>=0).AND.(iacat(MNA:MXA,:)==1)),d1*dat(MNA:MXA,:)%logwr*dat(MNA:MXA,:)%logwsp,mom,cnt,var)
            WRITE(name(im),'("log-wage-spXr")') 
            weights(im)=0.0_dp !wwage1
            im=im+1


            !headloc(ihead)=im ; headstr(ihead)='everyone misc' ; ihead=ihead+1
            headloc(ihead)=im; headstr(ihead)='Mean wage by gender/ed - FULL TIME (all ages)';ihead=ihead+1
            do g=1,2
                do j=1,2    
                    CALL condmom(im,((dat(MNA:MXA,:)%co==co).AND.(dat(MNA:MXA,:)%sexr==g).AND.(dat(MNA:MXA,:)%edr==j).AND.(dat(MNA:MXA,:)%hhr>=H_FULLTIME).AND.(dat(MNA:MXA,:)%logwr>=0).AND.(iacat(MNA:MXA,:)==1)),d1*dat(MNA:MXA,:)%logwr,mom,cnt,var)
                    WRITE(name(im),'("logw/gender/ed ",2I4)') g,j
                    weights(im)=0.0_dp !wwage0
                    im=im+1
                end do 
            end do 

            headloc(ihead)=im; headstr(ihead)='mean log wage and log wage sq by gender and by hours';ihead=ihead+1
            do g=1,2
                CALL condmom(im,((dat(MNA:MXA,:)%co==co).AND.(dat(MNA:MXA,:)%sexr==g)   .AND.(dat(MNA:MXA,:)%hhr>=H_FULLTIME).AND.(dat(MNA:MXA,:)%logwr>=0).AND.(iacat(MNA:MXA,:)==1)),d1*(dat(MNA:MXA,:)%logwr),mom,cnt,var)
                WRITE(name(im),'("mean-log-wage/full-time/gender ",I4)') g
                weights(im)=0.0_dp !wwage0
                im=im+1
                CALL condmom(im,((dat(MNA:MXA,:)%co==co).AND.(dat(MNA:MXA,:)%sexr==g)   .AND.(dat(MNA:MXA,:)%hhr>=H_FULLTIME).AND.(dat(MNA:MXA,:)%logwr>=0).AND.(iacat(MNA:MXA,:)==1)),d1*(dat(MNA:MXA,:)%logwr)**2,mom,cnt,var)
                WRITE(name(im),'("mean-log-wgsq/full-time/gender ",I4)') g
                weights(im)=0.0_dp !wwage0
                im=im+1
            end do 
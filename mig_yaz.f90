module myaz 
	! placer=0 if being called before getdec_s for sep
	! placer=1 if being called from checknb for within marriage decisions
	! placer=2 if being called from checknb for whether to get married when meet
	! placer=3 if being called before getdec_s for sin
	use params
	use share

	implicit none 

contains										
	subroutine yaz_decs(dd,vcho) !qstar,v0,v_cond_u,vbar,val)
	integer(i4b), intent(in) :: dd(:)  !,qstar
	real(dp), intent(in) :: vcho(ncs)   ! v0,v_cond_u,vbar,val(:,:)
	integer(i4b) :: ia,trueindex,q,x,z,q0,g,j,iepsmove
	ia=dd(1) 
	trueindex=dd(2) 
	q=dd(3) 
	x=dd(4) 
	z=dd(5) 
	q0=dd(6)		
	g=dd(7) 
	j=dd(8) 
    iepsmove=dd(11)
	!kid=x2kid(x)
	!altq=dd(9)
	write(100,'(" ************************************************************************************************************************** ")' ) 	
	if (iter==1)  	write(100,'("in solve/getdec_s iteration 1")' ) 
	if (iter==2)  	write(100,'("in solve/getdec_s iteration 2")' ) 
	if (g==1)  	write(100,'("male: age: ",i4)' ) 	ia
	if (g==2)  	write(100,'("female: age: ",i4)' ) ia
    write(100,*)
	write(100,'(" ************************************************************************************************************************** ")' ) 	
	write(100,'(tr1,"alphakid nokid",tr2,"alphakid kid", tr2,"alphaed noed",tr4,"alphaed ed",tr9,"uhome",tr5,"move cost")')
	write(100,'(6f14.2)') alphakid(g,1),alphakid(g,2),alphaed(g,1),alphaed(g,2),uhome(g),movecost(q0,x,trueindex)
	write(100,'(9(tr10,"uloc"))')
	write(100,'(9f14.2)') uloc(:)
	write(100,'(1x,"moveshock(ieps)")')
	if (g==1)  	write(100,'(1x,f14.2)') moveshock_m(iepsmove)
	if (g==2)  	write(100,'(1x,f14.2)') moveshock_f(iepsmove)
    write(100,*)
    write(100,'(" ************************************************************************************************************************** ")' ) 	
	write(100,*)
    write(100,'(14x,tr2,"qs",tr2,tr2,"ws",tr2,"ls")') 
	write(100,'(6x,"current ",i6,2i4)') q0,q2w(q0),q2l(q0)
	write(100,'(6x,"draw    ",i6,2i4)') q,q2w(q),q2l(q)
	!write(100,'(6x,"optimal ",i6,2i4)') qstar,q2w(qstar),q2l(qstar)
    write(100,*)    
	!write(100,'(tr12,"v0",tr6,"v_cond_u",tr10,"vbar",6(tr4,"v(:,ofloc)") ) ')  
	!write(100,'(3f14.2,6f14.2)') v0,v_cond_u,vbar,val( 1:np1, q2l(q) )
	write(100,*)
    !write(100,'(tr12,"v0",tr6,"v_cond_u",tr10,"vbar",9(tr6,"v(np1,:)") ) ')  
	!write(100,'(12f14.2)') v0,v_cond_u,vbar,val( np1, 1:nl )
	write(100,'(tr12,"vcho")')  
	write(100,'(11f10.2)') vcho(:)
    write(100,'(" ************************************************************************************************************************** ",3/)' ) 	
	end subroutine yaz_decs


	subroutine yaz_checknb(dd,vec,transfers,def)         !     (dd,vec,transfers,cornersol,vsum,pc,def)
	integer(i4b), intent(in) :: dd(:)
	real(dp), intent(in) :: vec(5),transfers(2)
	logical, intent(in) :: def
    logical :: pc(2),criter(3),pc_alt(2),criter_alt(3),haveenoughtotransfer,haveenoughtotransfer_alt
	real(dp) :: dumv(np2,nl),utility(4),wage(4),vdif(2),asum,surplus
	integer(i4b) :: ia,index,trueindex,q,x,z,q0,g,j,i
	integer(i4b) :: a,b,c,qbar(2),iepsmove,fnum
	ia=dd(1) 
	trueindex=dd(2) 
	q=dd(3) 
	x=dd(4) 
	z=dd(5) 
	q0=dd(6)
	g=dd(7) 
	j=dd(8) ! alternative that is being evaluated in this call to checknb
	i=dd(9) ! altq i.e. the q that alternative j corresponds to. if being called by mar market, then this is just q as there is no choice to be made. 
	iepsmove=dd(11)
    if (groups) then 
        index=1
    else 
        index=trueindex
    end if
    !kid=maxval(xx2kid(:,x))
	if (whereami==1) write(200,'(" in solve/getdec_c/checknb ........ j  ")' ) 	
	if (whereami==5) write(200,'(" in solve/marriage market ")' ) 	
    if (whereami==4) write(400,'(" in simulation/getdec_c/checknb ........ j  ")' ) 	
    !if (j==1.or.whereamI==1.or.whereami==2) then 			
	!	write(200,'(" ************************************************************************************************************************** ")' ) 	
	!	write(200,'(tr3,"alpha nokid",tr5,"alpha kid",tr9,"uhome",tr10,"umar",tr5,"move cost",9(tr5,"uloc noed"))')
	!	write(200,'(14f14.2)') alpha(1,1),alpha(1,2),uhome(1),mg(z),movecost(xx2x(1,x),trueindex),uloc(:,1)
	!	write(200,'(14f14.2)') alpha(2,1),alpha(2,2),uhome(2),mg(z),movecost(xx2x(2,x),trueindex),uloc(:,1)
	!	write(200,'(" ************************************************************************************************************************** ")' ) 	
	!end if     
    if (whereamI==1) then 
        fnum=200
    else if (whereamI==5) then 
        fnum=200
    else if (whereamI==4) then 
        fnum=400
    else 
        print*, 'no whereamI so fnum not assigned. being called by somewhere not meant to be!'
        stop
    end if 
    write(fnum,*)
	write(fnum,'(6x,tr5,"age",tr6,"q0",tr7,"q",tr7,"z",tr7,"j",tr4,"altq",tr1,"iepsmve",tr1,"trueind",tr1,"sex",TR3,"X",tr3,"def")' )
	write(fnum,'(6x,8i8,i4,I4,L6)')	ia,q0,q,z,j,i,iepsmove,trueindex,g,X,def    ! altq is just the q that altrnative j corresponds to
	write(fnum,*)
    write(fnum,'(14x         ,tr7,"q",2(tr6,"qs",tr6,"ws",tr6,"ls")     )') 
	if (WhereamI==1) then
		write(fnum,'(6x,"current ",7i8)') q0,qq2q(1,q0),qq2w(1,q0),qq2l(1,q0),qq2q(2,q0),qq2w(2,q0),qq2l(2,q0)
	end if 
	write(fnum,'(6x,"draw    ",7i8)') q,qq2q(1,q),qq2w(1,q),qq2l(1,q),qq2q(2,q),qq2w(2,q),qq2l(2,q)
	if (whereamI==1.or.whereamI==4) then 
        write(fnum,'(6x,"alt     ",7i8)') i,qq2q(1,i),qq2w(1,i),qq2l(1,i),qq2q(2,i),qq2w(2,i),qq2l(2,i)
        ! males and females: outside options
		a=qq2q(1,q) 
		b=xx2x(1,x) 
		c=qq2q(1,q0)
		qbar(1) = decm0_s(iepsmove,a,b,c,ia,index)
		a=qq2q(2,q) 
		b=xx2x(2,x) 
		c=qq2q(2,q0)		
		qbar(2) = decf0_s(iepsmove,a,b,c,ia,index)
		write(fnum,'(3x,"qbar",2i8)') qbar(1),qbar(2)
        write(fnum,'(3x,"outside opt",7i8)') q2qq(qbar(1),qbar(2)),qbar(1),q2w(qbar(1)),q2l(qbar(1)),qbar(2),q2w(qbar(2)),q2l(qbar(2))
        !if (def) tmp(1) = (  ( vsum(1) - vec(1) )**0.5_dp ) * (  ( vsum(2) - vec(2) )**0.5_dp ) 
		!if (def) tmp(2) = vsum(1) - vec(1) + vsum(2)- vec(2) 
		! males: utility,wage,value func
		a=qq2q(1,i) 
		b=xx2x(1,x) 
		utility(1) = utils(1,a,b,trueindex) !utilm_s()
		utility(3) = utilc(1,i,x,trueindex) !utilm_c(i,x)
		wage(1)=ws(1,a,b,trueindex)    !wm_s(a,b)
        wage(3)=wc(1,i,x,trueindex)
		! females: utility,wage,value func
		a=qq2q(2,i) 
		b=xx2x(2,x) 
		utility(2) = utils(2,a,b,trueindex)
		utility(4) = utilc(2,i,x,trueindex)
		wage(2)=ws(2,a,b,trueindex)    !wf_s(a,b)
        wage(4)=wc(2,i,x,trueindex)
        if (wage(3)+wage(4) - vec(5) > 0.0_dp ) then
            print*, 'sum of wage3 and wage4 not equal to vec5!'
            stop
        end if
        !vdif = vec(3:4) + mg( dd(5) ) - vec(1:2) 	! dd(5) is z
		vdif = vec(3:4) - vec(1:2) 	
	    surplus=vec(3)-vec(1)+vec(4)-vec(2)+vec(5)
        pc(1:2)	= ( vdif + eps >= 0.0_dp )	!pc(1:2)    = ( vec(3:4) - vec(1:2) >= 0.0_dp )						        
        pc_alt(1:2)=( vdif >= 0.0_dp )	
        asum = sum(  one(.not. pc)  *   abs( vdif )   )
		criter(1) = ( vec(5) + eps >= vec(1) - vec(3) ) 
		criter(2) = ( vec(5) + eps >= vec(2) - vec(4) ) 
		criter(3) = ( vec(5) + eps >= vec(1) - vec(3) + vec(2) - vec(4) )     		
		criter_alt(1) = ( vec(5) >= vec(1) - vec(3) ) 
		criter_alt(2) = ( vec(5) >= vec(2) - vec(4) ) 
		criter_alt(3) = ( vec(5) >= vec(1) - vec(3) + vec(2) - vec(4) )
        haveenoughtotransfer=(  vec(5) + eps - asum  >= 0.0_dp  )
        haveenoughtotransfer_alt=(  vec(5) - asum  >= 0.0_dp  )
        if (nonneg) then 
            if ( maxval(one(pc)-one(pc_alt)) >0) then 
                write(fnum,'("pc with eps no match to pc without eps")')
                write(*,'("pc with eps no match to pc without eps")')
                stop 
            end if
            if ( maxval(one(criter)-one(criter_alt)) >0) then 
                write(fnum,'("criter with eps no match to criter without eps")')
                write(*,'("criter with eps no match to criter without eps")')
                stop 
            end if
            if ( one(haveenoughtotransfer).ne.one(haveenoughtotransfer_alt) ) then 
                write(fnum,'("vec(5)-asum with eps is no match to without eps ")')
                write(*,'("vec(5)-asum with eps is no match to without eps ")')
                stop         
            end if 
        
            !write(fnum,'("these utilities are not indexed now so be careful")')
		    !write(fnum,'(6x,"alt spec values:") ') 
            write(fnum,'("PC before any transfers (pc):",2L6)') pc
            write(fnum,'("Abs Val of Total Transfers Needed (asum):",F10.2)') asum
            !ahu030822 write(fnum,'("Total Wages (wc):",2L6)') vec(5)
            write(fnum,'("Total Wages (wc):",F10.2)') vec(5) !ahu030822 replacing 2L6 with F10.2
            write(fnum,'("haveenoughtotransfer? ie vec5+eps-asum>=0.0dp?:",L6)') haveenoughtotransfer !ahu030822
            write(fnum,'("haveenoughtotransfer_alt? ie vec5-asum>=0.0dp?:",L6)') haveenoughtotransfer_alt !ahu030822   
            write(fnum,'("interior? ie vec5 vs. Vdif1-Vdif2:",2F14.4)') vec(5),Vdif(1)-Vdif(2) !ahumarch1122   
            if ( (.not.pc(1)).or.(.not.pc(2)) ) then 
                write(fnum,'("One or both PC do not hold")')
                if ( haveenoughtotransfer  ) then !haveenoughtotransfer=(  vec(5) + eps - asum  >= 0.0_dp  )
                    write(fnum,'("BUT they HAVE enough wages to make transfers!")')
                else 
                    write(fnum,'("AND they DO NOT HAVE enough wages to make transfers!")')
                end if 
            else
                write(fnum,'("Both PCs hold")')
            end if 
            if (surplus+eps>=0.0_dp.and. (.not.haveenoughtotransfer) ) then 
                write(fnum,'("Positive surplus but they are not able to make transfers!")' )
            end if 
            write(fnum,*)
            !write(fnum,'(6x,tr10,"asum",tr4,"w+eps-asum>=0.0dp",3(tr2,"crit")  )' ) 
		    !write(fnum,'(6x,f14.2,11x,L6,3L6)') asum,haveenoughtotransfer,criter(1:3)     !def  = ( vec(5) + eps - asum  >= 0.0_dp )
            write(fnum,'(6x,2(tr9,"trans"),tr3,"def" )')         
            write(fnum,'(6x,2F14.2,L6)') transfers,def                 
            if (     ( minval(transfers) + eps2) >= 0.0_dp  .and. ( minval(transfers) + epstest ) >= 0.0_dp  ) then 
                write(fnum,*) "ohere minval"
            end if 
            if (      maxval(transfers) <=  ( vec(5)+eps2 ) .and.  maxval(transfers) <= ( vec(5)+epstest ) ) then 
                write(fnum,*) "ohere maxval"
            end if 
        end if !NONNEG
        write(fnum,'(6x,2(tr10,"vbar"),2(tr12,"vc"),tr5,"wcsum",2(tr4,"pc"),2(tr8,"us"),2(tr8,"uc"),2(tr8,"ws"),2(tr8,"wc"))' )    
		write(fnum,'(6x,4f14.2,F10.2,2L6,8F10.2)') vec,pc,utility(1:4),wage(1:4)
        write(fnum,*)
        write(fnum,'(6x,tr7,"surplus",tr2,"sur+eps>=0.0dp",tr2,"sur>=0.0dp",tr3,"def")') 
		write(fnum,'(6x,f14.2,10x,L6,6x,2L6)') surplus,(surplus+eps>=0.0_dp),(surplus>=0.0_dp),def
        if (def) then 
            write(fnum,'(6x,tr10,"vec1",tr1,"vec1+0.5*surpls",tr10,"vec2",tr1,"vec2+0.5*surpls")')
		    write(fnum,'(6x,f14.2,2x,f14.2,f14.2,2x,f14.2)') vec(1),vec(1)+0.5_dp*surplus,vec(2),vec(2)+0.5_dp*surplus
        end if 
    end if !whereamI=1 or 4
    
    if (whereamI==5) then 
        !NOTE that when this is being called from the marriage market, the index that corresponds to altq (i.e. i) is just set to q since there is no altq in marriage market
		! males: utility,wage,value func
		a=qq2q(1,i) 
		b=xx2x(1,x) 
		utility(1) = utils(1,a,b,trueindex) !utilm_s()
		utility(3) = utilc(1,i,x,trueindex) !utilm_c(i,x)
		wage(1)=ws(1,a,b,trueindex)    !wm_s(a,b)
        wage(3)=wc(1,i,x,trueindex)
        ! females: utility,wage,value func
		a=qq2q(2,i) 
		b=xx2x(2,x) 
		utility(2) = utils(2,a,b,trueindex)
		utility(4) = utilc(2,i,x,trueindex)
		wage(2)=ws(2,a,b,trueindex)    !wf_s(a,b)
        wage(4)=wc(2,i,x,trueindex)        
        !vdif = vec(3:4) + mg( dd(5) ) - vec(1:2) 	! dd(5) is z
		vdif = vec(3:4) - vec(1:2) 	
	    surplus=vec(3)-vec(1)+vec(4)-vec(2)+vec(5)
        pc(1:2)	= ( vdif + eps >= 0.0_dp )	!pc(1:2)    = ( vec(3:4) - vec(1:2) >= 0.0_dp )		
        write(fnum,'("Here is q,x",2I4)') q,x
        write(fnum,'(6x,2(tr10,"vbar"),2(tr12,"vc"),tr5,"wcsum",2(tr4,"pc"),2(tr8,"us"),2(tr8,"uc"),2(tr8,"ws"),2(tr8,"wc") )' )    
		write(fnum,'(6x,4f14.2,F10.2,2L6,8F10.2)') vec,pc,utility(1:4),wage(1:4)        
        write(fnum,*)
		write(fnum,'(6x,tr7,"surplus",tr2,"sur+eps>=0.0dp",tr2,"sur>=0.0dp",tr7,"surplus")') 
		write(fnum,'(6x,f14.2,10x,L6,6x,L6,F14.2)') surplus,(surplus+eps>=0.0_dp),(surplus>=0.0_dp),surplus
        write(fnum,'(6x,2(tr9,"trans"),tr3,"def" )')         
        write(fnum,'(6x,2F14.2,L6)') transfers,def
        if (     ( minval(transfers) + eps2) >= 0.0_dp  .and. ( minval(transfers) + epstest ) >= 0.0_dp  ) then 
            write(fnum,*) "ohere minval"
        end if 
        if (      maxval(transfers) <=  ( vec(5)+eps2 ) .and.  maxval(transfers) <= ( vec(5)+epstest ) ) then 
            write(fnum,*) "ohere maxval"
        end if 
        !ahu 041118: the below stuff was trying to figure out the eps2 problem. see notes for feb, march and april 2018 more details. 
                        !if (  minval(transfers) + epstest4 >= 0.0_dp .and. maxval(transfers) <= vec(5)+epstest4 ) then
                            !if (  minval(transfers) + epstest5 < 0.0_dp .or. maxval(transfers) > vec(5)+epstest5 ) then
                                    !write(fnum,*) "ohere45"                    
                            !end if 
                        !end if 
    end if !whereamI=5
	!move(:) = fnmove(kid) 
	!move( qq2l(1,q0) )= 0.0_dp
	!if (groups) then 
	!	a=myco
	!	b=mytyp
	!	c=myhome		
	!else 
	!	call index2cotyphome(k,a,b,c)	! don't need this if groups, since myco,mytyp,myhome gives co,typ,home already			
	!end if 
	!write(200,'(6x,"outside options:") ') 
	!ya=xx2x(1,x)
	!dumv(:,:) = reshape( vm0_s(:,ya) , (/ np2, nl /) )				! turn q index into its components (w,l)
	!ya = qq2w(1,q0) ; yb=qq2l(1,q0) ; yc=qq2l(1,q)				! w0,l0,l
	!write(200,'(6x,"best q,w,l    ",6i8)') qbar(1),q2w(qbar(1)),q2l(qbar(1))
	!write(200,'(6x,tr10,"vs0",tr10,"vsu",6(tr6,"vs(:,l)") )' )
	!write(200,'(8f14.2)') dumv(ya,yb),maxval(dumv(np1,1:nl) + move ),dumv(1:np1,yc) + move(yc)
	!write(200,*)
	!ya=xx2x(2,x)
	!dumv(:,:) = reshape( vf0_s(:,ya) , (/ np2, nl /) )				! turn q index into its components (w,l)
	!ya = qq2w(2,q0) ; yb=qq2l(2,q0) ; yc=qq2l(2,q)				! w0,l0,l
	!write(200,'(tr10,"vfs0",tr10,"vfsu",6(tr6,"vfs(:,l)") )' )
	!write(200,'(8f14.2)') dumv(ya,yb),maxval(dumv(np1,1:nl) + move ),dumv(1:np1,yc) + move(yc)
	!do ya=1,2	! gender
	!	yb=qq2q(ya,altq) 
	!	if (qq2w(ya,altq) /= q2w(yb)  )			ier(2) = 1 ! qq and q are not equal in terms of their w 
	!	if (qq2l(ya,altq) /= q2l(yb)  )			ier(3) = 1 ! qq and q are not equal in terms of their l 
	!end do 	
	!if ( decsingle(1)<0 .or. decsingle(4)<0 )		ier(8) = 1				! single dec rule qsbest is < 0 
	!if ( abs( w_c(altq,x) - vec(5) ) > 0.0_dp )		ier(9) = 1				! in write_checknb: wc is not equal to vec5 ",ii 
	!if ( def == -2)					ier(10) = 1
	!if ( def .and. ia==mxa .and. abs(tmp(2)) > 2.0_dp*eps .and. abs(mg(z)) < eps ) 
	!	write(*,*) tmp(2),mg(z)
	!	write(*,'(2l6)') ( abs(tmp(2)) > eps ), ( abs(mg(z)) < eps  ) 
	!end if 
	end subroutine yaz_checknb

	subroutine yaz_decision(dd,vmax)
	integer(i4b), intent(in) :: dd(:)
    real(8), intent(in) :: vmax(2)
	integer(i4b) :: q0,jmax,qmax,relmax,fnum
    if (whereamI==1) then           !in solve/getdec_c/checknb
        fnum=200
    else if (whereamI==5) then      !in solve/marriage market
        fnum=200
    else if (whereamI==4) then      !in simulation/getdec_c/checknb 
        fnum=400
    else 
        print*, 'no whereamI so fnum not assigned. being called by somewhere not meant to be!'
        stop
    end if
    write(fnum,*)
    write(fnum,'(1x,"****************************************************************************************************************** ")' ) 	
	if (whereami==1.or.whereami==4) then
		q0=dd(6)
		jmax=dd(8)
		qmax=dd(9)	
		relmax=dd(10)
		if (relmax>0) then 
			if ( qq2l(1,qmax) == qq2l(1,q0) ) then 
				write(fnum,'(1x,"decision: stay married and stay put. best j,q is: ",2i8,2F14.2,4i4 )') jmax,qmax !,vmax(1:2),q0,dd(3:5)
			else 
				write(fnum,'(1x,"decision: stay married and move. best j,q is: ",2i8,2F14.2,4i4 )') jmax,qmax !,vmax(1:2),q0,dd(3:5)
			end if 
		else if (relmax==0) then 
			write(fnum,'(1x,"decision: get divorced",4i4)') !q0,dd(3:5)
		end if 
	else if (whereami==5) then
		relmax=dd(10)
		if (relmax==1) then 
			write(fnum,'(1x,"decision: get married")') 
		else 
			write(fnum,'(1x,"decision: stay single and continue searching")') 
		end if 
	end if 
	write(fnum,'(1x,"****************************************************************************************************************** ",3/)' ) 	
	write(fnum,*)
    end subroutine yaz_decision

	subroutine yaz_getmom(dat,ndat)
	integer(i4b), intent(in) :: ndat						! number of observations in dat array    
	type(statevar), dimension(mnad:mxa,ndat), intent(in) :: dat ! data set. first entry is ia index, second observation number
	integer :: ia,j
	write(12,*)
	!if (ndat==ndata) then
!		write(12,'("actual data")')			
	!else 
	if (ndat==ndat) then 		
		!		write(12,'("simulated data")')	
				write(12,*)
				do j=1,2
					write(12,'(tr6,"id",tr1,"age",tr1,"sexr",tr1,"exp",tr1,"hhr",tr4,"logwr",tr1,"kid",tr1,"edr",&
					& tr1,"rel",&
					& tr1,"loc",tr1,"mxa",tr1,"mis",tr1,"hme",tr1,"hhr",tr1,"hsp",tr4,"logwr",tr3,"logwsp",tr1,"lsp" )')	
					write(12,*)
					do ia=mnad,mxa
						write(12,'(i8,4i4,f9.2,9i4,2f9.2,i4)') j,&
						& ia,dat(ia,j)%sexr,dat(ia,j)%expr,dat(ia,j)%hhr,dat(ia,j)%logwr,dat(ia,j)%kidr,dat(ia,j)%edr,dat(ia,j)%rel,&
						& dat(ia,j)%l,dat(ia,j)%endage,dat(ia,j)%nomiss,dat(ia,j)%hme,&
						& dat(ia,j)%hhr,dat(ia,j)%hhsp,dat(ia,j)%logwr,dat(ia,j)%logwsp,dat(ia,j)%lsp
					write(12,*)
					write(12,*)
					end do
				end do
	end if
	if (ndat==nsim) then 		
!		write(12,'("simulated data")')	
		write(12,*)
		do j=1,1000
			write(12,'(tr6,"id",tr1,"age",tr1,"sexr",tr1,"exp",tr1,"hhr",tr4,"logwr",tr1,"kid",tr1,"edr",&
			& tr1,"rel",&
			& tr1,"loc",tr1,"mxa",tr1,"mis",tr1,"hme",tr1,"hhr",tr1,"hsp",tr4,"logwr",tr3,"logwsp",tr1,"lsp" )')	
			write(12,*)
			do ia=mnad,mxa
				write(12,'(i8,4i4,f9.2,9i4,2f9.2,i4)') j,&
				& ia,dat(ia,j)%sexr,dat(ia,j)%expr,dat(ia,j)%hhr,dat(ia,j)%logwr,dat(ia,j)%kidr,dat(ia,j)%edr,dat(ia,j)%rel,&
				& dat(ia,j)%l,dat(ia,j)%endage,dat(ia,j)%nomiss,dat(ia,j)%hme,&
				& dat(ia,j)%hhr,dat(ia,j)%hhsp,dat(ia,j)%logwr,dat(ia,j)%logwsp,dat(ia,j)%lsp
			write(12,*)
			write(12,*)
			end do
		end do
	end if 
	end subroutine yaz_getmom

	subroutine yaz_sim(gender,rel,q,x)
	integer(i4b), intent(in) :: gender,rel,q,x
		if (rel==0) then
			write(400,'(2x,tr3,"g",tr2,"qs",tr2,"ws",tr2,"ls",tr4,"xs",tr3,"eds",tr2,"exps",tr2,"kids")' ) 
			write(400,'(2x,4I4,4I6  )') gender,q,q2w(q),q2l(q),x,x2e(x),x2r(x),x2kid(x)      
		else if (rel==1) then 
			write(400,'(2x,tr6,"qc",tr2,"qm",tr2,"wm",tr2,"lm", tr6,"xc",tr4,"xm",tr3,"edm",tr2,"expm",tr2,"kidm"    )') 
			write(400,'(2x,I8,3I4,I8,4I6)') q,qq2q(1,q),qq2w(1,q),qq2l(1,q),x,xx2x(1,x),xx2e(1,x),xx2r(1,x),xx2kid(1,x)
			write(400,'(2x,tr6,"qc",tr2,"qf",tr2,"wf",tr2,"lf", tr6,"xc",tr4,"xf",tr3,"edf",tr2,"expf",tr2,"kidf"    )') 
			write(400,'(2x,I8,3I4,I8,4I6)') q,qq2q(2,q),qq2w(2,q),qq2l(2,q),x,xx2x(2,x),xx2e(2,x),xx2r(2,x),xx2kid(2,x)
        end if 
	end subroutine yaz_sim
	subroutine yaz_simmatch(meet,qmatch,xmatch,z)
	logical, intent(in) :: meet
	integer(i4b), intent(in) :: qmatch,xmatch,z
		if (meet) then
			write(400,'(2x,tr7,"z",tr2,"qs",tr2,"ws",tr2,"ls",tr2,"xs",tr3,"eds",tr2,"exps",tr2,"kids")' )
			write(400,'(2x,I8,3I4,I4,3I6)') z,qmatch,q2w(qmatch),q2l(qmatch),xmatch,x2e(xmatch),x2r(xmatch),x2kid(xmatch)  
		else 
			write(400,'("Did Not Meet Anyone")')            
		end if 
	end subroutine yaz_simmatch

    subroutine yaz_simdecmar(rel)
    integer(i4b), intent(in) :: rel
        if (rel==0) then 
            write(400,'("Relnext: Single")') 
        else if (rel==1) then 
            write(400,'("Relnext: Married")')
        else 
            write(400,'("Relnext: NOTHING! SOMETHING WRONG!")')
            print*, 'Rel is not 0 nor 1, something wrong!'
            stop
        end if 
    end subroutine yaz_simdecmar
    
	subroutine yaz_simpath(ia,nn,mm,r,dat)
	integer(i4b), intent(in) :: ia,nn,mm,r
	type(statevar) :: dat
	write(500,'(tr6,"nn",tr2,"mm",tr7,"r",tr2,"ia",tr2,"co",tr2,"sx",tr2,"hm",tr2,"ed",tr2,"ea",tr2,"rl",tr2,"kd",tr3,"l",tr2,"hr",tr2,"hp",tr1,"len",tr2,"ms",tr12,"wr",tr11,"wsp",tr1,"exp",tr2,"wr",tr1,"wsp")') 
	write(500,'(I8,I4,I8,13I4,2F14.2,I4,2F14.2)') nn,mm,r,ia,dat%co,dat%sexr,dat%hme,dat%edr,dat%endage,dat%rel,dat%kidr,dat%l,dat%hhr,dat%hhsp,dat%rellen,& 
    & dat%nomiss,dat%logwr,dat%logwsp,dat%expr,dat%wr,dat%wsp    !,dat%nn,dat%mm,dat%r
	end subroutine 
					
	
	subroutine yaz0
	integer(i4b) :: n,i,j,g,q0,x0,q,x
	real(dp) :: lb,ub,cdf(2),expv,ppcqx_sum_overqx
	if (skriv) then 
		print*, "here i am in skriv!" 
		write(50,*) 
		write(50,*) 
		do g=1,2
			if (g==1) then 
				write(50,'(" ---------------- male   wage grid ---------------- ")') 
			else 
				write(50,'(" ---------------- female wage grid ---------------- ")') 
			end if 
			write(50,'(tr3,"num",tr10,"wage",tr8,"weight",tr4,"sum=1?")')
			write(50,*)
			do i=1,np
				write(50,'(i6,2f14.3,f10.2)') i,wg(i,g),wgt(i),sum(wgt)
				write(50,*)
			end do
		end do 
		if (abs( sum(wgt) - 1.0_dp ) > 1.0d-6  ) then 
			print*, "sum of weight is not 1! " 
			stop 
		end if 
		write(50,'(2/," ---------------- mar grid ---------------- ")') 
		if (skriv) print*, "mar grid mg(:,1) ", mg(:,1)
		if (skriv) print*, "mar grid mg(:,2) ", mg(:,2)
		if (skriv) print*, "mar grid mg(:,3) ", mg(:,3)
		if (skriv) print*, "mar grid mg(:,4) ", mg(:,4)
		do i=1,nz
			write(50,'(1x,tr6,"z",tr9,"mg(z,.)",tr8,"mgt(z)",tr6,"sum=1.0?")')
			write(50,'(i8,3f14.3)') i,mg(i,1),mgt(i),sum(mgt)
			write(50,'(i8,3f14.3)') i,mg(i,2),mgt(i),sum(mgt)
			write(50,'(i8,3f14.3)') i,mg(i,3),mgt(i),sum(mgt)
			write(50,'(i8,3f14.3)') i,mg(i,4),mgt(i),sum(mgt)
        end do
		write(50,'(3/," ---------------- moveshock grid ---------------- ")') 
		!if (skriv) print*, "moveshock grid ", moveshock_m(:)
        !if (skriv) print*, "moveshock grid ", moveshock_f(:)
        do i=1,nepsmove
			write(50,'(1x,tr6,"i",tr5,"moveshock(i)",tr7,"ppso(i)",tr6,"sum=1.0?")')
			write(50,'(i8,3f14.3)') i,moveshock_m(i),ppso(i),sum(ppso)
            write(50,'(i8,3f14.3)') i,moveshock_f(i),ppso(i),sum(ppso)
        end do
		write(50,'(4/," ----------------  kidshock grid ---------------- ")') 
		if (skriv) print*, " kidshock grid ", kidshock(:)
		do i=1,nepskid
			write(50,'(1x,tr6,"i",tr5," kidshock(i)",tr7,"ppsk(i)",tr6,"sum=1.0?")')
			write(50,'(i8,f14.3)') i,kidshock(i)
		end do
        write(50,*) 
		write(50,*) "sum(wgt) - free" 
		write(50,*)  sum(wgt)
		write(50,*) 
		write(50,*) "1.0,1.0_sp,1.0_dp - free "
		write(50,*) 1.0,1.0_sp,1.0_dp
		lb=-3.0_dp ; ub=3.0_dp
		call checkgauss(np,lb,ub,cdf,expv)
		write(50,'(2/," ---------------- checkgauss(np,lb,ub,cdf,expv) ---------------- ")') 
		write(50,*) 
		write(50, '(1x,tr11,"lb",tr12,"ub",tr8,"cdf(1)",tr8,"cdf(2)",tr8,"expval" )' ) 
		write(50,'(5f14.3)') lb,ub,expv,cdf(1),cdf(2)
		write(50,'(3/," ---------------- fnprof(w,ed,g) ---------------- ")') 
		write(50,'(" raw parameters: ")') 
		write(50, '(1x,"wrking (w<=np)",5x,2(tr10,"m/hs"),2(tr9,"m/col"),2(tr10,"f/hs"),2(tr9,"f/col") )') 
		write(50,'(19x,8f14.3 )') psio(1:8)	
		write(50, '(1x,"not wrking (w=np1)",tr10,"m/hs",14x,tr9,"m/col",14x,tr10,"f/hs",14x,tr9,"f/col",14x )') 
		write(50,'(19x,f14.3,14x,f14.3,14x,f14.3,14x,f14.3,14x )') psio(9:12)	
		write(50,'(/," probabilities: ")') 
		do n=1,neduc		!ed
			do g=1,2	!sex
				write(50, '(1x,tr4,"sex",tr6,"ed",4x,"wrking (w<=np)",4x,tr5,"pr(offer)",tr4,"pr(layoff)",tr3,"pr(nothing)" )' ) 
				write(50,'(2i8,22x,3(f14.3) )') g,n,fnprof(np,n,g)	
				write(50, '(20x,"not wrking (w=np1)",2x,tr5,"pr(offer)",tr4,"pr(layoff)",tr3,"pr(nothing)" )' ) 
				write(50,'(38x,3(f14.3) )') fnprof(np1,n,g)	
				write(50,*) 
			end do 
		end do 
		write(50,'(/," ---------------- fnprloc parameters ---------------- ")') 
		!write(50,'(/,"iter psil(1) psil(2) psil(3) ",i6,f14.3)') iter,psil(:)
		write(50, '(/,6x,tr12,"ne",tr12,"ma",tr11,"enc",tr11,"wnc",tr12,"sa",tr11,"esc",tr11,"wsc",tr9,"mount",tr11,"pac" )' ) 
		write(50,'(6x,9f14.3)') popsize(1:nl)
		write(50,'(2/,4x,"ne",9(f14.3))') fnprloc(1)		
		write(50,'(6x,9(2x,i8,4x))') distance(:,1)
		write(50,'(/,4x,"ma",9(f14.3) )') fnprloc(2)
		write(50,'(6x,9(3x,i8,3x))') distance(:,2)
		write(50,'(3/,3x,"enc",9(f14.3) )') fnprloc(3)
		write(50,'(6x,9(3x,i8,3x))') distance(:,3)
		write(50,'(/,3x,"wnc",9(f14.3) )') fnprloc(4)
		write(50,'(6x,9(3x,i8,3x))') distance(:,4)
		write(50,'(/,4x,"sa",9(f14.3) )') fnprloc(5)
		write(50,'(6x,9(3x,i8,3x))') distance(:,5)
		write(50,'(/,3x,"esc",9(f14.3) )') fnprloc(6)
		write(50,'(6x,9(3x,i8,3x))') distance(:,6)
		write(50,'(/,3x,"wsc",9(f14.3) )') fnprloc(7)
		write(50,'(6x,9(3x,i8,3x))') distance(:,7)
		write(50,'(/,1x,"mount",9(f14.3) )') fnprloc(8)
		write(50,'(6x,9(3x,i8,3x))') distance(:,8)
		write(50,'(/,3x,"pac",9(f14.3) )') fnprloc(9)
		write(50,'(6x,9(3x,i8,3x))') distance(:,9)
		write(50,'(3/," ---------------- fnprhc(exp,w) ---------------- ")') 
		write(50, '(1x,tr7,"w",tr4,"exp",4(tr8,"fnprhc") )' ) 
        do n=1,np1
        do i=1,nexp
				if (n==1.or.n==np1) then
					write(50,'(2i8,4f14.3)') n,i,fnprhc(i,n)
					write(50,*) 
				end if 
			end do 
		end do 
		i=-1 ; n=-1 ; j=-1

		write(50,'(2/," ---------------- wage profiles given typ=1,loc=nloc ---------------- ")') 
		write(50, '(1x,tr3,"ed",tr3,"exp",2(tr11,"eps",tr10,"wage") )')
		do j=1,neduc
			do i=1,nexp
				do n=1,np
					write(50,'(2i6,4f14.3)') j,i,wg(n,1),fnwge(1,1,nl,wg(n,1),j,i),wg(n,2),fnwge(2,1,nl,wg(n,2),j,i)				!fnwge(dg,dtyp,dl,dw,de,dr)						
				end do 
			end do 
		end do 
		!write(50,'(" ---------------- ppsq(q,q0,x0,g) - males and females ---------------- ")') 
		!ppmq_check=0.0_dp
		!ppfq_check=0.0_dp
		!do n=1,nx
		!	do i=1,nq
		!		do j=1,nq				
		!			if ( pc0(i) ) then 	
		!				x0=xx2x(1,n)
		!				q0=qq2q(1,i)
		!				q=qq2q(1,j)
		!				ppmq_check(q,q0,x0) = ppmq_check(q,q0,x0) + wgt( qq2w(2,i) ) * ppcq(j,i,n)
		!				x0=xx2x(2,n)
		!				q0=qq2q(2,i)
		!				q=qq2q(2,j)
		!				ppfq_check(q,q0,x0) = ppfq_check(q,q0,x0) + wgt( qq2w(1,i) ) * ppcq(j,i,n)
		!			end if 
		!		end do 
		!	end do 
		!end do 
		do x0=1,nxs			
			do q0=1,nqs		
				!do q=1,nqs	
				!	write(50, '(1x,tr5,"x0",tr5,"ed0",tr4,"exp0")' ) 
				!	write(50,'(3i8)')	x0,x2e(x0),x2r(x0)	
				!	write(50, '(1x,tr5,"q0",tr6,"w0",tr6,"l0")' )
				!	write(50,'(3i8)')	q0,q2w(q0),q2l(q0)
				!	write(50, '(1x,tr6,"q" ,tr7,"w" ,tr7,"l",2(6x,"ppsq(q,q0,x0,g)"),2(tr6,"sum=1.0?") )') 
				!	write(50,'(3i8,2(7x,f14.3),2f14.3)') q,q2w(q),q2l(q),ppsq(q,q0,x0,1),ppsq(q,q0,x0,2),sum(ppsq(:,q0,x0,1)),sum(ppsq(:,q0,x0,2))
				!	write(50,'(3i8,2(7x,f14.3),2f14.3)') q,q2w(q),q2l(q),ppmq_check(q,q0,x0),ppfq_check(q,q0,x0),sum(ppmq_check(:,q0,x0)),sum(ppfq_check(:,q0,x0))
				!	write(50,*)
				!	write(50,*)
				!end do 
			   ! do g=1,2
				!	if ( q2w(q0)<=np1) then  
				!		if ( abs(  sum(ppsq(:,q0,x0,g)) - 1.0_dp) >  1.0d-6) then  ; 
				!			print*, "error in check: ppsq for w=1:np1 does not add up to 1 " , q0,x0,g,sum(ppsq(:,q0,x0,g)) 
				!			stop 
				!		end if 
				!	else if ( q2w(q0)>np1) then					
				!		if ( sum(ppsq(:,q0,x0,g)) > 1.0d-6 ) then	! q2w(q0)=np2 is not a state variable that we take in to account, so there's no proability attached to it 
				!			print*, "error in check: ppsq for w=np2 does not add up to 0 " , q0,x0,g,sum(ppsq(:,q0,x0,g)) 
				!			stop 
				!		end if 
				!	end if 
				!end do 
			end do 
		end do 
		write(50,'(2/," ---------------- ppcq(q,q0,x0) ---------------- ")') 
		do x0=1,nx		
			do q0=1,nq		
				!do q=1,nx	
				!	write(50,'(1x,tr5,"x0",2(tr5,"xs0",tr4,"eds0",tr3,"exps0")     )') 
				!	write(50,'(7i8)') x0,xx2x(1,x0),xx2e(1,x0),xx2r(1,x0),xx2x(2,x0),xx2e(2,x0),xx2r(2,x0)
				!	write(50,*)
				!	write(50,'(1x,tr5,"q0",2(tr5,"qs0",tr5,"ws0",tr5,"ls0")     )') 
				!	write(50,'(7i8)') q0,qq2q(1,q0),qq2w(1,q0),qq2l(1,q0),qq2q(2,q0),qq2w(2,q0),qq2l(2,q0)
				!	write(50,*)
				!	write(50,'(1x,tr6,"q",2(tr6,"qs",tr6,"ws",tr6,"ls"),6x,"ppcq(q,q0,x0)",tr6,"sum=1.0?" )') 
				!	write(50,'(7i8,5x,2f14.3)') q,qq2q(1,q),qq2w(1,q),qq2l(1,q),qq2q(2,q),qq2w(2,q),qq2l(2,q),ppcq(q,q0,x0),sum(ppcq(:,q0,x0)) 
				!	write(50,*)
				!	write(50,*)
				!end do 
				if ( qq2l(1,q0) /= qq2l(2,q0) ) then 
					print*, "error in check: locs not equal " 
					stop 
				end if 			 
				if (  qq2w(1,q0)<=np1 .and.  qq2w(2,q0)<=np1 ) then 
					if ( abs(  sum(ppcq(:,q0,x0)) - 1.0_dp  ) >  1.0d-6) then  
						print*, "error in check: ppcq for w=1:np1 does not add up to 1 " , q0,x0 , sum(ppcq(:,q0,x0)) 
						stop 
					end if 
				end if 
				if ( qq2w(1,q0)>np1 .or.  qq2w(2,q0)>np1 ) then 
					if ( sum(ppcq(:,q0,x0)) > 1.0d-6 ) then 
						print*, "error in check: ppcq for w=np2 does not add up to 0 " , q0,x0,sum(ppcq(:,q0,x0)) 
						stop 
					end if 
				end if
			end do 
		end do 

		write(50,'(2/," ---------------- ppsx(x,q0,x0) ---------------- ")') 
		do x0=1,nxs				
			do q0=1,nqs		
				!do x=1,nxs	
				!	write(50, '(1x,tr5,"x0",tr5,"ed0",tr4,"exp0")' ) 
				!	write(50,'(3i8)')	x0,x2e(x0),x2r(x0)	
				!	write(50, '(1x,tr5,"q0",tr6,"w0",tr6,"l0")' )
				!	write(50,'(3i8)')	q0,q2w(q0),q2l(q0)
				!	write(50, '(1x,tr6,"x" ,tr6,"ed" ,tr5,"exp" ,6x,tr1,"ppsx(x,q0,x0)",tr6,"sum=1.0?")') 
				!	write(50,'(3i8,6x,2f14.3)') x,x2e(x),x2r(x),ppsx( x, q0, x0 ),sum(ppsx(:,q0,x0)) 
				!	write(50,*)
				!	write(50,*)
				!end do 
				if (  q2w(q0)<=np1 ) then 
					if ( abs(  sum(ppsx(:,q0,x0)) - 1.0_dp  )   >   1.0d-6) then  
						print*, "error check_all: ppsx for w=1:np1 does not add up to 1 ", q0,x0,sum(ppsx(:,q0,x0)) 
						stop 
					end if 
				else if (  q2w(q0)>np1 ) then 
					if ( sum(ppsx(:,q0,x0)) > 1.0d-6 ) then 
						print*, "error in check_all: ppsx for w=np2 does not add up to 0 " , q0,x0,sum(ppsx(:,q0,x0)) 
						stop 
					end if 
				end if 
			end do 
		end do 
		write(50,'(2/," ---------------- ppcx(x,q0,x0) ---------------- ")') 
		do x0=1,nx		
			do q0=1,nq		
				!do x=1,nxs	
				!	write(50,'(1x,tr5,"x0",2(tr5,"xs0",tr4,"eds0",tr3,"exps0")     )') 
				!	write(50,'(7i8)') x0,xx2x(1,x0),xx2e(1,x0),xx2r(1,x0),xx2x(2,x0),xx2e(2,x0),xx2r(2,x0)
				!	write(50,*)
				!	write(50,'(1x,tr5,"q0",2(tr5,"qs0",tr5,"ws0",tr5,"ls0")     )') 
				!	write(50,'(7i8)') q0,qq2q(1,q0),qq2w(1,q0),qq2l(1,q0),qq2q(2,q0),qq2w(2,q0),qq2l(2,q0)
				!	write(50,*)
				!	write(50,'(1x,tr6,"x",2(tr6,"xs",tr5,"eds",tr4,"exps"),5x,tr1,"ppcx(x,q0,x0)",tr6,"sum=1.0?" )') 
				!	write(50,'(7i8,5x,2f14.3)') x,xx2x(1,x),xx2e(1,x),xx2r(1,x),xx2x(2,x),xx2e(2,x),xx2r(2,x),ppcx(x,q0,x0),sum(ppcx(:,q0,x0)) 
				!	write(50,*)
				!	write(50,*)
				!end do 
				if (  qq2w(1,q0)<=np1 .and.  qq2w(2,q0)<=np1 ) then 
					if ( abs( sum(ppcx(:,q0,x0)) - 1.0_dp  ) >  1.0d-6) then  
						print*, "error in getppcx: ppcx for w0=1:np1 does not add up ", q0 , sum(ppcx(:,q0,x0)) 
						stop 
					end if 
				end if 
				if ( qq2w(1,q0)>np1 .or.  qq2w(2,q0)>np1 ) then 
					if ( sum(ppcx(:,q0,x0)) > 1.0d-6 ) then 
						print*, "error in getppx: ppcx for w0=np2 is not 0 " , q0,x0,sum(ppcx(:,q0,x0)) 
						stop 
					end if 
				end if 
			end do 
		end do 
		
		if (icheck_probs) then 
			do x0=1,nx				
				do q0=1,nq			
					if ( maxval(qq2w(:,q0)) <= np1 ) then 
						ppcqx_sum_overqx=0.0_dp
						do x=1,nx
							do q=1,nq
								ppcqx_sum_overqx = ppcqx_sum_overqx + ppcq(q,q0,x0) * ppcx(x,q0,x0)      ! ppcq*ppcx is the joint probability of drawing q,x given q0,x0. ppcq(q,q0,x0) has no x in it because it's independent of x and ppcx doesn't have q in i tbecause is indeoendent of q .
							end do 
						end do 
						if ( abs(ppcqx_sum_overqx-1.0_dp) > eps ) then 
							print*, " pp does not sum to 1 ", ppcqx_sum_overqx
							stop 
						end if 
					end if 
				end do 
			end do 
		end if 
		
	end if 
	end subroutine yaz0

	subroutine  yaz1(index,age)
		integer(i4b), intent(in) :: index,age
		integer(i4b) :: q0,q,x,j,iepsmove
		real(dp) :: emp_m(2),emp_f(2) 
!		do x=1,nx
!			do q=1,nq				
!				if (age<mxa) then ; ageprime=age+1 ; else if (age==mxa) then ; ageprime=age ; end if 
!				write(300,'(tr4,"in",tr3,"age",tr5,"q",tr5,"x",tr7,"utilm_c",tr1,"emaxm_c(ia+1)",tr9,"vm0_c",tr3,"emaxm_c(ia)")')
!				write(300,'(4i6,4f14.2)') index,age,q,x,utilm_c(q,x),emaxm_c(q,x,ageprime),vm0_c(q,x,age),emaxm_c(q,x,age)
!				write(300,*)
!
!			end do 
!		end do 

		!if (age<mxa) then 
		!   j=age+1
		!    do x=1,nxs
		!	    do q=1,nqs				
		!		    write(300,'(tr4,"in",tr3,"age",tr5,"q",tr5,"x",tr8,"emax_s")')
		!		    write(300,'(4i6,2f14.2)') index,age,q,x,emaxm_s(q,x,j),emaxf_s(q,x,j)
		!		    write(300,*)
		!	    end do 
		!   end do 
		!end if 
		
		!print*, "hrtr id mr ",index,myindex

		emp_m=0.0_dp 
		emp_f=0.0_dp
		do x=1,nxs
		do q0=1,nqs	
			do q=1,nqs	
                do iepsmove=1,nepsmove
				if (q2w(q0)<=np1) then 
					j=decm_s(iepsmove,q,x,q0,age,index)
					emp_m(1)=emp_m(1)+1 
					emp_m(2)=emp_m(2)+one(j<=np)
					!write(300,'(tr2,"in",tr1,"age",tr3,"x",tr4,"q0",tr5,"q",tr3,"dec",tr4,"w0",tr5,"w",tr3,"dec",tr4,"l0",tr5,"l",tr3,"dec",tr6,"emax")')
					!write(300,'(3i4,9i6,f10.2)') index,age,x,q0,q,j,q2w(q0),q2w(q),q2w(j),q2l(q0),q2l(q),q2l(j),emaxm_s(q,x,age)
					j=decf_s(iepsmove,q,x,q0,age,index)
					emp_f(1)=emp_f(1)+1 
					emp_f(2)=emp_f(2)+one(j<=np)
					!write(300,'(3i4,9i6,f10.2)') index,age,x,q0,q,j,q2w(q0),q2w(q),q2w(j),q2l(q0),q2l(q),q2l(j),emaxf_s(q,x,age)
					if (icheck_eqvmvf .and. decm_s(iepsmove,q,x,q0,age,index)/=decf_s(iepsmove,q,x,q0,age,index)) then 
						print*, "decsingle not equal! " 
						stop 
					end if 
					!write(300,*)
				end if 
                end do !iepsmove
			end do 
		end do 
		end do 
		
		if (age==30) then 
		if (emp_m(1)>0.0_dp ) then
		 !   write(*,'("m: age,emp ",i4,f10.2,i4)') age,emp_m(2)/emp_m(1),iter
			empdec(age,1)=emp_m(2)/emp_m(1)
		end if 
		if (emp_f(1)>0.0_dp ) then 
		 !   write(*,'("f: age,emp ",i4,f10.2,i4)') age,emp_f(2)/emp_f(1),iter
			empdec(age,2)=emp_f(2)/emp_f(1)
		end if 
		end if 
	end subroutine yaz1

	subroutine check_eqvmvf(dd,vm_c,vf_c)
	integer(i4b), intent(in) :: dd(:)
	real(dp), dimension(nq,nx), intent(in) :: vm_c,vf_c
	real(dp) :: vec(4)
	integer(i4b) :: ia,k,q,x,z,q0,g,j,i
	integer(i4b) :: qr,xr,qp,xp
	ia=dd(1) 
	k=dd(2) 
	q=dd(3) 
	x=dd(4) 
	z=dd(5) 
	q0=dd(6)
	g=dd(7) 
	!j=dd(8) 
	if ( qq2q(1,q0)==qq2q(2,q0)  ) then 
	do xr=1,nxs
		do qr=1,nqs
			do xp=1,nxs
				do qp=1,nqs
					vec=pen
					q=q2qq(qp,qr) ; x=x2xx(xp,xr)
					if (q>0 .and. x>0) then 
						vec(1)=vm_c(q,x) 
						vec(2)=vf_c(q,x)
					end if 								
					i=q2qq(qr,qp) ; j=x2xx(xr,xp)
					if (i>0 .and. j>0) then 
						vec(3)=vm_c(i,j) 
						vec(4)=vf_c(i,j)
					end if 
					if ( abs(vec(1)-vec(4))>eps .and. abs(vec(2)-vec(3))>eps  ) then 
						write(*,'("not good",I4,4F14.2)') ia,vec(1:4)
						write(*,'(4i4)') qq2w(1,q0),qq2l(1,q0),qq2w(2,q0),qq2l(2,q0)
						write(*,'(4i4)') qq2w(1,q),qq2l(1,q),qq2w(2,q),qq2l(2,q)
						write(*,'(4i4)') q2w(qp),q2l(qp),q2w(qr),q2l(qr)
						write(*,'(6i4,2f14.2,4i4)') qr,xr,qp,xp,q,x,vec(1:2)
						write(*,*) ppc(q,q0),maxval(ppcq(q,q0,:))
						write(*,*)
						write(*,'(4i4)') qq2w(1,q0),qq2l(1,q0),qq2w(2,q0),qq2l(2,q0)
						write(*,'(4i4)') qq2w(1,i),qq2l(1,i),qq2w(2,i),qq2l(2,i)
						write(*,'(4i4)') q2w(qr),q2l(qr),q2w(qp),q2l(qp)
						write(*,'(6i4,2f14.2,4i4)') qp,xp,qr,xr,i,j,vec(3:4)
						write(*,*) ppc(i,q0),maxval(ppcq(i,q0,:))
						stop
					end if 
								
				end do 
			end do 
		end do 
	end do 
	end if 
	end subroutine
	
	SUBROUTINE yazsimtemp(dd)    
	integer(i4b), intent(in) :: dd(:)
	integer(i4b) :: ia,k,q,x,z,q0,g,j,altq,i
	integer(i4b) :: a,b,c,qbar(2),kid
	ia=dd(1) 
	k=dd(2) ! index
	q=dd(3) 
	x=dd(4) 
	z=dd(5) 
	q0=dd(6)
	g=dd(7) 
	j=dd(8) ! jmax
	i=dd(9) ! qmax
	kid=0
	if (whereamI==1) then 
		write(200,'(" -------------- SOL: DIV DIV DIV DIV DIV DIV DIV : JCHO,ICHO ",2I8, "--------------")') J,I
	else if (whereamI==4) then
		write(200,'(" -------------- SIM: DIV DIV DIV DIV DIV DIV DIV : JCHO,ICHO ",2I8, "--------------")') J,I
	end if 
	!if (insol) then 
	!	write(200,'(" -------------- SOL: MAR MAR MAR MAR MAR MAR MAR : JCHO,ICHO ",2I8, "--------------")') J,I
	!else 	
	!	write(200,'(" -------------- SIM: MAR MAR MAR MAR MAR MAR MAR : JCHO,ICHO ",2I8, "--------------")') J,I
	!end if 	
	10 FORMAT (/1x,TR5,'Q0',TR6,'X0',TR6,'IA',TR6,'IN',TR7,'Q',TR7,'X', /)
	20 FORMAT (/1x,TR5,'Q0',TR6,'X0',TR6,'IA',TR6,'IN',2(TR6,'W0'),2(TR6,'L0'),TR9,'EMAXC',TR9,'EMAXS' /)
	END SUBROUTINE yazsimtemp

	
end module myaz

	!else if (tag==4) then 
	!	write(200,'(" *************************************** check emax 1 ********************************************* ")') 
	!	write(200,10) 
	!	write(200,'(6i8,5f14.4)') q0,x0,ia,in,q,x,ppcq(q,q0,x0),ppcx(x,q0,x0)  !,val(1),ppcq(q,q0,x0)*ppcx(x,q0,x0)*val(1)
	!	write(200,*) 	
	!else if (tag==5) then 
	!	write(200,'(" *************************************** check emax 2 ********************************************* ")') 
	!	do q0=1,70
	!		x0=2
	!		write(200,20) 
	!		g=1 ; i=qq2q(g,q0) ; n=xx2x(g,x0) 
	!		write(200,'(8i8,2f14.2)') q0,x0,ia,in,qq2w(:,q0),qq2l(:,q0),emaxm_c(q0,x0,ia,in),emaxm_s(i,n,ia,in)
	!		g=2 ; i=qq2q(g,q0) ; n=xx2x(g,x0) 
	!		write(200,'(8i8,2f14.2)') q0,x0,ia,in,qq2w(:,q0),qq2l(:,q0),emaxf_c(q0,x0,ia,in),emaxf_s(i,n,ia,in)
	!		write(200,*) 
	!	end do 



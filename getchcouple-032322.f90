subroutine getch_couple
    ! choice vectors that indicate whether an alternative is feasible given the status quo and given the offer situation
    integer(i4b) :: i,j,a,b,cc,qc(2),is(2),js(2),l(2), d(nc),k

    !First nl options are being unemployed in each location
    cc=0
    do i=1,nl							    ! choice index
        cc=cc+1
        w(1:2)=np1							! w=np+1 denotes unemployment, i is location option, and j is the combo
        j=wl2q(np1,i)						! get what choosing unemployment (w=np+1) at location l corresponds to in terms of the composite index q
        ch(cc,:,:)=qq2q(j,j)                ! get the couple q as the combined index of the q of each partner (i.e. their j in this case). moving anywhere unemployed is always an option		
        !chs(c,:,:)=j									
        !dum(1,c)=w
        !dum(2,c)=i
    end do 	
    if (skriv) then 		
        if ( cc /= nl ) then ; print*, "cc is not equal to nl! ",cc,nl ; stop ; end if 		
    end if 
                
    
    
    qq0: do i=1,nq
        if (  qq2w(1,i)<=np1 .or.  qq2w(2,i)<=np1 )  then   !State space part of qspace
            qq: do j=1,nq			
            w0(1:2)=qq2w(1:2,i)         !get each partner's w0 from q0 (i) 
            l0=qq2l(1,i)                !get l0 from q0 (i)
            w(1:2)=qq2w(1:2,j)          !get each partner's w from q (j)
            l=qq2l(1,j)                 !get l from q (j)

            !if (w0(1)<=np.and.w0(2)<=np) then                   !hub employed and wife employed
            if (l==l0) then                                 !curloc draw
                !if (w(1)<=np.and.w(2)<=np) then             !SCENARIO 1: BOTH GET WAGE OFFER FROM CURRENT LOCATION L0
                                                            !ahuapr22: it should not make a difference whether you do this last if or not
                
                ch(nl+1,j,i) = i                        !hub keeps current w0   / wife keeps current w0 !status quo (both keep their current w0 i.e. reject their offers)
                
                qcm=wl2q(w0(1),l) ; qcf=wl2q(w(2),l)    !hub keeps current w0   / wife accepts offer
                ch(nl+2,j,i) = q2qq(qcm,qcf)            !hub keeps current w0   / wife accepts offer
                qcm=wl2q(w0(1),l) ; qcf=wl2q(np1,l)     !hub keeps current w0   / wife goes into u
                ch(nl+3,j,i) = q2qq(qcm,qcf)            !hub keeps current w0   / wife goes into u


                qcm=wl2q(w(1),l) ; qcf=wl2q(w0(2),l)    !hub accepts offer / wife keeps current w0
                ch(nl+4,j,i) = q2qq(qcm,qcf)            !hub accepts offer / wife keeps current w0
                qcm=wl2q(w(1),l) ; qcf=wl2q(w(2),l)     !hub accepts offer / wife accepts offer
                ch(nl+5,j,i) = q2qq(qcm,qcf)            !hub accepts offer / wife accepts offer
                qcm=wl2q(w(1),l) ; qcf=wl2q(np1,l)      !hub accepts offer / wife goes into u 
                ch(nl+6,j,i) = q2qq(qcm,qcf)            !hub accepts offer / wife goes into u

                qcm=wl2q(np1,l) ; qcf=wl2q(w0(2),l)     !hub goes into u / wife keeps current w0
                ch(nl+7,j,i) = q2qq(qcm,qcf)            !hub goes into u / wife keeps current w0
                qcm=wl2q(np1,l) ; qcf=wl2q(w(2),l)      !hub goes into u / wife accepts offer
                ch(nl+8,j,i) = q2qq(qcm,qcf)            !hub goes into u / wife accepts offer
                qcm=wl2q(np1,l) ; qcf=wl2q(np1,l)       !hub goes into u / wife goes into u
                ch(nl+9,j,i) = q2qq(qcm,qcf)            !hub goes into u / wife goes into u
    
                if (w0(1)==np1) then    !if hub is currently not working then certain options are not in the choice set 
                    ch(nl+1,j,i)=0      !ahuapril22: check if all is the same even when these are not set to 0 since it should not really matter
                    ch(nl+2,j,i)=0 
                    ch(nl+3,j,i)=0
                end if
                if (w0(2)==np1) then    !if wife is currently not working then certain optiosn are not int eh choice set 
                    ch(nl+1,j,i)=0      !ahuapril22: check if all is the same even when these are not set to 0 since it should not really matter
                    ch(nl+4,j,i)=0 
                    ch(nl+7,j,i)=0
                end if 
                if (w(1)==np2) then     !if hub gets laid off then keeping current w0 is not an option in th choice set nor are the accepting offer options
                    ch(nl+1:nl+6,j,i)=0 
                end if 
                if (w(2)==np2) then     !if wife gets laid off then keeping current w0is not an option in th choice set nor are the accepting offer options
                    ch(nl+1:nl+2,j,i)=0 
                    ch(nl+4:nl+5,j,i)=0 
                    ch(nl+7:nl+8,j,i)=0 
                end if

            else if (l/=l0) then                            !ofloc draw
                !STAY
                ch(nl+1,j,i) = i                        !hub keeps current w0   / wife keeps current w0 !status quo (both keep their current w0 i.e. reject their offers)
                
                qcm=wl2q(w0(1),l) ; qcf=wl2q(w(2),l)    !n/a:hub keeps current w0   / wife accepts offer
                ch(nl+2,j,i) = q2qq(qcm,qcf)            !n/a:hub keeps current w0   / wife accepts offer
                qcm=wl2q(w0(1),l0) ; qcf=wl2q(np1,l0)   !hub keeps current w0   / wife goes into u STAY
                ch(nl+3,j,i) = q2qq(qcm,qcf)            !hub keeps current w0   / wife goes into u STAY

                qcm=wl2q(w(1),l) ; qcf=wl2q(w0(2),l)    !n/a: hub accepts offer / wife keeps current w0 MOVE
                ch(nl+4,j,i) = 0                        !n/a: hub accepts offer / wife keeps current w0 MOVE
                qcm=wl2q(w(1),l) ; qcf=wl2q(w(2),l)     !hub accepts offer / wife accepts offer MOVE 
                ch(nl+5,j,i) = q2qq(qcm,qcf)            !hub accepts offer / wife accepts offer MOVE 
                qcm=wl2q(w(1),l) ; qcf=wl2q(np1,l)      !hub accepts offer / wife goes into u MOVE 
                ch(nl+6,j,i) = q2qq(qcm,qcf)            !hub accepts offer / wife goes into u MOVE
                
                qcm=wl2q(np1,l0) ; qcf=wl2q(w0(2),l0)   !hub goes into u / wife keeps current w0 STAY
                ch(nl+7,j,i) = q2qq(qcm,qcf)            !hub goes into u / wife keeps current w0 STAY
                qcm=wl2q(np1,l) ; qcf=wl2q(w(2),l)      !hub goes into u / wife accepts offer MOVE
                ch(nl+8,j,i) = q2qq(qcm,qcf)            !hub goes into u / wife accepts offer MOVE
                qcm=wl2q(np1,l0) ; qcf=wl2q(np1,l0)     !hub goes into u / wife goes into u STAY
                ch(nl+9,j,i) = q2qq(qcm,qcf)            !hub goes into u / wife goes into u STAY

                if (w0(1)==np1) then    !if hub is currently not working then certain options are not in the choice set 
                    ch(nl+1,j,i)=0      !ahuapril22: check if all is the same even when these are not set to 0 since it should not really matter
                    ch(nl+2,j,i)=0 
                    ch(nl+3,j,i)=0
                end if
                if (w0(2)==np1) then    !if wife is currently not working then certain optiosn are not int eh choice set 
                    ch(nl+1,j,i)=0      !ahuapril22: check if all is the same even when these are not set to 0 since it should not really matter
                    ch(nl+4,j,i)=0 
                    ch(nl+7,j,i)=0
                end if 
                if (w(1)==np2) then     !if hub gets laid off then keeping current w0 is not an option in th choice set nor are the accepting offer options
                    ch(nl+1:nl+6,j,i)=0 
                end if 
                if (w(2)==np2) then     !if wife gets laid off then keeping current w0is not an option in th choice set nor are the accepting offer options
                    ch(nl+1:nl+2,j,i)=0 
                    ch(nl+4:nl+5,j,i)=0 
                    ch(nl+7:nl+8,j,i)=0 
                end if

            end if !l and l0 comparison
        
        end if !state space part of q0
        end do !q
    end do !q0

end subroutine getch_couple
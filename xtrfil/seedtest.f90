program seedtest
    
    implicit none 

    real(dp) :: rand
    integer :: i
    !seedo=98765
	!seedo=seed_init    
	
    call random_seed (size=p)
    !p=12
    call random_seed (put = (/(k,k=1,p)/))
    !seed_init=3
    !call random_seed (put = seed_init )
	
    print*, 'Here is p',p
    print*, 'Here is k',(/(k,k=1,p)/)
	
    
    do i=1,10
    call random_number(rand)
    print*, "Here is rand", rand
    end do 
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

    end program seedtest 

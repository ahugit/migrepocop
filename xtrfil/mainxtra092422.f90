




    !pars(70)=ptypehs1  !set to 0 in getpars and kept fixed
    !pars(71)=ptypecol1 !set to 0 in getpars and kept fixed
    !pars(72)=alf1t1 
    !pars(73)=alf1t1
    !pars(74)=cst1  !is 0 (set to 0 in getpars) and should be kept fixed
    !pars(75)=0.0015_dp !mumar1 !should this be set to 0? 
    
    !pars(76)=0.0_dp !ptypehs2 
    !pars(77)=0.0_dp !ptypecol2
    !pars(78)=alf1t2 
    !pars(79)=alf1t2
    pars(80)=4000.0_dp !cst2 
    pars(81)=0.0015_dp !mumar2
    
    pars(82)=0.0_dp !ptypehs3 
    pars(83)=0.0_dp !ptypecol3
    !pars(84)=alf1t 3
    !pars(85)=alf1t 3
    pars(86)=4500.0_dp !cst3
    pars(87)=0.0015_dp !mumar3

    pars(88)=0.0_dp !ptypehs4 
    pars(89)=0.0_dp !ptypecol4
    !pars(90)=alf1t 4
    !pars(91)=alf1t 4
    pars(92)=5000.0_dp !cst4 
    pars(93)=0.0015_dp !mumar4
    call getpars(pars,realpars)
    call objfunc(pars,qval) ; realpars=realpartemp     


*************************

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

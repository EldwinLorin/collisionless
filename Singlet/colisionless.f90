program collisionless
    implicit NONE


!Let's declare many variables here later
    integer i, j,  nproducts, nTS, nbins, Jmax, DeltaJ, nJ, k
    character*30  input_energy_unit 


real(kind=16) Emax, DeltaE, T, Qfrag, kb, integral, h, rateconst, qts, E_independant_factor,&
integral_CCH2,integral_HCCH,rateconst_CCH2,rateconst_HCCH, kcal_to_invcm 

integer, Allocatable ::  angularmomentum(:), lines_to_add(:)


character(len=30), Allocatable ::  ts_dos(:)
character(len=100) :: filename_kJ, filename_results, filename_allk

real(kind=16), Allocatable ::  product_energy(:), ts_energy(:), energy(:),&
 denominator(:), frac(:),kTJ(:), from_CCH2(:), from_HCCH(:),&
kTJ_CCH2(:),kTJ_HCCH(:)
real(kind=16), Allocatable :: Energybin(:,:)
real, Allocatable :: sumstates(:,:,:)

real(kind=16) min_E,min_prod,min_TS!Those variables are used to determine lowest energy value

!Conversion factor
kcal_to_invcm=349.757
kb=0.00831446261815324*83.593
h=83.593*6.02214076*6.62607015d-14  ! cm-1.s-1 6.62607015d-34 

!Read the input file
    open (unit=1, file='collisionless.inp',status='old')

    read(1,*) Emax, DeltaE  !Max energy and grain size (cm-1)
    read(1,*) Jmax, DeltaJ  !Max J and J step
    read(1,*) T, Qfrag      !Temperature and down energy transfer value (cm-3)

    read(1,*) input_energy_unit  !Probably we'll only support kcal/mol
    read(1,*) nproducts   !Number of bimolecular products

ALLOCATE(product_energy(nproducts))


    do i=1, nproducts !Reading products data
        read(1,*) product_energy(i)  !product_index(i), product_name(i), 
    end do

    read(1,*) nTS !Reading TS data
    ALLOCATE(ts_dos(nTS), &
    ts_energy(nTS))

    do i=1, nTS
        read(1,*) ts_dos(i), &
         ts_energy(i)
    end do

    close(1)

    !Reading main input file over. Now we process it
    !Creating energy scale
    nbins=int(Emax/DeltaE)
    ALLOCATE(energy(nbins))
    do i=1, nbins
        energy(i)=(i-1)*DeltaE
    end do
    !Creating J scale
    nJ= Jmax/DeltaJ
    Allocate(angularmomentum(nJ))
    do i=1, nJ
        angularmomentum(i)=i*DeltaJ
    end do
    Allocate(kTJ(nJ))
    kTJ=0

    !Let's find the lowest energy point and define all others relative to it
    min_prod=MINVAL(product_energy)
    min_TS=MINVAL(ts_energy)
    min_E=MIN(min_prod,min_TS)

    product_energy=product_energy-min_E
    ts_energy=ts_energy-min_E

    !Now we convert all energies to cm-1
    product_energy=product_energy*kcal_to_invcm
    ts_energy=ts_energy*kcal_to_invcm

    !I think we want all eneregies to be a multiple of our deltaE, so if a TS is a E=114cm-1 and deltaE=10cm-1, we want that TS to be actually at 110cm-1
    do i=1, nproducts
        product_energy(i)=real(deltaE*int(product_energy(i)/deltaE),8)
    end do
    do i=1, nTS
        ts_energy(i)=real(deltaE*int(ts_energy(i)/deltaE),8)
    end do


    !Now we want to read the DOS files
    Allocate(Energybin(nts,nbins),  sumstates(nts,nbins,nJ))
    do i = 1, nTS

        open(unit=20, file=ts_dos(i), status='old')
        do j = 1, nbins
            read(20,*)  Energybin(i,j), (sumstates(i,j,k), k=1,nJ)
            Energybin(i,j)=Energybin(i,j)+ts_energy(i)
        end do
        close(20)
    end do
    do j=1,nbins
    end do
do j=1,nbins
end do
    !Now all my DOS files start from the 0 of energy of each specie, we changed them to start from the energy of the species with respect to the bottom
    !of the system. We have to make them all match by adding E=0, deltaE, 2*delta E etc with sum of states being 0 lines at the start, and cutting the
    !extra lines in the bottom
    !So we count the number of lines to add at the start (which is the same as the number of lines we'll have to trim at the bottom)
    Allocate(lines_to_add(nTS))
    lines_to_add=0
    do i=1, nTS
        do j = nbins, 1, -1
            if (Energybin(i,j) .ge. Emax) THEN
                lines_to_add(i) = lines_to_add(i)+1
            endif
        end do
    end do
    !Then we move the lines to the bottom and add extra lines at the start
    do i=1, nTS
        if (lines_to_add(i).gt.0) then
            do j = nbins, lines_to_add(i)+1, -1
            energybin(i,j)=energybin(i,j)-lines_to_add(i)*deltaE
            sumstates(i,j,:)=sumstates(i,j-lines_to_add(i),:)
        end do

        do j=1, lines_to_add(i)
            Energybin(i,j)=(j-1)*deltaE
            sumstates(i,j,:)=0
        end do
    end if
    end do

    !Now I am not sure how to generalize that so this will be only tailored to my system. We will use J dependant DoS, and model is collisionless
    Allocate(denominator(nbins),frac(nbins), from_CCH2(nbins), from_HCCH(nbins))
    Allocate(kTJ_CCH2(nJ),kTJ_HCCH(nJ))
 Qts=0
 E_independant_factor=kb*T/(Qfrag*h)
 
 write(filename_kJ, '(A,F6.1,A)') 'kJ_', T, '.out'
 write(filename_results, '(A,F6.1,A)') 'results_', T, '.out'
 write(filename_allk, '(A,F6.1,A)') 'all_k_', T, '.out'
 
 open(unit=100, file=filename_kJ, status='replace')
 open(unit=200, file=filename_results, status='replace')
 open(unit=300, file=filename_allk, status='replace')

 write(100,'(A5, A15, A15, A15, A15)') 'J', 'k(T,J)', 'kCCH2(T,J)', 'kHCCH(T,J)', 'kback(T,J)'

    do k=1, nJ

        denominator=sumstates(3,:,k)*(sumstates(1,:,k)+sumstates(4,:,k)+sumstates(2,:,k))+sumstates(2,:,k)&
        *(sumstates(1,:,k)+sumstates(4,:,k))
        from_CCH2=sumstates(4,:,k)*(sumstates(3,:,k)+sumstates(2,:,k))*sumstates(1,:,k)/denominator
        from_HCCH=sumstates(3,:,k)*sumstates(2,:,k)*sumstates(1,:,k)/denominator

        frac=sumstates(1,:,k)

    do i=1, nbins
        if (energy(i) .ge. product_energy(1)) then
        frac(i)=frac(i)*deltaE*exp(-(energy(i)-product_energy(1))/(kb*T))/(kb*T)
        from_CCH2(i)=from_CCH2(i)*deltaE*exp(-(energy(i)-product_energy(1))/(kb*T))/(kb*T)
        from_HCCH(i)=from_HCCH(i)*deltaE*exp(-(energy(i)-product_energy(1))/(kb*T))/(kb*T)
        
        else
        frac(i)=0
        from_CCH2(i)=0
        from_HCCH(i)=0
        endif
    end do

    integral=0
    integral_CCH2=0
    integral_HCCH=0
    do i=1, nbins
            integral=integral+frac(i)*(2*k+1)
            integral_CCH2=integral_CCH2+from_CCH2(i)*(2*k+1)
            integral_HCCH=integral_HCCH+from_HCCH(i)*(2*k+1)
    end do

    Qts=Qts+integral
    kTJ_CCH2(k)=integral_CCH2*E_independant_factor
    kTJ_HCCH(k)=integral_HCCH*E_independant_factor
    kTJ(k) = integral*E_independant_factor
    write(100,'(I5, 1PE15.5, 1PE15.5, 1PE15.5, 1PE15.5)') k, kTJ(k), kTJ_CCH2(k), kTJ_HCCH(k), kTJ(k)-(kTJ_CCH2(k)&
    +kTJ_HCCH(k))
end do

rateconst_CCH2=sum(kTJ_CCH2)
rateconst_HCCH=sum(kTJ_HCCH)
rateconst=sum(kTJ)

write(200,'(A, T50, F15.7)') 'Results for T=', T
write(200,'(A, T50, 1PE15.7)') 'k(T)=', rateconst
write(200,'(A, T50, 1PE15.7)') 'kprod(T)=', rateconst_CCH2+rateconst_HCCH
write(200,'(A, T50, 1PE15.7)') 'kCCH2(T)=', rateconst_CCH2
write(200,'(A, T50, 1PE15.7)') 'kHCCH(T)=', rateconst_HCCH
write(200,'(A, T50, 1PE15.7)') 'kback(T)=', rateconst - (rateconst_CCH2 + rateconst_HCCH)
write(200,'(A, T50, F15.7)') 'Branching ratio for products=', (rateconst_CCH2 + rateconst_HCCH)/rateconst
write(200,'(A, T50, F15.7)') 'Branching ratio for backdissociation=', 1-((rateconst_CCH2 + rateconst_HCCH)/rateconst)
write(200,'(A, T50, 1PE15.7)') 'Final entrance TS partition function value=', Qts

write(300,'(F15.7)') T
write(300,'(1PE15.7)') rateconst
write(300,'(1PE15.7)') rateconst_CCH2+rateconst_HCCH
write(300,'(1PE15.7)') rateconst_CCH2
write(300,'(1PE15.7)') rateconst_HCCH
write(300,'(1PE15.7)') rateconst - (rateconst_CCH2 + rateconst_HCCH)
write(300,'(F15.7)') (rateconst_CCH2 + rateconst_HCCH)/rateconst
write(300,'(F15.7)') 1-((rateconst_CCH2 + rateconst_HCCH)/rateconst)
write(300,'(1PE15.7)') Qts

    deALLOCATE(kTJ)
    DEALLOCATE(lines_to_add,Energybin, sumstates, angularmomentum,energy)
    DEALLOCATE(ts_dos,ts_energy)
    deALLOCATE(product_energy)
    !Maybe some variables are not properly deallocated (denominator at least), but we don't care because program ends. Probably not good practice though

end program




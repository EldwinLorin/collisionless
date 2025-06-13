program collisionless
    use, intrinsic :: ieee_arithmetic
    implicit NONE

!Let's declare many variables here
    integer i, k,  nproducts, nTS, nbins, Jmax, DeltaJ, nJ, j, l, nreact, nwell, info,&
    fr, too,nchan,rho_int
    character*30  input_energy_unit 
    character(len=50) :: fmtstr


real(kind=16) Emax, DeltaE, T, Qfrag, kb, integral_cap, h, k_cap, qts, E_independant_factor,&
 kcal_to_invcm,sm_integral_prod,kTJ, sum_entrance_TS,f_func

integer, Allocatable ::  angularmomentum(:), lines_to_add(:), ipiv(:), from(:), to(:),&
rho(:),rho_ts(:)


character(len=30), Allocatable ::  ts_dos(:)
character(len=100) :: filename_kJ, filename_results, filename_allk

real(kind=16), Allocatable ::  product_energy(:), ts_energy(:), energy(:),&
sm_integral_cap,kTJ_CCH2(:),kTJ_HCCH(:),k_prod(:),&
A(:,:), B(:), X(:), integral_prod(:), ktJ_prod(:)
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
    read(1,*) Jmax, DeltaJ  !Max k and k step
    read(1,*) T, Qfrag      !Temperature and down energy transfer value (cm-3)

    read(1,*) input_energy_unit  !Probably we'll only support kcal/mol
    read(1,*) nreact, nwell, nproducts   !Number of reactants, wells, bimolecular products
    !reactants are numbered from 1 to nreactants, wells are numbered from nreactants + 1 to nreactants + nwells
    !products are numbered from nreactants + nwells + 1 to nreactants + nwells + nproducts
    !reactants are considered as source term
    !products are irreversible sinks
    !so reactants should also be included in products


ALLOCATE(product_energy(nproducts),rho(nwell))
rho=0

    do i=1, nproducts !Reading products data
        if (i.eq.1) then
            read(1,*) product_energy(i), rho_int
        else
        read(1,*) product_energy(i)  !product_index(i), product_name(i), 
        endif
    end do

    do i=1, nwell
        read(1,*) rho(i)
    end do

    read(1,*) nTS !Reading TS data
    ALLOCATE(ts_dos(nTS), ts_energy(nTS), from(nTS), to(nTS),rho_ts(nTS))

    do i=1, nTS
        read(1,*) ts_dos(i), ts_energy(i), from(i), to (i), rho_TS(i)
    end do

    close(1)
    !write(*,*) rho_int, rho, rho_TS

    !Reading main input file over. Now we process it
    !Creating energy scale
    nbins=int(Emax/DeltaE)
    ALLOCATE(energy(nbins))
    do i=1, nbins
        energy(i)=(i-1)*DeltaE
    end do
    !Creating k scale
    nJ= Jmax/DeltaJ
    Allocate(angularmomentum(nJ))
    do i=1, nJ
        angularmomentum(i)=i*DeltaJ
    end do
!    Allocate(kTJ(nJ))
!    kTJ=0

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
        do k = 1, nbins
            read(20,*)  Energybin(i,k), (sumstates(i,k,j), j=1,nJ)
            Energybin(i,k)=Energybin(i,k)+ts_energy(i)
        end do
        close(20)
    end do
    do k=1,nbins
    end do
do k=1,nbins
end do
    !Now all my DOS files start from the 0 of energy of each specie, we changed them to start from the energy of the species with respect to the bottom
    !of the system. We have to make them all match by adding E=0, deltaE, 2*delta E etc with sum of states being 0 lines at the start, and cutting the
    !extra lines in the bottom
    !So we count the number of lines to add at the start (which is the same as the number of lines we'll have to trim at the bottom)
    Allocate(lines_to_add(nTS))
    lines_to_add=0
    do i=1, nTS
        do k = nbins, 1, -1
            if (Energybin(i,k) .ge. Emax) THEN
                lines_to_add(i) = lines_to_add(i)+1
            endif
        end do
    end do

    !Then we move the lines to the bottom and add extra lines at the start
    do i=1, nTS
        if (lines_to_add(i).gt.0) then
            do k = nbins, lines_to_add(i)+1, -1
            energybin(i,k)=energybin(i,k)-lines_to_add(i)*deltaE
            sumstates(i,k,:)=sumstates(i,k-lines_to_add(i),:)
        end do

        do k=1, lines_to_add(i)
            Energybin(i,k)=(k-1)*deltaE
            sumstates(i,k,:)=0
        end do
    end if
    end do

    nchan=0
    do l=1, nTS  !Counting the number of exit channels
        if(to(l).gt.nreact+nwell) then
            nchan=nchan+1
        endif
    end do
write(fmtstr,'("(I5,",I0,"(1PE15.5))")') nchan + 1 !Some formating for later

    Allocate(kTJ_CCH2(nJ),kTJ_HCCH(nJ), A(nwell,nwell),B(nwell),X(nwell), ipiv(nwell),&
    integral_prod(nchan),ktJ_prod(nchan), k_prod(nchan))

 !Some output filename formating
 write(filename_kJ, '(A,F0.0,A)') 'kJ_', T, 'out'
 write(filename_results, '(A,F0.0,A)') 'results_', T, 'out'
 write(filename_allk, '(A,F0.0,A)') 'all_k_', T, 'out'
 
 open(unit=100, file=filename_kJ, status='replace')
 open(unit=200, file=filename_results, status='replace')
 open(unit=300, file=filename_allk, status='replace')

 write(100,'(A5, A15, A15, A15)') 'J', 'k(T,J)', 'k_channels(T,J)'

k_prod=0
Qts=0
E_independant_factor=kb*T/(Qfrag*h)

        sum_entrance_TS=0
do j=1, nJ
    do i = 1, nbins
        sum_entrance_TS=sum_entrance_TS+(2*J+1)*(sumstates(1,i,j)*deltaE*exp(-energy(i)/(kb*T)))
    end do
end do
    do j=1, nJ
        ktJ=0
        ktJ_prod=0
        integral_prod=0
        integral_cap=0

        do i=1, nbins
            if (Energy(i) .ge. product_energy(1)) then
            A=0
            B=0
            X=0
            f_func=(2*j+1)*sumstates(1,i,j)*exp(-energy(i)/(kb*T))/sum_entrance_TS
            do l=1, nTS  !Building matrix
                if ((from(l) .gt. nreact .and. from(l) .le. nreact + nwell) .and. &
                (to(l) .gt. nreact .and. to(l) .le. nreact + nwell)) then
                    fr=from(l)-nreact
                    too=to(l)-nreact

                    A(fr,too)=sumstates(l,i,j)*rho(too)/rho_TS(l)
                    A(too,fr)=sumstates(l,i,j)*rho(fr)/rho_TS(l)

                    A(fr,fr)=A(fr,fr)-A(too,fr)
                    A(too,too)=A(too,too)-A(fr,too)
                else if (from(l) .le. nreact) then
                    too=to(l)-nreact
                    B(too)=B(too)-(f_func*rho_int/rho_TS(l))
                else if(to(l) .gt.nreact+nwell) then                   
                    fr=from(l)-nreact
                    A(fr,fr)=A(fr,fr)-sumstates(l,i,j)*rho(fr)/rho_TS(l)
                end if
            end do
                X = B
                call DGESV(nwell, 1, A, nwell, ipiv, X, nwell, info) !Solve linear system

            nchan=0
            !Compute k_cap(E,J)
            sm_integral_cap=sumstates(1,i,j)*deltaE*exp(-(energy(i)-product_energy(1))/(kb*T))/(kb*T)*rho_int/rho_TS(1)
            integral_cap=integral_cap+sm_integral_cap*(2*j+1)    !Sum over E
            do l=1, nTS  !Compute k(E,J) for each product channel
                if(to(l).gt.nreact+nwell) then
                    nchan=nchan+1
                    fr=from(l)-nreact
                    sm_integral_prod=X(fr)*sumstates(l,i,j)*deltaE*exp(-(energy(i)-product_energy(1))/(kb*T))/(kb*T)
                    sm_integral_prod=sm_integral_prod*(2*j+1)*rho(fr)/rho_TS(l)

                    integral_prod(nchan)=integral_prod(nchan)+sm_integral_prod  !Sum over E
                end if
            end do
        end if
        end do


    Qts=Qts+integral_cap

    kTJ_prod=integral_prod*E_independant_factor  !Compute k(J) for each product channel
    kTJ = integral_cap*E_independant_factor      !Compute k_cap(J)

    k_cap=k_cap+kTJ                              !Compute k_cap(T)
    k_prod=k_prod+kTJ_prod                       !Compute k(T) for each product channel

    write(100,fmtstr) j, kTJ, kTJ_prod(1:nchan)
end do


!Only outputs from now on
write(200,'(A, T50, F15.7)') 'Results for T=', T
write(200,'(A, T50, 1PE15.7)') 'k_capture(T)=', k_cap
do i = 1, nchan
    write(200,'(A,I0,A, T50, 1PE15.7)') 'kprod_', i, '(T)=', (k_prod(i)/sum(k_prod))*k_cap
end do

do i = 1, nchan
write(200,'(A,I0,A, T50, F15.7)') 'Branching ratio for products_',i,'=', k_prod(i)/sum(k_prod)
end do
write(200,'(A, T50, 1PE15.7)') 'TS partition function=', Qts


write(300,'(F15.7)') T
write(300,'(1PE15.7)') k_cap
do i = 1, nchan
    write(300,'(1PE15.7)') k_prod(i)
end do
do i = 1, nchan
write(300,'(F15.7)') k_prod(i)/k_cap
end do
write(300,'(1PE15.7)') Qts

    DEALLOCATE(lines_to_add,Energybin, sumstates, angularmomentum,energy)
    DEALLOCATE(ts_dos,ts_energy)
    deALLOCATE(product_energy)
    !Maybe some variables are not properly deallocated but we don't care because program ends. Probably not good practice though

end program




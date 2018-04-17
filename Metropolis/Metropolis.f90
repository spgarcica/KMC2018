module MetropolisHastings
        use mtmod
        use energiesmod
        contains
        subroutine Metropolis(N_particles,From,To,Temperature,Neighbours,Geomat,OccupiedVec,PotentialVec,PolarizationVec, & 
                        InterMat,TotInteraction,TotRan,TotPol,TotEnergy,Switch,EPot,Real_To)
                implicit none
                integer, intent(in) :: N_particles, Neighbours, From, To, Real_To
                integer, dimension(neighbours,N_particles), intent(in) :: Geomat
                integer, dimension(N_particles), intent(inout) :: OccupiedVec
                real(8), dimension(N_particles), intent(in) :: PotentialVec, PolarizationVec
                real(8), dimension(N_particles), intent(inout) :: InterMat
                real(8), intent(in) :: Temperature, EPot
                real(8), intent(inout) :: TotInteraction, TotRan, TotPol, TotEnergy
                integer :: ii, Energychange, Check
                real(8), dimension(4) :: Pivot
                integer, dimension(N_particles) :: OccupiedPivot
                real(8), dimension(N_particles) :: DummyMatrix
                Logical, intent(out) :: Switch

                Switch = .false.
                
       
                ! El cambio de energía se ha de multiplicar por dos para que todo cuadre
                ! debido a que cambia el signo
                OccupiedPivot(:) = OccupiedVec(:)
                OccupiedPivot(To) = 1
                OccupiedPivot(From) = 0
                ! Debido a que las partículas suelen saltar a un lado recalculo la energía en cada paso !
                ! sé que esto no es lo óptimo pero tengo poco tiempo de cálculo                         !
                call InteractionPotentialCalc(N_Particles,Neighbours,OccupiedPivot,Geomat,DummyMatrix,Pivot(1))
                Pivot(2) = (PotentialVec(To) - PotentialVec(From))
                Pivot(3) = (PolarizationVec(To) - PolarizationVec(From))
                Pivot(4) = Real_To*EPot-PolarizationVec(From)
                Energychange = Pivot(1) + Pivot(2) - Pivot(4) - TotInteraction

                if (energychange <= 0) then
                        switch = .true.
                else
                        check = exp(-Energychange/(Temperature)) + grnd()
                        if (check == 1) then
                                switch = .true.
                        end if
                end if

                ! Aquí hacemos un update en caso de ser válido el cambio !
                if (switch) then
                        TotInteraction = Pivot(1)
                        TotRan = TotRan + Pivot(2)
                        TotPol = TotPol + Pivot(3)
                        TotEnergy = TotInteraction + TotRan + TotPol
                        OccupiedVec(:) = OccupiedPivot(:)
                end if
        end subroutine                        
end module

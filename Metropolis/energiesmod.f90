module EnergiesMod
        use mtmod
        use Translate
        contains
                subroutine SitesInit(L_Box,OccupiedVec)
                        integer, intent(in) :: L_Box
                        integer, dimension(:), allocatable, intent(out) :: OccupiedVec
                        integer :: N_Sites, ii, jj, NChanges, SiteSelect
                        logical :: HalfCheck, SecondCheck

                        ! Primero reservamos la memoria para el vector de ocupación !
                        N_Sites = L_Box**2
                        allocate(OccupiedVec(N_Sites))

                        ! Ahora localizamos al azar los espines ! 
                        do ii=1, N_Sites
                                if (grnd() >=0.5) then
                                        OccupiedVec(ii) = 1
                                else
                                        OccupiedVec(ii) = 0
                                end if
                         end do
                         ! Ahora comprobamos si se cumple la condición de tener la mitad de huecos ocupados !
                         HalfCheck = .true.
                         do while (HalfCheck)
                                 if (sum(OccupiedVec) == N_Sites/2) then
                                         HalfCheck = .false.
                                 else if (sum(OccupiedVec) > N_Sites/2) then
                                         ! Elegimos un sitio al azar y quitamos una partícula !
                                         NChanges = sum(OccupiedVec) - N_Sites/2
                                         do ii=1, NChanges
                                                 SecondCheck = .true.
                                                 do while (SecondCheck)
                                                         SiteSelect = N_Sites
                                                         do while (SiteSelect == N_Sites)
                                                                 SiteSelect = (grnd() * N_Sites) + 1
                                                         end do
                                                         if (OccupiedVec(SiteSelect) == 1) then
                                                                 OccupiedVec(SiteSelect) = 0
                                                                 SecondCheck = .false.
                                                         end if
                                                 end do
                                         end do
                                 else if (sum(OccupiedVec) < N_Sites/2) then
                                         ! Elegimos un sitio al azar y añadimos una partícula !
                                         NChanges = N_Sites/2 - sum(OccupiedVec)
                                         do ii=1, NChanges
                                                 SecondCheck = .true.
                                                 do while (SecondCheck)
                                                         SiteSelect = N_Sites
                                                         do while (SiteSelect == N_Sites)
                                                                 SiteSelect = (grnd() * N_Sites) + 1
                                                         end do
                                                         if (OccupiedVec(SiteSelect) == 0) then
                                                                 OccupiedVec(SiteSelect) = 1
                                                                 SecondCheck = .false.
                                                         end if
                                                 end do
                                         end do
                                 end if
                         end do
                end subroutine

                ! Esta subrutina calcula un vector con los valores del potencial aleatorio en cada sitio !
                subroutine RanPotentialInit(N_Sites,WPotential,OccupiedVec,PotentialVec,TotRan)
                        implicit none
                        integer, intent(in) :: N_Sites
                        integer, dimension(N_Sites), intent(in) :: OccupiedVec
                        real(8), intent(in) :: WPotential
                        real(8), intent(out) :: TotRan
                        real(8), dimension(:), allocatable, intent(out) :: PotentialVec
                        integer :: ii

                        allocate(PotentialVec(N_Sites))
                        do ii=1, N_Sites
                                PotentialVec(ii) = grnd() * WPotential
                                if (grnd() >= 0.5) then 
                                        PotentialVec(ii) = PotentialVec(ii) * (-1)
                                end if
                                TotRan = PotentialVec(ii)
                        end do
                end subroutine

                ! Esta subrutina calcula un vector con los valores del potencial E en cada sitio !
                subroutine PolarizationCalc(L_Box,N_Sites,OccupiedVec,EPotential,PolarizationVec,TotPol,Polarization)
                        implicit none
                        integer, intent(in) :: L_Box, N_Sites
                        real(8), intent(in) :: EPotential
                        real(8), intent(inout) :: TotPol, Polarization
                        integer, dimension(N_Sites), intent(in) :: OccupiedVec
                        real(8), dimension(N_Sites), intent(out) :: PolarizationVec
                        integer, dimension(2) :: CartesianCoords
                        integer :: ii
                         
                        Polarization = 0
                        TotPol = 0
                        do ii=1, N_Sites
                                CartesianCoords = Index_To_XY(L_Box,ii)
                                PolarizationVec(ii) = CartesianCoords(1) * EPotential
                                TotPol = TotPol + (PolarizationVec(ii) * OccupiedVec(ii))
                                Polarization = Polarization + OccupiedVec(ii)*CartesianCoords(1)
                        end do
                end subroutine

                ! Esta subrutina calcula un vector con los valores del potencial de interacción en cada sitio !
                subroutine InteractionPotentialCalc(N_Sites,Neighbours,OccupiedVec,Geomat,InterMat,TotInteraction)
                        implicit none
                        integer, intent(in) :: N_Sites, Neighbours
                        integer, dimension(N_Sites), intent(in) :: OccupiedVec
                        integer, dimension(Neighbours,N_Sites), intent(in) :: Geomat
                        real(8), dimension(N_Sites), intent(inout) :: InterMat
                        real(8), intent(inout) :: TotInteraction
                        integer :: ii, jj

                        InterMat = 0
                        TotInteraction = 0
                        do ii = 1, N_Sites
                                do jj = 1, Neighbours
                                        InterMat(ii) = InterMat(ii) + (OccupiedVec(Geomat(jj,ii))-0.5) &
                                                * (OccupiedVec(ii)-0.5)
                                end do
                                TotInteraction = TotInteraction + InterMat(ii)
                        end do
                end subroutine
end module

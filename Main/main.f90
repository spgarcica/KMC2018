program test
      use mtmod
      use Translate
      use TheTower
      use EnergiesMod
      use GenGeom
      use MetropolisHastings
      implicit none
      integer ::  ii, jj, kk, Tower_Sites, Return_Index, Real_To, L_Box, N_Sites, Neighbours
      ! Varbiales de MonteCarlo !
      integer :: From, To, Seed, Snap
      integer(8) :: NSteps, MonteCount, SSteps, SnapCount, NSnap
      integer, dimension(2) :: Cartesian
      integer, dimension(:), allocatable :: OccupiedVec
      integer, dimension(:,:), allocatable :: Geomat, Tower_Mat
      real(8) :: Gamma_T, TotPol, TotInteraction, TotRan, Temperature, TotalEnergy
      real(8) :: EPotential, WPotential, Time, Polarization
      real(8), dimension(:), allocatable :: Tower_Prob_Vec
      real(8), dimension(:), allocatable :: InterMat, PotentialVec, PolarizationVec
      character(3) :: Simulation
      logical :: Change
      
      ! Leemos el archivo de parámetros !
      call ReadFile(L_Box,EPotential,WPotential,Temperature,SSteps,NSteps,Snap,Seed)
      ! Inicializamos el generador de números aleatorios !
      call sgrnd(seed)
      ! Calculamos el número de sititos del sistema !
      N_Sites = L_Box**2
      ! Inicializamos las energías y las matrices de todo el sistema !
      call Initialize(L_Box,N_Sites,Neighbours,WPotential,EPotential,Geomat,OccupiedVec,InterMat,PotentialVec, &
                      PolarizationVec,TotInteraction,TotPol,TotRan,Tower_Mat,Tower_Prob_Vec,Gamma_T,Polarization)

      ! Calculamos además la energía inicial !
      TotalEnergy = TotInteraction + TotRan - TotPol

      ! Realizamos los pasos de equilibrio !
      do MonteCount = 1, SSTeps
              call MonteCarlo()
      end do

      ! Para calcular los Snapshots real(8)izaremos dos Do !
      ! así ahorramos comprobaciones                    !
      if (Snap == 0) then
              NSnap = 1
      else if (mod(NSteps,Snap) .ne. 0) then
              print *, 'ERROR: Snapshot frequency not compatible with the number of MonteCarlo steps'
              STOP
      else
              NSnap = NSteps/Snap
              NSteps = NSteps/NSnap
      end if
      
              
      ! Inicializamos el tiempo antes de empezar !
      Time = 0
      open(unit=25, file='traj.xyz',form='formatted',status='unknown')
      open(unit=26, file='ener.dat',form='formatted',status='unknown')
      do SnapCount = 1, NSnap
              do MonteCount = 1, NSteps
                      Time = Time + (-log(grnd()))/Gamma_T
                      call MonteCarlo()
                      write(26,*) Time, TotInteraction, TotRan, TotPol, TotalEnergy, Polarization
              end do
                 
              ! Definimos ahora un archivo .xyz para poder visualizar la simulación !
              ! Los siguientes writes determinana el formato de este archivo !
              write(25,*) N_Sites
              ! Esta es la linea de comentarios, dónde aprovecharemos para poner el tiempo !
              ! y la energía en ese paso. !
              write(25,*) 'Time= ', Time, 'TotEnergy= ', TotalEnergy
              do kk=1, N_Sites
                      ! Utilizamos dos átomos para poder representarlo en .xyz !
                      if (OccupiedVec(kk) == 1) then
                              Simulation = 'O'
                      else
                              Simulation = 'N'
                      end if
                      write(25,*) Simulation, Index_To_XY(L_Box,kk) ,1.
              end do
      end do
      close(25)
      close(26)
      
        contains
        subroutine Initialize(L_Box,N_Sites,Neighbours,WPotential,EPotential,Geomat,OccupiedVec,InterMat,PotentialVec, &
                        PolarizationVec,TotInteraction,TotPol,TotRan,Tower_Mat,Tower_Prob_Vec,Gamma_T,Polarization)
                implicit none
                integer, intent(out) :: L_Box, Neighbours
                real(8), intent(in) :: WPotential, EPotential
                integer, intent(in) :: N_Sites
                integer, dimension(:,:), allocatable, intent(out) :: Geomat, Tower_Mat
                integer, dimension(:), allocatable, intent(out) :: OccupiedVec
                real(8), intent(out) :: TotInteraction, TotRan, TotPol, Gamma_T, Polarization
                real(8), dimension(:), allocatable, intent(out) :: InterMat, PotentialVec, PolarizationVec, Tower_Prob_Vec

                ! Generamos la geometría 2D del sistema !
                call d2geom(L_Box,Geomat,Neighbours)
                ! Inicializamos los sitios ocupados !
                call SitesInit(L_Box,OccupiedVec)
                ! Inicializamos el potencial aleatorio !
                call RanPotentialInit(N_Sites,WPotential,OccupiedVec,PotentialVec,TotRan)
                ! Reservamos la memoria para algunas matrices !
                allocate(PolarizationVec(N_Sites))
                allocate(InterMat(N_Sites))
                ! Calculamos el valor de la polarización de cada posición respecto a E !
                call PolarizationCalc(L_Box,N_Sites,OccupiedVec,EPotential,PolarizationVec,TotPol,Polarization)
                ! Calculamos el valor de la energía de las interacciones inicial !
                call InteractionPotentialCalc(N_Sites,Neighbours,OccupiedVec,Geomat,InterMat,TotInteraction)
                ! Creamos la torre !
                call Create_Tower(L_Box,Gamma_T,Tower_Sites,Tower_Mat,Tower_Prob_Vec)
        end subroutine
        subroutine MonteCarlo()
                ! Seleccionamos un hueco ocupado !
                From = grnd() * N_Sites
                if (From == 0) then
                        From = N_Sites
                end if
                do while (OccupiedVec(From).eq.0)
                        From = grnd() * N_Sites
                        if (From == 0) then
                                From = N_Sites
                        end if
                end do 

                ! Escogemos una nueva posición To !
                call Tower(L_Box,Tower_Sites,Tower_Mat,Tower_Prob_Vec,From,To,Real_To)
              
                ! Si esta posición no está ocupada entonces llamamos a Metropolis !
                if (OccupiedVec(To) == 0) then
                        call Metropolis(N_Sites,From,To,Temperature,Neighbours,Geomat,OccupiedVec,PotentialVec,PolarizationVec, & 
                                  InterMat,TotInteraction,TotRan,TotPol,TotalEnergy,Change,EPotential,Real_To)
                        if (Change) then
                                Cartesian(:) = Index_To_XY(L_Box,From)
                                Polarization = Polarization + (Real_To - Cartesian(1))
                        end if
                end if
        end subroutine
        subroutine ReadFile(L_Box,EPotential,WPotential,Temperature,SSteps,NSteps,Snap,Seed)
                implicit none
                integer, intent(out) :: L_Box, Seed, Snap
                integer(8), intent(out) :: SSteps, NSteps
                real(8), intent(out) :: EPotential, WPotential, Temperature

                open(unit=24, file='param.dat',form='formatted',status='old')
                read(24,*) L_Box
                read(24,*) EPotential
                read(24,*) WPotential
                read(24,*) Temperature
                read(24,*) SSteps
                read(24,*) NSteps
                read(24,*) Snap
                read(24,*) Seed
                close(24)
        end subroutine
end program

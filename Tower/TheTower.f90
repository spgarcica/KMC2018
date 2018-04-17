module TheTower
        use Translate
        use mtmod
        use PBCmod
        contains
                subroutine Tower(L_Box,Tower_Sites,Tower_Mat,Tower_Prob_Vec,NIndex,Return_Index,Real_X)
                        implicit none
                        integer, intent(in) :: L_Box, Tower_Sites, NIndex
                        integer, dimension(Tower_Sites,2), intent(in) :: Tower_Mat
                        real(8), dimension(Tower_Sites), intent(in) :: Tower_Prob_Vec
                        integer, intent(out) :: Return_Index, Real_X
                        integer :: ii, jj, TowerStep, Case_Sel, IntegerRan
                        integer, dimension(2) :: CartesianCoords
                        real(8) :: RanParam
                        logical :: CheckTower

                        ! Definimos las cartesian para facilitar el algoritmo !
                        CartesianCoords = Index_To_XY(L_Box,NIndex)
                        ! Con un número al azar decidimos la distancia recorrida !
                        RanParam = grnd()
                        CheckTower = .true.
                        TowerStep = 1
                        ! Escalando la torre !
                        do while (CheckTower)
                                if (RanParam <= Tower_Prob_Vec(TowerStep)) then
                                        CheckTower = .false.
                                else
                                        TowerStep = TowerStep + 1
                                end if
                        end do

                        ! Comprobamos si hay degeneración vertical y horizontal !
                        if (Tower_Mat(TowerStep,2) == 1) then
                                IntegerRan = 4
                                ! Evitamos el valor de 4 en RanParam !
                                do while (IntegerRan == 4)
                                        IntegerRan = int(grnd()*4)
                                end do
                                select case(IntegerRan)
                                        case(0)
                                                CartesianCoords(1) = CartesianCoords(1) + Tower_Mat(TowerStep,1) - 1
                                        case(1)
                                                CartesianCoords(2) = CartesianCoords(2) - Tower_Mat(TowerStep,1) + 1
                                        case(2)
                                                CartesianCoords(1) = CartesianCoords(1) - Tower_Mat(TowerStep,1) + 1
                                        case(3) 
                                                CartesianCoords(2) = CartesianCoords(2) + Tower_Mat(TowerStep,1) - 1
                                end select
                        ! Comprobamos si hay degeneración diagonal !
                        else if (Tower_Mat(TowerStep,1) == Tower_Mat(TowerStep,2)) then
                                IntegerRan = 4
                                ! Evitamos el valor de 4 en RanParam !
                                do while (IntegerRan == 4)
                                        IntegerRan = int(grnd()*4)
                                end do
                                !print *, IntegerRan
                                select case(IntegerRan)
                                        case(0)
                                                CartesianCoords(:) = CartesianCoords(:) + Tower_Mat(TowerStep,:) - 1
                                        case(1)
                                                CartesianCoords(1) = CartesianCoords(1) + Tower_Mat(TowerStep,1) - 1
                                                CartesianCoords(2) = CartesianCoords(2) - Tower_Mat(TowerStep,2) + 1
                                        case(2)
                                                CartesianCoords(:) = CartesianCoords(:) - Tower_Mat(TowerStep,:) + 1
                                        case(3) 
                                                CartesianCoords(1) = CartesianCoords(1) - Tower_Mat(TowerStep,1) + 1
                                                CartesianCoords(2) = CartesianCoords(2) + Tower_Mat(TowerStep,2) - 1
                                end select
                        ! En caso de que no exista degeneración !
                        else
                                IntegerRan = 8
                                ! Evitamos el valor de 8 en RanParam !
                                do while (IntegerRan == 8)
                                        IntegerRan = int(grnd()*8)
                                end do
                                select case(IntegerRan)
                                        case(0)
                                                CartesianCoords(:) = CartesianCoords(:) + Tower_Mat(TowerStep,:) - 1
                                        case(1)
                                                CartesianCoords(1) = CartesianCoords(1) + Tower_Mat(TowerStep,1) - 1
                                                CartesianCoords(2) = CartesianCoords(2) - Tower_Mat(TowerStep,2) + 1
                                        case(2)
                                                CartesianCoords(1) = CartesianCoords(1) + Tower_Mat(TowerStep,2) - 1
                                                CartesianCoords(2) = CartesianCoords(2) - Tower_Mat(TowerStep,1) + 1
                                        case(3)
                                                CartesianCoords(1) = CartesianCoords(1) - Tower_Mat(TowerStep,2) + 1
                                                CartesianCoords(2) = CartesianCoords(2) - Tower_Mat(TowerStep,1) + 1
                                        case(4)
                                                CartesianCoords(1) = CartesianCoords(1) - Tower_Mat(TowerStep,1) + 1
                                                CartesianCoords(2) = CartesianCoords(2) - Tower_Mat(TowerStep,2) + 1
                                        case(5) 
                                                CartesianCoords(1) = CartesianCoords(1) - Tower_Mat(TowerStep,1) + 1
                                                CartesianCoords(2) = CartesianCoords(2) + Tower_Mat(TowerStep,2) - 1
                                        case(6)
                                                CartesianCoords(1) = CartesianCoords(1) - Tower_Mat(TowerStep,2) + 1
                                                CartesianCoords(2) = CartesianCoords(2) + Tower_Mat(TowerStep,1) - 1
                                        case(7)
                                                CartesianCoords(1) = CartesianCoords(1) + Tower_Mat(TowerStep,2) - 1
                                                CartesianCoords(2) = CartesianCoords(2) + Tower_Mat(TowerStep,1) - 1
                                end select
                        end if
                        Real_X = CartesianCoords(1)
                        call PBC(L_Box,CartesianCoords)
                        Return_Index = XY_To_Index(L_Box,CartesianCoords(1),CartesianCoords(2))

                end subroutine
                ! Esta subrutina crea dos arrays que asocian un cambio de posición a una probabilidad de movimiento !
                subroutine Create_Tower(L_Box,Gamma_T,Tower_Sites,Tower_Mat,Tower_Prob_Vec)
                        implicit none
                        integer, intent(in) :: L_Box
                        integer, intent(out) :: Tower_Sites
                        integer, dimension(:,:), allocatable, intent(out) :: Tower_Mat
                        real(8), intent(out) :: Gamma_T
                        real(8), dimension(:), allocatable, intent(out) :: Tower_Prob_Vec
                        integer :: ii, jj, Half_Range, Accumulator
                        real(8) :: Pivot, Gamma_T_Intern
                        real(8), dimension(2) :: Pivot_Vec
                        logical :: Check_Sort

                        ! Se define la forma de la matriz !
                        Half_Range = L_Box/2
                        Accumulator = 0
                        do ii=2, Half_Range+1
                                Accumulator = Accumulator + ii
                        end do
                        Tower_Sites = Accumulator
                        allocate(Tower_Mat(Tower_Sites,2))
                        allocate(Tower_Prob_Vec(Tower_Sites))
                        
                        ! Se calculan los valores de la matriz !
                        Accumulator = 0
                        Gamma_T = 0
                        do ii=2, Half_Range+1 
                                do jj=1, ii
                                        Accumulator = Accumulator + 1
                                        Tower_Mat(Accumulator,1) = ii
                                        Tower_Mat(Accumulator,2) = jj
                                        Tower_Prob_Vec(Accumulator) = exp(-2.*(sqrt((real(ii)**2) + real(jj)**2)))

                                        Gamma_T_Intern = Gamma_T_Intern + Tower_Prob_Vec(Accumulator)
                                        ! Tenemos en cuenta la degeneración en gamma y que está repetida !
                                        ! Esta Gamma es la de salida !
                                        if (jj == 1 .or. ii == jj) then
                                                Gamma_T = Gamma_T + Tower_Prob_Vec(Accumulator)*4.
                                        else
                                                Gamma_T = Gamma_T + Tower_Prob_Vec(Accumulator)*8.
                                        end if
                                end do
                        end do 

                        ! Renormalizamos los valores de probabilidad !
                        do ii=1, Tower_Sites
                                Tower_Prob_Vec(ii) = Tower_Prob_Vec(ii) / Gamma_T_Intern
                        end do

                        ! Ordenamos la matriz !
                        Check_Sort = .true.
                        do while (Check_Sort)
                                Check_Sort = .false.
                                do ii=1, Tower_Sites-1
                                        if (Tower_Prob_Vec(ii+1) > Tower_Prob_Vec(ii)) then
                                                Pivot = Tower_Prob_Vec(ii)
                                                Tower_Prob_Vec(ii) = Tower_Prob_Vec(ii+1)
                                                Tower_Prob_Vec(ii+1) = Pivot

                                                Pivot_Vec(:) = Tower_Mat(ii,:)
                                                Tower_Mat(ii,:) = Tower_Mat(ii+1,:)
                                                Tower_Mat(ii+1,:) = Pivot_Vec(:)

                                                Check_Sort = .true.
                                        end if
                                end do
                        end do

                        ! Hacemos la suma !
                        do ii=2, Tower_Sites
                                Tower_Prob_Vec(ii) = Tower_Prob_Vec(ii) + Tower_Prob_Vec(ii-1)
                        end do
                end subroutine
end module

module PBCmod
        contains
                subroutine PBC(L_Box,CartesianCoord)
                        implicit none
                        integer, intent(in) :: L_Box
                        integer, dimension(2), intent(inout) :: CartesianCoord
                        integer :: Selector

                        do Selector=1, 2
                                if (CartesianCoord(Selector) > L_Box) then
                                        CartesianCoord(Selector) = CartesianCoord(Selector) - L_Box
                                else if (CartesianCoord(Selector) < 1) then
                                        CartesianCoord(Selector) = CartesianCoord(Selector) + L_Box
                                end if
                        end do
                end subroutine
end module

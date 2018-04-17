module gengeom

contains
        subroutine d3geom(L,geomat,neighbours)
        implicit none
        integer, intent(in) :: L
        integer, dimension(:,:), allocatable, intent(out) :: geomat
        integer, intent(out) :: neighbours
        integer, dimension(2,L) :: init
        integer :: i, iz, iy, ix, L2

        neighbours = 6
        allocate(geomat(neighbours,L**2))
        
        do i = 1, L
                init(1,i)=i-1
                init(2,i)=i+1
        end do

        init(1,1)=L
        init(2,L)=1
        i=0
        L2 = L**2

        do iz = 1, L
                do iy = 1, L
                        do ix = 1, L
                                i=i+1
                                geomat(1,i)=init(2,ix)+L*(iy-1)+L2*(iz-1)
                                geomat(2,i)=init(1,ix)+L*(iy-1)+L2*(iz-1)
                                geomat(3,i)=ix+L*(init(2,iy)-1)+L2*(iz-1)
                                geomat(4,i)=ix+L*(init(1,iy)-1)+L2*(iz-1)
                                geomat(5,i)=ix+L*(iy-1)+L2*(init(2,iz)-1)
                                geomat(6,i)=ix+L*(iy-1)+L2*(init(1,iz)-1)
                        end do
                end do
        end do
        end subroutine

        subroutine d2geom(L,geomat,neighbours)
        implicit none
        integer, intent(in) :: L
        integer, dimension(:,:),allocatable, intent(out) :: geomat
        integer, intent(out) :: neighbours
        integer, dimension(2,L) :: init
        integer :: i, iy, ix

        neighbours = 4
        allocate(geomat(neighbours,L**2))

        do i = 1, L
                init(1,i)=i-1
                init(2,i)=i+1
        end do

        init(1,1)=L
        init(2,L)=1
        i=0

        do iy = 1, L
                do ix = 1, L
                        i=i+1
                        geomat(1,i)=init(2,ix)+L*(iy-1)
                        geomat(2,i)=init(1,ix)+L*(iy-1)
                        geomat(3,i)=ix+L*(init(2,iy)-1)
                        geomat(4,i)=ix+L*(init(1,iy)-1)
                end do
        end do
        end subroutine
end module

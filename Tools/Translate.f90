module Translate
      contains
              integer function XY_To_Index(L_Box,xx,yy)
                      integer, intent(in) :: L_Box ,xx, yy

                      XY_To_Index = L_Box * (xx-1) + yy
              end function XY_To_Index

              function Index_To_XY(L_Box,NIndex)
                      integer, intent(in) :: NIndex
                      integer, dimension(2) :: Index_To_XY

                      if (mod(NIndex,L_Box) == 0) then
                              Index_To_XY(1) = NIndex/L_Box
                              Index_To_XY(2) = L_Box
                      else
                              Index_To_XY(1) = NIndex/L_Box + 1
                              Index_To_XY(2) = NIndex - (Index_To_XY(1) -1)*L_Box 
                      end if
              end function Index_To_XY
end module

module mexinterface
  use, intrinsic :: ISO_C_BINDING
    integer, parameter :: mwPointer = 8
    integer, parameter :: mwSize = 8
    interface

    integer*8 function mxgetpr(pm) bind(c, name = 'MXGETPR')
      integer*8 :: pm
    end function mxgetpr
    !
    integer*8 function mxgetpi(pm) bind(c, name = 'MXGETPI')
      integer*8 :: pm
    end function mxgetpi
    !
    integer*8 function mxgetm(pm) bind(c, name = 'MXGETM')
      integer*8 :: pm
    end function mxgetm
    !
    integer*8 function mxgetn(pm) bind(c, name = 'MXGETN')
      integer*8 pm
    end function mxgetn
    !
    integer*8 function mxcreatedoublematrix(m,n,type) bind(c, name = 'MXCREATEDOUBLEMATRIX')
      integer*8 m, n
      integer*4 type
    end function mxcreatedoublematrix
    !
    subroutine mexerrmsgtxt(s) bind(C, name = 'mexErrMsgTxt')
      import c_char
      character(c_char) s(*)
    end subroutine mexerrmsgtxt
    !
    integer*4 function mexprintf(s) bind(C, name = 'mexPrintf')
      import c_char
      character(c_char) s(*)
    end function mexprintf
    !
 logical function mxiscomplex(pm) bind(C, name = 'MXISCOMPLEX')
      integer*8 :: pm
    end function mxiscomplex
SUBROUTINE  mxcopyptrtocomplex16( x_r, x_i, x, n) bind(C, name = 'MXCOPYPTRTOCOMPLEX16')
      integer*8 :: n
      integer*8::  x_i,  x_r
      COMPLEX*16:: x

   END SUBROUTINE mxCopyPtrToComplex16
SUBROUTINE  mxsetpi( pm,pi) bind(C, name = 'MXSETPI')
      integer*8 :: pm
integer*8 :: pi

   END SUBROUTINE mxsetpi
SUBROUTINE  mxcopyreal8toptr(y, px, n) bind(C, name = 'MXCOPYREAL8TOPTR')

      integer*8 :: n
      integer*8::  px
      real*8 :: y

   END SUBROUTINE mxcopyreal8toptr
   SUBROUTINE  mxcopycomplex16toptr( x, x_r, x_i,  n) bind(C, name = 'MXCOPYCOMPLEX16TOPTR')

      integer*8 :: n
      integer*8 ::  x_i,  x_r
      COMPLEX*16 :: x

   END SUBROUTINE mxcopycomplex16toptr

   SUBROUTINE  mxcopyptrtoreal8( x_r, x, n) bind(C, name = 'MXCOPYPTRTOREAL8')

      integer*8 :: n
      integer*8 ::   x_r
      real*8:: x

   END SUBROUTINE mxcopyptrtoreal8

  end interface
end module mexinterface

! EISCOR_ROOTS
!
! INSTALL INSTRUCTIONS
!
!   To compile this mex file you need to install a copy of eiscor before. 
!   You can download and install eiscor with the following commands: 
!   
!    $ git clone git://github.com/eiscor/eiscor
!    $ cd eiscor && make && make install
!
!   If the compilation and installation succeeds, then you can compile this
!   mex file running the command 
! 
!    >> mex eiscor_roots.F90 ~/eiscor/lib/libeiscor.so.0.2.0
! 
!   at the MATLAB prompt. 
!
! LICENSE
! 
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include "fintrf.h"

subroutine mexFunction(nlhs, plhs, nrhs, prhs)	
  implicit none

  ! Argument declarations
  integer*8 :: nlhs, nrhs
  mwPointer :: plhs(*), prhs(*)

  ! Local variables
  integer :: degree, i, info
  mwSize :: mxGetNumberOfElements
  mwPointer :: mxGetPr, mxGetPi, pr, pi
  mwPointer :: mxCreateDoubleMatrix
  complex*16, allocatable :: p(:), r(:), res(:), qq(:,:), zz(:,:), ss(:,:), tt(:,:)
  logical :: qz, wantv

  ! The input should be a vector containing the coefficients
  ! of a polynomial, ordered according to the MATLAB use.
  degree = mxGetNumberOfElements(prhs(1)) - 1
  allocate(p(degree+1), r(degree), res(degree))

  ! Eigenvectors are required only if the user requests the Schur form
  wantv = (nlhs > 1)

  if (nrhs .gt. 1) then
     pr = mxGetPr(prhs(2))
     call mxCopyPtrToReal8(pr, p, 1)
     qz = (real(p(1)) .ne. 0.d0)
  else
     qz = .false.
  end if

  pr = mxGetPr(prhs(1))
  pi = mxGetPi(prhs(1))
  call mxCopyPtrToComplex16(pr, pi, p, degree + 1)

  ! Create the storage for the Schur form and the unitary matrices
  if (wantv) then
     allocate(qq(degree,degree), zz(degree,degree), & 
          ss(degree,degree), tt(degree,degree))
  else
     ! Only the space for the eigenvalues is required
     allocate(ss(degree,1), tt(degree,1))
  end if

  ! Compute the roots
  if (qz) then
     call mat_z_poly_roots_qz(degree, p, r, info, ss, tt, qq, zz, wantv)
  else
     call mat_z_poly_roots_qr(degree, p, r, info, ss, qq, wantv)
  end if

  ! Create a complex output
  plhs(1) = mxCreateDoubleMatrix(degree, 1, 1)
  pr = mxGetPr(plhs(1))
  pi = mxGetPi(plhs(1))

  ! Store the result
  call mxCopyComplex16ToPtr(r, pr, pi, degree)

  ! Do the same for the matrices, if the user requested it
  if (wantv) then
     do i = 2, 5
        plhs(i) = mxCreateDoubleMatrix(degree, degree, 1)
     end do

     pr = mxGetPr(plhs(2))
     pi = mxGetPi(plhs(2))
     call mxCopyComplex16ToPtr(qq, pr, pi, degree * degree)

     pr = mxGetPr(plhs(3))
     pi = mxGetPi(plhs(3))
     if (.not. qz) then 
        zz = qq
     end if
     call mxCopyComplex16ToPtr(zz, pr, pi, degree * degree)

     pr = mxGetPr(plhs(4))
     pi = mxGetPi(plhs(4))
     call mxCopyComplex16ToPtr(ss, pr, pi, degree * degree)

     pr = mxGetPr(plhs(5))
     pi = mxGetPi(plhs(5))
     if (.not. qz) then
        tt = 0.d0;
        do i = 1, degree
           tt(i,i) = 1.d0
        end do
     end if
     call mxCopyComplex16ToPtr(tt, pr, pi, degree * degree)

     deallocate(qq, zz)
  end if

  deallocate(p, r, res, ss, tt)

end subroutine mexFunction

subroutine mat_z_poly_roots_qz(N,COEFFS,ROOTS,INFO,SS,TT,QQ,ZZ,wantv)

  implicit none

  ! input variables
  integer, intent(in) :: N
  integer, intent(inout) :: INFO
  complex(8), intent(in) :: COEFFS(N+1)
  complex(8), intent(inout) :: ROOTS(N)
  complex(8), intent(out) :: SS(N,N), TT(N,N), QQ(N,N), ZZ(N,N)
  logical :: qz, wantv

  ! compute variables
  integer :: ii
  real(8) :: scl
  logical, allocatable :: P(:)
  integer, allocatable :: ITS(:)
  real(8), allocatable :: Q(:),D1(:),C1(:),B1(:)
  real(8), allocatable :: D2(:),C2(:),B2(:)
  complex(8), allocatable :: V(:),W(:),T(:,:)
  interface
     function l_upr1fact_hess(m,flags)
       logical :: l_upr1fact_hess
       integer, intent(in) :: m
       logical, dimension(m-2), intent(in) :: flags
     end function l_upr1fact_hess
  end interface

  ! allocate memory
  allocate(P(N-2),ITS(N-1),Q(3*(N-1)),D1(2*(N+1)),C1(3*N),B1(3*N))   
  allocate(V(N),W(N),D2(2*(N+1)),C2(3*N),B2(3*N))    

  ! initialize INFO
  INFO = 0

  ! fill P
  P = .FALSE.

  ! fill V and W
  W = cmplx(0d0,0d0,kind=8)
  scl = maxval(abs(COEFFS))
  V(N) = ((-1d0)**(N))*COEFFS(N+1)/scl
  do ii=1,(N-1)
     V(ii) = -1.0d0 * COEFFS(N+1-ii)/scl
  end do
  W(N) = COEFFS(1)/scl

  ! factor companion matrix
  call z_comppen_compress(N,P,V,W,Q,D1,C1,B1,D2,C2,B2)

  ! call z_upr1fpen_qz
  call z_upr1fpen_qz(wantv,.true.,l_upr1fact_hess,N,P,Q,D1,C1,B1,D2,C2,B2,N,qq,zz,ITS,INFO)

  ! extract roots
  call z_upr1utri_decompress(.not. wantv,N,D1,C1,B1,ss)
  call z_upr1utri_decompress(.not. wantv,N,D2,C2,B2,tt)
  do ii=1,N
     if (wantv) then
        ROOTS(ii) = ss(ii,ii)/tt(ii,ii)
     else
        ROOTS(ii) = ss(ii,1) / tt(ii,1)
     end if
  end do

  ! free memory
  deallocate(P,ITS,Q,D1,C1,B1,D2,C2,B2,V,W)

end subroutine mat_z_poly_roots_qz

subroutine mat_z_poly_roots_qr(N,COEFFS,ROOTS,INFO,SS,QQ,wantv)

  implicit none

  ! input variables
  integer, intent(in) :: N
  integer, intent(inout) :: INFO
  complex(8), intent(in) :: COEFFS(N+1)
  complex(8), intent(inout) :: ROOTS(N)
  complex(8), intent(out) :: SS(N,N), QQ(N,N)
  logical :: wantv

  ! compute variables
  integer :: ii
  real(8) :: scl
  logical, allocatable :: P(:)
  integer, allocatable :: ITS(:)
  real(8), allocatable :: Q(:),D1(:),C1(:),B1(:)
  complex(8), allocatable :: V(:),T(:,:)
  interface
     function l_upr1fact_hess(m,flags)
       logical :: l_upr1fact_hess
       integer, intent(in) :: m
       logical, dimension(m-2), intent(in) :: flags
     end function l_upr1fact_hess
  end interface

  ! allocate memory
  allocate(P(N-2),ITS(N-1),Q(3*(N-1)),D1(2*(N+1)),C1(3*N),B1(3*N))   
  allocate(V(N))

  ! initialize INFO
  INFO = 0

  ! fill P
  P = .FALSE.

  ! fill V
  V(N) = ((-1d0)**(N))*COEFFS(N+1)/COEFFS(1)
  do ii=1,(N-1)
     V(ii) = -COEFFS(N+1-ii)/COEFFS(1)
  end do

  ! factor companion matrix
  call z_compmat_compress(N,P,V,Q,D1,C1,B1)

  ! call z_upr1fpen_qr
  call z_upr1fact_qr(wantv,.true.,l_upr1fact_hess,N,P,Q,D1,C1,B1,N,qq,ITS,INFO)

  ! extract roots
  call z_upr1utri_decompress(.not. wantv,N,D1,C1,B1,ss)
  do ii=1,N
     if (wantv) then
        ROOTS(ii) = ss(ii,ii)
     else
        ROOTS(ii) = ss(ii,1)
     end if
  end do

  ! free memory
  deallocate(P,ITS,Q,D1,C1,B1,V)

end subroutine mat_z_poly_roots_qr


!This file is part of Notus 0.4.0

!Copyright Bordeaux-INP, Université de Bordeaux, CNRS
!
!Contributors:
!Antoine Lemoine, 14-06-2017, antoine.lemoine@bordeaux-inp.fr

!This software is a computer program whose purpose is to simulate fluid flows.

!This software is governed by the CeCILL license under French law and
!abiding by the rules of distribution of free software.  You can  use,
!modify and/ or redistribute the software under the terms of the CeCILL
!license as circulated by CEA, CNRS and INRIA at the following URL
!"http://www.cecill.info".

!As a counterpart to the access to the source code and  rights to copy,
!modify and redistribute granted by the license, users are provided only
!with a limited warranty  and the software's author,  the holder of the
!economic rights,  and the successive licensors  have only  limited
!liability.

!In this respect, the user's attention is drawn to the risks associated
!with loading,  using,  modifying and/or developing or reproducing the
!software by the user in light of its specific status of free software,
!that may mean  that it is complicated to manipulate,  and  that  also
!therefore means  that it is reserved for developers  and  experienced
!professionals having in-depth computer knowledge. Users are therefore
!encouraged to load and test the software's suitability as regards their
!requirements in conditions enabling the security of their systems and/or
!data to be ensured and,  more generally, to use and operate it in the
!same conditions as regards security.

!The fact that you are presently reading this means that you have had
!knowledge of the CeCILL license and that you accept its terms.

! Contains the minimization algorithm to find the MOF reconstruction in 3D
module mod_mof3d_bfgs
   use mod_cg3_polyhedron
   use mod_mof3d_gradient
   implicit none

   private

   public :: mof3d_bfgs

contains

   !> Broyden-Fletcher-Goldfarb-Shanno method to minimize the objective function of MOF.
   !!
   !! Reference:
   !!  - R. Fletcher, Practical methods of optimization, John Wiley & Sons, 1987.
   !!
   !! @param[in]     polyhedron:    Convex polyhedron.
   !! @param[in]     ref_centroid1: Coordinates of the reference centroid of the material 1.
   !! @param[in]     ref_centroid2: Coordinates of the reference centroid of the material 2.
   !! @param[in]     ref_volume:    Reference volume of the material 1.
   !! @param[inout]  angles:        Initial guess and final angles.
   !! @param[out]    normal:        Unit normal of the reconstructed interface (material 1 -> material 2).
   !! @param[out]    stat:          Number of calls of the gradient. stat(1): calls in BFGS, stat(2): calls in line-search.
   !! @param[out]    residual:      Residuals. residual(1): derivative, residual(2): objective function.
   !! @ingroup moment_of_fluid
   pure subroutine mof3d_bfgs(polyhedron, ref_centroid1, ref_centroid2, ref_volume, angles, normal, stat, residual)
      use mod_cg3_flood_polyhedron
      use variables_mof, only : mof3d_max_iter, mof3d_tol_derivative
      type(t_polyhedron), intent(in) :: polyhedron
      double precision, dimension(3), intent(in) :: ref_centroid1, ref_centroid2
      double precision, intent(in) :: ref_volume
      double precision, dimension(2), intent(inout) :: angles
      double precision, dimension(3), intent(out) :: normal
      integer, dimension(2), intent(out) :: stat
      double precision, dimension(2), intent(out) :: residual

      integer :: i
      double precision :: cell_volume, determinant, objective, objective_next, delta_objective, dot_delta, dot_hessian
      double precision, dimension(2) :: angles_next, gradient, gradient_next, direction, step, delta_gradient
      double precision, dimension(3) :: cell_centroid
      double precision, dimension(2,2) :: hessian, hessian_next

      double precision, parameter :: step_max = 2d0*acos(0d0)/40d0 ! pi/40

      ! Compute the volume and the centroid of the polyhedron
      call cg3_polyhedron_compute_centroid(polyhedron, cell_volume, cell_centroid)

      ! Initialize the statistics
      stat = 0
      residual = huge(residual)

      ! Initialize the Hessian approximation to the identity
      hessian(1,1) = 1d0
      hessian(1,2) = 0d0
      hessian(2,1) = 0d0
      hessian(2,2) = 1d0

      ! Compute the initial gradient
      call mof3d_compute_gradient(polyhedron, cell_volume, cell_centroid, angles, ref_volume, &
         &                        ref_centroid1, ref_centroid2, objective, gradient           )
      residual(2) = objective

      ! Update statistics
      stat(1) = stat(1) + 1

      delta_objective = step_max

      do i = 1, mof3d_max_iter
         residual(1) = norm2(gradient)

         ! Check for convergence
         if (residual(1) < mof3d_tol_derivative) exit

         ! Compute the opposite of the product of the Hessian approximation by the gradient
         determinant = hessian(1,1)*hessian(2,2) - hessian(1,2)*hessian(2,1)

         ! Singular Hessian, no more convergence
         if (abs(determinant) <= tiny(determinant)) exit

         direction = [hessian(1,2)*gradient(2) - hessian(2,2)*gradient(1), &
            &         hessian(2,1)*gradient(1) - hessian(1,1)*gradient(2)]/determinant

         ! Compute the step from the direction (linear search)
         call mof3d_line_search(polyhedron, cell_volume, cell_centroid, ref_centroid1, ref_centroid2, ref_volume, angles, &
            &                   direction, objective, gradient, delta_objective, mof3d_tol_derivative, step, stat         )

         ! Update the angles
         angles_next = angles + step

         ! Compute the next reconstruction and gradient
         call mof3d_compute_gradient(polyhedron, cell_volume, cell_centroid, angles_next, ref_volume, &
            &                        ref_centroid1, ref_centroid2, objective_next, gradient_next      )

         ! Compute the difference between the next gradient and the previous gradient
         delta_gradient = gradient_next - gradient

         ! Update statistics
         stat(1) = stat(1) + 1

         ! Compute dot product between the gradient increment and the step
         dot_delta = dot_product(delta_gradient, step)

         ! Check for potential singularity
         if (abs(dot_delta) <= tiny(dot_delta)) then
            angles = angles_next
            residual(2) = objective_next
            exit
         end if

         ! Compute dot product between pseudo direction and step
         dot_hessian = dot_product(matmul(hessian, step), step)

         ! Check for potential singularity
         if (abs(dot_hessian) <= tiny(dot_hessian)) then
            angles = angles_next
            residual(2) = objective_next
            exit
         end if

         ! Update the Hessian approximation
         hessian_next = hessian + mof3d_kronecker(delta_gradient, delta_gradient)/dot_delta &
            &         - matmul(matmul(hessian, mof3d_kronecker(step, step)), hessian)/dot_hessian

         delta_objective = objective - objective_next

         if (abs(delta_objective)/max(abs(objective), abs(objective_next), 1d0) <= epsilon(delta_objective)) then
            angles = angles_next
            exit
         end if

         delta_objective = min(delta_objective, step_max)

         angles      = angles_next
         gradient    = gradient_next
         hessian     = hessian_next
         objective   = objective_next
         residual(2) = objective_next
      end do

      normal = [cos(angles(1))*sin(angles(2)), &
         &      sin(angles(1))*sin(angles(2)), &
         &      cos(angles(2))]
   end subroutine mof3d_bfgs

   pure function mof3d_kronecker(u, v) result(m)
      double precision, dimension(2), intent(in) :: u, v
      double precision, dimension(2,2) :: m

      m(1,1) = u(1)*v(1)
      m(1,2) = u(1)*v(2)
      m(2,1) = u(2)*v(1)
      m(2,2) = u(2)*v(2)
   end function mof3d_kronecker

   pure subroutine mof3d_minimize_cubic(polyhedron, cell_volume, cell_centroid, ref_centroid1, ref_centroid2, &
      &                                 ref_volume, angles, direction, a, b, alpha, stat)
      type(t_polyhedron), intent(in) :: polyhedron
      double precision, intent(in) :: cell_volume
      double precision, dimension(3), intent(in) :: cell_centroid
      double precision, dimension(3), intent(in) :: ref_centroid1, ref_centroid2
      double precision, intent(in) :: ref_volume
      double precision, dimension(2), intent(in) :: angles, direction
      double precision, intent(in) :: a, b
      double precision, intent(out) :: alpha
      integer, dimension(2), intent(inout) :: stat

      double precision, dimension(4) :: cubic
      double precision, dimension(2) :: gradient
      double precision :: fa, fb, fpa, fpb
      double precision :: cubic3, cubic4, delta, z3, z4, z_min, c_min

      ! Evaluate the objective function at the beginning of the interval
      call mof3d_compute_gradient(polyhedron, cell_volume, cell_centroid, angles + a*direction, &
         &                        ref_volume, ref_centroid1, ref_centroid2, fa, gradient        )
      fpa = (b - a)*dot_product(gradient, direction)

      ! Evaluate the derivative at the beginning of the interval
      call mof3d_compute_gradient(polyhedron, cell_volume, cell_centroid, angles + b*direction, &
         &                        ref_volume, ref_centroid1, ref_centroid2, fb, gradient        )
      fpb = (b - a)*dot_product(gradient, direction)

      ! Update the statistics
      stat(2) = stat(2) + 2

      ! Extreme points of the interval
      if (fb < fa) then
         c_min = fb
         z_min = 1d0
      else
         c_min = fa
         z_min = 0d0
      end if

      cubic(1) = fa
      cubic(2) = fpa
      cubic(3) = 3d0*(fb - fa) - 2d0*fpa - fpb
      cubic(4) = fpa + fpb - 2d0*(fb - fa)

      delta = cubic(3)**2 - 3d0*cubic(2)*cubic(4)

      if (abs(cubic(4)) >= tiny(delta)) then
         ! Cubic approximation
         if (delta >= 0d0) then
            z3 = (-cubic(3) - sqrt(delta))/(3d0*cubic(4))
            z4 = (-cubic(3) + sqrt(delta))/(3d0*cubic(4))

            cubic3 = cubic(1) + z3*(cubic(2) + z3*(cubic(3) + z3*cubic(4)))
            cubic4 = cubic(1) + z4*(cubic(2) + z4*(cubic(3) + z4*cubic(4)))

            if (cubic3 < c_min) then
               if (z3 > 0d0 .and. z3 < 1d0) then
                  c_min = cubic3
                  z_min = z3
               end if
            end if

            if (cubic4 < c_min .and. z4 > 0d0 .and. z4 < 1d0) z_min = z4
         end if
      else
         if (abs(cubic(3)) >= tiny(delta)) then
            ! Parabolic approximation
            z3 = -cubic(2)/(2d0*cubic(3))

            cubic3 = cubic(1) + z3*(cubic(2) + z3*cubic(3))

            if (cubic3 < c_min .and. z3 > 0d0 .and. z3 < 1d0) z_min = z3
         end if
      end if

      alpha = a + z_min*(b - a)
   end subroutine mof3d_minimize_cubic

   pure subroutine mof3d_line_search(polyhedron, cell_volume, cell_centroid, ref_centroid1, ref_centroid2, ref_volume, angles, &
      &                              direction, objective0, gradient0, step_max, tolerance, step, stat)
      type(t_polyhedron), intent(in) :: polyhedron
      double precision, intent(in) :: cell_volume
      double precision, dimension(3), intent(in) :: cell_centroid
      double precision, dimension(3), intent(in) :: ref_centroid1, ref_centroid2
      double precision, intent(in) :: ref_volume
      double precision, dimension(2), intent(in) :: angles, direction
      double precision, intent(in) :: objective0
      double precision, dimension(2), intent(in) :: gradient0
      double precision, intent(in) :: step_max, tolerance
      double precision, dimension(2), intent(out) :: step
      integer, dimension(2), intent(inout) :: stat

      double precision, parameter :: rho = 0.25d0
      double precision, parameter :: sigma = 0.5d0
      double precision, parameter :: tau1 = 3d0
      double precision, parameter :: tau2 = 0.1d0
      double precision, parameter :: tau3 = 0.5d0
      double precision, parameter :: f_bar = 0d0 ! The objective function >= 0
      integer, parameter :: MAX_ITER = 10

      double precision, dimension(2) :: gradient
      double precision :: alpha, alpha_next, alpha_next_next, mu, f0, fp0, f, f_next, fp, fp_next, a, b, fa, fpa, a_j
      integer :: i

      f0 = objective0
      fp0 = dot_product(gradient0, direction)

      mu = (f_bar - f0)/(rho*fp0)

      ! Bracketting: start with α ∈ [0,μ]
      alpha = 0d0
      alpha_next = min(-2d0*max(step_max, 10d0*tolerance)/fp0, 1d0)

      f = f0
      fp = fp0
      fpa = fp0

      ! Bracket the minimum
      i = 1
      do
         step = alpha_next*direction

         ! Evaluate the gradient at angles + alpha_next*direction
         call mof3d_compute_gradient(polyhedron, cell_volume, cell_centroid, angles + step,      &
            &                        ref_volume, ref_centroid1, ref_centroid2, f_next, gradient)
         fp_next = dot_product(gradient, direction)

         ! Update the statistics
         stat(2) = stat(2) + 1

         if (f_next <= f_bar) return

         !if (f_next > f0 + alpha_next*fp0 .or. f_next >= f) then
         if (f_next >= f) then
            a = alpha
            b = alpha_next
            a_j = a
            fa = f
            fpa = fp
            exit
         end if

         if (abs(fp_next) <= -sigma*fp0) return

         ! Prevent infinite loop
         if (i > MAX_ITER) return

         if (fp_next >= 0d0) then
            a = alpha_next
            b = alpha
            a_j = a
            fa = f_next
            fpa = fp_next
            exit
         end if

         if (mu <= 2d0*alpha_next - alpha) then
            alpha = alpha_next
            alpha_next = mu
         else
            call mof3d_minimize_cubic(polyhedron, cell_volume, cell_centroid, ref_centroid1, ref_centroid2, ref_volume, angles, &
               &                      direction, 2d0*alpha_next - alpha, min(mu, alpha_next + tau1*(alpha_next - alpha)),       &
               &                      alpha_next_next, stat)
            alpha = alpha_next
            alpha_next = alpha_next_next
         end if

         fp = fp_next
         f = f_next
         i = i + 1
      end do

      ! Shorten the bracket
      do i = 1, MAX_ITER
         call mof3d_minimize_cubic(polyhedron, cell_volume, cell_centroid, ref_centroid1, ref_centroid2, ref_volume, &
            &                      angles, direction, a + tau2*(b - a), b - tau3*(b - a), alpha, stat)
         step = alpha*direction

         if ((a_j - alpha)*fpa <= tolerance) return

         ! Evaluate the objective function at angles + step
         call mof3d_compute_gradient(polyhedron, cell_volume, cell_centroid, angles + step, &
            &                        ref_volume, ref_centroid1, ref_centroid2, f, gradient  )

         ! Update the statistics
         stat(2) = stat(2) + 1

         if (f > f0 + rho*alpha*fp0 .or. f >= fa) then
            b = alpha
         else
            ! Evaluate the gradient at angles + step
            fp = dot_product(gradient, direction)

            if (abs(fp) <= -sigma*fp0) return

            if ((b - a)*fp >= 0d0) b = a

            a = alpha
         end if

         if (abs(b - a) < tolerance) exit
      end do
   end subroutine mof3d_line_search

end module mod_mof3d_bfgs

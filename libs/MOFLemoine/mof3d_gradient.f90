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

module mod_mof3d_gradient
   implicit none

   private

   public :: mof3d_compute_objective_function, mof3d_compute_gradient, mof3d_debug_geometric_gradient

contains

   pure subroutine mof3d_compute_objective_function(centroid, ref_volume, ref_centroid1, ref_centroid2, &
      &                                             cell_volume, cell_centroid, objective               )
      use variables_mof, only: mof_use_symmetric_reconstruction
      double precision, dimension(3), intent(in) :: centroid
      double precision, intent(in) :: ref_volume
      double precision, dimension(3), intent(in) :: ref_centroid1, ref_centroid2
      double precision, intent(in) :: cell_volume
      double precision, dimension(3), intent(in) :: cell_centroid
      double precision, intent(out) :: objective

      double precision, dimension(3) :: dual_centroid, diff

      objective = dot_product(centroid - ref_centroid1, centroid - ref_centroid1)

      if (mof_use_symmetric_reconstruction) then
         ! Compute the dual centroid
         dual_centroid = (cell_volume*cell_centroid - ref_volume*centroid)/(cell_volume - ref_volume)

         ! Add the contribution of the symmetric part
         diff = dual_centroid - ref_centroid2
         objective = objective + dot_product(diff, diff)
      end if
   end subroutine mof3d_compute_objective_function

   !> Compute the gradient of the objective function of MOF.
   !!
   !! @param[in]     polyhedron:    Convex polyhedron.
   !! @param[in]     cell_volume:   Volume of the cell.
   !! @param[in]     cell_centroid: Coordinates of the centroid of the cell.
   !! @param[in]     angles:        Spherical angles.
   !! @param[in]     ref_volume:    Reference volume of the material 1.
   !! @param[in]     ref_centroid1: Coordinates of the reference centroid of the material 1.
   !! @param[in]     ref_centroid2: Coordinates of the reference centroid of the material 2.
   !! @param[out]    objective:     Objective function.
   !! @param[out]    gradient:      Coordinates of the gradient.
   !! @ingroup moment_of_fluid
   pure subroutine mof3d_compute_gradient(polyhedron, cell_volume, cell_centroid, angles, ref_volume, &
      &                                   ref_centroid1, ref_centroid2, objective, gradient           )
      use mod_cg3_complete_polyhedron_structure
      use mod_cg3_flood_polyhedron
      use mod_cg3_polyhedron
      use mod_mof3d_analytic_centroid
      use variables_mof, only: mof3d_internal_is_analytic_gradient_enabled
      use variables_mof, only: mof_use_symmetric_reconstruction
      type(t_polyhedron), intent(in) :: polyhedron
      double precision, intent(in) :: cell_volume
      double precision, dimension(3), intent(in) :: cell_centroid
      double precision, dimension(2), intent(in) :: angles
      double precision, intent(in) :: ref_volume
      double precision, dimension(3), intent(in) :: ref_centroid1, ref_centroid2
      double precision, intent(out) :: objective
      double precision, dimension(2), intent(out) :: gradient

      double precision, dimension(3) :: c

      if (mof3d_internal_is_analytic_gradient_enabled) then
         c = polyhedron%point(:,8)

         if (mof_use_symmetric_reconstruction) then
            call mof3d_compute_analytic_gradient_symmetric(angles, ref_centroid1, ref_centroid2, ref_volume, c, objective, gradient)
         else
            call mof3d_compute_analytic_gradient(angles, ref_centroid1, ref_volume, c, objective, gradient)
         end if
      else
         call mof3d_compute_gradient_geometric(polyhedron, cell_volume, cell_centroid, angles, ref_volume, &
            &                                  ref_centroid1, ref_centroid2, objective, gradient)
      end if
   end subroutine mof3d_compute_gradient

   !> Compute the gradient of the objective function of MOF using geometric approaches.
   !!
   !! @param[in]     polyhedron:    Convex polyhedron.
   !! @param[in]     cell_volume:   Volume of the cell.
   !! @param[in]     cell_centroid: Coordinates of the centroid of the cell.
   !! @param[in]     angles:        Spherical angles.
   !! @param[in]     ref_volume:    Reference volume of the material 1.
   !! @param[in]     ref_centroid1: Coordinates of the reference centroid of the material 1.
   !! @param[in]     ref_centroid2: Coordinates of the reference centroid of the material 2.
   !! @param[out]    objective:     Objective function.
   !! @param[out]    gradient:      Coordinates of the gradient.
   !! @ingroup moment_of_fluid
   pure subroutine mof3d_compute_gradient_geometric(polyhedron, cell_volume, cell_centroid, angles, ref_volume, &
      &                                             ref_centroid1, ref_centroid2, objective, gradient           )
      use mod_cg3_complete_polyhedron_structure
      use mod_cg3_flood_polyhedron
      use mod_cg3_polyhedron
      use mod_mof3d_analytic_centroid
      use variables_mof, only: mof3d_use_analytical_derivatives, mof3d_use_optimized_centroid
      type(t_polyhedron), intent(in) :: polyhedron
      double precision, intent(in) :: cell_volume
      double precision, dimension(3), intent(in) :: cell_centroid
      double precision, dimension(2), intent(in) :: angles
      double precision, intent(in) :: ref_volume
      double precision, dimension(3), intent(in) :: ref_centroid1, ref_centroid2
      double precision, intent(out) :: objective
      double precision, dimension(2), intent(out) :: gradient

      type(t_polyhedron) :: polyhedron_full
      type(t_chained_polygon) :: polygon
      double precision, dimension(3) :: direction, centroid
      double precision :: volume, sin_phi
      integer :: error_id

      ! Geometric methods
      sin_phi = sin(angles(2))

      ! Central reconstruction
      direction = [cos(angles(1))*sin_phi, &
         &         sin(angles(1))*sin_phi, &
         &         cos(angles(2))]

      ! Compute the centroid
      if (mof3d_use_optimized_centroid) then
         call cg3_flood_polyhedron_centroid(polyhedron, direction, ref_volume, centroid, polygon)
      else if (mof3d_use_analytical_derivatives) then
         call cg3_flood_polyhedron(polyhedron, direction, ref_volume, polyhedron_full, polygon=polygon)
         call cg3_complete_polyhedron_structure(polyhedron_full, error_id)
         call cg3_polyhedron_compute_centroid(polyhedron_full, volume, centroid)
      else
         call cg3_flood_polyhedron(polyhedron, direction, ref_volume, polyhedron_full)
         call cg3_complete_polyhedron_structure(polyhedron_full, error_id)
         call cg3_polyhedron_compute_centroid(polyhedron_full, volume, centroid)
      end if

      ! Compute the derivatives
      if (mof3d_use_optimized_centroid) then
         call mof3d_compute_geometric_gradient(polygon, cell_volume, cell_centroid, angles, ref_volume,   &
            &                                  centroid, ref_centroid1, ref_centroid2, objective, gradient)
      else if (mof3d_use_analytical_derivatives) then
         call mof3d_compute_geometric_gradient(polygon, cell_volume, cell_centroid, angles, ref_volume,   &
            &                                  centroid, ref_centroid1, ref_centroid2, objective, gradient)
      else
         call mof3d_compute_gradient_finite_differences(polyhedron, cell_volume, cell_centroid, angles, ref_volume, &
            &                                           centroid, ref_centroid1, ref_centroid2, objective, gradient )
      end if
   end subroutine mof3d_compute_gradient_geometric

   pure subroutine mof3d_compute_gradient_finite_differences(polyhedron, cell_volume, cell_centroid, angles, ref_volume, &
      &                                                      centroid, ref_centroid1, ref_centroid2, objective, gradient )
      use mod_cg3_complete_polyhedron_structure
      use mod_cg3_flood_polyhedron
      use mod_cg3_polyhedron
      use variables_mof, only: mof_use_symmetric_reconstruction
      type(t_polyhedron), intent(in) :: polyhedron
      double precision, intent(in) :: cell_volume
      double precision, dimension(3), intent(in) :: cell_centroid
      double precision, dimension(2), intent(in) :: angles
      double precision, intent(in) :: ref_volume
      double precision, dimension(3), intent(in) :: centroid
      double precision, dimension(3), intent(in) :: ref_centroid1, ref_centroid2
      double precision, intent(out) :: objective
      double precision, dimension(2), intent(out) :: gradient

      double precision, parameter :: delta_angle = 1d-9

      type(t_polyhedron) :: polyhedron_full
      double precision, dimension(3) :: direction
      double precision, dimension(3) :: centroid_prev, centroid_next
      double precision, dimension(3) :: derivative_theta, derivative_phi
      double precision, dimension(3) :: dual_centroid, diff
      double precision :: volume
      integer :: error_id

      ! Compute d/dθ
      direction = [cos(angles(1) - delta_angle)*sin(angles(2)), &
         &         sin(angles(1) - delta_angle)*sin(angles(2)), &
         &         cos(angles(2))]
      call cg3_flood_polyhedron(polyhedron, direction, ref_volume, polyhedron_full)
      call cg3_complete_polyhedron_structure(polyhedron_full, error_id)
      call cg3_polyhedron_compute_centroid(polyhedron_full, volume, centroid_prev)

      direction = [cos(angles(1) + delta_angle)*sin(angles(2)), &
         &         sin(angles(1) + delta_angle)*sin(angles(2)), &
         &         cos(angles(2))]
      call cg3_flood_polyhedron(polyhedron, direction, ref_volume, polyhedron_full)
      call cg3_complete_polyhedron_structure(polyhedron_full, error_id)
      call cg3_polyhedron_compute_centroid(polyhedron_full, volume, centroid_next)

      derivative_theta = (centroid_next - centroid_prev)/(2d0*delta_angle)

      ! Compute d/dφ
      direction = [cos(angles(1))*sin(angles(2) - delta_angle), &
         &         sin(angles(1))*sin(angles(2) - delta_angle), &
         &         cos(angles(2) - delta_angle)]
      call cg3_flood_polyhedron(polyhedron, direction, ref_volume, polyhedron_full)
      call cg3_complete_polyhedron_structure(polyhedron_full, error_id)
      call cg3_polyhedron_compute_centroid(polyhedron_full, volume, centroid_prev)

      direction = [cos(angles(1))*sin(angles(2) + delta_angle), &
         &         sin(angles(1))*sin(angles(2) + delta_angle), &
         &         cos(angles(2) + delta_angle)]
      call cg3_flood_polyhedron(polyhedron, direction, ref_volume, polyhedron_full)
      call cg3_complete_polyhedron_structure(polyhedron_full, error_id)
      call cg3_polyhedron_compute_centroid(polyhedron_full, volume, centroid_next)

      derivative_phi = (centroid_next - centroid_prev)/(2d0*delta_angle)

      ! Compute the difference between the centroid and the reference centroid.
      diff = centroid - ref_centroid1

      ! Compute gradient
      gradient = 2d0*[dot_product(diff, derivative_theta), &
         &            dot_product(diff, derivative_phi)    ]

      objective = dot_product(diff, diff)

      ! Symmetric contribution
      if (mof_use_symmetric_reconstruction) then
         ! Compute the dual centroid
         dual_centroid = (cell_volume*cell_centroid - ref_volume*centroid)/(cell_volume - ref_volume)

         ! Compute the difference between the dual centroid and the second reference centroid.
         diff = dual_centroid - ref_centroid2

         ! Add the contribution of the symmetric part to the gradient.
         gradient = gradient + 2d0*(ref_volume/(ref_volume - cell_volume)) &
            &     * [dot_product(diff, derivative_theta),                  &
            &        dot_product(diff, derivative_phi)    ]

         ! Add the contribution of the symmetric part to the objective function.
         objective = objective + dot_product(diff, diff)
      end if
   end subroutine mof3d_compute_gradient_finite_differences

   ! Gradient computed using the method of Chen & Zhang
   ! Ref.: An improved 3d mof method based on analytical partial derivatives, Journal of Computational Physics 326, 156–170 (2016)
   pure subroutine mof3d_compute_geometric_gradient(polygon, cell_volume, cell_centroid, angles, ref_volume,   &
      &                                             centroid, ref_centroid1, ref_centroid2, objective, gradient)
      use mod_cg3_flood_polyhedron
      use mod_cg3_points
      use variables_mof, only: mof_use_symmetric_reconstruction
      type(t_chained_polygon), intent(in) :: polygon
      double precision, intent(in) :: cell_volume
      double precision, dimension(3), intent(in) :: cell_centroid
      double precision, dimension(2), intent(in) :: angles
      double precision, intent(in) :: ref_volume
      double precision, dimension(3), intent(in) :: centroid
      double precision, dimension(3), intent(in) :: ref_centroid1, ref_centroid2
      double precision, intent(out) :: objective
      double precision, dimension(2), intent(out) :: gradient

      double precision, dimension(3) :: i_phi, j_phi, k_phi
      double precision, dimension(3) :: origin, pref, translated_point
      double precision, dimension(2) :: c_diff, derivative_theta, derivative_phi
      double precision, dimension(:,:), allocatable :: p
      double precision :: cos_theta, sin_theta, cos_phi, sin_phi
      double precision :: integral_xy, integral_xx, integral_yy
      double precision :: area, triangle_area, x1, x2, x3, y1, y2, y3
      integer :: i

      ! Compute the origin of the coordinate system (the centroid of the polygon)
      area = 0d0
      origin = 0d0
      pref = polygon%point(1)%point

      do i = 2, polygon%nb_points - 1
         triangle_area = norm2(cg3_cross_product(polygon%point(i)%point - pref, polygon%point(i+1)%point - pref))
         area = area + triangle_area
         origin = origin + (polygon%point(i)%point - pref + polygon%point(i+1)%point - pref)*triangle_area
      end do

      area = area/2d0

      if (area > 0d0) then
         origin = pref + origin/(6d0*area)
      end if

      ! Compute the basis
      cos_theta = cos(angles(1))
      sin_theta = sin(angles(1))
      cos_phi = cos(angles(2))
      sin_phi = sin(angles(2))

      i_phi = [-sin_theta, cos_theta, 0d0]
      j_phi = [cos_theta*cos_phi, sin_theta*cos_phi, -sin_phi]
      k_phi = [cos_theta*sin_phi, sin_theta*sin_phi, cos_phi]

      allocate(p(3,polygon%nb_points))

      ! Project all the points
      do i = 1, polygon%nb_points
         translated_point = polygon%point(i)%point - origin
         p(1,i) = dot_product(translated_point, i_phi)
         p(2,i) = dot_product(translated_point, j_phi)
         p(3,i) = dot_product(translated_point, k_phi)
      end do

      ! Compute integrals
      integral_xx = 0d0
      integral_xy = 0d0
      integral_yy = 0d0

      x3 = (p(1,1) + p(1,2))/2d0
      y3 = (p(2,1) + p(2,2))/2d0
      do i = 2, polygon%nb_points - 1
         x1 = x3
         x2 = (p(1,i  ) + p(1,i+1))/2d0
         x3 = (p(1,i+1) + p(1,1  ))/2d0
         y1 = y3
         y2 = (p(2,i  ) + p(2,i+1))/2d0
         y3 = (p(2,i+1) + p(2,1  ))/2d0
         triangle_area = norm2(cg3_cross_product(p(:,i) - p(:,1), p(:,i+1) - p(:,1)))/2d0
         integral_xx = integral_xx + triangle_area*(x1**2 + x2**2 + x3**2)/3d0
         integral_xy = integral_xy + triangle_area*(x1*y1 + x2*y2 + x3*y3)/3d0
         integral_yy = integral_yy + triangle_area*(y1**2 + y2**2 + y3**2)/3d0
      end do

      ! Compute the partial derivatives
      derivative_theta = -sin_phi/ref_volume*[integral_xx, integral_xy]
      derivative_phi   = -1d0/ref_volume*[integral_xy, integral_yy]

      ! Project centroid difference
      translated_point = centroid - ref_centroid1
      c_diff(1) = dot_product(translated_point, i_phi)
      c_diff(2) = dot_product(translated_point, j_phi)

      ! Compute gradient
      gradient = 2d0*[dot_product(c_diff, derivative_theta), dot_product(c_diff, derivative_phi)]

      ! Compute the objective function
      objective = dot_product(translated_point, translated_point)

      ! Symmetric contribution
      if (mof_use_symmetric_reconstruction) then
         ! Project the difference between the dual centroid and the reference centroid of the symmetric part
         translated_point = (cell_volume*cell_centroid - ref_volume*centroid)/(cell_volume - ref_volume) - ref_centroid2
         c_diff(1) = dot_product(translated_point, i_phi)
         c_diff(2) = dot_product(translated_point, j_phi)

         ! Add the contribution of the symmetric part to the gradient.
         gradient = gradient + 2d0*(ref_volume/(ref_volume - cell_volume)) &
            &     * [dot_product(c_diff, derivative_theta), dot_product(c_diff, derivative_phi)]

         ! Add the contribution of the symmetric part to the objective function.
         objective = objective + dot_product(translated_point, translated_point)
      end if
   end subroutine mof3d_compute_geometric_gradient

   pure subroutine mof3d_debug_geometric_gradient(polyhedron, angles, ref_volume, centroid, derivative)
      use mod_cg3_polyhedron
      use mod_cg3_flood_polyhedron
      use mod_cg3_points
      type(t_polyhedron), intent(in) :: polyhedron
      double precision, dimension(2), intent(in) :: angles
      double precision, intent(in) :: ref_volume
      double precision, dimension(3), intent(out) :: centroid
      double precision, dimension(3,2), intent(out) :: derivative

      double precision, dimension(3) :: i_phi, j_phi, k_phi
      double precision, dimension(3) :: origin, pref, translated_point
      double precision, dimension(2) :: derivative_theta, derivative_phi
      double precision, dimension(:,:), allocatable :: p
      double precision :: cos_theta, sin_theta, cos_phi, sin_phi
      double precision :: integral_xy, integral_xx, integral_yy
      double precision :: area, triangle_area, x1, x2, x3, y1, y2, y3
      integer :: i
      double precision, dimension(3) :: direction
      type(t_chained_polygon) :: polygon

      ! Geometric methods
      sin_phi = sin(angles(2))

      ! Central reconstruction
      direction = [cos(angles(1))*sin_phi, &
         &         sin(angles(1))*sin_phi, &
         &         cos(angles(2))]

      call cg3_flood_polyhedron_centroid(polyhedron, direction, ref_volume, centroid, polygon)

      ! Compute the origin of the coordinate system (the centroid of the polygon)
      area = 0d0
      origin = 0d0
      pref = polygon%point(1)%point

      do i = 2, polygon%nb_points - 1
         triangle_area = norm2(cg3_cross_product(polygon%point(i)%point - pref, polygon%point(i+1)%point - pref))
         area = area + triangle_area
         origin = origin + (polygon%point(i)%point - pref + polygon%point(i+1)%point - pref)*triangle_area
      end do

      area = area/2d0

      if (area > 0d0) then
         origin = pref + origin/(6d0*area)
      end if

      ! Compute the basis
      cos_theta = cos(angles(1))
      sin_theta = sin(angles(1))
      cos_phi = cos(angles(2))
      sin_phi = sin(angles(2))

      i_phi = [-sin_theta, cos_theta, 0d0]
      j_phi = [cos_theta*cos_phi, sin_theta*cos_phi, -sin_phi]
      k_phi = [cos_theta*sin_phi, sin_theta*sin_phi, cos_phi]

      allocate(p(3,polygon%nb_points))

      ! Project all the points
      do i = 1, polygon%nb_points
         translated_point = polygon%point(i)%point - origin
         p(1,i) = dot_product(translated_point, i_phi)
         p(2,i) = dot_product(translated_point, j_phi)
         p(3,i) = dot_product(translated_point, k_phi)
      end do

      ! Compute integrals
      integral_xx = 0d0
      integral_xy = 0d0
      integral_yy = 0d0

      x3 = (p(1,1) + p(1,2))/2d0
      y3 = (p(2,1) + p(2,2))/2d0
      do i = 2, polygon%nb_points - 1
         x1 = x3
         x2 = (p(1,i  ) + p(1,i+1))/2d0
         x3 = (p(1,i+1) + p(1,1  ))/2d0
         y1 = y3
         y2 = (p(2,i  ) + p(2,i+1))/2d0
         y3 = (p(2,i+1) + p(2,1  ))/2d0
         triangle_area = norm2(cg3_cross_product(p(:,i) - p(:,1), p(:,i+1) - p(:,1)))/2d0
         integral_xx = integral_xx + triangle_area*(x1**2 + x2**2 + x3**2)/3d0
         integral_xy = integral_xy + triangle_area*(x1*y1 + x2*y2 + x3*y3)/3d0
         integral_yy = integral_yy + triangle_area*(y1**2 + y2**2 + y3**2)/3d0
      end do

      ! Compute the partial derivatives
      derivative_theta = -sin_phi/ref_volume*[integral_xx, integral_xy]
      derivative_phi   = -1d0/ref_volume*[integral_xy, integral_yy]

      derivative(:,1) = derivative_theta(1)*i_phi + derivative_theta(2)*j_phi
      derivative(:,2) = derivative_phi(1)*i_phi + derivative_phi(2)*j_phi
   end subroutine mof3d_debug_geometric_gradient

end module mod_mof3d_gradient

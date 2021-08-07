!This file is part of Notus 0.4.0

!Copyright Bordeaux-INP, Université de Bordeaux, CNRS
!
!Contributors:
!Antoine Lemoine, 30-06-2016, antoine.lemoine@bordeaux-inp.fr

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

!> @defgroup point_polyhedron Point/polyhedron relations
!! @brief Geometric tools relative to points and polyhedron
!! @ingroup computational_geometry_3d

module mod_cg3_point_polyhedron
   use mod_cg3_polyhedron
   implicit none

contains

   !> Find the closest point on a polyhedron face to a given point.
   !!
   !! @param[in]  polyhedron: any polyhedron
   !! @param[in]  face: face index of the @p polyhedron
   !! @param[in]  point: any point coordinates
   !! @param[out] closest_point: coordinates of the closest point on the given @p face
   !! @ingroup point_polyhedron
   pure subroutine cg3_point_face_closest_point(polyhedron, face, point, closest_point)
      type(t_polyhedron), intent(in) :: polyhedron
      integer, intent(in) :: face
      double precision, dimension(3), intent(in) :: point
      double precision, dimension(3), intent(out) :: closest_point

      double precision, dimension(3) :: projection
      double precision :: alpha, distance, current_distance
      integer :: i, edge

      ! Compute the orthogonal projection
      closest_point = point - dot_product(point - polyhedron%point(:,polyhedron%face(face)%id(1)), polyhedron%normal(:,face)) &
         &          * polyhedron%normal(:,face)

      ! If the orthogonal projection belongs to the face, return
      if (cg3_face_winding_number(polyhedron, face, closest_point) /= 0) return

      ! Find the closest point on the edges
      distance = huge(1d0)

      do i = 1, polyhedron%face_to_edge(face)%size
         edge = polyhedron%face_to_edge(face)%id(i)

         alpha = dot_product(point - polyhedron%point(:,polyhedron%edge(1,edge)), polyhedron%tangent(:,edge))

         if (alpha < 0d0) then
            projection = polyhedron%point(:,polyhedron%edge(1,edge))
         else if (alpha > norm2(polyhedron%point(:,polyhedron%edge(2,edge)) - polyhedron%point(:,polyhedron%edge(1,edge)))) then
            projection = polyhedron%point(:,polyhedron%edge(2,edge))
         else
            projection = polyhedron%point(:,polyhedron%edge(1,edge)) + alpha*polyhedron%tangent(:,edge)
         end if

         current_distance = norm2(point - projection)

         if (current_distance < distance) then
            distance = current_distance
            closest_point = projection
         end if
      end do
   end subroutine cg3_point_face_closest_point

   !> Compute the winding number of a point for a face of a polyhedron
   !!
   !! @param[in] polyhedron: any polyhedron
   !! @param[in] face: index of any face of the polyhedron
   !! @param[in] point: coordinates of any point that lies on the same plane as the face
   !! @ingroup point_polyhedron
   integer pure function cg3_face_winding_number(polyhedron, face, point) result(winding_number)
      use mod_cg2_points
      use mod_cg_transformation
      type(t_polyhedron), intent(in) :: polyhedron
      integer, intent(in) :: face
      double precision, dimension(3), intent(in) :: point

      double precision, parameter :: PI = 2d0*acos(0d0)

      type(cg_transformation) :: transformation
      double precision, dimension(2, polyhedron%face(face)%size) :: points
      double precision, dimension(3) :: axis, p3d
      double precision, dimension(2) :: p
      double precision :: cosine, sine
      integer :: i

      ! Rotation of the reference
      cosine = polyhedron%normal(3, face)
      sine = norm2([polyhedron%normal(1, face), polyhedron%normal(2, face)])

      ! Initialize the transformation matrix
      call cg_transformation_initialize(transformation, 3)

      ! Translate near the origin to limit numerical errors
      call cg_transformation_add_translation(transformation, -polyhedron%point(:,polyhedron%face(face)%id(1)))

      if (sine > epsilon(1d0)) then
         ! Axis of rotation
         axis = [polyhedron%normal(2, face)/sine, -polyhedron%normal(1, face)/sine, 0d0]

         ! Compute the transformation matrix
         call cg_transformation_add_rotation_cos_sin(transformation, axis, cosine, sine)
      else if (cosine < 0d0) then
         ! If the normal is equal to [0,0,-1], rotate the polyhedron of π around the x-axis
         call cg_transformation_add_rotation_x(transformation, PI)
      end if

      ! Apply the transformation
      do i = 1, polyhedron%face(face)%size
         p3d = cg_transform_point(transformation, polyhedron%point(:,polyhedron%face(face)%id(i)))
         points(:,i) = p3d(1:2)
      end do

      p3d = cg_transform_point(transformation, point)
      p = p3d(1:2)

      ! Compute winding number
      winding_number = 0

      do i = 1, polyhedron%face(face)%size - 1
         if (points(2,i) <= p(2)) then
            if (points(2,i+1) <= p(2)) cycle
            ! Check if the point is on the left of the segment
            if (.not. cg2_is_point_left_of_line(points(:,i), points(:,i+1), p)) cycle

            winding_number = winding_number + 1
         else
            if (points(2,i+1) > p(2)) cycle
            ! Check if the point is on the right of the segment
            if (cg2_is_point_left_of_line(points(:,i), points(:,i+1), p)) cycle

            winding_number = winding_number - 1
         end if
      end do

      ! Last segment [n+1, 1]
      if (points(2,polyhedron%face(face)%size) <= p(2)) then
         if (points(2,1) <= p(2)) return
         ! Check if the point is on the left of the segment
         if (.not. cg2_is_point_left_of_line(points(:,polyhedron%face(face)%size), points(:,1), p)) return

         winding_number = winding_number + 1
      else
         if (points(2,1) > p(2)) return
         ! Check if the point is on the right of the segment
         if (cg2_is_point_left_of_line(points(:,polyhedron%face(face)%size), points(:,1), p)) return

         winding_number = winding_number - 1
      end if
   end function cg3_face_winding_number

   !> Check if a sphere intersects a face polyhedron
   !!
   !! @param[in] polyhedron: any polyhedron
   !! @param[in] face: index of any face of the polyhedron
   !! @param[in] center: coordinates of the sphere center
   !! @param[in] radius: radius of the center
   !! @ingroup point_polyhedron
   logical pure function cg3_does_sphere_intersect_face(polyhedron, face, center, radius) result(is_intersection)
      type(t_polyhedron), intent(in) :: polyhedron
      integer, intent(in) :: face
      double precision, dimension(3), intent(in) :: center
      double precision, intent(in) :: radius

      double precision, dimension(3) :: face_center
      double precision :: face_radius, distance, center_center_distance
      integer :: i

      is_intersection = .true.

      ! Test with the first point
      if (norm2(center - polyhedron%point(:, polyhedron%face(face)%id(1))) <= radius) return

      ! Approximate face center
      face_center = polyhedron%point(:, polyhedron%face(face)%id(1))
      do i = 2, polyhedron%face(face)%size
         face_center = face_center + polyhedron%point(:, polyhedron%face(face)%id(i))
      end do

      face_center = face_center/dble(polyhedron%face(face)%size)

      center_center_distance = norm2(face_center - center)
      if (center_center_distance <= radius) return

      ! Compute face radius
      face_radius = norm2(face_center - polyhedron%point(:, polyhedron%face(face)%id(1)))
      do i = 2, polyhedron%face(face)%size
         distance = norm2(face_center - polyhedron%point(:, polyhedron%face(face)%id(i)))
         if (distance > face_radius) face_radius = distance
      end do

      if (center_center_distance > radius + face_radius) is_intersection = .false.
   end function cg3_does_sphere_intersect_face

end module mod_cg3_point_polyhedron

!This file is part of Notus 0.4.0

!Copyright Bordeaux-INP, Université de Bordeaux, CNRS
!
!Contributors:
!Antoine Lemoine, 04-07-2016, antoine.lemoine@bordeaux-inp.fr

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

!> @defgroup line_polyhedron Line/polyhedron intersection
!! @brief Computational geometry tools relative to lines and polyhedron intersection
!! @ingroup computational_geometry_3d

module mod_cg3_line_polyhedron
   use mod_cg3_point_polyhedron
   use mod_cg3_polyhedron
   implicit none

contains

   !> Determine if a ray (origin + direction) intersect a face of a given polyhedron.
   !!
   !! @param[in]  polyhedron: any polyhedron
   !! @param[in]  face: face index of the @p polyhedron
   !! @param[in]  origin: origin of the ray
   !! @param[in]  direction: direction of the ray
   !! @param[out] is_intersected: true if the ray intersects the face
   !! @param[out] is_origin_on_boundary: true if the origin of the ray belong to the face
   !! @ingroup line_polyhedron
   pure subroutine cg3_does_direction_intersect_face(polyhedron, face, origin, direction, is_intersected, is_origin_on_face)
      type(t_polyhedron), intent(in) :: polyhedron
      integer, intent(in) :: face
      double precision, dimension(3), intent(in) :: origin
      integer, intent(in) :: direction
      logical, intent(out) :: is_intersected
      logical, intent(out), optional :: is_origin_on_face

      double precision :: distance
      double precision, dimension(3) :: intersection, closest_point

      is_intersected = .false.
      if (present(is_origin_on_face)) is_origin_on_face = .false.

      ! Check if the direction of the ray is orthogonal to the face.
      if (abs(polyhedron%normal(direction, face)) <= epsilon(1d0)) then
         if (.not. present(is_origin_on_face)) return

         ! Check if the origin of the ray belongs to the face.
         call cg3_point_face_closest_point(polyhedron, face, origin, closest_point)

         ! If the point is too far from the face, it is considered outside the shape.
         if (dot_product(origin - closest_point, polyhedron%normal(:,face)) > tiny(1d0)) return

         ! Check if the closest point belongs to the face.
         if (cg3_face_winding_number(polyhedron, face, closest_point) /= 0) is_intersected = .true.

         ! Set is_origin_on_face flag.
         is_origin_on_face = is_intersected

         return
      end if

      ! Compute the distance from the origin of the ray to the the face along the direction of the ray.
      distance = dot_product(polyhedron%point(:,polyhedron%face(face)%id(1)) - origin, &
         &                   polyhedron%normal(:,face))                                &
         &     / polyhedron%normal(direction, face)

      if (distance < 0d0) return

      intersection = origin
      intersection(direction) = intersection(direction) + distance

      if (cg3_face_winding_number(polyhedron, face, intersection) /= 0) is_intersected = .true.

      if (present(is_origin_on_face)) then
         is_origin_on_face = is_intersected .and. (abs(distance) <= tiny(distance))
      end if
   end subroutine cg3_does_direction_intersect_face

   !> Intersect a face of a given polyhedron with a ray.
   !!
   !! @param[in]  polyhedron: any polyhedron
   !! @param[in]  face: face index of the @p polyhedron
   !! @param[in]  ray: ray
   !! @param[in]  has_normal: flag to determine if the normal vector must be computed
   !! @param[out] is_intersected: true if the ray intersects the face
   !! @param[out] intersection: intersection
   !! @ingroup line_polyhedron
   pure subroutine cg3_compute_ray_face_intersection(polyhedron, face, ray, has_normal, is_intersected, intersection)
      use mod_ray_tracing
      type(t_polyhedron), intent(in) :: polyhedron
      integer, intent(in) :: face
      type(t_ray), intent(in) :: ray
      logical, intent(in) :: has_normal
      logical, intent(out) :: is_intersected
      type(t_intersection), intent(out) :: intersection

      double precision :: distance, direction_normal
      double precision, dimension(3) :: point

      is_intersected = .false.

      ! Dot product between the normal vector to the face and the direction of the ray.
      direction_normal = dot_product(ray%direction, polyhedron%normal(:,face))

      ! Ignore the intersection if the ray is parallel to the face.
      if (abs(direction_normal) <= tiny(1d0)) return

      ! Compute the distance from the origin of the ray to the the face.
      distance = dot_product(polyhedron%point(:,polyhedron%face(face)%id(1)) - ray%origin, polyhedron%normal(:,face)) &
         &     / direction_normal

      ! Compute the intersection point.
      point = ray%origin + ray%direction*distance

      ! Check whether the intersection point belongs to the face.
      if (cg3_face_winding_number(polyhedron, face, point) == 0) return

      ! Set intersection flag.
      is_intersected = .true.

      ! Set the intersection depending on the has_normal flag.
      if (has_normal) then
         intersection = t_intersection(point, distance, polyhedron%normal(:,face))
      else
         intersection = t_intersection(point, distance)
      end if
   end subroutine cg3_compute_ray_face_intersection

end module mod_cg3_line_polyhedron

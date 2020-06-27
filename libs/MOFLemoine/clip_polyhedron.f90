!This file is part of Notus 0.4.0

!Copyright Bordeaux-INP, Université de Bordeaux, CNRS
!
!Contributors:
!Antoine Lemoine, 22-11-2017, antoine.lemoine@bordeaux-inp.fr

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

module mod_cg3_clip_polyhedron
   use mod_cg3_polyhedron
   implicit none
   private

   public :: cg3_polyhedron_clipping, cg3_plane_polyhedron_clipping

   type t_incidence_graph
      !> Each line contains the list of vertices connected to the vertex corresponding to the line
      integer, dimension(:,:), allocatable :: id
      !> Number of connected element in a given line
      integer, dimension(:), allocatable :: row_size
   end type t_incidence_graph

   !> @defgroup polyhedron_clipping Polyhedron clipping
   !! @brief Clip a polyhedron with a half-space or a convex polyhedron
   !! @ingroup computational_geometry_3d
   !!
   !! References:
   !!  - Sugihara, K. (1994). A robust and consistent algorithm for intersecting convex polyhedra. In: Computer Graphics Forum.
   !!    Blackwell Science Ltd. p. 45-54. doi:[10.1111/1467-8659.1330045]
   !!
   !! [10.1111/1467-8659.1330045]: https://doi.org/10.1111/1467-8659.1330045 "10.1111/1467-8659.1330045"

contains

   !> Clip a convex polyhedron with a convex polyhedron
   !! @param[in]  polyhedron:         convex polyhedron
   !! @param[in]  clip_polyhedron:    convex clipping polyhedron
   !! @param[out] clipped_polyhedron: clipped polyhedron (convex by construction)
   !! @param[out] is_empty:           true if all the points are outside the clipped polyhedron
   !! @param[out] is_clipped:         true if polyhedron is clipped
   !! @ingroup polyhedron_clipping
   pure subroutine cg3_polyhedron_clipping(polyhedron, clip_polyhedron, clipped_polyhedron, is_empty, is_clipped)
      use mod_cg3_complete_polyhedron_structure
      type(t_polyhedron), intent(in) :: polyhedron, clip_polyhedron
      type(t_polyhedron), allocatable, intent(out) :: clipped_polyhedron
      logical, intent(out) :: is_empty
      logical, intent(out) :: is_clipped

      type(t_polyhedron), allocatable :: tmp_polyhedron, tmp_clipped_polyhedron
      integer :: face
      logical :: is_clip_empty, is_cut

      is_clipped = .false.
      is_empty   = .false.

      ! Loop on all the clip polyhedron faces
      do face = 1, clip_polyhedron%nb_faces
         ! Clip the remaining part of the polyedron with the next plane
         if (is_clipped) then
            call cg3_plane_polyhedron_clipping(tmp_clipped_polyhedron, clip_polyhedron%normal(:,face),                    &
               &                               clip_polyhedron%point(:,clip_polyhedron%face(face)%id(1)), tmp_polyhedron, &
               &                               is_clip_empty, is_cut)
         else
            call cg3_plane_polyhedron_clipping(polyhedron, clip_polyhedron%normal(:,face),                                &
               &                               clip_polyhedron%point(:,clip_polyhedron%face(face)%id(1)), tmp_polyhedron, &
               &                               is_clip_empty, is_cut)
         end if

         if (is_cut) then
            call move_alloc(tmp_polyhedron, tmp_clipped_polyhedron)
            is_clipped = .true.
         else if (is_clip_empty) then
            ! The clipped polyhedron is empty. Return now
            if (allocated(tmp_clipped_polyhedron)) deallocate(tmp_clipped_polyhedron)
            is_empty = .true.
            is_clipped = .false.
            return
         end if
      end do

      if (is_clipped) call move_alloc(tmp_clipped_polyhedron, clipped_polyhedron)
   end subroutine cg3_polyhedron_clipping

   !> Clip a convex polyhedron with a plane (half-space) using Sugihara's algorithm
   !! @param[in]  polyhedron:         convex polyhedron
   !! @param[in]  normal:             outgoing normal of the half-space
   !! @param[in]  origin:             point of the clipping plane
   !! @param[out] clipped_polyhedron: clipped polyhedron (convex by construction)
   !! @param[out] is_empty:           true if all the points are above the clipping plane
   !! @param[out] is_cut:             true if the plane cut the polyhedron
   !! @ingroup polyhedron_clipping
   pure subroutine cg3_plane_polyhedron_clipping(polyhedron, normal, origin, clipped_polyhedron, is_empty, is_cut)
      use mod_cg3_points
      type(t_polyhedron), intent(in) :: polyhedron
      double precision, dimension(3), intent(in) :: normal
      double precision, dimension(3), intent(in) :: origin
      type(t_polyhedron), allocatable, intent(out) :: clipped_polyhedron
      logical, intent(out) :: is_empty
      logical, intent(out) :: is_cut

      type(t_incidence_graph) :: vertex_graph
      logical, dimension(polyhedron%nb_points) :: is_point_below
      integer :: i, j, max_degree

      is_empty = .false.
      is_cut   = .true.

      ! Compute the vertex partition. 'is_point_below(i)' is true if the point is below the clipping plane.
      do i = 1, polyhedron%nb_points
         is_point_below(i) = .not. (dot_product(polyhedron%point(:,i) - origin, normal) > 0d0)
      end do

      ! Check if the plane actually cut the polyhedron
      if (is_point_below(1)) then
         ! If all the points are below, set is_cut to false and return (clipped_polyhedron is not affected)
         if (all(is_point_below)) then
            is_empty = .false.
            is_cut   = .false.
            return
         end if
      else
         ! If all the points are above, set is_cut to false and return (clipped_polyhedron is not affected)
         if (.not. any(is_point_below)) then
            is_empty = .true.
            is_cut   = .false.
            return
         end if
      end if

      ! Check if the above/below sub-graphs are connected

      ! Find the vertex of the polyhedron of maximum degree
      max_degree = polyhedron%point_to_edge(1)%size
      do i = 2, polyhedron%nb_points
         if (polyhedron%point_to_edge(i)%size > max_degree) max_degree = polyhedron%point_to_edge(i)%size
      end do

      ! Allocate the vertex-edge graph represented by a vertex-vertex connection graph
      allocate(vertex_graph%id(max_degree, polyhedron%nb_points), source=0)
      allocate(vertex_graph%row_size(polyhedron%nb_points), source=0)

      ! Construct the graph
      do i = 1, polyhedron%nb_points
         vertex_graph%row_size(i) = polyhedron%point_to_edge(i)%size

         do j = 1, polyhedron%point_to_edge(i)%size
            if (polyhedron%edge(1,polyhedron%point_to_edge(i)%id(j)) == i) then
               vertex_graph%id(j,i) = polyhedron%edge(2,polyhedron%point_to_edge(i)%id(j))
            else
               vertex_graph%id(j,i) = polyhedron%edge(1,polyhedron%point_to_edge(i)%id(j))
            end if
         end do
      end do

      ! Check if the part below and the part above the clipping plane are connected
      if (.not. is_connected(vertex_graph, is_point_below, .true.)) then
         call combinatorial_partition(vertex_graph, is_point_below)
      else if (.not. is_connected(vertex_graph, is_point_below, .false.)) then
         call combinatorial_partition(vertex_graph, is_point_below)
      end if

      ! Compute the clipped polyhedron
      call partition_polyhedron(polyhedron, normal, origin, is_point_below, clipped_polyhedron)

      ! Check if we have a clipped polyhedron and set flags accordingly
      if (.not. allocated(clipped_polyhedron)) then
         is_empty = .true.
         is_cut = .false.
      end if
   end subroutine cg3_plane_polyhedron_clipping

   !> Reconstruct the part of the polyhedron that belongs to the half-space
   !! @param[in]  polyhedron:         any convex polyhedron
   !! @param[in]  normal:             outgoing normal of the half-space
   !! @param[in]  origin:             point of the clipping plane
   !! @param[in]  partition:          partition of vertices as a boolean array
   !! @param[out] clipped_polyhedron: clipped polyhedron (convex by construction)
   !! @ingroup polyhedron_clipping
   pure subroutine partition_polyhedron(polyhedron, normal, origin, partition, clipped_polyhedron)
      use mod_cg3_complete_polyhedron_structure
      type(t_polyhedron), intent(in) :: polyhedron
      double precision, dimension(3), intent(in) :: normal
      double precision, dimension(3), intent(in) :: origin
      logical, dimension(:), intent(in) :: partition
      type(t_polyhedron), allocatable, intent(out) :: clipped_polyhedron

      type(t_incidence_matrix), dimension(:), allocatable :: tmp_face
      double precision, dimension(3,polyhedron%nb_faces) :: cap_point
      double precision, dimension(3) :: intersection
      integer, dimension(polyhedron%nb_faces) :: cap_point_id
      integer, dimension(polyhedron%nb_points) :: new_id
      integer, dimension(:), allocatable :: point_buffer
      logical, dimension(polyhedron%nb_faces) :: clipped_face
      integer :: nb_faces, nb_cap_point, nb_new_point, max_face_point, nb_face_point, error_id
      integer :: edge, face, first_face, next_face, inside_point, outside_point
      integer :: i, j, curp, prevp, in_id, in_end, out_id, out_end, first_id, first_end
      logical :: is_inside

      ! Allocate the clipped polyhedron
      allocate(clipped_polyhedron)

      ! Create a new index of points in the clipped polyhedron
      clipped_polyhedron%nb_points = 0
      do i = 1, polyhedron%nb_points
         if (partition(i)) then
            clipped_polyhedron%nb_points = clipped_polyhedron%nb_points + 1
            new_id(i) = clipped_polyhedron%nb_points
         else
            new_id(i) = 0
         end if
      end do

      ! Count the number of faces of the clipped polyhedron
      clipped_face = .false.
      nb_faces = 0
      max_face_point = 0
      do i = 1, polyhedron%nb_faces
         ! Find the face with the maximum number of point
         if (polyhedron%face(i)%size > max_face_point) max_face_point = polyhedron%face(i)%size

         if (any(partition(polyhedron%face(i)%id))) then
            clipped_face(i) = .true.
            nb_faces = nb_faces + 1
         end if
      end do

      ! Add the cap polygon
      nb_faces = nb_faces + 1

      ! Allocate memory for the faces
      allocate(tmp_face(nb_faces))

      ! The face/plane intersection can create an extra point on a face
      max_face_point = max_face_point + 1

      ! Allocate the point buffer
      allocate(point_buffer(max_face_point))

      ! Find the first edge with two colors
      edge = 0
      do i = 1, polyhedron%nb_edges
         if (partition(polyhedron%edge(1,i)) .eqv. partition(polyhedron%edge(2,i))) cycle
         edge = i
         exit
      end do

      ! Find the point of the edge that is inside the half-space
      if (partition(polyhedron%edge(1,edge))) then
         inside_point  = polyhedron%edge(1,edge)
         outside_point = polyhedron%edge(2,edge)
      else
         inside_point  = polyhedron%edge(2,edge)
         outside_point = polyhedron%edge(1,edge)
      end if

      ! Find the first face
      first_face = polyhedron%edge_to_face(1,edge)
      if (polyhedron%face(first_face)%id(1) == inside_point) then
         if (polyhedron%face(first_face)%id(polyhedron%face(first_face)%size) /= outside_point) then
            first_face = polyhedron%edge_to_face(2,edge)
         end if
      else
         do i = 2, polyhedron%face(first_face)%size
            if (polyhedron%face(first_face)%id(i) /= inside_point) cycle
            if (polyhedron%face(first_face)%id(i-1) /= outside_point) first_face = polyhedron%edge_to_face(2,edge)
            exit
         end do
      end if

      ! Initialize the counter of new point generated by the cap polygon
      nb_new_point = 0

      ! Compute the coordinates of the first point of the cap polygon
      nb_cap_point = 1
      call edge_plane_intersection(polyhedron, inside_point, outside_point, normal, origin, out_end, cap_point(:,nb_cap_point))

      if (out_end == inside_point) then
         cap_point_id(nb_cap_point) = new_id(out_end)
      else
         ! Increment the number of new point in the cap polygon
         nb_new_point = nb_new_point + 1
         cap_point_id(nb_cap_point) = clipped_polyhedron%nb_points + nb_new_point
      end if

      out_id    = cap_point_id(nb_cap_point)
      first_id  = out_id
      first_end = out_end

      ! Loop on every cut faces until we reach the first face
      face = first_face
      next_face = face
      clipped_polyhedron%nb_faces = 0
      do
         in_id = out_id
         in_end = out_end
         out_id = -1
         out_end = -1

         ! Remove the face of the clipped face list
         clipped_face(face) = .false.

         ! Loop on every point of the face
         prevp = polyhedron%face(face)%size
         is_inside = partition(polyhedron%face(face)%id(prevp))
         nb_face_point = 0
         do curp = 1, polyhedron%face(face)%size
            associate(prevp_id => polyhedron%face(face)%id(prevp), curp_id => polyhedron%face(face)%id(curp))
               ! Detect if the current point is inside the clipped polyhedron
               if (partition(curp_id)) then
                  ! Detect transition outside → inside
                  if (is_inside) then
                     ! No transition
                     ! -------------

                     if (nb_face_point > 0) then
                        if (point_buffer(1) /= new_id(curp_id)) then
                           ! Add the current point to the face
                           nb_face_point = nb_face_point + 1
                           point_buffer(nb_face_point) = new_id(curp_id)
                        end if
                     else
                        ! Add the current point to the face
                        nb_face_point = nb_face_point + 1
                        point_buffer(nb_face_point) = new_id(curp_id)
                     end if
                  else
                     ! Transition outside → inside
                     ! ---------------------------

                     ! Set inside flag
                     is_inside = .true.

                     ! Add intersection point
                     if (nb_face_point > 0) then
                        if (point_buffer(nb_face_point) /= in_id) then
                           ! Add the entry point to the face
                           nb_face_point = nb_face_point + 1
                           point_buffer(nb_face_point) = in_id
                        end if
                     else
                        ! Add the entry point to the face
                        nb_face_point = nb_face_point + 1
                        point_buffer(nb_face_point) = in_id
                     end if

                     ! Avoid to add the current point twice
                     if (point_buffer(nb_face_point) /= new_id(curp_id) .and. new_id(curp_id) /= out_id) then
                        ! Add the current point to the face
                        nb_face_point = nb_face_point + 1
                        point_buffer(nb_face_point) = new_id(curp_id)
                     end if
                  end if
               else if (is_inside) then
                  ! Transition inside → outside
                  ! ---------------------------

                  ! Find the transition edge
                  ! The current point and the previous point belong to the edge common to the next face
                  do i = 1, polyhedron%point_to_edge(prevp_id)%size
                     if (polyhedron%edge(1,polyhedron%point_to_edge(prevp_id)%id(i)) == prevp_id) then
                        if (polyhedron%edge(2,polyhedron%point_to_edge(prevp_id)%id(i)) /= curp_id) cycle
                        ! We found the edge
                        edge = polyhedron%point_to_edge(prevp_id)%id(i)
                        exit
                     else
                        if (polyhedron%edge(1,polyhedron%point_to_edge(prevp_id)%id(i)) /= curp_id) cycle
                        ! We found the edge
                        edge = polyhedron%point_to_edge(prevp_id)%id(i)
                        exit
                     end if
                  end do

                  ! Find the next face
                  ! The next face is in the opposite side of the edge
                  if (polyhedron%edge_to_face(1,edge) == face) then
                     next_face = polyhedron%edge_to_face(2,edge)
                  else
                     next_face = polyhedron%edge_to_face(1,edge)
                  end if

                  ! Compute the intersection point
                  if (next_face /= first_face) then
                     ! Compute edge/plane intersection
                     call edge_plane_intersection(polyhedron, curp_id, prevp_id, normal, origin, out_end, intersection)

                     ! Check whether the exit point is an edge extremity
                     if (out_end == 0) then
                        ! The point is not at the extremity of the edge

                        ! Add the point to the cap point
                        nb_cap_point = nb_cap_point + 1
                        cap_point(:,nb_cap_point) = intersection

                        ! Create a new id for this point
                        nb_new_point = nb_new_point + 1
                        cap_point_id(nb_cap_point) = clipped_polyhedron%nb_points + nb_new_point

                        ! Set out id
                        out_id = cap_point_id(nb_cap_point)

                        ! Add the point to the current face
                        nb_face_point = nb_face_point + 1
                        point_buffer(nb_face_point) = out_id
                     else if (out_end /= first_end) then
                        ! The end point is not the same as the end point of the first face

                        ! Check whether we does not generate duplicate points when the entry point corresponds to the exit point
                        if (out_end /= in_end) then
                           ! Add the point to the cap point
                           nb_cap_point = nb_cap_point + 1
                           ! Determine if we need to generate a new point
                           if (out_end == prevp_id) then
                              ! The out point is inside the clipped polyhedron: no need to generate a new point
                              cap_point_id(nb_cap_point) = new_id(out_end)
                           else
                              ! Generate a new point
                              nb_new_point = nb_new_point + 1
                              cap_point_id(nb_cap_point) = clipped_polyhedron%nb_points + nb_new_point
                           end if

                           ! Set out id
                           out_id = cap_point_id(nb_cap_point)

                           ! Add intersection point to the cap
                           cap_point(:,nb_cap_point) = intersection

                           if (nb_face_point > 0) then
                              if (point_buffer(nb_face_point) /= out_id) then
                                 ! Add the point to the current face
                                 nb_face_point = nb_face_point + 1
                                 point_buffer(nb_face_point) = out_id
                              end if
                           else
                              ! Add the point to the current face
                              nb_face_point = nb_face_point + 1
                              point_buffer(nb_face_point) = out_id
                           end if
                        else
                           ! Set out id
                           out_id = in_id
                        end if
                     else
                        ! The extremity of the edge corresponds to the entry point of the first face

                        if (first_id /= in_id) then
                           if (nb_face_point > 0) then
                              if (point_buffer(nb_face_point) /= first_id) then
                                 ! Add the point to the current face
                                 nb_face_point = nb_face_point + 1
                                 point_buffer(nb_face_point) = first_id
                              end if
                           else
                              ! Add the point to the current face
                              nb_face_point = nb_face_point + 1
                              point_buffer(nb_face_point) = first_id
                           end if
                        end if

                        ! Set out id and end
                        out_id = first_id
                        out_end = first_end
                     end if
                  else if (in_id /= first_id) then
                     ! The entry point end of the last face is different from the entry point end of the first face
                     if (nb_face_point >= 1) then
                        if (point_buffer(nb_face_point) /= first_id) then
                           ! Add the first point to the face
                           nb_face_point = nb_face_point + 1
                           point_buffer(nb_face_point) = first_id
                           ! Set out id and end
                           out_id = first_id
                           out_end = first_end
                        end if
                     else
                        ! Add the first point to the face
                        nb_face_point = nb_face_point + 1
                        point_buffer(nb_face_point) = first_id
                        ! Set out id and end
                        out_id = first_id
                        out_end = first_end
                     end if
                  end if

                  is_inside = .false.
               end if

               prevp = curp
            end associate
         end do

         if (nb_face_point > 2) then
            ! Increment the number of faces
            clipped_polyhedron%nb_faces = clipped_polyhedron%nb_faces + 1

            ! Register the new face
            allocate(tmp_face(clipped_polyhedron%nb_faces)%id(nb_face_point), source=point_buffer(:nb_face_point))
            tmp_face(clipped_polyhedron%nb_faces)%size = nb_face_point
         end if

         ! Switch to next face
         face = next_face

         if (face == first_face) exit
      end do

      ! Add the cap face if it is not singular
      if (nb_cap_point > 2) then
         ! Increment the number of faces
         clipped_polyhedron%nb_faces = clipped_polyhedron%nb_faces + 1
         ! Register the cap polygon
         allocate(tmp_face(clipped_polyhedron%nb_faces)%id(nb_cap_point), source=cap_point_id(:nb_cap_point))
         tmp_face(clipped_polyhedron%nb_faces)%size = nb_cap_point
      end if

      ! Register the remaining faces (those that are inside and not cut by the plane)
      do i = 1, polyhedron%nb_faces
         if (clipped_face(i)) then
            ! Increment the number of faces
            clipped_polyhedron%nb_faces = clipped_polyhedron%nb_faces + 1
            ! Register the face
            allocate(tmp_face(clipped_polyhedron%nb_faces)%id(polyhedron%face(i)%size), source=new_id(polyhedron%face(i)%id))
            ! Add the size of the face
            tmp_face(clipped_polyhedron%nb_faces)%size = polyhedron%face(i)%size
         end if
      end do

      if (clipped_polyhedron%nb_faces == size(tmp_face)) then
         call move_alloc(tmp_face, clipped_polyhedron%face)
      else
         allocate(clipped_polyhedron%face(clipped_polyhedron%nb_faces))
         do i = 1, clipped_polyhedron%nb_faces
            call move_alloc(tmp_face(i)%id, clipped_polyhedron%face(i)%id)
            clipped_polyhedron%face(i)%size = tmp_face(i)%size
         end do
      end if

      ! Allocate memory for point coordinates
      allocate(clipped_polyhedron%point(3,clipped_polyhedron%nb_points + nb_new_point))

      ! Add old points
      j = 0
      do i = 1, polyhedron%nb_points
         if (new_id(i) == 0) cycle
         j = j + 1
         clipped_polyhedron%point(:,j) = polyhedron%point(:,i)
      end do

      ! Add new points
      do i = 1, nb_cap_point
         if (cap_point_id(i) <= clipped_polyhedron%nb_points) cycle
         j = j + 1
         clipped_polyhedron%point(:,j) = cap_point(:,i)
      end do

      ! Set the number of point in the clipped polyhedron structure
      clipped_polyhedron%nb_points = clipped_polyhedron%nb_points + nb_new_point

      ! If the polyhedron has less than 4 points, consider there is no intersection.
      if (clipped_polyhedron%nb_points < 4) then
         deallocate(clipped_polyhedron)
         return
      end if

      ! Complete the polyhedron structure
      call cg3_complete_polyhedron_structure(clipped_polyhedron, error_id)

      if (error_id == 1 .or. error_id == 3) deallocate(clipped_polyhedron)
   end subroutine partition_polyhedron

   !> Compute the intersection of a plane and an edge
   pure subroutine edge_plane_intersection(polyhedron, p1, p2, normal, origin, end_id, intersection_point)
      type(t_polyhedron), intent(in) :: polyhedron
      integer, intent(in) :: p1, p2
      double precision, dimension(3), intent(in) :: normal
      double precision, dimension(3), intent(in) :: origin
      integer, intent(out) :: end_id
      double precision, dimension(3), intent(out) :: intersection_point

      double precision, parameter :: THRESHOLD = 1d4*epsilon(1d0)
      double precision, dimension(3) :: u, v
      double precision :: alpha, u_n, v_n

      end_id = 0

      u = polyhedron%point(:,p2) - polyhedron%point(:,p1)
      v = origin - polyhedron%point(:,p1)

      u_n = dot_product(u, normal)

      ! Check if the edge is parallel to the plane
      if (abs(u_n) < tiny(1d0)) then
         end_id = p1
         intersection_point = polyhedron%point(:,p1)
         return
      end if

      v_n = dot_product(v, normal)
      alpha = v_n/u_n

      if (alpha < THRESHOLD) then
         end_id = p1
         intersection_point = polyhedron%point(:,p1)
         return
      else if (1d0 - alpha < THRESHOLD) then
         end_id = p2
         intersection_point = polyhedron%point(:,p2)
         return
      end if

      intersection_point = (1d0 - alpha)*polyhedron%point(:,p1) + alpha*polyhedron%point(:,p2)

      ! Avoid to create short edges. A posteriori correction
      if (norm2(intersection_point - polyhedron%point(:,p1)) < THRESHOLD) then
         end_id = p1
         intersection_point = polyhedron%point(:,p1)
      else if (norm2(intersection_point - polyhedron%point(:,p2)) < THRESHOLD) then
         end_id = p2
         intersection_point = polyhedron%point(:,p2)
      end if
   end subroutine edge_plane_intersection

   !> Partition the vertex-graph of a polyhedron with the Sugihara's algorithm
   pure subroutine combinatorial_partition(graph, partition)
      type(t_incidence_graph), intent(in) :: graph
      logical, dimension(:), intent(inout) :: partition

      logical, dimension(size(partition)) :: new_partition, is_visited
      logical :: inside_flag
      integer :: i, first

      is_visited = .false.

      ! Find the biggest part to build the smallest partition
      if (count(partition) > size(partition)/2) then
         inside_flag = .false.
      else
         inside_flag = .true.
      end if

      ! Initialize the new partition
      new_partition = .not. inside_flag

      ! Find the first vertex inside
      do i = 1, size(graph%row_size)
         ! If the current vertex is not below the cut plane, cycle.
         if (partition(i) .neqv. inside_flag) cycle
         ! Retain the index of the first element
         first = i
         exit
      end do

      ! Add the first point to the list
      new_partition(first) = inside_flag

      ! Mark the first point as visited
      is_visited(first) = .true.

      ! Use a depth-first search algorithm to construct the inside connected component
      call depth_first_search_partition(graph, partition, is_visited, new_partition, inside_flag, first)

      ! Return the new partition
      partition = new_partition
   end subroutine combinatorial_partition

   pure recursive subroutine depth_first_search_partition(graph, partition, is_visited, new_partition, inside_flag, i)
      type(t_incidence_graph), intent(in) :: graph
      logical, dimension(:), intent(in) :: partition
      logical, dimension(:), intent(inout) :: is_visited
      logical, dimension(:), intent(inout) :: new_partition
      logical, intent(in) :: inside_flag
      integer, intent(in) :: i

      integer :: j

      ! Look around each point connected to the current vertex
      do j = 1, graph%row_size(i)
         ! If the vertex has been already visited, do not visit it
         if (is_visited(graph%id(j,i))) cycle

         ! Try outside (no need to change the value, it is already initialized as outside)
         if (is_connected(graph, new_partition, .not. inside_flag)) then
            ! Try inside
            new_partition(graph%id(j,i)) = inside_flag
            ! Check if the outside part is always connected
            if (is_connected(graph, new_partition, .not. inside_flag)) then
               ! It is not possible to decide if the point is inside or outside ⇒ use the geometric partition
               new_partition(graph%id(j,i)) = partition(graph%id(j,i))
            else
               ! Impossible to apply inside. It must be outside.
               new_partition(graph%id(j,i)) = .not. inside_flag
            end if
         else
            ! Impossible to apply outside. It must be inside.
            new_partition(graph%id(j,i)) = inside_flag
         end if

         ! Mark as visited
         is_visited(graph%id(j,i)) = .true.

         ! Do not visit the neighbor vertex if the current vertex is outside
         if (new_partition(graph%id(j,i)) .neqv. inside_flag) cycle

         ! Visit the vertices connected to the current vertex
         call depth_first_search_partition(graph, partition, is_visited, new_partition, inside_flag, graph%id(j,i))
      end do
   end subroutine depth_first_search_partition

   logical pure function is_connected(graph, sub_graph, inside_flag)
      type(t_incidence_graph), intent(in) :: graph
      logical, dimension(:), intent(in) :: sub_graph
      logical, intent(in) :: inside_flag

      integer :: i, first
      logical, dimension(size(sub_graph)) :: is_visited

      is_connected = .true.
      is_visited = .false.

      ! Find the first element that belong to the sub-graph
      do i = 1, size(graph%row_size)
         ! If the current vertex does not belong to the sub-graph, cycle.
         if (sub_graph(i) .neqv. inside_flag) cycle
         ! Retain the index of the first element
         first = i
         exit
      end do

      ! Mark the first element as visited
      is_visited(first) = .true.

      ! Use a depth-first search algorithm to find the first connected component
      call depth_first_search(graph, sub_graph, is_visited, inside_flag, first)

      ! Check if all the remaining vertices has been visited by the depth-first search algorithm
      do i = first + 1, size(graph%row_size)
         ! Do not visit vertices that does not belong to the sub-graph
         if (sub_graph(i) .neqv. inside_flag) cycle
         ! Check if the point is visited
         if (is_visited(i)) cycle
         ! The point is not visited ⇒ the sub-graph is not connected ⇒ return
         is_connected = .false.
         return
      end do
   end function is_connected

   pure recursive subroutine depth_first_search(graph, sub_graph, is_visited, inside_flag, i)
      type(t_incidence_graph), intent(in) :: graph
      logical, dimension(:), intent(in) :: sub_graph
      logical, dimension(:), intent(inout) :: is_visited
      logical, intent(in) :: inside_flag
      integer, intent(in) :: i

      integer :: j

      ! Look around each point connected to the current vertex
      do j = 1, graph%row_size(i)
         ! Do not visit vertices that not belong to the sub-graph
         if (sub_graph(graph%id(j,i)) .neqv. inside_flag) cycle
         ! Do not visit vertices that has been already visited
         if (is_visited(graph%id(j,i))) cycle
         ! Mark this vertex as visited
         is_visited(graph%id(j,i)) = .true.
         ! Visit the vertices connected to the current vertex
         call depth_first_search(graph, sub_graph, is_visited, inside_flag, graph%id(j,i))
      end do
   end subroutine depth_first_search

end module mod_cg3_clip_polyhedron

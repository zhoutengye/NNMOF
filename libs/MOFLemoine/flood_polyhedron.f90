!This file is part of Notus 0.4.0

!Copyright Bordeaux-INP, Université de Bordeaux, CNRS
!
!Contributors:
!Antoine Lemoine, 16-02-2016, antoine.lemoine@bordeaux-inp.fr

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

!> @defgroup flood_polyhedron Polyhedron flooding
!! @brief Computational geometry tool to fill a convex polyhedron with a given volume fraction in a given direction
!! @ingroup computational_geometry_3d

module mod_cg3_flood_polyhedron
   use mod_cg3_points
   use mod_cg3_polyhedron
   use mod_cg_transformation
   implicit none
   private

   !> Single item of the structure that track the state of the current algorithm.
   !! This structure contains a polygon stored as a list of point and edges.
   !! It also embbed information about neighbor edge, neighbor face, and if the current
   !! point is a point of the polyhedron.
   !! @ingroup flood_polyhedron
   type t_chained_polygon_item
      !> Coordinates of the point
      double precision, dimension(3) :: point = 0
      !> Set of edges that go up
      integer, dimension(:), allocatable :: edge_set
      !> Opposite point of the edge
      integer, dimension(:), allocatable :: edge_point
      !> Tangents of the edges that go up
      double precision, dimension(:,:), allocatable :: tangent
      !> Corresponding point number in the polyhedron structure. 0 if the point is not on the polyhedron.
      integer :: point_id = 0
      !> If point_id == 0, contains the edge number in the polyhedron structure. edge = 0 if point_id ≠ 0.
      integer :: edge = 0
      !> Contains the face number in the polyhedron structure that lies between this item and the next.
      integer :: face = 0
   end type t_chained_polygon_item

   !> Structure to track the state of the current algorithm.
   !! It stores the status of the cap polygon.
   !!
   !! @verbatim
   !!          5    f4    4              Legend:
   !!           *--------*  f3           ------
   !!       f5 /          \
   !!         /            * 3           ×   point of the polyhedron
   !!      6 *            /              *   point on an edge
   !!     f6 |           / f2            n   point n
   !!        ×----------*                fn  face between point n and n+1
   !!      1      f1     2
   !! @endverbatim
   !! @ingroup flood_polyhedron
   type t_chained_polygon
      !> Pointer to the last point
      type(t_chained_polygon_item), dimension(:), allocatable :: point
      !> Number of points in the chained polygon
      integer :: nb_points = 0
   contains
      procedure :: add_point => cg3_chained_polygon_add_point
   end type t_chained_polygon

   public :: t_chained_polygon_item, t_chained_polygon
   public :: cg3_flood_polyhedron, cg3_flood_polyhedron_centroid

contains

   !> Flood algorithm based on Diot & Francois 2016 publication
   !!
   !! Reference:
   !!  - Diot, S., & François, M. M. (2016). An interface reconstruction method based on an analytical formula for 3D arbitrary
   !!    convex cells. Journal of Computational Physics, 305, 63-74. doi:[10.1016/j.jcp.2015.10.011]
   !!
   !! @param[in]  polyhedron: any convex polyhedron. Must be fully initialized.
   !! @param[in]  normal: flood direction
   !! @param[in]  volume: prescribed volume to flood
   !! @param[out] polyhedron_full: part of @p polyhedron flooded with prescribed @p volume
   !! @param[out] polyhedron_empty: empty part of @p polyhedron
   !!
   !! [10.1016/j.jcp.2015.10.011]: https://doi.org/10.1016/j.jcp.2015.10.011 "10.1016/j.jcp.2015.10.011"
   !! @ingroup flood_polyhedron
   pure subroutine cg3_flood_polyhedron(polyhedron, normal, volume, polyhedron_full, polyhedron_empty, polygon)
      type(t_polyhedron), intent(in) :: polyhedron
      double precision, dimension(3), intent(in) :: normal
      double precision, intent(in) :: volume
      type(t_polyhedron), intent(out) :: polyhedron_full
      type(t_polyhedron), intent(out), optional :: polyhedron_empty
      type(t_chained_polygon), intent(out), optional :: polygon

      double precision, parameter :: DISTANCE_THRESHOLD = 1d8*epsilon(1d0)
      double precision, parameter :: NORMAL_THRESHOLD = 1d-7
      double precision, parameter :: PI = 2d0*acos(0d0)
      double precision, parameter :: NEWTON_EPSILON = 10d0*epsilon(1d0)
      integer, parameter :: ITER_NEWTON = 100

      ! Variables for ordering points along the n-axis
      type(t_incidence_matrix), dimension(:), allocatable :: point_to_edge
      double precision, dimension(polyhedron%nb_points) :: n_distance
      integer, dimension(polyhedron%nb_points) :: n_point, inv_n_point, inv_n_group
      integer, dimension(polyhedron%nb_points+1) :: n_group
      integer :: nb_planes

      ! Variables for the flood algorithm
      type(t_chained_polygon) :: bottom_polygon, top_polygon, cap_polygon
      double precision, dimension(3) :: coords, surface_vector
      double precision :: height, height_old, total_volume
      double precision :: A, B, C, volume_prismatoid
      integer :: nb_points_top

      ! Variables for creating full and empty polyhedron.
      ! 'cut_face' is true if the face is cut by the cap polygon.
      ! 'full_face' is true if the face is striclty below the cap polygon.
      logical, dimension(polyhedron%nb_faces) :: cut_face, full_face
      logical :: last_up, cap_is_cut, cap_is_bottom
      integer :: first, last, first_empty, last_empty, nb_new_points_cap, nb_faces_empty, nb_points_empty, empty_shift, full_padding

      ! Other variables
      integer :: i, j, k, point, edge, face, nb_edges, nb_faces, nb_points, curp, prevp, cur_full, cur_empty

      ! Algorithm:
      ! =========
      ! The complete algorithm contains 7 steps:
      ! 1 - Compute the distances of the planes
      ! 2 - Initialization of the first polygon
      ! 3 - Bracket the solution
      ! 4 - Compute cap polygon
      ! 5 - Compute full and empty polyhedron
      !-------------------------

      !-------------------------
      ! 1 - Compute the distances of the planes

      ! Summary: This step is the most sensible part of the algorithm. It is the only step that compares real numbers.
      ! =======  In this step, the point are sorted along the n-axis. If many points have the same n-distance, they
      !          are stored in the same groupolyp. nb_planes groups are generated during this step.
      !
      !          Three arrays are created during this step:
      !          1 - n_distance [double precision (nb_points)]:
      !                 contains the z-coordinate of every points.
      !          2 - n_point [integer (nb_points)]:
      !                 maps the sorted point list to the polyhedron point list.
      !          3 - n_group [integer (nb_planes+1)]:
      !                 n_group(k) contains the index of the first sorted point that belongs to the plane k.
      !
      !          To facilitate the inverse correspondance, two other inverse arrays are created:
      !          1 - inv_n_point [integer (nb_points)]:
      !                 associates a point of the polyhedron list to the index of the sorted list.
      !          2 - inv_n_group [integer (nb_points)]:
      !                 associates a point of the polyhedron list to its group number.

      ! Copy 'point_to_edge' incidence matrix
      allocate(point_to_edge, source=polyhedron%point_to_edge)

      ! Compute the distance of the points along the n-axis
      do i = 1, polyhedron%nb_points
         n_point(i) = i
         n_distance(i) = dot_product(polyhedron%point(:,i), normal)
      end do

      ! Sort the distance of the points along the n-axis
      call cg3_quicksort_n_distance(n_distance, n_point, 1, polyhedron%nb_points)

      ! Find the set of points with same distance
      ! The first index in n_point of the set of point k is equal to n_group(k)
      ! The last index in n_point of the set of point k is equal to n_group(k+1)-1
      ! This representation is inspired from CSR matrix storage
      nb_planes = 1
      n_group = 0
      n_group(nb_planes) = 1
      do i = 2, polyhedron%nb_points
         if (abs(n_distance(i) - n_distance(i-1)) > DISTANCE_THRESHOLD) then
            nb_planes = nb_planes + 1
            n_group(nb_planes) = i
         end if
      end do
      n_group(nb_planes+1) = polyhedron%nb_points + 1

      ! Test if the first three points belong to the same face
      ! If the face normal is approximately collinear to the flood direction,
      ! re-structure the order of the points.
      face = cg3_polyhedron_find_common_face_3_points(polyhedron, point_to_edge, n_point(1), n_point(2), n_point(3))
      if (face > 0) then
         if (norm2(cg3_cross_product(polyhedron%normal(:,face), normal)) < NORMAL_THRESHOLD) then
            ! Group all the points of the bottom face in the group 1
            n_group(2) = polyhedron%face(face)%size + 1
            n_point(1:n_group(2)-1) = polyhedron%face(face)%id

            nb_planes = 1
            do i = n_group(2), polyhedron%nb_points
               if (abs(n_distance(i) - n_distance(i-1)) > DISTANCE_THRESHOLD) then
                  nb_planes = nb_planes + 1
                  n_group(nb_planes) = i
               end if
            end do
            n_group(nb_planes+1) = polyhedron%nb_points + 1
         end if
      end if

      ! Create inverse tables
      do i = 1, nb_planes
         do j = n_group(i), n_group(i+1) - 1
            inv_n_point(n_point(j)) = j
            inv_n_group(n_point(j)) = i
         end do
      end do

      !-------------------------
      ! 2 - Initialization of the first polygon

      ! Summary: This step initialize the polygon defined as the intersection of the first plane and the polyhedron.
      ! =======  The polygon is stored as a list. Each item of the list contains information about
      !          neighbor edges and faces. Each item corresponds to a point of the polygon. The coordinates of this
      !          point is stored in the structure. There are two kinds of points: there are points that belongs to
      !          the polyhedron -- in such a way, the corresponding index is stored in the structure -- and there are
      !          points that are the intersection of the plane and an edge of the polyhedron -- the corresponding
      !          edge index is stored in the structure. The face above the edge joining the point of the previous
      !          item to the point of the current item is stored in the structure. The set of edges that start from
      !          the point of the current item and that go up is stored in the structure.

      !          To build the first polygon, there are 3 configurations to consider:
      !
      !          1 - The first group contains only one point. The polygon is a single point. No faces are associated.
      !          2 - The first group contains two points. The polygon is a single edge. Two faces are stored.
      !          3 - The first group contains >2 points. The polygon is a face of the polyhedron.
      !
      !          During this step, we begin to track the faces of the polygon that are fully or partially flooded.
      !          Two arrays of logical are dedicated to this task: full_face and cut_face. They are indexed with
      !          the numbering of the polyhedron faces. Note that when cut_face is true, full_face is also true for
      !          a given face.
      !
      !          The first polygon is stored in the top_polygon variable.

      full_face = .false.
      cut_face = .false.

      nb_points_top = n_group(2) - n_group(1)

      ! Allocate memory for top_polygon and top_polygon
      ! The maximal size for top_polygon is equal to the number of faces (worst case: quadrilateral section of a tetrahedron).
      allocate(top_polygon%point(polyhedron%nb_faces))
      allocate(bottom_polygon%point(polyhedron%nb_faces))

      if (nb_points_top == 1) then
         ! The top polygon is a point, there is no face associated
         call top_polygon%add_point(polyhedron%point(:,n_point(1)), n_point(1), 0, 0)
      else if (nb_points_top == 2) then
         ! The top polygon is a edge

         ! Find the common edge
         edge = cg3_polyhedron_find_common_edge(polyhedron, point_to_edge, n_point(1), n_point(2))

         call top_polygon%add_point(polyhedron%point(:,n_point(1)), n_point(1), 0, polyhedron%edge_to_face(1,edge))
         call top_polygon%add_point(polyhedron%point(:,n_point(2)), n_point(2), 0, polyhedron%edge_to_face(2,edge))

         full_face(polyhedron%edge_to_face(1,edge)) = .true.
         full_face(polyhedron%edge_to_face(2,edge)) = .true.
      else
         ! The top polygon is a face
         ! Find the common face
         face = cg3_polyhedron_find_common_face_3_points(polyhedron, point_to_edge, n_point(1), n_point(2), n_point(3))

         full_face(face) = .true.

         ! Copy the points of the face to ensure the right order
         do i = 1, nb_points_top
            n_point(i) = polyhedron%face(face)%id(i)
         end do

         ! Create the polygon with the point of the face in reverse order to get a polygon with its points
         ! oriented in counter-clockwise order when watched from above.
         do i = nb_points_top, 1, -1
            point = polyhedron%face(face)%id(i)

            if (i > 1) then
               edge = cg3_polyhedron_find_common_edge(polyhedron, point_to_edge, n_point(i), n_point(i-1))
            else
               edge = cg3_polyhedron_find_common_edge(polyhedron, point_to_edge, n_point(1), n_point(nb_points_top))
            end if

            if (polyhedron%edge_to_face(1,edge) == face) then
               call top_polygon%add_point(polyhedron%point(:,point), point, 0, polyhedron%edge_to_face(2,edge))
               full_face(polyhedron%edge_to_face(2,edge)) = .true.
            else
               call top_polygon%add_point(polyhedron%point(:,point), point, 0, polyhedron%edge_to_face(1,edge))
               full_face(polyhedron%edge_to_face(1,edge)) = .true.
            end if
         end do
      end if

      ! Create and order the set of edges that go up passing through a point
      do i = 1, top_polygon%nb_points
         call cg3_create_edge_set(polyhedron, point_to_edge, top_polygon%point(i), inv_n_group, 1)
      end do

      ! Register all the faces around the first point as full face
      do curp = 1, top_polygon%nb_points
         do i = 1, size(top_polygon%point(curp)%edge_set)
            full_face(polyhedron%edge_to_face(1,top_polygon%point(curp)%edge_set(i))) = .true.
            full_face(polyhedron%edge_to_face(2,top_polygon%point(curp)%edge_set(i))) = .true.
         end do
      end do

      !-------------------------
      ! 3 - Bracket the solution

      ! Summary: Find the plane such that the volume below it is equal to the prescribed volume.
      ! =======  This step consists in bracketting this plane between two planes of the sorted points.
      !          The shape contained between two planes is known as a prismatoid.
      !
      !          The two polygons bracketting the plane are stored in top_polygon and top_polygon.
      !
      !          Before constructing the top polygon, we check if the volume of the prismatoid above
      !          the top polygon + the volume below the top polygon is greater than the prescribed
      !          volume. Note that we do not need to construct the top polygon to compute the volume!
      !          If the volume exceed the prescribed volume, the final plane is found using the fact that
      !          the evolution of the volume in a prismatoid is a cubic function of the n-distance. A Newton-
      !          Raphson algorithm is used to find the correct root. If the volume is smaller that the
      !          prescribed volume, the top polygon is constructed and we start again this step from
      !          the top polygon.

      total_volume = 0d0
      height = 0d0
      A = 0d0; B = 0d0; C = 0d0
      cap_is_cut = .false.
      cap_is_bottom = .false.

      ! Loop over all the planes
      do k = 2, nb_planes
         ! Compute the height of the prismatoid
         height = n_distance(n_group(k)) - n_distance(n_group(k-1))

         ! Compute the volume of the prismatoid
         associate(top_point => top_polygon%point)
            A = 0d0; B = 0d0; C = 0d0
            surface_vector = 0d0

            prevp = top_polygon%nb_points
            do curp = 1, top_polygon%nb_points
               ! Tetrahedron contributions
               do i = 1, size(top_point(curp)%edge_set) - 1
                  A = A + dot_product(cg3_cross_product(top_point(curp)%tangent(:,i), top_point(curp)%tangent(:,i+1)), normal) &
                     &  / (dot_product(top_point(curp)%tangent(:,i), normal)*dot_product(top_point(curp)%tangent(:,i+1), normal))
               end do

               ! Compute the area of the top polygon
               if (curp > 2) then
                  surface_vector = surface_vector + cg3_cross_product(top_point(prevp)%point - top_point(1)%point, &
                     &                                                top_point(curp)%point  - top_point(1)%point  )
               end if

               prevp = curp
            end do

            C = norm2(surface_vector)/2d0

            ! Wedge contribution
            if (top_polygon%nb_points > 1) then
               prevp = top_polygon%nb_points
               do curp = 1, top_polygon%nb_points
                  B = B - norm2(top_point(curp)%point - top_point(prevp)%point)           &
                     &  * dot_product(polyhedron%normal(:,top_point(prevp)%face), normal) &
                     &  / norm2(cg3_cross_product(normal, polyhedron%normal(:,top_point(prevp)%face)))

                  nb_edges = size(top_point(prevp)%edge_set)
                  A = A + dot_product(cg3_cross_product(top_point(prevp)%tangent(:,nb_edges),  &
                     &                                  top_point(curp)%tangent(:,1)), normal) &
                     &  / (dot_product(top_point(prevp)%tangent(:,nb_edges), normal)           &
                     &  * dot_product(top_point(curp)%tangent(:,1), normal))

                  prevp = curp
               end do
            else
               ! Remaining tetrahedron
               A = A + dot_product(cg3_cross_product(top_point(1)%tangent(:,size(top_point(1)%edge_set)), &
                  &                                  top_point(1)%tangent(:,1)), normal)                  &
                  &  / (dot_product(top_point(1)%tangent(:,size(top_point(1)%edge_set)), normal)          &
                  &  * dot_product(top_point(1)%tangent(:,1), normal))
            end if
         end associate

         volume_prismatoid = height*(C + height*(B/2d0 + height*A/6d0))

         ! Check if we have bracketed the volume between two planes
         if (total_volume + volume_prismatoid >= volume) then
            ! Newton-Raphson algorithm to find the final height.
            ! The algorithm starts from the middle of the prismatoid
            height_old = 0.5d0*height
            do i = 1, ITER_NEWTON
               ! Compute the derivative
               height = (C + height_old*(B + height_old*A/2d0))

               ! Check for zero derivative
               if (abs(height) < NEWTON_EPSILON) exit

               height = height_old - (total_volume - volume + height_old*(C + height_old*(B/2d0 + height_old*A/6d0))) &
                  &   / height

               ! Check for convergence
               if (abs(height - height_old) < NEWTON_EPSILON) exit

               height_old = height
            end do

            if (abs(height) <= NEWTON_EPSILON) then
               ! The bottom polygon becomes the cap polygon
               call cg3_move_alloc_chained_polygon(top_polygon, cap_polygon)
               ! Set bottom flag
               cap_is_bottom = .true.
               ! Set flag
               cap_is_cut = .true.
               ! Exit the loop
               exit
            else if (abs(n_distance(n_group(k)) - n_distance(n_group(k-1)) - height) <= NEWTON_EPSILON) then
               ! Set flag
               cap_is_cut = .true.
               ! Compute the height of the prismatoid
               height = n_distance(n_group(k)) - n_distance(n_group(k-1))
               ! Construct the next plane
            else
               ! Exit the loop
               exit
            end if
         end if

         ! The top polygon becomes the bottom polygon!
         call cg3_move_alloc_chained_polygon(top_polygon, bottom_polygon)

         total_volume = total_volume + volume_prismatoid

         ! Compute the next polygon
         ! Loop on every point of the bottom polygon
         p_loop: do curp = 1, bottom_polygon%nb_points
            ! Check if the current point is a point of the polyhedron
            if (bottom_polygon%point(curp)%edge > 0) then
               ! The current point is not a polyhedron point => we are on an edge of the polyhedron
               ! Check if the end of the edge belongs to the next plane
               point = polyhedron%edge(1,bottom_polygon%point(curp)%edge)
               if (inv_n_group(point) == k) then
                  call top_polygon%add_point(polyhedron%point(:,point), point, 0, 0)
                  cycle p_loop
               end if

               point = polyhedron%edge(2,bottom_polygon%point(curp)%edge)
               if (inv_n_group(point) == k) then
                  call top_polygon%add_point(polyhedron%point(:,point), point, 0, 0)
                  cycle p_loop
               end if

               ! The end of the edge does not belong the the next plane
               ! => compute the coordinates of the plane/edge intersection
               coords = bottom_polygon%point(curp)%point                                                   &
                  &   + height/dot_product(polyhedron%tangent(:, bottom_polygon%point(curp)%edge), normal) &
                  &   * polyhedron%tangent(:, bottom_polygon%point(curp)%edge)

               ! Add the intersection to the top polygon
               call top_polygon%add_point(coords, 0, bottom_polygon%point(curp)%edge, bottom_polygon%point(curp)%face)

               ! Add edge information
               allocate(top_polygon%point(top_polygon%nb_points)%edge_set(1))
               allocate(top_polygon%point(top_polygon%nb_points)%edge_point(1))
               allocate(top_polygon%point(top_polygon%nb_points)%tangent(3,1))

               top_polygon%point(top_polygon%nb_points)%edge_set = bottom_polygon%point(curp)%edge
               top_polygon%point(top_polygon%nb_points)%edge_point = 0
               if (dot_product(polyhedron%tangent(:,bottom_polygon%point(curp)%edge), normal) > 0) then
                  top_polygon%point(top_polygon%nb_points)%tangent(:,1) = polyhedron%tangent(:,bottom_polygon%point(curp)%edge)
               else
                  top_polygon%point(top_polygon%nb_points)%tangent(:,1) = -polyhedron%tangent(:,bottom_polygon%point(curp)%edge)
               end if
            else
               ! The current point is a polyhedron point
               point = bottom_polygon%point(curp)%point_id

               ! Loop on every edges of the current point that go up
               nb_edges = size(bottom_polygon%point(curp)%edge_set)
               j_loop: do j = 1, nb_edges
                  ! Check if the extreme point of the edge belongs to the next plane
                  if (inv_n_group(bottom_polygon%point(curp)%edge_point(j)) == k) then
                     call top_polygon%add_point(polyhedron%point(:,bottom_polygon%point(curp)%edge_point(j)), &
                        &                       bottom_polygon%point(curp)%edge_point(j), 0, 0)
                     cycle j_loop
                  end if

                  ! The extreme point of the edge does not belong to the next plane
                  ! => compute the coordinates of the plane/edge intersection
                  coords = bottom_polygon%point(curp)%point                                    &
                     &   + height/dot_product(bottom_polygon%point(curp)%tangent(:,j), normal) &
                     &   * bottom_polygon%point(curp)%tangent(:,j)

                  ! Add the intersection to the top polygon
                  call top_polygon%add_point(coords, 0, bottom_polygon%point(curp)%edge_set(j), 0)

                  ! Add edge information
                  allocate(top_polygon%point(top_polygon%nb_points)%edge_set(1))
                  allocate(top_polygon%point(top_polygon%nb_points)%edge_point(1))
                  allocate(top_polygon%point(top_polygon%nb_points)%tangent(3,1))

                  top_polygon%point(top_polygon%nb_points)%edge_set = bottom_polygon%point(curp)%edge_set(j)
                  top_polygon%point(top_polygon%nb_points)%edge_point = 0
                  if (dot_product(polyhedron%tangent(:,bottom_polygon%point(curp)%edge_set(j)), normal) > 0) then
                     top_polygon%point(top_polygon%nb_points)%tangent(:,1) = &
                        &  polyhedron%tangent(:,bottom_polygon%point(curp)%edge_set(j))
                  else
                     top_polygon%point(top_polygon%nb_points)%tangent(:,1) = &
                        &  -polyhedron%tangent(:,bottom_polygon%point(curp)%edge_set(j))
                  end if
               end do j_loop
            end if
         end do p_loop

         associate(top_point => top_polygon%point)
            ! Order the edge that go up passing through the points of the next planes
            do i = 1, top_polygon%nb_points
               if (top_point(i)%edge == 0) call cg3_create_edge_set(polyhedron, point_to_edge, top_point(i), inv_n_group, k)
            end do

            ! Find the faces between the points
            prevp = top_polygon%nb_points
            do curp = 1, top_polygon%nb_points
               ! If multiple edge goes up from one point, mark all faces between as full faces
               if (top_point(curp)%point_id > 0) then
                  do i = 1, size(top_point(curp)%edge_set) - 1
                     full_face(cg3_polyhedron_find_common_face_2_edges(polyhedron, top_point(curp)%edge_set(i), &
                        &                                              top_point(curp)%edge_set(i+1))) = .true.
                  end do
               end if

               if (top_point(prevp)%face > 0) then
                  full_face(top_point(prevp)%face) = .true.
                  prevp = curp
                  cycle
               end if

               ! Book-keeping algorithm… I am sure that we can optimize this part that travels connectivities to find a face.
               if (top_point(prevp)%edge > 0) then
                  if (top_point(curp)%edge > 0) then
                     top_point(prevp)%face = cg3_polyhedron_find_common_face_2_edges(polyhedron, top_point(prevp)%edge, &
                        &                                                            top_point(curp)%edge               )
                  else
                     top_point(prevp)%face = cg3_polyhedron_find_common_face_edge_point(polyhedron, point_to_edge, &
                        &                                                               top_point(prevp)%edge,     &
                        &                                                               top_point(curp)%point_id   )
                  end if
               else if (top_point(curp)%edge > 0) then
                  top_point(prevp)%face = cg3_polyhedron_find_common_face_edge_point(polyhedron, point_to_edge, &
                     &                                                               top_point(curp)%edge,      &
                     &                                                               top_point(prevp)%point_id  )
               else ! top_point(prevp)%point_id /= top_point(curp)%point_id
                  top_point(prevp)%face = cg3_polyhedron_find_common_face_2_points_up(polyhedron, point_to_edge,              &
                     &                                                                top_point(prevp)%point_id,              &
                     &                                                                top_point(curp)%point_id, inv_n_group, k)
               end if

               full_face(top_point(prevp)%face) = .true.

               prevp = curp
            end do
         end associate

         if (cap_is_cut) then
            ! The bottom polygon becomes the cap polygon
            call cg3_move_alloc_chained_polygon(top_polygon, cap_polygon)
            exit
         end if
      end do

      if (cap_is_cut .and. cap_is_bottom) k = k - 1

      !-------------------------
      ! 4 - Compute cap polygon

      ! Summary: This step consists in constructing the cap polygon that corresponds to the intersection of the
      ! =======  final plane and the polyhedron. This steps re-uses some code used in the previous step, but it
      !          is adaptated to the current situtation.
      !
      !          This step tracks the faces that are cut by the final plane.

      if (.not. cap_is_cut) then
         allocate(cap_polygon%point(polyhedron%nb_faces))

         do curp = 1, top_polygon%nb_points
            ! Check if the current point is a point of the polyhedron
            if (top_polygon%point(curp)%edge > 0) then
               ! The current point is not a point of the polyhedron
               ! compute the coordinates of the plane/edge intersection
               coords = top_polygon%point(curp)%point                                                   &
                  &   + height/dot_product(polyhedron%tangent(:, top_polygon%point(curp)%edge), normal) &
                  &   * polyhedron%tangent(:, top_polygon%point(curp)%edge)

               ! Add the intersection to the top polygon
               call cap_polygon%add_point(coords, 0, top_polygon%point(curp)%edge, top_polygon%point(curp)%face)

               ! Add edge information
               allocate(cap_polygon%point(cap_polygon%nb_points)%edge_set(1))
               allocate(cap_polygon%point(cap_polygon%nb_points)%edge_point(1))
               allocate(cap_polygon%point(cap_polygon%nb_points)%tangent(3,1))

               cap_polygon%point(cap_polygon%nb_points)%edge_set = top_polygon%point(curp)%edge
               cap_polygon%point(cap_polygon%nb_points)%edge_point = 0
               if (dot_product(polyhedron%tangent(:,top_polygon%point(curp)%edge), normal) > 0) then
                  cap_polygon%point(cap_polygon%nb_points)%tangent(:,1) = polyhedron%tangent(:,top_polygon%point(curp)%edge)
               else
                  cap_polygon%point(cap_polygon%nb_points)%tangent(:,1) = -polyhedron%tangent(:,top_polygon%point(curp)%edge)
               end if
            else
               ! The current point is a point of the polyhedron
               point = top_polygon%point(curp)%point_id

               ! Loop on every edges of the current point that go up
               nb_edges = size(top_polygon%point(curp)%edge_set)
               do j = 1, nb_edges
                  ! The extreme point of the edge does not belong to the next plane
                  ! => compute the coordinates of the plane/edge intersection
                  coords = top_polygon%point(curp)%point                                    &
                     &   + height/dot_product(top_polygon%point(curp)%tangent(:,j), normal) &
                     &   * top_polygon%point(curp)%tangent(:,j)

                  ! Add the intersection to the top polygon
                  call cg3_chained_polygon_add_point(cap_polygon, coords, 0, top_polygon%point(curp)%edge_set(j), 0)

                  ! Add edge information
                  allocate(cap_polygon%point(cap_polygon%nb_points)%edge_set(1))
                  allocate(cap_polygon%point(cap_polygon%nb_points)%edge_point(1))
                  allocate(cap_polygon%point(cap_polygon%nb_points)%tangent(3,1))

                  cap_polygon%point(cap_polygon%nb_points)%edge_set = top_polygon%point(curp)%edge_set(j)
                  cap_polygon%point(cap_polygon%nb_points)%edge_point = 0
                  if (dot_product(polyhedron%tangent(:,top_polygon%point(curp)%edge_set(j)), normal) > 0) then
                     cap_polygon%point(cap_polygon%nb_points)%tangent(:,1) = &
                        & polyhedron%tangent(:,top_polygon%point(curp)%edge_set(j))
                  else
                     cap_polygon%point(cap_polygon%nb_points)%tangent(:,1) = &
                        & -polyhedron%tangent(:,top_polygon%point(curp)%edge_set(j))
                  end if
               end do
            end if
         end do

         ! Find the faces between the points of the cap polygon
         ! Detect also the cut faces
         ! Give a number to point without id
         associate(cap_point => cap_polygon%point)
            prevp = cap_polygon%nb_points
            do curp = 1, cap_polygon%nb_points
               cap_point(curp)%point_id = -curp

               ! If multiple edge goes up from one point, mark all faces between as full faces
               if (cap_point(curp)%point_id > 0) then
                  do i = 1, size(cap_point(curp)%edge_set) - 1
                     full_face(cg3_polyhedron_find_common_face_2_edges(polyhedron, cap_point(curp)%edge_set(i), &
                        &                                              cap_point(curp)%edge_set(i+1))) = .true.
                  end do
               end if

               if (cap_point(prevp)%face > 0) then
                  full_face(cap_point(prevp)%face) = .true.
                  cut_face(cap_point(prevp)%face) = .true.
                  prevp = curp
                  cycle
               end if

               if (cap_point(prevp)%edge > 0) then
                  if (cap_point(curp)%edge > 0) then
                     cap_point(prevp)%face = cg3_polyhedron_find_common_face_2_edges(polyhedron, cap_point(prevp)%edge, &
                        &                                                            cap_point(curp)%edge)
                  else
                     cap_point(prevp)%face = cg3_polyhedron_find_common_face_edge_point(polyhedron, point_to_edge, &
                        &                                                               cap_point(prevp)%edge,     &
                        &                                                               cap_point(curp)%point_id   )
                  end if
               else if (cap_point(curp)%edge > 0) then
                  cap_point(prevp)%face = cg3_polyhedron_find_common_face_edge_point(polyhedron, point_to_edge, &
                     &                                                               cap_point(curp)%edge,      &
                     &                                                               cap_point(prevp)%point_id  )
               else ! cap_point(prevp)%point_id /= cap_point(curp)%point_id
                  cap_point(prevp)%face = cg3_polyhedron_find_common_face_2_points_up(polyhedron, point_to_edge, &
                     &                                                                cap_point(prevp)%point_id, &
                     &                                                                cap_point(curp)%point_id, inv_n_group, k-1)
               end if

               full_face(cap_point(prevp)%face) = .true.
               cut_face(cap_point(prevp)%face) = .true.

               prevp = curp
            end do
         end associate

         ! The number of new points in the polyhedron is the number of points of the cap polygon.
         nb_new_points_cap = cap_polygon%nb_points
      else
         ! Detect cut faces and give a number to point without id
         prevp = cap_polygon%nb_points
         point = 1
         do curp = 1, cap_polygon%nb_points
            if (cap_polygon%point(curp)%point_id == 0) then
               cut_face(cap_polygon%point(prevp)%face) = .true.
               cap_polygon%point(curp)%point_id = -point
               point = point + 1
               prevp = curp
               cycle
            else
               ! If multiple edge goes up from one point, do not mark the faces beetween as cut faces or full face
               do i = 1, size(cap_polygon%point(curp)%edge_set) - 1
                  face = cg3_polyhedron_find_common_face_2_edges(polyhedron, cap_polygon%point(curp)%edge_set(i), &
                     &                                           cap_polygon%point(curp)%edge_set(i+1))
                  full_face(face) = .false.
                  cut_face(face)  = .false.
               end do
            end if

            if (cap_polygon%point(prevp)%point_id == 0) then
               cut_face(cap_polygon%point(prevp)%face) = .true.
               prevp = curp
               cycle
            end if

            ! Detect if there is a common edge shared by the two points
            edge = cg3_polyhedron_find_common_edge(polyhedron, point_to_edge, cap_polygon%point(curp)%point_id, &
               &                                   cap_polygon%point(prevp)%point_id)

            ! If there is a common edge, the face is not cut
            if (edge == 0) then
               cut_face(cap_polygon%point(prevp)%face) = .true.
            else
               full_face(cap_polygon%point(prevp)%face) = .false.
               cut_face(cap_polygon%point(prevp)%face)  = .false.
            end if

            prevp = curp
         end do

         ! The number of new points in the polyhedron is the number of points of the cap polygon minus the size of group k.
         nb_new_points_cap = cap_polygon%nb_points - n_group(k+1) + n_group(k)
      end if

      ! Output cap polygon if required
      if (present(polygon)) polygon = cap_polygon

      !-------------------------
      ! 5 - Compute full and empty polyhedron

      ! Summary: This step splits the polyherdon in two polyhedrons: a polyhedron flooded and an empty polyhedron.
      ! =======

      ! This variable corresponds to the number of points in the empty polyhedron without the new points of the cap
      empty_shift = polyhedron%nb_points - n_group(k) + 1

      ! Last index of the full polyhedron that contains a point of the original polyhedron
      full_padding = n_group(k) - 1 + cap_polygon%nb_points - nb_new_points_cap

      ! Allocate memory for full and empty polyhedron
      polyhedron_full%nb_points = n_group(k) - 1 + cap_polygon%nb_points
      allocate(polyhedron_full%point(3,polyhedron_full%nb_points))

      if (present(polyhedron_empty)) then
         polyhedron_empty%nb_points = empty_shift + nb_new_points_cap
         allocate(polyhedron_empty%point(3,polyhedron_empty%nb_points))
      end if

      ! Copy the points that belongs to the full polyhedron
      do i = 1, full_padding
         polyhedron_full%point(:,i) = polyhedron%point(:,n_point(i))
      end do

      ! Copy the points that belongs to the empty polyhedron
      if (present(polyhedron_empty)) then
         do i = n_group(k), polyhedron%nb_points
            polyhedron_empty%point(:,i-n_group(k)+1) = polyhedron%point(:,n_point(i))
         end do
      end if

      ! Add cap points to the full and empty polyhedrons
      ! Counter for full polyhedron
      cur_full = full_padding + 1
      ! Counter for empty polyhedron
      cur_empty = empty_shift + 1
      do i = 1, cap_polygon%nb_points
         if (cap_polygon%point(i)%point_id > 0) cycle
         ! Full polyhedron
         polyhedron_full%point(:,cur_full) = cap_polygon%point(i)%point
         cur_full = cur_full + 1
         ! Empty polyhedron
         if (present(polyhedron_empty)) then
            polyhedron_empty%point(:,cur_empty) = cap_polygon%point(i)%point
            cur_empty = cur_empty + 1
         end if
      end do

      ! Count the number of faces (+1 for the cap)
      polyhedron_full%nb_faces = count(full_face) + 1
      if (present(polyhedron_empty)) polyhedron_empty%nb_faces = count(.not. full_face .or. cut_face) + 1

      ! Allocate memory for faces
      allocate(polyhedron_full%face(polyhedron_full%nb_faces))
      if (present(polyhedron_empty)) allocate(polyhedron_empty%face(polyhedron_empty%nb_faces))

      ! Fill faces of the full polyhedron
      nb_faces = 0
      nb_faces_empty = 0
      do i = 1, polyhedron%nb_faces
         if (.not. full_face(i)) then ! Face that belong to empty polyhedron
            if (present(polyhedron_empty)) then
               ! The current face is a non-cut face of the empty polyhedron
               nb_faces_empty = nb_faces_empty + 1

               allocate(polyhedron_empty%face(nb_faces_empty)%id(polyhedron%face(i)%size))
               polyhedron_empty%face(nb_faces_empty)%size = polyhedron%face(i)%size

               do j = 1, polyhedron_empty%face(nb_faces_empty)%size
                  ! This relation is just magic!
                  polyhedron_empty%face(nb_faces_empty)%id(j) = inv_n_point(polyhedron%face(i)%id(j)) - n_group(k) + 1
               end do
            end if
         else if (cut_face(i)) then ! Cut face
            nb_faces = nb_faces + 1
            nb_faces_empty = nb_faces_empty + 1

            ! Count the number of point
            ! Find the first point below the plane such that the next point is above (last)
            ! Find the first point below the plane such that the previous point is above (first)
            nb_points = 2
            first = 0
            last = 0
            last_up = inv_n_group(polyhedron%face(i)%id(polyhedron%face(i)%size)) > k-1
            do j = 1, polyhedron%face(i)%size
               if (inv_n_group(polyhedron%face(i)%id(j)) > k-1) then
                  if (last == 0 .and. .not. last_up) last = j - 1

                  last_up = .true.

                  cycle
               end if

               nb_points = nb_points + 1

               if (first == 0 .and. last_up) first = j
               last_up = .false.
            end do

            if (last == 0) last = polyhedron%face(i)%size

            ! Deduce the number of points for the face of the empty polyhedron (temporary, see corrections below)
            nb_points_empty = polyhedron%face(i)%size - nb_points + 4

            ! Allocate memory for the points of the current fulle face
            allocate(polyhedron_full%face(nb_faces)%id(nb_points))
            polyhedron_full%face(nb_faces)%size = nb_points

            ! Find the edge of the cap polygon which corresponds to the current face
            do point = 1, cap_polygon%nb_points
               if (cap_polygon%point(point)%face == i) exit
            end do

            ! Add the first point to the full face
            nb_points = 1
            if (cap_polygon%point(point)%point_id < 0) then
               polyhedron_full%face(nb_faces)%id(nb_points) = full_padding - cap_polygon%point(point)%point_id
               last_empty = first - 1
            else
               polyhedron_full%face(nb_faces)%id(nb_points) = inv_n_point(cap_polygon%point(point)%point_id)
               ! Correction of the number of points for the empty cell
               nb_points_empty = nb_points_empty - 1
               last_empty = first - 1
               if (last_empty < 1) last_empty = polyhedron%face(i)%size
               last_empty = last_empty - 1
            end if

            if (last_empty < 1) last_empty = polyhedron%face(i)%size

            ! Add points of the current face to the full face
            if (first <= last) then
               do j = first, last
                  nb_points = nb_points + 1
                  polyhedron_full%face(nb_faces)%id(nb_points) = inv_n_point(polyhedron%face(i)%id(j))
               end do
            else
               do j = first, polyhedron%face(i)%size
                  nb_points = nb_points + 1
                  polyhedron_full%face(nb_faces)%id(nb_points) = inv_n_point(polyhedron%face(i)%id(j))
               end do

               do j = 1, last
                  nb_points = nb_points + 1
                  ! Correction of the number of points for the empty cell
                  polyhedron_full%face(nb_faces)%id(nb_points) = inv_n_point(polyhedron%face(i)%id(j))
               end do
            end if

            ! Add the last point to the full face
            nb_points = nb_points + 1
            point = point + 1
            if (point > cap_polygon%nb_points) point = 1
            if (cap_polygon%point(point)%point_id < 0) then
               polyhedron_full%face(nb_faces)%id(nb_points) = full_padding - cap_polygon%point(point)%point_id
               first_empty = last + 1
            else
               polyhedron_full%face(nb_faces)%id(nb_points) = inv_n_point(cap_polygon%point(point)%point_id)
               nb_points_empty = nb_points_empty - 1
               first_empty = last + 1
               if (first_empty > polyhedron%face(i)%size) first_empty = 1
               first_empty = first_empty + 1
            end if

            if (first_empty > polyhedron%face(i)%size) first_empty = 1

            if (present(polyhedron_empty)) then
               ! Allocate memory for the points of the current empty face
               allocate(polyhedron_empty%face(nb_faces_empty)%id(nb_points_empty))
               polyhedron_empty%face(nb_faces_empty)%size = nb_points_empty

               ! Add the first point for the empty face
               nb_points_empty = 1
               if (cap_polygon%point(point)%point_id < 0) then
                  polyhedron_empty%face(nb_faces_empty)%id(nb_points_empty) = empty_shift - cap_polygon%point(point)%point_id
               else
                  polyhedron_empty%face(nb_faces_empty)%id(nb_points_empty) = inv_n_point(cap_polygon%point(point)%point_id) &
                     &                                                      - n_group(k) + 1
               end if

               ! Add points of the current face to the empty face
               if (first_empty <= last_empty) then
                  do j = first_empty, last_empty
                     nb_points_empty = nb_points_empty + 1
                     polyhedron_empty%face(nb_faces_empty)%id(nb_points_empty) = inv_n_point(polyhedron%face(i)%id(j))-n_group(k)+1
                  end do
               else
                  do j = first_empty, polyhedron%face(i)%size
                     nb_points_empty = nb_points_empty + 1
                     polyhedron_empty%face(nb_faces_empty)%id(nb_points_empty) = inv_n_point(polyhedron%face(i)%id(j))-n_group(k)+1
                  end do

                  do j = 1, last_empty
                     nb_points_empty = nb_points_empty + 1
                     polyhedron_empty%face(nb_faces_empty)%id(nb_points_empty) = inv_n_point(polyhedron%face(i)%id(j))-n_group(k)+1
                  end do
               end if
            end if

            ! Add the last point to the empty face
            point = point - 1
            if (point < 1) point = cap_polygon%nb_points
            if (present(polyhedron_empty)) then
               nb_points_empty = nb_points_empty + 1
               if (cap_polygon%point(point)%point_id < 0) then
                  polyhedron_empty%face(nb_faces_empty)%id(nb_points_empty) = empty_shift - cap_polygon%point(point)%point_id
               else
                  polyhedron_empty%face(nb_faces_empty)%id(nb_points_empty) = inv_n_point(cap_polygon%point(point)%point_id) &
                     &                                                      - n_group(k) + 1
               end if
            end if
         else ! Face that belong to full polyhedron
            nb_faces = nb_faces + 1

            allocate(polyhedron_full%face(nb_faces)%id(polyhedron%face(i)%size))
            polyhedron_full%face(nb_faces)%size = polyhedron%face(i)%size

            do j = 1, polyhedron_full%face(nb_faces)%size
               polyhedron_full%face(nb_faces)%id(j) = inv_n_point(polyhedron%face(i)%id(j))
            end do
         end if
      end do

      ! Fill the cap
      nb_faces = nb_faces + 1
      nb_faces_empty = nb_faces_empty + 1

      allocate(polyhedron_full%face(nb_faces)%id(cap_polygon%nb_points))
      polyhedron_full%face(nb_faces)%size = cap_polygon%nb_points
      if (present(polyhedron_empty)) then
         allocate(polyhedron_empty%face(nb_faces_empty)%id(cap_polygon%nb_points))
         polyhedron_empty%face(nb_faces_empty)%size = cap_polygon%nb_points
      end if

      do i = 1, cap_polygon%nb_points
         if (cap_polygon%point(i)%point_id < 0) then
            polyhedron_full%face(nb_faces)%id(i) = full_padding - cap_polygon%point(i)%point_id
            if (present(polyhedron_empty)) then
               ! Reverse order for the empty face to ensure an outgoing normal
               polyhedron_empty%face(nb_faces_empty)%id(cap_polygon%nb_points - i + 1) = empty_shift - cap_polygon%point(i)%point_id
            end if
         else
            polyhedron_full%face(nb_faces)%id(i) = inv_n_point(cap_polygon%point(i)%point_id)
            if (present(polyhedron_empty)) then
               ! Reverse order for the empty face to ensure an outgoing normal
               polyhedron_empty%face(nb_faces_empty)%id(cap_polygon%nb_points - i + 1) = inv_n_point(cap_polygon%point(i)%point_id)&
                  &                                                                    - n_group(k) + 1
            end if
         end if
      end do

      ! Thanks for reading this long code. Here is an ASCII potato:
      !              ***
      !             ****
      !            ****
      !            **
   end subroutine cg3_flood_polyhedron

   pure subroutine cg3_flood_polyhedron_centroid(polyhedron, normal, volume, centroid_full, polygon)
      type(t_polyhedron), intent(in) :: polyhedron
      double precision, dimension(3), intent(in) :: normal
      double precision, intent(in) :: volume
      double precision, dimension(3), intent(out) :: centroid_full
      type(t_chained_polygon), intent(out), optional :: polygon

      double precision, parameter :: DISTANCE_THRESHOLD = 1d8*epsilon(1d0)
      double precision, parameter :: NORMAL_THRESHOLD = 1d-7
      double precision, parameter :: PI = 2d0*acos(0d0)
      double precision, parameter :: NEWTON_EPSILON = 10d0*epsilon(1d0)
      integer, parameter :: ITER_NEWTON = 100

      ! Variables for ordering points along the n-axis
      type(t_incidence_matrix), dimension(:), allocatable :: point_to_edge
      double precision, dimension(polyhedron%nb_points) :: n_distance
      integer, dimension(polyhedron%nb_points) :: n_point, inv_n_point, inv_n_group
      integer, dimension(polyhedron%nb_points+1) :: n_group
      integer :: nb_planes

      ! Variables for the flood algorithm
      type(t_chained_polygon) :: bottom_polygon, top_polygon, cap_polygon
      double precision, dimension(3) :: coords, prism_momentum, edge_prism
      double precision :: height, area, polygon_area, height_old
      double precision :: A, B, C, volume_prismatoid, total_volume, volume_tetra, volume_prism, Bf, S
      integer :: nb_points_top, nb_new_points_cap

      ! Variables for computing the momentum
      logical :: cap_is_cut, cap_is_bottom
      double precision, dimension(:,:), allocatable :: momentum

      ! Other variables
      integer :: i, j, k, point, edge, face, nb_edges, curp, prevp

      ! Algorithm:
      ! =========
      ! The complete algorithm contains 7 steps:
      ! 1 - Compute the distances of the planes
      ! 2 - Initialization of the first polygon
      ! 3 - Bracket the solution
      !-------------------------

      !-------------------------
      ! 1 - Compute the distances of the planes

      ! Summary: This step is the most sensible part of the algorithm. It is the only step that compares real numbers.
      ! =======  In this step, the point are sorted along the n-axis. If many points have the same n-distance, they
      !          are stored in the same groupolyp. nb_planes groups are generated during this step.
      !
      !          Three arrays are created during this step:
      !          1 - n_distance [double precision (nb_points)]:
      !                 contains the z-coordinate of every points.
      !          2 - n_point [integer (nb_points)]:
      !                 maps the sorted point list to the polyhedron point list.
      !          3 - n_group [integer (nb_planes+1)]:
      !                 n_group(k) contains the index of the first sorted point that belongs to the plane k.
      !
      !          To facilitate the inverse correspondance, two other inverse arrays are created:
      !          1 - inv_n_point [integer (nb_points)]:
      !                 associates a point of the polyhedron list to the index of the sorted list.
      !          2 - inv_n_group [integer (nb_points)]:
      !                 associates a point of the polyhedron list to its group number.

      ! Copy 'point_to_edge' incidence matrix
      allocate(point_to_edge, source=polyhedron%point_to_edge)

      ! Compute the distance of the points along the n-axis
      do i = 1, polyhedron%nb_points
         n_point(i) = i
         n_distance(i) = dot_product(polyhedron%point(:,i), normal)
      end do

      ! Sort the distance of the points along the n-axis
      call cg3_quicksort_n_distance(n_distance, n_point, 1, polyhedron%nb_points)

      ! Find the set of points with same distance
      ! The first index in n_point of the set of point k is equal to n_group(k)
      ! The last index in n_point of the set of point k is equal to n_group(k+1)-1
      ! This representation is inspired from CSR matrix storage
      nb_planes = 1
      n_group = 0
      n_group(nb_planes) = 1
      do i = 2, polyhedron%nb_points
         if (abs(n_distance(i) - n_distance(i-1)) > DISTANCE_THRESHOLD) then
            nb_planes = nb_planes + 1
            n_group(nb_planes) = i
         end if
      end do
      n_group(nb_planes+1) = polyhedron%nb_points + 1

      ! Test if the first three points belong to the same face
      ! If the face normal is approximately collinear to the flood direction,
      ! re-structure the order of the points.
      face = cg3_polyhedron_find_common_face_3_points(polyhedron, point_to_edge, n_point(1), n_point(2), n_point(3))
      if (face > 0) then
         if (norm2(cg3_cross_product(polyhedron%normal(:,face), normal)) < NORMAL_THRESHOLD) then
            ! Group all the points of the bottom face in the group 1
            n_group(2) = polyhedron%face(face)%size + 1
            n_point(1:n_group(2)-1) = polyhedron%face(face)%id

            nb_planes = 1
            do i = n_group(2), polyhedron%nb_points
               if (abs(n_distance(i) - n_distance(i-1)) > DISTANCE_THRESHOLD) then
                  nb_planes = nb_planes + 1
                  n_group(nb_planes) = i
               end if
            end do
            n_group(nb_planes+1) = polyhedron%nb_points + 1
         end if
      end if

      ! Create inverse tables
      do i = 1, nb_planes
         do j = n_group(i), n_group(i+1) - 1
            inv_n_point(n_point(j)) = j
            inv_n_group(n_point(j)) = i
         end do
      end do

      ! Allocate the momentum arrary
      allocate(momentum(3,nb_planes))

      !-------------------------
      ! 2 - Initialization of the first polygon

      ! Summary: This step initialize the polygon defined as the intersection of the first plane and the polyhedron.
      ! =======  The polygon is stored as a list. Each item of the list contains information about
      !          neighbor edges and faces. Each item corresponds to a point of the polygon. The coordinates of this
      !          point is stored in the structure. There are two kinds of points: there are points that belongs to
      !          the polyhedron -- in such a way, the corresponding index is stored in the structure -- and there are
      !          points that are the intersection of the plane and an edge of the polyhedron -- the corresponding
      !          edge index is stored in the structure. The face above the edge joining the point of the previous
      !          item to the point of the current item is stored in the structure. The set of edges that start from
      !          the point of the current item and that go up is stored in the structure.

      !          To build the first polygon, there are 3 configurations to consider:
      !
      !          1 - The first group contains only one point. The polygon is a single point. No faces are associated.
      !          2 - The first group contains two points. The polygon is a single edge. Two faces are stored.
      !          3 - The first group contains >2 points. The polygon is a face of the polyhedron.
      !
      !          The first polygon is stored in the top_polygon variable.

      nb_points_top = n_group(2) - n_group(1)

      ! Allocate memory for top_polygon and top_polygon
      ! The maximal size for top_polygon is equal to the number of faces (worst case: quadrilateral section of a tetrahedron).
      allocate(top_polygon%point(polyhedron%nb_faces))
      allocate(bottom_polygon%point(polyhedron%nb_faces))

      if (nb_points_top == 1) then
         ! The top polygon is a point, there is no face associated
         call top_polygon%add_point(polyhedron%point(:,n_point(1)), n_point(1), 0, 0)
      else if (nb_points_top == 2) then
         ! The top polygon is a edge

         ! Find the common edge
         edge = cg3_polyhedron_find_common_edge(polyhedron, point_to_edge, n_point(1), n_point(2))

         call top_polygon%add_point(polyhedron%point(:,n_point(1)), n_point(1), 0, polyhedron%edge_to_face(1,edge))
         call top_polygon%add_point(polyhedron%point(:,n_point(2)), n_point(2), 0, polyhedron%edge_to_face(2,edge))
      else
         ! The top polygon is a face
         ! Find the common face
         face = cg3_polyhedron_find_common_face_3_points(polyhedron, point_to_edge, n_point(1), n_point(2), n_point(3))

         ! Copy the points of the face to ensure the right order
         do i = 1, nb_points_top
            n_point(i) = polyhedron%face(face)%id(i)
         end do

         ! Create the polygon with the point of the face in reverse order to get a polygon with its points
         ! oriented in counter-clockwise order when watched from above.
         do i = nb_points_top, 1, -1
            point = polyhedron%face(face)%id(i)

            if (i > 1) then
               edge = cg3_polyhedron_find_common_edge(polyhedron, point_to_edge, n_point(i), n_point(i-1))
            else
               edge = cg3_polyhedron_find_common_edge(polyhedron, point_to_edge, n_point(1), n_point(nb_points_top))
            end if

            if (polyhedron%edge_to_face(1,edge) == face) then
               call top_polygon%add_point(polyhedron%point(:,point), point, 0, polyhedron%edge_to_face(2,edge))
            else
               call top_polygon%add_point(polyhedron%point(:,point), point, 0, polyhedron%edge_to_face(1,edge))
            end if
         end do
      end if

      ! Create and order the set of edges that go up passing through a point
      do i = 1, top_polygon%nb_points
         call cg3_create_edge_set(polyhedron, point_to_edge, top_polygon%point(i), inv_n_group, 1)
      end do

      !-------------------------
      ! 3 - Bracket the solution

      ! Summary: Find the plane such that the volume below it is equal to the prescribed volume.
      ! =======  This step consists in bracketting this plane between two planes of the sorted points.
      !          The shape contained between two planes is known as a prismatoid.
      !
      !          The two polygons bracketting the plane are stored in top_polygon and top_polygon.
      !
      !          Before constructing the top polygon, we check if the volume of the prismatoid above
      !          the top polygon + the volume below the top polygon is greater than the prescribed
      !          volume. Note that we do not need to construct the top polygon to compute the volume!
      !          If the volume exceed the prescribed volume, the final plane is found using the fact that
      !          the evolution of the volume in a prismatoid is a cubic function of the n-distance. A Newton-
      !          Raphson algorithm is used to find the correct root. If the volume is smaller that the
      !          prescribed volume, the top polygon is constructed and we start again this step from
      !          the top polygon.

      momentum = 0d0
      total_volume = 0d0
      height = 0d0
      A = 0d0; B = 0d0; C = 0d0
      cap_is_cut = .false.
      cap_is_bottom = .false.

      ! Loop over all the planes
      do k = 2, nb_planes
         ! Compute the height of the prismatoid
         height = n_distance(n_group(k)) - n_distance(n_group(k-1))

         ! Compute the volume of the prismatoid
         associate(top_point => top_polygon%point)
            A = 0d0; B = 0d0; C = 0d0
            polygon_area = 0d0
            prism_momentum = 0d0
            area = 0d0

            prevp = top_polygon%nb_points
            do curp = 1, top_polygon%nb_points
               ! Tetrahedron contributions
               do i = 1, size(top_point(curp)%edge_set) - 1
                  volume_tetra = dot_product( &
                     &              cg3_cross_product(top_point(curp)%tangent(:,i), top_point(curp)%tangent(:,i+1)), &
                     &              normal &
                     &           ) &
                     &         / (dot_product(top_point(curp)%tangent(:,i),   normal) &
                     &         *  dot_product(top_point(curp)%tangent(:,i+1), normal))

                  A = A + volume_tetra

                  momentum(:,k) = momentum(:,k) + volume_tetra*height**3/6d0*(top_point(curp)%point                   &
                     &          + height/4d0*(normal                                                                  &
                     &          + top_point(curp)%tangent(:,i)  /dot_product(normal, top_point(curp)%tangent(:,i))    &
                     &          + top_point(curp)%tangent(:,i+1)/dot_product(normal, top_point(curp)%tangent(:,i+1))))
               end do

               ! Compute the area the top polygon
               if (curp > 2) then
                  area = norm2(cg3_cross_product(top_point(prevp)%point - top_point(1)%point,&
                     &                           top_point(curp)%point  - top_point(1)%point))
                  polygon_area = polygon_area + area
                  ! Compute the centroid of the top polygon
                  prism_momentum = prism_momentum + area*(top_point(prevp)%point + top_point(curp)%point - 2d0*top_point(1)%point)
               end if

               prevp = curp
            end do

            C = polygon_area/2d0

            ! Add the contribution of the prism momentum
            if (top_polygon%nb_points > 2) then
               momentum(:,k) = momentum(:,k) + C*height*((top_point(1)%point + prism_momentum/(6d0*C)) + height*normal/2d0)
            end if

            ! Wedge contribution
            if (top_polygon%nb_points > 1) then
               prevp = top_polygon%nb_points
               do curp = 1, top_polygon%nb_points
                  Bf = norm2(top_point(curp)%point - top_point(prevp)%point)
                  S = norm2(cg3_cross_product(normal, polyhedron%normal(:,top_point(prevp)%face)))
                  volume_prism = - Bf*dot_product(polyhedron%normal(:,top_point(prevp)%face), normal)/S

                  B = B + volume_prism

                  edge_prism = cg3_cross_product(polyhedron%normal(:,top_point(prevp)%face),        &
                     &                           top_point(curp)%point - top_point(prevp)%point)/Bf

                  ! May happen if nb_point_top = 2 (I don't know how to fix this without this test)
                  if (dot_product(edge_prism, normal) < 0d0) edge_prism = -edge_prism

                  momentum(:,k) = momentum(:,k) + volume_prism*height**2/2d0*(height/3d0*(normal + edge_prism/S) &
                     &          + (top_point(curp)%point + top_point(prevp)%point)/2d0)

                  nb_edges = size(top_point(prevp)%edge_set)

                  volume_tetra = dot_product(cg3_cross_product(top_point(prevp)%tangent(:,nb_edges), edge_prism), normal) &
                     &         / (dot_product(top_point(prevp)%tangent(:,nb_edges), normal)*dot_product(edge_prism, normal))

                  A = A + volume_tetra

                  momentum(:,k) = momentum(:,k) + volume_tetra*height**3/6d0*(top_point(prevp)%point                             &
                     &          + height/4d0*(normal                                                                             &
                     &          + top_point(prevp)%tangent(:,nb_edges)/dot_product(normal, top_point(prevp)%tangent(:,nb_edges)) &
                     &          + edge_prism/dot_product(normal, edge_prism)))

                  volume_tetra = dot_product(cg3_cross_product(edge_prism, top_point(curp)%tangent(:,1)), normal)  &
                     &         / (dot_product(edge_prism, normal)*dot_product(top_point(curp)%tangent(:,1), normal))

                  A = A + volume_tetra

                  momentum(:,k) = momentum(:,k) + volume_tetra*height**3/6d0*(top_point(curp)%point              &
                     &          + height/4d0*(normal                                                             &
                     &          + top_point(curp)%tangent(:,1)/dot_product(normal, top_point(curp)%tangent(:,1)) &
                     &          + edge_prism/dot_product(normal, edge_prism)))

                  prevp = curp
               end do
            else
               ! Remaining tetrahedron
               nb_edges = size(top_point(1)%edge_set)
               volume_tetra = dot_product(cg3_cross_product(top_point(1)%tangent(:,nb_edges),   &
                  &                                         top_point(1)%tangent(:,1)), normal) &
                  &         / (dot_product(top_point(1)%tangent(:,nb_edges), normal)            &
                  &         * dot_product(top_point(1)%tangent(:,1), normal))
               A = A + volume_tetra
               momentum(:,k) = momentum(:,k) + volume_tetra*height**3/6d0*(top_point(1)%point                         &
                  &          + height/4d0*(normal                                                                     &
                  &          + top_point(1)%tangent(:,1)/dot_product(normal, top_point(1)%tangent(:,1))               &
                  &          + top_point(1)%tangent(:,nb_edges)/dot_product(normal, top_point(1)%tangent(:,nb_edges))))
            end if
         end associate

         volume_prismatoid = height*(C + height*(B/2d0 + height*A/6d0))

         ! Check if we have bracketed the volume between two planes
         if (total_volume + volume_prismatoid >= volume) then
            ! Newton-Raphson algorithm to find the final height.
            ! The algorithm starts from the middle of the prismatoid
            height_old = 0.5d0*height
            do i = 1, ITER_NEWTON
               ! Compute the derivative
               height = (C + height_old*(B + height_old*A/2d0))

               ! Check for zero derivative
               if (abs(height) < NEWTON_EPSILON) exit

               height = height_old - (total_volume - volume + height_old*(C + height_old*(B/2d0 + height_old*A/6d0)))/height

               ! Check for convergence
               if (abs(height - height_old) < NEWTON_EPSILON) exit

               height_old = height
            end do

            if (abs(height) <= NEWTON_EPSILON) then
               ! The bottom polygon becomes the cap polygon
               call cg3_move_alloc_chained_polygon(top_polygon, cap_polygon)
               ! Set bottom flag
               cap_is_bottom = .true.
               ! Set flag
               cap_is_cut = .true.
               ! Exit the loop
               exit
            else if (abs(n_distance(n_group(k)) - n_distance(n_group(k-1)) - height) <= NEWTON_EPSILON) then
               ! Set flag
               cap_is_cut = .true.
               ! Compute the height of the prismatoid
               height = n_distance(n_group(k)) - n_distance(n_group(k-1))
               ! Construct the next plane
            else
               ! Exit the loop
               exit
            end if
         end if

         ! The top polygon becomes the bottom polygon!
         call cg3_move_alloc_chained_polygon(top_polygon, bottom_polygon)

         total_volume = total_volume + volume_prismatoid

         ! Compute the next polygon
         ! Loop on every point of the bottom polygon
         p_loop: do curp = 1, bottom_polygon%nb_points
            ! Check if the current point is a point of the polyhedron
            if (bottom_polygon%point(curp)%edge > 0) then
               ! The current point is not a polyhedron point => we are on an edge of the polyhedron
               ! Check if the end of the edge belongs to the next plane
               point = polyhedron%edge(1,bottom_polygon%point(curp)%edge)
               if (inv_n_group(point) == k) then
                  call top_polygon%add_point(polyhedron%point(:,point), point, 0, 0)
                  cycle p_loop
               end if

               point = polyhedron%edge(2,bottom_polygon%point(curp)%edge)
               if (inv_n_group(point) == k) then
                  call top_polygon%add_point(polyhedron%point(:,point), point, 0, 0)
                  cycle p_loop
               end if

               ! The end of the edge does not belong the the next plane
               ! => compute the coordinates of the plane/edge intersection
               coords = bottom_polygon%point(curp)%point                                                   &
                  &   + height/dot_product(polyhedron%tangent(:, bottom_polygon%point(curp)%edge), normal) &
                  &   * polyhedron%tangent(:, bottom_polygon%point(curp)%edge)

               ! Add the intersection to the top polygon
               call top_polygon%add_point(coords, 0, bottom_polygon%point(curp)%edge, bottom_polygon%point(curp)%face)

               ! Add edge information
               allocate(top_polygon%point(top_polygon%nb_points)%edge_set(1))
               allocate(top_polygon%point(top_polygon%nb_points)%edge_point(1))
               allocate(top_polygon%point(top_polygon%nb_points)%tangent(3,1))

               top_polygon%point(top_polygon%nb_points)%edge_set = bottom_polygon%point(curp)%edge
               top_polygon%point(top_polygon%nb_points)%edge_point = 0
               if (dot_product(polyhedron%tangent(:,bottom_polygon%point(curp)%edge), normal) > 0) then
                  top_polygon%point(top_polygon%nb_points)%tangent(:,1) = polyhedron%tangent(:,bottom_polygon%point(curp)%edge)
               else
                  top_polygon%point(top_polygon%nb_points)%tangent(:,1) = -polyhedron%tangent(:,bottom_polygon%point(curp)%edge)
               end if
            else
               ! The current point is a polyhedron point
               point = bottom_polygon%point(curp)%point_id

               ! Loop on every edges of the current point that go up
               nb_edges = size(bottom_polygon%point(curp)%edge_set)
               j_loop: do j = 1, nb_edges
                  ! Check if the extreme point of the edge belongs to the next plane
                  if (inv_n_group(bottom_polygon%point(curp)%edge_point(j)) == k) then
                     call top_polygon%add_point(polyhedron%point(:,bottom_polygon%point(curp)%edge_point(j)), &
                        &                       bottom_polygon%point(curp)%edge_point(j), 0, 0)
                     cycle j_loop
                  end if

                  ! The extreme point of the edge does not belong to the next plane
                  ! => compute the coordinates of the plane/edge intersection
                  coords = bottom_polygon%point(curp)%point                                    &
                     &   + height/dot_product(bottom_polygon%point(curp)%tangent(:,j), normal) &
                     &   * bottom_polygon%point(curp)%tangent(:,j)

                  ! Add the intersection to the top polygon
                  call top_polygon%add_point(coords, 0, bottom_polygon%point(curp)%edge_set(j), 0)

                  ! Add edge information
                  allocate(top_polygon%point(top_polygon%nb_points)%edge_set(1))
                  allocate(top_polygon%point(top_polygon%nb_points)%edge_point(1))
                  allocate(top_polygon%point(top_polygon%nb_points)%tangent(3,1))

                  top_polygon%point(top_polygon%nb_points)%edge_set = bottom_polygon%point(curp)%edge_set(j)
                  top_polygon%point(top_polygon%nb_points)%edge_point = 0
                  if (dot_product(polyhedron%tangent(:,bottom_polygon%point(curp)%edge_set(j)), normal) > 0) then
                     top_polygon%point(top_polygon%nb_points)%tangent(:,1) = &
                        &  polyhedron%tangent(:,bottom_polygon%point(curp)%edge_set(j))
                  else
                     top_polygon%point(top_polygon%nb_points)%tangent(:,1) = &
                        &  -polyhedron%tangent(:,bottom_polygon%point(curp)%edge_set(j))
                  end if
               end do j_loop
            end if
         end do p_loop

         associate(top_point => top_polygon%point)
            ! Order the edge that go up passing through the points of the next planes
            do i = 1, top_polygon%nb_points
               if (top_point(i)%edge == 0) call cg3_create_edge_set(polyhedron, point_to_edge, top_point(i), inv_n_group, k)
            end do

            ! Find the faces between the points
            prevp = top_polygon%nb_points
            do curp = 1, top_polygon%nb_points
               ! Book-keeping algorithm… I am sure that we can optimize this part that travels connectivities to find a face.
               if (top_point(prevp)%edge > 0) then
                  if (top_point(curp)%edge > 0) then
                     top_point(prevp)%face = cg3_polyhedron_find_common_face_2_edges(polyhedron, top_point(prevp)%edge, &
                        &                                                            top_point(curp)%edge               )
                  else
                     top_point(prevp)%face = cg3_polyhedron_find_common_face_edge_point(polyhedron, point_to_edge, &
                        &                                                               top_point(prevp)%edge,     &
                        &                                                               top_point(curp)%point_id   )
                  end if
               else if (top_point(curp)%edge > 0) then
                  top_point(prevp)%face = cg3_polyhedron_find_common_face_edge_point(polyhedron, point_to_edge, &
                     &                                                               top_point(curp)%edge,      &
                     &                                                               top_point(prevp)%point_id  )
               else ! top_point(prevp)%point_id /= top_point(curp)%point_id
                  top_point(prevp)%face = cg3_polyhedron_find_common_face_2_points_up(polyhedron, point_to_edge,              &
                     &                                                                top_point(prevp)%point_id,              &
                     &                                                                top_point(curp)%point_id, inv_n_group, k)
               end if

               prevp = curp
            end do
         end associate

         if (cap_is_cut) then
            ! The bottom polygon becomes the cap polygon
            call cg3_move_alloc_chained_polygon(top_polygon, cap_polygon)
            exit
         end if
      end do

      ! Re-compute the momentum of the last slice
      if (.not. cap_is_cut) then
         ! Compute the volume of the prismatoid
         associate(top_point => top_polygon%point)
            A = 0d0; B = 0d0; C = 0d0
            polygon_area = 0d0
            prism_momentum = 0d0
            area = 0d0

            ! Reset momentum
            momentum(:,k) = 0d0

            prevp = top_polygon%nb_points
            do curp = 1, top_polygon%nb_points
               ! Tetrahedron contributions
               do i = 1, size(top_point(curp)%edge_set) - 1
                  volume_tetra = dot_product( &
                     &              cg3_cross_product(top_point(curp)%tangent(:,i), top_point(curp)%tangent(:,i+1)), &
                     &              normal &
                     &           ) &
                     &         / (dot_product(top_point(curp)%tangent(:,i),   normal) &
                     &         *  dot_product(top_point(curp)%tangent(:,i+1), normal))

                  A = A + volume_tetra

                  momentum(:,k) = momentum(:,k) + volume_tetra*height**3/6d0*(top_point(curp)%point                   &
                     &          + height/4d0*(normal                                                                  &
                     &          + top_point(curp)%tangent(:,i)  /dot_product(normal, top_point(curp)%tangent(:,i))    &
                     &          + top_point(curp)%tangent(:,i+1)/dot_product(normal, top_point(curp)%tangent(:,i+1))))
               end do

               ! Compute the area the top polygon
               if (curp > 2) then
                  area = norm2(cg3_cross_product(top_point(prevp)%point - top_point(1)%point,&
                     &                           top_point(curp)%point  - top_point(1)%point))
                  polygon_area = polygon_area + area
                  ! Compute the centroid of the top polygon
                  prism_momentum = prism_momentum + area*(top_point(prevp)%point + top_point(curp)%point - 2d0*top_point(1)%point)
               end if

               prevp = curp
            end do

            C = polygon_area/2d0

            ! Add the contribution of the prism momentum
            if (top_polygon%nb_points > 2) then
               momentum(:,k) = momentum(:,k) + C*height*((top_point(1)%point + prism_momentum/(6d0*C)) + height*normal/2d0)
            end if

            ! Wedge contribution
            if (top_polygon%nb_points > 1) then
               prevp = top_polygon%nb_points
               do curp = 1, top_polygon%nb_points
                  Bf = norm2(top_point(curp)%point - top_point(prevp)%point)
                  S = norm2(cg3_cross_product(normal, polyhedron%normal(:,top_point(prevp)%face)))
                  volume_prism = - Bf*dot_product(polyhedron%normal(:,top_point(prevp)%face), normal)/S

                  B = B + volume_prism

                  edge_prism = cg3_cross_product(polyhedron%normal(:,top_point(prevp)%face),        &
                     &                           top_point(curp)%point - top_point(prevp)%point)/Bf

                  ! May happen if nb_point_top = 2 (I don't know how to fix this without this test)
                  if (dot_product(edge_prism, normal) < 0d0) edge_prism = -edge_prism

                  momentum(:,k) = momentum(:,k) + volume_prism*height**2/2d0*(height/3d0*(normal + edge_prism/S) &
                     &          + (top_point(curp)%point + top_point(prevp)%point)/2d0)

                  nb_edges = size(top_point(prevp)%edge_set)

                  volume_tetra = dot_product(cg3_cross_product(top_point(prevp)%tangent(:,nb_edges), edge_prism), normal) &
                     &         / (dot_product(top_point(prevp)%tangent(:,nb_edges), normal)*dot_product(edge_prism, normal))

                  A = A + volume_tetra

                  momentum(:,k) = momentum(:,k) + volume_tetra*height**3/6d0*(top_point(prevp)%point                             &
                     &          + height/4d0*(normal                                                                             &
                     &          + top_point(prevp)%tangent(:,nb_edges)/dot_product(normal, top_point(prevp)%tangent(:,nb_edges)) &
                     &          + edge_prism/dot_product(normal, edge_prism)))

                  volume_tetra = dot_product(cg3_cross_product(edge_prism, top_point(curp)%tangent(:,1)), normal)  &
                     &         / (dot_product(edge_prism, normal)*dot_product(top_point(curp)%tangent(:,1), normal))

                  A = A + volume_tetra

                  momentum(:,k) = momentum(:,k) + volume_tetra*height**3/6d0*(top_point(curp)%point              &
                     &          + height/4d0*(normal                                                             &
                     &          + top_point(curp)%tangent(:,1)/dot_product(normal, top_point(curp)%tangent(:,1)) &
                     &          + edge_prism/dot_product(normal, edge_prism)))

                  prevp = curp
               end do
            else
               ! Remaining tetrahedron
               nb_edges = size(top_point(1)%edge_set)
               volume_tetra = dot_product(cg3_cross_product(top_point(1)%tangent(:,nb_edges),   &
                  &                                         top_point(1)%tangent(:,1)), normal) &
                  &         / (dot_product(top_point(1)%tangent(:,nb_edges), normal)            &
                  &         * dot_product(top_point(1)%tangent(:,1), normal))
               A = A + volume_tetra
               momentum(:,k) = momentum(:,k) + volume_tetra*height**3/6d0*(top_point(1)%point                         &
                  &          + height/4d0*(normal                                                                     &
                  &          + top_point(1)%tangent(:,1)/dot_product(normal, top_point(1)%tangent(:,1))               &
                  &          + top_point(1)%tangent(:,nb_edges)/dot_product(normal, top_point(1)%tangent(:,nb_edges))))
            end if
         end associate
      end if

      if (cap_is_cut .and. cap_is_bottom) k = k - 1

      centroid_full = 0d0
      do i = 2, k
         centroid_full = centroid_full + momentum(:,i)
      end do
      centroid_full = centroid_full/volume

      ! Exit here if we do not want the cap polygon
      if (.not. present(polygon)) return

      !-------------------------
      ! 4 - Compute cap polygon

      ! Summary: This step consists in constructing the cap polygon that corresponds to the intersection of the
      ! =======  final plane and the polyhedron. This steps re-uses some code used in the previous step, but it
      !          is adaptated to the current situtation.
      !
      !          This step tracks the faces that are cut by the final plane.

      if (.not. cap_is_cut) then
         allocate(cap_polygon%point(polyhedron%nb_faces))

         do curp = 1, top_polygon%nb_points
            ! Check if the current point is a point of the polyhedron
            if (top_polygon%point(curp)%edge > 0) then
               ! The current point is not a point of the polyhedron
               ! compute the coordinates of the plane/edge intersection
               coords = top_polygon%point(curp)%point                                                   &
                  &   + height/dot_product(polyhedron%tangent(:, top_polygon%point(curp)%edge), normal) &
                  &   * polyhedron%tangent(:, top_polygon%point(curp)%edge)

               ! Add the intersection to the top polygon
               call cap_polygon%add_point(coords, 0, top_polygon%point(curp)%edge, top_polygon%point(curp)%face)

               ! Add edge information
               allocate(cap_polygon%point(cap_polygon%nb_points)%edge_set(1))
               allocate(cap_polygon%point(cap_polygon%nb_points)%edge_point(1))
               allocate(cap_polygon%point(cap_polygon%nb_points)%tangent(3,1))

               cap_polygon%point(cap_polygon%nb_points)%edge_set = top_polygon%point(curp)%edge
               cap_polygon%point(cap_polygon%nb_points)%edge_point = 0
               if (dot_product(polyhedron%tangent(:,top_polygon%point(curp)%edge), normal) > 0) then
                  cap_polygon%point(cap_polygon%nb_points)%tangent(:,1) = polyhedron%tangent(:,top_polygon%point(curp)%edge)
               else
                  cap_polygon%point(cap_polygon%nb_points)%tangent(:,1) = -polyhedron%tangent(:,top_polygon%point(curp)%edge)
               end if
            else
               ! The current point is a point of the polyhedron
               point = top_polygon%point(curp)%point_id

               ! Loop on every edges of the current point that go up
               nb_edges = size(top_polygon%point(curp)%edge_set)
               do j = 1, nb_edges
                  ! The extreme point of the edge does not belong to the next plane
                  ! => compute the coordinates of the plane/edge intersection
                  coords = top_polygon%point(curp)%point                                    &
                     &   + height/dot_product(top_polygon%point(curp)%tangent(:,j), normal) &
                     &   * top_polygon%point(curp)%tangent(:,j)

                  ! Add the intersection to the top polygon
                  call cg3_chained_polygon_add_point(cap_polygon, coords, 0, top_polygon%point(curp)%edge_set(j), 0)

                  ! Add edge information
                  allocate(cap_polygon%point(cap_polygon%nb_points)%edge_set(1))
                  allocate(cap_polygon%point(cap_polygon%nb_points)%edge_point(1))
                  allocate(cap_polygon%point(cap_polygon%nb_points)%tangent(3,1))

                  cap_polygon%point(cap_polygon%nb_points)%edge_set = top_polygon%point(curp)%edge_set(j)
                  cap_polygon%point(cap_polygon%nb_points)%edge_point = 0
                  if (dot_product(polyhedron%tangent(:,top_polygon%point(curp)%edge_set(j)), normal) > 0) then
                     cap_polygon%point(cap_polygon%nb_points)%tangent(:,1) = &
                        & polyhedron%tangent(:,top_polygon%point(curp)%edge_set(j))
                  else
                     cap_polygon%point(cap_polygon%nb_points)%tangent(:,1) = &
                        & -polyhedron%tangent(:,top_polygon%point(curp)%edge_set(j))
                  end if
               end do
            end if
         end do

         ! Find the faces between the points of the cap polygon
         ! Detect also the cut faces
         ! Give a number to point without id
         associate(cap_point => cap_polygon%point)
            prevp = cap_polygon%nb_points
            do curp = 1, cap_polygon%nb_points
               cap_point(curp)%point_id = -curp

               if (cap_point(prevp)%face > 0) then
                  prevp = curp
                  cycle
               end if

               if (cap_point(prevp)%edge > 0) then
                  if (cap_point(curp)%edge > 0) then
                     cap_point(prevp)%face = cg3_polyhedron_find_common_face_2_edges(polyhedron, cap_point(prevp)%edge, &
                        &                                                            cap_point(curp)%edge)
                  else
                     cap_point(prevp)%face = cg3_polyhedron_find_common_face_edge_point(polyhedron, point_to_edge, &
                        &                                                               cap_point(prevp)%edge,     &
                        &                                                               cap_point(curp)%point_id   )
                  end if
               else if (cap_point(curp)%edge > 0) then
                  cap_point(prevp)%face = cg3_polyhedron_find_common_face_edge_point(polyhedron, point_to_edge, &
                     &                                                               cap_point(curp)%edge,      &
                     &                                                               cap_point(prevp)%point_id  )
               else ! cap_point(prevp)%point_id /= cap_point(curp)%point_id
                  cap_point(prevp)%face = cg3_polyhedron_find_common_face_2_points_up(polyhedron, point_to_edge, &
                     &                                                                cap_point(prevp)%point_id, &
                     &                                                                cap_point(curp)%point_id, inv_n_group, k-1)
               end if

               prevp = curp
            end do
         end associate

         ! The number of new points in the polyhedron is the number of points of the cap polygon.
         nb_new_points_cap = cap_polygon%nb_points
      else
         ! The number of new points in the polyhedron is the number of points of the cap polygon minus the size of group k.
         nb_new_points_cap = cap_polygon%nb_points - n_group(k+1) + n_group(k)
      end if

      ! Output cap polygon if required
      if (present(polygon)) polygon = cap_polygon
   end subroutine cg3_flood_polyhedron_centroid

   !==========================================================================================================
   ! Private routines

   !> Add a point to the chained polygon
   pure subroutine cg3_chained_polygon_add_point(polygon, point, point_id, edge, face)
      class(t_chained_polygon), intent(inout) :: polygon
      double precision, dimension(3), intent(in) :: point
      integer, intent(in) :: point_id, edge, face

      integer :: i

      ! Avoid duplicate points
      if (point_id > 0) then
         if (polygon%nb_points > 0) then
            ! If the last point corresponds to the new point, do not add new point and free memory
            if (polygon%point(polygon%nb_points)%point_id == point_id) then
               if (allocated(polygon%point(polygon%nb_points)%edge_set))   deallocate(polygon%point(polygon%nb_points)%edge_set)
               if (allocated(polygon%point(polygon%nb_points)%edge_point)) deallocate(polygon%point(polygon%nb_points)%edge_point)
               if (allocated(polygon%point(polygon%nb_points)%tangent))    deallocate(polygon%point(polygon%nb_points)%tangent)
               return
            end if

            ! If the first point corresponds to the new point, shift everything on the left
            if (polygon%point(1)%point_id == point_id) then
               do i = 1, polygon%nb_points - 1
                  polygon%point(i)%point = polygon%point(i+1)%point
                  call move_alloc(polygon%point(i+1)%edge_set, polygon%point(i)%edge_set)
                  call move_alloc(polygon%point(i+1)%edge_point, polygon%point(i)%edge_point)
                  call move_alloc(polygon%point(i+1)%tangent, polygon%point(i)%tangent)
                  polygon%point(i)%point_id = polygon%point(i+1)%point_id
                  polygon%point(i)%edge     = polygon%point(i+1)%edge
                  polygon%point(i)%face     = polygon%point(i+1)%face
               end do

               polygon%nb_points = polygon%nb_points - 1
            end if
         end if
      end if

      polygon%nb_points = polygon%nb_points + 1

      polygon%point(polygon%nb_points)%point    = point
      polygon%point(polygon%nb_points)%point_id = point_id
      polygon%point(polygon%nb_points)%edge     = edge
      polygon%point(polygon%nb_points)%face     = face
   end subroutine cg3_chained_polygon_add_point

   !> Reset a chained polygon
   pure subroutine cg3_reset_chained_polygon(polygon)
      type(t_chained_polygon), intent(inout) :: polygon

      integer :: i

      do i = 1, polygon%nb_points
         if (allocated(polygon%point(i)%edge_set))   deallocate(polygon%point(i)%edge_set)
         if (allocated(polygon%point(i)%edge_point)) deallocate(polygon%point(i)%edge_point)
         if (allocated(polygon%point(i)%tangent))    deallocate(polygon%point(i)%tangent)
      end do

      polygon%nb_points = 0
   end subroutine cg3_reset_chained_polygon

   !> move_alloc of a chained_polygon. Avoid unnecessary copies.
   pure subroutine cg3_move_alloc_chained_polygon(polygon1, polygon2)
      type(t_chained_polygon), intent(inout) :: polygon1, polygon2

      integer :: i

      if (.not. allocated(polygon2%point)) then
         allocate(polygon2%point(size(polygon1%point)))
      else if (size(polygon2%point) /= size(polygon1%point)) then
         deallocate(polygon2%point)
         allocate(polygon2%point(size(polygon1%point)))
      end if

      ! Reset polygon2
      call cg3_reset_chained_polygon(polygon2)

      do i = 1, polygon1%nb_points
         polygon2%point(i)%point = polygon1%point(i)%point
         call move_alloc(polygon1%point(i)%edge_set, polygon2%point(i)%edge_set)
         call move_alloc(polygon1%point(i)%edge_point, polygon2%point(i)%edge_point)
         call move_alloc(polygon1%point(i)%tangent, polygon2%point(i)%tangent)
         polygon2%point(i)%point_id = polygon1%point(i)%point_id
         polygon2%point(i)%edge     = polygon1%point(i)%edge
         polygon2%point(i)%face     = polygon1%point(i)%face
      end do

      polygon2%nb_points = polygon1%nb_points

      ! Reset polygon1
      call cg3_reset_chained_polygon(polygon1)
   end subroutine cg3_move_alloc_chained_polygon

   !> Create and order the set of edges passing through a point
   !!
   !! Seen from the inside of the polyhedron, edges are ordered counterclockwise around the point.
   !! Place the first edge that goes up in first position
   !!
   !! From Diot & François 2016
   pure subroutine cg3_create_edge_set(polyhedron, point_to_edge, item, inv_n_group, group_max)
      type(t_polyhedron), intent(in) :: polyhedron
      type(t_incidence_matrix), dimension(:), intent(inout) :: point_to_edge
      type(t_chained_polygon_item), intent(inout) :: item
      integer, dimension(:), intent(in) :: inv_n_group
      integer, intent(in) :: group_max

      double precision :: triple_product
      integer :: point, edge, face, i, j, first, nb_edges
      logical :: last_down, goes_up

      point = item%point_id
      edge = point_to_edge(point)%id(1)

      triple_product = dot_product(cg3_cross_product(polyhedron%normal(:,polyhedron%edge_to_face(1,edge)),  &
         &                                           polyhedron%normal(:,polyhedron%edge_to_face(2,edge))), &
         &                         polyhedron%tangent(:,edge))

      if (polyhedron%edge(2,edge) == point) triple_product = -triple_product

      if (triple_product > 0) then
         face = polyhedron%edge_to_face(2,edge)
      else
         face = polyhedron%edge_to_face(1,edge)
      end if

      ! Sort the edges
      do i = 1, point_to_edge(point)%size - 1
         do j = i + 1, point_to_edge(point)%size
            edge = point_to_edge(point)%id(j)

            if (face == polyhedron%edge_to_face(1,edge)) then
               point_to_edge(point)%id(j) = point_to_edge(point)%id(i+1)
               point_to_edge(point)%id(i+1) = edge
               face = polyhedron%edge_to_face(2,edge)
               exit
            else if (face == polyhedron%edge_to_face(2,edge)) then
               point_to_edge(point)%id(j) = point_to_edge(point)%id(i+1)
               point_to_edge(point)%id(i+1) = edge
               face = polyhedron%edge_to_face(1,edge)
               exit
            end if
         end do
      end do

      ! Find the first edge that goes up and such that the previous edge goes down
      ! and count the number of edges that goes up
      edge = point_to_edge(point)%id(point_to_edge(point)%size)
      if (polyhedron%edge(1, edge) == point) then
         last_down = inv_n_group(polyhedron%edge(2,edge)) <= group_max
      else
         last_down = inv_n_group(polyhedron%edge(1,edge)) <= group_max
      end if

      first = 0
      nb_edges = 0
      do i = 1, point_to_edge(point)%size
         edge = point_to_edge(point)%id(i)

         if (polyhedron%edge(1, edge) == point) then
            goes_up = inv_n_group(polyhedron%edge(2,edge)) > group_max
         else
            goes_up = inv_n_group(polyhedron%edge(1,edge)) > group_max
         end if

         if (goes_up) nb_edges = nb_edges + 1

         if (first > 0) cycle

         if (.not. goes_up) then
            last_down = .true.
         else if (last_down) then
            first = i
         end if
      end do

      ! Circular shift
      if (first > 1) point_to_edge(point)%id = cshift(point_to_edge(point)%id, first-1)

      ! Allocate memory
      allocate(item%edge_set(nb_edges))
      allocate(item%edge_point(nb_edges))
      allocate(item%tangent(3,nb_edges))

      ! Fill the array of edges that go up and their item%tangent
      do i = 1, nb_edges
         edge = point_to_edge(point)%id(i)

         item%edge_set(i) = edge

         ! If the first point of the edge is not the current point, reverse the item%tangent
         if (polyhedron%edge(1, edge) /= point) then
            item%edge_point(i) = polyhedron%edge(1,edge)
            item%tangent(:,i) = -polyhedron%tangent(:,edge)
         else
            item%edge_point(i) = polyhedron%edge(2,edge)
            item%tangent(:,i) = polyhedron%tangent(:,edge)
         end if
      end do
   end subroutine cg3_create_edge_set

   !> Quick sort algorithm for the distance along the n-axis
   pure recursive subroutine cg3_quicksort_n_distance(n_distance, n_point, p, r)
      double precision, dimension(:), intent(inout) :: n_distance
      integer, dimension(:), intent(inout) :: n_point
      integer, intent(in) :: p, r

      integer :: q

      if (p >= r) return

      call cg3_partition(n_distance, n_point, p, r, q)
      call cg3_quicksort_n_distance(n_distance, n_point, p, q)
      call cg3_quicksort_n_distance(n_distance, n_point, q+1, r)
   end subroutine cg3_quicksort_n_distance

   !> Partition function for the quick sort algorithm
   pure subroutine cg3_partition(n_distance, n_point, p, r, partition)
      double precision, dimension(:), intent(inout) :: n_distance
      integer, dimension(:), intent(inout) :: n_point
      integer, intent(in) :: p, r
      integer, intent(out) :: partition

      integer :: i, j, p_tmp
      double precision :: d_tmp, pivot

      pivot = n_distance(p)
      i = p - 1
      j = r + 1

      do
         do
            j = j - 1
            if (n_distance(j) <= pivot) exit
         end do

         do
            i = i + 1
            if (pivot <= n_distance(i)) exit
         end do

         if (i >= j) exit

         d_tmp = n_distance(i)
         p_tmp = n_point(i)
         n_distance(i) = n_distance(j)
         n_point(i)    = n_point(j)
         n_distance(j) = d_tmp
         n_point(j)    = p_tmp
      end do

      partition = j
   end subroutine cg3_partition

   !> Find an edge common to two points. Return 0 if no edge is found.
   pure function cg3_polyhedron_find_common_edge(polyhedron, point_to_edge, p1, p2) result(edge)
      type(t_polyhedron), intent(in) :: polyhedron
      type(t_incidence_matrix), dimension(:), intent(in) :: point_to_edge

      integer, intent(in) :: p1, p2
      integer :: edge

      integer :: i, e

      edge = 0

      do i = 1, point_to_edge(p1)%size
         e = point_to_edge(p1)%id(i)

         if (polyhedron%edge(1,e) == p1) then
            if (polyhedron%edge(2,e) == p2) then
               edge = e
               return
            end if
         else
            if (polyhedron%edge(1,e) == p2) then
               edge = e
               return
            end if
         end if
      end do
   end function cg3_polyhedron_find_common_edge

   !> Find a face common to three points. Return 0 if no face is found.
   pure function cg3_polyhedron_find_common_face_3_points(polyhedron, point_to_edge, p1, p2, p3) result(face)
      type(t_polyhedron), intent(in) :: polyhedron
      type(t_incidence_matrix), dimension(:), intent(in) :: point_to_edge
      integer, intent(in) :: p1, p2, p3
      integer :: face

      integer :: i, j, k, e1, e2, e3, f1, f2

      face = 0

      do i = 1, point_to_edge(p1)%size
         e1 = point_to_edge(p1)%id(i)

         do f1 = 1, 2
            face = polyhedron%edge_to_face(f1,e1)

            do j = 1, point_to_edge(p2)%size
               e2 = point_to_edge(p2)%id(j)

               do f2 = 1, 2
                  if (face /= polyhedron%edge_to_face(f2,e2)) cycle

                  do k = 1, point_to_edge(p3)%size
                     e3 = point_to_edge(p3)%id(k)

                     if (face == polyhedron%edge_to_face(1,e3)) return
                     if (face == polyhedron%edge_to_face(2,e3)) return
                  end do
               end do
            end do
         end do
      end do
   end function cg3_polyhedron_find_common_face_3_points

   !> Find a face common to two edges. Return 0 if no face is found.
   pure function cg3_polyhedron_find_common_face_2_edges(polyhedron, e1, e2) result(face)
      type(t_polyhedron), intent(in) :: polyhedron
      integer, intent(in) :: e1, e2
      integer :: face

      if (polyhedron%edge_to_face(1,e1) == polyhedron%edge_to_face(1,e2)) then
         face = polyhedron%edge_to_face(1,e1)
      else if (polyhedron%edge_to_face(1,e1) == polyhedron%edge_to_face(2,e2)) then
         face = polyhedron%edge_to_face(1,e1)
      else if (polyhedron%edge_to_face(2,e1) == polyhedron%edge_to_face(2,e2)) then
         face = polyhedron%edge_to_face(2,e1)
      else if (polyhedron%edge_to_face(2,e1) == polyhedron%edge_to_face(1,e2)) then
         face = polyhedron%edge_to_face(2,e1)
      else
         face = 0
      end if
   end function cg3_polyhedron_find_common_face_2_edges

   !> Find a face common to an edge and a point. Return 0 if no face is found.
   pure function cg3_polyhedron_find_common_face_edge_point(polyhedron, point_to_edge, e, p) result(face)
      type(t_polyhedron), intent(in) :: polyhedron
      type(t_incidence_matrix), dimension(:), intent(in) :: point_to_edge
      integer, intent(in) :: e, p
      integer :: face

      integer :: i, edge

      face = 0

      do i = 1, point_to_edge(p)%size
         edge = point_to_edge(p)%id(i)

         face = cg3_polyhedron_find_common_face_2_edges(polyhedron, e, edge)

         if (face > 0) return
      end do
   end function cg3_polyhedron_find_common_face_edge_point

   !> Find a face common to two points such that the face goes up. Return 0 if no face is found.
   !! Note that this algorithm requires inv_n_group to avoid real number comparison.
   pure function cg3_polyhedron_find_common_face_2_points_up(polyhedron, point_to_edge, p1, p2, inv_n_group, max_group) result(face)
      type(t_polyhedron), intent(in) :: polyhedron
      type(t_incidence_matrix), dimension(:), intent(in) :: point_to_edge
      integer, intent(in) :: p1, p2
      integer, dimension(:), intent(in) :: inv_n_group
      integer, intent(in) :: max_group
      integer :: face

      integer :: i, j, e1, e2

      face = 0

      do i = 1, point_to_edge(p1)%size
         e1 = point_to_edge(p1)%id(i)
         if (polyhedron%edge(1,e1) == p1 .and. inv_n_group(polyhedron%edge(2,e1)) <= max_group) cycle
         if (polyhedron%edge(2,e1) == p1 .and. inv_n_group(polyhedron%edge(1,e1)) <= max_group) cycle

         do j = 1, point_to_edge(p2)%size
            e2 = point_to_edge(p2)%id(j)
            if (polyhedron%edge(1,e2) == p2 .and. inv_n_group(polyhedron%edge(2,e2)) <= max_group) cycle
            if (polyhedron%edge(2,e2) == p2 .and. inv_n_group(polyhedron%edge(1,e2)) <= max_group) cycle

            face = cg3_polyhedron_find_common_face_2_edges(polyhedron, e1, e2)

            if (face > 0) return
         end do
      end do
   end function cg3_polyhedron_find_common_face_2_points_up

end module mod_cg3_flood_polyhedron

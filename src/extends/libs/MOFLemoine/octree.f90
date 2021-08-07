!This file is part of Notus 0.4.0

!Copyright Bordeaux-INP, Université de Bordeaux, CNRS
!
!Contributors:
!Antoine Lemoine, 01-07-2016, antoine.lemoine@bordeaux-inp.fr

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

!> @defgroup octree_group Sparse octree
!! @brief Sparse octree
!!
!! Sparse octree
!! =============
!!
!! Example of a valid sparse octree:
!! @verbatim
!!     ┌───────┬───┬───┐
!!     │       │   │   │     ┌┐
!!     │       ├───┼─┬─┤     ││: full nodes
!!     │       │   ├─┼─┤     └┘
!!     ├───────╆━━━┷━┷━┪     ┏┓
!!     │       ┃       ┃     ┃┃: empty nodes
!!     │       ┡━┳━┳━━━┩     ┗┛
!!     │       ┢━┛ ┃   │
!!     └───────┺━━━┹───┘
!! @endverbatim
!! Empty leaves are pruned from the tree
!!
!! @ingroup computational_geometry_3d

module mod_cg3_octree
   implicit none

   private

   !> Definition of an octree node
   !! @ingroup octree_group
   type, public :: t_octree_node
      !> Pointer to the child nodes
      type(t_octree_node), dimension(:), allocatable :: child
      !> Bounding box corner with minimal coordinates
      double precision, dimension(3) :: corner_min = -huge(1d0)
      !> Bounding box corner with maximal coordinates
      double precision, dimension(3) :: corner_max =  huge(1d0)
      !> List of objects contained in the leaves
      integer, dimension(:), allocatable :: object
   end type t_octree_node

   !> Definition of a sparse octree
   !! @ingroup octree_group
   type, public :: t_octree
      !> Pointer to the root node
      type(t_octree_node), allocatable :: root
      !> Total number of objects tracked by the octree
      integer :: nb_objects = 0
      !> Maximal depth
      integer :: max_depth = 0
      !> Maximal number of object per leaf
      integer :: max_object = huge(1)
   end type t_octree

   interface unalloc
      module procedure cg3_finalize_octree
   end interface unalloc

   public :: cg3_create_octree_from_polyhedron
   public :: cg3_octree_write_vtk_file
   public :: cg3_octree_direction_intersection
   public :: cg3_octree_ray_intersection
   public :: cg3_is_point_inside_octree
   public :: cg3_octree_get_closest_objects
   public :: unalloc

contains

   !> Finalization routine for an octree
   !!
   !! @param[inout] octree: octree to finalize
   !! @ingroup octree_group
   subroutine cg3_finalize_octree(octree)
      type(t_octree), intent(inout) :: octree

      if (allocated(octree%root)) deallocate(octree%root)
      octree%nb_objects = 0
      octree%max_depth = 0
      octree%max_object = huge(1)
   end subroutine cg3_finalize_octree

   !> Create an octree from a polyhedron
   !!
   !! @param[out] octree: octree
   !! @param[in]  polyhedron: any polyhedron
   !! @ingroup octree_group
   pure subroutine cg3_create_octree_from_polyhedron(octree, polyhedron)
      use mod_cg3_polyhedron
      type(t_octree), intent(out) :: octree
      type(t_polyhedron), intent(in) :: polyhedron

      integer :: i

      octree%nb_objects = polyhedron%nb_faces

      ! Good rule of thumb for max_depth and max_object
      octree%max_depth = int(log(dble(polyhedron%nb_faces))/log(8d0)) + 3
      octree%max_object = 15

      ! Compute bounding box
      allocate(octree%root)
      octree%root%corner_min = minval(polyhedron%point, 2)
      octree%root%corner_max = maxval(polyhedron%point, 2)

      allocate(octree%root%object(polyhedron%nb_faces))
      octree%root%object = [(i, i=1, polyhedron%nb_faces)]

      if (octree%nb_objects <= octree%max_object) return

      call distribute_faces_in_childs(octree, octree%root, polyhedron, 0)
   end subroutine cg3_create_octree_from_polyhedron

   !> Determine if a point lies inside the octree
   !!
   !! @param[in] octree: octree
   !! @param[in] point: any point
   !! @return return .true. if the point is inside the octree
   !! @ingroup octree_group
   logical pure function cg3_is_point_inside_octree(octree, point) result(is_inside)
      type(t_octree), intent(in) :: octree
      double precision, dimension(3), intent(in) :: point

      is_inside = .false.

      if (.not. allocated(octree%root)) return

      if (is_point_in_bounding_box(octree%root, point)) is_inside = .true.
   end function cg3_is_point_inside_octree

   !> Intersection of an octree with a directional ray
   !!
   !! Get a list of object that lie in intersected leaves.
   !! The ray point to positive values of the given direction.
   !!
   !! @param[in]  octree: octree
   !! @param[in]  origin: origin of the directional ray
   !! @param[in]  direction: axis of the ray
   !! @param[out] nb_objects: number of objects in intersected leaves
   !! @param[out] object_list: list of objects in intersected leaves
   !! @ingroup octree_group
   pure subroutine cg3_octree_direction_intersection(octree, origin, direction, nb_objects, object_list)
      type(t_octree), intent(in) :: octree
      double precision, dimension(3), intent(in) :: origin
      integer, intent(in) :: direction
      integer, intent(out) :: nb_objects
      integer, dimension(:), allocatable, intent(out) :: object_list

      logical, dimension(octree%nb_objects) :: is_intersected
      integer :: i, n

      nb_objects = 0

      if (.not. allocated(octree%root)) return

      is_intersected = .false.

      call octree_node_direction_intersection(octree%root, origin, direction, is_intersected)

      nb_objects = count(is_intersected)

      if (nb_objects == 0) return

      allocate(object_list(nb_objects))

      n = 1
      do i = 1, octree%nb_objects
         if (.not. is_intersected(i)) cycle

         object_list(n) = i
         n = n + 1
      end do
   end subroutine cg3_octree_direction_intersection

   !> Intersection of an octree with a ray
   !!
   !! Get a list of object that lie in intersected leaves.
   !!
   !! @param[in]  octree: octree
   !! @param[in]  ray: ray
   !! @param[out] nb_objects: number of objects in intersected leaves
   !! @param[out] object_list: list of objects in intersected leaves
   !! @ingroup octree_group
   pure subroutine cg3_octree_ray_intersection(octree, ray, nb_objects, object_list)
      use mod_ray_tracing
      type(t_octree), intent(in) :: octree
      type(t_ray), intent(in) :: ray
      integer, intent(out) :: nb_objects
      integer, dimension(:), allocatable, intent(out) :: object_list

      logical, dimension(octree%nb_objects) :: is_intersected
      integer :: i, n

      nb_objects = 0

      if (.not. allocated(octree%root)) return

      is_intersected = .false.

      call octree_node_ray_intersection(octree%root, ray, is_intersected)

      nb_objects = count(is_intersected)

      if (nb_objects == 0) return

      allocate(object_list(nb_objects))

      n = 1
      do i = 1, octree%nb_objects
         if (.not. is_intersected(i)) cycle

         object_list(n) = i
         n = n + 1
      end do
   end subroutine cg3_octree_ray_intersection

   !> Find the closest objects to a given point
   !!
   !! Get a list of object that lie in closest leaves.
   !!
   !! @param[in]  octree: octree
   !! @param[in]  point: any point
   !! @param[out] nb_objects: number of objects in closest leaves
   !! @param[out] object_list: list of objects in closest leaves
   !! @ingroup octree_group
   pure subroutine cg3_octree_get_closest_objects(octree, point, nb_objects, object_list)
      type(t_octree), intent(in) :: octree
      double precision, dimension(3), intent(in) :: point
      integer, intent(out) :: nb_objects
      integer, dimension(:), allocatable, intent(out) :: object_list

      logical, dimension(octree%nb_objects) :: is_close_enough
      double precision :: search_distance
      integer :: i, n

      nb_objects = 0

      if (.not. allocated(octree%root)) return

      search_distance = huge(1d0)

      call octree_node_get_minimal_search_distance(octree%root, point, search_distance)

      is_close_enough = .false.

      call octree_node_get_objects_in_search_distance(octree%root, point, search_distance, is_close_enough)

      nb_objects = count(is_close_enough)

      if (nb_objects == 0) return

      allocate(object_list(nb_objects))

      n = 1
      do i = 1, octree%nb_objects
         if (.not. is_close_enough(i)) cycle

         object_list(n) = i
         n = n + 1
      end do
   end subroutine cg3_octree_get_closest_objects

   !> Write the sparse octree to a VTK file
   !!
   !! @param[in]  octree: octree
   !! @param[in]  filename: VTK file name (with extension)
   !! @ingroup octree_group
   subroutine cg3_octree_write_vtk_file(octree, filename)
      type(t_octree), intent(in) :: octree
      character(len=*), intent(in) :: filename

      integer :: tree_unit, nb_points, nb_faces, nb_face_points, offset

      open(newunit=tree_unit, file=trim(filename), status="unknown")

      write(tree_unit,'("# vtk DataFile Version 3.0")')
      write(tree_unit,'(a)') trim(filename)
      write(tree_unit,'("ASCII")')
      write(tree_unit,'("DATASET POLYDATA")')
      write(tree_unit,*)

      nb_points = 0
      nb_faces = 0
      call octree_node_count_points_and_faces(octree%root, nb_points, nb_faces)

      write(tree_unit,'("POINTS ",i8," float")') nb_points

      ! Write the list of points
      call octree_write_points(octree%root, tree_unit)

      ! Count the number of integer for the faces
      nb_face_points = nb_faces*5

      write(tree_unit,'("POLYGONS ",i8," ",i8)') nb_faces, nb_face_points

      offset = 0
      call octree_write_faces(octree%root, tree_unit, offset)

      close(tree_unit)
   end subroutine cg3_octree_write_vtk_file

   !----------------------------
   ! Private routines

   pure recursive subroutine distribute_faces_in_childs(octree, node, polyhedron, depth)
      use mod_cg3_polyhedron
      type(t_octree), intent(inout) :: octree
      type(t_octree_node), intent(inout) :: node
      type(t_polyhedron), intent(in) :: polyhedron
      integer, intent(in) :: depth

      type integer_array
         integer, dimension(:), allocatable :: object
      end type integer_array

      double precision, dimension(3) :: middle
      double precision, dimension(3,3) :: grid
      double precision, dimension(3,8) :: corner_min, corner_max
      type(integer_array), dimension(8) :: childs
      logical, dimension(:), allocatable :: is_inside
      integer, dimension(8) :: nb_objects
      integer :: nb_childs, face, i, j, k, n, p

      ! Compute the bounding box of the current child node
      middle = 0.5d0*(node%corner_min + node%corner_max)

      grid(1,:) = node%corner_min
      grid(2,:) = middle
      grid(3,:) = node%corner_max

      do k = 1, 2
         do j = 1, 2
            do i = 1, 2
               corner_min(:,i+2*j+4*k-6) = [grid(i  ,1), grid(j  ,2), grid(k  ,3)]
               corner_max(:,i+2*j+4*k-6) = [grid(i+1,1), grid(j+1,2), grid(k+1,3)]
            end do
         end do
      end do

      allocate(is_inside(size(node%object)))

      ! Distribute faces in childs
      do p = 1, 8
         ! Look for faces inside the new node
         is_inside = .false.
         do face = 1, size(node%object)
            if (face_in_bounding_box(polyhedron, node%object(face), corner_min(:,p), corner_max(:,p))) is_inside(face) = .true.
         end do

         ! Count the number of faces inside the new node
         nb_objects(p) = count(is_inside)

         if (nb_objects(p) == 0) cycle

         ! Add the faces in the new node
         allocate(childs(p)%object(nb_objects(p)))

         ! Select object that are inside the new node
         ! Note: PACK is slower here because it has to count the number of elements with MASK to .true.
         n = 1
         do face = 1, size(node%object)
            if (.not. is_inside(face)) cycle

            childs(p)%object(n) = node%object(face)
            n = n + 1
         end do
      end do

      ! Count and allocate the childs
      nb_childs = count(nb_objects /= 0) ! Should be at least equal to 1
      allocate(node%child(nb_childs))

      ! Fill the childs
      n = 1
      do p = 1, 8
         if (nb_objects(p) == 0) cycle

         node%child(n)%corner_min = corner_min(:,p)
         node%child(n)%corner_max = corner_max(:,p)
         call move_alloc(childs(p)%object, node%child(n)%object)
         n = n + 1
      end do

      deallocate(node%object)
      deallocate(is_inside)

      if (depth >= octree%max_depth) return

      ! Subdivide if required
      do p = 1, size(node%child)
         if (size(node%child(p)%object) <= octree%max_object) cycle

         call distribute_faces_in_childs(octree, node%child(p), polyhedron, depth + 1)
      end do
   end subroutine distribute_faces_in_childs

   logical pure function face_in_bounding_box(polyhedron, face, corner_min, corner_max) result(is_inside)
      use mod_cg3_polyhedron
      type(t_polyhedron), intent(in) :: polyhedron
      integer, intent(in) :: face
      double precision, dimension(3), intent(in) :: corner_min, corner_max

      double precision, dimension(3) :: inf, sup
      integer :: i

      is_inside = .false.

      inf = polyhedron%point(:,polyhedron%face(face)%id(1))
      sup = polyhedron%point(:,polyhedron%face(face)%id(1))
      do i = 2, polyhedron%face(face)%size
         if (polyhedron%point(1,polyhedron%face(face)%id(i)) < inf(1)) inf(1) = polyhedron%point(1,polyhedron%face(face)%id(i))
         if (polyhedron%point(2,polyhedron%face(face)%id(i)) < inf(2)) inf(2) = polyhedron%point(2,polyhedron%face(face)%id(i))
         if (polyhedron%point(3,polyhedron%face(face)%id(i)) < inf(3)) inf(3) = polyhedron%point(3,polyhedron%face(face)%id(i))
         if (polyhedron%point(1,polyhedron%face(face)%id(i)) > sup(1)) sup(1) = polyhedron%point(1,polyhedron%face(face)%id(i))
         if (polyhedron%point(2,polyhedron%face(face)%id(i)) > sup(2)) sup(2) = polyhedron%point(2,polyhedron%face(face)%id(i))
         if (polyhedron%point(3,polyhedron%face(face)%id(i)) > sup(3)) sup(3) = polyhedron%point(3,polyhedron%face(face)%id(i))
      end do

      if (sup(1) < corner_min(1)) return
      if (sup(2) < corner_min(2)) return
      if (sup(3) < corner_min(3)) return
      if (inf(1) > corner_max(1)) return
      if (inf(2) > corner_max(2)) return
      if (inf(3) > corner_max(3)) return

      is_inside = .true.
   end function face_in_bounding_box

   logical pure function is_point_in_bounding_box(node, point) result(is_inside)
      type(t_octree_node), intent(in) :: node
      double precision, dimension(3), intent(in) :: point

      is_inside = .false.

      if (point(1) < node%corner_min(1)) return
      if (point(2) < node%corner_min(2)) return
      if (point(3) < node%corner_min(3)) return
      if (point(1) > node%corner_max(1)) return
      if (point(2) > node%corner_max(2)) return
      if (point(3) > node%corner_max(3)) return

      is_inside = .true.
   end function is_point_in_bounding_box

   pure recursive subroutine octree_node_direction_intersection(node, origin, direction, is_intersected)
      type(t_octree_node), intent(in) :: node
      double precision, dimension(3), intent(in) :: origin
      integer, intent(in) :: direction
      logical, dimension(:), intent(inout) :: is_intersected

      integer :: d2, d3, i

      select case (direction)
      case(1)
         d2 = 2; d3 = 3
      case(2)
         d2 = 3; d3 = 1
      case default ! 3
         d2 = 1; d3 = 2
      end select

      ! Check the current node
      if (origin(direction) > node%corner_max(direction)) return
      if (origin(d2) < node%corner_min(d2)) return
      if (origin(d3) < node%corner_min(d3)) return
      if (origin(d2) > node%corner_max(d2)) return
      if (origin(d3) > node%corner_max(d3)) return

      ! If the node has childs, explore the childs. Otherwise, add all the objects of the node to the list.
      if (allocated(node%child)) then
         do i = 1, size(node%child)
            call octree_node_direction_intersection(node%child(i), origin, direction, is_intersected)
         end do
      else
         do i = 1, size(node%object)
            is_intersected(node%object(i)) = .true.
         end do
      end if
   end subroutine octree_node_direction_intersection

   pure recursive subroutine octree_node_ray_intersection(node, ray, is_intersected)
      use mod_ray_tracing
      type(t_octree_node), intent(in) :: node
      type(t_ray), intent(in) :: ray
      logical, dimension(:), intent(inout) :: is_intersected

      double precision, dimension(3,2) :: bounding_box
      integer :: i

      ! Generate the bounding-box.
      bounding_box(:,1) = node%corner_min
      bounding_box(:,2) = node%corner_min

      ! Intersect the bounding-box with the ray.
      if (.not. rt_does_ray_hit_aligned_axis_bounding_box(ray, bounding_box)) return

      ! If the node has childs, explore the childs. Otherwise, add all the objects of the node to the list.
      if (allocated(node%child)) then
         do i = 1, size(node%child)
            call octree_node_ray_intersection(node%child(i), ray, is_intersected)
         end do
      else
         do i = 1, size(node%object)
            is_intersected(node%object(i)) = .true.
         end do
      end if
   end subroutine octree_node_ray_intersection

   pure subroutine octree_node_get_maximal_distance(node, point, distance)
      type(t_octree_node), intent(in) :: node
      double precision, dimension(3), intent(in) :: point
      double precision, intent(out) :: distance

      double precision :: current_distance

      distance = norm2(point - node%corner_min)

      current_distance = norm2(point - [node%corner_min(1), node%corner_min(2), node%corner_max(3)])
      if (current_distance > distance) distance = current_distance

      current_distance = norm2(point - [node%corner_min(1), node%corner_max(2), node%corner_min(3)])
      if (current_distance > distance) distance = current_distance

      current_distance = norm2(point - [node%corner_min(1), node%corner_max(2), node%corner_max(3)])
      if (current_distance > distance) distance = current_distance

      current_distance = norm2(point - [node%corner_max(1), node%corner_min(2), node%corner_min(3)])
      if (current_distance > distance) distance = current_distance

      current_distance = norm2(point - [node%corner_max(1), node%corner_min(2), node%corner_max(3)])
      if (current_distance > distance) distance = current_distance

      current_distance = norm2(point - [node%corner_max(1), node%corner_max(2), node%corner_min(3)])
      if (current_distance > distance) distance = current_distance

      current_distance = norm2(point - [node%corner_max(1), node%corner_max(2), node%corner_max(3)])
      if (current_distance > distance) distance = current_distance
   end subroutine octree_node_get_maximal_distance

   pure recursive subroutine octree_node_get_minimal_search_distance(node, point, distance)
      type(t_octree_node), intent(in) :: node
      double precision, dimension(3), intent(in) :: point
      double precision, intent(inout) :: distance

      double precision :: child_distance
      integer :: node_min, i

      ! Two configurations:
      !  - the point is in the current node
      !  - the point is outise the current node
      if (is_point_in_bounding_box(node, point)) then
         ! The minimal search distance corresponds to the maximal distance from the point to the node bounding box
         call octree_node_get_maximal_distance(node, point, distance)

         ! Explore the childs to find a shorter minimal search distance
         if (allocated(node%child)) then
            do i = 1, size(node%child)
               call octree_node_get_minimal_search_distance(node%child(i), point, distance)
            end do
         end if
      else
         ! If there are no child to explore, we can exit here
         if (.not. allocated(node%child)) return

         ! Find the child with a search distance shorter than the current minimal distance
         node_min = 0
         do i = 1, size(node%child)
            ! Compute the maximal distance from the point to the bounding box of the child node i
            call octree_node_get_maximal_distance(node%child(i), point, child_distance)

            if (child_distance > distance) cycle

            distance = child_distance
            node_min = i
         end do

         ! Explore the child if it is necessary
         if (node_min > 0) call octree_node_get_minimal_search_distance(node%child(node_min), point, distance)
      end if
   end subroutine octree_node_get_minimal_search_distance

   logical pure function does_octree_node_intersect_sphere(node, center, radius) result(is_intersected)
      type(t_octree_node), intent(in) :: node
      double precision, dimension(3), intent(in) :: center
      double precision, intent(in) :: radius

      double precision, dimension(3) :: relative_position

      is_intersected = .false.

      if (is_point_in_bounding_box(node, center)) then
         is_intersected = .true.
         return
      end if

      relative_position(1) = max(node%corner_min(1) - center(1), 0d0, center(1) - node%corner_max(1))
      relative_position(2) = max(node%corner_min(2) - center(2), 0d0, center(2) - node%corner_max(2))
      relative_position(3) = max(node%corner_min(3) - center(3), 0d0, center(3) - node%corner_max(3))

      if (norm2(relative_position) > radius) return

      is_intersected = .true.
   end function does_octree_node_intersect_sphere

   recursive pure subroutine octree_node_get_objects_in_search_distance(node, point, distance, is_close_enough)
      type(t_octree_node), intent(in) :: node
      double precision, dimension(3), intent(in) :: point
      double precision, intent(in) :: distance
      logical, dimension(:), intent(inout) :: is_close_enough

      integer :: i

      if (.not. does_octree_node_intersect_sphere(node, point, distance)) return

      if (allocated(node%object)) then
         do i = 1, size(node%object)
            is_close_enough(node%object(i)) = .true.
         end do
      end if

      if (.not. allocated(node%child)) return

      do i = 1, size(node%child)
         call octree_node_get_objects_in_search_distance(node%child(i), point, distance, is_close_enough)
      end do
   end subroutine octree_node_get_objects_in_search_distance

   pure recursive subroutine octree_node_count_points_and_faces(node, nb_points, nb_faces)
      type(t_octree_node), intent(in) :: node
      integer, intent(inout) :: nb_points, nb_faces

      integer :: i

      nb_points = nb_points + 8
      nb_faces  = nb_faces  + 6

      if (allocated(node%child)) then
         do i = 1, size(node%child)
            call octree_node_count_points_and_faces(node%child(i), nb_points, nb_faces)
         end do
      end if
   end subroutine octree_node_count_points_and_faces

   recursive subroutine octree_write_points(node, tree_unit)
      type(t_octree_node), intent(in) :: node
      integer, intent(in) :: tree_unit

      integer :: i

      write(tree_unit,*) node%corner_min(1), node%corner_min(2), node%corner_min(3)
      write(tree_unit,*) node%corner_max(1), node%corner_min(2), node%corner_min(3)
      write(tree_unit,*) node%corner_max(1), node%corner_max(2), node%corner_min(3)
      write(tree_unit,*) node%corner_min(1), node%corner_max(2), node%corner_min(3)
      write(tree_unit,*) node%corner_min(1), node%corner_max(2), node%corner_max(3)
      write(tree_unit,*) node%corner_max(1), node%corner_max(2), node%corner_max(3)
      write(tree_unit,*) node%corner_max(1), node%corner_min(2), node%corner_max(3)
      write(tree_unit,*) node%corner_min(1), node%corner_min(2), node%corner_max(3)

      if (allocated(node%child)) then
         do i = 1, size(node%child)
            call octree_write_points(node%child(i), tree_unit)
         end do
      end if
   end subroutine octree_write_points

   recursive subroutine octree_write_faces(node, tree_unit, offset)
      type(t_octree_node), intent(in) :: node
      integer, intent(in) :: tree_unit
      integer, intent(inout) :: offset

      integer :: i

      write(tree_unit,'(5i8)') 4, offset + 0, offset + 1, offset + 2, offset + 3
      write(tree_unit,'(5i8)') 4, offset + 1, offset + 2, offset + 5, offset + 6
      write(tree_unit,'(5i8)') 4, offset + 4, offset + 7, offset + 6, offset + 5
      write(tree_unit,'(5i8)') 4, offset + 0, offset + 7, offset + 4, offset + 3
      write(tree_unit,'(5i8)') 4, offset + 0, offset + 1, offset + 6, offset + 7
      write(tree_unit,'(5i8)') 4, offset + 2, offset + 3, offset + 4, offset + 5

      if (allocated(node%child)) then
         do i = 1, size(node%child)
            offset = offset + 8
            call octree_write_faces(node%child(i), tree_unit, offset)
         end do
      end if
   end subroutine octree_write_faces

end module mod_cg3_octree

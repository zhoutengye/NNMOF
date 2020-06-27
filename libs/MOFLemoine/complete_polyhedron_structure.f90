!This file is part of Notus 0.4.0

!Copyright Bordeaux-INP, Université de Bordeaux, CNRS
!
!Contributors:
!Antoine Lemoine, 10-02-2016, antoine.lemoine@bordeaux-inp.fr

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

!> @defgroup complete_polyhedron_structure Complete polyhedron structure
!! @brief Complete the polyhedron structure
!! @ingroup computational_geometry_3d

module mod_cg3_complete_polyhedron_structure
   use mod_cg3_polyhedron
   implicit none

   private

   ! Structure used to generate edges
   type t_edge
      integer, dimension(2) :: p ! End-points of the edge p(1) < p(2)
      integer :: face1 = 0       ! First neighbor face
      integer :: face2 = 0       ! Second neighbor face
      integer :: k = 0           ! Index of the edge
      logical :: mask = .true.   ! Flag used in the edge sorting algorithm
   end type t_edge

   public :: cg3_complete_polyhedron_structure

contains

   !> Complete the polyhedron structure
   !!
   !! Generate the polyhedron connectivies from face and vertices connectivities.
   !!
   !! - Check if the polyhedron is a 2-manifold with no boundary (error_id = 1)
   !! - Check if every faces has non-zero area (error_id = 2)
   !! - Check if every edges has non-zero length (error_id = 3)
   !!
   !! @param[inout] polyhedron: any polyhedron with incomplete structure
   !! @param[out]   error_id:   > 0 if an error occurs. 0 otherwise.
   !! @ingroup complete_polyhedron_structure
   pure subroutine cg3_complete_polyhedron_structure(polyhedron, error_id)
      type(t_polyhedron), intent(inout) :: polyhedron
      integer, intent(out) :: error_id

      type(t_edge), dimension(:), allocatable :: edges
      integer, dimension(:), allocatable :: f2e, p2e
      integer :: nb_edges, i, j, k, tmp, face_error, edge_error

      error_id = 0

      ! Count the number of edges
      ! Note: for a given face, the number of edges is equal to the number of vertices.
      nb_edges = 0
      do i = 1, polyhedron%nb_faces
         nb_edges = nb_edges + polyhedron%face(i)%size
      end do

      ! Generate edges
      allocate(edges(nb_edges))

      k = 1
      do i = 1, polyhedron%nb_faces
         edges(k)%p(1) = polyhedron%face(i)%id(polyhedron%face(i)%size)
         edges(k)%p(2) = polyhedron%face(i)%id(1)
         edges(k)%face1 = i
         k = k + 1

         do j = 2, polyhedron%face(i)%size
            edges(k)%p(1) = polyhedron%face(i)%id(j-1)
            edges(k)%p(2) = polyhedron%face(i)%id(j)
            edges(k)%face1 = i
            k = k + 1
         end do
      end do

      ! Sort edge points
      do i = 1, nb_edges
         if (edges(i)%p(1) < edges(i)%p(2)) cycle
         tmp = edges(i)%p(1)
         edges(i)%p(1) = edges(i)%p(2)
         edges(i)%p(2) = tmp
      end do

      ! Sort edges
      call quicksort(edges, 1, nb_edges)

      ! Count the number of unique edges
      edges(1)%k = 1
      polyhedron%nb_edges = 1
      do i = 2, nb_edges
         if (.not. all(edges(i)%p == edges(i-1)%p)) then
            polyhedron%nb_edges = polyhedron%nb_edges + 1
            edges(i)%k = polyhedron%nb_edges
            cycle
         end if

         edges(i)%mask = .false.
         edges(i-1)%face2 = edges(i)%face1
         edges(i)%face2 = edges(i-1)%face1
         edges(i)%k = polyhedron%nb_edges
      end do

      ! Allocate the real edges
      allocate(polyhedron%edge(2,polyhedron%nb_edges))

      ! Fill the edges array
      j = 1
      do i = 1, nb_edges
         if (.not. edges(i)%mask) cycle
         polyhedron%edge(:,j) = edges(i)%p
         j = j + 1
      end do

      ! Allocate the edges -> faces incidence matrix
      allocate(polyhedron%edge_to_face(2,polyhedron%nb_edges))
      allocate(f2e(polyhedron%nb_faces))

      ! Fill the edges -> faces incidence matrix and compute the number of element in the faces -> edges incidence matrix
      f2e = 0
      do i = 1, nb_edges
         if (.not. edges(i)%mask) cycle

         if (edges(i)%face2 < 1) then
            error_id = 1
            return
         end if

         polyhedron%edge_to_face(1,edges(i)%k) = edges(i)%face1
         polyhedron%edge_to_face(2,edges(i)%k) = edges(i)%face2
         f2e(edges(i)%face1) = f2e(edges(i)%face1) + 1
         f2e(edges(i)%face2) = f2e(edges(i)%face2) + 1
      end do

      ! Allocate the faces -> edges incidence matrix
      allocate(polyhedron%face_to_edge(polyhedron%nb_faces))

      do i = 1, polyhedron%nb_faces
         allocate(polyhedron%face_to_edge(i)%id(f2e(i)))
         polyhedron%face_to_edge(i)%size = f2e(i)
      end do

      deallocate(edges)

      ! Compute the faces -> edges incidence matrix
      f2e = 0
      do i = 1, polyhedron%nb_edges
         f2e(polyhedron%edge_to_face(1,i)) = f2e(polyhedron%edge_to_face(1,i)) + 1
         polyhedron%face_to_edge(polyhedron%edge_to_face(1,i))%id(f2e(polyhedron%edge_to_face(1,i))) = i
         f2e(polyhedron%edge_to_face(2,i)) = f2e(polyhedron%edge_to_face(2,i)) + 1
         polyhedron%face_to_edge(polyhedron%edge_to_face(2,i))%id(f2e(polyhedron%edge_to_face(2,i))) = i
      end do

      deallocate(f2e)

      call sort_face_to_edge(polyhedron)

      ! Compute the number of element in the points -> edges incidence matrix
      allocate(p2e(polyhedron%nb_points))
      p2e = 0
      do i = 1, polyhedron%nb_edges
         p2e(polyhedron%edge(1,i)) = p2e(polyhedron%edge(1,i)) + 1
         p2e(polyhedron%edge(2,i)) = p2e(polyhedron%edge(2,i)) + 1
      end do

      ! Allocate points -> edges incidence matrix
      allocate(polyhedron%point_to_edge(polyhedron%nb_points))

      do i = 1, polyhedron%nb_points
         allocate(polyhedron%point_to_edge(i)%id(p2e(i)))
         polyhedron%point_to_edge(i)%size = p2e(i)
      end do

      ! Fill the points -> edges incidence matrix
      p2e = 0
      do i = 1, polyhedron%nb_edges
         p2e(polyhedron%edge(1,i)) = p2e(polyhedron%edge(1,i)) + 1
         polyhedron%point_to_edge(polyhedron%edge(1,i))%id(p2e(polyhedron%edge(1,i))) = i
         p2e(polyhedron%edge(2,i)) = p2e(polyhedron%edge(2,i)) + 1
         polyhedron%point_to_edge(polyhedron%edge(2,i))%id(p2e(polyhedron%edge(2,i))) = i
      end do

      deallocate(p2e)

      call cg3_polyhedron_compute_normals(polyhedron, face_error)
      if (face_error > 0) error_id = 2
      call cg3_polyhedron_compute_tangents(polyhedron, edge_error)
      if (edge_error > 0) then
         error_id = 3
         return
      end if
   end subroutine cg3_complete_polyhedron_structure

   !====================================================================================
   ! Private routines

   !> Sort the face edges in the same order than face points
   !! @verbatim
   !!  e1 e2 e3 e4
   !! +--+--+--+--+ -> p1
   !! p1 p2 p3 p4 p5
   !! @endverbatim
   pure subroutine sort_face_to_edge(polyhedron)
      type(t_polyhedron), intent(inout) :: polyhedron

      integer :: i, j, k, tmp, cur
      integer, dimension(:), allocatable :: tab
      logical :: is_reversed_order

      allocate(tab(3))

      do i = 1, polyhedron%nb_faces
         if (size(tab) /= polyhedron%face_to_edge(i)%size) then
            deallocate(tab)
            allocate(tab(polyhedron%face_to_edge(i)%size))
         end if

         tab = polyhedron%face_to_edge(i)%id

         ! Orient the first edge (we suppose that the face points are already ordered)
         is_reversed_order = .false.
         do j = 1, polyhedron%face(i)%size-1
            if (polyhedron%edge(2,tab(1)) /= polyhedron%face(i)%id(j)) cycle
            if (polyhedron%edge(1,tab(1)) /= polyhedron%face(i)%id(j+1)) exit
            is_reversed_order = .true.
            exit
         end do
         if (polyhedron%edge(2,tab(1)) == polyhedron%face(i)%id(polyhedron%face(i)%size)) then
            if (polyhedron%edge(1,tab(1)) == polyhedron%face(i)%id(1)) is_reversed_order = .true.
         end if

         ! Sort edges
         cur = 2
         do j = 1, size(tab)-1
            do k = j+1, size(tab)
               if (polyhedron%edge(cur,tab(j)) == polyhedron%edge(1,tab(k))) then
                  tmp = tab(j+1)
                  tab(j+1) = tab(k)
                  tab(k) = tmp
                  cur = 2
                  exit
               else if (polyhedron%edge(cur,tab(j)) == polyhedron%edge(2,tab(k))) then
                  tmp = tab(j+1)
                  tab(j+1) = tab(k)
                  tab(k) = tmp
                  cur = 1
                  exit
               end if
            end do
         end do

         if (is_reversed_order) then
            polyhedron%face_to_edge(i)%id = tab(size(tab):1:-1)
            cycle
         end if

         polyhedron%face_to_edge(i)%id = tab
      end do

      deallocate(tab)
   end subroutine sort_face_to_edge

   !> Quick sort algorithm for edges
   pure recursive subroutine quicksort(tab, p, r)
      integer, intent(in) :: p, r
      type(t_edge), dimension(:), intent(inout) :: tab

      integer :: q

      if (p >= r) return

      call partition(tab, p, r, q)
      call quicksort(tab, p, q)
      call quicksort(tab, q+1, r)
   end subroutine quicksort

   !> Partition function of the quick sort algorithm
   pure subroutine partition(tab, p, r, q)
      integer, intent(in) :: p, r
      type(t_edge), dimension(:), intent(inout) :: tab
      integer, intent(out) :: q

      integer :: i, j
      type(t_edge) :: temp, pivot

      pivot = tab(p)
      i = p-1
      j = r+1

      do
         do
            j = j - 1
            if (inf(tab(j), pivot)) exit
         end do

         do
            i = i + 1
            if (inf(pivot, tab(i))) exit
         end do

         if (i >= j) exit

         temp = tab(i)
         tab(i) = tab(j)
         tab(j) = temp
      end do

      q = j
   end subroutine partition

   !> Comparison function used in the quick sort algorithm
   logical pure function inf(p1, p2)
      type(t_edge), intent(in) :: p1, p2

      inf = .true.

      if (p1%p(1) > p2%p(1)) then
         inf = .false.
         return
      end if

      if (p1%p(1) /= p2%p(1)) return
      if (p1%p(2) <= p2%p(2)) return

      inf = .false.
   end function inf

end module mod_cg3_complete_polyhedron_structure

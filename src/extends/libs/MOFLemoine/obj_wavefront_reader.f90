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

!> @defgroup obj_wavefront_reader OBJ Wavefront reader
!! @brief OBJ Wavefront reader
!! @ingroup computational_geometry_3d

module mod_cg3_obj_wavefront_reader
   use mod_parser
   use mod_cg3_polyhedron
   implicit none
   private

   enum, bind(c)
      ! Vertex data
      enumerator :: kw_obj_v = 1      ! geometric vertices
      enumerator :: kw_obj_vn         ! vertex normals
      enumerator :: kw_obj_vp         ! parameter space vertices
      enumerator :: kw_obj_vt         ! texture vertices
      ! Free-form curve/surface attributes
      enumerator :: kw_obj_cstype     ! rational or non-rational forms of curve or surface type: basis matrix, Bezier, B-spline, etc.
      enumerator :: kw_obj_deg        ! degree
      enumerator :: kw_obj_bmat       ! basis matrix
      enumerator :: kw_obj_step       ! step size
      ! Elements
      enumerator :: kw_obj_p          ! point
      enumerator :: kw_obj_l          ! line
      enumerator :: kw_obj_f          ! face
      enumerator :: kw_obj_curv       ! curve
      enumerator :: kw_obj_curv2      ! 2D curve
      enumerator :: kw_obj_surf       ! surface
      ! Free-form curve/surface body statements
      enumerator :: kw_obj_parm       ! parameter values
      enumerator :: kw_obj_trim       ! outer trimming loop
      enumerator :: kw_obj_hole       ! inner trimming loop
      enumerator :: kw_obj_rat        ! rational form flag
      enumerator :: kw_obj_scrv       ! special curve
      enumerator :: kw_obj_sp         ! special point
      enumerator :: kw_obj_end        ! end statement
      ! Connectivity between free-form surfaces
      enumerator :: kw_obj_con        ! connect
      ! Grouping
      enumerator :: kw_obj_g          ! group name
      enumerator :: kw_obj_s          ! smoothing group
      enumerator :: kw_obj_mg         ! merging group
      enumerator :: kw_obj_o          ! object name
      ! Display/render attributes
      enumerator :: kw_obj_bevel      ! bevel interpolation
      enumerator :: kw_obj_c_interp   ! color interpolation
      enumerator :: kw_obj_d_interp   ! dissolve interpolation
      enumerator :: kw_obj_lod        ! level of detail
      enumerator :: kw_obj_usemtl     ! material name
      enumerator :: kw_obj_mtllib     ! material library
      enumerator :: kw_obj_shadow_obj ! shadow casting
      enumerator :: kw_obj_trace_obj  ! ray tracing
      enumerator :: kw_obj_ctech      ! curve approximation technique
      enumerator :: kw_obj_stech      ! surface approximation technique
      enumerator :: kw_obj_zzzzzzzzz  ! DO NOT REMOVE!
   end enum

   type t_integer_list
      integer, dimension(:), allocatable :: i
      integer :: size = 0
   end type t_integer_list

   type t_vertex_list_item
      type(t_vertex_list_item), pointer :: next => null()
      double precision, dimension(3) :: coord
   end type t_vertex_list_item

   type t_vertex_list
      type(t_vertex_list_item), pointer :: first => null()
      type(t_vertex_list_item), pointer :: last => null()
      integer :: size = 0
   end type t_vertex_list

   type t_face_list_item
      type(t_face_list_item), pointer :: next => null()
      type(t_integer_list) :: vertices
   end type t_face_list_item

   type t_face_list
      type(t_face_list_item), pointer :: first => null()
      type(t_face_list_item), pointer :: last => null()
      integer :: size = 0
   end type t_face_list

   integer, parameter :: INTEGER_LIST_BLOCK_SIZE = 3

   public :: read_obj_polyhedron

contains

   !> Read a OBJ Wavefront file (.obj) in a polyhedron structure
   !!
   !! This reader ignores normal and texture information.
   !! Curve and surfaces (NURBS) are also ignored.
   !!
   !! @warning Negative references for face description are not supported.
   !!
   !! @param[in]    filename: .obj file name
   !! @param[inout] polyhedron: resulting polyhedron
   !! @ingroup obj_wavefront_reader
   subroutine read_obj_polyhedron(filename, polyhedron)
      character(len=*), intent(in) :: filename
      type(t_polyhedron), intent(inout) :: polyhedron

      type(t_parser) :: parser
      type(t_token) :: tok
      type(t_vertex_list) :: vertices
      type(t_vertex_list_item), pointer :: curv
      type(t_face_list) :: faces
      type(t_face_list_item), pointer :: curf
      integer :: i

      call initialize_obj_reader(parser, filename)

      do
         call parser%get(tok)

         select case(tok%kind)
         case(tk_eof)
            exit
         case(tk_backslash)
            cycle
         case(tk_newline)
            cycle
         case(tk_keyword)
            select case(tok%keyword_id)
            case(kw_obj_v)
               call obj_reader_read_vertex(parser, vertices)
            case(kw_obj_f)
               call obj_reader_read_face(parser, faces)
            case(kw_obj_g, kw_obj_s, kw_obj_mg, kw_obj_o, kw_obj_vn, kw_obj_vp, kw_obj_vt, kw_obj_usemtl, kw_obj_mtllib)
               call obj_reader_ignore_line(parser)
            case default
               call parser%throw_error(tok, "unsupported keyword '"//parser%lex%keywords(tok%keyword_id)%name//&
                  &                         "' in Notus Wavefront OBJ reader")
            end select
         case(tk_identifier)
            call parser%throw_error(tok, "invalid keyword '"//trim(tok%string)//"' in Notus Wavefront OBJ reader")
         case default
            call parser%throw_error(tok, "invalid Wavefront OBJ file")
         end select
      end do

      call finalize_obj_reader(parser)

      polyhedron%nb_points = vertices%size
      allocate(polyhedron%point(3, vertices%size))
      curv => vertices%first
      do i = 1, vertices%size
         polyhedron%point(:,i) = curv%coord
         curv => curv%next
      end do
      call finalize_vertex_list(vertices)

      polyhedron%nb_faces = faces%size
      allocate(polyhedron%face(faces%size))
      curf => faces%first
      do i = 1, faces%size
         polyhedron%face(i)%size = curf%vertices%size
         allocate(polyhedron%face(i)%id(curf%vertices%size))
         polyhedron%face(i)%id = curf%vertices%i(:curf%vertices%size)
         curf => curf%next
      end do
      call finalize_face_list(faces)
   end subroutine read_obj_polyhedron

   subroutine initialize_integer_list(list)
      type(t_integer_list), intent(inout) :: list

      if (allocated(list%i)) deallocate(list%i)
      allocate(list%i(INTEGER_LIST_BLOCK_SIZE))
      list%size = 0
   end subroutine initialize_integer_list

   subroutine add_integer_to_list(list, i)
      type(t_integer_list), intent(inout) :: list
      integer, intent(in) :: i

      integer, dimension(:), allocatable :: tmp

      ! Re-allocate if necessary
      if (list%size == size(list%i)) then
         call move_alloc(list%i, tmp)
         allocate(list%i(list%size+INTEGER_LIST_BLOCK_SIZE))
         list%i(1:list%size) = tmp
      end if

      list%size = list%size + 1
      list%i(list%size) = i
   end subroutine add_integer_to_list

   subroutine finalize_integer_list(list)
      type(t_integer_list), intent(inout) :: list

      if (allocated(list%i)) deallocate(list%i)
   end subroutine finalize_integer_list

   subroutine add_vertex_to_list(list, vertex)
      type(t_vertex_list), intent(inout) :: list
      double precision, dimension(3), intent(in) :: vertex

      type(t_vertex_list_item), pointer :: new_item

      allocate(new_item)
      new_item%coord = vertex

      if (list%size == 0) then
         list%first => new_item
         list%last => new_item
         list%size = list%size + 1
         return
      end if

      list%last%next => new_item
      list%last => new_item
      list%size = list%size + 1
   end subroutine add_vertex_to_list

   subroutine finalize_vertex_list(list)
      type(t_vertex_list), intent(inout) :: list

      type(t_vertex_list_item), pointer :: cur, tmp

      cur => list%first
      do while(associated(cur))
         tmp => cur%next
         deallocate(cur)
         cur => tmp
      end do

      list%first => null()
      list%last => null()
      list%size = 0
   end subroutine finalize_vertex_list

   subroutine add_item_to_face_list(list, item)
      type(t_face_list), intent(inout) :: list
      type(t_face_list_item), pointer :: item

      if (list%size == 0) then
         list%first => item
         list%last => item
         list%size = list%size + 1
         return
      end if

      list%last%next => item
      list%last => item
      list%size = list%size + 1
   end subroutine add_item_to_face_list

   subroutine finalize_face_list(list)
      type(t_face_list), intent(inout) :: list

      type(t_face_list_item), pointer :: cur, tmp

      cur => list%first
      do while(associated(cur))
         tmp => cur%next
         call finalize_integer_list(cur%vertices)
         deallocate(cur)
         cur => tmp
      end do

      list%first => null()
      list%last => null()
      list%size = 0
   end subroutine finalize_face_list

   subroutine initialize_obj_reader(parser, filename)
      type(t_parser) :: parser
      character(len=*), intent(in) :: filename

      type(t_keyword_name), dimension(kw_obj_zzzzzzzzz-1) :: keywords

      ! Vertex data
      keywords(kw_obj_v         )%name = "v"
      keywords(kw_obj_vn        )%name = "vn"
      keywords(kw_obj_vp        )%name = "vp"
      keywords(kw_obj_vt        )%name = "vt"
      ! Free-form curve/surface attributes
      keywords(kw_obj_cstype    )%name = "cstype"
      keywords(kw_obj_deg       )%name = "deg"
      keywords(kw_obj_bmat      )%name = "bmat"
      keywords(kw_obj_step      )%name = "step"
      ! Elements
      keywords(kw_obj_p         )%name = "p"
      keywords(kw_obj_l         )%name = "l"
      keywords(kw_obj_f         )%name = "f"
      keywords(kw_obj_curv      )%name = "curv"
      keywords(kw_obj_curv2     )%name = "curv2"
      keywords(kw_obj_surf      )%name = "surf"
      ! Free-form curve/surface body statements
      keywords(kw_obj_parm      )%name = "parm"
      keywords(kw_obj_trim      )%name = "trim"
      keywords(kw_obj_hole      )%name = "hole"
      keywords(kw_obj_rat       )%name = "rat"
      keywords(kw_obj_scrv      )%name = "scrv"
      keywords(kw_obj_sp        )%name = "sp"
      keywords(kw_obj_end       )%name = "end"
      ! Connectivity between free-form surfaces
      keywords(kw_obj_con       )%name = "con"
      ! Grouping
      keywords(kw_obj_g         )%name = "g"
      keywords(kw_obj_s         )%name = "s"
      keywords(kw_obj_mg        )%name = "mg"
      keywords(kw_obj_o         )%name = "o"
      ! Display/render attributes
      keywords(kw_obj_bevel     )%name = "bevel"
      keywords(kw_obj_c_interp  )%name = "c_interp"
      keywords(kw_obj_d_interp  )%name = "d_interp"
      keywords(kw_obj_lod       )%name = "lod"
      keywords(kw_obj_usemtl    )%name = "usemtl"
      keywords(kw_obj_mtllib    )%name = "mtllib"
      keywords(kw_obj_shadow_obj)%name = "shadow_obj"
      keywords(kw_obj_trace_obj )%name = "trace_obj"
      keywords(kw_obj_ctech     )%name = "ctech"
      keywords(kw_obj_stech     )%name = "stech"

      call parser_initialize(parser, filename, keywords)

      ! Do not ignore new lines
      parser%lex%ignore_newline = .false.
   end subroutine initialize_obj_reader

   subroutine finalize_obj_reader(parser)
      type(t_parser) :: parser

      call parser_finalize(parser)
   end subroutine finalize_obj_reader

   subroutine obj_reader_read_vertex(parser, vertices)
      type(t_parser) :: parser
      type(t_vertex_list) :: vertices

      type(t_token) :: tok
      integer :: nb_coord
      double precision :: coord
      double precision, dimension(3) :: coords

      coord = 0d0
      nb_coord = 0

      do
         call parser%get(tok)

         select case(tok%kind)
         case(tk_real)
            coord = tok%double_value
         case(tk_integer)
            coord = tok%integer_value
         case(tk_op_minus)
            call parser%get(tok)

            if (tok%kind == tk_integer) then
               coord = -tok%integer_value
            else if (tok%kind == tk_real)  then
               coord = -tok%double_value
            end if
         case(tk_op_plus)
            cycle
         case(tk_backslash)
            call parser%expect(tk_newline)
            cycle
         case(tk_newline)
            if (nb_coord < 3) call parser%throw_error(tok, "unexepected end of line: coordinate missing")
            exit
         case(tk_eof)
            call parser%throw_error(tok, "unexpected end of file")
         case(tk_identifier)
            call parser%throw_error(tok, "invalid token '"//trim(tok%string)//"' in Notus Wavefront OBJ reader")
         case default
            call parser%throw_error(tok, "invalid token '"//trim(token_id_to_string(tok%kind))// &
               &                         "' in Notus Wavefront OBJ reader")
         end select

         nb_coord = nb_coord + 1

         if (nb_coord > 3) cycle

         coords(nb_coord) = coord
      end do

      call add_vertex_to_list(vertices, coords)
   end subroutine obj_reader_read_vertex

   subroutine obj_reader_read_face(parser, faces)
      type(t_parser) :: parser
      type(t_face_list) :: faces

      type(t_token) :: tok
      type(t_face_list_item), pointer :: face => null()
      integer :: nb_vertex

      allocate(face)
      call initialize_integer_list(face%vertices)
      nb_vertex = 0

      do
         call parser%get(tok)

         select case(tok%kind)
         case(tk_integer)
            nb_vertex = nb_vertex + 1
            call add_integer_to_list(face%vertices, tok%integer_value)

            ! Ignore texture vertex and normal vertex information
            if (parser%next_token(tk_op_divide)) then
               if (.not. parser%next_token(tk_integer)) then
                  call parser%expect(tk_op_divide, .false.)
                  call parser%expect(tk_integer, .false.)
               else if (parser%next_token(tk_op_divide)) then
                  call parser%expect(tk_integer, .false.)
               end if
            end if
         case(tk_backslash)
            call parser%expect(tk_newline)
            cycle
         case(tk_op_minus)
            call parser%throw_error(tok, "negative reference number not supported")
         case(tk_newline)
            if (nb_vertex < 3) call parser%throw_error(tok, "incomplete face")
            exit
         case(tk_eof)
            if (nb_vertex < 3) call parser%throw_error(tok, "incomplete face")
            exit
         case(tk_identifier)
            call parser%throw_error(tok, "invalid keyword '"//trim(tok%string)//"' in Notus Wavefront OBJ reader")
         case default
            call parser%throw_error(tok, "invalid token '"//trim(token_id_to_string(tok%kind))// &
               &                         "' in Notus Wavefront OBJ reader")
         end select
      end do

      call add_item_to_face_list(faces, face)
   end subroutine obj_reader_read_face

   subroutine obj_reader_ignore_line(parser)
      type(t_parser) :: parser

      character :: c
      logical :: ignore_newline

      ignore_newline = .false.

      do
         c = source_file_read_character(parser%lex%source_file(parser%lex%nb_file))

         if (parser%lex%source_file(parser%lex%nb_file)%eof) return

         if (c == achar(92)) ignore_newline = .true. ! Backslash '\'

         if (c == new_line('a')) then
            if (ignore_newline) then
               ignore_newline = .false.
               cycle
            end if
            call source_file_unread_character(parser%lex%source_file(parser%lex%nb_file))
            return
         end if
      end do
   end subroutine obj_reader_ignore_line

end module mod_cg3_obj_wavefront_reader

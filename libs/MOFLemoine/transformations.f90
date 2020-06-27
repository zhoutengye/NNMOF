!This file is part of Notus 0.4.0

!Copyright Bordeaux-INP, Université de Bordeaux, CNRS
!
!Contributors:
!Antoine Lemoine, 28-10-2015, antoine.lemoine@bordeaux-inp.fr

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

!> @defgroup transformations Transformations
!!
!! @brief Transformations using homogeneous coordinates
!!
!! This module contains routines to apply transformations on points and
!! directions in 2D and 3D. Available transformations are scaling, rotation,
!! and translation.
!!
!! # Mathematical formulation
!!
!! This implementation is based on the homogeneous coordinate system of Möbius.
!! Points and directions, or vectors, are represented by column-matrices
!! composed of 3 rows in 2D and 4 rows in 3D; transformations are represented
!! by 3×3 matrices in 2D and 4×4 matrices in 3D. Their structures are:
!! @verbatim
!!       Rotation   Translation
!!  ┌───────┴───────┐ ┌─┴─┐
!! ┌─────┬─────┬─────┬─────┐             ┌────┐                ┌────┐
!! │ Sxx │ Rxy │ Rxz │  Tx │             │ Px │                │ Dx │
!! ├─────┼─────┼─────┼─────┤             ├────┤                ├────┤
!! │ Ryx │ Syy │ Ryz │  Ty │             │ Py │                │ Dy │
!! ├─────┼─────┼─────┼─────┤     Points: ├────┤    Directions: ├────┤
!! │ Rzx │ Rzy │ Szz │  Tz │             │ Pz │                │ Dz │
!! ├─────┼─────┼─────┼─────┤             ├────┤                ├────┤
!! │  0  │  0  │  0  │  1  │             │  1 │                │  0 │
!! └─────┴─────┴─────┴─────┘             └────┘                └────┘
!! @endverbatim
!! where `R` coefficients are linked to rotations, `T` coefficients are linked
!! to translations, and `S` coefficients are linked to rotation and scaling.
!! Points and directions are differentiated by their extra coefficient: 1 for
!! a point and 0 for a direction.
!!
!! Transformations are applied to points and directions through classical
!! matrix product. Composition of transformations is also applied through
!! matrix product. The reciprocal transformation is obtained by matrix
!! inversion.
!!
!! # How to use transformations?
!!
!! First, use the module mod_cg_transformation:
!!
!! ~~~{.f90}
!! use mod_cg_transformation
!! ~~~
!! Then, declare a transformation:
!!
!! ~~~{.f90}
!! type(cg_transformation) :: transformation
!! ~~~
!!
!! Before any use of a transformation, you need to initialize a transformation.
!! To do so, you need to provide the spatial dimension to the initialization routine:
!!
!! ~~~{.f90}
!! call initialize(transformation, 3)
!! ~~~
!! This operation initializes the transformation matrices to the identity.
!! After that, you can add transformations. For instance, you can add a rotation
!! of pi/4 around the z axis:
!!
!! ~~~{.f90}
!! call cg_transformation_add_rotation_z(transformation, pi/4d0)
!! ~~~
!!
!! You can compose this transformation with a translation with the following code:
!!
!! ~~~{.f90}
!! call cg_transformation_add_translation(transformation, [1d0, 0d0, 0d0])
!! ~~~
!!
!! Now consider two variables representing a point and a direction, namely:
!!
!! ~~~{.f90}
!! double precision, dimension(3) :: point, transformed_point
!! double precision, dimension(3) :: direction, transformed_direction
!! ~~~
!!
!! To apply a transformation on these objects, simply use:
!!
!! ~~~{.f90}
!! transformed_point = cg_transform_point(transformation, point)
!! transformed_direction = cg_transform_direction(transformation, direction)
!! ~~~
!!
!! You can apply the inverse transformation to the result:
!!
!! ~~~{.f90}
!! point = cg_inverse_transform_point(transformation, transformed_point)
!! direction = cg_inverse_transform_direction(transformation, transformed_direction)
!! ~~~
!!
!! @ingroup computational_geometry

module mod_cg_transformation
   implicit none

   !> Representation of a transformation
   !! @ingroup transformations
   type cg_transformation
      !> Transformation matrix
      double precision, dimension(4,4) :: matrix = reshape([1d0, 0d0, 0d0, 0d0, &
         &                                                  0d0, 1d0, 0d0, 0d0, &
         &                                                  0d0, 0d0, 1d0, 0d0, &
         &                                                  0d0, 0d0, 0d0, 1d0], [4,4])
      !> Inverse of the transformation matrix
      double precision, dimension(4,4) :: inverse_matrix = reshape([1d0, 0d0, 0d0, 0d0, &
         &                                                          0d0, 1d0, 0d0, 0d0, &
         &                                                          0d0, 0d0, 1d0, 0d0, &
         &                                                          0d0, 0d0, 0d0, 1d0], [4,4])
      !> Dimension of the transformation matrix
      integer :: dimension = 0
   contains
      ! Initialize the transformation to the identity.
      procedure :: initialize => cg_transformation_initialize
      ! Compose two transformations.
      procedure :: compose => cg_transformation_compose
      ! Add a rotation around the X axis.
      procedure :: add_rotation_x => cg_transformation_add_rotation_x
      ! Add a rotation around the Y axis.
      procedure :: add_rotation_y => cg_transformation_add_rotation_y
      ! Add a rotation around the Z axis.
      procedure :: add_rotation_z => cg_transformation_add_rotation_z
      ! Add a rotation around a given axis.
      procedure :: add_rotation => cg_transformation_add_rotation
      ! Add a rotation around a given axis.
      procedure :: add_rotation_cos_sin => cg_transformation_add_rotation_cos_sin
      ! Add a translation.
      procedure :: add_translation => cg_transformation_add_translation
      ! Add a scaling factor in direction X.
      procedure :: add_scale_x => cg_transformation_add_scale_x
      ! Add a scaling factor in direction Y.
      procedure :: add_scale_y => cg_transformation_add_scale_y
      ! Add a scaling factor in direction Z.
      procedure :: add_scale_z => cg_transformation_add_scale_z
      ! Add a scaling factor in every directions.
      procedure :: add_scale => cg_transformation_add_scale
      ! Apply the transformation to a point.
      procedure :: transform_point => cg_transform_point
      ! Apply the transformation to a direction.
      procedure :: transform_direction => cg_transform_direction
      ! Apply the inverse transformation to a point.
      procedure :: inverse_transform_point => cg_inverse_transform_point
      ! Apply the inverse transformation to a direction.
      procedure :: inverse_transform_direction => cg_inverse_transform_direction
   end type cg_transformation

   interface initialize
      module procedure cg_transformation_initialize
   end interface initialize

contains

   !> Initialize the transformation to the identity
   !!
   !! The initialization results in the following matrix:
   !! @verbatim
   !! ┌─────┬─────┬─────┬─────┐
   !! │  1  │  0  │  0  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  1  │  0  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  0  │  1  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  0  │  0  │  1  │
   !! └─────┴─────┴─────┴─────┘
   !! @endverbatim
   !!
   !! @param[inout] transformation: transformation
   !! @param[in]    dimension: spatial dimension
   !! @ingroup transformations
   pure subroutine cg_transformation_initialize(transformation, dimension)
      class(cg_transformation), intent(inout) :: transformation
      integer, intent(in) :: dimension

      integer :: i

      transformation%matrix = 0d0

      do i = 1, dimension+1
         transformation%matrix(i,i) = 1d0
      end do

      transformation%inverse_matrix = transformation%matrix

      transformation%dimension = dimension
   end subroutine cg_transformation_initialize

   !> Compose two transformations.
   !!
   !! Add a transformation to an other. Perform a matrix-matrix product.
   !!
   !! @param[in]    transformation: transformation to add
   !! @param[inout] target: transformation that receive the transformation
   !! @ingroup transformations
   pure subroutine cg_transformation_compose(transformation, target)
      class(cg_transformation), intent(in) :: transformation
      type(cg_transformation), intent(inout) :: target

      target%matrix = matmul(transformation%matrix, target%matrix)
      target%inverse_matrix = matmul(target%inverse_matrix, transformation%inverse_matrix)
   end subroutine cg_transformation_compose

   !> Add a rotation around the x axis
   !!
   !! Compose by the following rotation matrix:
   !! @verbatim
   !! ┌─────┬─────┬─────┬─────┐
   !! │  1  │  0  │  0  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  c  │ -s  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  s  │  c  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  0  │  0  │  1  │
   !! └─────┴─────┴─────┴─────┘
   !! @endverbatim
   !! Where s and c are the sine and the cosine of the @p angle.
   !!
   !! @warning Do not apply this transformation in 2D
   !!
   !! @param[inout] transformation: transformation that receive the transformation
   !! @param[in]    angle: angle of rotation in radians
   !! @ingroup transformations
   pure subroutine cg_transformation_add_rotation_x(transformation, angle)
      class(cg_transformation), intent(inout) :: transformation
      double precision, intent(in) :: angle

      double precision, dimension(4,4) :: rotation_x

      rotation_x = 0d0

      rotation_x(1,1) = 1d0
      rotation_x(2,2) = cos(angle)
      rotation_x(2,3) = -sin(angle)
      rotation_x(3,2) = sin(angle)
      rotation_x(3,3) = cos(angle)
      rotation_x(4,4) = 1d0

      transformation%matrix = matmul(rotation_x, transformation%matrix)
      transformation%inverse_matrix = matmul(transformation%inverse_matrix, transpose(rotation_x))
   end subroutine cg_transformation_add_rotation_x

   !> Add a rotation around the y axis
   !!
   !! Compose by the following rotation matrix:
   !! @verbatim
   !! ┌─────┬─────┬─────┬─────┐
   !! │  c  │  0  │  s  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  1  │  0  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │ -s  │  0  │  c  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  0  │  0  │  1  │
   !! └─────┴─────┴─────┴─────┘
   !! @endverbatim
   !! Where s and c are the sine and the cosine of the @p angle.
   !!
   !! @warning Do not apply this transformation in 2D
   !!
   !! @param[inout] transformation: transformation that receive the transformation
   !! @param[in]    angle: angle of rotation in radians
   !! @ingroup transformations
   pure subroutine cg_transformation_add_rotation_y(transformation, angle)
      class(cg_transformation), intent(inout) :: transformation
      double precision, intent(in) :: angle

      double precision, dimension(4,4) :: rotation_y

      rotation_y = 0d0

      rotation_y(1,1) = cos(angle)
      rotation_y(1,3) = sin(angle)
      rotation_y(2,2) = 1d0
      rotation_y(3,1) = -sin(angle)
      rotation_y(3,3) = cos(angle)
      rotation_y(4,4) = 1d0

      transformation%matrix = matmul(rotation_y, transformation%matrix)
      transformation%inverse_matrix = matmul(transformation%inverse_matrix, transpose(rotation_y))
   end subroutine cg_transformation_add_rotation_y

   !> Add a rotation around the z axis
   !!
   !! Compose by the following rotation matrix:
   !! @verbatim
   !! ┌─────┬─────┬─────┬─────┐
   !! │  c  │ -s  │  0  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  s  │  c  │  0  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  0  │  1  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  0  │  0  │  1  │
   !! └─────┴─────┴─────┴─────┘
   !! @endverbatim
   !! Where s and c are the sine and the cosine of the @p angle.
   !!
   !! @param[inout] transformation: transformation that receive the transformation
   !! @param[in]    angle: angle of rotation in radians
   !! @ingroup transformations
   pure subroutine cg_transformation_add_rotation_z(transformation, angle)
      class(cg_transformation), intent(inout) :: transformation
      double precision, intent(in) :: angle

      double precision, dimension(4,4) :: rotation_z

      rotation_z = 0d0

      rotation_z(1,1) = cos(angle)
      rotation_z(1,2) = -sin(angle)
      rotation_z(2,1) = sin(angle)
      rotation_z(2,2) = cos(angle)
      rotation_z(3,3) = 1d0
      rotation_z(4,4) = 1d0

      transformation%matrix = matmul(rotation_z, transformation%matrix)
      transformation%inverse_matrix = matmul(transformation%inverse_matrix, transpose(rotation_z))
   end subroutine cg_transformation_add_rotation_z

   !> Add a rotation around a given axis
   !!
   !! Compose by a rotation matrix. Refer to the code to see the coefficients.
   !! Note that the @p axis can be a non-unit vector. It is normalized in the routine.
   !! A null axis results in a null rotation, that is the identity matrix.
   !!
   !! @param[inout] transformation: transformation that receive the transformation
   !! @param[in]    axis: rotation axis which will be normalized
   !! @param[in]    angle: angle of rotation in radians
   !! @ingroup transformations
   pure subroutine cg_transformation_add_rotation(transformation, axis, angle)
      class(cg_transformation), intent(inout) :: transformation
      double precision, dimension(3), intent(in) :: axis
      double precision, intent(in) :: angle

      double precision, dimension(4,4) :: rotation
      double precision, dimension(3) :: unit_axis
      double precision :: norm, c, s

      norm = norm2(axis)

      if (norm < epsilon(1d0)) return

      unit_axis = axis/norm

      c = cos(angle)
      s = sin(angle)

      rotation = 0d0

      rotation(1,1) = unit_axis(1)**2*(1d0 - c) + c
      rotation(2,2) = unit_axis(2)**2*(1d0 - c) + c
      rotation(3,3) = unit_axis(3)**2*(1d0 - c) + c
      rotation(4,4) = 1d0

      rotation(1,2) = unit_axis(1)*unit_axis(2)*(1d0 - c) - unit_axis(3)*s
      rotation(1,3) = unit_axis(3)*unit_axis(1)*(1d0 - c) + unit_axis(2)*s
      rotation(2,3) = unit_axis(2)*unit_axis(3)*(1d0 - c) - unit_axis(1)*s

      rotation(2,1) = unit_axis(1)*unit_axis(2)*(1d0 - c) + unit_axis(3)*s
      rotation(3,1) = unit_axis(3)*unit_axis(1)*(1d0 - c) - unit_axis(2)*s
      rotation(3,2) = unit_axis(2)*unit_axis(3)*(1d0 - c) + unit_axis(1)*s

      transformation%matrix = matmul(rotation, transformation%matrix)
      transformation%inverse_matrix = matmul(transformation%inverse_matrix, transpose(rotation))
   end subroutine cg_transformation_add_rotation

   !> Add a rotation around a given axis
   !!
   !! Compose by a rotation matrix. Refer to the code to see the coefficients.
   !! Note that the @p axis can be a non-unit vector. It is normalized in the routine.
   !! A null axis results in a null rotation, that is the identity matrix.
   !!
   !! @param[inout] transformation: transformation that receive the transformation
   !! @param[in]    axis: rotation axis which will be normalized
   !! @param[in]    c, s: cosine and sine of the angle of rotation
   !! @ingroup transformations
   pure subroutine cg_transformation_add_rotation_cos_sin(transformation, axis, c, s)
      class(cg_transformation), intent(inout) :: transformation
      double precision, dimension(3), intent(in) :: axis
      double precision, intent(in) :: c, s

      double precision, dimension(4,4) :: rotation
      double precision, dimension(3) :: unit_axis
      double precision :: norm

      norm = norm2(axis)

      if (norm < epsilon(1d0)) return

      unit_axis = axis/norm

      rotation = 0d0

      rotation(1,1) = unit_axis(1)**2*(1d0 - c) + c
      rotation(2,2) = unit_axis(2)**2*(1d0 - c) + c
      rotation(3,3) = unit_axis(3)**2*(1d0 - c) + c
      rotation(4,4) = 1d0

      rotation(1,2) = unit_axis(1)*unit_axis(2)*(1d0 - c) - unit_axis(3)*s
      rotation(1,3) = unit_axis(3)*unit_axis(1)*(1d0 - c) + unit_axis(2)*s
      rotation(2,3) = unit_axis(2)*unit_axis(3)*(1d0 - c) - unit_axis(1)*s

      rotation(2,1) = unit_axis(1)*unit_axis(2)*(1d0 - c) + unit_axis(3)*s
      rotation(3,1) = unit_axis(3)*unit_axis(1)*(1d0 - c) - unit_axis(2)*s
      rotation(3,2) = unit_axis(2)*unit_axis(3)*(1d0 - c) + unit_axis(1)*s

      transformation%matrix = matmul(rotation, transformation%matrix)
      transformation%inverse_matrix = matmul(transformation%inverse_matrix, transpose(rotation))
   end subroutine cg_transformation_add_rotation_cos_sin

   !> Add a translation
   !!
   !! Compose by the following translation matrix:
   !! @verbatim
   !! ┌─────┬─────┬─────┬─────┐
   !! │  1  │  0  │  0  │ tx  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  1  │  0  │ ty  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  0  │  1  │ tz  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  0  │  0  │  1  │
   !! └─────┴─────┴─────┴─────┘
   !! @endverbatim
   !! Where tx, ty, tz are the coordinates of the translation @p vector.
   !!
   !! @param[inout] transformation: transformation that receive the transformation
   !! @param[in]    vector: translation vector
   !! @ingroup transformations
   pure subroutine cg_transformation_add_translation(transformation, vector)
      class(cg_transformation), intent(inout) :: transformation
      double precision, dimension(:), intent(in) :: vector

      double precision, dimension(size(vector)+1) :: T
      integer :: i

      do i = 1, transformation%dimension
         transformation%matrix(i,transformation%dimension+1) = transformation%matrix(i,transformation%dimension+1) + vector(i)
      end do

      T = 1d0
      T(1:transformation%dimension) = -vector

      transformation%inverse_matrix(:transformation%dimension+1,transformation%dimension+1) = &
         matmul(transformation%inverse_matrix(:transformation%dimension+1,:transformation%dimension+1), T)
   end subroutine cg_transformation_add_translation

   !> Add a scaling factor in direction x
   !!
   !! Compose by the following translation matrix:
   !! @verbatim
   !! ┌─────┬─────┬─────┬─────┐
   !! │ sx  │  0  │  0  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  1  │  0  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  0  │  1  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  0  │  0  │  1  │
   !! └─────┴─────┴─────┴─────┘
   !! @endverbatim
   !! Where sx is the scale factor in the direction x.
   !!
   !! @param[inout] transformation: transformation that receive the transformation
   !! @param[in]    factor: scaling factor in direction x.
   !! @ingroup transformations
   pure subroutine cg_transformation_add_scale_x(transformation, factor)
      class(cg_transformation), intent(inout) :: transformation
      double precision, intent(in) :: factor

      transformation%matrix(1,:) = transformation%matrix(1,:)*factor
      transformation%inverse_matrix(:,1) = transformation%inverse_matrix(:,1)/factor
   end subroutine cg_transformation_add_scale_x

   !> Add a scaling factor in direction y
   !!
   !! Compose by the following translation matrix:
   !! @verbatim
   !! ┌─────┬─────┬─────┬─────┐
   !! │  1  │  0  │  0  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │ sy  │  0  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  0  │  1  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  0  │  0  │  1  │
   !! └─────┴─────┴─────┴─────┘
   !! @endverbatim
   !! Where sy is the scale factor in the direction y.
   !!
   !! @param[inout] transformation: transformation that receive the transformation
   !! @param[in]    factor: scaling factor in direction y.
   !! @ingroup transformations
   pure subroutine cg_transformation_add_scale_y(transformation, factor)
      class(cg_transformation), intent(inout) :: transformation
      double precision, intent(in) :: factor

      transformation%matrix(2,:) = transformation%matrix(2,:)*factor
      transformation%inverse_matrix(:,2) = transformation%inverse_matrix(:,2)/factor
   end subroutine cg_transformation_add_scale_y

   !> Add a scaling factor in direction z
   !!
   !! Compose by the following translation matrix:
   !! @verbatim
   !! ┌─────┬─────┬─────┬─────┐
   !! │  1  │  0  │  0  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  1  │  0  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  0  │ sz  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  0  │  0  │  1  │
   !! └─────┴─────┴─────┴─────┘
   !! @endverbatim
   !! Where sz is the scale factor in the direction z.
   !!
   !! @warning Do not apply this transformation in 2D
   !!
   !! @param[inout] transformation: transformation that receive the transformation
   !! @param[in]    factor: scaling factor in direction z.
   !! @ingroup transformations
   pure subroutine cg_transformation_add_scale_z(transformation, factor)
      class(cg_transformation), intent(inout) :: transformation
      double precision, intent(in) :: factor

      transformation%matrix(3,:) = transformation%matrix(3,:)*factor
      transformation%inverse_matrix(:,3) = transformation%inverse_matrix(:,3)/factor
   end subroutine cg_transformation_add_scale_z

   !> Add a scaling in every direction
   !!
   !! Compose by the following translation matrix:
   !! @verbatim
   !! ┌─────┬─────┬─────┬─────┐
   !! │ sx  │  0  │  0  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │ sy  │  0  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  0  │ sz  │  0  │
   !! ├─────┼─────┼─────┼─────┤
   !! │  0  │  0  │  0  │  1  │
   !! └─────┴─────┴─────┴─────┘
   !! @endverbatim
   !! Where sx, sy, and sz are the scaling factor in direction x, y, and z.
   !!
   !! @param[inout] transformation: transformation that receive the transformation
   !! @param[in]    vector: array containing the scaling factors per direction.
   !! @ingroup transformations
   pure subroutine cg_transformation_add_scale(transformation, vector)
      class(cg_transformation), intent(inout) :: transformation
      double precision, dimension(:), intent(in) :: vector

      integer :: i

      do i = 1, transformation%dimension
         transformation%matrix(i,:transformation%dimension+1) = transformation%matrix(i,:transformation%dimension+1)*vector(i)
         transformation%inverse_matrix(:transformation%dimension+1,i) = &
            transformation%inverse_matrix(:transformation%dimension+1,i)/vector(i)
      end do
   end subroutine cg_transformation_add_scale

   !> Apply a transformation to a point
   !!
   !! Return the transformed point.
   !!
   !! @param[in] transformation: transformation
   !! @param[in] point: point to transform
   !! @ingroup transformations
   pure function cg_transform_point(transformation, point) result(r)
      class(cg_transformation), intent(in) :: transformation
      double precision, dimension(:), intent(in) :: point
      double precision, dimension(size(point)) :: r

      associate(m => transformation%matrix)
         if (transformation%dimension == 3) then
            r(1) = m(1,1)*point(1) + m(1,2)*point(2) + m(1,3)*point(3) + m(1,4)
            r(2) = m(2,1)*point(1) + m(2,2)*point(2) + m(2,3)*point(3) + m(2,4)
            r(3) = m(3,1)*point(1) + m(3,2)*point(2) + m(3,3)*point(3) + m(3,4)
         else
            r(1) = m(1,1)*point(1) + m(1,2)*point(2) + m(1,3)
            r(2) = m(2,1)*point(1) + m(2,2)*point(2) + m(2,3)
            ! Avoid non-initialized values.
            if (size(r) > 2) r(3:) = 0d0
         end if
      end associate
   end function cg_transform_point

   !> Apply a transformation to a point
   !!
   !! Elemental version of cg_transform_point. Only in 3D.
   !!
   !! @param[in] transformation: transformation
   !! @param[in] x0,y0,z0: coordinates of the original point
   !! @param[out] x,y,z: coordinates of the transformed point
   !! @ingroup transformations
   elemental subroutine cg_transform_point_elemental(transformation, x0,y0,z0, x,y,z)
      type(cg_transformation), intent(in) :: transformation
      double precision, intent(in) :: x0,y0,z0
      double precision, intent(out) :: x,y,z

      x = transformation%matrix(1, 1) * x0 + transformation%matrix(1, 2) * y0 + transformation%matrix(1, 3) * z0
      y = transformation%matrix(2, 1) * x0 + transformation%matrix(2, 2) * y0 + transformation%matrix(2, 3) * z0
      z = transformation%matrix(3, 1) * x0 + transformation%matrix(3, 2) * y0 + transformation%matrix(3, 3) * z0

      x = x + transformation%matrix(1, 4)
      y = y + transformation%matrix(2, 4)
      z = z + transformation%matrix(3, 4)
   end subroutine cg_transform_point_elemental

   !> Apply a transformation to a direction
   !!
   !! Return the transformed direction.
   !!
   !! @param[in] transformation: transformation
   !! @param[in] direction: direction to transform
   !! @ingroup transformations
   pure function cg_transform_direction(transformation, direction) result(r)
      class(cg_transformation), intent(in) :: transformation
      double precision, dimension(:), intent(in) :: direction
      double precision, dimension(size(direction)) :: r

      associate(m => transformation%matrix)
         if (transformation%dimension == 3) then
            r(1) = m(1,1)*direction(1) + m(1,2)*direction(2) + m(1,3)*direction(3)
            r(2) = m(2,1)*direction(1) + m(2,2)*direction(2) + m(2,3)*direction(3)
            r(3) = m(3,1)*direction(1) + m(3,2)*direction(2) + m(3,3)*direction(3)
         else
            r(1) = m(1,1)*direction(1) + m(1,2)*direction(2)
            r(2) = m(2,1)*direction(1) + m(2,2)*direction(2)
            ! Avoid non-initialized values.
            if (size(r) > 2) r(3:) = 0d0
         end if
      end associate
   end function cg_transform_direction

   !> Apply a transformation to a direction
   !!
   !! Elemental version of cg_transform_direction. Only in 3D.
   !!
   !! @param[in] transformation: transformation
   !! @param[in] x0,y0,z0: coordinates of the original point
   !! @param[out] x,y,z: coordinates of the transformed point
   !! @ingroup transformations
   elemental subroutine cg_transform_direction_elemental(transformation, x0,y0,z0, x,y,z)
      type(cg_transformation), intent(in) :: transformation
      double precision, intent(in) :: x0,y0,z0
      double precision, intent(out) :: x,y,z

      x = transformation%matrix(1, 1)*x0 + transformation%matrix(1, 2)*y0 + transformation%matrix(1, 3)*z0
      y = transformation%matrix(2, 1)*x0 + transformation%matrix(2, 2)*y0 + transformation%matrix(2, 3)*z0
      z = transformation%matrix(3, 1)*x0 + transformation%matrix(3, 2)*y0 + transformation%matrix(3, 3)*z0
   end subroutine cg_transform_direction_elemental

   !> Apply an inverse transformation to a point
   !!
   !! Return the transformed point.
   !!
   !! @param[in] transformation: transformation
   !! @param[in] point: point to transform
   !! @ingroup transformations
   pure function cg_inverse_transform_point(transformation, point) result(r)
      class(cg_transformation), intent(in) :: transformation
      double precision, dimension(:), intent(in) :: point
      double precision, dimension(size(point)) :: r

      associate(m => transformation%inverse_matrix)
         if (transformation%dimension == 3) then
            r(1) = m(1,1)*point(1) + m(1,2)*point(2) + m(1,3)*point(3) + m(1,4)
            r(2) = m(2,1)*point(1) + m(2,2)*point(2) + m(2,3)*point(3) + m(2,4)
            r(3) = m(3,1)*point(1) + m(3,2)*point(2) + m(3,3)*point(3) + m(3,4)
         else
            r(1) = m(1,1)*point(1) + m(1,2)*point(2) + m(1,3)
            r(2) = m(2,1)*point(1) + m(2,2)*point(2) + m(2,3)
            ! Avoid non-initialized values.
            if (size(r) > 2) r(3:) = 0d0
         end if
      end associate
   end function cg_inverse_transform_point

   !> Apply an inverse transformation to a point
   !!
   !! Elemental version of cg_inverse_transform_point. Only in 3D.
   !!
   !! @param[in] transformation: transformation
   !! @param[in] x0,y0,z0: coordinates of the original point
   !! @param[out] x,y,z: coordinates of the transformed point
   !! @ingroup transformations
   elemental subroutine cg_inverse_transform_point_elemental(transformation, x0,y0,z0, x,y,z)
      type(cg_transformation), intent(in) :: transformation
      double precision, intent(in) :: x0,y0,z0
      double precision, intent(out) :: x,y,z

      x = transformation%inverse_matrix(1, 1)*x0 + transformation%inverse_matrix(1, 2)*y0 + transformation%inverse_matrix(1, 3)*z0
      y = transformation%inverse_matrix(2, 1)*x0 + transformation%inverse_matrix(2, 2)*y0 + transformation%inverse_matrix(2, 3)*z0
      z = transformation%inverse_matrix(3, 1)*x0 + transformation%inverse_matrix(3, 2)*y0 + transformation%inverse_matrix(3, 3)*z0

      x = x + transformation%inverse_matrix(1, 4)
      y = y + transformation%inverse_matrix(2, 4)
      z = z + transformation%inverse_matrix(3, 4)
   end subroutine cg_inverse_transform_point_elemental

   !> Apply an inverse transformation to a direction
   !!
   !! Return the transformed direction.
   !!
   !! @param[in] transformation: transformation
   !! @param[in] direction: direction to transform
   !! @ingroup transformations
   pure function cg_inverse_transform_direction(transformation, direction) result(r)
      class(cg_transformation), intent(in) :: transformation
      double precision, dimension(:), intent(in) :: direction
      double precision, dimension(size(direction)) :: r

      associate(m => transformation%inverse_matrix)
         if (transformation%dimension == 3) then
            r(1) = m(1,1)*direction(1) + m(1,2)*direction(2) + m(1,3)*direction(3)
            r(2) = m(2,1)*direction(1) + m(2,2)*direction(2) + m(2,3)*direction(3)
            r(3) = m(3,1)*direction(1) + m(3,2)*direction(2) + m(3,3)*direction(3)
         else
            r(1) = m(1,1)*direction(1) + m(1,2)*direction(2)
            r(2) = m(2,1)*direction(1) + m(2,2)*direction(2)
            ! Avoid non-initialized values.
            if (size(r) > 2) r(3:) = 0d0
         end if
      end associate
   end function cg_inverse_transform_direction

   !> Apply an inverse transformation to a direction
   !!
   !! Elemental version of cg_inverse_transform_direction. Only in 3D.
   !!
   !! @param[in] transformation: transformation
   !! @param[in] x0,y0,z0: coordinates of the original point
   !! @param[out] x,y,z: coordinates of the transformed point
   !! @ingroup transformations
   elemental subroutine cg_inverse_transform_direction_elemental(transformation, x0,y0,z0, x,y,z)
      type(cg_transformation), intent(in) :: transformation
      double precision, intent(in) :: x0,y0,z0
      double precision, intent(out) :: x,y,z

      x = transformation%inverse_matrix(1, 1)*x0 + transformation%inverse_matrix(1, 2)*y0 + transformation%inverse_matrix(1, 3)*z0
      y = transformation%inverse_matrix(2, 1)*x0 + transformation%inverse_matrix(2, 2)*y0 + transformation%inverse_matrix(2, 3)*z0
      z = transformation%inverse_matrix(3, 1)*x0 + transformation%inverse_matrix(3, 2)*y0 + transformation%inverse_matrix(3, 3)*z0
   end subroutine cg_inverse_transform_direction_elemental

end module mod_cg_transformation

!This file is part of Notus 0.4.0

!Copyright Bordeaux-INP, Université de Bordeaux, CNRS
!
!Contributors:
!Antoine Lemoine, 08-08-2015, antoine.lemoine@bordeaux-inp.fr

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

!> Module containing the parameters of the Moment-of-Fluid method
!!
!! @ingroup moment_of_fluid
module variables_mof
   !> Tolerance angle of the minimization algorithm in 2D
   double precision :: mof2d_tol_angle = 1d-8
   !> Tolerance of the derivative of the minimization algorithm in 2D
   double precision :: mof_tol_derivative = 0d0
   !> Tolerance of the derivative of the minimization algorithm in 3D
   double precision :: mof3d_tol_derivative = 1d-8
   !> Maximum number of iteration for minimization algorithm in 3D
   integer :: mof3d_max_iter = 100
   !> CFL for the MOF advection scheme
   double precision :: mof_cfl = 0.9d0
   !> Flag to enable the analytic reconstruction on rectangular and hexahedral cells
   logical :: mof_use_analytic_reconstruction = .true.
   !> Flag to enable the usage of the B-tree in the reconstruction (nb_phases >= 4)
   logical :: mof_use_b_tree = .false.
   !> Flag to enable symmetric reconstruction
   logical :: mof_use_symmetric_reconstruction = .true.
   !> Flag to enable the use of filaments
   logical :: mof_has_filaments = .false.
   !> Maximum number of filament parts per phase
   integer :: mof_max_filament_parts = 2
   !> Flag to use analytical derivatives (Chen & Zhang 2016)
   logical :: mof3d_use_analytical_derivatives = .false.
   !> Flag to use optimized centroid (Milcent & Lemoine 2019)
   logical :: mof3d_use_optimized_centroid = .true.
   !> Flag to use analytic formulas for centroid and gradient (Milcent & Lemoine 2019)
   logical :: mof3d_internal_is_analytic_gradient_enabled = .true.
   !> Flag to activate or deactivate the Navier boundary condition for MoF 
   logical :: mof2d_use_navier_boundary_condition = .false.
   !> Value of the Navier boundary coefficient (beta) for MoF
   double precision :: mof2d_navier_boundary_coefficient = 0.75d0

   !> Toggle filter
   logical :: mof_use_filter = .false.
   !> Filter the isolated structures
   double precision :: mof_filter_min_isolated_volume_fraction = 1d-2
   !> Smallest authorized volume fraction for filaments
   double precision :: mof_filament_min_volume_fraction = 1d-3

   !> Toggle volume difference diagnostic.
   logical :: mof_has_volume_difference = .false.

   !> Smallest authorized volume fraction
   double precision, parameter :: MOF_VOLUME_FRACTION_EPSILON = 100d0*epsilon(1d0)
end module variables_mof

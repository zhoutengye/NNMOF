Module Decision_Tree

    Implicit None
    Type Decision_Tree_Binary
    Integer :: node
    Integer :: children_left
    Integer :: children_right
    Integer :: feature
    Real(8) :: threshold
    Real(8) :: value(3)
    End Type Decision_Tree_Binary

    Integer :: node_count
    Integer :: n_features
    Integer :: n_outputs
    Integer :: max_depth

    Type(Decision_Tree_Binary), Allocatable :: DTree(:)

Contains

!---------------------------------------------------------------------------------
!
! Read the Binary tree
!
!---------------------------------------------------------------------------------
! 
! The binary tree from Python sklearn library is stored with pre-order traversal.
! More details see: https://en.wikipedia.org/wiki/Tree_traversal
! 
! For each node in decision tree, 5 parameters are stored. They are:
!   1. children_left  : integer, corresponds with the node number of its left child
!   2. children_right : integer, corresponds with the node number of its left child
!   4. feature        : integer, indicates which feature used for threshold
!   3. threshold      : float,   the creterion for choosing left or right. 
!                                When input(feature) < threshold, go LEFT
!                                When input(feature) >= threshold, go RIGHT
!   5. value(:)       : float, corresponds with the results at each leaf node.
!                       the size of value(:) at each node equals to the number of 
!                       output, however, FORTRAN cannot make allocatable variable 
!                       in type, so value is assign as the size of 10. If more than 
!                       10 outputs are needed, the size should be increased.
!
!---------------------------------------------------------------------------------
!
! The DATA/tree_parameters.namelist file is used to store parameters of the tree, 
! including:
!       node_count : number of total nodes in binary tree
!       n_features : number of input parameters
!       n_outputs  : number of output parameters
!       max_depth  : maximum levels of binary tree
!
! The binary tree data generated from Python code are stored in the following 
! data files:
!       DATA/tree_left.dat      : 
!       DATA/tree_right.dat     :
!       DATA/tree_feature.dat   :
!       DATA/tree_threshold.dat :
!       DATA/tree_value**.dat'  : the number of value files depends on the number of outpurs
!
!---------------------------------------------------------------------------------
!
! By Zhouteng Ye
! Update: August 23th, 2018 
!---------------------------------------------------------------------------------
    Subroutine Read_Decision_Tree

        Implicit None

        Namelist /tree_Parameters/ node_count, n_features, n_outputs, max_depth

        Integer :: kk
        integer :: i
        Character(50) :: ii

        ! Read the parameters for decision tree
        Open(61,file='DATA/tree_parameters.namelist')
        Read(61,nml=tree_Parameters)

        ! Allocate variables according to the tree parameters
        Allocate(DTree(node_count))

        ! Read tree data
        Open(62,file='DATA/tree_left.dat')
        Open(63,file='DATA/tree_right.dat')
        Open(64,file='DATA/tree_feature.dat')
        Open(65,file='DATA/tree_threshold.dat')
        kk = 66
        Do i = 1, n_outputs
            write(ii,'(I8)') i
            Open(kk,file='DATA/tree_value'//trim(adjustl(ii))//'.dat',status = 'old')
            kk = kk+1
        End Do

        ! Read decision tree data
        Read(62,*)  DTree%children_left
        Read(63,*)  DTree%children_right
        Read(64,*)  DTree%feature
        Read(65,*)  DTree%threshold
        kk = 66
        Do i = 1,3
            Read(kk,'(F12.8)')  DTree%value(i)
            kk = kk+1
        End Do

    End Subroutine Read_Decision_Tree


!---------------------------------------------------------------------------------
!
! Search the Binary tree
! 
!
!
!---------------------------------------------------------------------------------
!
! By Zhouteng Ye
! Update: August 23th, 2018 
!---------------------------------------------------------------------------------

    Function Predict_Decision_Tree(input) Result(output)

        Implicit None

        Real(8) :: input(:)
        Real(8) :: output(n_outputs)
        Integer :: node

        Integer :: i

        ! Find the predicted results with binary tree search
        node = 1
        Do i = 1,max_depth

            !! This line can be used to show the decision path
            ! write(*,*) node, DTree(node)%feature, DTree(node)%threshold, DTree(node)%children_left, DTree(node)%children_right
        
            If ( DTree(node)%feature .eq. -2) Exit

            If ( input(DTree(node)%feature+1) .le. DTree(node)%threshold ) Then
                node = DTree(node)%children_left +1
            Else
                node = DTree(node)%children_right +1
            EndIf

        EndDo

        output = DTree(node)%value(1:n_outputs)

    End Function Predict_Decision_Tree

!---------------------------------------------------------------------------------
!
! Search the Binary tree
! 
! return norm
!
!---------------------------------------------------------------------------------
!
! By Zhouteng Ye
! Update: September 3th, 2018 
!---------------------------------------------------------------------------------

    Function Predict_Norm_Decision_Tree(input) Result(output)

        Implicit None

        Real(8) :: input(:)
        Real(8) :: output(n_outputs-1)
        Integer :: node

        Integer :: i

        ! Find the predicted results with binary tree search
        node = 1
        Do i = 1,max_depth

            !! This line can be used to show the decision path
            ! write(*,*) node, DTree(node)%feature, DTree(node)%threshold, DTree(node)%children_left, DTree(node)%children_right
        
            If ( DTree(node)%feature .eq. -2) Exit

            If ( input(DTree(node)%feature+1) .le. DTree(node)%threshold ) Then
                node = DTree(node)%children_left +1
            Else
                node = DTree(node)%children_right +1
            EndIf

        EndDo

        output = DTree(node)%value(1:n_outputs-1)

    End Function Predict_Norm_Decision_Tree

!---------------------------------------------------------------------------------
!
! Search the Binary tree
! 
! return intercepts
!
!---------------------------------------------------------------------------------
!
! By Zhouteng Ye
! Update: September 3th, 2018 
!---------------------------------------------------------------------------------

    Function Predict_Intercepts_Decision_Tree(input) Result(output)

        Implicit None

        Real(8) :: input(:)
        Real(8) :: output
        Integer :: node

        Integer :: i

        ! Find the predicted results with binary tree search
        node = 1
        Do i = 1,max_depth

            !! This line can be used to show the decision path
            ! write(*,*) node, DTree(node)%feature, DTree(node)%threshold, DTree(node)%children_left, DTree(node)%children_right
        
            If ( DTree(node)%feature .eq. -2) Exit

            If ( input(DTree(node)%feature+1) .le. DTree(node)%threshold ) Then
                node = DTree(node)%children_left +1
            Else
                node = DTree(node)%children_right +1
            EndIf

        EndDo

        output = DTree(node)%value(n_outputs)

    End Function Predict_Intercepts_Decision_Tree

End Module Decision_Tree

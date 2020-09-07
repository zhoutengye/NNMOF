Program Main
  Use NueralNetwork
  Implicit None
  real(8)  :: inputs(4), outputs(2)
  Type(Neural_Network) :: ml

  Call ml%Initialization

  ! inputs = [0.125,0.125,0.125,1.0/6.0]
  inputs = [ 0.22438486,  0.02104238, -0.05832734,  0.5269191 ]


  print *, 'n_inputs', ml%n_inputs
  print *, 'n_outputs', ml%n_outputs

  print *, ml%Predict(inputs)
End Program Main
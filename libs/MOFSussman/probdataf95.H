      INTEGER_T fort_max_level
      INTEGER_T bfact_space_order(0:100)
      INTEGER_T bfact_time_order
 
      REAL_T adv_vel
      REAL_T rgasinlet
      REAL_T vinletgas
      REAL_T denfact
      REAL_T xblobarr(10)
      REAL_T yblobarr(10)
      REAL_T zblobarr(10)
      REAL_T radblobarr(10)
      REAL_T xblob,yblob,zblob,radblob
      REAL_T xblob2,yblob2,zblob2,radblob2
      REAL_T xblob3,yblob3,zblob3,radblob3
      REAL_T xblob4,yblob4,zblob4,radblob4
      REAL_T xblob5,yblob5,zblob5,radblob5
      REAL_T xblob6,yblob6,zblob6,radblob6
      REAL_T xblob7,yblob7,zblob7,radblob7
      REAL_T xblob8,yblob8,zblob8,radblob8
      REAL_T xblob9,yblob9,zblob9,radblob9
      REAL_T xblob10,yblob10,zblob10,radblob10
      REAL_T xactive,yactive,zactive
      REAL_T ractivex,ractivey,ractivez
      REAL_T velfact
      REAL_T advbot
      REAL_T inflow_pressure,outflow_pressure
      REAL_T period_time
      REAL_T twall
      REAL_T problox,probloy,probloz
      REAL_T probhix,probhiy,probhiz
      INTEGER_T fort_is_periodic(3)
      REAL_T fort_stop_time
      REAL_T fort_energyconst(100)
      REAL_T fort_energycutoff(100)
      REAL_T fort_drhodt(100)
      REAL_T fort_drhodz(100)
      REAL_T fort_tempconst(100)
      REAL_T fort_initial_temperature(100)
      REAL_T fort_tempcutoff(100)
      REAL_T fort_tempcutoffmax(100)
      REAL_T fort_denconst(100)
      REAL_T fort_viscconst(100)
      INTEGER_T fort_viscosity_state_model(100)
      REAL_T fort_heatviscconst(100)

      REAL_T fort_prerecalesce_heatviscconst(100)
      REAL_T fort_prerecalesce_viscconst(100)
      REAL_T fort_prerecalesce_stiffCP(100)
      
      REAL_T fort_speciesconst(1000)
      REAL_T fort_speciesviscconst(1000)
      REAL_T fort_tension(1000)
      REAL_T fort_tension_slope(1000)
      REAL_T fort_tension_T0(1000)
      REAL_T fort_tension_min(1000)
      REAL_T fort_prefreeze_tension(1000)

      REAL_T fort_alpha(1000)
      REAL_T fort_beta(1000)
      REAL_T fort_expansion_factor(1000)
      REAL_T fort_time_radblob(1000)
      REAL_T fort_stefan_number(1000)
      REAL_T fort_jacob_number(1000)

      REAL_T fort_stiffPINF(100)
      REAL_T fort_stiffCP(100)
      REAL_T fort_stiffGAMMA(100)

      REAL_T FORT_MUSHY_THICK  ! default==2.0
      REAL_T density_at_depth
      REAL_T global_pressure_scale
      REAL_T global_velocity_scale
      REAL_T gravity

      INTEGER_T num_species_var
      INTEGER_T configuration_tensor_var
      INTEGER_T elastic_tensor_var
      INTEGER_T num_state_material
      INTEGER_T num_state_base
      INTEGER_T ngeom_raw
      INTEGER_T ngeom_recon
      INTEGER_T num_materials
      INTEGER_T num_materials_vel
      INTEGER_T num_materials_scalar_solve
      INTEGER_T fort_max_num_materials
      INTEGER_T fort_max_num_eos
      INTEGER_T probtype, adv_dir, axis_dir
      INTEGER_T levelrz
      INTEGER_T prescribe_temperature_outflow
      INTEGER_T invert_solid_levelset
      INTEGER_T FSI_flag(100)
      INTEGER_T fort_material_type(100)
      INTEGER_T gravity_dir,invert_gravity
      INTEGER_T pres_homflag
      INTEGER_T vel_homflag
      INTEGER_T temp_homflag
      INTEGER_T species_homflag
      INTEGER_T ls_homflag
      INTEGER_T fort_finest_level
      INTEGER_T n_sites
      REAL_T nucleation_init_time
      REAL_T pos_sites(400)


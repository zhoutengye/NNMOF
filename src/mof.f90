Module ModMOF
  Use ModGlobal
  Use SussmanMOF
  Implicit None
Contains
  Subroutine MOF_Init
    Implicit None
    Call Initialize_SussmanMOF
  End Subroutine MOF_Init

  Subroutine NormSussmanMOF(f, c, norm)
    Use geometry_intersect_module
    Use MOF_routines_module
    Implicit None
    Real(sp) , Intent(In)   :: f
    Real(sp) , Intent(In)   :: c(3)
    Real(sp) , Intent(out)  :: norm(3)

    mofdata = 0.0_sp
    mofdata(1) = f
    mofdata(2) = c(1) - 0.5
    mofdata(3) = c(2) - 0.5
    mofdata(4) = c(3) - 0.5
    mofdata(5) = 1.0_sp

    mofdata(10)  = 1.0-f
    mofdata(11) = - f * c(1) / mofdata(10)
    mofdata(12) = - f * c(2) / mofdata(10)
    mofdata(13) = - f * c(3) / mofdata(10)
    mofdata(14) = 2.0_sp

    Call multimaterial_MOF( &
        bfact,dl,xsten0,nhalf0, &
        mof_verbose, &
        use_ls_data, &
        LS_stencil, &
        xtetlist_vof, &
        xtetlist_cen, &
        nmax, &
        mofdata, &
        multi_centroidA, &
        continuous_mof, &
        levelrz,nmat,sdim, &
        ngeom_recon, &
        caller_id)

    norm(1) = -mofdata(6)
    norm(2) = -mofdata(7)
    norm(3) = -mofdata(8)

  End Subroutine NormSussmanMOF

End Module MODMOF

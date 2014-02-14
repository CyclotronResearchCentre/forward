// Run with
// getdp eeg_forward.pro -solve Electrostatics
//

// Get Real Valued GetDP: http://geuz.org/getdp/bin/MacOSX/getdp-svn-MacOSX64r.tgz
//
// Solve with these flags to save memory:
// getdp eeg_forward.pro -msh TMS007_running.msh -bin -solve Electrostatics -ksp_type gmres -pc_type ilu -pc_factor_levels 2 -ksp_gmres_restart 1000 -ksp_rtol 1e-10
// use -v2 to save output in v2 mesh format to inspect each element

// Opened eeg_forward.geo
// Defined 2D Mesh
// Saved as .msh
// Merged with anatomy (3D mesh)
// Saved as new .msh
// Ran with
// getdp eeg_forward.pro -solve Electrostatics -msh eeg_forward_both.msh
//
// Or in parallel with
// openmpirun -np 4 getdp eeg_forward.pro -pre Electrostatics -cal -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps -msh TMS007_all.msh
//
// getdp eeg_forward.pro -post Electrostatics_PostProcessing

/// Can't use "Skin" because it's a function name!!!

scale = 0.0001;
Group {
  // The PhysicalVolumes
  WhiteMatter_Cerebellum = Region[1001];
  GrayMatter = Region[1002];
  CSF_Ventricles = Region[1003];
  Skull = Region[1004];
  Scalp = Region[1005];
  Sink = Region[5008];
  Elec2 = Region[5001];
  Source = Region[8000];

  Omega = Region[{WhiteMatter_Cerebellum, GrayMatter, CSF_Ventricles, Skull, Scalp}];
  
  // Find a way to specify Region[5000:5060]? with python?
  Electrodes = Region[{Sink, Source}];
  DipoleRegions = Region[{Sink, Source}];
}

// Define the conductivities
// Values are in Siemens per metre [S / mm]
//
// See Gullmar 2010:
// 0.337 S/m for gray matter
// 0.0042 S/m for the skull compartment
// 0.33 S/m for soft tissue/skin
// (Baumann et al., 1997; Haueisen et al., 2002; Rullmann et al., 2009; Wolters, 2002).
Function {
    sigma[GrayMatter]=0.33;

    // while GmshRead is commented out!
    sigma[WhiteMatter_Cerebellum] = 0.33; //sigma[WhiteMatter_Cerebellum] = TensorField[XYZ[]] #1 ? #1 : 0.33 ;//: 0.33*1e-3 ;

    sigma[CSF_Ventricles] = 1.79;
    sigma[Skull]=0.0042;
    sigma[Scalp]=0.33;

    sigma[Sink]=0.33; // On scalp
    sigma[Source]=0.33; // In Gray Matter
    DefineConstant[ Length = 1. ] ;

}


/* --------------------------------------------------------------------------*/

Constraint {

  { Name ElectricScalarPotential ;
    Case { /* A reference must be given for the scalar potential */
      {
        Region Sink ; Value 0. ;
      }
    }
  }

  /* ... or the current is fixed. Uncomment one of these two... */
  { Name GlobalElectricCurrentSource ;
    Case {
      { Region Source ; Value 1. ; }
    }
  }
  { Name GlobalElectricCurrentSink ;
    Case {
      { Region Sink ; Value -1. ; }
    }
  }


}

Jacobian {
  { Name Volume ;
    Case {
      { Region All ; Jacobian Vol ; }
    }
  }
}


Integration {
  { Name GradGrad ;
    Case {{
      Type Gauss ;
        Case {
         {GeoElement Triangle; NumberOfPoints 4;}
         {GeoElement Tetrahedron; NumberOfPoints 1;}
       }
    }}
}}


FunctionSpace {{
  Name Hgrad_vf_Ele;
  Type Form0;

  BasisFunction {
    {
    Name sn;
    NameOfCoef vn;
    Function BF_Node;
    Support Region[{Omega}];
    Entity NodesOf[All, Not Electrodes];
    }
    //
    {
    Name sck ;
    NameOfCoef vck_source ;
    Function BF_GroupOfNodes ;
    Support Region[{Omega}];
    Entity GroupsOfNodesOf[ Source ];
    }
    //
    {
    Name sck ;
    NameOfCoef vck_sink ;
    Function BF_GroupOfNodes ;
    Support Region[{Omega}];
    Entity GroupsOfNodesOf[ Sink ];
    }
  }
  GlobalQuantity {
    {
    Name V_source;
    Type AliasOf;
    NameOfCoef vck_source;
    }
    {
    Name V_sink;
    Type AliasOf;
    NameOfCoef vck_sink;
    }
    {
    Name I_source;
    Type AssociatedWith;
    NameOfCoef vck_source;
    }
    {
    Name I_sink;
    Type AssociatedWith;
    NameOfCoef vck_sink;
    }

  }
    Constraint {
      { NameOfCoef vn ;
        EntityType NodesOf ; NameOfConstraint ElectricScalarPotential ; }
      { NameOfCoef V_sink ;
        EntityType GroupsOfNodesOf ; NameOfConstraint GlobalElectricPotential ; }
      { NameOfCoef I_source ;
        EntityType GroupsOfNodesOf ; NameOfConstraint GlobalElectricCurrentSource ; }
      { NameOfCoef I_sink ;
        EntityType GroupsOfNodesOf ; NameOfConstraint GlobalElectricCurrentSink ; }
    }

}}


Formulation {{
  Name Electrostatics_Formulation;
  Type FemEquation;
  Quantity {
    {
    Name v;
    Type Local;
    NameOfSpace Hgrad_vf_Ele;
    }
    {
    Name I_source;
    Type Global;
    NameOfSpace Hgrad_vf_Ele [I_source];
    }
    {
    Name I_sink;
    Type Global;
    NameOfSpace Hgrad_vf_Ele [I_sink];
    }
    {
    Name V_source;
    Type Global;
    NameOfSpace Hgrad_vf_Ele [V_source];
    }
    {
    Name V_sink;
    Type Global;
    NameOfSpace Hgrad_vf_Ele [V_sink];
    }
  }
  Equation {
    Galerkin {
      [sigma[] * Dof{Grad v}, {Grad v}];
      In Omega;
      Jacobian Volume;
      Integration GradGrad;
    }
    GlobalTerm
    {
      [ Dof{I_source} / Length , {V_source} ];
      In Source ;
    }
    GlobalTerm
    {
      [ Dof{I_sink} / Length , {V_sink} ];
      In Sink ;
    }
  }
}}


Resolution {{
  Name Electrostatics ;
  System {{
    Name Electrostatic_System;
    NameOfFormulation Electrostatics_Formulation;
  }}
  Operation {
    SetGlobalSolverOptions["-ksp_type gmres -ksp_gmres_restart 1000 -ksp_rtol 1e-10"];
    SetGlobalSolverOptions["-pc_type ilu -pc_factor_levels 2 "];
    Generate[Electrostatic_System];
    Solve[Electrostatic_System];
    SaveSolution[Electrostatic_System];
    PostOperation[v_j_e];
  }
}}


PostProcessing {{
  Name Electrostatics_PostProcessing;
  NameOfFormulation Electrostatics_Formulation;
  NameOfSystem Electrostatic_System;
  PostQuantity {
    { Name v; Value {Term {[{v}]; In Omega; Jacobian Volume;}}}
    { Name j; Value {Term {[sigma[] * (-{Grad v}) ]; In Omega; Jacobian Volume;}}}
    { Name e; Value {Term {[-{Grad v}]; In Omega; Jacobian Volume;}}}

    { Name V_source; Value { Term { [ {V_source} ]; In Source; Jacobian Volume;} } }
    { Name I_source; Value { Term { [ {I_source} ]; In Source; Jacobian Volume;} } }
    { Name R_source; Value { Term { [ -{V_source}/{I_source} ]; In Source; Jacobian Volume;} } }

    { Name V_sink; Value { Term { [ {V_sink} ]; In Sink; Jacobian Volume;} } }
    { Name I_sink; Value { Term { [ {I_sink} ]; In Sink; Jacobian Volume;} } }
    { Name R_sink; Value { Term { [ -{V_sink}/{I_sink} ]; In Sink; Jacobian Volume;} } }
    // Sanity check for electrode potential
    {Name v_elec; Value {Term {[{v}]; In Electrodes; Jacobian Volume;}}}
    {Name e_brain; Value {Term {[-{Grad v}]; In GrayMatter; Jacobian Volume;}}}
  }
}}


PostOperation v_j_e UsingPost Electrostatics_PostProcessing {
  Print [v, OnElementsOf Omega, File "v.pos"];
  Print [j, OnElementsOf Omega, File "j.pos"];
  Print [e, OnElementsOf Omega, File "e.pos"];

  Print [V_source, OnElementsOf Source, File "V_source.pos"];
  Print [I_source, OnElementsOf Source, File "I_source.pos"];
  Print [R_source, OnElementsOf Source, File "R_source.pos"];

  Print [V_sink, OnElementsOf Sink, File "V_sink.pos"];
  Print [I_sink, OnElementsOf Sink, File "I_sink.pos"];
  Print [R_sink, OnElementsOf Sink, File "R_sink.pos"];
  
  Print [v_elec, OnElementsOf Electrodes, Depth 0, Format SimpleTable, File "v_elec.txt"];
  Print [e_brain, OnElementsOf GrayMatter,  Depth 0, Format SimpleTable, File "e_brain.txt" ];
}
//End of File
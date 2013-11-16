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
  Anode = Region[5008];
  Cathode = Region[5011];

  Omega = Region[{WhiteMatter_Cerebellum, GrayMatter, CSF_Ventricles, Skull, Scalp}];
  Electrodes = Region[{Anode, Cathode}];
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

    sigma[WhiteMatter_Cerebellum] = TensorField[XYZ[]] #1 ? #1 : 0.33 ;//: 0.33*1e-3 ;
    //sigma[WhiteMatter_Cerebellum] = 0.33;

    sigma[CSF_Ventricles] = 1.79;
    sigma[Skull]=0.0042;
    sigma[Scalp]=0.33;
    sigma[Anode]=0.33;
    sigma[Cathode]=0.33;
}


Constraint {{
  Name ElectricScalarPotential;
  Type Assign;
  Case {
    {Region Region[{Anode}]; Value  1.;}
    {Region Region[{Cathode}]; Value -1.;}
  }
}}



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
  BasisFunction {{
    Name sn;
    NameOfCoef vn;
    Function BF_Node;
    Support Region[{Omega}];
    Entity NodesOf[All];
  }}
  Constraint {{
    NameOfCoef vn;
    EntityType NodesOf;
    NameOfConstraint ElectricScalarPotential;
  }}
}}


Formulation {{
  Name Electrostatics_Formulation;
  Type FemEquation;
  Quantity {{
    Name v;
    Type Local;
    NameOfSpace Hgrad_vf_Ele;
  }}
  Equation {
    Galerkin {
      [sigma[] * Dof{Grad v}, {Grad v}];
      In Omega;
      Jacobian Volume;
      Integration GradGrad;
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
    SetGlobalSolverOptions["-ksp_type gmres -ksp_gmres_restart 1000 -ksp_rtol 1e-8"];
    SetGlobalSolverOptions["-pc_type ilu -pc_factor_levels 2 "];
    GmshRead["TMS007_running.msh"];
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
    {Name v; Value {Term {[{v}]; In Omega; Jacobian Volume;}}}
    {Name j; Value {Term {[sigma[] * (-{Grad v}) ]; In Omega; Jacobian Volume;}}}
    {Name e; Value {Term {[-{Grad v}]; In Omega; Jacobian Volume;}}}
    // Sanity check for electrode potential
    {Name v_elec; Value {Term {[{v}]; In Electrodes; Jacobian Volume;}}}

    // Sanity check for electrode potential
    {Name e_brain; Value {Term {[-{Grad v}]; In GrayMatter; Jacobian Volume;}}}
  }
}}


PostOperation v_j_e UsingPost Electrostatics_PostProcessing {
  Print [v, OnElementsOf Omega, File "v_eeg_forward.pos"];
  Print [v_elec, OnElementsOf Electrodes, File "v_elec.pos"];

  // Current density is in amperes per square metre [A / m^2]
  //Print [j, OnElementsOf Omega, File "j_eeg_forward.pos"];
  //Print [e, OnElementsOf Omega, File "e_eeg_forward.pos"];

  // Save the electric field elements in a text file
  // This is done for the reciprocity calculations later
  // Depth 0 gives returns values in the barycenter of each element
  // Electric field is in volts per metre [V / M]
  //Print [e, OnElementsOf Omega, Depth 0, Format Table, File "e_eeg_forward.txt" ];

  // Electric field is in volts per metre [V / M]
  //Print [e_brain, OnElementsOf GrayMatter, Depth 0, Format Table, File "e_eeg_forward.txt" ];
  Print [e_brain, OnElementsOf GrayMatter, File "e_eeg_forward.pos" ];

  //Print [v, OnElementsOf Omega, Format SimpleTable, File "v_eeg_forward.txt" ];

  // Sanity check for electrode potential
  // V is in Volts [V]
  Print [v_elec, OnElementsOf Electrodes, Format Table, File "v_electrodes.txt" ];

  // Sanity check for electrode potential
  //Print [v_elec, OnElementsOf Electrodes, Depth 0, Format SimpleTable, File "v_electrodes_d0.txt" ];

}

Group {
  // The PhysicalVolumes
  GrayMatter = Region[1001];
  CSF_Ventricles = Region[1002];
  Skull = Region[1003];
  Scalp = Region[1004];
  Anode = Region[5008];
  Cathode = Region[5011];

  Omega = Region[{GrayMatter, CSF_Ventricles, Skull, Scalp}];
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
    {Name v; Value {Term {[{v}]; In Omega; Jacobian Volume;}}}
    {Name j; Value {Term {[sigma[] * (-{Grad v}) ]; In Omega; Jacobian Volume;}}}
    {Name e; Value {Term {[-{Grad v}]; In Omega; Jacobian Volume;}}}
    // Sanity check for electrode potential
    //{Name v_elec; Value {Term {[{v}]; In Electrodes; Jacobian Volume;}}}
    //{Name e_sphere; Value {Term {[-{Grad v}]; In GrayMatter; Jacobian Volume;}}}
  }
}}


PostOperation v_j_e UsingPost Electrostatics_PostProcessing {
  Print [v, OnElementsOf Omega, File "v_sphere.pos"];
  Print [j, OnElementsOf Omega, File "j_sphere.pos"];
  Print [e, OnElementsOf Omega, File "e_sphere.pos"];
  //Print [v_elec, OnElementsOf Electrodes, File "v_elec_sphere.pos"];
  //Print [e_sphere, OnElementsOf GrayMatter, File "e_sphere.pos" ];
}
//End of File
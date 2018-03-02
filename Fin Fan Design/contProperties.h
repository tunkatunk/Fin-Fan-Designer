//
//  contProperties.h
//  Semi Fin Fan Design
//
//  Created by Jay Tunkie Saunders on 2/27/18.
//  Copyright © 2018 Jay Tunkie Saunders. All rights reserved.
//

#ifndef contProperties_h
#define contProperties_h

char state = 'c';

// Inlet Process Fluid Properties From Aspen HYSYS Semi Simulation
const float Tp_in = 271.8; // F
const float p_in = 122.4; // psia
//const float m = 2.249e+5; // lb/hr
//const float m = 2.249e+5; // lb/hr
const float Cp_in = 0.6201; // Btu/lb/F
const float rho_in = 0.4678; // lb/ft^3
const float mu_l_in = 0.5317; // lb/ft/h from HYSYS 0.2198 cP
const float mu_v_in = 0.01713; // lb/ft/h from HYSYS 0.708e-02 cP
const float vf_l_in = 2.989e-03; // dimless
const float vf_v_in = 0.997; // dimless
const float H_in = -217.6; // Btu/lb from HYSYS

// Outlet Process Fluid Properties
const float Tp_out = 104.0; // F
const float Cp_out = 0.5833; // Btu/lb/F
const float rho_out = 0.5993; // lb/ft^3
const float mu_l_out = 0.8999; // lb/ft/h : from HYSYS 0.4561 cP
const float mu_v_out = 0.026343; // lb/ft/h : from HYSYS 1.01e-002 cP
const float vf_l_out = 9.699e-003; // dimless
const float vf_v_out = 0.9903; // dimless
const float H_out = -391.8; // Btu/lb from HYSYS
/*
 // Average Process Fluid Properties
 const float Cp = (Cp_in + Cp_out) / 2;
 const float rho = (rho_in + rho_out) / 2;
 const float mu_in = vf_l_in*mu_l_in + vf_v_in*mu_v_in;
 const float mu_out = vf_l_out*mu_l_out + vf_v_out*mu_v_out;
 const float mu = (mu_in + mu_out) / 2;
 const float k = 0.073736; // Btu/hr/ft/F*/

// Process Fluid Properties At Mean Temperature
const float Cp = 0.6141; // Btu/lb/F
const float rho = 0.5825; // lb/ft^3
const float mu = 0.3033; // lb/ft/hr
const float k = 0.073736; // Btu/hr/ft/F

#endif /* contProperties_h */


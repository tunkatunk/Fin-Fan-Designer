//
//  Semi_properties.h
//  Semi Fin Fan Design
//
//  Created by Jay Tunkie Saunders on 2/27/18.
//  Copyright © 2018 Jay Tunkie Saunders. All rights reserved.
//

#ifndef Semi_properties_h
#define Semi_properties_h

char state = 's';

// Inlet Process Fluid Properties From Aspen HYSYS Semi Simulation
const float Tp_in = 287.0; // F
const float p_in = 236; // psia
//float m = 3.182e+5; // lb/hr
//float m = 3.182e+5 / 3; // lb/hr
const float Cp_in = 0.7112; // Btu/lb/F
const float rho_in = 0.5445; // lb/ft^3
const float mu_l_in = 0.4962; // lb/ft/h : from HYSYS 0.2051 cP
const float mu_v_in = 0.039359; // lb/ft/h : from HYSYS 1.627e-02 cP
const float vf_l_in = 1.163e-004; // dimless
const float vf_v_in = 0.9999; // dimless
const float H_in = -438.1; // Btu/lb from HYSYS

// Outlet Process Fluid Properties
const float Tp_out = 104.0; // F
const float Cp_out = 0.6732; // Btu/lb/F
const float rho_out = 0.7376; // lb/ft^3
const float mu_l_out = 0.8999; // lb/ft/h : from HYSYS 0.3720 cP
const float mu_v_out = 0.026343; // lb/ft/h : from HYSYS 1.089e-002 cP
const float vf_l_out = 8.312e-003; // dimless
const float vf_v_out = 0.9917; // dimless
const float H_out = -639.1; // Btu/lb from HYSYS

// Average Process Fluid Properties
//const float Cp = (Cp_in + Cp_out) / 2;
//const float rho = (rho_in + rho_out) / 2;
//const float mu_in = vf_l_in*mu_l_in + vf_v_in*mu_v_in;
//const float mu_out = vf_l_out*mu_l_out + vf_v_out*mu_v_out;
//const float mu = (mu_in + mu_out) / 2;
//const float k = 0.067786; // Btu/hr/ft/F

// Process Fluid Properties At Mean Temperature
const float Cp = 0.6983; // Btu/lb/F
const float rho = 0.6069; // lb/ft^3
const float mu = 0.128086; // lb/ft/hr
const float k = 0.067786; // Btu/hr/ft/F

#endif /* Semi_properties_h */


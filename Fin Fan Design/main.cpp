//
//  main.cpp
//  Fin Fan Design
//
//  Created by Jay Tunkie Saunders on 2/23/18.
//  Copyright © 2018 Jay Tunkie Saunders. All rights reserved.
//
//  Design Methodology Sourced From: "Process Heat Transfer: Principles, Applications, and Rules of Thumb"
//
//  Data Notes:
//  - Tube conductivity
//  - Fouling factor
//  - Air-side pressure needs to account for enclosures, etc ...
//

#include <iostream>
#include <cmath>
#include <math.h>
#include <vector>
#include <tuple>
#include <iomanip>
//#include "semiProperties.h" // semi-regenerative process properties
#include "contProperties.h" // continuous regenerative proces properties


// Air Properties
const float Ta_in = 95.2; // F
const float rho_std = 0.075; // lb/ft^3 : Standard air conditions
const float rho_avg = 0.05993; // lb/ft^3 : Average air condtions for Denver, Colorado
const float Cp_air = 0.241; // Btu/lb/F
const float mu_air = 0.05191; // lb/ft/h : 200 F at standard
const float k_air = 0.018; // Btu/h/ft/F : thermal conductivity of air at 200 F and 1 bar
const float Pr_air = Cp_air * mu_air / k; // dimless

// Miscellaneous Properties
const float k_al = 137.5; // Btu/h/ft/F : thermal conductivity of aluminum

// Calculate Heat Transferered
float energyBalance(float m)
{
    float q;
    q = m * abs(H_out - H_in);
    return q;
}

// Calculate Mass Flow Rate Of Air
float airFlow(float q, float Ta_out)
{
    float ma;
    ma = q / (Cp_air * abs(Ta_out - Ta_in));
    return ma;
}

// Calculate LMTD
float lmtd(float F, float Ta_out)
{
    float delT_p, delT_a, LMTD;
    delT_p = abs(Tp_out - Tp_in);
    delT_a = abs(Ta_out - Ta_in);
    LMTD = (delT_p - delT_a) / log(delT_p / delT_a);
    LMTD = F * LMTD;
    return LMTD;
}

// Calculate Heat Transfer Area
float area(float q, float U_d, float F, float LMTD)
{
    float A;
    A = q / (U_d * F * LMTD);
    return A;
}

// Calculate Tubing
std::tuple<int, int, float, float, float, float> tubular(float ma, float V_face, float A, int config)
{
    float A_face, ratio, table_ratio;
    A_face = ma / (rho_std * V_face);
    ratio = A / A_face;
    
    std::vector<float> ratios = {0, 0, 0, 0};
    std::vector<float> distances = {0, 0, 0, 0};
    
    // External surface area per unit bundle face areas taken from Table 12.1.
    // Each corresponds to a different standard configuration
    if (config == 0)
    {
        ratios[0] = 60.6;
        ratios[1] = 80.8;
        ratios[2] = 101.0;
        ratios[3] = 121.2;
    }
    else if (config == 1)
    {
        ratios[0] = 80.4;
        ratios[1] = 107.2;
        ratios[2] = 134.0;
        ratios[3] = 160.8;
    }
    
    
    for (int i = 0; i < ratios.size(); i++)
    {
        distances[i] = abs(ratio - ratios[i]);
    }
    
    float small;
    small = distances[0];
    for (int i = 1; i < distances.size(); i++)
    {
        if (distances[i] < small)
        {
            small = distances[i];
        }
    }
    
    int idx = -1, count = 0;
    for (int i = 0; i < distances.size(); i++)
    {
        if (distances[i] == small)
        {
            idx = i;
            count++;
            
        }
        else if (count == 0)
        {
            std::cout << "Error: Failed to find a matching table ratio!" << std::endl;
        }
    }
    
    int n_tube_rows;
    float W, L, n_tubes_tot, AL = 0, pitch = 0, V_face_std;
    n_tube_rows = idx + 4;
    table_ratio = ratios[idx];
    A_face = A / table_ratio;
    W = sqrt(A_face / 3);
    L = 3 * W;
    
    if (config == 0)
    {
        AL = 3.90; // Table 12.1
        pitch = 2.25;
    }
    else if (config == 1)
    {
        AL = 5.58; // Table 12.1
        pitch = 2.50;
    }
    
    n_tubes_tot = A / (AL * L);
    
    int n_t = (int)(roundf(n_tubes_tot / 4) * 4 + 0.5); // Round to nearest integer divisible by 4
    int n_t_row = n_t / 4;
    W = (pitch*n_t_row + 2) / 12;
    A_face = W * L;
    V_face_std = (ma / 60) / (rho_std * A_face);
    
    std::tuple<int, int, float, float, float> tube_output;
    {
        return std::make_tuple(n_t, n_t_row, W, L, A_face, V_face_std);
    }
}

// Calculate Tube Passes
float tubePasses(int n_t, float m)
{
    float V;
    V = ((m / 3600) * (1.0 / n_t)) / (rho * M_PI * pow((0.81 / 12), 2) / 4);
    return V;
}

// Calculate Reynolds Number
float reynolds(int n_p, int n_t, float m)
{
    float Re;
    Re = (4 * m * ((float)n_p / n_t)) / (M_PI * (0.81 / 12) * mu);
    return Re;
}

// Calculate Tube-side Pressure Drop
float tubeDrop(float Re, float n_p, float n_t, float L, float m)
{
    float f, G, s, delP, phi;
    f = 0.4137 * pow(Re,-0.2585);
    G = m * (n_p / n_t) / (M_PI * pow((0.81 / 12),2) / 4);
    s = rho / 62.43;
    phi = 1;
    delP = f * n_p * L * pow(G,2) / (7.50e+12 * (0.81 / 12) * s * phi);
    delP += ((1.334e-13) * (2 * n_p - 1.5) * (pow(G, 2) / s)); // Nozzle pressure drop is causing issues
    
    float G_n, Re_n;
    G_n = m / 0.1390; // Assuming 5 in, schedule 40 nozzles -> source for 0.1390 ft^2
    Re_n = ((5.047 / 12)*G_n) / mu;
    delP += (2.0e-13*pow(G_n, 2)) / s;
    
    return delP;
}

// Calculate Required Overall Coefficient
float requiredU(float q, float A, float F, float LMTD)
{
    float U_r = q / (A * F * LMTD);
    return U_r;
}

// Calculate Tube-side Convective Heat Transfer Coefficient
float hTube(float k, float Re, float Pr)
{
    float phi = 1;
    float h = (k / 0.0675) * 0.023 * pow(Re, 0.8) * pow(Pr, 1/3) * pow(phi, 0.14);
    return h;
}

// Calculate Air-side Convective Heat Transfer Coefficient
std::tuple<float, float, float> hAir(float V_face_std, float pitch, float D_r, float n_f, float b, float AtAo)
{
    float tau = 0.013; // in : from the SOURCE
    float V_face_avg = V_face_std * (rho_std / rho_avg);
    float V_max = (pitch * V_face_avg / (pitch - D_r - 2 * n_f * b * tau)) * 60;
    float Re = (D_r / 12) * V_max * rho_avg / mu_air;
    float Nu = 0.38 * pow(Re, 0.6) * pow(Pr_air, 1/3) * pow(AtAo, -0.15);
    float h_o = (k_air / (D_r / 12)) * Nu;
    
    std::tuple<float, float, float> air_output;
    {
        return std::make_tuple(h_o, Re, V_max);
    }
}

// Calculate Fin Efficiency
float finEff(float D_r, float b, float h_o, float n_f)
{
    float tau = 0.013; // fin thickness in inches
    float r_1 = D_r / 2;
    float r_2 = r_1 + b;
    float r_2c = r_2 + tau / 2;
    float psi = ((r_2c - r_1) * (1 + 0.35 * log(r_2c / r_1))) / 12.0;
    float m = pow(2 * h_o / (k_al * (tau / 12)), 0.5);
    float eff = tanh(m * psi) / (m * psi);
    float A_fins = 2 * n_f * M_PI * (pow(r_2c, 2) - pow(r_1, 2)); // in^2
    float A_prime = 2 * M_PI * r_1 * (1.0 - n_f * tau);
    float AfAt = A_fins / (A_fins + A_prime);
    float ApAt = 1 - AfAt;
    float eff_w = ApAt + eff * AfAt;
    return eff_w;
}

// Calculate Clean Overall Coefficient
float cleanU(float AtAi, float h_i, float AtAl, float D_r, float D_i, float eff_w, float h_o)
{
    float k_tube = k_al; // INVESTIGATE
    float U_c = pow((AtAi / h_i) + (AtAl * log(D_r / D_i) / (2 * M_PI * k_tube)) + (1 / (eff_w * h_o)), -1);
    return U_c;
}

// Calculate Fouling Allowance
float fouling(float R_di, float AtAi, float eff_w)
{
    float R_do = 0; // Air-side fouling taken to be zero
    float R_d = R_di * AtAi + R_do / eff_w;
    return R_d;
}

// Calculate Overall Design Coefficient
float designU(float U_C, float R_D)
{
    float U_D = pow(1 / U_C + R_D, -1);
    return U_D;
}

// Calculate Air-side Pressure Drop
float airDrop(float D_r, float b, float pitch, float n_f, float Re_air, float V_max, float rows)
{
    float D_f = D_r + 2 * b;
    float a = (pitch - D_f) / D_r;
    float tau = 0.013; // Suggested by SOURCE
    float l = 1 / n_f - tau;
    float Re_eff = Re_air * (l / b);
    float f = (1 + 2 * exp(-a / 4) / (1 + a)) * (0.021 + 27.2 / Re_eff + 0.29 / pow(Re_eff, 0.2));
    float G = rho_avg * V_max;
    float delP = 9.22e-10 * f * rows * pow(G, 2) / rho_avg;
    delP = delP * 1.30; // Account for losses from hail guards, etc ...
    return delP;
}

// Calculate Fan Sizing
std::tuple<float, float, float> fanSize(float A_face, float fbp, float delP, float ma)
{
    float D_fan_min = sqrt(0.4 * A_face * 4 / (fbp * M_PI));
    float FSP = delP; // Assumption from SOURCE
    float v_fan = ((1 / fbp) * (ma / 60)) / rho_avg;
    
    std::tuple<float, float, float> fan_output;
    {
        return std::make_tuple(D_fan_min, FSP, v_fan);
    }
}

// Calculate Motor Size
std::tuple<float, float, float> motorSize(float FSP, float v_fan, float D_fr, float alpha_fr, float rho_fr, float eff_fan)
{
    float g_c = 32.174; // acceleration constant
    float V_fr = (v_fan / (M_PI * pow(D_fr, 2) / 4)) / 60;
    float delP_fan = FSP + alpha_fr * rho_fr * pow(V_fr, 2) / (2 * g_c);
    
    float W_fan = delP_fan * v_fan / (6342.0 * eff_fan);
    float eff_sr = 0.95;
    float W_motor = W_fan / eff_sr;
    
    std::tuple<float, float, float> motor_output;
    {
        return std::make_tuple(delP_fan, W_fan, W_motor);
    }
}

// Summarize Design Results - WORK IN PROGRESS
void summary()
{
    float holder = 0; // place holder until function is complete
    std::cout << "Design Summary: " << std::endl;
    std::cout << std::endl;
    std::cout << "Number of fan bays: " << holder << std::endl;
    std::cout << "Number of tube bundles per bay: " << holder << std::endl;
    std::cout << "Number of fan per bays: " << holder << std::endl;
    std::cout << "Bundle width and length: " << holder << std::endl;
    std::cout << "Number of tube rows: " << holder << std::endl;
    std::cout << "Number of tube passes: " << holder << std::endl;
    std::cout << "Number of tubes: " << holder << std::endl;
    std::cout << "Number of fan bays: " << holder << std::endl;
    std::cout << "Tubing type: " << holder << std::endl;
    std::cout << "Tube size: " << holder << std::endl;
    std::cout << "Tube layout: " << holder << std::endl;
    std::cout << "Fins: " << holder << std::endl;
    std::cout << "Heat transfer surface area: " << holder << std::endl;
    std::cout << "Draft type: " << holder << std::endl;
    std::cout << "Fan diameter: " << holder << std::endl;
    std::cout << "Motor size: " << holder << std::endl;
    std::cout << "Tube-side nozzles: " << holder << std::endl;
    std::cout << "Headers: " << holder << std::endl;
    std::cout << "Materials: " << holder << std::endl;
}

int main()
{
    // Indicates whether properties corresponding to the semi-regen or continous-regen system are being utilized
    float m = 0; // initialize process fluid mass flow rate
    if (state == 's')
    {
        std::cout << "Preliminary Design Results For The Semi Process:" << std::endl << std::endl;
        m = 3.182e+5; // lb/hr
    }
    else if (state == 'c')
    {
        std::cout << "Preliminary Design Results For The Continuous Process:" << std::endl << std::endl;
        m = 2.249e+5; // lb/hr
    }
    else
    {
        std::cout << "Error";
    }
    
    // set output style
    std::cout << std::fixed;
    std::cout << std::setprecision(2);
    
    // Manipulate process fluid mass flow rate to achieve thermally, dimesionally, physicall feasible design ...
    int bundles = 3;
    m = m / bundles;
    std::cout << "Number of tube bundles: " << bundles << std::endl;
    const float q = energyBalance(m); // Btu/hr
    std::cout << "Quantity of heat removed from the process fluid: " << q << " Btu/hr" << std::endl;
    
    float ma, Ta_out;
    Ta_out = 200; // F : MANIPULATE / NOTE: was previously 200 F for calcs
    std::cout << "Air outlet temperature: " << Ta_out << " F" << std::endl;
    ma = airFlow(q, Ta_out); // lb/hr
    std::cout << "Mass flow rate of air required: " << ma << " lb/hr" << std::endl;
    
    float F, LMTD;
    F = 0.9; // dimless : Correction factor guess
    LMTD = lmtd(F, Ta_out); // F
    std::cout << "LMTD with correction factor: " << LMTD << " F" << std::endl;
    
    float U_d, A;
    U_d = (3.3 + 4.7) / 2; // Btu/h/ft^2/F : From Table X
    A = area(q, U_d, F, LMTD); // ft^2
    std::cout << "Heat transfer area: " << A << " ft^2" << std::endl;
    
    int config = 0; // set to 0 for first configuration or 1 for second configuration of Table 12.1
    std::cout << "Configuration: " << config << std::endl;
    float AtAl = 0, AtAi = 0, pitch = 0, b = 0, n_f = 0, D_r = 1.0, D_i = 0.81, AtAo = 0;
    
    if (config == 0)
    {
        AtAl = 3.90; // Table 12.1
        AtAi = 17.9; // Table 12.1 for 13 BWG
        pitch = 2.25; // in
        b = 0.5; // in : fin height
        n_f = 9; // fins per inch
        AtAo = 14.5; // Ratio of total area to outer area
    }
    else if (config == 1)
    {
        AtAl = 5.58; // Table 12.1
        AtAi = 26.3; // Table 12.1 for 13 BWG
        pitch = 2.50;
        b = 0.625; // in : fin height
        n_f = 10; // fins per inch
        AtAo = 21.4; // Ratio of total area to outer area
    }
    
    float V_face = 600.0; // ft/min : Chosen based on Example 12.1 INVESTIGATE
    std::tuple<int, int, float, float, float, float> tube_output = tubular(ma, V_face, A, config);
    float n_t, n_t_row, W, L, A_act, V_face_std;
    n_t = std::get<0>(tube_output);
    n_t_row = std::get<1>(tube_output);
    W = std::get<2>(tube_output);
    L = std::get<3>(tube_output);
    A_act = std::get<4>(tube_output);
    V_face_std = std::get<5>(tube_output);
    std::cout << "Number of tubes: " << n_t << std::endl;
    std::cout << "Number of tube per row: " << n_t_row << std::endl;
    std::cout << "Width of bundle: " << W << " ft" << std::endl;
    std::cout << "Length of bundle: " << L << " ft" << std::endl;
    std::cout << "Actual bundle face area: " << A_act << " ft^2" << std::endl;
    std::cout << "Actual face velocity (std): " << V_face_std << " ft/min" << std::endl;
    
    std::cout << std::endl << "Output for Rigorous Design Process:" << std::endl << std::endl;
    
    float V;
    V = tubePasses(n_t, m);
    for (int i = 1; i < 5; i++)
    {
        if (i == 1)
        {
            std::cout << "Tube-side velocity with 1 pass: " << V << " ft/s" << " with corresponding Reynolds number: " << reynolds(i, n_t, m) << " and pressure drop: " << tubeDrop(reynolds(i, n_t, m), i, n_t, L, m) << " psi" << std::endl;
        }
        else
        {
            std::cout << "Tube-side velocity with " << i << " passes: " << i*V << " ft/s" << " with corresponding Reynolds number: " << reynolds(i, n_t, m) << " and pressure drop: " << tubeDrop(reynolds(i, n_t, m), i, n_t, L, m) << " psi" << std::endl;
        }
    }
    
    float n_p = 1;
    F = 1; // correction factor for one-pass
    float U_r = requiredU(q, A, F, LMTD); // Btu/h/ft^2/F
    std::cout << "Required overall heat transfer coefficient: " << U_r << " Btu/h/ft^2/F" << std::endl;
    
    float Re = reynolds(n_p, n_t, m);
    float Pr = Cp * mu / k;
    float h_i = hTube(k, Re, Pr);
    std::cout << "Tube-side convective heat transfer coefficient: " << h_i << " Btu/h/ft^2/F" << std::endl;
    
    //float h_o = hAir(V_face_std, pitch, D_r, n_f, b, AtAo);
    std::tuple<float, float, float> air_output = hAir(V_face_std, pitch, D_r, n_f, b, AtAo);
    float h_o = std::get<0>(air_output);
    float Re_air = std::get<1>(air_output);
    float V_max = std::get<2>(air_output);
    std::cout << "Air-side convective heat transfer coefficient: " << h_o << " Btu/h/ft^2/F" << std::endl;
    
    float eff_w = finEff(D_r, b, h_o, n_f);
    std::cout << "Weighted fin efficiency: " << eff_w << std::endl;
    
    float U_c = cleanU(AtAi, h_i, AtAl, D_r, D_i, eff_w, h_o);
    std::cout << "Clean overall heat transfer coefficient: " << U_c << " Btu/h/ft^2/F" << std::endl;
    
    std::string ok;
    if (U_c > U_r)
    {
        ok = "Solution is valid -- Continue.";
    }
    else
    {
        ok = "Solution is not valid -- Do not continue!";
    }
    
    std::cout << std::endl << "STATUS: " << ok << std::endl << std::endl;
    
    float R_di = 0.001; // Btu/h/ft^2/F : Fouling given by algo NEED TO FIND VALUE FOR OUR SYSTEM
    float R_d = fouling(R_di, AtAi, eff_w);
    U_d = designU(U_c, R_d);
    std::cout << "Design overall heat transfer coefficient: " << U_d << " Btu/h/ft^2/F" << std::endl;
    
    float over_surf = (U_c / U_r - 1) * 100;
    float over_des = (U_d / U_r - 1) * 100;
    std::cout << "Over-surface: " << over_surf << "%" << std::endl;
    std::cout << "Over-design: " << over_des << "%" << std::endl;
    if (over_surf >= 0 && over_des >= 0)
    {
        std::cout << std::endl << "FEASIBILITY: " << "The design is thermally feasible!" << std::endl << std::endl;
    }
    else
    {
        std::cout << std::endl << "FEASIBILITY: " << "The design is not thermally feasible." << std::endl << std::endl;
    }
    
    float delP_tube = tubeDrop(Re, n_p, n_t, L, m);
    std::cout << "Tube-side pressure drop: " << delP_tube << " psi" << std::endl;
    
    float rows = n_t / n_t_row;
    float delP_air = airDrop(D_r, b, pitch, n_f, Re_air, V_max, rows);
    std::cout << "Air-side pressure drop: " << delP_air << " in. H2O" << std::endl;
    float fbp = 2; // Numbers of fans per bay
    std::tuple<float, float, float> fan_output = fanSize(A_act, fbp, delP_air, ma);
    float D_fan_min = std::get<0>(fan_output);
    float FSP = std::get<1>(fan_output);
    float v_fan = std::get<2>(fan_output);
    std::cout << "Minimum fan diameter: " << D_fan_min << " ft" << std::endl;
    std::cout << "Fan static pressure: " << FSP << " in. H2O" << std::endl;
    std::cout << "Volumetric flow rate per fan: " << v_fan << " acfm" << std::endl;
    
    std::cout << std::endl;
    
    float D_fr = 12.0;
    float alpha_fr = 1.08;
    float rho_fr = 0.04911;
    float eff_fan = 0.65;
    std::tuple<float, float, float> motor_output = motorSize(FSP, v_fan, D_fr, alpha_fr, rho_fr, eff_fan);
    float delP_fan = std::get<0>(motor_output);
    float W_fan = std::get<1>(motor_output);
    float W_motor = std::get<2>(motor_output);
    std::cout << "Total pressure difference across fan: " << delP_fan << " in. H2O" << std::endl;
    std::cout << "Brake power: " << W_fan << " hp" << std::endl;
    std::cout << "Motor size: " << W_motor << " hp" << std::endl;
    
    
    std::cout << std::endl;
    
    //summary(bundles); // Function to summarize final design results - WORK IN PROGRESS
    
    std::cout << std::endl;
    return 0;
}


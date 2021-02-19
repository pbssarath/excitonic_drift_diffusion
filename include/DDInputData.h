/*
  Copyright 2016-2017 Baskar Ganapathysubramanian

  This file is part of TALYFem.

  TALYFem is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 2.1 of the
  License, or (at your option) any later version.

  TALYFem is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TALYFem.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef DD_INPUT_DATA_HPP
#define DD_INPUT_DATA_HPP

#include <string>
#include <talyfem/talyfem.h>
struct DDInputData:public InputData {
    std::string output_extension;
    int ifPrintPltFiles;          ///< whether to print .plt files at start and end
    int shouldFail;
    int use_elemental_assembly_;  ///< 1 if elemental assembly is to be used

    // Universal constants
    double epsilon_0 = 8.854e-12; //< Permittivity of free space      --As m^-1 V^-1
    double k_B = 1.3806503e-23;   //< Boltzmann constant        --m^2 kg s^(-2) K^(-1)
    double T = 300;               //< Room temperature        --K
    double q = 1.60217646e-19;    //< Elementary charge       --As,coulombs
    double Vt = k_B * T / q;
    // double Lx = 1, Ly = 1;
    // Parameters, in M.K.S

    // Material Parameters
    double N_C = 2.5e25;          //< Electron effective density of states  --m^(-3)
    double N_V = 2.5e25;          //< Hole effective density of states    --m^(-3)
    double E_g = 1.1;             //< Band gap, (E_V - E_C)     --eV
    double mu_n = 1e-7;           //< Electron mobility when in majority          --m^2/(Vs)
    double mu_p = 1e-7;           //< Hole mobility when in majority              --m^2/(Vs)
    double G = 1e28;                     //< Generation rate         --m^(-3)/s
    double a = 1.8e-9;            //< e/h pair distance (most probable)   --m
    double Height;                //< Domain thickness, with additional materials if applicable         --m
    double realHeight;            //< real active layer thickness without adding additional material;
    double Geo2Phy;               //< factor to convert geometry distance contour to physical contour (in nanometer)
    double eps_D = 3 * 8.854e-12;
    double eps_A = 3 * 8.854e-12;
    double threshold = 0.5, transportLayer = 0.;
    double heightofTransportLayer;  //hight of added material layer;
    double gridTkness;
    double interfaceThk = 2.;
    bool isPeriodic = false, threephase = false;  //flags for periodic boundary and three phase model;
    double D_fraction = 0.5;     //if three phase model, then this is the donor volume fraction in the mixture phase.
    double muRatio = 0.01;
    double potentialFactor = 0.;
    double lamdaFactor = 3;
    int ngp = 1;

    // Structured Mesh parameters
    int elmnx = 0;                //< Number of elements in x-direction
    int elmny = 0;                //< Number of elements in y-direction
    int elmnz = 0;                //< Number of elements in z-direction
    double beta_1;                //< Clustering parameter

    // Calculated parameters for 1-D case
    double Voc;
    double vapp_step = 0.02;
    double vapp_ratio = 5;
    int new_bc = 0;

    // other parameters
    char **args;
    int couple = 3;
    // Recombination model number being used
    int rmodelno;

    double alpha0;                //< absorption coefficient      --α0, m^(-1)
    double Gamma0;                //< photon flux         --Γ0, m^(-2) s^(-1)
    double tau_x;                 //< average lifetime of an excitons   --τX, s
    double mu_x = 1.e-9;          //< exciton mobility        --µX, m^2 V s

    // Morphology parameter
    double midphi;
    double chi;

    // Input ms file name
    std::string MSfile;

    // Input mesh file if _readmesh_ is defined
    char *meshFile;
    // Input ms file [tecplot] node number in x-dir and y-dir
    int msDataNx, msDataNy, msDataNz;

    // Exciton generation rate profile
    // 1-constant 2-exponential 3-reflecting cathode
    int Gprofile;
    //----------------------------------------------------------------------------------------------

    double insulatingdomain;      //the size of insulating domain around high-k particles.


    double V_t() {
        return Vt;
    }                             // Thermal voltage [V] [V]=0.0258520269

    double n_int(double ms, double DT) {
        double returnVal;
        returnVal = sqrt(this->N_C * this->N_V) * exp(-this->E_g / (2 * V_t()));
        return returnVal;
    }

    // Scaling parameters [ref: S. Selberherr, 1984, page 141]
    double x0() {
        return this->Height;
    }

    double V0() {
        return this->Vt;
    }

    double C0() {
        return this->N_C;
    }

    double mu0() {
        return (this->mu_n > this->mu_p ? this->mu_n : this->mu_p);
    }

    double U0() {
        return mu0() * V0() * C0() / (x0() * x0());
    }

    // Perturbation parameter
    double lambda2_(double epsilong) {
        double eps;
        eps = epsilong; //*lamdaFactor;

        double returnVal = V0() * eps / (x0() * x0() * this->q * C0());
        return returnVal;
    }                             // Perturbation parameter

    // Scaled variables
    //----------------------------------------------------------------------------------------------
    double x_hat(double x) {
        return x / x0();
    }

    double V_hat(double V) {
        return V / V0();
    }

    double c_hat(double c) {
        return c / C0();
    }

    double mu_hat(double mu) {
        return mu / mu0();
    }

    double U_hat(double U) {
        return U / U0();
    }
    //----------------------------------------------------------------------------------------------

    // Unscaled variables
    //----------------------------------------------------------------------------------------------
    double x_unhat(double x_hat) {
        return x_hat * x0();
    }

    double V_unhat(double V_hat) {
        return V_hat * V0();
    }

    double c_unhat(double c_hat) {
        return c_hat * C0();
    }

    double mu_unhat(double mu_hat) {
        return mu_hat * mu0();
    }

    double U_unhat(double U_hat) {
        return U_hat * U0();
    }

    //----------------------------------------------------------------------------------------------
    // Intrinsic voltage [V]
    double Vbi(double ms) {
        return this->E_g;
    }


    DDInputData() : InputData(), ifPrintPltFiles(1), shouldFail(0), use_elemental_assembly_(0) {
    }

    bool ReadFromFile(const std::string &filename = std::string("config.txt")) {
        // Read config file and initialize basic fields
        InputData::ReadFromFile(filename);
        if (ReadValue("outputExtension", output_extension)) { }
        if (ReadValue("D_fraction", D_fraction)) { }
        if (ReadValue("threephase", threephase)) { }
        if (ReadValue("ifPrintPltFiles", ifPrintPltFiles)) { }
        if (ReadValue("shouldFail", shouldFail)) { }
        if (ReadValue("use_elemental_assembly", use_elemental_assembly_)) { }
        if (ReadValue("muRatio", muRatio)) {}
        if (ReadValue("mu_n", mu_n)) { }
        if (ReadValue("mu_p", mu_p)) { }
        if (ReadValue("mu_x", mu_x)) { }
        if (ReadValue("E_g", E_g)) { }
        if (ReadValue("potentialFactor", potentialFactor)) {}
        if (ReadValue("elmnx", elmnx)) { }
        if (ReadValue("Weak_bc", new_bc)) { }
        if (ReadValue("transportLayer", transportLayer)) { }

        if (ReadValue("elmny", elmny)) {
            L[0] = double(elmnx) / elmny / (1. + 2. * transportLayer);
            L[1] = 1;
        }
        if (ReadValue("elmnz", elmnz)) {
            L[2] = double(elmnz) / elmny / (1. + 0. * transportLayer);
        }

        if (ReadValue("threshold", threshold)) { }
        if (ReadValue("Height", Height)) {
            Geo2Phy = Height * 1.e9 / (this->elmny);
            realHeight = Height / (1. + 2. * transportLayer);
            heightofTransportLayer = (Height - realHeight) / 2. * 1.e9;
        }

        if (ReadValue("eps_A", eps_A)) {
            eps_A *= 8.854e-12;
        }

        if (ReadValue("eps_D", eps_D)) {
            eps_D *= 8.854e-12;
        }

        if (ReadValue("isPeriodic", isPeriodic)) { }
        if (ReadValue("Voc", Voc)) { }
        if (ReadValue("vapp_step", vapp_step)) { }
        if (ReadValue("vapp_ratio", vapp_ratio)) { }
        if (ReadValue("Gx", G)) { }
        if (ReadValue("tau_x", tau_x)) { }
        if (ReadValue("ngp", ngp)) { }
        if (ReadValue("a", a)) { }
        if (ReadValue("beta", beta_1)) { }
        if (ReadValue("MSfile", MSfile)) { }
        if (ReadValue("interfaceThk", interfaceThk)) { }
        if (ReadValue("T", T)) {
            //then update thermal voltage with new temperature.
            this->Vt = k_B * T / q;
        }
        return true;
    }

    bool CheckInputData() const {
        if (((typeOfIC == 0) && (inputFilenameGridField.empty())) ||
            ((typeOfIC != 0) && (!inputFilenameGridField.empty()))) {
            PrintWarning("IC not set properly check!", typeOfIC, " ", inputFilenameGridField);
            return false;
        }

        return InputData::CheckInputData();
    }

    std::ostream &print(std::ostream &oss) const {
        PrintLogStream(oss, "ifPrintPltFiles = ", ifPrintPltFiles);
        PrintLogStream(oss, "shouldFail = ", shouldFail);
        return InputData::print(oss);
    }

};

#endif

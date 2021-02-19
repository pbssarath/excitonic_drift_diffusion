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
// --- end license text --- //
#ifndef INCLUDE_DDEXCITON_H_
#define INCLUDE_DDEXCITON_H_

#include "DDNodeData.h"
#include "initial_guess.h"
#include "gen_diss.hh"

/**
 * This class implements exciton diffusion problem.
 *
 * The class can handle 1, 2, or 3 dimensions. Boundary conditions are
 * Direchlet (fixed value on the boundaries) with the boundary values set to be
 * equal to the experimental value at the boundaries.
 *
 */
class DDExciton:public CEquation < DDNodeData > {
public:
    // indices to timing arrays. These are locations in the timers_ array that
    // correspond to specific stages of the code that we wish to time.
    static const int kTimerSolve = 0; ///< index to timer for solve process
    static const int kTimerAssemble = 1;  ///< index to timer for assemble
    static const int kTimerKSPSolve = 2;  ///< index to timer for KSPsolve
    static const int kTimerUpdate = 3;  ///< index to timer for update process
    gen_diss *gendiss;
    float cont = 1;

    /**
     * Constucts the solver by setting up timers.
     */
    explicit DDExciton(DDInputData *input_data, AssemblyMethod assembly_method = kAssembleGaussPoints)
            : CEquation<DDNodeData>(assembly_method),
              input_data_
                      (input_data) {
        timers_[kTimerSolve].set_label("Solve");
        timers_[kTimerAssemble].set_label("Assemble");
        timers_[kTimerKSPSolve].set_label("KSPSolve");
        timers_[kTimerUpdate].set_label("Update");
    }

    /**
     * Destroys the object
     *
     * This also prints the average times required for each of the processes that
     * we times during execution.
     */
    virtual ~ DDExciton() {
        timers_[kTimerSolve].PrintGlobalAverageSeconds();
        timers_[kTimerAssemble].PrintGlobalAverageSeconds();
        timers_[kTimerKSPSolve].PrintGlobalAverageSeconds();
        timers_[kTimerUpdate].PrintGlobalAverageSeconds();
    }

    /**
     * Sets up the essential (Direchlet) boundary conditions
     *
     * The values of the nodes on each boundary are set to experimental values.
     */
    virtual void fillEssBC() {
        // Setting boundary conditions. for 2D problem
        this->initEssBC();
        for (int nodeID = 0; nodeID < this->p_grid_->n_nodes(); nodeID++) {
            /**
             * boundary IDs (second argument to BoNode):
                4
              _ _ _
             |     |
            1|     |2
             |_ _ _|
                3
            */
            if (this->p_grid_->BoNode(nodeID, 3)) {
                for (int i = 0; i < 1; i++) {
                    this->specifyValue(nodeID, i, 0);
                }

                this->p_data_->GetNodeData(nodeID).x[0] = 0;
            } else if (this->p_grid_->BoNode(nodeID, 4)) {
                for (int i = 0; i < 1; i++) {
                    this->specifyValue(nodeID, i, 0);
                }
                this->p_data_->GetNodeData(nodeID).x[0] = 0;
            }
        }
    }

    /**
     * Solves the system for a given time value.
     *
     * @param t the current time value
     * @param dt the time step between solve times
     */
    virtual void Solve(double dt = 0.0, double t = 0.0) {
        timers_[kTimerSolve].Start();  // we want to time the entire solve process
        this->t_ = t;               // store these in the object data for later use
        this->dt_ = dt;

        PrintInfo("  Filling essential boundary conditions...");
        fillEssBC();               // Set the boundary conditions
        PrintInfo("  Applying Ess BC to solution...");
        ApplyEssBCToSolution();    // apply boundary conditions to the solution vector
        timers_[kTimerAssemble].Start(); // we're timing just the assembly
        PrintInfo("  Assembling...");
        Assemble(false);           // assemble the system
        timers_[kTimerAssemble].Stop();
        PrintInfo("  Applying Ess BC...");
        ApplyEssBC();              // apply boundary conditions to the assembled system

        timers_[kTimerKSPSolve].Start(); // time the KSPSolve step
        PrintInfo("  Solving with KSP...");
        SolveKSP(solution_, 1, 0); // run the KSP solver (uses PETSc)
        timers_[kTimerKSPSolve].Stop();
        if (input_data_->isPeriodic) updatePeriodicSol();
        timers_[kTimerUpdate].Start(); // time the process of saving the solution
        // The result of the system solve is in the solution_ vector. We want to
        // store this data in our node data arrays where we keep track of the
        // current and previous heat values. This function call will put the data
        // from solution_ into location 0 of each NodeData object.
        PrintInfo("  Copying data from solution vector into node data...");
        for (int i = 0; i < p_grid_->n_nodes(); i++) {
            //double new_value = solution_(i*n_dof());  // the value for the ith node
            p_data_->GetNodeData(i).value(X_ID) = solution_(i * n_dof()); // store in node data object
            // p_data_->GetNodeData(i).value (25) = solution_(i*n_dof() + 1);  // store in node data object
            // p_data_->GetNodeData(i).value (26) = solution_(i*n_dof() + 2);  // store in node data object
            // p_data_->GetNodeData(i).value (27) = solution_(i*n_dof() + 3);  // store in node data object
            if (p_data_->GetNodeData(i).value(X_ID) < 0)
                p_data_->GetNodeData(i).value(X_ID) *= -1;  // store in node data object;  // store in node data object
            //p_data_->GetNodeData(i).value (25) = solution_(i*n_dof()+1);  // store in node data object
            //p_data_->GetNodeData(i).value (26) = solution_(i*n_dof()+2);  // store in node data object
            //p_data_->GetNodeData(i).value (27) = solution_(i*n_dof()+3);  // store in node data object
        }

        PrintInfo("  Done!");

        timers_[kTimerUpdate].Stop();
        timers_[kTimerSolve].Stop();
    }

    /**
     * Fills the Ae and be structures with data for a single Gauss point.
     *
     * @param fe the element we are assembling a Gauss point from
     * @param Ae the element matrix to put data in
     * @param be the element vector to put data in
     */
    virtual void Integrands(const FEMElm &fe, ZeroMatrix<double> &Ae, ZEROARRAY<double> &be) {
        FEMElm fe_linear(p_grid_, BASIS_ALL);
        //   if (fe.basis_function() == BASIS_HERMITE) {
        fe_linear.refill(fe.elem(), BASIS_LINEAR, 0);
        fe_linear.calc_at(fe.itg_pt());
        //  } else {
        //   fe_linear = fe;
        // }

        const int nsd = fe.nsd(); // # of dimensions: 1D, 2D, or 3D
        const int n_basis_functions = fe.nbf();  // # of basis functions
        const double detJxW = fe.detJxW(); // (determinant of J) cross W = this->p_data_->valueFEM (fe_linear,D2PHIBAR_ID);
        double Dphi = this->p_data_->valueFEM(fe_linear, PHIBAR_ID);
        double dDTdx = this->p_data_->valueFEM(fe_linear, 32);
        double dDTdy = this->p_data_->valueFEM(fe_linear, 33);
        double dDTdz = 0;
        double absdDT;
        if (nsd == 2) {
            absdDT = sqrt(dDTdx * dDTdx + dDTdy * dDTdy) + 1e-16;
        } else {
            absdDT = sqrt(dDTdx * dDTdx + dDTdy * dDTdy + dDTdz * dDTdz) + 1e-16;
        }
        double phi_x = this->p_data_->valueDerivativeFEM(fe_linear, PHI_ID, 0);
        double phi_y = this->p_data_->valueDerivativeFEM(fe_linear, PHI_ID, 1);
        double phi_z = 0;
        if (nsd == 3) {
            phi_z = this->p_data_->valueDerivativeFEM(fe_linear, PHI_ID, 2);
            dDTdz = this->p_data_->valueFEM(fe_linear, 48);
        }
        double dV_norm = ((phi_x - dDTdx / absdDT * Dphi) * (phi_x - dDTdx / absdDT * Dphi) +
                          (phi_y - dDTdy / absdDT * Dphi) *
                          (phi_y - dDTdy / absdDT * Dphi));  // norm of derivative of phi from previous iteration
        if (nsd == 3) {
            dV_norm += (phi_z - dDTdz / absdDT * Dphi) * (phi_z - dDTdz / absdDT * Dphi);
        }
        dV_norm = sqrt(dV_norm);
        double DT = this->p_data_->valueFEM(fe_linear,
                                            DT_ID); // Distance contour, showing distance from transition layer
        double mun_dim = this->p_data_->valueFEM(fe_linear, MUN_ID); // dimensional form of mu_n from config file
        double mup_dim = this->p_data_->valueFEM(fe_linear, MUP_ID); // dimensional form of mu_p from config file
        double eps_dim = this->p_data_->valueFEM(fe_linear, EPS_ID); // dimensional form of epsilon from config file
        double mu_dim = mun_dim + mup_dim;
        double mu_x = input_data_->mu_hat(input_data_->mu_x);  // non dimensional form of mu_x
        double tauX = input_data_->c_unhat(input_data_->U_hat(1. / input_data_->tau_x));
        double G = this->p_data_->valueFEM(fe_linear, G_ID); ///get exciton generation
        double Gx = input_data_->U_hat(G);
        const ZEROPTV &pt = fe.position();
        double n_elec_pre = this->p_data_->valueFEM(fe_linear, N_ID); // dimensional form of mu_n from config file
        double p_hole_pre = this->p_data_->valueFEM(fe_linear, P_ID); // dimensional form of mu_p from config file
        double y_h = pt(1);        // number of DOF of system (n, p, phi, x)
        double kdiss = input_data_->c_unhat(input_data_->U_hat(gendiss->calcKdiss(y_h, dV_norm, DT, eps_dim,
                                                                                  mu_dim)));
        // K dissociation, this will be used for calculation of dissociation
        double Rnp = gendiss->calcR(y_h, n_elec_pre, p_hole_pre, DT, eps_dim, mu_dim);
        Rnp = input_data_->U_hat(Rnp);
        if (y_h > 0.95) {
            Gx = 0;
            Rnp = 0;
        }
        if (y_h < 0.05) {
            Gx = 0;
            Rnp = 0;
        }
        // in order to assemble the gauss point, we loop over each pair of basis
        // functions and calculate the individual contributions to the 'Ae' matrix
        // and 'be' vector.
        for (int a = 0; a < n_basis_functions; a++) {
            for (int b = 0; b < n_basis_functions; b++) {
                double M = fe.N(a) * fe.N(b) * detJxW;
                double N = 0;
                for (int k = 0; k < nsd; k++) {
                    N += fe.dN(a, k) * fe.dN(b, k) * detJxW;
                }
                // Add term to the A element matrix
                Ae(a, b) += (tauX + kdiss) * M + mu_x * N;
            }
            // Add term to the b element vector
            be(a) += fe.N(a) * (Gx + Rnp) * detJxW;
        }
    }

private:
    MPITimer timers_[4];          ///< for timing several parts of the code
    DDInputData *input_data_;     ///< pointer to input data
};

#endif // INCLUDE_DDEXCITON_H_

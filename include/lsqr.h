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
#ifndef INCLUDE_DDLSQR_H_
#define INCLUDE_DDLSQR_H_

#include "DDNodeData.h"
#include <vector>
typedef vector <float>vec;

/**
 * This class implements Least square method.
 *
 * The class can handle 1, 2, or 3 dimensions.daries as well as different number of G.Ps
 *
 */
class lsqr:public CEquation < DDNodeData > {
public:
    // indices to timing arrays. These are locations in the timers_ array that
    // correspond to specific stages of the code that we wish to time.
    static const int kTimerSolve = 0; ///< index to timer for solve process
    static const int kTimerAssemble = 1;  ///< index to timer for assemble
    static const int kTimerKSPSolve = 2;  ///< index to timer for KSPsolve
    static const int kTimerUpdate = 3;  ///< index to timer for update process
    ZEROMATRIX<double> J_n;
    int kk = 0;

    /**
     * Constucts the solver by setting up timers.
     */
    explicit lsqr(DDInputData *input_data_, AssemblyMethod assembly_method = kAssembleGaussPoints)
            : CEquation<DDNodeData>(assembly_method),
              input_data_
                      (input_data_) {
      timers_[kTimerSolve].set_label("Solve");
      timers_[kTimerAssemble].set_label("Assemble");
      timers_[kTimerKSPSolve].set_label("KSPSolve");
      timers_[kTimerUpdate].set_label("Update");

      // SetPreallocator(new PreallocatorPerfect<HTNodeData>(this));
    }

    /**
     * Destroys the object
     *
     * This also prints the average times required for each of the processes that
     * we times during execution.
     */
    virtual ~ lsqr() {
      timers_[kTimerSolve].PrintGlobalAverageSeconds();
      timers_[kTimerAssemble].PrintGlobalAverageSeconds();
      timers_[kTimerKSPSolve].PrintGlobalAverageSeconds();
      timers_[kTimerUpdate].PrintGlobalAverageSeconds();
    }

    /**
     * Sets up the essential (Direchlet) boundary conditions
     *
     * No boundary condition availabe for this methood.
     */
    virtual void fillEssBC() {
      // Setting boundary conditions.
      this->initEssBC();
    }

    /**
     * Solves the system for a given time value.
     *
     * @param t the current time value
     * @param jj and kk: ID of variable we want to store the solution in it
     */
    virtual void Solve(double jj, double t = 0.0) {
      int kk = int(jj);

      PrintInfo("  Filling essential boundary conditions...");
      fillEssBC();               // Set the boundary conditions
      PrintInfo("  Applying Ess BC to solution...");
      ApplyEssBCToSolution();    // apply boundary conditions to the solution vector

      PrintInfo("  Assembling...");
      Assemble(false);           // assemble the system

      PrintInfo("  Applying Ess BC...");
      ApplyEssBC();              // apply boundary conditions to the assembled system

      PetscViewer view;
      PetscViewerASCIIOpen(PETSC_COMM_WORLD, "mat.m", &view);
      PetscViewerPushFormat(view, PETSC_VIEWER_ASCII_MATLAB);
      MatView(Ag_, view);
      PetscViewerDestroy(&view);

      PetscViewerASCIIOpen(PETSC_COMM_WORLD, "vec.m", &view);
      PetscViewerPushFormat(view, PETSC_VIEWER_ASCII_MATLAB);
      VecView(bg_, view);
      PetscViewerDestroy(&view);

      PrintInfo("  Solving with KSP...");
      SolveKSP(solution_, 1, 0); // run the KSP solver (uses PETSc)

      // The result of the system solve is in the solution_ vector. We want to
      // store this data in our node data arrays where we keep track of the
      // current and previous heat values. This function call will put the data
      // from solution_ into location 0 of each NodeData object.
      PrintInfo("  Copying data from solution vector into node data...");
      for (int i = 0; i < p_grid_->n_nodes(); i++) {
        //double new_value = solution_(i*n_dof());  // the value for the ith node
        if (kk == 6 || kk == 7) {
          p_data_->GetNodeData(i).phi[kk - 5] = solution_(i * n_dof()); // store in node data object
        } else if (kk == 8 || kk == 9) {
          p_data_->GetNodeData(i).n[kk - 7] = solution_(i * n_dof()); // store in node data object
        } else if (kk == 10 || kk == 11) {
          p_data_->GetNodeData(i).p[kk - 9] = solution_(i * n_dof()); // store in node data object
        } else if (kk == 12 || kk == 13 || kk == 14) {
          p_data_->GetNodeData(i).j[kk - 6] = solution_(i * n_dof()); // store in node data object
        } else {
          p_data_->GetNodeData(i).j[0 + kk] = 0;
          for (int c = 0; c < input_data_->ngp; c++) {
            p_data_->GetNodeData(i).j[0 + kk] += solution_(i * n_dof() + c); // store in node data object
          }
          p_data_->GetNodeData(i).j[0 + kk] /= input_data_->ngp;
          if (p_data_->GetNodeData(i).j[0] < 0)
            p_data_->GetNodeData(i).j[0 + kk] *= -1;  // store in node data object;  // store in node data object
        }

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
      fe_linear.refill(fe.elem(), BASIS_LINEAR, 0);
      fe_linear.calc_at(fe.itg_pt());
      const int nbf = fe.nbf();  // # of basis functions
      const double detJxW = fe.detJxW(); // (determinant of J) cross W
      const int elm_id = fe.elem()->elm_id();
      const int lsq_ndof = input_data_->ngp;

      for (int a = 0; a < nbf; a++) {
        for (int b = 0; b < nbf; b++) {
          for (int ndofi = 0; ndofi < lsq_ndof; ndofi++) {
            for (int ndofj = 0; ndofj < lsq_ndof; ndofj++) {

              if (ndofi == ndofj) {
                Ae(lsq_ndof * a + ndofi, lsq_ndof * b + ndofj) += fe.N(a) * fe.N(b) * detJxW;
              } else {
                Ae(lsq_ndof * a + ndofi, lsq_ndof * b + ndofj) += 0;
              }
              if (b == 0) {
                if (ndofj == 0) {
                  be(lsq_ndof * a + ndofi) += J_n(elm_id, kk) * fe.N(a) *
                                              detJxW;//be(lsq_ndof * a +ndofi) += J_n(elm_id,ndofi) * fe.N(a)*detJxW;
                }
              }
            }
          }
        }
      }
    }

private:
    MPITimer timers_[4];          ///< for timing several parts of the code
    DDInputData *input_data_;     ///< pointer to input data
};

#endif // INCLUDE_DDEXCITON_H_

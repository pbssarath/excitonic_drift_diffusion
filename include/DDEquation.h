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
#ifndef DD_EQUATION_HPP
#define DD_EQUATION_HPP

#include "DDNodeData.h"
#include "initial_guess.h"
#include "gen_diss.hh"

using namespace std;

MPITimer timers[4];             ///< for timing several parts of the code

// indices to timing array
const int kTimerSolve = 0;
const int kTimerAssemble = 1;
const int kTimerSNESSolve = 2;
const int kTimerUpdate = 3;

// need to track the first assemble and last update time
// see destructor for details.
double g_first_assemble;
double g_last_update;

/**
 * Class calculates and returns value of A and b matrix for drift diffusion equations.
 * The boundary conditions are set up in DDEquation::fillEssBC()
 * The non-linear solver SNES parameters are set in DDEquation::Solve()
 * The solution is passed to GridField in each itteration in: UpdateGridField
 *
 */
class DDEquation:public CEquation < DDNodeData > {
public:
    gen_diss *gendiss;           ///< context for calculating k_diss
    float Vapp = 0;
    int cont = 8;                 ///< max iterations
    SNES snes;                    ///< nonlinear solver context
    int counter;                  ///< number of function calls
    DDInputData *idata;           ///< pointer to inputdata
    InitialGuess *init_;          ///< pointer to initial guess function + valueFEM function
    TezduyarUpwindFE PG_n;
    TezduyarUpwindFE PG_p;

    DDEquation(DDInputData *idata, InitialGuess *init_,
               bool has_uniform_mesh = false, AssemblyMethod assembly_method = kAssembleGaussPoints);

    virtual ~ DDEquation();

    virtual void fillEssBC();

    virtual void Solve(double dt, double t);

    virtual void Integrands(const FEMElm &fe, ZeroMatrix<double> &Ae, ZEROARRAY<double> &be);

    virtual void Integrands4side(const FEMElm &fe, int sideInd,
                                 ZeroMatrix<double> &Ae,
                                 ZEROARRAY<double> &be);
};

PetscErrorCode UpdateGridField (Vec _xg, DDEquation * ceqn, int flag = 0);
PetscErrorCode FormFunction (SNES snes, Vec x, Vec f, void *ceqn);
PetscErrorCode FormJacobian (SNES snes, Vec _xg, Mat jac, Mat B, void *ceqn_);

DDEquation::DDEquation (DDInputData * idata, InitialGuess * init_, bool has_uniform_mesh, AssemblyMethod assembly_method):CEquation < DDNodeData > (has_uniform_mesh, assembly_method),
idata (idata),
init_ (init_) {
    counter = 0;
    timers[kTimerSolve].set_label("Solve");
    timers[kTimerAssemble].set_label("Assemble");
    timers[kTimerSNESSolve].set_label("SNESSolve");
    timers[kTimerUpdate].set_label("Update");
    g_first_assemble = 0.0;
    g_last_update = 0.0;
}

DDEquation::~DDEquation () {
    timers[kTimerSolve].PrintGlobalTotalSeconds();
    timers[kTimerAssemble].PrintGlobalTotalSeconds();
    double assemble_mod = (timers[kTimerAssemble].GetTotalTimeSeconds() - g_first_assemble) * -1.0;
    double update_mod = (timers[kTimerUpdate].GetTotalTimeSeconds() - g_last_update) * -1.0;
    timers[kTimerSNESSolve].AddToTotalTime(assemble_mod);
    timers[kTimerSNESSolve].AddToTotalTime(update_mod);
    timers[kTimerSNESSolve].PrintGlobalTotalSeconds();
    timers[kTimerUpdate].PrintGlobalTotalSeconds();
}

void DDEquation::fillEssBC () {
    // Setting boundary conditions. for 2D problem
    this->initEssBC();
    double Vt = idata->Vt;        // scale factor for phi (0.0256)
    int n_boundaries;
    if (this->idata->nsd == 2) {
        if (this->idata->basisFunction == BASIS_HERMITE) {
            n_boundaries = N_COUPLED * 4; // number of boundaries for hermite basis
        } else {
            n_boundaries = N_COUPLED; // number of boundaries for linear basis
        }
    } else if (this->idata->nsd == 3) {
        if (this->idata->basisFunction == BASIS_HERMITE) {
            n_boundaries = 32;        // number of boundaries for hermite basis
        } else {
            n_boundaries = N_COUPLED; // number of boundaries for linear basis
        }
    }


    for (int nodeID = 0; nodeID < this->p_grid_->n_nodes(); nodeID++) {
        /**
         boundary IDs (second argument to BoNode):
            4
          _ _ _
         |     |
        1|     |2
         |_ _ _|
            3
            */
        if (this->p_grid_->BoNode(nodeID, 3)) {
            for (int i = 0; i < 3; i++) {
                if (idata->new_bc == 0)
                    this->specifyValue(nodeID, i, 0);
                if (idata->new_bc == 1 && i != 2)
                    this->specifyValue(nodeID, i, 0);
            }
            this->p_data_->GetNodeData(nodeID).phi[0] =
                    (idata->E_g - this->Vapp) / 2. / Vt; // this number can be changed for computing j-v curve
            this->p_data_->GetNodeData(nodeID).n[0] = 1; // this number should be 1 if N_c= N_v (N_c/N_c)
            this->p_data_->GetNodeData(
                    nodeID).p[0] = 1e-19; // this number is a small number. number of holes in cathode side
            this->p_data_->GetNodeData(nodeID).x[0] = 0;
        } else if (this->p_grid_->BoNode(nodeID, 4)) {
            for (int i = 0; i < 3; i++) {
                if (idata->new_bc == 0)
                    this->specifyValue(nodeID, i, 0);
                if (idata->new_bc == 1 && i != 1)
                    this->specifyValue(nodeID, i, 0);
            }
            this->p_data_->GetNodeData(nodeID).phi[0] =
                    -(idata->E_g - this->Vapp) / 2. / Vt;  // this number can be changed for computing j-v curve
            this->p_data_->GetNodeData(
                    nodeID).n[0] = 1e-19; // this number should be small number. number of electrons in anode side
            this->p_data_->GetNodeData(nodeID).p[0] = 1; // this number should be 1 (N_V/N_V)
            this->p_data_->GetNodeData(nodeID).x[0] = 0;
        }
    }
}

void DDEquation::Solve (double dt = 0.0, double t = 0.0) {
    timers[kTimerSolve].Start();
    PetscErrorCode ierr;

    // register dirichlet boundary conditions
    this->fillEssBC();

    // Assemble the Jacobian matrix and Function vector
    timers[kTimerAssemble].Start();
    if (idata->new_bc == 1) {
        this->Assemble(true /* surface integration */ );
    } else if (idata->new_bc == 0) {
        this->Assemble(false /* no surface integration */ );
    }
    timers[kTimerAssemble].Stop();

    g_first_assemble = timers[kTimerAssemble].GetLastTimeSeconds();

    // apply dirichlet boundary conditions
    this->ApplyEssBC();

    // Create nonlinear solver context
    ierr = SNESCreate(PETSC_COMM_WORLD, &snes);
    CHKERRABORT (PETSC_COMM_WORLD, ierr);

    // Set initial guess while checking if mesh
    // partitioning is used i.e. parallel_type_ == kWithDomainDecomp
    if (this->p_grid_->parallel_type_ == kWithDomainDecomp) {
        // array holds the local DOF values we want to put into the global
        // solution vector as our initial guess
        PetscScalar *array;
        ierr = PetscMalloc (this->n_total_dof() * sizeof(PetscScalar), &array);
        CHKERRABORT (PETSC_COMM_WORLD, ierr);

        // The global solution vector is stored as a distributed PETSc vector.
        // That means the storage for it is split up among all processors in the
        // system. Since the vector is distributed, we can't just write our local
        // DOF values to it when running in parallel.
        // However, PETSc has a "vector scatter" system that can be used to
        // "scatter" values from a local vector into a distributed vector.

        // So, we  create initVector (and use 'array' as the storage backing it).
        // initVector (aka 'array') will hold the DOF values stored on *this*
        // processor that we want to put into the global vector (this->xg_).
        // The mapping indices have been pre-computed by the CEquation base class
        // (which we inherit from) and are available as this->to() and this->from().
        Vec initVec;
        VecScatter scatter;
        ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, this->n_total_dof(), array, &initVec);
        CHKERRABORT (PETSC_COMM_WORLD, ierr);

        // fill 'array' with our local DOF values from p_data_ (the GridField)
        for (int nodeID = 0; nodeID < this->p_grid_->n_nodes(); nodeID++) {
            for (int i = 0; i < n_dof(); i++) {
                if (i % N_COUPLED == 0) {
                    array[n_dof() * nodeID + i] = this->p_data_->GetNodeData(nodeID).phi[int(i / N_COUPLED)];
                }
                if (i % N_COUPLED == 1) {
                    array[n_dof() * nodeID + i] = this->p_data_->GetNodeData(nodeID).n[int(i / N_COUPLED)];
                }
                if (i % N_COUPLED == 2) {
                    array[n_dof() * nodeID + i] = this->p_data_->GetNodeData(nodeID).p[int(i / N_COUPLED)];
                }
                if (i % N_COUPLED == 3) {
                    array[n_dof() * nodeID + i] = this->p_data_->GetNodeData(nodeID).x[int(i / N_COUPLED)];
                }
            }
        }

        // Use the PETSc scatter system to copy initVec into xg_
        ierr = VecScatterCreate(initVec, this->to(), this->xg_, this->from(), &scatter);
        CHKERRABORT (PETSC_COMM_WORLD, ierr);
        ierr = VecScatterBegin(scatter, initVec, this->xg_, INSERT_VALUES, SCATTER_FORWARD);
        CHKERRABORT (PETSC_COMM_WORLD, ierr);
        ierr = VecScatterEnd(scatter, initVec, this->xg_, INSERT_VALUES, SCATTER_FORWARD);
        CHKERRABORT (PETSC_COMM_WORLD, ierr);

        // free everything
        ierr = VecDestroy(&initVec);
        CHKERRABORT (PETSC_COMM_WORLD, ierr);
        ierr = VecScatterDestroy(&scatter);
        CHKERRABORT (PETSC_COMM_WORLD, ierr);
        ierr = PetscFree (array);
        CHKERRABORT (PETSC_COMM_WORLD, ierr);
    } else {
        // when we aren't running in distributed mode,
        // we can just directly copy the values from p_data_
        // into xg_ with VecSetValues (in chunks of size n_dof())
        double *val = new double[this->n_dof()];
        PetscInt *index = new PetscInt[this->n_dof()];
        for (int nodeID = 0; nodeID < this->p_grid_->n_nodes(); nodeID++) {
            for (int i = 0; i < this->n_dof(); i++) {
                if (i % N_COUPLED == 0) {
                    val[i] = this->p_data_->GetNodeData(nodeID).phi[int(i / N_COUPLED)];
                }
                if (i % N_COUPLED == 1) {
                    val[i] = this->p_data_->GetNodeData(nodeID).n[int(i / N_COUPLED)];
                }
                if (i % N_COUPLED == 2) {
                    val[i] = this->p_data_->GetNodeData(nodeID).p[int(i / N_COUPLED)];
                }
                if (i % N_COUPLED == 3) {
                    val[i] = this->p_data_->GetNodeData(nodeID).x[int(i / N_COUPLED)];
                }
                index[i] = this->n_dof() * nodeID + i;
            }

            ierr = VecSetValues(this->xg_, this->n_dof(), index, val, INSERT_VALUES);
            CHKERRABORT (PETSC_COMM_WORLD, ierr);
        }
        delete[]val;
        delete[]index;
        ierr = VecAssemblyBegin(this->xg_);
        CHKERRABORT (PETSC_COMM_WORLD, ierr);
        ierr = VecAssemblyEnd(this->xg_);
        CHKERRABORT (PETSC_COMM_WORLD, ierr);
    }

    // Set function evaluation routine and vector.
    ierr = SNESSetFunction(snes, this->bg_, FormFunction, this);
    CHKERRABORT (PETSC_COMM_WORLD, ierr);

    // Set jacobian matrix
    ierr = SNESSetJacobian(snes, this->Ag_, this->Ag_, FormJacobian, this);
    CHKERRABORT (PETSC_COMM_WORLD, ierr);

    // Set non-linear solver tolerances
    double rtol = 1e-14;
    int maxit = cont, maxf = 100000;
    counter = 0;                  // Resetting function counter
    SNESSetTolerances(snes, PETSC_DEFAULT, rtol, PETSC_DEFAULT, maxit, maxf);
    SNESSetFromOptions(snes);

    // Do the solve
    timers[kTimerSNESSolve].Start();
    ierr = SNESSolve(snes, PETSC_NULL, this->xg_);
    CHKERRABORT (PETSC_COMM_WORLD, ierr);
    timers[kTimerSNESSolve].Stop();

    // Get number of iterations used for logging
    PetscInt its;
    SNESGetIterationNumber(snes, &its);
    PetscPrintf(PETSC_COMM_WORLD, "number of iterations = %d\n", its);

    // Store solution to p_data_
    timers[kTimerUpdate].Start();
    UpdateGridField(this->xg_, this);
    timers[kTimerUpdate].Stop();
    g_last_update = timers[kTimerUpdate].GetLastTimeSeconds();

    // free the SNES object
    ierr = SNESDestroy(&snes);
    CHKERRABORT (PETSC_COMM_WORLD, ierr);

    timers[kTimerSolve].Stop();
}

/**
 * Returns the values of Ae[j,k] and be[j] for drift diffusion equations.
 * This code is written for solving 4 equations of drift diffusion (exiton and poison equations). 
 * The variables are named n, p, phi and x.
 * By using linear basis one can access and find these 4 unknowns.
 * By using Hermite one can find these variables + their derivative + jn and jp + grad of jn and jp in the whole domain.
 *
 * There are 4 loops in this integrands:
 * The 1st one is looping over the basis functions (j), for linear 0<=j<4 and for hermite 0<=j<12
 * The 2nd one is looping over the basis functions (k), for linear 0<=k<4 and for hermite 0<=k<12
 * The 3rd one is looping over degrees of freedom (ndof_i), for this system 0<=ndof_i<3 (n,p,phi)
 * The 4th one is looping over degrees of freedom (ndof_j), for this system 0<=ndof_j<3 (n,p,phi)
 * non_dim_fact: non dimensionalized factore for dissosiation.
 */
void DDEquation::Integrands (const FEMElm & fe, ZeroMatrix < double >&Ae, ZEROARRAY < double >&be) {
    /***
    * General parameters of finite element system:
    */
    FEMElm fe_linear(p_grid_, BASIS_ALL);
    if (fe.basis_function() == BASIS_HERMITE) {
        fe_linear.refill(fe.elem(), BASIS_LINEAR, 0);
        fe_linear.calc_at(fe.itg_pt());
    } else {
        fe_linear = fe;
    }
    const int nsd = fe.nsd();    // number of spatial dimensions of system (1, 2, or 3)
    const int nbf = fe.nbf();    // number of basis functions of system (4, or 16) for 2D
    const double detJxW = fe.detJxW(); ///< Jacobian of the system
    const int nndof = N_COUPLED;  // number of DOF of system (n, p, phi, x)
    /***
    * Computations regarding to previous iteration:
    */
    ZEROPTV dphi_pre;             // derivative of "phi" from last iteration in x, y and z direction
    ZEROPTV dn_elec_pre;          // derivative of "n" from last iteration in x, y and z direction
    ZEROPTV dp_hole_pre;          // derivative of "p" from last iteration in x, y and z direction
    double n_elec_pre = 0;        // previous "n" from last iteration
    double p_hole_pre = 0;        // previous "p" from last iteration
    double X_pre = 0;             // previous "x" from last iteration
    double G = this->p_data_->valueFEM(fe_linear, G_ID); ///get exciton generation
    double D2phi = this->p_data_->valueFEM(fe_linear, D2PHIBAR_ID);
    double Dphi = this->p_data_->valueFEM(fe_linear, PHIBAR_ID);
    double dDTdx = this->p_data_->valueFEM(fe_linear, 32);
    double dDTdy = this->p_data_->valueFEM(fe_linear, 33);
    double dDTdz = this->p_data_->valueFEM(fe_linear, 48);
    double absdDT;
    if (nsd == 2) {
        absdDT = sqrt(dDTdx * dDTdx + dDTdy * dDTdy) + 1e-16;
    } else {
        absdDT = sqrt(dDTdx * dDTdx + dDTdy * dDTdy + dDTdz * dDTdz) + 1e-16;
    }
    const ZEROPTV &pt = fe.position();
    double y_h = pt(1);
    // variable 0: phi
    // variable 1: n
    // variable 2: p

    if (fe.basis_function() == BASIS_HERMITE) {
        n_elec_pre = this->init_->valueFEM(fe, N_ID, nndof); // Electron density from last iteration
        p_hole_pre = this->init_->valueFEM(fe, P_ID, nndof); // hole density from last iteration
        X_pre = this->p_data_->valueFEM(fe, X_ID); // Exciton density from last iteration
        for (int i = 0; i < nsd; i++) {
            dphi_pre(i) = this->init_->valueDerivativeFEM(fe, PHI_ID, i,
                                                          nndof); // derivative of Electric Field from last iteration
            dn_elec_pre(i) = this->init_->valueDerivativeFEM(fe, N_ID, i,
                                                             nndof); // derivative of  Electron from last iteration
            dp_hole_pre(i) = this->init_->valueDerivativeFEM(fe, P_ID, i,
                                                             nndof); // derivative of  hole from last iteration
        }
    } else {
        n_elec_pre = this->p_data_->valueFEM(fe, N_ID); // Electron density from last iteration
        p_hole_pre = this->p_data_->valueFEM(fe, P_ID); // hole density from last iteration
        X_pre = this->p_data_->valueFEM(fe_linear, X_ID); // Exciton density from last iteration
        for (int i = 0; i < nsd; i++) {
            dphi_pre(i) = this->p_data_->valueDerivativeFEM(fe, PHI_ID,
                                                            i); // derivative of Electric Field from last iteration
            dn_elec_pre(i) = this->p_data_->valueDerivativeFEM(fe, N_ID,
                                                               i); // derivative of Electron from last iteration
            dp_hole_pre(i) = this->p_data_->valueDerivativeFEM(fe, P_ID, i); // derivative of hole from last iteration
        }
    }
    /***
    * Calculation of k dissociation, Lambda, Recombination and SUPG weights:
    */
    double dV_norm = ((dphi_pre(0) - dDTdx / absdDT * Dphi) * (dphi_pre(0) - dDTdx / absdDT * Dphi)
                      + (dphi_pre(1) - dDTdy / absdDT * Dphi) *
                        (dphi_pre(1) - dDTdy / absdDT * Dphi));  // norm of derivative of phi from previous iteration
    if (idata->nsd == 3) dV_norm += (dphi_pre(2) - dDTdz / absdDT * Dphi) * (dphi_pre(2) - dDTdz / absdDT * Dphi);
    dV_norm = sqrt(dV_norm);
    double DT = this->p_data_->valueFEM(fe_linear, DT_ID); // Distance contour, showing distance from transition layer
    double mun_dim = this->p_data_->valueFEM(fe_linear, MUN_ID); // dimensional form of mu_n from config file
    double mup_dim = this->p_data_->valueFEM(fe_linear, MUP_ID); // dimensional form of mu_p from config file
    double eps_dim = this->p_data_->valueFEM(fe_linear, EPS_ID); // dimensional form of epsilon from config file
    double mu_dim = mun_dim + mup_dim;
    double kdiss = idata->c_unhat(idata->U_hat(gendiss->calcKdiss(y_h, dV_norm, DT, eps_dim,
                                                                  mu_dim)));  // K dissosiation, non dimensional form. this will be used for calculation of dissosiation
    double lambda_2 = idata->lambda2_(eps_dim);  // Lambda sq.
    double Rnp = gendiss->calcR(y_h, n_elec_pre, p_hole_pre, DT, eps_dim, mu_dim); // Recombination calculation
    double D = kdiss * X_pre; // Dissociation calculation
    Rnp = idata->U_hat(Rnp); // Non dimensional form of recombination
    double gamma = 0;
    if (n_elec_pre * p_hole_pre != 0)
        gamma = Rnp / (n_elec_pre * p_hole_pre); // non dimensional form of gamma, used in Jacobian
    ZEROPTV U_n; // Advection term for electron
    ZEROPTV U_p; // Advection term for hole
    U_n(0) = dphi_pre(0); // x component term of advection for electron
    U_n(1) = dphi_pre(1); // y component term of advection for electron
    U_n(2) = 0;
    U_p(0) = -dphi_pre(0); // x component term of advection for hole
    U_p(1) = -dphi_pre(1); // y component term of advection for hole
    U_p(2) = 0;
    if (idata->nsd == 3) {
        U_n(2) = dphi_pre(2);
        U_p(2) = -dphi_pre(2);
    }
    PG_n.calcSUPGWeightingFunction(fe, U_n, 1);
    PG_p.calcSUPGWeightingFunction(fe, U_p, 1);
    /***
    * Computations of non dimensional form of input parameters:
    */
    double mu_n = idata->mu_hat(mun_dim);  // non dimensional form of mu_n
    double mu_p = idata->mu_hat(mup_dim);  // non dimensional form of mu_p

    // Looping over basis functions
    for (int j = 0; j < nbf; j++) { // Looping over basis functions
        // in 2D: nbf for linear is 4, Quad is 9, Cubic and Hermite are 16
        for (int k = 0; k < nbf; k++) { // Looping of over coupled variables (phi, n, p)
            // ndofi = 0 --> phi
            // ndofi = 1 --> n
            // ndofi = 2 --> p
            for (int n_dofi = 0; n_dofi < nndof; n_dofi++) { // Looping of over coupled variables (phi, n, p)
                // ndofj = 0 --> phi
                // ndofj = 1 --> n
                // ndofj = 2 --> p
                for (int n_dofj = 0; n_dofj < nndof; n_dofj++) {
                    /**
                     Equation for potential (Poison's Eq.)-------> \nabla.(\lambda^{2} \varphi) - n + p = 0
                     @ return --> A11, A12, A13 and b1
                    */
                    if (n_dofi ==
                        0) { // This if loop returns A11 and b1. Details of these calculation can be found in formula sheet.
                        // A11 in Latex form is: -(\nabla \omega, \lambda^{2} \delta \varphi) = 0 ----> -(dN(i,k), \lambda^{2} dN(j,k))
                        if (n_dofj == 0) {
                            double N = 0;
                            for (int index = 0; index < nsd; index++) {
                                N += fe.dN(j, index) * fe.dN(k, index);
                            }
                            Ae(nndof * j + n_dofi, nndof * k + n_dofj) -=
                                    N * lambda_2 * detJxW; ///< A[11] in Eq. 4 in formula sheet
                            // This if loop returns b1 and make sure that b1 is computed by looping over basis functions and number of coupled variables (just one time).
                            //b1 in Latex form is: (dN(i,k), \lambda^{2} d\varphi_{,k}) - (N(i) , n_{pre}) - (N(i) , p_{pre}) = 0;
                            if (k == 0) {
                                double dphi_dN = 0;
                                for (int index = 0; index < nsd; index++) {
                                    dphi_dN += dphi_pre(index) * fe.dN(j, index);
                                }
                                be(nndof * j + n_dofi) += (-dphi_dN * lambda_2 - (fe.N(j) + PG_n.SUPG(j)) *
                                                                                 (n_elec_pre - p_hole_pre -
                                                                                  D2phi * eps_dim)) *
                                                          detJxW;  ///< b[1] in Eq. 5 in formula sheet
                            }
                        }
                            // This if loop returns A12
                            // A12 in Latex form is: - (\omega, \delta n) = 0 ----> -(N(i), N(j))
                        else if (n_dofj == 1) {
                            Ae(nndof * j + n_dofi, nndof * k + n_dofj) -= (fe.N(j) + PG_n.SUPG(j)) * fe.N(k) *
                                                                          detJxW;  ///< A[12] in Eq. 4 in formula sheet
                        }
                            // This if loop returns A13
                            // A13 in Latex form is: (\omega, \delta p) = 0 ----> (N(i), N(j))
                        else if (n_dofj == 2) {
                            Ae(nndof * j + n_dofi, nndof * k + n_dofj) +=
                                    (fe.N(j) + PG_p.SUPG(j)) * fe.N(k) * detJxW;  ///< A[13] in Eq. 4 in formula sheet
                        }
                    }
                        /**
                         Equation for electron (Drift Diffusion Eq.)-------> \nabla.(- \mu_{n} n \nabla \varphi + \mu_{n} \nabla n) - (-D + R) = 0
                         @ return --> A21, A22, A23 and b2
                        */
                    else if (n_dofi ==
                             1) { // This if loop returns A21 and b2. Details of these calculation can be found in formula sheet.
                        // A21 in Latex form is: (\nabla \omega, \mu_{n} n \nabla \delta \varphi) ----> (dN(i,k), \mu_{n} n dN(j,k))
                        if (n_dofj == 0) {
                            double N = 0;
                            double NN = 0;
                            for (int index = 0; index < nsd; index++) {
                                N += fe.dN(j, index) * fe.dN(k, index);
                            }
                            Ae(nndof * j + n_dofi, nndof * k + n_dofj) +=
                                    mu_n * n_elec_pre * N * detJxW;  ///< A[21] in Eq. 10 in formula sheet
                            // This if loop returns b1 and make sure that b2 is computed by looping over basis functions and number of coupled variables (just one time).
                            // b2 in Latex form is: -(dN(i,k),-\mu_{n} n^{pre} d\varphi_{,k} + \mu_{n} dn^{pre}_{,k}) + (N(i),K_{diss}X^{pre}) - (N(i),\gamma n^{pre} p^{pre})
                            if (k == 0) {
                                double dn_dN = 0;
                                double dphi_n_dN = 0;
                                for (int index = 0; index < nsd; index++) {
                                    dn_dN += dn_elec_pre(index) * fe.dN(j, index);
                                    dphi_n_dN += dphi_pre(index) * fe.dN(j, index);
                                }
                                be(nndof * j + n_dofi) +=
                                        (mu_n * n_elec_pre * dphi_n_dN - mu_n * dn_dN + D * fe.N(j) - Rnp * fe.N(j) -
                                         (-D + Rnp) * PG_n.SUPG(j)) * detJxW;
                                ///< b[2] in Eq. 11 in formula sheet
                            }
                        }
                            // This if loop returns A22
                            // (\nabla \omega, \mu_{n} \nabla \varphi \delta n - \mu_{n} \nabla \delta n)) + (\omega, - \gamma \delta n p) = 0 ---->
                            // (dN(i,k), \mu_{n} \varphi_{,k} N(j) - \mu_{n} dN(j,k)) + (N(i), - \gamma p N(j))
                        else if (n_dofj == 1) {
                            double N = 0;
                            for (int index = 0; index < nsd; index++) {
                                N += fe.dN(j, index) * fe.dN(k, index);
                            }
                            double K_phi = 0;
                            for (int index = 0; index < nsd; index++) {
                                K_phi += fe.dN(j, index) * fe.N(k) * dphi_pre(index);
                            }
                            Ae(nndof * j + n_dofi, nndof * k + n_dofj) -=
                                    (mu_n * (N - K_phi) + gamma * p_hole_pre * (fe.N(j) + PG_n.SUPG(j)) * fe.N(k)) *
                                    detJxW; ///< A[22] in Eq. 10 in formula sheet
                        }
                            // This if loop returns A23
                            // (\omega, - \gamma n \delta p) = 0 ----> (N(i), - \gamma n N(j))
                        else if (n_dofj == 2) {
                            Ae(nndof * j + n_dofi, nndof * k + n_dofj) -=
                                    gamma * n_elec_pre * (fe.N(j) + PG_n.SUPG(j)) * fe.N(k) *
                                    detJxW; ///< A[23] in Eq. 10 in formula sheet
                        }
                    }
                        /**
                         Equation for hole (Drift Diffusion Eq.)-------> \nabla.(- \mu_{p} p \nabla \varphi - \mu_{p} \nabla p) + (-D + R) = 0
                         @ return --> A31, A32, A33 and b2
                        */
                    else if (n_dofi == 2) { // This if loop returns A31 and b3
                        //A31 in Latex form is: -(\nabla \omega, \mu_{p} p \nabla \delta \varphi) = 0 ----> (dN(i,k), \mu_{p} p dN(j,k))
                        if (n_dofj == 0) {
                            double N = 0;
                            for (int index = 0; index < nsd; index++) {
                                N += fe.dN(j, index) * fe.dN(k, index);
                            }
                            Ae(nndof * j + n_dofi, nndof * k + n_dofj) +=
                                    mu_p * p_hole_pre * N * detJxW;  ///< A[31] in Eq. 13 in formula sheet
                            // This if loop returns b1 and make sure that b3 is computed by looping over basis functions and number of coupled variables (just one time).
                            if (k == 0) {
                                double dp_dN = 0;
                                double dphi_p_dN = 0;
                                for (int index = 0; index < nsd; index++) {
                                    dp_dN += dp_hole_pre(index) * fe.dN(j, index);
                                    dphi_p_dN += dphi_pre(index) * fe.dN(j, index);
                                }
                                be(nndof * j + n_dofi) +=
                                        (mu_p * (dp_dN + p_hole_pre * dphi_p_dN) - D * fe.N(j) + Rnp * fe.N(j) +
                                         (-D + Rnp) * PG_p.SUPG(j)) * detJxW;
                                ///< b[3] in Eq. 14 in formula sheet
                            }
                        }
                            // This if loop returns A32
                            // A32 in Latex form is: (\omega, - \gamma \delta n p) = 0 ---> (N(i), \gamma p N(j))
                        else if (n_dofj == 1) {
                            Ae(nndof * j + n_dofi, nndof * k + n_dofj) +=
                                    gamma * p_hole_pre * (fe.N(j) + PG_p.SUPG(j)) * fe.N(k) *
                                    detJxW; ///< A[32] in Eq. 13 in formula sheet
                        }
                            // This if loop returns A33
                            // A33 in Latex form is: - (\nabla \omega, \mu_{p} \nabla \varphi \delta p + \mu_{p} \nabla \delta p)) + (\omega, - \gamma n \delta p) = 0 ---->
                            // (dN(i,k), \mu_{p} \varphi_{,k} N(j) + \mu_{p} dN(j,k)) + (N(i), n \gamma N(j))
                        else if (n_dofj == 2) {
                            double N = 0;
                            double K_phi = 0;
                            for (int index = 0; index < nsd; index++) {
                                N += fe.dN(j, index) * fe.dN(k, index);
                                K_phi += fe.dN(j, index) * fe.N(k) * dphi_pre(index);
                            }
                            Ae(nndof * j + n_dofi, nndof * k + n_dofj) +=
                                    (mu_p * (N + K_phi) + gamma * n_elec_pre * (fe.N(j) + PG_p.SUPG(j)) * fe.N(k)) *
                                    detJxW; ///< A[33] in Eq. 13 in formula sheet
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i < be.size(); i++) {
        assert(!std::isnan(be(i)));
        //if (std::isnan(be(i)))
        //  std::cout << "===== be(" << i << ") == NaN!\n";
    }
}

/**
 Defining Neumann boundary condition of the form \surf_int w*(grad (u))dA
 which is a surface integral, where A is the surface area
 Only if the boundary is Neumann
 */
void DDEquation::Integrands4side(const FEMElm& fe, int sideInd,
                               ZeroMatrix<double>& Ae,
                               ZEROARRAY<double>& be) {
        // sideInd == 1 means left boundary
        // sideInd == 2 means right boundary
        // sideInd == 3 means bottom boundary
        // sideInd == 4 means top boundary

        ZEROPTV n_hat;
        if (sideInd == 3 || sideInd == 4) {
            if (sideInd == 3) {
                n_hat(0) = 0;
                n_hat(1) = 1;
            } else if (sideInd == 4) {
                n_hat(0) = 0;
                n_hat(1) = -1;
            }

            FEMElm fe_linear(p_grid_, BASIS_ALL);
            if (fe.basis_function() == BASIS_HERMITE) {
                fe_linear.refill(fe.elem(), BASIS_LINEAR, 0);
                fe_linear.calc_at(fe.itg_pt());
            } else {
                fe_linear = fe;
            }
            const ZEROPTV &normal = fe.surface()->normal(); // Normal of the boundary side
            const int nsd = fe.nsd();    // number of spatial dimensions of system (1, 2, or 3)
            const int nbf = fe.nbf();    // number of basis functions of system (4, or 16) for 2D
            const double detJxW = fe.detJxW(); ///< Jacobian of the system
            const int nndof = N_COUPLED;  // number of DOF of system (n, p, phi, x)
            ZEROPTV dphi_pre;             // derivative of "phi" from last iteration in x, y and z direction
            ZEROPTV dn_elec_pre;          // derivative of "n" from last iteration in x, y and z direction
            ZEROPTV dp_hole_pre;          // derivative of "p" from last iteration in x, y and z direction
            double n_elec_pre = 0;        // previous "n" from last iteration
            double p_hole_pre = 0;        // previous "p" from last iteration
            double phi_pre = 0;
            // variable 0: phi
            // variable 1: n
            // variable 2: p
            // variable 3: x
            n_elec_pre = this->p_data_->valueFEM(fe, 1);
            p_hole_pre = this->p_data_->valueFEM(fe, 2);
            phi_pre = this->p_data_->valueFEM(fe, 0);
            for (int i = 0; i < nsd; i++) {
                dphi_pre(i) = this->p_data_->valueDerivativeFEM(fe, 0, i);
                dn_elec_pre(i) = this->p_data_->valueDerivativeFEM(fe, 1, i);
                dp_hole_pre(i) = this->p_data_->valueDerivativeFEM(fe, 2, i);
            }
            double mu_n, mu_p;  // non dimensional form of mu_n, mu_p
            double u_n = 0;
            double u_p = 0;
            ZEROPTV dN_hat;
            ZEROPTV dN_hat1;
            double fact = 1;
            // h in penalty term
            ZEROPTV h_elm(idata->L[0] / idata->Nelem[0],
                          idata->L[1] / idata->Nelem[1],
                          idata->L[2] / idata->Nelem[2]);
            double h = sqrt(h_elm(0) * h_elm(1));
            if (sideInd == 3) {
                u_n = 1.00000000001;
                u_p = 1e-19;
                mu_n = 1.;
                mu_p = 0.01;
            } else {
                u_n = 1e-19;
                u_p = 1.00000000001;
                mu_n = 0.01;
                mu_p = 1.;
            }
            for (int j = 0; j < nbf; j++) {
                dN_hat(0) = fe.dN(j, 0);
                dN_hat(1) = fe.dN(j, 1);
                for (int k = 0; k < nbf; k++) {
                    dN_hat1(0) = fe.dN(k, 0);
                    dN_hat1(1) = fe.dN(k, 1);
                    for (int n_dofi = 0; n_dofi < 3; n_dofi++) {
                        for (int n_dofj = 0; n_dofj < 3; n_dofj++) {
                            if (n_dofi == 1) {
                                if (n_dofj == 0) {
                                    if (k == 0) {
                                        // Consistency term for n
                                        be(3 * j + n_dofi) += mu_n * fe.N(j) *
                                                              (-n_elec_pre * dphi_pre.innerProduct(normal) +
                                                               dn_elec_pre.innerProduct(normal)) * detJxW;
                                        // Adjoint term
                                        be(3 * j + n_dofi) += mu_n * (n_elec_pre - u_n) * (dN_hat.innerProduct(normal) -
                                                                                           fe.N(j) *
                                                                                           dphi_pre.innerProduct(
                                                                                                   normal)) * detJxW;
                                        // Penalty term
                                        be(3 * j + n_dofi) += 16 / h * mu_n * (n_elec_pre - u_n) * fe.N(j) * detJxW;
                                    }
                                    Ae(3 * j + n_dofi, 3 * k + n_dofj) +=
                                            mu_n * fe.N(j) * (-n_elec_pre * dN_hat1.innerProduct(normal)) * detJxW;
                                }
                                if (n_dofj == 1) {
                                    Ae(3 * j + n_dofi, 3 * k + n_dofj) += mu_n * fe.N(j) *
                                                                          (-fe.N(k) * dphi_pre.innerProduct(normal) +
                                                                           dN_hat1.innerProduct(normal)) * detJxW;
                                    // Adjoint term
                                    Ae(3 * j + n_dofi, 3 * k + n_dofj) += mu_n * (fe.N(k) - 0) *
                                                                          (dN_hat.innerProduct(normal) -
                                                                           fe.N(j) * dphi_pre.innerProduct(normal)) *
                                                                          detJxW;
                                    // Penalty term
                                    Ae(3 * j + n_dofi, 3 * k + n_dofj) +=
                                            mu_n * 16 / h * (fe.N(k) - 0) * fe.N(j) * detJxW;
                                }
                                if (n_dofj == 2) {
                                }
                            } else if (n_dofi == 2) {
                                if (n_dofj == 0) {
                                    if (k == 0) {
                                        // Consistency term for p
                                        be(3 * j + n_dofi) += mu_p * fe.N(j) *
                                                              (p_hole_pre * dphi_pre.innerProduct(normal) +
                                                               dp_hole_pre.innerProduct(normal)) * detJxW;
                                        // Adjoint term
                                        be(3 * j + n_dofi) += mu_p * (p_hole_pre - u_p) * (dN_hat.innerProduct(normal) -
                                                                                           fe.N(j) *
                                                                                           dphi_pre.innerProduct(
                                                                                                   normal)) * detJxW *
                                                              fact;
                                        // Penalty term
                                        be(3 * j + n_dofi) +=
                                                16 / h * mu_p * (p_hole_pre - u_p) * fe.N(j) * detJxW * fact;
                                    }
                                    Ae(3 * j + n_dofi, 3 * k + n_dofj) +=
                                            mu_p * fe.N(j) * (-p_hole_pre * dN_hat1.innerProduct(normal)) * detJxW *
                                            fact;
                                }
                                if (n_dofj == 1) {
                                }
                                if (n_dofj == 2) {
                                    Ae(3 * j + n_dofi, 3 * k + n_dofj) += mu_p * fe.N(j) *
                                                                          (fe.N(k) * dphi_pre.innerProduct(normal) +
                                                                           dN_hat1.innerProduct(normal)) * detJxW *
                                                                          fact;
                                    // Adjoint term
                                    Ae(3 * j + n_dofi, 3 * k + n_dofj) += mu_p * (fe.N(k) - 0) *
                                                                          (dN_hat.innerProduct(normal) -
                                                                           fe.N(j) * dphi_pre.innerProduct(normal)) *
                                                                          detJxW * fact;
                                    // Penalty term
                                    Ae(3 * j + n_dofi, 3 * k + n_dofj) +=
                                            16 / h * mu_p * (fe.N(k) - 0) * fe.N(j) * detJxW * fact;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

   
/***
 * Copies data from the global PETSc solution vector (_xg) into
 * ceqn->solution_, and then copies ceqn->solution_ into the
 * GridField (ceqn->p_data_).
 * @param _xg global solution vector
 * @param ceqn DDEquation instance to pull grid values from
 */
PetscErrorCode UpdateGridField (Vec _xg, DDEquation * ceqn, int flag) {
    PetscErrorCode ierr;
    if (ceqn->p_grid_->parallel_type_ == kWithDomainDecomp) {
        // Create a PETSc vector backed by ceqn->solution_.data() so we can
        // scatter from the global vector (_xg) into local ceqn->solution_ vector.
        VecScatter scatter;
        Vec SolutionVec;
        ierr = VecCreateSeqWithArray(PETSC_COMM_SELF, 1, ceqn->n_total_dof(), ceqn->solution_.data(), &SolutionVec);
        CHKERRQ (ierr);

        // do the scatter
        ierr = VecScatterCreate(_xg, ceqn->from(), SolutionVec, ceqn->to(), &scatter);
        CHKERRQ (ierr);
        ierr = VecScatterBegin(scatter, _xg, SolutionVec, INSERT_VALUES, SCATTER_FORWARD);
        CHKERRQ (ierr);
        ierr = VecScatterEnd(scatter, _xg, SolutionVec, INSERT_VALUES, SCATTER_FORWARD);
        CHKERRQ (ierr);

        // clean up
        ierr = VecDestroy(&SolutionVec);
        CHKERRQ (ierr);
        ierr = VecScatterDestroy(&scatter);
        CHKERRQ (ierr);
    } else {
        // send all of _xg to every processor into solution_vec
        // VecScatterCreateToAll also creates the solution_vec vector
        VecScatter scatter;
        Vec solution_vec;
        ierr = VecScatterCreateToAll(_xg, &scatter, &solution_vec);
        CHKERRQ (ierr);

        // do the scatter
        ierr = VecScatterBegin(scatter, _xg, solution_vec, INSERT_VALUES, SCATTER_FORWARD);
        CHKERRQ (ierr);
        ierr = VecScatterEnd(scatter, _xg, solution_vec, INSERT_VALUES, SCATTER_FORWARD);
        CHKERRQ (ierr);

        // copy the solution_vec vector (containing the new data)
        // into ceqn->solution_
        double *array;
        ierr = VecGetArray(solution_vec, &array);
        CHKERRQ (ierr);
        memcpy(ceqn->solution_.data(), array, sizeof(double) * (ceqn->n_total_dof()));
        ierr = VecRestoreArray(solution_vec, &array);
        CHKERRQ (ierr);

        // free the PETSc scatter and temporary vector
        ierr = VecScatterDestroy(&scatter);
        CHKERRQ (ierr);
        ierr = VecDestroy(&solution_vec);
        CHKERRQ (ierr);
    }
    if (ceqn->idata->isPeriodic) ceqn->updatePeriodicSol();
    // Copy ceqn->solution_'s values into p_data_ (GridField)
    for (int i = 0; i < ceqn->p_grid_->n_nodes(); i++) {
        for (int j = 0; j < ceqn->n_dof(); j++) {
            if (j % N_COUPLED == 0) {
                ceqn->p_data_->GetNodeData(i).phi[int(j / N_COUPLED)] = ceqn->solution_(i * ceqn->n_dof() + j);
            }
            if (j % N_COUPLED == 1) {
                ceqn->p_data_->GetNodeData(i).n[int(j / N_COUPLED)] = ceqn->solution_(i * ceqn->n_dof() + j);
            }
            if (j % N_COUPLED == 2) {
                ceqn->p_data_->GetNodeData(i).p[int(j / N_COUPLED)] = ceqn->solution_(i * ceqn->n_dof() + j);
            }
            if (j % N_COUPLED == 3) {
                ceqn->p_data_->GetNodeData(i).x[int(j / N_COUPLED)] = ceqn->solution_(i * ceqn->n_dof() + j);
            }
        }

        // ceqn->p_data_->GetNodeData (i).value (3) = 0.001;
    }
    return ierr;
}

/***
 * Computes the global matrix and vector.
 */
PetscErrorCode FormFunction (SNES snes, Vec _xg, Vec _f, void *ceqn_) {
    PetscErrorCode ier;
    DDEquation *ceqn = (DDEquation *) ceqn_;
    ceqn->counter++;

    // copy guess into the GridField so we can use its values
    // when calculating residual
    timers[kTimerUpdate].Start();
    UpdateGridField(_xg, ceqn);
    timers[kTimerUpdate].Stop();

    // the normal solve step (set dirichlet, assemble, apply dirichlet)
    ceqn->fillEssBC();
    timers[kTimerAssemble].Start();
    if (ceqn->idata->new_bc == 1) {
        ceqn->Assemble(true /* no surface integration */ );
    } else if (ceqn->idata->new_bc == 0) {
        ceqn->Assemble(false /* no surface integration */ );
    }
    timers[kTimerAssemble].Stop();
    ceqn->ApplyEssBC();

    // Assemble() fills ceqn->bg_, but we need to store the result in f
    ier = VecCopy(ceqn->bg(), _f);
    CHKERRQ (ier);
    return 0;
}

/***
 * The Jacobian is already calculated by FormFunction's Assemble() call.
 * We do both at once so we can reuse the basis function values at each element
 * (which take a significant amount of time to calculate).
 * So, we do nothing here.
 */
PetscErrorCode FormJacobian (SNES snes, Vec _xg, Mat jac, Mat B, void *ceqn_) {
    return 0;
}

#endif

/* ------------------------------------------------------------------------
 *  Created: May 8, 2017
 *
 *  Authors: Ramin Noruzi, Pengfei Du, Baskar Ganapathysubramanian
 *  Copyright (c) 2017 Baskar Ganapathysubramanian
 *  See accompanying LICENSE.
  ------------------------------------------------------------------------- */
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
#include <talyfem/talyfem.h>
#include <math.h>
#include <string>
#include "DDInputData.h"
#include "DDNodeData.h"
#include "DDGridField.h"
#include "lsqr.h"
#include "DDEquation.h"
#include "DDExciton.h"
#include "../include/DT.hh"
#include "initial_guess.h"
#include "gen_diss.hh"

#define BUILDINFO_DD_IMPL
#include "build_info_DD.h"

using namespace TALYFEMLIB;
static char help[] = "Solves a steady state Excitonic Drift-Diffusion problem!";


inline bool SetIC (DDGridField & data, DDInputData & idata) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    switch (idata.typeOfIC) {
        case 0:
            load_gf(&data, &idata);
            data.UpdateDataStructures();
            return true;
        case 1:
            data.SetIC(idata.nsd);
            return true;
        default:
            if (rank == 0)
                std::cerr << "IC not set up " << std::endl;
            return false;
    }
}

/*
  This code solves STEADY STATE drift-diffusion equations in organic solar cell.
  It solves electron(n) hole(h) exciton(X) densities and electric field (V) distribution.

  n,p, phi are coupled together and X is solved separately. This is iterated to reduce the residual
  to get the converged result. After having n,p,phi,x values we find J values in G.Ps (in post process function).
  By usage of least-square method nodal values of J are computed. 
*/
int DDNodeData::nsd = 0;
int main (int argc, char **args) {
    PetscInitialize(&argc, &args, (char *) nullptr, help);

    PetscBool write_repro = PETSC_FALSE;
    Repro r(argc, args);
    r.add_build_info<BuildInfo_DD>();
    r.write();

    try {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        PetscMPIInt mpi_rank, mpi_size;
        MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);
        MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank);
        //define the solver related parameters;
        DDInputData inputData;
        GRID *pGrid = nullptr;
        {
            if ((argc == 2) && (strcmp(args[1], "--printParameters") == 0)) {
                if (rank == 0) {
                    inputData.printAll(std::cout);
                }
                PetscFinalize();
                return -1;
            }
            // Read input data from file "config.txt"
            inputData.ReadFromFile(std::string("config.txt"));
            // Loging inputdata to std::clog or file
            if (rank == 0)
                std::clog << inputData;
            //  Check if inputdata is complete
            if (!inputData.CheckInputData()) {
                throw TALYException() << "[ERR] Problem with input data, check the config file!";
            }

            DDNodeData::nsd = inputData.nsd;
            // Based on inputdata create Grid
            CreateGrid(pGrid, &inputData);
            // when loading a mesh, fix the NSD (it is always 3 for Gmsh meshes)
            if (!inputData.ifBoxGrid)
                pGrid->set_nsd(inputData.nsd);
            DDGridField data(&inputData);

            // set up periodic boundary object -- always periodic.
            PeriodicBounds *pbc = nullptr;
            if (inputData.nsd == 2) {
                int *pbc_indices = new int[1];
                int pbc_count = 1;  //number of periodic boundary PAIRS (1(left)-2(right) for 2D)
                pbc_indices[0] = 2;  // define the second boundary in the pair
                pbc = new PeriodicBounds(pGrid, pbc_indices, pbc_count);
                delete[] pbc_indices;
            } else if (inputData.nsd == 3) {
                int *pbc_indices = new int[2];
                int pbc_count = 2;  //number of periodic boundary PAIRS (1(left)-2(right), 3(front)-4(back))
                pbc_indices[0] = 2;
                pbc_indices[1] = 6;
                pbc = new PeriodicBounds(pGrid, pbc_indices, pbc_count);
                delete[] pbc_indices;
            }
            PeriodicData periodic_dd(pbc, 3);
            periodic_dd.SetVarIndexPeriodic(0);  //set indicators to each periodic variable.
            periodic_dd.SetVarIndexPeriodic(1);
            periodic_dd.SetVarIndexPeriodic(2);
            PeriodicData periodic_exc(pbc, 1);
            periodic_exc.SetVarIndexPeriodic(0);

            // Construct gridfield based on Grid
            data.redimGrid(pGrid);

            //cluster mesh
            if (inputData.beta_1 > 0) {
                data.mesh_clustering(&inputData);
            }
            data.redimNodeData();

            // Read in a morphology from file and add transport layer
            std::size_t found = inputData.MSfile.find(".bits");
            if (found != std::string::npos) {
                PrintStatus("The input file is \'.txt\' format!");
                data.SetCoarseData_txt(inputData.MSfile);
            } else {
                PrintStatus("The input file is \'.plt\' format!");
                data.SetCoarseData(inputData.MSfile);
            }

            // Calculate the signed distance contour function using level set method.
            data.SetIndicators(inputData);
            auto *distance = new DistContour();
            int transportLayer = int(inputData.elmny * inputData.transportLayer);
            int extendNy = inputData.elmny + 1 + transportLayer * 2;
            if (inputData.nsd == 2) {
                // distance->setthreshold(inputData.threshold);
                distance->DistContour2d(data.morph2d, inputData.elmnx + 1, extendNy, 1, inputData.isPeriodic,
                                        inputData.Geo2Phy);
                data.interpolate(distance->dt, distance->dtdx, distance->dtdy);  // interpolation function
            } else if (inputData.nsd == 3) {
                // distance->setthreshold(inputData.threshold);
                distance->DistContour3d(data.morph3d, inputData.elmnx + 1, extendNy, inputData.elmnz + 1,
                                        inputData.isPeriodic, inputData.Geo2Phy);
                data.interpolate(distance->dt, distance->dtdx, distance->dtdy,
                                 distance->dtdz);  // interpolation function
            }
#ifdef NDEBUG
            save_gf(&data, &inputData, "debug_interpolated.plt", "Interpolated");
#endif
            if (rank == 0) distance->printMsg();
            InitialGuess init_(inputData.nsd, pGrid, &data); //initial guess pointer.
            data.SetInitialGuess(&init_);  // set initial guess
            // Set Initial Conditions
            if (!SetIC(data, inputData)) {
                PrintResults("Failed to set initial conditions", false, inputData.shouldFail);
                delete pGrid;
                throw (std::string("Problem with IC, not loaded"));
            }
            /*
            for (int nodeID = 0; nodeID < pGrid->n_nodes(); nodeID++) {
                 const ZEROPTV& p = pGrid->GetNode(nodeID)->location();
                 DDNodeData *pData = &(data.GetNodeData (nodeID));
                 if (p.y() > 0.99 && pData->dt > -1.) {
                    pGrid->GetNode(nodeID)->addIndicatorNum(change);
                 }
                 if (p.y() < 0.01 && pData->dt < 1.) {
                 pGrid->GetNode(nodeID)->addIndicatorNum(change);
                 }
            }
            pGrid->SetCaredSurfaceIndicator();
            pGrid->GenElmSurfaceIndicator();
             */
            // Set Solver parameters
            int nOfDofPerNode = inputData.couple; // number of degree of freedom per node
            int nOfDofPerNode_ksp = 1;  // number of degree of freedom per node for ksp
            int nOfGPs = inputData.ngp;  // number of degree of freedom per node
            int total_n_post = 15; // total number of column in post process vector; this number has been designed for both 2D and 3D
            int total_n_j = 25; // total number of column in post process vector that return J; this number has been designed for both 2D and 3D
            bool do_accelerate = (inputData.use_elemental_assembly_ == 1);
            AssemblyMethod assembly_method = kAssembleGaussPoints;
            if (inputData.use_elemental_assembly_ == 1) {
                assembly_method = kAssembleElements;
            }
            DDEquation ddEq(&inputData, &init_, do_accelerate, assembly_method); // the equation solver
            if (inputData.isPeriodic) {
                ddEq.redimSolver(pGrid, nOfDofPerNode, inputData.basisFunction, 0, &periodic_dd);
                ddEq.SetPreallocator();
                ddEq.PresetPeriodicData(1, 3, 3);
            } else {
                ddEq.redimSolver(pGrid, nOfDofPerNode, inputData.basisFunction,
                                 0); // Setting up solver for coupled DD Eq.
            }
            ddEq.setData(&data);
            DDExciton ddEx(&inputData, assembly_method); // the exciton equation solver
            if (inputData.isPeriodic) {
                ddEx.redimSolver(pGrid, nOfDofPerNode_ksp, inputData.basisFunction, 0, &periodic_exc);
                ddEx.SetPreallocator();
                ddEx.PresetPeriodicData(1, 1, 1);
            } else {
                ddEx.redimSolver(pGrid, nOfDofPerNode_ksp, inputData.basisFunction,
                                 0); // Setting up ksp solver for exciton
            }
            ddEx.setData(&data);
            lsqr ddls(&inputData,
                      assembly_method); // The least square solver. This function finds nodal values from G.P values
            ddls.redimSolver(pGrid, nOfGPs, inputData.basisFunction == BASIS_HERMITE,
                             0); // Setting up ksp solver for least square problem
            ddls.setData(&data);
            PrintStatus("solver initialized!", rank);
            // Class/Function to calculate dissociation and recombination rate
            gen_diss gendiss; // This function calculates kdiss and R_recombination
            gendiss.inputData = &inputData;
            ddEq.gendiss = &gendiss;
            ddEx.gendiss = &gendiss;
            data.gendiss = &gendiss;
            // The following line output the initial setup, it is only used for testing
            save_gf(&data, &inputData, "data_init.plt", "Initial Condition");
            // Loop to get the Jv curve;
            bool loopend = true;      //
            ddEq.Vapp = 0;
            ofstream Jvcurvefile;     // Create a file to save JV curve;
            char jvcurve[100];
            sprintf(jvcurve, "JVcurve_%s", inputData.MSfile.c_str());
            Jvcurvefile.open(jvcurve);
            Jvcurvefile << "Applied_voltage (eV)     Jpy     Jny    J(mA/cm^-2)" << endl;
            char resultfile[150];
            std::vector<float> nc;
            double error = 1.;        // L2 Error computed from results of current Exciton and previous exciton
            double converge_number = 1e-4;
            int index_continuation = 4;
            while (loopend){
                // This loop is for computing JV curve. You can compute JV curve by setting Voc in config file.

                if (rank == 0)
                    PrintStatus("=================Current applied voltage =", ddEq.Vapp, "eV=================");
                error = 0.;
                while (error > converge_number || index_continuation < 2) {
                    // This loop determines how many time Exciton equation and drift diffusion equation should be solve to converge. ||index_continuation<3
                    ZEROMATRIX<double> vector_old = data.Error(
                            &inputData);  // results of electron, hole density from previous exciton
                    ddEx.Solve();        // Solving Exciton Equation
                    ddEq.Solve();        // Solving drift diffusion Equation
                    ZEROMATRIX<double> vector_new = data.Error(
                            &inputData);  // results of electron, hole density from current exciton
                    /* L2 error computation
                       @ param error refers to summation of (n_new - n_old)
                       @ param norm refers to norm of n_old
                     */
                    double error_n = 0;
                    double norm_n = 0;
                    double error_p = 0;
                    double norm_p = 0;
                    double error_x = 0;
                    double norm_x = 0;
                    for (int i = 0; i < data.p_grid_->n_nodes(); i++) {
                        error_n +=
                                fabs(vector_new(i, 0) - vector_old(i, 0)) * fabs(vector_new(i, 0) - vector_old(i, 0));
                        norm_n += vector_old(i, 0) * vector_old(i, 0);
                        error_p +=
                                fabs(vector_new(i, 1) - vector_old(i, 1)) * fabs(vector_new(i, 1) - vector_old(i, 1));
                        norm_p += vector_old(i, 1) * vector_old(i, 1);
                        error_x +=
                                fabs(vector_new(i, 2) - vector_old(i, 2)) * fabs(vector_new(i, 2) - vector_old(i, 2));
                        norm_x += vector_new(i, 2) * vector_new(i, 2);
                    }
                    double totalnorm_n = 1;
                    double totalnorm_p = 1;
                    double totalnorm_x = 1;
                    MPI_Allreduce (&norm_n, &totalnorm_n, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
                    MPI_Allreduce (&norm_p, &totalnorm_p, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
                    MPI_Allreduce (&norm_x, &totalnorm_x, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
                    double totalerror_n = 1;
                    double totalerror_p = 1;
                    double totalerror_x = 1;
                    MPI_Allreduce (&error_n, &totalerror_n, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
                    MPI_Allreduce (&error_p, &totalerror_p, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
                    MPI_Allreduce (&error_x, &totalerror_x, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
                    totalnorm_n = sqrt(totalnorm_n);
                    error = sqrt(totalerror_n) / totalnorm_n;
                    double error_1 = sqrt(totalerror_p) / sqrt(totalnorm_p);
                    double error_2 = totalerror_x / sqrt(totalnorm_x);
                    PrintStatus(" condition error_n is = ", error, " condition error_p is = ", error_1,
                                " condition error_x is = ", error_2);
                    if (index_continuation < 1) {
                        data.increasePHIBAR();
                    }
                    index_continuation += 1;
                    PrintStatus("index_continuation", index_continuation);
                }

                // Check if test passed and print result
                PrintStatus("Checking solution...");
                // For usage of hermite basis we need to have post process
                if (inputData.basisFunction == BASIS_HERMITE ||
                    inputData.basisFunction == BASIS_LINEAR) {
                    data.CalcJ(data);
                    data.findGrad(data);
                    if (inputData.basisFunction == BASIS_HERMITE)
                        data.post_process(&inputData);
                } else {
                    /*
                     Looping Over all of post process variables obtained in G.Ps.
                     This variables of the solotion vector are:
                     Jny,Jnx,Jnz,Jpy,Jpx,Jpz;
                     phi_x; phi_y; n_x; n_y; p_x;p_y;
                     grad_jn; grad_jp; grad_jc
                     in the following loop we just do least square method for J.
                     if one need to have access to all of derivative can change/ comment
                     the if statement.
                     */
                    for (int kk = 0; kk < 4; kk++) {
                        if (kk < total_n_j) {
                            ddls.J_n = data.CalcJQuadCube(data);
                            ddls.kk = kk;
                            ddls.Solve(kk + 0.01);
                        }
                    }
                }

                data.checkJ(&inputData);

                if (inputData.ifPrintPltFiles) {
                    data.Calc_DR();
                    data.dimlize();
                    sprintf(resultfile, "result_%4.2f_%s%s", ddEq.Vapp, inputData.MSfile.c_str(),
                            inputData.output_extension.c_str());
                    save_gf(&data, &inputData, resultfile, 0.0);
                    data.nondimlize();
                }
                /*
                 * We have two functions that calculate Jc in boundaries.
                 * The first method has been prefered.
                 * One can use the 2nd one with un-commenting that.
                 */
                // Finding integrate of Jc (collected J in boundaries)
                std::vector<float> J = data.integrateJ(data);
                double maxerrorcheck = data.FindMaxGradJ(data);
                Jvcurvefile << ddEq.Vapp << "          " << J.at(0) << "          " << J.at(1) << "          "
                            << J.at(2)
                            << "          " << maxerrorcheck << endl;
                // We are doing this integration only if vapp = 0;
                // These are some information to check if our result is reasonable.
                if (ddEq.Vapp < 0.00001) {
                    std::vector<double> DRG = data.domain_integrations(data);
                    ofstream DRGfile;
                    DRGfile.open("drg.txt");
                    double Height = inputData.Height;
                    double conversion_factor = 1. * Height * Height * (0.1 / (inputData.L[0] * Height)) * inputData.q;
                    DRGfile << DRG.at(0) * conversion_factor << " " << DRG.at(1) * conversion_factor << " "
                            << DRG.at(2) * conversion_factor << endl;
                }

                //Increase applied voltage by 0.1eV;
                if (ddEq.Vapp < 2 * (inputData.Voc) / 3)
                    ddEq.Vapp += (inputData.vapp_step) * inputData.vapp_ratio;
                else
                    ddEq.Vapp += inputData.vapp_step;
                //Terminate the calculation at certain Vapp.
                if (ddEq.Vapp > inputData.Voc)
                    loopend = false;
            }
            Jvcurvefile.close();
        }
        PrintStatus("End of program");
        DestroyGrid(pGrid);
    }
    catch (const std::string &s) {
        PetscFinalize();
        std::cerr << s << std::endl;
        return -1;
    }
    catch (std::bad_alloc &e) {
        PetscFinalize();
        std::cerr << "Problem with memory allocation " << e.what();
        return -1;
    }
    catch (const TALYException &e) {
        e.print();
        PetscFinalize();
        return -1;
    }
    PetscFinalize();
    return 0;
}

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
#ifndef BT_GRIDFIELD_HPP
#define BT_GRIDFIELD_HPP
#include "../include/linearInterpolation.hxx"
#include "../include/interpolator.h"
#include "../include/gen_diss.hh"
#include <iostream>
#include <fstream>
#include "DDNodeData.h"
#include "DDInputData.h"
#include "initial_guess.h"
#include <vector>
using std::vector;
typedef vector < float >vec;
typedef vector < vec > Mat2d;
typedef vector < vector < vec >> Mat3d;

using namespace std;
using namespace linearinterp;
class DDGridField:public GridField < DDNodeData > {
public:
    gen_diss *gendiss;           ///< context for calculating dissociation and recombination
    DDInputData *input_;
    Mat2d morph2d;                //2d array/matrix storing the morphology;
    Mat3d morph3d;                //3d ....................................
    Mat2d DT2d;                   //2d array/matrix storing the distance contour;
    Mat2d DT2ddx;                 //2d array/matrix storing the x derivative of distance contour
    Mat2d DT2ddy;                 //............................y...............................
    Mat3d DT3ddx, DT3ddy, DT3ddz;
    Mat3d DT3d;                   //3d ....................................
    InitialGuess *initial_guess_;

    DDGridField(DDInputData *input) : input_(input), initial_guess_(NULL) {
    }

    virtual ~ DDGridField() {
    }

    void increasePHIBAR() {
        for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
            DDNodeData *pData = &(GetNodeData(nodeID));
            if (pData->PHIBAR != 0 || pData->D2PHIBAR != 0) {
                pData->D2PHIBAR *= 22.;
                pData->PHIBAR *= 22.;
            }
        }
    }

    void SetInitialGuess(InitialGuess *ig) {
        initial_guess_ = ig;
    }

    /*
     void increaseInterfacialPotential(double times) {
        for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++)
            DDNodeData *pData = &(GetNodeData(nodeID));
        pData->D2PHIBAR *= times;
        pData->PHIBAR *= times;
    }
    */


    void SetIC(int nsd) {

        double mu_n = input_->mu_n;
        double mu_p = input_->mu_p;
        double coeff = input_->muRatio;
        double DT_D = input_->interfaceThk / 2., DT_A = input_->interfaceThk / 2;
        double y_h = 0;
        for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
            DDNodeData *pData = &(GetNodeData(nodeID));
            y_h = p_grid_->GetCoord(nodeID, 1);
            for (int i = 0; i < 8; i++) {
                pData->phi[i] = 0;      // initial value of phi

                pData->n[i] = 0;        // initial value of n

                pData->p[i] = 0;        // initial value of p

                pData->x[i] = 0;        // initial value of x
            }
            pData->phi[0] = initial_guess_->ValueAt(nodeID, 0);  // Guess for phi (linear)
            pData->n[0] = initial_guess_->ValueAt(nodeID, 1);  // Guess for n (exponential)
            pData->p[0] = initial_guess_->ValueAt(nodeID, 2);  // Guess for p (exponential)
            double DT = pData->dt;
            double MS = pData->ms;
            if ((DT < DT_D) && (DT > -DT_A)) {
                pData->epsn = 0.5 * (input_->eps_A + input_->eps_D);
                pData->mun = (mu_n) * 0.5;
                pData->mup = (mu_p) * 0.5;
                if (y_h > 0.92)
                    pData->mun *= 0.001;
                if (y_h < 0.08)
                    pData->mup *= 0.001;
            } else {
                pData->mun = ((DT >= DT_D) ? mu_n : coeff * mu_n);
                pData->mup = ((DT >= DT_A) ? coeff * mu_p : mu_p);
                pData->epsn = DT > 0 ? input_->eps_A : input_->eps_D;
            }
            // Set the second derivative of interfacial potential = -1/3.14*arctan(DT)
            double temp = 1. + DT * DT;
            if (DT > 2 || DT < -2) {
                pData->D2PHIBAR = 0;
            } else {
                pData->D2PHIBAR = 2. * DT / temp / temp * (1. / 3.14);
            }
            // multiplied by a factor, which means currently the interfacial potential is part of real one.
            pData->D2PHIBAR *= input_->potentialFactor * 1e18 / (input_->q * input_->C0());

            // set the first derivative of interfacial potential = -1/3.14*arctan(DT)
            if (DT > 2) {
                //pData->PHIBAR = -1./3.14*atan(2.)/input_->V_t();
                pData->PHIBAR = 0;  //away from the inerface, the gradient of interfacial potential is close to 0.
            } else if (DT < -2) {
                //pData->PHIBAR = -1./3.14*atan(-2.)/input_->V_t();
                pData->PHIBAR = 0;
            } else {
                //pData->PHIBAR = -1./3.14*atan(DT)/input_->V_t();
                pData->PHIBAR = -1. / 3.14 / temp * input_->potentialFactor;
            }
            //set exciton generation curves
            if (MS < 0.5 && DT > -100) {
                double heightinnm = input_->realHeight * 1.e9;
                double relativeHeight = y_h * input_->x0() * 1.e9 -
                                        input_->heightofTransportLayer; //distance away from the bottom of true active layer in other word
                //relative to the top surface of bottom added material/layer.
                if (input_->G > 1) {    //constant generation
                    pData->GenX = input_->G;
                } else if (heightinnm < 140 && heightinnm > 130) { //in this range, use curver for thickness=135nm

                } else if (heightinnm <= 130 && heightinnm > 120) { //use curve for thickness=125nm
                    double x2 = relativeHeight * relativeHeight;
                    double x3 = x2 * relativeHeight;
                    double x4 = x2 * x2;
                    double x5 = x2 * x3;
                    double x6 = x3 * x3;
                    pData->GenX =
                            -2.838348315e10 * x6 + 7.383988429e12 * x5 - 1.514751248e14 * x4 - 7.492291583e16 * x3 +
                            3.334406435e18 * x2 + 1.31360722e20 * relativeHeight + 1.48625748e21;
                    pData->GenX *= 2 * 1e6;
                } else if (heightinnm <= 120 && heightinnm > 110) { //use curve for thickness=115nm
                    double x2 = relativeHeight * relativeHeight;
                    double x3 = x2 * relativeHeight;
                    double x4 = x2 * x2;
                    double x5 = x2 * x3;
                    double x6 = x3 * x3;
                    pData->GenX =
                            -3.638749855e10 * x6 + 9.965731197e12 * x5 - 4.368840167e14 * x4 - 6.283771826e16 * x3 +
                            3.099825813e18 * x2 + 1.42636502e20 * relativeHeight + 1.555488201e21;
                    pData->GenX *= 2 * 1e6;
                } else if (heightinnm <= 110 && heightinnm >= 100) { //use curve for thickness=105nm
                    double x2 = relativeHeight * relativeHeight;
                    double x3 = x2 * relativeHeight;
                    double x4 = x2 * x2;
                    double x5 = x2 * x3;
                    double x6 = x3 * x3;
                    pData->GenX =
                            -4.701565876e10 * x6 + 1.289735413e13 * x5 - 6.870417035e14 * x4 - 5.643274827e16 * x3 +
                            2.871177515e18 * x2 + 1.644511982e20 * relativeHeight + 1.848587415e21;
                    pData->GenX *= 2 * 1e6;
                } else if (heightinnm <= 300 && heightinnm >= 240) {
                    double x = relativeHeight;
                    double x2 = x * x;
                    double x3 = x2 * relativeHeight;
                    double x4 = x2 * x2;
                    double x5 = x2 * x3;
                    double x6 = x3 * x3;
                    double x7 = x3 * x4;
                    double x8 = x4 * x4;
                    pData->GenX =
                            -209897.2991 * x8 + 232740851.5 * x7 - 1.010498938 * 1e11 * x6 + 2.158052245 * 1e13 * x5 -
                            2.322604418 * 1e15 * x4 + 1.145233677 * 1e17 * x3 - 2.049681128 * 1e18 * x2 -
                            1.625666414 * 1e18 * x + 5.618839707 * 1e21;
                    pData->GenX *= 2 * 1e6;
                } else if (heightinnm <= 550 && heightinnm >= 480) {
                    double x = relativeHeight;
                    double x2 = x * x;
                    double x3 = x2 * relativeHeight;
                    double x4 = x2 * x2;
                    double x5 = x2 * x3;
                    double x6 = x3 * x3;
                    double x7 = x3 * x4;
                    double x8 = x4 * x4;
                    double x9 = x4 * x5;
                    pData->GenX = 7.912187168 * x9 - 18120.2018 * x8 + 17269971.83 * x7 - 8877985912 * x6 +
                                  2.667433909 * 1e12 * x5 - 4.741357483 * 1e14 * x4 + 4.786531782 * 1e16 * x3 -
                                  2.365249186 * 1e18 * x2 + 7.11392015 * 1e18 * x + 6.299113872 * 1e21;
                    pData->GenX *= 2 * 1e6;
                } else if (heightinnm <= 1100 && heightinnm >= 950) {
                    double x = relativeHeight;
                    double x2 = x * x;
                    double x3 = x2 * relativeHeight;
                    double x4 = x2 * x2;
                    double x5 = x2 * x3;
                    double x6 = x3 * x3;
                    double x7 = x3 * x4;
                    double x8 = x4 * x4;
                    double x9 = x4 * x5;
                    double x10 = x5 * x5;
                    pData->GenX =
                            -1.640795468 * 1e-5 * x10 + 7.698748822 * 1e-2 * x9 - 152.9974546 * x8 + 167583.6276 * x7 -
                            110113971.6 * x6 + 4.404467363 * 1e10 * x5 - 1.017112188 * 1e13 * x4 +
                            1.019111146 * 1e15 * x3 + 9.657183117 * 1e16 * x2 - 4.763205275 * 1e19 * x +
                            6.542031817 * 1e21;
                    pData->GenX *= 2 * 1e6;
                }

            }
            if (MS > 1) {  //This is for 3phase model with the mixture phase being MS==3;
                pData->GenX = input_->G * input_->D_fraction;
            }


        }

        PrintStatusStream(std::cerr, "IC set ");
    }

    /**
     * Calculate the current density for Hermite Basis Functions and linear basis function. if Hermite is the basis
     * calculated J is not the final result because the existance of bad values near to cathode and anode.
     * This values happen near to first and last 10% of domain. Because we are using directly values of
     * gradiant in y direction and There is no B.C for them so it causes bad results near to cathode and anode.
     * To remove this Effect later you need to call postprocess function.
     * If we are using linear Basis function calculated J is the final result,
     * Because we are calculating J in G.Ps then map it to nodes.
     */

    void CalcJ(DDGridField &data) {
        if (input_->basisFunction == BASIS_HERMITE) {
            double hx = input_->Height / input_->Nelem[0];
            double hy = (input_->nsd >= 2) ? input_->Height / input_->Nelem[1] : 0.0;
            double q = input_->q;
            double Vt = input_->Vt;
            for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
                DDNodeData *pData = &(GetNodeData(nodeID));
                double mu_n = pData->mun;
                double mu_p = pData->mup;
                pData->j[0] = mu_n * q * Vt / hy * (-pData->n[0] * pData->phi[2] + pData->n[2]);
                pData->j[0] = input_->c_unhat(pData->j[0]) / 10.;  //Divided by 10, so the the unit is : mA/cm^-2
                pData->j[1] = mu_n * q * Vt / hx * (-pData->n[0] * pData->phi[1] + pData->n[1]);
                pData->j[1] = input_->c_unhat(pData->j[1]) / 10.;
                pData->j[2] = mu_p * q * Vt / hy * (-pData->p[0] * pData->phi[2] - pData->p[2]);
                pData->j[2] = input_->c_unhat(pData->j[2]) / 10.;
                pData->j[3] = mu_p * q * Vt / hx * (-pData->p[0] * pData->phi[1] - pData->p[1]);
                pData->j[3] = input_->c_unhat(pData->j[3]) / 10.;
            }
            for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
                DDNodeData *pData = &(GetNodeData(nodeID));
                double DT = pData->dt;
                pData->j[4] = pData->j[0];
                if (DT < -1)
                    pData->j[4] = 0;
            }
            for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
                DDNodeData *pData = &(GetNodeData(nodeID));
                double DT = pData->dt;
                pData->j[5] = pData->j[2];
                if (DT > 1)
                    pData->j[5] = 0;
            }
        } else if (input_->basisFunction == BASIS_LINEAR) {
            ofstream file;
            FEMElm fe(p_grid_, BASIS_FIRST_DERIVATIVE | BASIS_POSITION);
            FEMElm fe_linear(p_grid_, BASIS_ALL);
            ZEROPTV Jn_x;               // Value of Jn_x in G.Ps
            ZEROPTV Jn_y;               // Value of Jn_y in G.Ps
            ZEROPTV Jn_z;
            ZEROPTV Jp_x;               // Value of Jp_x in G.Ps
            ZEROPTV Jp_y;               // Value of Jp_y in G.Ps
            ZEROPTV Jp_z;
            ZEROPTV phi_y;               // Value of Jp_x in G.Ps    ZEROPTV Jp_x;               // Value of Jp_x in G.Ps
            ZEROPTV phi_x;               // Value of Jp_x in G.Ps
            ZEROPTV phi_z;
            double phi_1 = 0;
            double phi_2 = 0;
            double phi_3 = 0;
            double n_0 = 0;
            double n_1 = 0;
            double n_2 = 0;
            double n_3 = 0;
            double p_0 = 0;
            double p_1 = 0;
            double p_2 = 0;
            double p_3 = 0;
            double q = input_->q;
            double Vt = input_->Vt;
            const double n_elements = p_grid_->n_elements();
            std::vector<double> elemental_Jny(n_elements);
            std::vector<double> elemental_Jnx(n_elements);
            std::vector<double> elemental_Jpy(n_elements);
            std::vector<double> elemental_Jpx(n_elements);
            std::vector<double> elemental_phix(n_elements);
            std::vector<double> elemental_phiy(n_elements);
            std::vector<double> elemental_phiz(n_elements);
            std::vector<double> elemental_Jnz(n_elements);
            std::vector<double> elemental_Jpz(n_elements);
            int nx = (input_->Nelem[0] + 1);
            // Looping over elements
            for (int elm_ID = 0; elm_ID < n_elements; elm_ID++) {
                fe.refill(elm_ID, BASIS_LINEAR, -1);
                Jn_y(0) = 0;
                Jp_y(0) = 0;
                Jn_x(0) = 0;
                Jp_x(0) = 0;
                Jn_y(1) = 0;
                Jn_z(0) = 0;
                n_0 = 0;
                int number = 0;
                // Loping over G.ps
                while (fe.next_itg_pt()) {
                    ZEROPTV pt = fe.position();
                    double y_1 = pt(1);
                    double DT = this->initial_guess_->valueFEM(fe, DT_ID, N_COUPLED);
                    double mu_n = this->initial_guess_->valueFEM(fe, MUN_ID, N_COUPLED);
                    double mu_p = this->initial_guess_->valueFEM(fe, MUP_ID, N_COUPLED);

                    phi_1 = this->initial_guess_->valueDerivativeFEM(fe, PHI_ID, 0, N_COUPLED);
                    phi_2 = this->initial_guess_->valueDerivativeFEM(fe, PHI_ID, 1, N_COUPLED);
                    phi_3 = (input_->nsd == 3) ? this->initial_guess_->valueDerivativeFEM(fe, PHI_ID, 2, N_COUPLED) : 0;
                    n_1 = this->initial_guess_->valueDerivativeFEM(fe, N_ID, 0, N_COUPLED);
                    n_2 = this->initial_guess_->valueDerivativeFEM(fe, N_ID, 1, N_COUPLED);
                    n_3 = (input_->nsd == 3) ? this->initial_guess_->valueDerivativeFEM(fe, N_ID, 2, N_COUPLED) : 0;
                    p_1 = this->initial_guess_->valueDerivativeFEM(fe, P_ID, 0, N_COUPLED);
                    p_2 = this->initial_guess_->valueDerivativeFEM(fe, P_ID, 1, N_COUPLED);
                    p_3 = (input_->nsd == 3) ? this->initial_guess_->valueDerivativeFEM(fe, P_ID, 2, N_COUPLED) : 0;

                    // Average value of n in G.Ps
                    n_0 = this->initial_guess_->valueFEM(fe, N_ID, N_COUPLED);
                    p_0 = this->initial_guess_->valueFEM(fe, P_ID, N_COUPLED);
                    Jn_y(0) += mu_n * q * Vt * (-n_0 * phi_2 + n_2);
                    Jn_x(0) += mu_n * q * Vt * (-n_0 * phi_1 + n_1);
                    Jp_y(0) += mu_p * q * Vt * (-p_0 * phi_2 - p_2);
                    Jp_x(0) += mu_p * q * Vt * (-p_0 * phi_1 - p_1);
                    if (input_->nsd == 3) {
                        Jn_z(0) += mu_n * q * Vt * (-p_0 * phi_3 + n_3);
                        Jp_z(0) += mu_p * q * Vt * (-p_0 * phi_3 - p_3);
                    }
                    number++;
                }
                Jn_y(0) = input_->c_unhat(Jn_y(0)) / input_->x0();  //Changing J to dimensional form
                Jn_x(0) = input_->c_unhat(Jn_x(0)) / input_->x0();  //Changing J to dimensional form
                Jp_y(0) = input_->c_unhat(Jp_y(0)) / input_->x0();  //Changing J to dimensional form
                Jp_x(0) = input_->c_unhat(Jp_x(0)) / input_->x0();  //Changing J to dimensional form

                phi_y(0) = phi_2;
                phi_x(0) = phi_1;
                if (input_->nsd == 3) {
                    Jn_z(0) = input_->c_unhat(Jn_z(0)) / input_->x0();
                    Jp_z(0) = input_->c_unhat(Jp_z(0)) / input_->x0();
                    phi_z(0) = phi_3;
                }
                elemental_Jny.at(elm_ID) = 0.1 * Jn_y(0) / number;  //Divided by 10, so the the unit is : mA/cm^-2
                elemental_Jnx.at(elm_ID) = 0.1 * Jn_x(0) / number;  //Divided by 10, so the the unit is : mA/cm^-2
                elemental_Jpy.at(elm_ID) = 0.1 * Jp_y(0) / number;  //Divided by 10, so the the unit is : mA/cm^-2
                elemental_Jpx.at(elm_ID) = 0.1 * Jp_x(0) / number;  //Divided by 10, so the the unit is : mA/cm^-2
                elemental_phiy.at(elm_ID) = phi_y(0) / number;  //Divided by 10, so the the unit is : mA/cm^-2
                elemental_phix.at(elm_ID) = phi_x(0) / number;  //Divided by 10, so the the unit is : mA/cm^-2
                if (input_->nsd == 3) {
                    elemental_Jnz.at(elm_ID) = 0.1 * Jn_z(0) / number;
                    elemental_Jpz.at(elm_ID) = 0.1 * Jp_z(0) / number;
                    elemental_phiz.at(elm_ID) = phi_z(0) / number;
                }
            }

            // build a map that says which elements each node belongs to
            // node_to_elems[0] -> the list of elements that node 0 belongs to
            std::vector<std::vector<int >> node_to_elems(p_grid_->n_nodes());
            for (int elm_ID = 0; elm_ID < n_elements; elm_ID++) {
                ELEM *elem = p_grid_->GetElm(elm_ID);
                for (int lcl_node = 0; lcl_node < elem->n_nodes(); lcl_node++) {
                    int node_ID = elem->ElemToLocalNodeID(lcl_node);
                    node_to_elems.at(node_ID).push_back(elm_ID);
                }
            }
            for (int node_ID = 0; node_ID < p_grid_->n_nodes(); node_ID++) {
                // calculate Jn averaged across all the elements containing this node
                double val = 0.0;
                double val_1 = 0.0;
                double val_2 = 0.0;
                double val_3 = 0.0;
                double val_4 = 0.0;
                double val_5 = 0.0;
                double val_6 = 0.0;
                std::vector<int> &elms = node_to_elems.at(node_ID);
                for (int e = 0; e < elms.size(); e++) {
                    int elm_ID = elms.at(e);
                    val += elemental_Jny.at(elm_ID);
                }
                val /= elms.size();
                std::vector<int> &elms_1 = node_to_elems.at(node_ID);
                for (int e = 0; e < elms_1.size(); e++) {
                    int elm_ID = elms_1.at(e);
                    val_1 += elemental_Jpy.at(elm_ID);
                }
                val_1 /= elms_1.size();
                std::vector<int> &elms_2 = node_to_elems.at(node_ID);
                for (int e = 0; e < elms_2.size(); e++) {
                    int elm_ID = elms_2.at(e);
                    val_2 += elemental_Jnx.at(elm_ID);
                }
                val_2 /= elms_2.size();
                std::vector<int> &elms_3 = node_to_elems.at(node_ID);
                for (int e = 0; e < elms_3.size(); e++) {
                    int elm_ID = elms_3.at(e);
                    val_3 += elemental_Jpx.at(elm_ID);
                }
                val_3 /= elms.size();
                std::vector<int> &elms_4 = node_to_elems.at(node_ID);
                for (int e = 0; e < elms_4.size(); e++) {
                    int elm_ID = elms_4.at(e);
                    val_4 += elemental_phix.at(elm_ID);
                }
                val_4 /= elms.size();
                std::vector<int> &elms_5 = node_to_elems.at(node_ID);
                for (int e = 0; e < elms_5.size(); e++) {
                    int elm_ID = elms_5.at(e);
                    val_5 += elemental_phiy.at(elm_ID);
                }
                val_5 /= elms.size();
                std::vector<int> &elms_6 = node_to_elems.at(node_ID);
                for (int e = 0; e < elms_6.size(); e++) {
                    int elm_ID = elms_6.at(e);
                    val_6 += elemental_phiz.at(elm_ID);
                }
                val_6 /= elms.size();

                // save on DDNodeData
                DDNodeData &pData = GetNodeData(node_ID);
                pData.j[0] = val;   //j_0: jny
                pData.j[1] = val_2; //j_1: jnx
                pData.j[2] = val_1; //j_2: jpy
                pData.j[3] = val_3; //j_3: jpx
                pData.j[6] = val_4; //j6: dphi_dx
                pData.j[7] = val_5; //j7: dphi_dy
                pData.j[5] = val_6; //j5: dphi_dz
            }
            for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
                DDNodeData &pData = GetNodeData(nodeID);
                pData.j[0] = fabs(pData.j[0]);
                pData.j[2] = fabs(pData.j[2]);
                pData.j[1] = fabs(pData.j[1]);
                pData.j[3] = fabs(pData.j[3]);
                pData.j[4] = pData.j[0];
                pData.j[5] = pData.j[2];
            }
            /*for (int nodeID = 0; nodeID < p_grid_->n_nodes (); nodeID++)
    {
      DDNodeData & pData = GetNodeData (nodeID);
      double y_1 = p_grid_->GetCoord (nodeID, 1);
      if (y_1 >= 1.-input_->transportLayer/(1.+2*input_->transportLayer)-0.01 && y_1 <= 1.-input_->transportLayer/(1.+2*input_->transportLayer)+0.01)
        {
        DDNodeData & pData1 = GetNodeData (nodeID - 1 * nx);  // node data of down node
        DDNodeData & pData2 = GetNodeData (nodeID - 2 * nx);  // node data of 2nd down node
        DDNodeData & pData3 = GetNodeData (nodeID - 3 * nx);  // node data of 3rd down node
        DDNodeData & pData4 = GetNodeData (nodeID - 4 * nx);  // node data of 4th down node
        double Jpy_1 = pData1.j[2];  //Jp of down node
        double Jpy_2 = pData2.j[2];  //Jp of 2nd down node
        double Jpy_3 = pData3.j[2];  //Jp of 3rd down node
        double Jpy_4 = pData4.j[2];  //Jp of 4th down node
        double Jpy_i = 12 / 25. * (4. * Jpy_1 - 3. * Jpy_2 + 4 / 3. * Jpy_3 - 1 / 4. * Jpy_4);  // 4th order forward difference
        pData.j[2] = fabs(Jpy_i) + 1/600. * fabs (pData.j[3]);
        }
    }
    for (int nodeID = p_grid_->n_nodes () - 1; nodeID >= 0; nodeID--)
    {
      DDNodeData & pData = GetNodeData (nodeID);
      double y_1 = p_grid_->GetCoord (nodeID, 1);
      if (y_1 <= input_->transportLayer/(1.+2*input_->transportLayer)+0.01 && y_1 >= input_->transportLayer/(1.+2*input_->transportLayer)-0.01)
        {
        DDNodeData & pData1 = GetNodeData (nodeID + 1 * nx);  // node data of down node
        DDNodeData & pData2 = GetNodeData (nodeID + 2 * nx);  // node data of 2nd down node
        DDNodeData & pData3 = GetNodeData (nodeID + 3 * nx);  // node data of 3rd down node
        DDNodeData & pData4 = GetNodeData (nodeID + 4 * nx);  // node data of 4th down node
        double Jny_1 = pData1.j[0];  //Jny of above node
        double Jny_2 = pData2.j[0];  //Jny of 2nd above node
        double Jny_3 = pData3.j[0];  //Jny of 3rd above node
        double Jny_4 = pData4.j[0];  //Jny of 4th above node
        double Jny_i = 12 / 25. * (4. * Jny_1 - 3. * Jny_2 + 4 / 3. * Jny_3 - 1 / 4. * Jny_4);  // 4th order backward difference
        pData.j[0] = fabs(Jny_i) + 1/600. * fabs (pData.j[1]);
        }
    }*/
        }
    }

    /**
     * NOTE: This function is used for calculatation of J and all other gradients of the problem.
     * (We will use this function for Quad-cubic Basis functions). After calling this function
     * in main we need to call "lsq.h" function to map variables to nodes
     */
    std::vector<std::vector<ElemNodeID> > small_elem_connectivities(ElemType type, kBasisFunction bf) {
        if (type == kElem2dBox && bf == BASIS_CUBIC) {
            return {
                    {0,  4,  12, 11},
                    {4,  5,  13, 12},
                    {5,  1,  6,  13},
                    {11, 12, 15, 10},
                    {12, 13, 14, 15},
                    {13, 6,  7,  14},
                    {10, 15, 9,  3},
                    {15, 14, 8,  9},
                    {14, 7,  2,  8},
            };
        } else if (type == kElem2dBox && bf == BASIS_QUADRATIC) {
            return {
                    {0, 4, 8, 7},
                    {4, 1, 5, 8},
                    {7, 8, 6, 3},
                    {8, 5, 2, 6},
            };
        } else if (type == kElem2dBox && bf == BASIS_LINEAR) {
            return {
                    {0, 1, 2, 3},
            };
        }

        throw NotImplementedException() << "Small elements not added for type '" << type << "'";
    }

    ZEROMATRIX<double> CalcJQuadCube(DDGridField &data) {

        FEMElm fe(p_grid_, BASIS_FIRST_DERIVATIVE | BASIS_POSITION);
        FEMElm fe_linear(p_grid_, BASIS_ALL);
        ZEROPTV Jn_y;               // Declaration of  Value of Jn_y in G.Ps
        ZEROPTV Jp_y;               // Declaration of  Value of Jp_y in G.Ps
        ZEROPTV Jn_x;               // Declaration of  Value of Jn_x in G.Ps
        ZEROPTV Jp_x;               // Declaration of  Value of Jp_x in G.Ps
        ZEROPTV phi_y;               // Declaration of  Value of Jp_x in G.Ps
        ZEROPTV phi_x;               // Declaration of  Value of Jp_x in G.Ps
        ZEROPTV n_y;               // Declaration of  Value of Jp_x in G.Ps
        ZEROPTV n_x;               // Declaration of  Value of Jp_x in G.Ps
        ZEROPTV p_y;               // Declaration of  Value of Jp_x in G.Ps
        ZEROPTV p_x;               // Declaration of  Value of Jp_x in G.Ps
        double phi_2 = 0; // Declaration of Derivative of phi respecto y
        double n_2 = 0; //  Declaration of Derivative of n respecto y
        double phi_1 = 0; // Declaration of  Derivative of phi respecto x
        double n_1 = 0; // Declaration of  Derivative of n respecto x
        double p_1 = 0; // Declaration of  Derivative of p respecto x
        double p_2 = 0; // Declaration of  Derivative of p respecto y
        double q = input_->q; // 1.6e-19
        double Vt = input_->Vt; //0.0259
        double mu_n = 0;
        double mu_p = 0;
        int n_var = 12; // number of variables this function should return
        const double n_elements = p_grid_->n_elements();
        std::vector<double> elemental_Jny(n_elements); // Storage vector for jny
        std::vector<double> elemental_Jnx(n_elements); // Storage vector for jnx
        std::vector<double> elemental_Jpy(n_elements); // Storage vector for jpy
        std::vector<double> elemental_Jpx(n_elements); // Storage vector for jpx
        std::vector<double> elemental_phix(n_elements); // Storage vector for phi_x
        std::vector<double> elemental_phiy(n_elements); // Storage vector for phi_y
        std::vector<double> elemental_nx(n_elements); // Storage vector for nx
        std::vector<double> elemental_ny(n_elements); // Storage vector for ny
        std::vector<double> elemental_px(n_elements); // Storage vector for px
        std::vector<double> elemental_py(n_elements); // Storage vector for py
        ZEROMATRIX<double> J;
        J.redim(n_elements, n_var);
        double y_1;
        // Looping over elements
        ofstream file;
        for (int elm_ID = 0; elm_ID < n_elements; elm_ID++) {
            Jn_y(0) = 0;
            Jp_y(0) = 0;
            Jn_x(0) = 0;
            Jp_x(0) = 0;
            Jn_y(1) = 0;
            //  n_0 = 0;
            int max_super = 0;
            if (input_->basisFunction == BASIS_QUADRATIC) {
                fe.refill(elm_ID, BASIS_QUADRATIC, -1);
                max_super = 4;
            } else if (input_->basisFunction == BASIS_CUBIC) {
                fe.refill(elm_ID, BASIS_CUBIC, 1);
                max_super = 16;
            }
            int number = 0;
            double node_left = -sqrt(35 + 8 * sqrt(7)) / (sqrt(3) * 5);
            double node_mid_left = -sqrt(35 - 8 * sqrt(7)) / (sqrt(3) * 5);
            double node_center = 0;
            double node_mid_right = sqrt(35 - 8 * sqrt(7)) / (sqrt(3) * 5);
            double node_right = sqrt(35 + 8 * sqrt(7)) / (sqrt(3) * 5);
            std::vector<ZEROPTV> points_cube = {
                    ZEROPTV(node_left, node_left, 0),
                    ZEROPTV(node_mid_left, node_left, 0),
                    ZEROPTV(node_center, node_left, 0),
                    ZEROPTV(node_mid_right, node_left, 0),
                    ZEROPTV(node_right, node_left, 0),

                    ZEROPTV(node_left, node_mid_left, 0),
                    ZEROPTV(node_mid_left, node_mid_left, 0),
                    ZEROPTV(node_center, node_mid_left, 0),
                    ZEROPTV(node_mid_right, node_mid_left, 0),
                    ZEROPTV(node_right, node_mid_left, 0),

                    ZEROPTV(node_left, node_mid_right, 0),
                    ZEROPTV(node_mid_left, node_mid_right, 0),
                    ZEROPTV(node_center, node_mid_right, 0),
                    ZEROPTV(node_mid_right, node_mid_right, 0),
                    ZEROPTV(node_right, node_mid_right, 0),

                    ZEROPTV(node_left, node_right, 0),
                    ZEROPTV(node_mid_left, node_right, 0),
                    ZEROPTV(node_center, node_right, 0),
                    ZEROPTV(node_mid_right, node_right, 0),
                    ZEROPTV(node_right, node_right, 0),

                    ZEROPTV(node_left, node_center, 0),
                    ZEROPTV(node_mid_left, node_center, 0),
                    ZEROPTV(node_center, node_center, 0),
                    ZEROPTV(node_mid_right, node_center, 0),
                    ZEROPTV(node_right, node_center, 0),
            };
            std::vector<ZEROPTV> points_quad = {
                    ZEROPTV(-1.0 / 1.7320508075688772935, -1.0 / 1.7320508075688772935, 0),
                    ZEROPTV(1.0 / 1.7320508075688772935, -1.0 / 1.7320508075688772935, 0),
                    ZEROPTV(-1.0 / 1.7320508075688772935, 1.0 / 1.7320508075688772935, 0),
                    ZEROPTV(1.0 / 1.7320508075688772935, 1.0 / 1.7320508075688772935, 0),
            };

            std::vector<ZEROPTV> points;
            if (input_->basisFunction == BASIS_QUADRATIC)
                points = points_quad;
            if (input_->basisFunction == BASIS_CUBIC)
                points = points_cube;
            // Loping over G.ps
            for (int i = 0; i < points.size(); i++) {
                fe.calc_at(points[i]);

                mu_n = this->initial_guess_->valueFEM(fe, MUN_ID, N_COUPLED); // mu_n value
                mu_p = this->initial_guess_->valueFEM(fe, MUP_ID, N_COUPLED); // mu_p value
                if (mu_n < mu_p) // Check if ratio of mun/mup is correct
                    mu_n *= 0.1;
                if (mu_p < mu_n) // Check if ratio of mup/mun is correct
                    mu_p *= 0.1;
                // Finding parameters to calculate Jn and Jp at gauss point
                phi_2 = this->initial_guess_->valueDerivativeFEM(fe, PHI_ID, 1,
                                                                 N_COUPLED); // Derivative of phi respecto y
                n_2 = this->initial_guess_->valueDerivativeFEM(fe, N_ID, 1, N_COUPLED); // Derivative of n respecto y
                phi_1 = this->initial_guess_->valueDerivativeFEM(fe, PHI_ID, 0,
                                                                 N_COUPLED); // Derivative of phi respecto x
                n_1 = this->initial_guess_->valueDerivativeFEM(fe, N_ID, 0, N_COUPLED); // Derivative of n respecto x
                p_1 = this->initial_guess_->valueDerivativeFEM(fe, P_ID, 0, N_COUPLED); // Derivative of p respecto x
                p_2 = this->initial_guess_->valueDerivativeFEM(fe, P_ID, 1, N_COUPLED); // Derivative of p respecto y
                // Average value of n in G.Ps
                double n_0 = this->initial_guess_->valueFEM(fe, N_ID, N_COUPLED); // amount of n in G.P
                double p_0 = this->initial_guess_->valueFEM(fe, P_ID, N_COUPLED); // amount of p in G.P
                //  if (number!=1 || number!=4 || number!=7)
                Jn_y(0) += mu_n * q * Vt * (-n_0 * phi_2 + n_2); //Jny in G.P
                Jn_x(0) += mu_n * q * Vt * (-n_0 * phi_1 + n_1); //Jnx in G.P
                Jp_y(0) += mu_p * q * Vt * (-p_0 * phi_2 - p_2); //Jpy in G.P
                Jp_x(0) += mu_p * q * Vt * (-p_0 * phi_1 - p_1); //Jpx in G.P
                number++;
            }
            //PrintError(number);
            Jn_y(0) = input_->c_unhat(Jn_y(0)) / input_->x0();  //Changing J to dimensional form
            Jn_x(0) = input_->c_unhat(Jn_x(0)) / input_->x0();  //Changing J to dimensional form
            Jp_y(0) = input_->c_unhat(Jp_y(0)) / input_->x0();  //Changing J to dimensional form
            Jp_x(0) = input_->c_unhat(Jp_x(0)) / input_->x0();  //Changing J to dimensional form

            phi_y(0) = phi_2; // Storing phi_y
            phi_x(0) = phi_1; // Storing phi_x
            n_y(0) = n_2; // Storing n_y
            n_x(0) = n_1; // Storing n_x
            p_y(0) = p_2; // Storing p_y
            p_x(0) = p_1; // Storing p_x
            elemental_Jny.at(elm_ID) = 0.1 * Jn_y(0) / number;  //Divided by 10, so the the unit is : mA/cm^-2
            elemental_Jnx.at(elm_ID) = 0.1 * Jn_x(0) / number;  //Divided by 10, so the the unit is : mA/cm^-2
            elemental_Jpy.at(elm_ID) = 0.1 * Jp_y(0) / number;  //Divided by 10, so the the unit is : mA/cm^-2
            elemental_Jpx.at(elm_ID) = 0.1 * Jp_x(0) / number;  //Divided by 10, so the the unit is : mA/cm^-2
            elemental_phiy.at(elm_ID) = phi_y(0) / number;  //Divided by 10, so the the unit is : mA/cm^-2
            elemental_phix.at(elm_ID) = phi_x(0) / number;  //Divided by 10, so the the unit is : mA/cm^-2
            elemental_ny.at(elm_ID) = n_y(0) / number;  //Divided by 10, so the the unit is : mA/cm^-2
            elemental_nx.at(elm_ID) = n_x(0) / number;  //Divided by 10, so the the unit is : mA/cm^-2
            elemental_py.at(elm_ID) = p_y(0) / number;  //Divided by 10, so the the unit is : mA/cm^-2
            elemental_px.at(elm_ID) = p_x(0) / number;  //Divided by 10, so the the unit is : mA/cm^-2

            //  Storing all of the parametrs in J matrix

            J(elm_ID, 0) = elemental_Jny.at(elm_ID); //Jny
            J(elm_ID, 1) = elemental_Jnx.at(elm_ID); //Jnx
            J(elm_ID, 2) = elemental_Jpy.at(elm_ID); //Jpy
            J(elm_ID, 3) = elemental_Jpx.at(elm_ID); //Jpx
            //J(elm_ID,4) = J(elm_ID,0); //Jny post process
            //J(elm_ID,5) = J(elm_ID,2); //Jpy post process
            J(elm_ID, 6) = elemental_phix.at(elm_ID); //phi_x
            J(elm_ID, 7) = elemental_phiy.at(elm_ID); //phi_y
            J(elm_ID, 8) = elemental_nx.at(elm_ID); //n_x
            J(elm_ID, 9) = elemental_ny.at(elm_ID); //n_y
            J(elm_ID, 10) = elemental_px.at(elm_ID); //p_x
            J(elm_ID, 11) = elemental_py.at(elm_ID);  //p_y

        }                           // end elements loop
        return (J);
    }

    /**
     * This function checks if there is any negative values for jn and jp.
     * Because negative values physically are not meaningful so it will detect and eliminate them.
    */
    void checkJ(const DDInputData *input_) {
        for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
            DDNodeData *pData = &(GetNodeData(nodeID));
            if (pData->j[0] < 0)
                pData->j[0] *= -1;
            if (pData->j[2] < 0)
                pData->j[2] *= -1;
        }
    }

    /**
     * This function calculates gradient of jn and gradient of Jp and finally
     * gradient of Jc. Calculated values are in G.Ps so we need to map them
     * to nodes. For using this map one need to call LSQ function. index of
     * J(elm_ID, 12) -> grad of Jn
     * J(elm_ID, 13) -> grad of Jp
     * J(elm_ID, 14) -> grad of Jc
     */
    void findGrad(DDGridField &data) {
        FEMElm fe(p_grid_, BASIS_FIRST_DERIVATIVE | BASIS_POSITION);
        FEMElm fe_linear(p_grid_, BASIS_ALL);
        ZEROPTV g_j;               // Value of J in G.Ps
        const double n_elements = p_grid_->n_elements();
        std::vector<double> elemental_grad_Jn(n_elements); // Storage vector for grad jn
        std::vector<double> elemental_grad_Jp(n_elements); // Storage vector for grad jp
        std::vector<double> elemental_grad_J(n_elements); // Storage vector for grad jc
        double grad_jn; // Declaration of grad_jn
        double grad_jp; // Declaration of grad_jp
        // Looping over elements
        double norm = 0;
        for (int elm_ID = 0; elm_ID < n_elements; elm_ID++) {
            grad_jn = 0;
            grad_jp = 0;
            fe.refill(elm_ID, BASIS_LINEAR, -1);
            int number = 0;
            // Loping over G.ps
            while (fe.next_itg_pt()) {
                double dy_jny = this->initial_guess_->valueDerivativeFEM(fe, JNY_ID, 1,
                                                                         N_COUPLED); // Derivative of jny respecto y
                double dx_jnx = this->initial_guess_->valueDerivativeFEM(fe, JNX_ID, 0,
                                                                         N_COUPLED); // Derivative of jnx respecto x
                grad_jn += dy_jny + dx_jnx;
                double dy_jpy = this->initial_guess_->valueDerivativeFEM(fe, JPY_ID, 1,
                                                                         N_COUPLED); // Derivative of jpy respecto y
                double dx_jpx = this->initial_guess_->valueDerivativeFEM(fe, JPX_ID, 0,
                                                                         N_COUPLED); // Derivative of jpx respecto x
                grad_jp += dy_jpy + dx_jpx;
                number++;
            }
            g_j(0) = grad_jn;
            g_j(1) = grad_jp;
            g_j(2) = grad_jn + grad_jp;
            elemental_grad_J.at(elm_ID) = 0.01 * g_j(2) / number;  //Storing grad_j
            norm += g_j(0) * g_j(0);
        }                           // end elements loop
        norm = sqrt(norm);
        std::vector<std::vector<int >> node_to_elems(p_grid_->n_nodes());
        for (int elm_ID = 0; elm_ID < n_elements; elm_ID++) {
            ELEM *elem = p_grid_->GetElm(elm_ID);
            for (int lcl_node = 0; lcl_node < elem->n_nodes(); lcl_node++) {
                int node_ID = elem->ElemToLocalNodeID(lcl_node);
                node_to_elems.at(node_ID).push_back(elm_ID);
            }
        }
        for (int node_ID = 0; node_ID < p_grid_->n_nodes(); node_ID++) {
            double val_4 = 0;
            std::vector<int> &elms = node_to_elems.at(node_ID);
            for (int e = 0; e < elms.size(); e++) {
                int elm_ID = elms.at(e);
                val_4 += elemental_grad_J.at(elm_ID);
            }
            val_4 /= (norm * elms.size());

            // save on DDNodeData
            DDNodeData &pData = GetNodeData(node_ID);
            pData.j[8] = val_4;
        }
    }

    /**
     * Compute and stores new values of jn and jp by doing post process.
     * Generally, in this method first we need to find the node that Jn and Jp start to give bad results
     * After finding that point we will do 4th order forward and backward finite difference to
     * substitute Jn and Jp values
     * At the end of the function values of grad of jn and grad of Jp are computed as well to
     * test accuracy of results
     */
    void post_process(const DDInputData *input_) {
        FEMElm fe(p_grid_, BASIS_FIRST_DERIVATIVE | BASIS_POSITION);
        double delta_jnx = 0;
        double delta_jpx = 0;
        int accuracy_order = 4;
        int n_tot = (input_->Nelem[0] + 1) * (input_->Nelem[1] + 1);
        int nx = (input_->Nelem[0] + 1);
        double hx = input_->Height / input_->Nelem[0];
        double hy = (input_->nsd >= 2) ? input_->Height / input_->Nelem[1] : 0.0;
        int jn_point = 0;
        int jp_point = 0;
        // Finding the node that Jn starts to have bad results
        for (int nodeID = p_grid_->n_nodes() - 1; nodeID >= 0; nodeID--) {
            DDNodeData *pData = &(GetNodeData(nodeID));
            double DT = pData->dt;
            if (nodeID <= 1 * nx || nodeID > n_tot - 2 * nx)
                continue;
            DDNodeData *pData1 = &(GetNodeData(nodeID + nx));  // data of above node
            DDNodeData *pData2 = &(GetNodeData(nodeID + 1)); // data of right node
            double y_1 = p_grid_->GetCoord(nodeID, 1);
            if (DT < -1)              // if we are donor phase
                pData->j[0] = 0;
            double delta_1 = fabs((pData1->j[0] - pData->j[0]) / hy + (pData2->j[1] - pData->j[1]) / hx);
            if (delta_1 > 20 && y_1 < 0.2)  // if grad of Jp is a large number (not reasonable)
            {
                jn_point = nodeID;
                std::cout << " At y = " << y_1 << " and node_ID = " << nodeID << " Jp post process started"
                          << std::endl;
                break;
            }
        }
        // Post process of Jn. From the node we found in previous step we do 4th Forward difference to substitute bad
        // results
        for (int nodeID = jn_point; nodeID >= 0; nodeID--) {
            DDNodeData *pData = &(GetNodeData(nodeID));
            double DT = pData->dt;
            int nex = (input_->Nelem[0]);
            int ney = (input_->Nelem[1]);
            float factor = 0;
            if (nex < ney) {
                factor = float(nex) / ney;
            } else {
                factor = float(ney) / nex;
            }
            factor *= 0.1;
            if (DT > -1.) {                         // if we are in acceptor phase
                if (nodeID > n_tot - accuracy_order * nx)
                    continue;
                DDNodeData *pData1 = &(GetNodeData(nodeID + (accuracy_order - 3) * nx)); // node data of above node
                DDNodeData *pData2 = &(GetNodeData(nodeID + (accuracy_order - 2) * nx)); // node data of 2nd above node
                DDNodeData *pData3 = &(GetNodeData(nodeID + (accuracy_order - 1) * nx)); // node data of 3rd above node
                DDNodeData *pData4 = &(GetNodeData(nodeID + accuracy_order * nx)); // node data of 4th above node
                double Jny_1 = pData1->j[4];  //Jny of above node
                double Jny_2 = pData2->j[4];  //Jny of 2nd above node
                double Jny_3 = pData3->j[4];  //Jny of 3rd above node
                double Jny_4 = pData4->j[4];  //Jny of 4th above node
                double Jny_i = 12 / 25. * (4. * Jny_1 - 3. * Jny_2 + 4 / 3. * Jny_3 -
                                           1 / 4. * Jny_4);  // 4th order backward difference
                if (nodeID % nx == 0) {
                    DDNodeData *pData5 = &(GetNodeData(nodeID + 1));
                    delta_jnx = pData5->j[1] - pData->j[1];
                } else {
                    DDNodeData *pData6 = &(GetNodeData(nodeID - 1));
                    delta_jnx = pData->j[1] - pData6->j[1];
                }
                if (nodeID < 3 * nx) {
                    pData->j[4] = Jny_i;  // substituting wrong Jny by finite difference Jny
                } else {
                    pData->j[4] = Jny_i + factor * fabs(delta_jnx);  // substituting Jny by finite difference Jny
                }
            } else {
                pData->j[4] = 0;
            }
        }
        /**
         * Finding the node that Jp starts to have bad results
         */
        for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
            DDNodeData *pData = &(GetNodeData(nodeID));
            double DT = pData->dt;

            if (nodeID > n_tot - 2 * nx)
                continue;
            if (DT > 2.)              // If we are in acceptor phase
                pData->j[2] = 0;
            DDNodeData *pData1 = &(GetNodeData(nodeID + nx));  // data of above node
            DDNodeData *pData2 = &(GetNodeData(nodeID + 1)); // data of right node
            double y_1 = p_grid_->GetCoord(nodeID, 1);
            double delta_1 = fabs((pData1->j[2] - pData->j[2]) / hy + (pData2->j[3] - pData->j[3]) / hx);
            if (delta_1 > 20 && y_1 > 0.8)  // if grad of Jp is a large number (not reasonable)
            {
                jp_point = nodeID;
                std::cout << " At y = " << y_1 << " and node_ID = " << nodeID << " Jp post process started"
                          << std::endl;
                break;
            }
        }
        /**
         * Post process of Jp. From the node we found in previous step we do 4th backward difference
         * to substitute bad results
         */
        for (int nodeID = jp_point; nodeID < n_tot; nodeID++) {
            DDNodeData *pData = &(GetNodeData(nodeID));
            double DT = pData->dt;
            int nex = (input_->Nelem[0]);
            int ney = (input_->Nelem[1]);
            float factor = 0;
            if (nex < ney) {
                factor = float(nex) / ney;
            } else {
                factor = float(ney) / nex;
            }
            factor *= 0.1;
            if (DT < 1.) {
                DDNodeData *pData1 = &(GetNodeData(nodeID - 1 * nx));  // node data of down node
                DDNodeData *pData2 = &(GetNodeData(nodeID - 2 * nx));  // node data of 2nd down node
                DDNodeData *pData3 = &(GetNodeData(nodeID - 3 * nx));  // node data of 3rd down node
                DDNodeData *pData4 = &(GetNodeData(nodeID - 4 * nx));  // node data of 4th down node
                double Jpy_1 = pData1->j[5];  //Jp of down node
                double Jpy_2 = pData2->j[5];  //Jp of 2nd down node
                double Jpy_3 = pData3->j[5];  //Jp of 3rd down node
                double Jpy_4 = pData4->j[5];  //Jp of 4th down node
                double Jpy_i = 12 / 25. * (4. * Jpy_1 - 3. * Jpy_2 + 4 / 3. * Jpy_3 -
                                           1 / 4. * Jpy_4);  // 4th order forward difference
                if (nodeID % nx == 0) {
                    //   DDNodeData *pData5 = &(GetNodeData (nodeID + 1));
                } else {
                    DDNodeData *pData6 = &(GetNodeData(nodeID - 1));
                    delta_jpx = pData->j[3] - pData6->j[3];
                }
                if (nodeID >= n_tot - 4 * nx) {
                    if (DT > -3 && DT < 1) {
                        pData->j[5] = pData1->j[5] +
                                      0.001 * hy * fabs(delta_jpx); // substituting Jny by finite difference Jny
                    } else {
                        pData->j[5] =
                                Jpy_i + 0.1 * hy * fabs(delta_jpx);
                    }
                } else {
                    pData->j[5] = Jpy_i + factor * fabs(delta_jpx);
                }

            } else {
                pData->j[5] = 0;
            }
        }
        for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
            DDNodeData *pData = &(GetNodeData(nodeID));
            if (pData->j[4] < 0)
                pData->j[4] *= -1;
            if (pData->j[5] < 0)
                pData->j[5] *= -1;
        }
    }

    /**
     * Reading data from bits file
     */
    void SetCoarseData_txt(const std::string &microstructure_filename) {
        if (input_->nsd < 3) {
            int Nx_coarse = input_->elmnx + 1;
            int Ny_coarse = input_->elmny + 1;
            int transLayer = int(input_->elmny * input_->transportLayer);
            int extendNy = input_->elmny + 1 + transLayer * 2;
            morph2d.resize(extendNy);
            DT2d.resize(extendNy);
            DT2ddx.resize(extendNy);
            DT2ddy.resize(extendNy);
            for (int i = 0; i < extendNy; i++) {
                morph2d[i].resize(Nx_coarse);
                DT2d[i].resize(Nx_coarse);
                DT2ddx[i].resize(Nx_coarse);
                DT2ddy[i].resize(Nx_coarse);
            }
            PrintStatus("Reading file : ", microstructure_filename.c_str(), "\n");
            FILE *fp = fopen(microstructure_filename.c_str(), "r");
            static char value[40000];
            char *position[40000];
            if (fp != nullptr) {
                getLine(fp, value);
                divideParameters(value, position, "(,)[]\t ");
                for (int j = 0; j < Ny_coarse; j++) {
                    for (int i = 0; i < Nx_coarse; i++) {
                        morph2d[j + transLayer][i] = atof(position[i + j * Nx_coarse]);
                        if (morph2d[j][i] < input_->threshold) morph2d[j][i] = 0;
                        else morph2d[j][i] = 1;
                    }
                }

                for (int j = transLayer - 1; j > -1; j--) {
                    for (int i = 0; i < Nx_coarse; i++) {
                        morph2d[j][i] = morph2d[j + 1][i];
                    }
                }
                for (int j = extendNy - 1 - transLayer; j < extendNy - 1; j++) {
                    for (int i = 0; i < Nx_coarse; i++) {
                        morph2d[j][i] = morph2d[j - 1][i];
                    }
                }

                fclose(fp);

            }
        } else {
            int Nx_coarse = input_->elmnx + 1;
            int Ny_coarse = input_->elmny + 1;
            int Nz_coarse = input_->elmnz + 1;
            int transLayer = int(input_->elmny * input_->transportLayer);
            int extendNy = input_->elmny + 1 + transLayer * 2;
            morph3d.resize(Nz_coarse);
            DT3d.resize(Nz_coarse);
            DT3ddx.resize(Nz_coarse);
            DT3ddy.resize(Nz_coarse);
            DT3ddz.resize(Nz_coarse);
            for (int i = 0; i < Nz_coarse; i++) {
                morph3d[i].resize(extendNy);
                DT3d[i].resize(extendNy);
                DT3ddx[i].resize(extendNy);
                DT3ddy[i].resize(extendNy);
                DT3ddz[i].resize(extendNy);

                for (int j = 0; j < extendNy; j++) {
                    morph3d[i][j].resize(Nx_coarse);
                    DT3d[i][j].resize(Nx_coarse);
                    DT3ddx[i][j].resize(Nx_coarse);
                    DT3ddy[i][j].resize(Nx_coarse);
                    DT3ddz[i][j].resize(Nx_coarse);
                }
            }

            PrintStatus("Reading file ...\n", microstructure_filename.c_str(), "\n");
            //FILE *fp = fopen (microstructure_filename.c_str (), "r");
            std::ifstream fin(microstructure_filename);
            //skip the first line
            //getLine (fp, value);
            //read the second line
            std::string line;
            getline(fin, line);
            std::stringstream ss(line);
            int morph = 0;
            //one row format
            /*
             for (int k = 0; k < Nz_coarse; k++) {
                 for (int i = 0; i < Ny_coarse; i++) {
                     for (int j = 0; j < Nx_coarse; j++) {
                         ss >> morph;
                         int id = j + i * Nx_coarse + k * (Nx_coarse * extendNy);
                         morph3d[k][i + transLayer][j] = morph;
                         if (morph3d[k][i + transLayer][j] < input_->threshold) morph3d[k][i + transLayer][j] = 0;
                         else morph3d[k][i + transLayer][j] = 1;
                     }
                 }
             }
            */

            for (int i = 0; i < Ny_coarse; i++) {
                for (int k = 0; k < Nz_coarse; k++) {
                    for (int j = 0; j < Nx_coarse; j++) {
                        fin >> morph;
                        morph3d[k][i + transLayer][j] = morph;
                        if (morph3d[k][i + transLayer][j] < input_->threshold) morph3d[k][i + transLayer][j] = 0;
                        else morph3d[k][i + transLayer][j] = 1;

                        /*
                        morph3d[i+transLayer][k][j] = morph;
                        if(morph3d[i+transLayer][k][j]<input_->threshold) morph3d[i+transLayer][k][j]=0;
                        else morph3d[i+transLayer][k][j] = 1;
                        */
                    }
                }
            }
            for (int i = transLayer - 1; i > -1; i--) {
                for (int k = 0; k < Nz_coarse; k++) {
                    for (int j = 0; j < Nx_coarse; j++) {
                        morph3d[k][i][j] = morph3d[k][i + 1][j];
                    }
                }
            }

            for (int i = extendNy - 1 - transLayer; i < extendNy - 1; i++) {
                for (int k = 0; k < Nz_coarse; k++) {
                    for (int j = 0; j < Nx_coarse; j++) {
                        morph3d[k][i][j] = morph3d[k][i - 1][j];
                    }
                }
            }

        }
        PrintStatus("Finish morphology reading");
    }

    /**
     * Computing surface(2d)/volume(3d) integral of diss and recombination and exciton generation.
     * @return value: 1- Dissociation 2- Recombination 3- Generation in vector form;
     */
    std::vector<double> domain_integrations(DDGridField &data) {
        int e;
        int nel = p_grid_->n_elements();
        std::vector<double> DRG(3);
        double intgD = 0, intgR = 0, intgG = 0;
        for (e = 0; e < nel; e++) {
            ELEM *pElm;
            pElm = this->p_grid_->GetElm(e);
            FEMElm fe(p_grid_, BASIS_ALL);
            fe.refill(e, BASIS_LINEAR, 0);
            while (fe.next_itg_pt()) {
                intgD += this->valueFEM(fe, 32) * fe.detJxW();
                intgR += this->valueFEM(fe, 33) * fe.detJxW();
                intgG += this->valueFEM(fe, 44) * fe.detJxW();
            }
        }
        DRG.at(0) = intgD;
        DRG.at(1) = intgR;
        DRG.at(2) = intgG;

        return DRG;
    }

    /**
    * Computing Jc, output current density.
    * @return value: current density of Jp, Jn and J in a vector form;
    */
    std::vector<float> integrateJ(DDGridField &data) {
        int j, e;
        int nel = p_grid_->n_elements();

        int boundary_bottom = 3;
        int boundary_top = 4;
        int flag, flag_bottom;

        double bfarea = 0, tparea = 0;  // for area integration

        double Jbot = 0;
        double Jtop = 0;
        double n_grad_phi(0.0), p_grad_phi(0.0), grad_n(0.0), grad_p(0.0);
        //loop over all the elements
        for (e = 0; e < nel; e++) {
            // initialize flags
            flag = 0;
            flag_bottom = 0;
            //STEP1:  Check if this element lies on either boundary
            // find the number of nodes in this element
            int nne = this->p_grid_->GetElm(e)->n_nodes();
            //Check to see if this node lies on boundary
            for (j = 0; j < nne; j++) {
                //Check to see if the element is a boundary element
                int nodeID = this->p_grid_->GetLocalNodeID(e, j);
                if (this->p_grid_->BoNode(nodeID, boundary_bottom)) {
                    flag = 1;
                    flag_bottom = 1;
                    break;
                } else if (this->p_grid_->BoNode(nodeID, boundary_top)) {
                    flag = 1;
                    break;
                }
            }
            // STEP 2: Find the line integral
            //-----------------------------------------------------
            if (flag) {
                ELEM *pElm;
                pElm = this->p_grid_->GetElm(e);
                FEMElm fe(p_grid_, BASIS_ALL);
                //set the indicator
                int interfaceInd = (flag_bottom > 0) ? boundary_bottom : boundary_top;
                //perform surface/line integration
                //loop over all the surfaces in the element
                for (ELEM::SurfaceList_type::iterator it = pElm->surface_indicator_.begin();
                     it != pElm->surface_indicator_.end(); it++) {
                    double area = 0;
                    //check to see if this surface is the correct surface
                    //cout << "indicators: " << it->indicators() << endl;
                    if (it->has_indicator(interfaceInd)) {
                        // loop over the surface gauss points
                        fe.refill_surface(e, &*it, BASIS_LINEAR, 0);
                        while (fe.next_itg_pt()) {
                            double detJxW = fe.detJxW();
                            double ms = data.valueFEM(fe, 38);
                            double Jny = data.valueFEM(fe, JNY_ID); //ms>0.5?data.valueFEM (fe, JNY_ID):0;
                            double Jpy = data.valueFEM(fe, JPY_ID); //ms<0.5?data.valueFEM (fe, JPY_ID):0;
                            if (interfaceInd == boundary_bottom) {
                                Jbot += Jny * detJxW;
                                grad_n += data.valueDerivativeFEM(fe, N_ID, 1) * detJxW;
                                n_grad_phi += data.valueFEM(fe, N_ID) * data.valueDerivativeFEM(fe, PHI_ID, 1) * detJxW;
                            }
                            else {
                                Jtop += Jpy * detJxW;
                                grad_p += data.valueDerivativeFEM(fe, P_ID, 1) * detJxW;
                                p_grad_phi += data.valueFEM(fe, P_ID) * data.valueDerivativeFEM(fe, PHI_ID, 1) * detJxW;
                            }
                            area += detJxW;
                        }
                        if (flag_bottom) {
                            bfarea += area;
                        } else {
                            tparea += area;
                        }
                    }
                }
            }
        }
        // send current density integral to proc 0
        double rJtop;
        double rJbot;
        MPI_Allreduce (&Jtop, &rJtop, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        MPI_Allreduce (&Jbot, &rJbot, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        // send area to proc 0
        double rtop;
        double rbot;
        MPI_Allreduce (&tparea, &rtop, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        MPI_Allreduce (&bfarea, &rbot, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        tparea = rtop;
        bfarea = rbot;
        rJtop /= tparea;
        rJbot /= bfarea;
        std::vector<float> J;
        J.push_back(-rJtop);
        J.push_back(-rJbot);
        rJtop > rJbot ? J.push_back(-rJbot) : J.push_back(-rJtop);
        grad_p /= tparea;
        grad_n /= bfarea;
        n_grad_phi /= bfarea;
        p_grad_phi /= tparea;
        PrintStatus("Calculated values Non-dimensional : ",
                    "\n\tgrad_p: ", grad_p,
                    "\n\tgrad_n: ", grad_n,
                    "\n\tn_grad_phi: ", n_grad_phi,
                    "\n\tp_grad_phi: ", p_grad_phi
        );
        PrintStatus("Calculated values dimensional : ",
                    "\n\tgrad_p: ", 1/input_->x0() * input_->c_unhat(grad_p),
                    "\n\tgrad_n: ", 1/input_->x0() * input_->c_unhat(grad_n),
                    "\n\tn_grad_phi: ", 1/input_->x0() * input_->c_unhat(n_grad_phi) * input_->V0(),
                    "\n\tp_grad_phi: ", 1/input_->x0() * input_->c_unhat(p_grad_phi) * input_->V0()
        );
        return J;
    }

    /**
     * This function changes regular mesh to clustered mesh.
     * It means we will have more mesh density near to cathode and anode.
     * Goal: By using this function you can obtain more accurate results.
     */
    void mesh_clustering(const DDInputData *input_data) {
        // Grid clustering
        double y, y_new, eta, x_hat = 1.;
        double beta = 1 + input_data->beta_1;
        for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
            y = p_grid_->GetCoord(nodeID, 1);
            eta = y / x_hat;
            y_new =
                    x_hat * ((1 + beta) * pow((beta + 1) / (beta - 1), 2 * eta - 1)
                             + 1 - beta) / (2 * (1 + pow((beta + 1) / (beta - 1), 2 * eta - 1)));
            p_grid_->node_array_[nodeID]->setCoor(1, y_new);
        }
        PrintStatus("finished Interpolate initialize...\n");
    }

    /**
    * Computing Jc, output current density.
    */
    std::vector<float> FindJc(const DDInputData *input_) {
        std::vector<float> J;
        int nx = (input_->elmnx + 1);
        int n_tot = p_grid_->n_nodes();
        float jc1 = 0;
        float jc2 = 0;
        for (int nodeID = 0; nodeID < nx; nodeID++) {
            DDNodeData *pData = &(GetNodeData(nodeID));
            jc1 += pData->j[0];
        }
        for (int nodeID = n_tot - 1; nodeID >= n_tot - nx; nodeID--) {
            DDNodeData *pData = &(GetNodeData(nodeID));
            jc2 += pData->j[2];
        }
        J.push_back(-jc1 / nx);
        J.push_back(-jc2 / nx);
        J.push_back((jc1 < jc2 ? -jc1 : -jc2) / nx);
        return J;
    }

    /**
     * Result checkting
     */
    double FindMaxGradJ(DDGridField &data) {
        double max_j = 0;
        for (int node_ID = 0; node_ID < p_grid_->n_nodes(); node_ID++) {
            DDNodeData &pData = GetNodeData(node_ID);
            double y_1 = p_grid_->GetCoord(node_ID, 1);
            double DT = pData.dt;
            if (DT < -2 || (DT > 2 && y_1 > 0.1 && y_1 < 0.9)) {
                if (pData.j[8] > max_j)
                    max_j = pData.j[8];
            }
        }
        return max_j;
    }


    /**
    * Reading morphology from plt input file and extending it to incorporate transport layer
    */
    void SetCoarseData(const std::string &MSfile) {
        PrintStatus("Reading .plt file! ");
        //< NOTE: Set length scales correctly in Input parameters
        if (input_->nsd == 2) {
            //< Fill the microstructure array
            int Nx_coarse = input_->elmnx + 1;
            int Ny_coarse = input_->elmny + 1;
            int transLayer = int(input_->elmny * input_->transportLayer);
            int extendNy = input_->elmny + 1 + transLayer * 2;
            morph2d.resize(extendNy);
            DT2d.resize(extendNy);
            DT2ddx.resize(extendNy);
            DT2ddy.resize(extendNy);
            for (int i = 0; i < extendNy; i++) {
                morph2d[i].resize(Nx_coarse);
                DT2d[i].resize(Nx_coarse);
                DT2ddx[i].resize(Nx_coarse);
                DT2ddy[i].resize(Nx_coarse);
            }
            //< Read in the coarse data from the file
            PrintStatus("Reading file MSfile ", MSfile.c_str(), "\n");
            FILE *fp = fopen(MSfile.c_str(), "r");
            static char value[1024];
            char *position[Nx_coarse];
            if (fp != nullptr) {
                // skip first three lines
                for (int j = 1; j <= 3; j++)
                    getLine(fp, value);

                for (int j = 0; j < Ny_coarse; j++) {
                    for (int i = 0; i < Nx_coarse; i++) {
                        getLine(fp, value);
                        divideParameters(value, position, "(,)[]\t ");
                        morph2d[j + transLayer][i] = atof(position[2]);
                    }
                }

                // set bottom transport layer
                for (int j = transLayer - 1; j >= 0 ; j--) {
                    for (int i = 0; i < Nx_coarse; i++) {
                        morph2d[j][i] = 1.0;
                    }
                }

                // set top transport layer
                for (int j = extendNy - 1; j > extendNy + Ny_coarse - 1; j--) {
                    for (int i = 0; i < Nx_coarse; i++) {
                        morph2d[j][i] = 0.0;
                    }
                }

                getLine(fp, value);
                fclose(fp);
#ifdef NDEBUG
                if(!GetMPIRank()) {
                    PrintStatus("Number of elements after adding transport layer, coarse morph:\n");
                    std::cout << "\t Nx: " << Nx_coarse << ", Ny: " << extendNy << std::endl;
                    std::fstream fp_morph;
                    fp_morph.open("debug_coarse_data.plt", std::fstream::out);
                    for(int j = 0; j < extendNy; j++){
                        for(int i = 0; i < Nx_coarse; i++){
                            fp_morph << i << "\t" << j << "\t" << morph2d[j][i] << "\n";
                        }
                    }
                    fp_morph.close();
                }
#endif
            } else {
                throw TALYException()  << "Fatal error: no morphology file exists;";
            }
        } else if (input_->nsd == 3) {
            int Nx_coarse = input_->elmnx + 1;
            int Ny_coarse = input_->elmny + 1;
            int Nz_coarse = input_->elmnz + 1;
            int transLayer = int(input_->elmny * input_->transportLayer);
            int extendNy = input_->elmny + 1 + transLayer * 2;
            morph3d.resize(Nx_coarse);
            DT3d.resize(Nx_coarse);
            DT3ddx.resize(Nx_coarse);
            DT3ddy.resize(Nx_coarse);
            DT3ddz.resize(Nx_coarse);
            //int nodeno_coarse = Nx_coarse*Ny_coarse*Nz_coarse;
            //< Read in the coarse data from the file
            for (int i = 0; i < Nx_coarse; i++) {
                morph3d[i].resize(extendNy);
                DT3d[i].resize(extendNy);
                DT3ddx[i].resize(extendNy);
                DT3ddy[i].resize(extendNy);
                DT3ddz[i].resize(extendNy);
                for (int j = 0; j < extendNy; j++) {
                    morph3d[i][j].resize(Nz_coarse);
                    DT3d[i][j].resize(Nz_coarse);
                    DT3ddx[i][j].resize(Nz_coarse);
                    DT3ddy[i][j].resize(Nz_coarse);
                    DT3ddz[i][j].resize(Nz_coarse);
                }
            }
            PrintStatus("reading file MSfile ", MSfile.c_str(), "\n");
            FILE *fp = fopen(MSfile.c_str(), "r");
            static char value[1024];
            char *position[Nx_coarse];
            //int itr_ = 1;
            // skip first three lines
            for (int j = 1; j <= 3; j++) getLine(fp, value);
            for (int k = 0; k < Nz_coarse; k++) {
                for (int j = 0; j < Ny_coarse; j++) {
                    for (int i = 0; i < Nx_coarse; i++) {
                        getLine(fp, value);
                        divideParameters(value, position, "(,)[]\t ");
                        morph3d[i][j + transLayer][k] = atof(position[3]);
                        if (morph3d[i][j + transLayer][k] < input_->threshold) morph3d[i][j + transLayer][k] = 0;
                        else morph3d[i][j + transLayer][k] = 1;
                    }
                }
            }
            for (int k = 0; k < Nz_coarse; k++) {
                for (int j = transLayer - 1; j > -1; j--) {
                    for (int i = 0; i < Nx_coarse; i++) {
                        morph3d[i][j][k] = morph3d[i][j + 1][k];
                    }
                }
            }
            for (int k = 0; k < Nz_coarse; k++) {
                for (int j = extendNy - 1 - transLayer; j < extendNy; j++) {
                    for (int i = 0; i < Nx_coarse; i++) {
                        morph3d[i][j][k] = morph3d[i][j - 1][k];
                    }
                }
            }
        }

        PrintStatus("Finished SetCoarseData\n");
    }

    // set surface indicator for each boundary
    void SetIndicators(DDInputData &inputData) {
        double tol = 1e-6;
        int transLayer = int(input_->elmny * input_->transportLayer);
        int extendNy = input_->elmny + 1 + transLayer * 2;
        morph2d.resize(extendNy);
        DT2d.resize(extendNy);
        DT2ddx.resize(extendNy);
        DT2ddy.resize(extendNy);
        for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
            NodeIndicator indicators = 0;
            double x = p_grid_->GetCoord(nodeID, 0);
            double y = p_grid_->GetCoord(nodeID, 1);

            if (fabs(x) < tol) {
                indicators |= INDICATOR_NUM(1);
            }
            if (fabs(x - inputData.L[0]) < tol) {
                indicators |= INDICATOR_NUM(2);
            }
            if (fabs(y) < tol) {
                indicators |= INDICATOR_NUM(3);
            }
            if (fabs(y - 1.0) < tol) {
                indicators |= INDICATOR_NUM(4);
            }
            if (input_->nsd == 3) { //if 3d then set the front, back surface id
                double z = p_grid_->GetCoord(nodeID, 2);
                if (fabs(z) < tol) {
                    indicators |= INDICATOR_NUM(5);
                }
                if (fabs(z - inputData.L[2]) < tol) {
                    indicators |= INDICATOR_NUM(6);
                }
            }
            p_grid_->GetNode(nodeID)->setIndicators(indicators);
        }
        p_grid_->cared_surface_indicator_.appendData(1);
        p_grid_->cared_surface_indicator_.appendData(2);
        p_grid_->cared_surface_indicator_.appendData(3);
        p_grid_->cared_surface_indicator_.appendData(4);
        if (input_->nsd == 3) {
            p_grid_->cared_surface_indicator_.appendData(5);
            p_grid_->cared_surface_indicator_.appendData(6);
        }
        p_grid_->GenElmSurfaceIndicator();
        PrintStatus("finished SetIndicators");
    }


    void interpolate(double *dt, double *dtdx, double *dtdy, double *dtdz = NULL) {
        double valms, valdt, valdtdx, valdtdy, valdtdz;
        //first map the 1d dt(distance contour) to a 3d structure;
        int Nx_coarse = input_->elmnx + 1;
        int Ny_coarse = input_->elmny + 1;
        int transLayer = int(input_->elmny * input_->transportLayer);
        int extendNy = input_->elmny + 1 + transLayer * 2;
        int Nz_coarse = input_->elmnz + 1;
        int id = 0;
        if (input_->nsd == 2) {
            for (int jj = 0; jj < extendNy; jj++) {
                for (int ii = 0; ii < Nx_coarse; ii++) {
                    id = ii + jj * Nx_coarse;
                    DT2d[jj][ii] = dt[id];
                    DT2ddx[jj][ii] = dtdx[id];
                    DT2ddy[jj][ii] = dtdy[id];
                }
            }
        } else if (input_->nsd == 3) {
            for (int k = 0; k < Nz_coarse; k++) {
                for (int j = 0; j < extendNy; j++) {
                    for (int i = 0; i < Nx_coarse; i++) {
                        int id = i + j * Nx_coarse + k * (Nx_coarse * extendNy);
                        DT3d[i][j][k] = dt[id];
                        DT3ddx[i][j][k] = dtdx[id];
                        DT3ddy[i][j][k] = dtdy[id];
                        DT3ddz[i][j][k] = dtdz[id];
                    }
                }
            }
        }

        int nodeNum = p_grid_->n_nodes();
        // cell size in each direction on the coarse mesh;
        double dy = 1. / input_->elmny * 1.0001;
        double dx = dy;             //the nature of the coarse mesh;
        if (input_->nsd == 2) {
            for (int nodeID = 0; nodeID < nodeNum; nodeID++) {
                //get the coordinates of node point.
                ZEROPTV pt_gb;
                pt_gb(0) = p_grid_->GetCoord(nodeID, 0);
                pt_gb(1) = p_grid_->GetCoord(nodeID, 1);
                //determine which cell the node belongs to on the coarse grid: TODO : check this definition
                int x_index = floor(pt_gb(0) / dx);
                int y_index = floor(pt_gb(1) / dy);
                // interpolate with bilinear interpolation;
                // First map the small cell to  1by1 domain and then do the interpolation;
                double local_x = (pt_gb(0) / dx - x_index);
                double local_y = (pt_gb(1) / dy - y_index);
#ifdef NDEBUG
                assert(x_index < input_->Nelem[0]);
                assert(y_index < input_->Nelem[1]);
                assert(local_x < dx);
                assert(local_y < dy);
#endif
                valms =
                        morph2d[y_index][x_index] * (1. - local_x) * (1. - local_y) +
                        morph2d[y_index][x_index + 1] * local_x * (1. - local_y) +
                        morph2d[y_index + 1][x_index] * (1. - local_x) * local_y +
                        morph2d[y_index + 1][x_index + 1] * local_x * local_y;
                valdt =
                        DT2d[y_index][x_index] * (1. - local_x) * (1. - local_y) +
                        DT2d[y_index][x_index + 1] * local_x * (1. - local_y) +
                        DT2d[y_index + 1][x_index] * (1. - local_x) * local_y +
                        DT2d[y_index + 1][x_index + 1] * local_x * local_y;
                valdtdx =
                        DT2ddx[y_index][x_index] * (1. - local_x) * (1. - local_y) +
                        DT2ddx[y_index][x_index + 1] * local_x * (1. - local_y) +
                        DT2ddx[y_index + 1][x_index] * (1. - local_x) * local_y +
                        DT2ddx[y_index + 1][x_index + 1] * local_x * local_y;
                valdtdy =
                        DT2ddy[y_index][x_index] * (1. - local_x) * (1. - local_y) +
                        DT2ddy[y_index][x_index + 1] * local_x * (1. - local_y) +
                        DT2ddy[y_index + 1][x_index] * (1. - local_x) * local_y +
                        DT2ddy[y_index + 1][x_index + 1] * local_x * local_y;
                //set the interpolated value to the fine mesh;
                DDNodeData *pData = &(GetNodeData(nodeID));
                pData->dt = valdt;
                pData->Dis = valdtdx;
                pData->R = valdtdy;
                pData->ms = valms;
                // acceptor layer
                if (pt_gb(1) < float(transLayer) / float(extendNy)) {
                    pData->ms = 0.99;
                    pData->dt = 2000;
                }
                // donor layer
                if (pt_gb(1) > 1. - 1. * float(transLayer) / float(extendNy)) {
                    pData->ms = 0.01;
                    pData->dt = -2000;
                }

            }
        } else if (input_->nsd == 3) {
            double dz = dx;
            for (int nodeID = 0; nodeID < nodeNum; nodeID++) {
                ZEROPTV pt_gb;
                pt_gb(0) = p_grid_->GetCoord(nodeID, 0);
                pt_gb(1) = p_grid_->GetCoord(nodeID, 1);
                pt_gb(2) = p_grid_->GetCoord(nodeID, 2);
                int ruler_x = floor(pt_gb(0) / dx);
                int ruler_y = floor(pt_gb(1) / dy);
                int ruler_z = floor(pt_gb(2) / dz);
                float pt_loc[3];
                pt_loc[0] = (pt_gb(0) / dx - ruler_x);
                pt_loc[1] = (pt_gb(1) / dy - ruler_y);
                pt_loc[2] = (pt_gb(2) / dz - ruler_z);
                //do Tri-linear interplation mapping on to the calculation mesh
                if (pt_gb(1) < 1.1 * float(transLayer) / float(extendNy)) {
                    DDNodeData *pData = &(GetNodeData(nodeID));
                    pData->dt = 0;
                    pData->Dis = 0;
                    pData->ms = 0.99;
                    pData->dt = 2000;
                    pData->j[8] = 0;
                } else if (pt_gb(1) > 1. - 1.1 * float(transLayer) / float(extendNy)) {
                    DDNodeData *pData = &(GetNodeData(nodeID));
                    pData->dt = 0;
                    pData->Dis = 0;
                    pData->ms = 0.01;
                    pData->dt = -2000;
                    pData->j[8] = 0;
                } else {
                    LinearInterpolator<double> interp;
                    double v[8];
                    v[0] = morph3d[ruler_x][ruler_y][ruler_z];
                    v[1] = morph3d[ruler_x + 1][ruler_y][ruler_z];
                    v[2] = morph3d[ruler_x][ruler_y + 1][ruler_z];
                    v[3] = morph3d[ruler_x + 1][ruler_y + 1][ruler_z];
                    v[4] = morph3d[ruler_x][ruler_y][ruler_z + 1];
                    v[5] = morph3d[ruler_x + 1][ruler_y][ruler_z + 1];
                    v[6] = morph3d[ruler_x][ruler_y + 1][ruler_z + 1];
                    v[7] = morph3d[ruler_x + 1][ruler_y + 1][ruler_z + 1];
                    valms = interp.Trilinear(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], pt_loc[0], pt_loc[1],
                                             pt_loc[2]);
                    v[0] = DT3d[ruler_x][ruler_y][ruler_z];
                    v[1] = DT3d[ruler_x + 1][ruler_y][ruler_z];
                    v[2] = DT3d[ruler_x][ruler_y + 1][ruler_z];
                    v[3] = DT3d[ruler_x + 1][ruler_y + 1][ruler_z];
                    v[4] = DT3d[ruler_x][ruler_y][ruler_z + 1];
                    v[5] = DT3d[ruler_x + 1][ruler_y][ruler_z + 1];
                    v[6] = DT3d[ruler_x][ruler_y + 1][ruler_z + 1];
                    v[7] = DT3d[ruler_x + 1][ruler_y + 1][ruler_z + 1];
                    valdt = interp.Trilinear(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], pt_loc[0], pt_loc[1],
                                             pt_loc[2]);
                    v[0] = DT3ddx[ruler_x][ruler_y][ruler_z];
                    v[1] = DT3ddx[ruler_x + 1][ruler_y][ruler_z];
                    v[2] = DT3ddx[ruler_x][ruler_y + 1][ruler_z];
                    v[3] = DT3ddx[ruler_x + 1][ruler_y + 1][ruler_z];
                    v[4] = DT3ddx[ruler_x][ruler_y][ruler_z + 1];
                    v[5] = DT3ddx[ruler_x + 1][ruler_y][ruler_z + 1];
                    v[6] = DT3ddx[ruler_x][ruler_y + 1][ruler_z + 1];
                    v[7] = DT3ddx[ruler_x + 1][ruler_y + 1][ruler_z + 1];
                    valdtdx = interp.Trilinear(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], pt_loc[0], pt_loc[1],
                                               pt_loc[2]);
                    v[0] = DT3ddy[ruler_x][ruler_y][ruler_z];
                    v[1] = DT3ddy[ruler_x + 1][ruler_y][ruler_z];
                    v[2] = DT3ddy[ruler_x][ruler_y + 1][ruler_z];
                    v[3] = DT3ddy[ruler_x + 1][ruler_y + 1][ruler_z];
                    v[4] = DT3ddy[ruler_x][ruler_y][ruler_z + 1];
                    v[5] = DT3ddy[ruler_x + 1][ruler_y][ruler_z + 1];
                    v[6] = DT3ddy[ruler_x][ruler_y + 1][ruler_z + 1];
                    v[7] = DT3ddy[ruler_x + 1][ruler_y + 1][ruler_z + 1];
                    valdtdy = interp.Trilinear(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], pt_loc[0], pt_loc[1],
                                               pt_loc[2]);
                    v[0] = DT3ddz[ruler_x][ruler_y][ruler_z];
                    v[1] = DT3ddz[ruler_x + 1][ruler_y][ruler_z];
                    v[2] = DT3ddz[ruler_x][ruler_y + 1][ruler_z];
                    v[3] = DT3ddz[ruler_x + 1][ruler_y + 1][ruler_z];
                    v[4] = DT3ddz[ruler_x][ruler_y][ruler_z + 1];
                    v[5] = DT3ddz[ruler_x + 1][ruler_y][ruler_z + 1];
                    v[6] = DT3ddz[ruler_x][ruler_y + 1][ruler_z + 1];
                    v[7] = DT3ddz[ruler_x + 1][ruler_y + 1][ruler_z + 1];
                    valdtdz = interp.Trilinear(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], pt_loc[0], pt_loc[1],
                                               pt_loc[2]);
                    //set the interpolated value to the fine mesh;
                    DDNodeData *pData = &(GetNodeData(nodeID));
                    pData->dt = valdt;
                    pData->Dis = valdtdx;
                    //because Dis, R and j[8] are only used/calculated at post processing, I can safely store distance contour derivatives here.
                    pData->R = valdtdy;
                    pData->j[8] = valdtdz;
                    pData->ms = valms;
                }
            }
        }


    }

    /**
     * Calculate (postprocesing) dissociation and recombination;
     */
    void Calc_DR() {
        for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
            DDNodeData *pData = &(GetNodeData(nodeID));
            double dphidx = pData->j[6];
            double dphidy = pData->j[7];
            double dphidz = pData->j[5];
            double n = pData->n[0];
            double p = pData->p[0];
            double Dphi = pData->PHIBAR;
            double dDTdx = pData->Dis;  //derivatives of distance contour was stored in Dis and R.
            double dDTdy = pData->R;
            double dDTdz = 0;
            double absdDT = sqrt(dDTdx * dDTdx + dDTdy * dDTdy);
            if (p_grid_->nsd() == 3) {
                dDTdz = pData->j[8];
                absdDT += (dDTdz * dDTdz);
            }
            absdDT = sqrt(absdDT);
            double dVnorm = ((dphidx - 0 * dDTdx / absdDT * Dphi) * (dphidx - 0 * dDTdx / absdDT * Dphi)
                             + (dphidy - 0 * dDTdy / absdDT * Dphi) * (dphidy - 0 * dDTdy / absdDT * Dphi));
            // TODO : check this calculation with Baskar
            //if (p_grid_->nsd()==3) dVnorm += ((dphidz-dDTdz/absdDT*Dphi) * (dphidz-dDTdz/absdDT*Dphi));
            dVnorm = sqrt(dVnorm);
            double DT = pData->dt;
            double eps = pData->epsn;
            double mu_dim = pData->mup + pData->mun;
            double X = input_->c_unhat(pData->x[0]);
            double y_1 = p_grid_->GetCoord(nodeID, 1);
            pData->Dis = X * gendiss->calcKdiss(y_1, dVnorm, DT, eps, mu_dim);
            pData->R = gendiss->calcR(y_1, n, p, DT, eps, mu_dim);  // input_->q / eps * mu_dim * n * p;
        }
    }

    /**
     * Convert non-dimensional variable to dimensional form;
     */
    void dimlize() {
        PrintStatus("Dimensionalization is applied only the the primary unknowns, i.e., n,p,V,X.");
        for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
            DDNodeData *pData = &(GetNodeData(nodeID));
            pData->n[0] = input_->c_unhat(pData->n[0]);  // electron
            pData->p[0] = input_->c_unhat(pData->p[0]);  // hole
            pData->phi[0] = input_->V_unhat(pData->phi[0]);  // potential
            pData->x[0] = input_->c_unhat(pData->x[0]);  // exciton
            /*
             pData->phi[1] = input_->c_unhat(pData->phi[1])/input_->x0();
             pData->n[1] = input_->c_unhat(pData->n[1])/input_->x0();

             pData->p[1] = input_->V_unhat(pData->p[1])/input_->x0();
             pData->u[7] = input_->c_unhat(pData->u[7])/input_->x0();
             pData->phi[2] = input_->c_unhat(pData->phi[2])/input_->x0();
             pData->n[2] = input_->V_unhat(pData->n[2])/input_->x0()/input_->x0();
             pData->p[2] = input_->c_unhat(pData->p[2])/input_->x0()/input_->x0();
             pData->u[11] = input_->c_unhat(pData->u[11])/input_->x0()/input_->x0();
             */
        }
    }

    /**
     Convert dimensional variable to non-dimensional form;
    */
    void nondimlize() {
        for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
            DDNodeData *pData = &(GetNodeData(nodeID));
            pData->n[0] = input_->c_hat(pData->n[0]);
            pData->p[0] = input_->c_hat(pData->p[0]);
            pData->phi[0] = input_->V_hat(pData->phi[0]);
            pData->x[0] = input_->c_hat(pData->x[0]);
            /*
             pData->phi[1] = input_->c_hat(pData->phi[1])*input_->x0();
             pData->n[1] = input_->c_hat(pData->n[1])*input_->x0();
             pData->p[1] = input_->V_hat(pData->p[1])*input_->x0();
             pData->u[7] = input_->c_hat(pData->u[7])*input_->x0();
             pData->phi[2] = input_->c_hat(pData->phi[2])*input_->x0();
             pData->n[2] = input_->V_hat(pData->n[2])*input_->x0()*input_->x0();
             pData->p[2] = input_->c_hat(pData->p[2])*input_->x0()*input_->x0();
             pData->u[11] = input_->c_hat(pData->u[11])*input_->x0()*input_->x0();
             */
        }
    }

    /**
     * This functions determines the cretia for stopping the convergence loop over the functions.
     * Each step computes L2Error of previous srep and current step. In main by setting:
     * converge_number we will compare L2Error of each step with current step.
     * if L2Error is less than converge_number loop will be stopped.
     */
    ZEROMATRIX<double> Error(const DDInputData *input_) {
        int n_tot = p_grid_->n_nodes();
        ZEROMATRIX<double> nc;
        nc.redim(n_tot, 3);
        for (int nodeID = 0; nodeID < n_tot; nodeID++) {
            DDNodeData *pData = &(GetNodeData(nodeID));
            nc(nodeID, 0) = pData->n[0];
        }
        for (int nodeID = 0; nodeID < n_tot; nodeID++) {
            DDNodeData *pData = &(GetNodeData(nodeID));
            nc(nodeID, 1) = pData->p[0];
        }
        for (int nodeID = 0; nodeID < n_tot; nodeID++) {
            DDNodeData *pData = &(GetNodeData(nodeID));
            nc(nodeID, 2) = pData->x[0];
        }
        return nc;
    }

};

#endif

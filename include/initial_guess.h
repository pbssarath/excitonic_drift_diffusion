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
#ifndef INCLUDE_INITIALGUESS_H_
#define INCLUDE_INITIALGUESS_H_

#include <string>               // for error string on construction


/**
 * Class calculates the initial value at a node.
 *
 * The values are calculated as the following:
 * 1 -> v = -Y + 0.5
 * 2 -> n = exp(-60 Y)
 * 3 -> p = exp(60(Y-1))
 * 4 -> x = 0
 * * Usage:
 * InitialGuess sol(n_dimensions, p_grid);
 * double intial_value = sol.ValueAt(node_id);  
 *
 * Class is also responsible too find valueFEM and ValueDerivativeFEM for Hermite Basis function
 *
 * Usage:
 * InitialGuess sol(n_dimensions, p_grid);
 * double intial_value = sol.ValueAt(node_id);  
 */
class InitialGuess {
public:
    InitialGuess(int n_dimensions, GRID *p_grid, GridField<DDNodeData> *p_data) : n_dimensions_(n_dimensions),
                                                                                  p_grid_(p_grid), p_data_(p_data) {
        if (n_dimensions_ < 1 || n_dimensions_ > 3) {
            throw (std::string("Invalid number of dimensions in initial_guess."));
        }
    }

    ~InitialGuess() {
    }

    /**
     * Returns the value at the given node at time zero.
     *
     * @param node_id ID of node where the value will be calculated.
     * @return value of initial guess for the given node
     */
    double ValueAt(int node_id, int n_dof) const {
        const ZEROPTV &pNode = p_grid_->GetNode(node_id)->location();  // node we care about

        double value = 1.0;
        double scale = 0.0258520269359094;
        switch (n_dof % 12) {
            case 0: {
                value *= 1. / scale * (0.5 - pNode(1));
                break;
            }
            case 1: {
                value *= exp(-60 * pNode(1));
                break;
            }
            case 2: {
                value *= exp(60 * (pNode(1) - 1));
                break;
            }
            case 3: {
                value *= 0;
                break;
            }
        }
        return value;
    }

    /**
     * The library has a built-in function GridField::valueFEM that is used to calculate the interpolated
     * node value at the gauss point. However, the library's implementation assumes that the hermite
     * derivatives are stored contiguously (at index + 1, index + 2, index +3, etc.) This code stores
     * the hermite derivatives in a stride (at index + 4*1, index + 4*2, index + 3*4, etc.), so we have
     * our own valueFEM function to interpolate values when using hermite derivatives.
     * @param index index of the variable to interpolate
     * @param dof the total number of degrees of freedom
     * @returns variable[index] at fe.position()
     */
    double valueFEM(const FEMElm &fe, int index, int dof) const {
        double sum = 0.0;
        if (fe.basis_function() == BASIS_HERMITE) {
            // Hermite is a special case - some basis functions are derivatives.
            // This assumes the node data is stored such that u = index
            // and du_x = index + 1, du_y = index + 2, etc.
            // 1D: -, x
            // 2D: -, x, y, xy
            // 3D: -, x, y, z, xy, xz, yz, xyz
            const int nbf = fe.nbf();
            const int nbf_per_node = fe.nbf_per_node();
            for (int a = 0; a < nbf; a++) {
                const int node_idx = fe.elem()->ElemToLocalNodeID(a / nbf_per_node);
                const int deriv = a % nbf_per_node;

                // DIFFERENT FROM LIBRARY: deriv is multiplied by dof!
                sum += fe.N(a) * this->p_data_->GetNodeData(node_idx).value(index + deriv * dof);
            }
        } else {
            const int nbf = fe.nbf();
            for (ElemNodeID a = 0; a < nbf; a++) {
                const int node_index = fe.elem()->ElemToLocalNodeID(a);
                sum += fe.N(a) * this->p_data_->GetNodeData(node_index).value(index);
            }
        }
        return sum;
    }

    /**
     * The library has a built-in function GridField::valueDerivativeFEM that is used to calculate the interpolated
     * node derivative values at the gauss point. However, the library's implementation assumes that the hermite
     * derivatives are stored contiguously (at index + 1, index + 2, index +3, etc.) This code stores
     * the hermite derivatives in a stride (at index + 4*1, index + 4*2, index + 3*4, etc.), so we have
     * our own valueDerivativeFEM function to interpolate values when using hermite derivatives.
     * @param index index of the variable to interpolate
     * @param dof the total number of degrees of freedom
     * @returns derivative of variable[index] at fe.position()
     */
    double valueDerivativeFEM(const FEMElm &fe, int index, int dir, int dof) const {
        double sum = 0.0;

        if (fe.basis_function() == BASIS_HERMITE) {
            const int nbf = fe.nbf();
            const int nbf_per_node = fe.nbf_per_node();
            for (int a = 0; a < nbf; a++) {
                const int node_idx = fe.elem()->ElemToLocalNodeID(a / nbf_per_node);
                const int deriv = a % nbf_per_node;
                sum += fe.dN(a, dir) * this->p_data_->GetNodeData(node_idx).value(index + deriv * dof);
            }
        } else {
            const int nbf = fe.nbf();
            for (ElemNodeID a = 0; a < nbf; a++) {
                const int node_index = fe.elem()->ElemToLocalNodeID(a);
                sum += fe.dN(a, dir) * this->p_data_->GetNodeData(node_index).value(index);
            }
        }
        return sum;
    }

    double dot(ZEROPTV &idx_1, ZEROPTV &idx_2, int nsd = 2) const {
        double sum = 0.0;
        if (nsd == 2) {
            sum += idx_1(0) * idx_2(0) + idx_1(1) * idx_2(1);
        } else if (nsd == 3) {
            sum += idx_1(0) * idx_2(0) + idx_1(1) * idx_2(1) + idx_1(2) * idx_2(2);
        }
        return (sum);
    }


private:
    int n_dimensions_;            ///< number of spatial dimensions of system (1, 2, or 3)
    GRID *p_grid_;                ///< pointer to grid
    GridField<DDNodeData> *p_data_;  ///< pointer to grid field data
};

#endif // INCLUDE_HTANALYTICSOLUTION_H_

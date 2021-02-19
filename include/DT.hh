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
#include "math.h"
#include "DDGridField.h"
/**
* This class solves level-set equation to get the distance contour.
* Simple finite difference is used.
* If you know nothing about level-set method, check out this link:
* http://profs.etsmtl.ca/hlombaert/levelset/
*/
class DistContour {

private:
    double dtdn[3];
    int dims[3];
    float threshold = 0.5;

    inline double pos_sq(double x);

    inline double neg_sq(double x);

    inline double max(double a, double b);

    inline int index(int i, int j, int k, int *dims, bool isPeriodic);

    double norm_grad_2d(double *lsf, double *sgn, int i, int j, int *dims, bool isPeriodic);

    double norm_grad_3d(double *lsf, double *sgn, int i, int j, int k, int *dims, bool isPeriodic);

    void gradient(int *dims, int ndims, double band, bool isPeriodic);

public:
    double *dt;
    double *dtdx;
    double *dtdy;
    double *dtdz;

    void set_threshold(float);

    void printMsg();

    DistContour();

    int construct(double *lsf, int *dims, int ndims, double band, double tolerance, bool isPeriodic);

    void DistContour2d(Mat2d &coarseMorph, int i, int j, int k, bool isPeriodic, double Geo2Phy);

    void DistContour3d(Mat3d &coarseMorph, int i, int j, int k, bool isPeriodic, float Geo2Phy);
};

void DistContour::set_threshold(float threshold){
  this->threshold = threshold;
}

void DistContour::printMsg () {
  std::cout << "Distance contour has been calculated" << endl;
}

DistContour::DistContour () {
}

/**
 * Set the 2d distance contour.
 *
 * @param coarseMorph grid field structure
 * @param i number of nodes along x direction -- usually Nelemx + 1
 * @param j number of nodes in y direction -- using Nelemy + 1
 * @param k number of nodes in z direction -- this will be '1' for 2D problem
 * @param isPeriodic is the problem periodic or not
 * @param Geo2Phy coefficient to convert distance contour to real physical distance : constant
 */
void DistContour::DistContour2d (Mat2d & coarseMorph, int i, int j, int k, bool isPeriodic, double Geo2Phy) {
    int nodeID;
    dims[0] = i;                  //node number along each direction; mesh#+1 !
    dims[1] = j;
    dims[2] = k;
    int ndims = 2;
    int nodeNum = i * j * k;
    this->dt = new double[nodeNum];
    this->dtdx = new double[nodeNum];
    this->dtdy = new double[nodeNum];
    for (int jindex = 0; jindex < j; jindex++) {
        for (int iindex = 0; iindex < i; iindex++) {
            nodeID = iindex + jindex * i;
            double msNode = coarseMorph[jindex][iindex];
            if (msNode == threshold) {
                this->dt[nodeID] = 0;
            } else if (msNode >= threshold) {
                this->dt[nodeID] = 1;
            } else {
                this->dt[nodeID] = -1;
            }
        }
    }

    construct(this->dt, this->dims, ndims, 50, 1e-6, isPeriodic);
#ifdef NDEBUG
    if(!GetMPIRank()) {
        std::fstream fp_dt;
        fp_dt.open("debug_dt.plt", std::fstream::out);
        for(int jj = 0; jj < this->dims[1]; jj++){
            for(int ii = 0; ii < this->dims[0]; ii++){
                fp_dt << ii << "\t" << jj << "\t" << this->dt[jj * i + ii] << "\n";
            }
        }
        fp_dt.close();
        std::cout << "DT written to debug_dt.plt" << std::endl;
    }
#endif
    gradient(this->dims, ndims, 50, isPeriodic);
    for (nodeID = 0; nodeID < nodeNum; nodeID++) {
        this->dt[nodeID] *= Geo2Phy;
        //this->dtdx[nodeID] /= Geo2Phy;
        //this->dtdy[nodeID] /= Geo2Phy;
    }
    PrintStatus("Distance contour calculated");
}

/**
* Set the 3d distance contour.
*
* @params: 1, grid field structure; 2-4, node# along x, y, z direction; where the problem is periodic;
* coefficient, which is a constant number, to convert distance contour to real physical distance.
*/

void DistContour::DistContour3d (Mat3d & coarseMorph, int i, int j, int k, bool isPeriodic, float Geo2Phy) {
    int nodeID;
    dims[0] = i;                  //node number along each direction; mesh#+1 !
    dims[1] = j;
    dims[2] = k;
    int ndims = 3;
    int nodeNum = i * j * k;
    this->dt = new double[nodeNum];
    this->dtdx = new double[nodeNum];
    this->dtdy = new double[nodeNum];
    this->dtdz = new double[nodeNum];
    for (int kindex = 0; kindex < k; kindex++) {
        for (int jindex = 0; jindex < j; jindex++) {
            for (int iindex = 0; iindex < i; iindex++) {
                nodeID = iindex + jindex * i + kindex * (i * j);
                float msNode = coarseMorph[iindex][jindex][kindex];
                if (msNode == threshold) {
                    this->dt[nodeID] = 0;
                } else if (msNode > threshold) {
                    this->dt[nodeID] = 1;
                } else {
                    this->dt[nodeID] = -1;
                }
            }
        }
    }
    construct(this->dt, this->dims, ndims, 50, 5e-3, isPeriodic);
    gradient(this->dims, ndims, 50, isPeriodic);
    for (nodeID = 0; nodeID < nodeNum; nodeID++) {
        this->dt[nodeID] *= Geo2Phy;
    }
    PrintStatus("Distance contour calculated");
}

inline double DistContour::pos_sq (double x) {
    return x > 0 ? x * x : 0;
}

inline double DistContour::neg_sq (double x) {
    return x < 0 ? x * x : 0;
}

inline double DistContour::max (double a, double b) {
    return a >= b ? a : b;
}

inline int DistContour::index (int i, int j, int k, int *dims, bool isPeriodic) {
    if (!isPeriodic) {
        if (i == dims[0])
            i = i - 1;
        if (i == -1)
            i = 0;
    } else {
        if (i == dims[0])
            i = 0;
        if (i == -1)
            i = dims[0] - 1;
    }
    if (j == dims[1])
        j = j - 1;                  //j = 0;
    if (j == -1)
        j = 0;                      // dims[1] - 1;
    //if (k > 1)
    {                             //only needed for 3d.
        if (k == dims[2])
            k = k - 1;                //k = 0;
        if (k == -1)
            k = 0;                    // dims[2] - 1;
    }
    return k * dims[1] * dims[0] + j * dims[0] + i;
}

// Calculate the gradient in 3D.
double DistContour::norm_grad_3d (double *lsf, double *sgn, int i, int j, int k, int *dims, bool isPeriodic) {
    unsigned p = index(i, j, k, dims, isPeriodic);
    double a = lsf[p] - lsf[index(i - 1, j, k, dims, isPeriodic)];
    double b = lsf[index(i + 1, j, k, dims, isPeriodic)] - lsf[p];
    double c = lsf[p] - lsf[index(i, j - 1, k, dims, isPeriodic)];
    double d = lsf[index(i, j + 1, k, dims, isPeriodic)] - lsf[p];
    double e = lsf[p] - lsf[index(i, j, k - 1, dims, isPeriodic)];
    double f = lsf[index(i, j, k + 1, dims, isPeriodic)] - lsf[p];

    double grad_x, grad_y, grad_z;
    if (sgn[p] > 0) {
        grad_x = max(pos_sq(a), neg_sq(b));
        grad_y = max(pos_sq(c), neg_sq(d));
        grad_z = max(pos_sq(e), neg_sq(f));
    } else {
        grad_x = max(pos_sq(b), neg_sq(a));
        grad_y = max(pos_sq(d), neg_sq(c));
        grad_z = max(pos_sq(f), neg_sq(e));
    }

    return sqrt(grad_x + grad_y + grad_z);
}

// Calculate the gradient in 2D.
double DistContour::norm_grad_2d (double *lsf, double *sgn, int i, int j, int *dims, bool isPeriodic) {
    int k = 0;
    unsigned p = index(i, j, k, dims, isPeriodic);

    double a = lsf[p] - lsf[index(i - 1, j, k, dims, isPeriodic)];
    double b = lsf[index(i + 1, j, k, dims, isPeriodic)] - lsf[p];

    double c = lsf[p] - lsf[index(i, j - 1, k, dims, isPeriodic)];
    double d = lsf[index(i, j + 1, k, dims, isPeriodic)] - lsf[p];

    double grad_x, grad_y;

    if (sgn[p] > 0) {
        grad_x = max(pos_sq(a), neg_sq(b));
        grad_y = max(pos_sq(c), neg_sq(d));
    } else {
        grad_x = max(pos_sq(b), neg_sq(a));
        grad_y = max(pos_sq(d), neg_sq(c));
    }
    return sqrt(grad_x + grad_y);
}

void DistContour::gradient(int *dims, int ndims, double band, bool isPeriodic) {
    if (ndims == 2) {
        for (int i = 0; i < dims[0]; i++)
            for (int j = 0; j < dims[1]; j++) {
                int ndx = j * dims[0] + i;
                int ndxright = j * dims[0] + i + 1;
                int ndxleft = j * dims[0] + i - 1;
                int ndxtop = (j + 1) * dims[0] + i;
                if (j == 0 || (i == dims[0] - 1) || (j == dims[1] - 1)) {
                    dtdx[ndx] = 0;
                    dtdy[ndx] = (dt[ndxtop] - dt[ndx]);
                } else {
                    dtdx[ndx] = (dt[ndxright] - dt[ndx]);
                    dtdy[ndx] = (dt[ndxtop] - dt[ndx]);
                }
                if (i == dims[0] - 1) {  //set right boundary
                    if (isPeriodic) dtdx[ndx] = dt[0] - dt[ndx];
                    else dtdx[ndx] = dt[ndx] - dt[ndxleft];
                }
            }
    } else if (ndims == 3) {
        for (int i = 0; i < dims[0]; i++)
            for (int j = 0; j < dims[1]; j++)
                for (int k = 0; k < dims[2]; k++) {
                    int ndx = k * dims[1] * dims[0] + j * dims[0] + i;
                    int ndxright = k * dims[1] * dims[0] + j * dims[0] + i + 1;
                    int ndxtop = k * dims[1] * dims[0] + (j + 1) * dims[0] + i;
                    int ndxback = (k - 1) * dims[1] * dims[0] + j * dims[0] + i;
                    if (j == 0 || k == 0 || (i == dims[0] - 1) || (j == dims[1] - 1) || (k == dims[2] - 1)) {
                        dtdx[ndx] = 0;
                        dtdy[ndx] = 0;
                        dtdz[ndx] = 0;
                    } else {
                        dtdx[ndx] = (dt[ndxright] - dt[ndx]);
                        dtdy[ndx] = (dt[ndxtop] - dt[ndx]);
                        dtdz[ndx] = (dt[ndx] - dt[ndxback]);
                    }
                }

    }
}



int DistContour::construct (double *lsf, int *dims, int ndims, double band, double tolerance, bool isPeriodic) {
    double dt = 0.2;
    double max_dt = tolerance + 1;
    double *sign, *sdf_dt;
    int sdf_dims[3];

    sdf_dims[0] = dims[0];
    sdf_dims[1] = dims[1];
    if (ndims == 2)
        sdf_dims[2] = 1;
    else if (ndims == 3)
        sdf_dims[2] = dims[2];
    sign = new double[dims[0] * dims[1] * dims[2]];
    sdf_dt = new double[dims[0] * dims[1] * dims[2]];
    for (int i = 0; i < dims[0]; i++)
        for (int j = 0; j < dims[1]; j++)
            for (int k = 0; k < dims[2]; k++) {
                int ndx = k * dims[1] * dims[0] + j * dims[0] + i;
                if (lsf[ndx] > 0.9)
                    sign[ndx] = 1;
                else if (lsf[ndx] < 0.1)
                    sign[ndx] = -1;
                else
                    sign[ndx] = 0;
            }
    int iter = 0;
    while (max_dt > tolerance && iter < 1000) {
        max_dt = 0.0;
        for (int i = 0; i < dims[0]; i++)
            for (int j = 0; j < dims[1]; j++)
                for (int k = 0; k < dims[2]; k++) {
                    int ndx = k * dims[1] * dims[0] + j * dims[0] + i;
                    if (fabs(lsf[ndx]) >= band) {
                        sdf_dt[ndx] = 0;
                    } else {
                        double s, grad;

                        s = lsf[ndx] / sqrt(lsf[ndx] * lsf[ndx] + 1.0);

                        if (ndims == 2) {
                            grad = norm_grad_2d(lsf, sign, i, j, sdf_dims, isPeriodic);
                        } else if (ndims == 3) {
                            grad = norm_grad_3d(lsf, sign, i, j, k, sdf_dims, isPeriodic);
                        } else {
                            return -1;
                        }
                        sdf_dt[ndx] = s * (1.0 - grad);
                        if (fabs(sdf_dt[ndx]) > max_dt)
                            max_dt = fabs(sdf_dt[ndx]);
                    }
                }

        for (int i = 0; i < dims[0]; i++)
            for (int j = 0; j < dims[1]; j++)
                for (int k = 0; k < dims[2]; k++) {
                    int ndx = k * dims[1] * dims[0] + j * dims[0] + i;
                    lsf[ndx] += dt * sdf_dt[ndx];
                }
        iter++;
    }

    delete[]sign;
    delete[]sdf_dt;
    return 0;
}

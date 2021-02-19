/*
  Copyright 2014-2015 Baskar Ganapathysubramanian

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
#ifndef DD_NODEDATA_HPP
#define DD_NODEDATA_HPP

// Global variables define a meaningful name for the index.
enum NodeDataIndices:int {
    G_ID = 44,                          // Generation
    DT_ID = 39,                         // Distance contour
    MUN_ID = 40,                        // Electron mobility
    MUP_ID = 41,                        // Hole mobility
    EPS_ID = 45,                        // Epsilon
    PHIBAR_ID = 49,                     // Interpotential potential
    D2PHIBAR_ID = 50,                   // Interpotential potential derivative
    N_COUPLED = 3,                      // Number of variables coupled together
    PHI_ID = 0,                         // Electric field
    N_ID = 1,                           // Electron
    P_ID = 2,                           // Hole
    X_ID = 24,                          // Exciton
    JNY_ID = 34,                        // JNy
    JPY_ID = 36,                        // JPy
    JNX_ID = 35,                        // JNx
    JPX_ID = 37,                        // JPX
    GRADJ_ID = 14,                      // Gradient of Jn and Jp
    change = 10,
    dDTdzID = 48,
};
class DDNodeData:public NODEData {

public:
    static int nsd;
    static int couple;
    double phi[8];
    double n[8];
    double p[8];
    double x[8];
    double j[9];
    double ms;
    double dt;
    double GenX;
    double epsn;
    double mun;
    double mup;
    double Dis;
    double R;
    double PHIBAR;
    double D2PHIBAR;

    virtual double &value(int index) {
      if (index < 24) {
        switch (index % N_COUPLED) {
          case 0: {
            return phi[int(index / N_COUPLED)];
            break;
          }
          case 1: {
            return n[int(index / N_COUPLED)];
            break;
          }
          case 2: {
            return p[int(index / N_COUPLED)];
            break;
          }
        }
      } else if (index > 23 && index < 32) {
        return x[index - 24];
      } else if (index == 16 + 16) {
        return Dis;
      } else if (index == 17 + 16) {
        return R;
      } else if (index == 18 + 16) {
        return j[0];
      } else if (index == 19 + 16) {
        return j[1];
      } else if (index == 20 + 16) {
        return j[2];
      } else if (index == 21 + 16) {
        return j[3];
      } else if (index == 22 + 16) {
        return ms;
      } else if (index == 23 + 16) {
        return dt;
      } else if (index == 24 + 16) {
        return mun;
      } else if (index == 25 + 16) {
        return mup;
      } else if (index == 26 + 16) {
        return j[4];
      } else if (index == 27 + 16) {
        return j[5];
      } else if (index == 28 + 16) {
        return GenX;
      } else if (index == 29 + 16) {
        return epsn;
      } else if (index == 30 + 16) {
        return j[6];
      } else if (index == 31 + 16) {
        return j[7];
      } else if (index == 32 + 16) {
        return j[8];
      } else if (index == 33 + 16) {
        return PHIBAR;
      } else if (index == 34 + 16) {
        return D2PHIBAR;
      } else {
        throw TALYException() << "Invalid DDNodeData index";
      }
      return PHIBAR;
    }

    virtual const double &value(int index) const {
      if (index < 24) {
        switch (index % N_COUPLED) {
          case 0: {
            return phi[int(index / N_COUPLED)];
            break;
          }
          case 1: {
            return n[int(index / N_COUPLED)];
            break;
          }
          case 2: {
            return p[int(index / N_COUPLED)];
            break;
          }
        }
      } else if (index > 23 && index < 32) {
        return x[index - 24];
      } else if (index == 16 + 16) {
        return Dis;
      } else if (index == 17 + 16) {
        return R;
      } else if (index == 18 + 16) {
        return j[0];
      } else if (index == 19 + 16) {
        return j[1];
      } else if (index == 20 + 16) {
        return j[2];
      } else if (index == 21 + 16) {
        return j[3];
      } else if (index == 22 + 16) {
        return ms;
      } else if (index == 23 + 16) {
        return dt;
      } else if (index == 24 + 16) {
        return mun;
      } else if (index == 25 + 16) {
        return mup;
      } else if (index == 26 + 16) {
        return j[4];
      } else if (index == 27 + 16) {
        return j[5];
      } else if (index == 28 + 16) {
        return GenX;
      } else if (index == 29 + 16) {
        return epsn;
      } else if (index == 30 + 16) {
        return j[6];
      } else if (index == 31 + 16) {
        return j[7];
      } else if (index == 32 + 16) {
        return j[8];
      } else if (index == 33 + 16) {
        return PHIBAR;
      } else if (index == 34 + 16) {
        return D2PHIBAR;
      } else {
        throw TALYException() << "Invalid DDNodeData index";
      }
      return PHIBAR;
    }


    static char *name(int index) {
      static char str[256];
      if (index == PHI_ID) {
        snprintf(str, 256, "phi");
      } else if (index == N_ID) {
        snprintf(str, 256, "n");
      } else if (index == P_ID) {
        snprintf(str, 256, "p");
      } else if (index == 3) {
        snprintf(str, 256, "phi_x");
      } else if (index == 4) {
        snprintf(str, 256, "n_x");
      } else if (index == 5) {
        snprintf(str, 256, "p_x");
      } else if (index == 6) {
        snprintf(str, 256, "phi_y");
      } else if (index == 7) {
        snprintf(str, 256, "n_y");
      } else if (index == 8) {
        snprintf(str, 256, "p_y");
      } else if (index == 9) {
        snprintf(str, 256, "phi_xy");
      } else if (index == 10) {
        snprintf(str, 256, "n_xy");
      } else if (index == 11) {
        snprintf(str, 256, "p_xy");
      } else if (index == 9) {
        if (nsd == 3) {
          snprintf(str, 256, "phi_z");
        } else if (nsd == 2) {
          snprintf(str, 256, "phi_xy_2D");
        }
      } else if (index == 10) {
        if (nsd == 3) {
          snprintf(str, 256, "n_Z");
        } else if (nsd == 2) {
          snprintf(str, 256, "n_xy_2D");
        }
      } else if (index == 11) {
        if (nsd == 3) {
          snprintf(str, 256, "p_z");
        } else if (nsd == 2) {
          snprintf(str, 256, "p_xy_2D");
        }
      } else if (index == 12) {
        snprintf(str, 256, "phi_xy_3d");
      } else if (index == 13) {
        snprintf(str, 256, "n_xy_3d");
      } else if (index == 14) {
        snprintf(str, 256, "p_xy_3d");
      } else if (index == 15) {
        snprintf(str, 256, "phi_xz_3D");
      } else if (index == 16) {
        snprintf(str, 256, "n_xz_3D");
      } else if (index == 17) {
        snprintf(str, 256, "p_xz_3D");
      } else if (index == 18) {
        snprintf(str, 256, "phi_yz_3D");
      } else if (index == 19) {
        snprintf(str, 256, "n_yz_3D");
      } else if (index == 20) {
        snprintf(str, 256, "p_yz_3D");
      } else if (index == 21) {
        snprintf(str, 256, "phi_xyz_3D");
      } else if (index == 22) {
        snprintf(str, 256, "n_xyz_3D");
      } else if (index == 23) {
        snprintf(str, 256, "p_xyz_3D");
      } else if (index == 24) {
        snprintf(str, 256, "Exciton");
      } else if (index == 25) {
        snprintf(str, 256, "X_x");
      } else if (index == 26) {
        snprintf(str, 256, "X_y");
      } else if (index == 27) {
        snprintf(str, 256, "X_z");
      } else if (index == 28) {
        snprintf(str, 256, "X_xy");
      } else if (index == 29) {
        snprintf(str, 256, "X_xz");
      } else if (index == 30) {
        snprintf(str, 256, "X_yz");
      } else if (index == 31) {
        snprintf(str, 256, "X_xyz");
      } else if (index == 32) {
        snprintf(str, 256, "Dis");
      } else if (index == 33) {
        snprintf(str, 256, "R");
      } else if (index == 34) {
        snprintf(str, 256, "Jn_y");
      } else if (index == 35) {
        snprintf(str, 256, "Jn_x");
      } else if (index == 36) {
        snprintf(str, 256, "Jp_y");
      } else if (index == 37) {
        snprintf(str, 256, "Jp_x");
      } else if (index == 38) {
        snprintf(str, 256, "MS");
      } else if (index == 39) {
        snprintf(str, 256, "DT");
      } else if (index == 40) {
        snprintf(str, 256, "mu_n");
      } else if (index == 41) {
        snprintf(str, 256, "mu_p");
      } else if (index == 42) {
        snprintf(str, 256, "jn_post");
      } else if (index == 43) {
        snprintf(str, 256, "jp_post");
      } else if (index == 44) {
        snprintf(str, 256, "G");
      } else if (index == 45) {
        snprintf(str, 256, "epsilon");
      } else if (index == 46) {
        snprintf(str, 256, "|grad_jn|");
      } else if (index == 47) {
        snprintf(str, 256, "|grad_jp|");
      } else if (index == 48) {
        snprintf(str, 256, "relativeErrorOfGradJ");
      } else if (index == 49) {
        snprintf(str, 256, "PHIBAR");
      } else if (index == 50) {
        snprintf(str, 256, "D2PHIBAR");
      } else {
        throw TALYException() << "Unknown node variable '" << index << "'";
      }

      return str;
    }

    static int valueno() {
      return 51;
    }
};




#endif

#ifndef GEN_DISS_H
#define GEN_DISS_H

#include "DDInputData.h"
class gen_diss {
public:
    DDInputData *inputData;

    gen_diss() {
    };

    ~gen_diss() {
    };

    double calcKdiss(double y_h, double F_hat, double DT, double epsilong, double mu_dim);

    double calcR(double y_h, double n_hat, double p_hat, double DT, double epsilong, double mu_dim);
};

#endif

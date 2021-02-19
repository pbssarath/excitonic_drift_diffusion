#include "gen_diss.hh"

double gen_diss::calcKdiss (double y_hat, double F_hat, double DT, double epsilong, double mu_dim) {
    double kdiss = 0;
    if ((DT > inputData->interfaceThk / 2.) || (DT < -inputData->interfaceThk / 2.)) {
        // we are far from the interface => no dissociation i.e. kdiss = 0.0;
    } else {
        // mu_dim = inputData->mu_n + inputData->mu_p;
        double F = F_hat / inputData->x0() * inputData->V0();
        double gammar = inputData->q * mu_dim / epsilong; //inputData->epsilon(ms,DT,ycoord);
        double q = inputData->q;
        double a = inputData->a;
        double k_B = inputData->k_B;
        double T = inputData->T;
        double ia = (a < 1e-15 ? 0 : 1. / a);
        double E_B = q * q * ia / (4 * M_PI * epsilong);
        double b = q * q * q * F / (8 * M_PI * epsilong * k_B * k_B * T * T);
        kdiss = (3 * gammar * ia * ia * ia / (4 * M_PI)) * exp(-E_B / (k_B * T))
                * (1 + b + b * b / 3 + b * b * b / 18 + b * b * b * b / 180 + b * b * b * b * b / 2700);
    }

    if (kdiss < 0)
        kdiss = 0;                  //does not allow negative dissociation.
    if (y_hat > 0.92 || y_hat < 0.08)
        kdiss = 0;

    return kdiss;
}


//recombination rate:
double gen_diss::calcR (double y_hat, double n_hat, double p_hat, double DT, double epsilong, double mu_dim) {
    double R = 0;
    double q = inputData->q;
    double n = inputData->c_unhat(n_hat);
    double p = inputData->c_unhat(p_hat);
    double ni = inputData->n_int(1., DT);
    double gamma = q * (inputData->mu_n + inputData->mu_p) / epsilong;
    R = gamma * (n * p - 0 * ni * ni);
    if (R < 0)
        R = 0;
    if (y_hat > 0.92 || y_hat < 0.08)
        R = 0;
    return R;
}

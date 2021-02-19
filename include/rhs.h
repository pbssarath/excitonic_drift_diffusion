#ifndef INCLUDE_RHS_H_
#define INCLUDE_RHS_H_

#include <string>               // for error string on construction
/**
 * Class calculates and returns value of Dissosiation Recombination as well as their derivatives.
 *
 * The values of Dissosiation are calculated as the following:
 * D = K_diss * X
 * The values of Recombination are calculated as the following:
 * R = gamma * n * p
 *  
 */
class RhS
{
public:
  RhS (int n_dimensions):n_dimensions_ (n_dimensions)
  {
    if (n_dimensions_ < 1 || n_dimensions_ > 3)
    {
      throw (std::string ("Invalid number of dimensions in initial_guess."));
    }
  }

   ~RhS ()
  {
  }

  /**
   * Returns the value of the non dimensionalized form of dissosiation for be[j].
   *
   * @param non_dim_fact non dimensionalized factore for dissosiation.
   * @return value of dissosiation for the given j
   */
  double CalcDiss (const FEMElm & fe, double k_diss, double x_pre, int j) const
  {
    double value = k_diss * x_pre * fe.N (j);
      return value;
  }
  /**
   * Returns the value of dissosiation derivative for Ae[j,k].
   *
   * @param non_dim_fact non dimensionalized factore for dissosiation.
   * @return value of derivative of dissosiation for the given j and k
   */
  double CalcDissDer (const FEMElm & fe, double k_diss, int j, int k) const
  {
    double value = k_diss * fe.N (j);
      return value;
  }
 /**
   * Returns the value of the non dimensionalized form of Recombination for be[j].
   *
   * @param gamma recombination coefficient.
   * @return value of Recombination for the given j
   */
  double CalcRecom (const FEMElm & fe, double gamma, double n_elec_pre, double p_hole_pre, int j) const
  {
    double value = gamma * n_elec_pre * p_hole_pre;
      return value;
  }

 /**
   * Returns the value of the Recombination derivative for Ae[j,k].
   *
   * @param gamma recombination coefficient.
   * @return value of Recombination derivative for the given j and k
   */
/**  double CalcRecomDer(const FEMElm& fe, double gamma, int j) const {
    //double non_dim_fact = 1.287 * pow(10,-6);
    //double value = k_diss * non_dim_fact;
    //return value;
  }
*/
private:
  int n_dimensions_;            ///< number of spatial dimensions of system (1, 2, or 3)

};

#endif // INCLUDE_HTANALYTICSOLUTION_H_

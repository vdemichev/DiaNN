#ifndef OUTPUT_REGRESSIONMSE_H_
#define OUTPUT_REGRESSIONMSE_H_

#include <Eigen/Core>
#include <stdexcept>
#include "../Config.h"

namespace MiniDNN
{


///
/// \ingroup Outputs
///
/// Regression output layer using Mean Squared Error (MSE) criterion
///
class RegressionMSE: public Output
{
    private:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;

        Matrix m_din;  // Derivative of the input of this layer.
        // Note that input of this layer is also the output of previous layer

    public:
        void evaluate(const Matrix& prev_layer_data, const Matrix& target)
        {
            // Check dimension
            const int nobs = prev_layer_data.cols();
            const int nvar = prev_layer_data.rows();

            if ((target.cols() != nobs) || (target.rows() != nvar))
            {
                throw std::invalid_argument("[class RegressionMSE]: Target data have incorrect dimension");
            }

            // Compute the derivative of the input of this layer
            // L = 0.5 * ||yhat - y||^2
            // in = yhat
            // d(L) / d(in) = yhat - y
            m_din.resize(nvar, nobs);
            m_din.noalias() = prev_layer_data - target;
        }

        const Matrix& backprop_data() const
        {
            return m_din;
        }

        Scalar loss() const
        {
            // L = 0.5 * ||yhat - y||^2
            return m_din.squaredNorm() / m_din.cols() * Scalar(0.5);
        }

        std::string output_type() const
        {
            return "RegressionMSE";
        }
};


} // namespace MiniDNN


#endif /* OUTPUT_REGRESSIONMSE_H_ */

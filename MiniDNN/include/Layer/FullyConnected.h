#ifndef LAYER_FULLYCONNECTED_H_
#define LAYER_FULLYCONNECTED_H_

#include <Eigen/Core>
#include <vector>
#include <stdexcept>
#include "../Config.h"
#include "../Layer.h"
#include "../Utils/Random.h"
#include "../Utils/IO.h"
#include "../Utils/Enum.h"

namespace MiniDNN
{


///
/// \ingroup Layers
///
/// Fully connected hidden layer
///
template <typename Activation>
class FullyConnected: public Layer
{
    private:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
        typedef Vector::ConstAlignedMapType ConstAlignedMapVec;
        typedef Vector::AlignedMapType AlignedMapVec;
        typedef std::map<std::string, int> MetaInfo;

        Matrix m_weight;  // Weight parameters, W(in_size x out_size)
        Vector m_bias;    // Bias parameters, b(out_size x 1)
        Matrix m_dw;      // Derivative of weights
        Vector m_db;      // Derivative of bias
        Matrix m_z;       // Linear term, z = W' * in + b
        Matrix m_a;       // Output of this layer, a = act(z)
        Matrix m_din;     // Derivative of the input of this layer.
                          // Note that input of this layer is also the output of previous layer

    public:
        ///
        /// Constructor
        ///
        /// \param in_size  Number of input units.
        /// \param out_size Number of output units.
        ///
        FullyConnected(const int in_size, const int out_size) :
            Layer(in_size, out_size)
        {}

        void init(const Scalar& mu, const Scalar& sigma, RNG& rng)
        {
            // Set parameter dimension
            init();
            // Set random coefficients
            internal::set_normal_random(m_weight.data(), m_weight.size(), rng, mu, sigma);
            internal::set_normal_random(m_bias.data(), m_bias.size(), rng, mu, sigma);
        }

        void init()
        {
            // Set parameter dimension
            m_weight.resize(this->m_in_size, this->m_out_size);
            m_bias.resize(this->m_out_size);
            m_dw.resize(this->m_in_size, this->m_out_size);
            m_db.resize(this->m_out_size);

			m_z.resize(0, 0);
			m_a.resize(0, 0);
			m_din.resize(0, 0);
        }

        // prev_layer_data: in_size x nobs
        void forward(const Matrix& prev_layer_data)
        {
            const int nobs = prev_layer_data.cols();
            // Linear term z = W' * in + b
			m_z.resize(this->m_out_size, nobs);
            m_z.noalias() = m_weight.transpose() * prev_layer_data;
            m_z.colwise() += m_bias;
            // Apply activation function
            m_a.resize(this->m_out_size, nobs);
            Activation::activate(m_z, m_a);
        }

        const Matrix& output() const
        {
            return m_a;
        }

        // prev_layer_data: in_size x nobs
        // next_layer_data: out_size x nobs
        void backprop(const Matrix& prev_layer_data, const Matrix& next_layer_data)
        {
            const int nobs = prev_layer_data.cols();
            // After forward stage, m_z contains z = W' * in + b
            // Now we need to calculate d(L) / d(z) = [d(a) / d(z)] * [d(L) / d(a)]
            // d(L) / d(a) is computed in the next layer, contained in next_layer_data
            // The Jacobian matrix J = d(a) / d(z) is determined by the activation function
            Matrix& dLz = m_z;
            Activation::apply_jacobian(m_z, m_a, next_layer_data, dLz);
            // Now dLz contains d(L) / d(z)
            // Derivative for weights, d(L) / d(W) = [d(L) / d(z)] * in'
            m_dw.noalias() = prev_layer_data * dLz.transpose() / nobs;
            // Derivative for bias, d(L) / d(b) = d(L) / d(z)
            m_db.noalias() = dLz.rowwise().mean();
            // Compute d(L) / d_in = W * [d(L) / d(z)]
            if (m_din.rows() != this->m_in_size && m_din.cols() != nobs) 
				m_din.resize(this->m_in_size, nobs);
            m_din.noalias() = m_weight * dLz;
        }

        const Matrix& backprop_data() const
        {
            return m_din;
        }

        void update(Optimizer& opt)
        {
            AlignedMapVec dw(m_dw.data(), m_dw.size());
            AlignedMapVec db(m_db.data(), m_db.size());
            AlignedMapVec      w(m_weight.data(), m_weight.size());
            AlignedMapVec      b(m_bias.data(), m_bias.size());
            opt.update(dw, w);
            opt.update(db, b);
        }

        std::vector<Scalar> get_parameters() const
        {
            std::vector<Scalar> res(m_weight.size() + m_bias.size());
            // Copy the data of weights and bias to a long vector
            std::copy(m_weight.data(), m_weight.data() + m_weight.size(), res.begin());
            std::copy(m_bias.data(), m_bias.data() + m_bias.size(),
                      res.begin() + m_weight.size());
            return res;
        }

        void set_parameters(const std::vector<Scalar>& param)
        {
            if (static_cast<int>(param.size()) != m_weight.size() + m_bias.size())
            {
                throw std::invalid_argument("[class FullyConnected]: Parameter size does not match");
            }

            std::copy(param.begin(), param.begin() + m_weight.size(), m_weight.data());
            std::copy(param.begin() + m_weight.size(), param.end(), m_bias.data());
        }

        std::vector<Scalar> get_derivatives() const
        {
            std::vector<Scalar> res(m_dw.size() + m_db.size());
            // Copy the data of weights and bias to a long vector
            std::copy(m_dw.data(), m_dw.data() + m_dw.size(), res.begin());
            std::copy(m_db.data(), m_db.data() + m_db.size(), res.begin() + m_dw.size());
            return res;
        }

        std::string layer_type() const
        {
            return "FullyConnected";
        }

        std::string activation_type() const
        {
            return Activation::return_type();
        }

        void fill_meta_info(MetaInfo& map, int index) const
        {
            std::string ind = internal::to_string(index);
            map.insert(std::make_pair("Layer" + ind, internal::layer_id(layer_type())));
            map.insert(std::make_pair("Activation" + ind, internal::activation_id(activation_type())));
            map.insert(std::make_pair("in_size" + ind, in_size()));
            map.insert(std::make_pair("out_size" + ind, out_size()));
        }
};


} // namespace MiniDNN


#endif /* LAYER_FULLYCONNECTED_H_ */

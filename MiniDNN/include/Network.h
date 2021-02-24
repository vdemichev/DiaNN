#ifndef NETWORK_H_
#define NETWORK_H_

#include <Eigen/Core>
#include <vector>
#include <map>
#include <stdexcept>
#include "Config.h"
#include "RNG.h"
#include "Layer.h"
#include "Output.h"
#include "Callback.h"
#include "Utils/Random.h"
#include "Utils/IO.h"
#include "Utils/Factory.h"
#include "math.h"

namespace MiniDNN
{


///
/// \defgroup Network Neural Network Model
///

///
/// \ingroup Network
///
/// This class represents a neural network model that typically consists of a
/// number of hidden layers and an output layer. It provides functions for
/// network building, model fitting, and prediction, etc.
///
class Network
{
    private:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
        typedef Eigen::RowVectorXi IntegerVector;
        typedef std::map<std::string, int> MetaInfo;

        RNG                 m_default_rng;      // Built-in RNG
        RNG&                m_rng;              // Reference to the RNG provided by the user,
                                                // otherwise reference to m_default_rng
        std::vector<Layer*> m_layers;           // Pointers to hidden layers
        Output*             m_output;           // The output layer
        Callback            m_default_callback; // Default callback function
        Callback*           m_callback;         // Points to user-provided callback function,
                                                // otherwise points to m_default_callback
		int dim;

        // Check dimensions of layers
        void check_unit_sizes() const
        {
            const int nlayer = num_layers();

            if (nlayer <= 1)
            {
                return;
            }

            for (int i = 1; i < nlayer; i++)
            {
                if (m_layers[i]->in_size() != m_layers[i - 1]->out_size())
                {
                    throw std::invalid_argument("[class Network]: Unit sizes do not match");
                }
            }
        }

        // Let each layer compute its output
        void forward(const Matrix& input)
        {
            const int nlayer = num_layers();

            if (nlayer <= 0)
            {
                return;
            }

            // First layer
            if (input.rows() != m_layers[0]->in_size())
            {
                throw std::invalid_argument("[class Network]: Input data have incorrect dimension");
            }

            m_layers[0]->forward(input);

            // The following layers
            for (int i = 1; i < nlayer; i++)
            {
                m_layers[i]->forward(m_layers[i - 1]->output());
            }
        }

        // Let each layer compute its gradients of the parameters
        // target has two versions: Matrix and RowVectorXi
        // The RowVectorXi version is used in classification problems where each
        // element is a class label
        template <typename TargetType>
        void backprop(const Matrix& input, const TargetType& target)
        {
            const int nlayer = num_layers();

            if (nlayer <= 0)
            {
                return;
            }

            Layer* first_layer = m_layers[0];
            Layer* last_layer = m_layers[nlayer - 1];
            // Let output layer compute back-propagation data
            m_output->check_target_data(target);
            m_output->evaluate(last_layer->output(), target);

            // If there is only one hidden layer, "prev_layer_data" will be the input data
            if (nlayer == 1)
            {
                first_layer->backprop(input, m_output->backprop_data());
                return;
            }

            // Compute gradients for the last hidden layer
            last_layer->backprop(m_layers[nlayer - 2]->output(), m_output->backprop_data());

            // Compute gradients for all the hidden layers except for the first one and the last one
            for (int i = nlayer - 2; i > 0; i--)
            {
                m_layers[i]->backprop(m_layers[i - 1]->output(),
                                      m_layers[i + 1]->backprop_data());
            }

            // Compute gradients for the first layer
            first_layer->backprop(input, m_layers[1]->backprop_data());
        }

        // Update parameters
        void update(Optimizer& opt)
        {
            const int nlayer = num_layers();

            if (nlayer <= 0)
            {
                return;
            }

            for (int i = 0; i < nlayer; i++)
            {
                m_layers[i]->update(opt);
            }
        }

        // Get the meta information of the network, used to export the NN model
        MetaInfo get_meta_info() const
        {
            const int nlayer = num_layers();
            MetaInfo map;
            map.insert(std::make_pair("Nlayers", nlayer));

            for (int i = 0; i < nlayer; i++)
            {
                m_layers[i]->fill_meta_info(map, i);
            }

            map.insert(std::make_pair("OutputLayer", internal::output_id(m_output->output_type())));
            return map;
        }

    public:
        ///
        /// Default constructor that creates an empty neural network
        ///
        Network() :
            m_default_rng(1),
            m_rng(m_default_rng),
            m_output(NULL),
            m_default_callback(),
            m_callback(&m_default_callback)
        {}

        ///
        /// Constructor with a user-provided random number generator
        ///
        /// \param rng A user-provided random number generator object that inherits
        ///            from the default RNG class.
        ///
        Network(RNG& rng) :
            m_default_rng(1),
            m_rng(rng),
            m_output(NULL),
            m_default_callback(),
            m_callback(&m_default_callback)
        {}

        ///
        /// Destructor that frees the added hidden layers and output layer
        ///
        ~Network()
        {
            const int nlayer = num_layers();

            for (int i = 0; i < nlayer; i++)
            {
                delete m_layers[i];
            }

            if (m_output)
            {
                delete m_output;
            }
        }

        ///
        /// Add a hidden layer to the neural network
        ///
        /// \param layer A pointer to a Layer object, typically constructed from
        ///              layer classes such as FullyConnected and Convolutional.
        ///              **NOTE**: the pointer will be handled and freed by the
        ///              network object, so do not delete it manually.
        ///
        void add_layer(Layer* layer)
        {
            m_layers.push_back(layer);
        }

        ///
        /// Set the output layer of the neural network
        ///
        /// \param output A pointer to an Output object, typically constructed from
        ///               output layer classes such as RegressionMSE and MultiClassEntropy.
        ///               **NOTE**: the pointer will be handled and freed by the
        ///               network object, so do not delete it manually.
        ///
        void set_output(Output* output)
        {
            if (m_output)
            {
                delete m_output;
            }

            m_output = output;
        }

        ///
        /// Number of hidden layers in the network
        ///
        int num_layers() const
        {
            return m_layers.size();
        }

        ///
        /// Get the list of hidden layers of the network
        ///
        std::vector<const Layer*> get_layers() const
        {
            const int nlayer = num_layers();
            std::vector<const Layer*> layers(nlayer);
            std::copy(m_layers.begin(), m_layers.end(), layers.begin());
            return layers;
        }

        ///
        /// Get the output layer
        ///
        const Output* get_output() const
        {
            return m_output;
        }

        ///
        /// Set the callback function that can be called during model fitting
        ///
        /// \param callback A user-provided callback function object that inherits
        ///                 from the default Callback class.
        ///
        void set_callback(Callback& callback)
        {
            m_callback = &callback;
        }
        ///
        /// Set the default silent callback function
        ///
        void set_default_callback()
        {
            m_callback = &m_default_callback;
        }

        ///
        /// Initialize layer parameters in the network using normal distribution
        ///
        /// \param mu    Mean of the normal distribution.
        /// \param sigma Standard deviation of the normal distribution.
        /// \param seed  Set the random seed of the %RNG if `seed > 0`, otherwise
        ///              use the current random state.
        ///
        void init(int _dim, const Scalar& mu = Scalar(0), const Scalar& sigma = Scalar(1.0),
                  int seed = -1)
        {
			dim = _dim;
            check_unit_sizes();

            if (seed > 0)
            {
                m_rng.seed(seed);
            }

            const int nlayer = num_layers();

            for (int i = 0; i < nlayer; i++)
            {
                m_layers[i]->init(mu, sigma / sqrt(0.5 * (double)(m_layers[i]->in_size())), m_rng);
            }
        }

        ///
        /// Get the serialized layer parameters
        ///
        std::vector< std::vector<Scalar> > get_parameters() const
        {
            const int nlayer = num_layers();
            std::vector< std::vector<Scalar> > res;
            res.reserve(nlayer);

            for (int i = 0; i < nlayer; i++)
            {
                res.push_back(m_layers[i]->get_parameters());
            }

            return res;
        }

        ///
        /// Set the layer parameters
        ///
        /// \param param Serialized layer parameters
        ///
        void set_parameters(const std::vector< std::vector<Scalar> >& param)
        {
            const int nlayer = num_layers();

            if (static_cast<int>(param.size()) != nlayer)
            {
                throw std::invalid_argument("[class Network]: Parameter size does not match");
            }

            for (int i = 0; i < nlayer; i++)
            {
                m_layers[i]->set_parameters(param[i]);
            }
        }

        ///
        /// Get the serialized derivatives of layer parameters
        ///
        std::vector< std::vector<Scalar> > get_derivatives() const
        {
            const int nlayer = num_layers();
            std::vector< std::vector<Scalar> > res;
            res.reserve(nlayer);

            for (int i = 0; i < nlayer; i++)
            {
                res.push_back(m_layers[i]->get_derivatives());
            }

            return res;
        }

        ///
        /// Debugging tool to check parameter gradients
        ///
        template <typename TargetType>
        void check_gradient(const Matrix& input, const TargetType& target, int npoints,
                            int seed = -1)
        {
            if (seed > 0)
            {
                m_rng.seed(seed);
            }

            this->forward(input);
            this->backprop(input, target);
            std::vector< std::vector<Scalar> > param = this->get_parameters();
            std::vector< std::vector<Scalar> > deriv = this->get_derivatives();
            const Scalar eps = 1e-5;
            const int nlayer = deriv.size();

            for (int i = 0; i < npoints; i++)
            {
                // Randomly select a layer
                const int layer_id = int(m_rng.rand() * nlayer);
                // Randomly pick a parameter, note that some layers may have no parameters
                const int nparam = deriv[layer_id].size();

                if (nparam < 1)
                {
                    continue;
                }

                const int param_id = int(m_rng.rand() * nparam);
                // Turbulate the parameter a little bit
                const Scalar old = param[layer_id][param_id];
                param[layer_id][param_id] -= eps;
                this->set_parameters(param);
                this->forward(input);
                this->backprop(input, target);
                const Scalar loss_pre = m_output->loss();
                param[layer_id][param_id] += eps * 2;
                this->set_parameters(param);
                this->forward(input);
                this->backprop(input, target);
                const Scalar loss_post = m_output->loss();
                const Scalar deriv_est = (loss_post - loss_pre) / eps / 2;
                std::cout << "[layer " << layer_id << ", param " << param_id <<
                          "] deriv = " << deriv[layer_id][param_id] << ", est = " << deriv_est <<
                          ", diff = " << deriv_est - deriv[layer_id][param_id] << std::endl;
                param[layer_id][param_id] = old;
            }

            // Restore original parameters
            this->set_parameters(param);
        }

        ///
        /// Fit the model based on the given data
        ///
        /// \param opt        An object that inherits from the Optimizer class, indicating the optimization algorithm to use.
        /// \param x          The predictors. Each column is an observation.
        /// \param y          The response variable. Each column is an observation.
        /// \param batch_size Mini-batch size.
        /// \param epoch      Number of epochs of training.
        /// \param seed       Set the random seed of the %RNG if `seed > 0`, otherwise
        ///                   use the current random state.
        ///
        template <typename DerivedX, typename DerivedY>
        bool fit(Optimizer& opt, std::vector<float*> xx, std::vector<float*> yy,
				const Eigen::MatrixBase<DerivedX>& x, const Eigen::MatrixBase<DerivedY>& y, // dummys
				float *m, float *sd,
                int batch_size, int epoch, int seed = -1)
        {
            // We do not directly use PlainObjectX since it may be row-majored if x is passed as mat.transpose()
            // We want to force XType and YType to be column-majored
            typedef typename Eigen::MatrixBase<DerivedX>::PlainObject PlainObjectX;
            typedef typename Eigen::MatrixBase<DerivedY>::PlainObject PlainObjectY;
            typedef Eigen::Matrix<typename PlainObjectX::Scalar, PlainObjectX::RowsAtCompileTime, PlainObjectX::ColsAtCompileTime>
            XType;
            typedef Eigen::Matrix<typename PlainObjectY::Scalar, PlainObjectY::RowsAtCompileTime, PlainObjectY::ColsAtCompileTime>
            YType;
            const int nlayer = num_layers();

            if (nlayer <= 0)
            {
                return false;
            }

            // Reset optimizer
            opt.reset();

            // Create shuffled mini-batches
            if (seed > 0)
            {
                m_rng.seed(seed);
            }

            XType x_batch;
            YType y_batch;

			int N = xx.size();
			assert(N == yy.size());

			std::vector<int> ind(N);
            for (int k = 0; k < epoch; k++)
            {
				for (int i = 0; i < N; i++) ind[i] = i;
				std::mt19937_64 gen(seed + k * 12345);
				std::shuffle(ind.begin(), ind.end(), gen);

				int start = 0;
				for (int i = 0; i < N; i++) {
					int s = i % batch_size;
					if (s == batch_size - 1 || i == N - 1) {
						x_batch.resize(dim, s + 1);
						y_batch.resize(2, s + 1);
						for (int sample = start; sample <= i; sample++) {
							int pos = ind[sample];
							int shift = sample - start;
							for (int f = 0; f < dim; f++) x_batch(f, shift) = (xx[pos][f] - m[f]) / sd[f];
							y_batch(0, shift) = yy[pos][0], y_batch(1, shift) = yy[pos][1];
						}

						this->forward(x_batch);
						this->backprop(x_batch, y_batch);
						this->update(opt);

						start = i + 1;
					}
				}
            }

            return true;
        }

        ///
        /// Use the fitted model to make predictions
        ///
        /// \param x The predictors. Each column is an observation.
        ///
        Matrix predict(const Matrix& x)
        {
            const int nlayer = num_layers();

            if (nlayer <= 0)
            {
                return Matrix();
            }

            this->forward(x);
            return m_layers[nlayer - 1]->output();
        }

        ///
        /// Export the network to files.
        ///
        /// \param folder   The folder where the network is saved.
        /// \param fileName The filename for the network.
        ///
        void export_net(const std::string& folder, const std::string& filename) const
        {
            bool created = internal::create_directory(folder);
            if (!created)
                throw std::runtime_error("[class Network]: Folder creation failed");

            MetaInfo map = this->get_meta_info();
            internal::write_map(folder + "/" + filename, map);
            std::vector< std::vector<Scalar> > params = this->get_parameters();
            internal::write_parameters(folder, filename, params);
        }

        ///
        /// Read in a network from files.
        ///
        /// \param folder   The folder where the network is saved.
        /// \param fileName The filename for the network.
        ///
        void read_net(const std::string& folder, const std::string& filename)
        {
            MetaInfo map;
            internal::read_map(folder + "/" + filename, map);
            int nlayer = map.find("Nlayers")->second;
            std::vector< std::vector<Scalar> > params = internal::read_parameters(folder, filename, nlayer);
            m_layers.clear();

            for (int i = 0; i < nlayer; i++)
            {
                this->add_layer(internal::create_layer(map, i));
            }

            this->set_parameters(params);
            this->set_output(internal::create_output(map));
        }

		Network& operator=(const Network& other) {
			auto map = other.get_meta_info();
			auto par = other.get_parameters();

			int nlayer = par.size();
			m_layers.clear();
			for (int i = 0; i < nlayer; i++) this->add_layer(internal::create_layer(map, i));

			this->set_parameters(par);
			this->set_output(internal::create_output(map));

			return *this;
		}
};


} // namespace MiniDNN


#endif /* NETWORK_H_ */

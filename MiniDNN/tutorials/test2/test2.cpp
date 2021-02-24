#include <MiniDNN.h>
using namespace MiniDNN;

typedef Eigen::MatrixXd Matrix;
typedef Eigen::VectorXd Vector;


int main()
{
    // Set random seed and generate some data
    std::srand(123);
    // Predictors -- each column is an observation
    Matrix x = Matrix::Random(400, 100);
    // Response variables -- each column is an observation
    Matrix y = Matrix::Random(2, 100);
    // Construct a network object
    Network net;
    // Create three layers
    // Layer 1 -- convolutional, input size 20x20x1, 3 output channels, filter size 5x5
    Layer* layer1 = new Convolutional<ReLU>(20, 20, 1, 3, 5, 5);
    // Layer 2 -- max pooling, input size 16x16x3, pooling window size 3x3
    Layer* layer2 = new MaxPooling<ReLU>(16, 16, 3, 3, 3);
    // Layer 3 -- fully connected, input size 5x5x3, output size 2
    Layer* layer3 = new FullyConnected<Identity>(5 * 5 * 3, 2);
    // Add layers to the network object
    net.add_layer(layer1);
    net.add_layer(layer2);
    net.add_layer(layer3);
    // Set output layer
    net.set_output(new RegressionMSE());
    // Create optimizer object
    RMSProp opt;
    opt.m_lrate = 0.001;
    // (Optional) set callback function object
    VerboseCallback callback;
    net.set_callback(callback);
    // Initialize parameters with N(0, 0.01^2) using random seed 123
    net.init(0, 0.01, 123);
    // Fit the model with a batch size of 100, running 10 epochs with random seed 123
    net.fit(opt, x, y, 100, 10, 123);
    net.export_net("NetFolder", "NetFile");
    // Create a new network
    Network netFromFile;
    // Read structure and paramaters from file
    netFromFile.read_net("./NetFolder/", "NetFile");
    // Obtain prediction -- each column is an observation
    std::cout << net.predict(x) << std::endl;
    std::cout << netFromFile.predict(x) - net.predict(x) << std::endl;
    // Layer objects will be freed by the network object,
    // so do not manually delete them
    return 0;
}

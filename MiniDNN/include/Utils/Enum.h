#ifndef UTILS_ENUM_H_
#define UTILS_ENUM_H_

#include <string>
#include <stdexcept>
#include "../Config.h"

namespace MiniDNN
{

namespace internal
{


// Enumerations for hidden layers
enum LAYER_ENUM
{
    FULLY_CONNECTED = 0,
    CONVOLUTIONAL,
    MAX_POOLING
};

// Convert a hidden layer type string to an integer
inline int layer_id(const std::string& type)
{
    if (type == "FullyConnected")
        return FULLY_CONNECTED;
    if (type == "Convolutional")
        return CONVOLUTIONAL;
    if (type == "MaxPooling")
        return MAX_POOLING;

    throw std::invalid_argument("[function layer_id]: Layer is not of a known type");
    return -1;
}

// Enumerations for activation functions
enum ACTIVATION_ENUM
{
    IDENTITY = 0,
    RELU,
    SIGMOID,
    SOFTMAX,
    TANH,
    MISH
};

// Convert an activation type string to an integer
inline int activation_id(const std::string& type)
{
    if (type == "Identity")
        return IDENTITY;
    if (type == "ReLU")
        return RELU;
    if (type == "Sigmoid")
        return SIGMOID;
    if (type == "Softmax")
        return SOFTMAX;
    if (type == "Tanh")
        return TANH;
    if (type == "Mish")
        return MISH;

    throw std::invalid_argument("[function activation_id]: Activation is not of a known type");
    return -1;
}

// Enumerations for output layers
enum OUTPUT_ENUM
{
    REGRESSION_MSE = 0,
    BINARY_CLASS_ENTROPY,
    MULTI_CLASS_ENTROPY
};

// Convert an output layer type string to an integer
inline int output_id(const std::string& type)
{
    if (type == "RegressionMSE")
        return REGRESSION_MSE;
    if (type == "MultiClassEntropy")
        return BINARY_CLASS_ENTROPY;
    if (type == "BinaryClassEntropy")
        return MULTI_CLASS_ENTROPY;

    throw std::invalid_argument("[function output_id]: Output is not of a known type");
    return -1;
}


} // namespace internal

} // namespace MiniDNN


#endif /* UTILS_ENUM_H_ */

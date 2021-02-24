/// \file
/// Source code file of the EigenStream class, it contains the implementation of
/// several methods for input output operations.

#ifndef MiniDNNStream_H
#define MiniDNNStream_H

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <algorithm>
#include <fstream>
#include <iterator>
#include "Assert.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor>
#pragma GCC diagnostic pop
#define MAXBUFSIZE (static_cast<int> (1e6))


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/// Class for input-output manipulation
namespace MiniDNN
{

//--------------------------------------------------------------------------
/// @brief      function to convert a number to a string in plan c++98
///
/// @param[in]  num         The number
///
/// @tparam     numberType  type of number
///
/// @return     a std::string containing the number
///
template<class numberType>
std::string to_string(numberType num)
{
    std::ostringstream convert;
    convert << num;
    return convert.str();
}


//--------------------------------------------------------------------------
/// Export the matrices in numpy (type=python), matlab (type=matlab) format and txt (type=eigen) format
/* In this case the function is implemented for a second order matrix  */
///
/// @param[in] matrice Eigen::MatrixXd that you want to export.
/// @param[in] Name string to identify the name you want to use to save the file.
/// @param[in] type string to identify format to export the matrix if numpy (type="python"), if matlab (type="matlab") if txt (type="eigen").
/// @param[in] folder string to identify the folder where you want to save the file.
///
void export_matrix(Eigen::MatrixXd& matrix, std::string Name,
                   std::string type = "python", std::string folder = "./Model")
{
    mkdir(folder.c_str(), ACCESSPERMS);
    std::string est;

    if (type == "python")
    {
        est = ".py";
        std::string filename(folder + "/" + Name + "_mat" + est);
        std::ofstream str(filename.c_str());
        str << Name << "=np.array([";

        for (int i = 0; i < matrix.rows(); i++)
        {
            for (int j = 0; j < matrix.cols(); j++)
            {
                if (j == 0)
                {
                    str << "[" << matrix(i, j);
                }
                else
                {
                    str << "," << matrix(i, j);
                }
            }

            if (i != (matrix.rows() - 1))
            {
                str << "]," << std::endl;
            }
        }

        str << "]])" << std::endl;
    }

    if (type == "matlab")
    {
        est = ".m";
        std::string filename(folder + "/" + Name + "_mat" + est);
        std::ofstream str(filename.c_str());
        str << Name << "=[";

        for (int i = 0; i < matrix.rows(); i++)
        {
            for (int j = 0; j < matrix.cols(); j++)
            {
                str << " " << matrix(i, j);
            }

            if (i != (matrix.rows() - 1))
            {
                str << ";" << std::endl;
            }
        }

        str << "];" << std::endl;
    }

    if (type == "eigen")
    {
        std::ofstream ofs;
        std::string filename(folder + "/" + Name + "_mat.txt");
        ofs.open(filename.c_str());
        ofs << matrix << std::endl;
        ofs.close();
    }
}

//--------------------------------------------------------------------------
/// @brief      Saves a dense matrix to a binary format file
///
/// @param[in]  Matrix      The  Eigen dense matrix
/// @param[in]  folder      the folder where you want to save the matrix
/// @param[in]  MatrixName  The matrix name for the output file
///
/// @tparam     MatrixType  type of the matrix, i.e. double, float, ...
///
template <typename MatrixType>
void save_dense_matrix(MatrixType& Matrix, std::string folder,
                       std::string MatrixName)
{
    mkdir(folder.c_str(), ACCESSPERMS);
    std::ofstream out(folder + MatrixName,
                      std::ios::out | std::ios::binary | std::ios::trunc);
    typename MatrixType::Index rows = Matrix.rows(), cols = Matrix.cols();
    out.write(reinterpret_cast<char*> (&rows), sizeof(typename MatrixType::Index));
    out.write(reinterpret_cast<char*> (&cols), sizeof(typename MatrixType::Index));
    out.write(reinterpret_cast<char*> (Matrix.data()),
              rows * cols * sizeof(typename MatrixType::Scalar) );
    out.close();
}

//--------------------------------------------------------------------------
/// @brief      Reads a dense matrix from a binary format file
///
/// @param[in,out]  Matrix      The Eigen dense matrix
/// @param[in]      folder      The folder from where you want to read the matrix
/// @param[in]      MatrixName  The matrix name of the input file
///
/// @tparam     MatrixType  type of the matrix, i.e. double, float, ...
///
template <typename MatrixType>
void read_dense_matrix(MatrixType& Matrix, std::string folder,
                       std::string MatrixName)
{
    std::ifstream in;
    in.open((folder + MatrixName).c_str(), std::ios::in | std::ios::binary);
    std::string message(folder + MatrixName +
                        " file does not exist, Check if the file is existing");
    M_Assert(in.good(), message.c_str());

    if (in.is_open())
    {
        typename MatrixType::Index rows = 0, cols = 0;
        in.read(reinterpret_cast<char*> (&rows), sizeof(typename MatrixType::Index));
        in.read(reinterpret_cast<char*> (&cols), sizeof(typename MatrixType::Index));
        Matrix.resize(rows, cols);
        in.read( reinterpret_cast<char*>(Matrix.data()),
                 rows * cols * sizeof(typename MatrixType::Scalar) );
        in.close();
    }
}

//----------------------------------------------------------------------
/// @brief      Saves a dense tensor.
///
/// @param      Tensor      The tensor
/// @param[in]  folder      The folder
/// @param[in]  MatrixName  The matrix name
///
/// @tparam     TensorType  type of the tensor, i.e. double, float, ...
///
template <typename TensorType>
void save_dense_tensor(TensorType& Tensor, std::string folder,
                       std::string MatrixName)
{
    std::ofstream out(folder + MatrixName,
                      std::ios::out | std::ios::binary | std::ios::trunc);
    typename TensorType::Dimensions dim = Tensor.dimensions();
    int tot = 1;

    for (unsigned int k = 0; k < dim.size(); k++)
    {
        tot *= dim[k];
    }

    out.write(reinterpret_cast<char*> (&dim),
              sizeof(typename TensorType::Dimensions));
    out.write(reinterpret_cast<char*> (Tensor.data()),
              tot * sizeof(typename TensorType::Scalar) );
    out.close();
}

//----------------------------------------------------------------------
/// @brief      Reads a dense tensor.
///
/// @param      Tensor      The tensor
/// @param[in]  folder      The folder
/// @param[in]  MatrixName  The matrix name
///
/// @tparam     TensorType   type of the tensor, i.e. double, float, ...
///
template <typename TensorType>
void read_dense_tensor(TensorType& Tensor, std::string folder,
                       std::string MatrixName)
{
    std::ifstream in;
    in.open((folder + MatrixName).c_str(), std::ios::in | std::ios::binary);
    typename TensorType::Dimensions dim;
    in.read(reinterpret_cast<char*> (&dim),
            sizeof(typename TensorType::Dimensions));
    const typename TensorType::Dimensions& dims = Tensor.dimensions();
    M_Assert(dims.size() == dim.size(),
             "The rank of the tensor you want to fill does not coincide with the rank of the tensor you are reading");
    int tot = 1;

    for (unsigned int k = 0; k < dim.size(); k++)
    {
        tot *= dim[k];
    }

    Tensor.resize(dim);
    in.read( reinterpret_cast<char*>(Tensor.data()),
             tot * sizeof(typename TensorType::Scalar) );
    in.close();
}

///
/// @brief      Write a std::vector<Scalar> to file.
///
/// @param[in]  myVector  The vector you want to write
/// @param[in]  filename  The filename of the std::vector
///
void write_vector_to_file(const std::vector<Scalar>& myVector,
                          std::string filename)
{
    std::ofstream ofs(filename.c_str(), std::ios::out | std::ios::binary);
    std::ostream_iterator<char> osi(ofs);
    const char* beginByte = (char*)&myVector[0];
    const char* endByte = (char*)&myVector.back() + sizeof(Scalar);
    std::copy(beginByte, endByte, osi);
}
///
/// @brief      Reads a std::vector<Scalar> from file.
///
/// @param[in]  filename  The filename of the vector
///
/// @return     The vector
///
std::vector<Scalar> read_vector_from_file(std::string filename)
{
    std::vector<char> buffer;
    std::ifstream ifs(filename.c_str(), std::ios::in | std::ifstream::binary);
    std::istreambuf_iterator<char> iter(ifs);
    std::istreambuf_iterator<char> end;
    std::copy(iter, end, std::back_inserter(buffer));
    std::vector<Scalar> newVector(buffer.size() / sizeof(Scalar));
    memcpy(&newVector[0], &buffer[0], buffer.size());
    return newVector;
}

///
/// @brief      Writes parameters of the net from file.
///
/// @param[in]  folder    The folder where the parameter files are stored
/// @param[in]  fileName  The file name prefix of the parameter files
/// @param      params    The parameters of the net
///
void write_parameters(std::string folder, std::string fileName,
                      std::vector<std::vector< Scalar> >& params)
{
    for (int i = 0; i < params.size(); i++)
    {
        write_vector_to_file(params[i], folder + "/" + fileName + to_string(i));
    }
}

///
/// @brief      Reads parameters of the net from file
///
/// @param[in]  folder    The folder where the parameter files are stored
/// @param[in]  fileName  The file name prefix of the parameter files
/// @param[in]  Nlayers   The number of layers
///
/// @return     a vector which contains the parameters
///
std::vector<std::vector< Scalar> > read_parameters(std::string folder,
        std::string fileName, int Nlayers)
{
    std::vector<std::vector< Scalar> > params;

    for (int i = 0; i < Nlayers; i++)
    {
        params.push_back(read_vector_from_file(folder + "/" + fileName + to_string(i)));
    }

    return params;
}



///
/// @brief      Convert a layer type string to an integer
///
/// @param[in]  type  The type of layer (string)
///
/// @return     0 for Convolutional, 1 for MaxPooling, 2 for FullyConnected
///
int layer_type(std::string type)
{
    M_Assert(type == "Convolutional" || type == "MaxPooling" ||
             type == "FullyConnected", "Layer is not of a known type");
    int out;

    if (type == "Convolutional")
    {
        out = 0;
    }

    if (type == "MaxPooling")
    {
        out = 1;
    }

    if (type == "FullyConnected")
    {
        out = 2;
    }

    return out;
}
///
/// @brief      Convert an activation type string to an integer
///
/// @param[in]  type  The type of the activation (string)
///
/// @return     0 for Identity, 1 for ReLU, 2 for Sigmoid. 3 for Softmax, 4 for Mish
///
int activation_type(std::string type)
{
    M_Assert(type == "Identity" || type == "ReLU" ||
             type == "Sigmoid" || type == "Softmax" ||
             type == "Mish" || type == "Tanh", "Activation is not of a known type");
    int out;

    if (type == "Identity")
    {
        out = 0;
    }

    if (type == "ReLU")
    {
        out = 1;
    }

    if (type == "Sigmoid")
    {
        out = 2;
    }

    if (type == "Softmax")
    {
        out = 3;
    }

    if (type == "Mish")
    {
        out = 4;
    }

    if (type == "Tanh")
    {
        out = 5;
    }

    return out;
}
///
/// @brief      Convert an output layer type string to an integer
///
/// @param[in]  type  The type of the output layer in string
///
/// @return     0 for RegressionMSE, 1 for MultiClassEntropy, 2 for BinaryClassEntropy.
///
int output_type(std::string type)
{
    M_Assert(type == "RegressionMSE" || type == "MultiClassEntropy" ||
             type == "BinaryClassEntropy", "output is not of a known type");
    int out;

    if (type == "RegressionMSE")
    {
        out = 0;
    }

    if (type == "MultiClassEntropy")
    {
        out = 1;
    }

    if (type == "BinaryClassEntropy")
    {
        out = 2;
    }

    return out;
}



//--------------------------------------------------------------------------
/// @brief      Writes a map model.
///
/// @param[in]  fileName  The file name of the model you want to write
/// @param[in]  map       The model you want to write
///
/// @return     integer for success
///
int write_map (std::string fileName, std::map<std::string, int> map)
{
    int count = 0;

    if (map.empty())
    {
        return 0;
    }

    FILE* fp = fopen(fileName.c_str(), "w");

    if (! fp)
    {
        return - errno;
    }

    for (std::map<std::string, int>::iterator it = map.begin();
            it != map.end(); it++)
    {
        fprintf(fp, "%s=%d\n", it->first.c_str(), it->second);
        count++;
    }

    fclose(fp);
    return count;
}

//--------------------------------------------------------------------------
/// @brief      Read a map Model.
///
/// @param[in]  fileName  The file name of the model you want to read.
/// @param[in]  map       Map to store the model.
///
/// @return     integer for success
///
int read_map (std::string fname, std::map<std::string, int>& map)
{
    int count = 0;
    FILE* fp = fopen(fname.c_str(), "r");

    if (!fp)
    {
        return -errno;
    }

    map.clear();
    char* buf = 0;
    size_t buflen = 0;

    while (getline(&buf, &buflen, fp) > 0)
    {
        char* nl = strchr(buf, '\n');

        if (nl == NULL)
        {
            continue;
        }

        *nl = 0;
        char* sep = strchr(buf, '=');

        if (sep == NULL)
        {
            continue;
        }

        *sep = 0;
        sep++;
        std::string s1 = buf;
        int s2 = atoi(sep);
        (map)[s1] = s2;
        count++;
    }

    if (buf)
    {
        free(buf);
    }

    fclose(fp);
    return count;
}
};


std::ostream& operator<<(std::ostream& os,
                         std::map<std::string, int> const& myMap)
{
    for (std::map<std::string, int>::const_iterator it = myMap.begin();
            it != myMap.end();
            ++it)
    {
        os << it->first << " " << it->second << "\n";
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, std::vector<Scalar>& myVector)
{
    for (int i = 0; i < myVector.size(); i++)
    {
        os << myVector[i] << "\n";
    }

    return os;
}



#endif







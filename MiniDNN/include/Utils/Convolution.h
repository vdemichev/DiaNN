#ifndef UTILS_CONVOLUTION_H_
#define UTILS_CONVOLUTION_H_

#include <Eigen/Core>
#include "../Config.h"

namespace MiniDNN
{

namespace internal
{


// We assume the following memory layout:
// There are 'n_obs' images, each with 'in_channels' channels
// Each channel has 'channel_rows' rows and 'channel_cols' columns
// The data starts from 'src'
// If 'image_outer_loop == true', the data first iterates on channels,
// and then images:
/*
 * ###############################################################
 * #           #           #           #           #
 * # channel 1 # channel 2 # channel 3 # channel 1 # ...
 * #           #           #           #           #
 * ###############################################################
 * |<------------ image 1 ------------>|<------------ image 2 ----
 */
// If 'image_outer_loop == false', the layout looks like below:
/*
 * ###############################################################
 * #           #           #           #           #
 * #  image 1  #  image 2  #  image 3  #  image 1  # ...
 * # channel 1 # channel 1 # channel 1 # channel 2 #
 * #           #           #           #           #
 * ###############################################################
 * |<----------- channel 1 ----------->|<----------- channel 2 ----
 */
//
// Then we assume there are 'out_channels' output channels, so in total
// we have 'in_channels * out_channels' filters for each image
// Each filter has 'filter_rows' rows and 'filter_cols' columns
// Filters start from 'filter_data'. The layout looks like below, with each
// block consisting of 'filter_rows * filter_cols' elements
/*
 * #########################################################################
 * #               #               #               #               #
 * # out channel 1 # out channel 2 # out channel 3 # out channel 1 # ...
 * #               #               #               #               #
 * #########################################################################
 * |<---------------- in channel 1 --------------->|<---------------- in channel 2 ----
 */
//
// Convolution results from different input channels are summed up to produce the
// result for each output channel
// Convolution results for different output channels are concatenated to preoduce
// the result for each image
//
// The final result is written to the memory pointed by 'dest', with a similar
// layout to 'src' in the 'image_outer_loop == true' case
//
// Memory efficient convolution (MEC)
// Algorithm is based on https://arxiv.org/abs/1706.06873
//
// First define a simple structure to store the various dimensions of convolution
struct ConvDims
{
    // Input parameters
    const int in_channels;
    const int out_channels;
    const int channel_rows;
    const int channel_cols;
    const int filter_rows;
    const int filter_cols;
    // Image dimension -- one observation with all channels
    const int img_rows;
    const int img_cols;
    // Dimension of the convolution result for each output channel
    const int conv_rows;
    const int conv_cols;

    ConvDims(
        const int in_channels_, const int out_channels_,
        const int channel_rows_, const int channel_cols_,
        const int filter_rows_, const int filter_cols_
    ) :
        in_channels(in_channels_), out_channels(out_channels_),
        channel_rows(channel_rows_), channel_cols(channel_cols_),
        filter_rows(filter_rows_), filter_cols(filter_cols_),
        img_rows(channel_rows_), img_cols(in_channels_ * channel_cols_),
        conv_rows(channel_rows_ - filter_rows_ + 1),
        conv_cols(channel_cols_ - filter_cols_ + 1)
    {}
};
// Transform original matrix to "lower" form as described in the MEC paper
// I feel that it is better called the "flat" form
//
// Helper function to "flatten" source images
// 'flat_mat' will be overwritten
// We focus on one channel, and let 'stride' be the distance between two images
inline void flatten_mat(
    const ConvDims& dim, const Scalar* src, const int stride, const int n_obs,
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& flat_mat
)
{
    // Number of bytes in the segment that will be copied at one time
    const int& segment_size = dim.filter_rows;
    const std::size_t copy_bytes = sizeof(Scalar) * segment_size;
    Scalar* writer = flat_mat.data();
    const int channel_size = dim.channel_rows * dim.channel_cols;

    for (int i = 0; i < n_obs; i++, src += stride)
    {
        const Scalar* reader_row = src;
        const Scalar* const reader_row_end = src + dim.conv_rows;

        for (; reader_row < reader_row_end; reader_row++)
        {
            const Scalar* reader = reader_row;
            const Scalar* const reader_end = reader + channel_size;

            for (; reader < reader_end; reader += dim.channel_rows, writer += segment_size)
            {
                std::memcpy(writer, reader, copy_bytes);
            }
        }
    }
}
// A special matrix product. We select a window from 'mat1' and calculates its product with 'mat2',
// and progressively move the window to the right
inline void moving_product(
    const int step,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>&
    mat1,
    Eigen::Map< const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> >& mat2,
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& res
)
{
    const int row1 = mat1.rows();
    const int col1 = mat1.cols();
    const int row2 = mat2.rows();
    const int col2 = mat2.cols();
    const int col_end = col1 - row2;
    int res_start_col = 0;

    for (int left_end = 0; left_end <= col_end;
            left_end += step, res_start_col += col2)
    {
        res.block(0, res_start_col, row1, col2).noalias() += mat1.block(0, left_end,
                row1, row2) * mat2;
    }
}
// The main convolution function using the "valid" rule
inline void convolve_valid(
    const ConvDims& dim,
    const Scalar* src, const bool image_outer_loop, const int n_obs,
    const Scalar* filter_data,
    Scalar* dest)
{
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    RMatrix;
    typedef Eigen::Map<const Matrix> ConstMapMat;
    // Flat matrix
    const int flat_rows = dim.conv_rows * n_obs;
    const int flat_cols = dim.filter_rows * dim.channel_cols;
    const int channel_size = dim.channel_rows * dim.channel_cols;
    // Distance between two images
    const int img_stride = image_outer_loop ? (dim.img_rows * dim.img_cols) :
                           channel_size;
    // Distance between two channels
    const int channel_stride = image_outer_loop ? channel_size :
                               (channel_size * n_obs);
    RMatrix flat_mat(flat_rows, flat_cols);
    // Convolution results
    const int& res_rows = flat_rows;
    const int res_cols = dim.conv_cols * dim.out_channels;
    Matrix res = Matrix::Zero(res_rows, res_cols);
    const int& step = dim.filter_rows;
    const int filter_size = dim.filter_rows * dim.filter_cols;
    const int filter_stride = filter_size * dim.out_channels;

    for (int i = 0; i < dim.in_channels;
            i++, src += channel_stride, filter_data += filter_stride)
    {
        // Flatten source image
        flatten_mat(dim, src, img_stride, n_obs, flat_mat);
        // Compute the convolution result
        ConstMapMat filter(filter_data, filter_size, dim.out_channels);
        moving_product(step, flat_mat, filter, res);
    }

    // The layout of 'res' is very complicated
    /*
     * obs0_out0[0, 0] obs0_out1[0, 0] obs0_out2[0, 0] obs0_out0[0, 1] obs0_out1[0, 1] obs0_out2[0, 1] ...
     * obs0_out0[1, 0] obs0_out1[1, 0] obs0_out2[1, 0] obs0_out0[1, 1] obs0_out1[1, 1] obs0_out2[1, 1] ...
     * obs0_out0[2, 0] obs0_out1[2, 0] obs0_out2[2, 0] obs0_out0[2, 1] obs0_out1[2, 1] obs0_out2[2, 1] ...
     * obs1_out0[0, 0] obs1_out1[0, 0] obs1_out2[0, 0] obs1_out0[0, 1] obs1_out1[0, 1] obs1_out2[0, 1] ...
     * obs1_out0[1, 0] obs1_out1[1, 0] obs1_out2[1, 0] obs1_out0[1, 1] obs1_out1[1, 1] obs1_out2[1, 1] ...
     * obs1_out0[2, 0] obs1_out1[2, 0] obs1_out2[2, 0] obs1_out0[2, 1] obs1_out1[2, 1] obs1_out2[2, 1] ...
     * ...
     *
     */
    // obs<k>_out<l> means the convolution result of the k-th image on the l-th output channel
    // [i, j] gives the matrix indices
    // The destination has the layout
    /*
     * obs0_out0[0, 0] obs0_out0[0, 1] obs0_out0[0, 2] obs0_out1[0, 0] obs0_out1[0, 1] obs0_out1[0, 2] ...
     * obs0_out0[1, 0] obs0_out0[1, 1] obs0_out0[1, 2] obs0_out1[1, 0] obs0_out1[1, 1] obs0_out1[1, 2] ...
     * obs0_out0[2, 0] obs0_out0[2, 1] obs0_out0[2, 2] obs0_out1[2, 0] obs0_out1[2, 1] obs0_out1[2, 2] ...
     *
     */
    // which in a larger scale looks like
    // [obs0_out0 obs0_out1 obs0_out2 obs1_out0 obs1_out1 obs1_out2 obs2_out0 ...]
    // Copy data to destination
    // dest[a, b] corresponds to obs<k>_out<l>[i, j]
    // where k = b / (conv_cols * out_channels),
    //       l = (b % (conv_cols * out_channels)) / conv_cols
    //       i = a,
    //       j = b % conv_cols
    // and then obs<k>_out<l>[i, j] corresponds to res[c, d]
    // where c = k * conv_rows + i,
    //       d = j * out_channels + l
    const int dest_rows = dim.conv_rows;
    const int dest_cols = res_cols * n_obs;
    const Scalar* res_data = res.data();
    const std::size_t copy_bytes = sizeof(Scalar) * dest_rows;

    for (int b = 0; b < dest_cols; b++, dest += dest_rows)
    {
        const int k = b / res_cols;
        const int l = (b % res_cols) / dim.conv_cols;
        const int j = b % dim.conv_cols;
        const int d = j * dim.out_channels + l;
        const int res_col_head = d * res_rows;
        std::memcpy(dest, res_data + res_col_head + k * dim.conv_rows, copy_bytes);
    }
}



// The moving_product() function for the "full" rule
inline void moving_product(
    const int padding, const int step,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>&
    mat1,
    const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& mat2,
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& res
)
{
    const int row1 = mat1.rows();
    const int col1 = mat1.cols();
    const int row2 = mat2.rows();
    const int col2 = mat2.cols();
    int res_start_col = 0;
    // Left padding
    int left_end = -padding;
    int right_end = step;

    for (; left_end < 0
            && right_end <= col1;
            left_end += step, right_end += step, res_start_col += col2)
    {
        res.block(0, res_start_col, row1, col2).noalias() += mat1.leftCols(right_end) *
                mat2.bottomRows(right_end);
    }

    // Main part
    for (; right_end <= col1;
            left_end += step, right_end += step, res_start_col += col2)
    {
        res.block(0, res_start_col, row1, col2).noalias() += mat1.block(0, left_end,
                row1, row2) * mat2;
    }

    // Right padding
    for (; left_end < col1; left_end += step, res_start_col += col2)
    {
        if (left_end <= 0)
        {
            res.block(0, res_start_col, row1, col2).noalias() += mat1 * mat2.block(0,
                    -left_end, col1, row2);
        }
        else
        {
            const int overlap = col1 - left_end;
            res.block(0, res_start_col, row1, col2).noalias() += mat1.rightCols(overlap) *
                    mat2.topRows(overlap);
        }
    }
}
// The main convolution function for the "full" rule
inline void convolve_full(
    const ConvDims& dim,
    const Scalar* src, const int n_obs, const Scalar* filter_data,
    Scalar* dest)
{
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    RMatrix;
    typedef Eigen::Map<const Matrix> ConstMapMat;
    // Padding sizes
    const int padding_top = dim.filter_rows - 1;
    const int padding_left = dim.filter_cols - 1;
    // Dimension of convolution result using "full" rule
    const int conv_rows = dim.channel_rows + padding_top;
    const int conv_cols = dim.channel_cols + padding_left;
    // Add (top and bottom) padding to source images
    const int pad_rows = dim.img_rows + padding_top * 2;
    const int pad_cols = dim.img_cols * n_obs;
    Matrix pad_mat(pad_rows, pad_cols);
    ConstMapMat src_mat(src, dim.img_rows, pad_cols);
    pad_mat.topRows(padding_top).setZero();
    pad_mat.bottomRows(padding_top).setZero();
    pad_mat.block(padding_top, 0, dim.img_rows, pad_cols).noalias() = src_mat;
    src = pad_mat.data();
    ConvDims pad_dim(dim.in_channels, dim.out_channels, pad_rows, dim.channel_cols,
                     dim.filter_rows, dim.filter_cols);
    // Flat matrix
    const int flat_rows = conv_rows * n_obs;
    const int flat_cols = dim.filter_rows * dim.channel_cols;
    const int img_stride = pad_rows * dim.img_cols;
    const int channel_stride = pad_rows * dim.channel_cols;
    RMatrix flat_mat(flat_rows, flat_cols);
    // The processing of filters are different from the "valid" rule in two ways:
    // 1. The layout of input channels and output channels are switched
    // 2. The filters need to be rotated, which is equivalent to reversing the vector of each filter
    // We also separate filters that belong to different input channels
    std::vector<Matrix> filters_in(dim.in_channels);
    const int filter_size = dim.filter_rows * dim.filter_cols;
    const int nfilter = dim.in_channels * dim.out_channels;

    for (int i = 0; i < dim.in_channels; i++)
    {
        filters_in[i].resize(filter_size, dim.out_channels);
    }

    const Scalar* reader = filter_data;

    for (int i = 0; i < nfilter; i++, reader += filter_size)
    {
        Scalar* writer = filters_in[i % dim.in_channels].data() +
                         (i / dim.in_channels) * filter_size;
        std::reverse_copy(reader, reader + filter_size, writer);
    }

    // Convolution results
    const int& res_rows = flat_rows;
    const int res_cols = conv_cols * dim.out_channels;
    Matrix res = Matrix::Zero(res_rows, res_cols);
    const int& step = dim.filter_rows;
    const int filter_padding = padding_left * dim.filter_rows;

    for (int i = 0; i < dim.in_channels; i++, src += channel_stride)
    {
        // Flatten source image
        flatten_mat(pad_dim, src, img_stride, n_obs, flat_mat);
        // Compute the convolution result
        moving_product(filter_padding, step, flat_mat, filters_in[i], res);
    }

    // Copy results to destination
    const int& dest_rows = conv_rows;
    const int  dest_cols = res_cols * n_obs;
    const Scalar* res_data = res.data();
    const std::size_t copy_bytes = sizeof(Scalar) * dest_rows;

    for (int b = 0; b < dest_cols; b++, dest += dest_rows)
    {
        const int k = b / res_cols;
        const int l = (b % res_cols) / conv_cols;
        const int j = b % conv_cols;
        const int d = j * dim.out_channels + l;
        const int res_col_head = d * res_rows;
        std::memcpy(dest, res_data + res_col_head + k * conv_rows, copy_bytes);
    }
}


} // namespace internal

} // namespace MiniDNN


#endif /* UTILS_CONVOLUTION_H_ */

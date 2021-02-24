#ifndef UTILS_FINDMAX_H_
#define UTILS_FINDMAX_H_

#include "../Config.h"

namespace MiniDNN
{

namespace internal
{


// Find the location of the maximum element in x[0], x[1], ..., x[n-1]
// Special cases for small n using recursive template
// N is assumed to be >= 2
template <int N>
inline int find_max(const Scalar* x)
{
    const int loc = find_max < N - 1 > (x);
    return (x[N - 1] > x[loc]) ? (N - 1) : loc;
}

template <>
inline int find_max<2>(const Scalar* x)
{
    return int(x[1] > x[0]);
}

// n is assumed be >= 2
inline int find_max(const Scalar* x, const int n)
{
    switch (n)
    {
        case 2:
            return find_max<2>(x);

        case 3:
            return find_max<3>(x);

        case 4:
            return find_max<4>(x);

        case 5:
            return find_max<5>(x);
    }

    int loc = find_max<6>(x);

    for (int i = 6; i < n; i++)
    {
        loc = (x[i] > x[loc]) ? i : loc;
    }

    return loc;
}

// Find the maximum element in the block x[0:(nrow-1), 0:(ncol-1)]
// col_stride is the distance between x[0, 0] and x[0, 1]
// Special cases for small n
inline Scalar find_block_max(const Scalar* x, const int nrow, const int ncol,
                             const int col_stride, int& loc)
{
    // Max element in the first column
    loc = find_max(x, nrow);
    Scalar val = x[loc];
    // 2nd column
    x += col_stride;
    int loc_next = find_max(x, nrow);
    Scalar val_next = x[loc_next];

    if (val_next > val)
    {
        loc = col_stride + loc_next;
        val = val_next;
    }

    if (ncol == 2)
    {
        return val;
    }

    // 3rd column
    x += col_stride;
    loc_next = find_max(x, nrow);
    val_next = x[loc_next];

    if (val_next > val)
    {
        loc = 2 * col_stride + loc_next;
        val = val_next;
    }

    if (ncol == 3)
    {
        return val;
    }

    // 4th column
    x += col_stride;
    loc_next = find_max(x, nrow);
    val_next = x[loc_next];

    if (val_next > val)
    {
        loc = 3 * col_stride + loc_next;
        val = val_next;
    }

    if (ncol == 4)
    {
        return val;
    }

    // 5th column
    x += col_stride;
    loc_next = find_max(x, nrow);
    val_next = x[loc_next];

    if (val_next > val)
    {
        loc = 4 * col_stride + loc_next;
        val = val_next;
    }

    if (ncol == 5)
    {
        return val;
    }

    // Other columns
    for (int i = 5; i < ncol; i++)
    {
        x += col_stride;
        loc_next = find_max(x, nrow);
        val_next = x[loc_next];

        if (val_next > val)
        {
            loc = i * col_stride + loc_next;
            val = val_next;
        }
    }

    return val;
}


} // namespace internal

} // namespace MiniDNN


#endif /* UTILS_FINDMAX_H_ */

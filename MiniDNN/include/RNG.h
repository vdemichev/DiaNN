#ifndef RNG_H_
#define RNG_H_

namespace MiniDNN
{


///
/// The interface and default implementation of the random number generator (%RNG).
/// The %RNG is used to shuffle data and to initialize parameters of hidden layers.
/// This default implementation is based on the public domain code by Ray Gardner
/// <http://stjarnhimlen.se/snippets/rg_rand.c>.
///
class RNG
{
    private:
        const unsigned int m_a;     // multiplier
        const unsigned long m_max;  // 2^31 - 1
        long m_rand;

        inline long next_long_rand(long seed)
        {
            unsigned long lo, hi;
            lo = m_a * (long)(seed & 0xFFFF);
            hi = m_a * (long)((unsigned long)seed >> 16);
            lo += (hi & 0x7FFF) << 16;

            if (lo > m_max)
            {
                lo &= m_max;
                ++lo;
            }

            lo += hi >> 15;

            if (lo > m_max)
            {
                lo &= m_max;
                ++lo;
            }

            return (long)lo;
        }
    public:
        RNG(unsigned long init_seed) :
            m_a(16807),
            m_max(2147483647L),
            m_rand(init_seed ? (init_seed & m_max) : 1)
        {}

        virtual ~RNG() {}

        virtual void seed(unsigned long seed)
        {
            m_rand = (seed ? (seed & m_max) : 1);
        }

        virtual Scalar rand()
        {
            m_rand = next_long_rand(m_rand);
            return Scalar(m_rand) / Scalar(m_max);
        }
};


} // namespace MiniDNN


#endif /* RNG_H_ */

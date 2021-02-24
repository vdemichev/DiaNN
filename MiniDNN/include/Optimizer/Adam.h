#ifndef OPTIMIZER_ADAM_H_
#define OPTIMIZER_ADAM_H_

#include <Eigen/Core>
#include <map>
#include "../Config.h"
#include "../Optimizer.h"

#define _USE_MATH_DEFINES
#include "math.h"

#define Min(x, y) ((x) < (y) ? (x) : (y))
#define Max(x, y) ((x) > (y) ? (x) : (y))
#define Abs(x) ((x) >= 0 ? (x) : (-(x)))
#define Sgn(x) ((x) >= 0.0 ? (1) : (-1))
#define Sqr(x) ((x) * (x))
#define Cube(x) ((x) * Sqr(x))

namespace MiniDNN
{


///
/// \ingroup Optimizers
///
/// The Adam algorithm
///
class Adam: public Optimizer
{
    private:
        typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
        typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> Array;
        typedef Vector::ConstAlignedMapType ConstAlignedMapVec;
        typedef Vector::AlignedMapType AlignedMapVec;

        std::map<const Scalar*, Array> m_history_m;
        std::map<const Scalar*, Array> m_history_v;
        Scalar m_beta1t;
        Scalar m_beta2t;
		long long m_iter;

    public:
        Scalar m_lrate;
        Scalar m_eps;
        Scalar m_beta1;
        Scalar m_beta2;
		Scalar m_gradient_norm;
		Scalar m_gradient_clip;
		Scalar m_update_clip;
		Scalar m_reg_radius;
		Scalar m_reg_rate;
		bool m_cosine_annealing;

        Adam(const Scalar& lrate = Scalar(0.001), 
			 const Scalar& gradient_norm = Scalar(1.0), const Scalar& gradient_clip = Scalar(10.0), const Scalar& update_clip = Scalar(0.1),
			 const Scalar& reg_radius = Scalar(100.0), 
			 const Scalar& reg_rate = Scalar(0.0), // 0 = use learning rate
			 bool cosine_annealing = true,  
			 const Scalar& eps = Scalar(1e-6),
             const Scalar& beta1 = Scalar(0.9), const Scalar& beta2 = Scalar(0.999)) :
            m_beta1t(beta1), m_beta2t(beta2),
            m_lrate(lrate), m_eps(eps),
            m_beta1(beta1), m_beta2(beta2),
			m_gradient_norm(gradient_norm), m_gradient_clip(gradient_clip), m_update_clip(update_clip),
			m_reg_radius(reg_radius), m_reg_rate(reg_rate),
			m_cosine_annealing(cosine_annealing)
        {}

        void reset()
        {
            m_history_m.clear();
            m_history_v.clear();
            m_beta1t = m_beta1;
            m_beta2t = m_beta2;
			m_iter = 1;
        }

        // https://ruder.io/optimizing-gradient-descent/index.html
        void update(AlignedMapVec& dvec, AlignedMapVec& vec)
        {
            using std::sqrt;
            // Get the m and v vectors associated with this gradient
            Array& mvec = m_history_m[dvec.data()];
            Array& vvec = m_history_v[dvec.data()];

            // If length is zero, initialize it
            if (mvec.size() == 0)
            {
                mvec.resize(dvec.size());
                mvec.setZero();
            }

            if (vvec.size() == 0)
            {
                vvec.resize(dvec.size());
                vvec.setZero();
            }

			double norm = sqrt(dvec.array().square().mean());
			if (norm > m_gradient_norm) dvec.array() *= m_gradient_norm / norm;
			dvec.array() = dvec.array().unaryExpr([](double v) { return std::isfinite(v) ? v : 0.0; });
			dvec.array() = dvec.array().cwiseMax(-m_gradient_clip).cwiseMin(m_gradient_clip);

            mvec = m_beta1 * mvec + (Scalar(1) - m_beta1) * dvec.array();
            vvec = m_beta2 * vvec + (Scalar(1) - m_beta2) * dvec.array().square();
            // Correction coefficients
            const Scalar correct1 = Scalar(1) / (Scalar(1) - m_beta1t);
            const Scalar correct2 = Scalar(1) / sqrt(Scalar(1) - m_beta2t);

			// Cosine annealing
			double ca = 1.0;
			if (m_cosine_annealing) if (m_iter > 16) {
				double bin = (m_iter >> 1) << 1;
				double current = m_iter - bin;
				double phase = current / bin;
				ca = 0.01 + 0.99 * 0.5 * (1.0 + cos(phase * M_PI));
			}

            // Update parameters
            vec.array() -= ((ca * m_lrate * correct1) * mvec / (correct2 * vvec.sqrt() + m_eps)).cwiseMax(-m_update_clip).cwiseMin(m_update_clip);

			// Weight decay
			double rr = m_reg_rate > 0.00000001 ? m_reg_rate : m_lrate;
			rr /= Sqr(m_reg_radius);
			rr *= ca;

			vec.array() = vec.array().unaryExpr([&](double v) { 
				double a = Abs(v);
				double x = a * a;
				double reg = rr * x;
				return v > 0.0 ? Max(v * 0.9, v - reg) : Min(v * 0.9, v + reg); 
			});

            m_beta1t *= m_beta1;
            m_beta2t *= m_beta2;
			m_iter++;
        }
};


} // namespace MiniDNN


#endif /* OPTIMIZER_ADAM_H_ */

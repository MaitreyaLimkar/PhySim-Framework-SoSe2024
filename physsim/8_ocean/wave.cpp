#include "wave.hpp"

namespace physsim
{
    double Wave::waveNumber() const { return 2 * (float)EIGEN_PI / waveLength; }

    Eigen::Vector2d Wave::direction() const { return Eigen::Vector2d(std::cos(angle), std::sin(angle)); }

    Eigen::Vector2d Wave::waveVector() const { return direction() * waveNumber(); }

    double Wave::angularFrequency() const
    {
        double k = waveNumber();
        return k * dispersion.phaseSpeed(k);
    }

    Eigen::Vector3d Wave::offset(const Eigen::Vector2d& position, const double& t) const
    {
        // TODO: compute the offset for the displaced wave from the given position.
        float eta = amplitude * cos((waveVector()).dot(position) - angularFrequency() * t);
        Eigen::Vector2d h_disp = - steepness * (waveVector() / waveNumber()) * amplitude * sin((waveVector()).dot(position) - angularFrequency() * t);
        return Eigen::Vector3d(h_disp.x(), h_disp.y(), eta);
    }
}

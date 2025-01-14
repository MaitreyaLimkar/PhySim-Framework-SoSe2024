#include "fluid_solver.hpp"

#include <vislab/field/regular_field.hpp>
#include <vislab/graphics/colormap_texture.hpp>

namespace physsim
{
    FluidSolver::FluidSolver(const Eigen::AlignedBox2d& domain, const Eigen::Vector2i& resolution)
        : mDomain(domain)
        , mResolution(resolution)
        , mDensity(std::make_shared<vislab::RegularSteadyScalarField2f>())
        , mPressure(std::make_shared<vislab::RegularSteadyScalarField2f>())
        , mDivergence(std::make_shared<vislab::RegularSteadyScalarField2f>())
        , mVelocity_u(std::make_shared<vislab::RegularSteadyScalarField2f>())
        , mVelocity_v(std::make_shared<vislab::RegularSteadyScalarField2f>())
        , mAcceleration_u(std::make_shared<vislab::RegularSteadyScalarField2f>())
        , mAcceleration_v(std::make_shared<vislab::RegularSteadyScalarField2f>())
        , densityTexture(std::make_shared<vislab::ColormapTexture>())
        , accuracy(1E-3)
        , iterations(1000)
        , stepSize(0.005 * std::sqrt((resolution.x() + resolution.y()) * 0.5))
        , pressureSolver(EPressureSolver::Iterative)
        , boundary(EBoundary::Closed)
    {
        mSpacing = (mDomain.max() - mDomain.min()).cwiseQuotient(mResolution.cast<double>());

        auto gridCentered = std::make_shared<vislab::RegularGrid2d>();
        gridCentered->setResolution(mResolution);
        gridCentered->setDomain(Eigen::AlignedBox2d(
            mDomain.min() + mSpacing * 0.5,
            mDomain.max() - mSpacing * 0.5));
        mDensity->setGrid(gridCentered);
        mDensity->setArray(std::make_shared<vislab::Array1f>());
        mDensity->getArray()->setSize(gridCentered->getResolution().prod());
        mPressure->setGrid(gridCentered);
        mPressure->setArray(std::make_shared<vislab::Array1f>());
        mPressure->getArray()->setSize(gridCentered->getResolution().prod());
        mDivergence->setGrid(gridCentered);
        mDivergence->setArray(std::make_shared<vislab::Array1f>());
        mDivergence->getArray()->setSize(gridCentered->getResolution().prod());

        auto gridStaggeredX = std::make_shared<vislab::RegularGrid2d>();
        gridStaggeredX->setResolution(mResolution + Eigen::Vector2i(1, 0));
        gridStaggeredX->setDomain(Eigen::AlignedBox2d(
            mDomain.min() + Eigen::Vector2d(0, mSpacing.y() * 0.5),
            mDomain.max() - Eigen::Vector2d(0, mSpacing.y() * 0.5)));
        mVelocity_u->setGrid(gridStaggeredX);
        mVelocity_u->setArray(std::make_shared<vislab::Array1f>());
        mVelocity_u->getArray()->setSize(gridStaggeredX->getResolution().prod());
        mAcceleration_u->setGrid(gridStaggeredX);
        mAcceleration_u->setArray(std::make_shared<vislab::Array1f>());
        mAcceleration_u->getArray()->setSize(gridStaggeredX->getResolution().prod());

        auto gridStaggeredY = std::make_shared<vislab::RegularGrid2d>();
        gridStaggeredY->setResolution(mResolution + Eigen::Vector2i(0, 1));
        gridStaggeredY->setDomain(Eigen::AlignedBox2d(
            mDomain.min() + Eigen::Vector2d(mSpacing.x() * 0.5, 0),
            mDomain.max() - Eigen::Vector2d(mSpacing.x() * 0.5, 0)));
        mVelocity_v->setGrid(gridStaggeredY);
        mVelocity_v->setArray(std::make_shared<vislab::Array1f>());
        mVelocity_v->getArray()->setSize(gridStaggeredY->getResolution().prod());
        mAcceleration_v->setGrid(gridStaggeredY);
        mAcceleration_v->setArray(std::make_shared<vislab::Array1f>());
        mAcceleration_v->getArray()->setSize(gridStaggeredY->getResolution().prod());

        mDensity_tmp    = std::shared_ptr<vislab::RegularSteadyScalarField2f>(mDensity->clone());
        mVelocity_u_tmp = std::shared_ptr<vislab::RegularSteadyScalarField2f>(mVelocity_u->clone());
        mVelocity_v_tmp = std::shared_ptr<vislab::RegularSteadyScalarField2f>(mVelocity_v->clone());

        densityTexture->scalarField = mDensity;
    }

    void FluidSolver::advance()
    {
        // apply source in density field
        applyDensitySource();

        // accumulate forces
        addBuoyancy();

        // advect the flow under application of the forces to the next time step. the flow becomes compressible.
        solveAdvection();

        // apply pressure-projection to make fluid incompressible
        solvePressure();

        // advect density in incompressible flow
        advectDensity();

        // for debugging, compute divergence
        computeDivergence();

        // reset forces
        mAcceleration_u->getArray()->setZero();
        mAcceleration_v->getArray()->setZero();

        // inform that texture has changed
        densityTexture->scalarFieldChanged.notify(densityTexture.get());
    }

    void FluidSolver::reset()
    {
        mDensity->getArray()->setZero();
        applyDensitySource();
        mPressure->getArray()->setZero();
        mDivergence->getArray()->setZero();
        mVelocity_u->getArray()->setZero();
        mVelocity_v->getArray()->setZero();
        mAcceleration_u->getArray()->setZero();
        mAcceleration_v->getArray()->setZero();

        // precompute sparse Cholesky factorization for Laplace operator on 2D grid
        Eigen::SparseMatrix<double> A(mResolution.prod(), mResolution.prod());
        buildLaplace2d(A);
        mPoissonFactorization.compute(A);
    }

    void FluidSolver::applyDensitySource()
    {
        Eigen::AlignedBox2d domain = mDensity->getDomain();
        Eigen::Vector2i resolution = mDensity->getGrid()->getResolution();
        Eigen::Vector2i low        = (sourceRegion.min() - domain.min()).cwiseQuotient(domain.max() - domain.min()).cwiseProduct(resolution.cast<double>() - Eigen::Vector2d::Ones()).cast<int>().cwiseMin(resolution - Eigen::Vector2i::Ones()).cwiseMax(Eigen::Vector2i::Zero());
        Eigen::Vector2i high       = (sourceRegion.max() - domain.min()).cwiseQuotient(domain.max() - domain.min()).cwiseProduct(resolution.cast<double>() - Eigen::Vector2d::Ones()).cast<int>().cwiseMin(resolution - Eigen::Vector2i::Ones()).cwiseMax(Eigen::Vector2i::Zero());

        for (int y = low.y(); y < high.y(); y++)
        {
            for (int x = low.x(); x < high.x(); x++)
            {
                float rho_ij = 1.f;
                mDensity->setVertexDataAt({ x, y }, rho_ij);
            }
        }
    }

    void FluidSolver::addBuoyancy()
    {
        Eigen::Vector2i resolution = mDensity->getGrid()->getResolution();
        float scaling              = 64.f / resolution.x();

        // approximate buoyancy simply from density
        for (int j = 1; j < resolution.y(); ++j)
        {
            for (int i = 0; i < resolution.x(); ++i)
            {
                float acc_v         = mAcceleration_v->getVertexDataAt({ i, j }).x();
                float density_below = mDensity->getVertexDataAt({ i, j - 1 }).x();
                float density_above = mDensity->getVertexDataAt({ i, j }).x();
                float density       = (density_below + density_above) / 2.f;
                acc_v += 0.1f * density * scaling;
                mAcceleration_v->setVertexDataAt({ i, j }, acc_v);
            }
        }
    }

    void FluidSolver::addAcceleration()
    {
        Eigen::Vector2i resu = mVelocity_u->getGrid()->getResolution();
        for (int j = 0; j < resu.y(); ++j)
        {
            for (int i = 0; i < resu.x(); ++i)
            {
                float vel_u = mVelocity_u->getVertexDataAt({ i, j }).x();
                float acc_u = mAcceleration_u->getVertexDataAt({ i, j }).x();
                vel_u += stepSize * acc_u;
                mVelocity_u->setVertexDataAt({ i, j }, vel_u);
            }
        }

        Eigen::Vector2i resv = mVelocity_v->getGrid()->getResolution();
        for (int j = 0; j < resv.y(); ++j)
        {
            for (int i = 0; i < resv.x(); ++i)
            {
                float vel_v = mVelocity_v->getVertexDataAt({ i, j }).x();
                float acc_v = mAcceleration_v->getVertexDataAt({ i, j }).x();
                vel_v += stepSize * acc_v;
                mVelocity_v->setVertexDataAt({ i, j }, vel_v);
            }
        }
    }

    void FluidSolver::solveAdvection()
    {
        // accelerate the velocity field with the accumulated forces
        addAcceleration();

        // advect fluid -> the fluid becomes compressible
        advectVelocity();

        // set boundary conditions on velocity. it would be smarter if the advection procedure meets these automatically. In this form here, we lose energy.
        switch (boundary)
        {
        case EBoundary::Open:
            setNormalNeumann();
            break;
        case EBoundary::Closed:
            setNormalDirichlet();
            setTangentialNoSlip();
            break;
        }
    }

    void FluidSolver::solvePressure()
    {

        // compute divergence as input for the pressure solver
        computeDivergence();

        // compute pressure
        switch (pressureSolver)
        {
        case EPressureSolver::Iterative:
            solvePoissonIterative();
            break;
        case EPressureSolver::Linear:
            solvePoissonLinear();
            break;
        }

        // apply pressure to correct velocity
        correctVelocity();
    }

    void FluidSolver::setNormalNeumann()
    {
        // x-velocity
        Eigen::Vector2i resu = mVelocity_u->getGrid()->getResolution();
        for (int y = 0; y < resu.y(); ++y)
        {
            mVelocity_u->setVertexDataAt({ 0, y }, mVelocity_u->getVertexDataAt({ 1, y }));
            mVelocity_u->setVertexDataAt({ resu.x() - 1, y }, mVelocity_u->getVertexDataAt({ resu.x() - 2, y }));
        }

        // y-velocity
        Eigen::Vector2i resv = mVelocity_v->getGrid()->getResolution();
        for (int x = 0; x < resv.x(); ++x)
        {
            mVelocity_v->setVertexDataAt({ x, 0 }, mVelocity_v->getVertexDataAt({ x, 1 }));
            mVelocity_v->setVertexDataAt({ x, resv.y() - 1 }, mVelocity_v->getVertexDataAt({ x, resv.y() - 2 }));
        }
    }

    void FluidSolver::setNormalDirichlet()
    {
        // x-velocity
        Eigen::Vector2i resu = mVelocity_u->getGrid()->getResolution();
        for (int y = 0; y < resu.y(); ++y)
        {
            mVelocity_u->setVertexDataAt({ 0, y }, 0);
            mVelocity_u->setVertexDataAt({ resu.x() - 1, y }, 0);
        }

        // y-velocity
        Eigen::Vector2i resv = mVelocity_v->getGrid()->getResolution();
        for (int x = 0; x < resv.x(); ++x)
        {
            mVelocity_v->setVertexDataAt({ x, 0 }, 0);
            mVelocity_v->setVertexDataAt({ x, resv.y() - 1 }, 0);
        }
    }

    void FluidSolver::setTangentialNoSlip()
    {
        // x-velocity
        Eigen::Vector2i resu = mVelocity_u->getGrid()->getResolution();
        for (int x = 0; x < resu.x(); ++x)
        {
            mVelocity_u->setVertexDataAt({ x, 0 }, 0);
            mVelocity_u->setVertexDataAt({ x, resu.y() - 1 }, 0);
        }

        // y-velocity
        Eigen::Vector2i resv = mVelocity_v->getGrid()->getResolution();
        for (int y = 0; y < resv.y(); ++y)
        {
            mVelocity_v->setVertexDataAt({ 0, y }, 0);
            mVelocity_v->setVertexDataAt({ resv.x() - 1, y }, 0);
        }
    }

    void FluidSolver::computeDivergence()
    {
        // calculate divergence
        Eigen::Vector2i res = mDivergence->getGrid()->getResolution();
        for (int y = 0; y < res.y(); ++y)
        {
            for (int x = 0; x < res.x(); ++x)
            {
                float xComponent = (mVelocity_u->getVertexDataAt({ x + 1, y }) - mVelocity_u->getVertexDataAt({ x, y })).x() / mSpacing.x();
                float yComponent = (mVelocity_v->getVertexDataAt({ x, y + 1 }) - mVelocity_v->getVertexDataAt({ x, y })).x() / mSpacing.y();
                float divergence = xComponent + yComponent;
                mDivergence->setVertexDataAt({ x, y }, divergence);
            }
        }
    }

    void FluidSolver::solvePoissonIterative()
    {
        float dx2      = mSpacing.prod();
        float residual = accuracy + 1; // initial residual
        float rho      = 1;

        for (int it = 0; residual > accuracy && it < iterations; ++it)
        {
            switch (boundary)
            {
            case EBoundary::Open:
                for (int j = 0; j < mResolution.y(); ++j)
                {
                    for (int i = 0; i < mResolution.x(); ++i)
                    {
                        float b = mDivergence->getVertexDataAt({ i, j }).x() / stepSize * rho; // right-hand
                        // TODO: update the pressure values
                        mPressure->setVertexDataAt({ i, j }, (1 / 4.0) * (((i != 0) ? mPressure->getVertexDataAt({ i - 1, j }).x() : 0) +
                                                                          ((i != mResolution.x() - 1) ? mPressure->getVertexDataAt({ i + 1, j }).x() : 0) +
                                                                          ((j != 0) ? mPressure->getVertexDataAt({ i, j - 1 }).x() : 0) +
                                                                          ((j != mResolution.y() - 1) ? mPressure->getVertexDataAt({ i, j + 1 }).x() : 0) - b * dx2));
                    }
                }
                break;
            case EBoundary::Closed:
                for (int j = 0; j < mResolution.y(); ++j)
                {
                    for (int i = 0; i < mResolution.x(); ++i)
                    {
                        float b = mDivergence->getVertexDataAt({ i, j }).x() / stepSize * rho; // right-hand
                        float k = 4;
                        // TODO: update the pressure values
                        if (i == 0 || i == mResolution.x() - 1)
                        {
                            k--;
                        }
                        if (j == 0 || j == mResolution.y() - 1)
                        {
                            k--;
                        }
                        mPressure->setVertexDataAt({ i, j }, (1 / k) * (((i != 0) ? mPressure->getVertexDataAt({ i - 1, j }).x() : 0) +
                                                                        ((i != mResolution.x() - 1) ? mPressure->getVertexDataAt({ i + 1, j }).x() : 0) +
                                                                        ((j != 0) ? mPressure->getVertexDataAt({ i, j - 1 }).x() : 0) +
                                                                        ((j != mResolution.y() - 1) ? mPressure->getVertexDataAt({ i, j + 1 }).x() : 0) - b * dx2));
                    }
                }
                break;
            }

            // Compute the new residual, i.e. the sum of the squares of the individual residuals (squared L2-norm)
            residual = 0;
            for (int j = 1; j < mResolution.y() - 1; ++j)
            {
                for (int i = 1; i < mResolution.x() - 1; ++i)
                {
                    float b = mDivergence->getVertexDataAt({ i, j }).x() / stepSize * rho; // right-hand
                    // TODO: compute the cell residual
                    float abc = b - ((((-4 * mPressure->getVertexDataAt({ i, j })) + (mPressure->getVertexDataAt({ i - 1, j })) + (mPressure->getVertexDataAt({ i + 1, j })) + (mPressure->getVertexDataAt({ i, j - 1 })) + (mPressure->getVertexDataAt({ i, j + 1 }))).x()) / (dx2));
                    residual += abc * abc;
                }
            }

            // Get the L2-norm of the residual
            residual = sqrt(residual);

            // We assume the accuracy is meant for the average L2-norm per grid cell
            residual /= (mResolution.x() - 2) * (mResolution.y() - 2);
        }
    }

    void FluidSolver::solvePoissonLinear()
    {
        float dx2 = mSpacing.prod();
        float rho = 1;

        // TODO: right hand side
        Eigen::VectorXd b(mResolution.x() * mResolution.y());
        for (int j = 0; j < mResolution.y(); ++j)
            for (int i = 0; i < mResolution.x(); ++i)
                b(j * (size_t)mResolution.x() + i) = mDivergence->getVertexDataAt({ i, j }).x() * dx2 * rho / stepSize; // ... set correct value
               

        // solve sparse linear SPD system
        Eigen::VectorXd pout = mPoissonFactorization.solve(b);

        // copy result to output
        for (int j = 0; j < mResolution.y(); ++j)
            for (int i = 0; i < mResolution.x(); ++i)
                mPressure->setVertexDataAt({ i, j }, pout(j * (size_t)mResolution.x() + i));
    }

    void FluidSolver::correctVelocity()
    {
        //     the staggered grid indices look like this
        //     ------------------------------
        //     |                            |
        //     |u_i-0.5,j      p_i,j        |u_i+0.5,j
        //     |                            |
        //     ------------------------------

        //     velocity u is actually stored at
        //     ------------------------------
        //     |                            |
        //     |u_i,j          p_i,j        |u_i+1,j
        //     |                            |
        //     ------------------------------
        //     hence, updating u_i,j needs p_i,j and p_i-1,j

        // Note: velocity u_{i+1/2} is practically stored at i+1, hence xV_{i} -= dt * (p_{i} - p_{i-1}) / dx
        for (int j = 1; j < mResolution.y() - 1; ++j)
            for (int i = 1; i < mResolution.x(); ++i)
            {
                // TODO: update u
                float u_corr = mVelocity_u->getVertexDataAt({ i, j }).x() - (stepSize * (mPressure->getVertexDataAt({ i, j }) - mPressure->getVertexDataAt({ i - 1, j })).x() / mSpacing.x());
                mVelocity_u->setVertexDataAt({ i, j }, u_corr);
            }

        // Same for velocity v_{i+1/2}.
        for (int j = 1; j < mResolution.y(); ++j)
            for (int i = 1; i < mResolution.x() - 1; ++i)
            {
                // TODO: update v
                float v_corr = mVelocity_v->getVertexDataAt({ i, j }).x() - (stepSize * (mPressure->getVertexDataAt({ i, j }) - mPressure->getVertexDataAt({ i, j - 1 })).x() / mSpacing.y());
                mVelocity_v->setVertexDataAt({ i, j }, v_corr);
            }
    }

    void FluidSolver::advectVelocity()
    {
        // Velocities live on the MAC grid
        // velocity_u resolution is: [0,m_res_x] x [0,m_res_y-1]
        // velocity_v resolution is: [0,m_res_x-1] x [0,m_res_y]

        //     the staggered grid indices look like this (with half offsets added already)
        //        v_i-1,j+1      v_i,j+1
        //     ------O-------------O--------
        //     |            |
        //     |            O uij  X pij
        //     |            |
        //     ------O-------------O--------
        //        v_i-1,j        v_i,j

        // Velocities (u), MAC grid
        for (int j = 0; j < mResolution.y(); ++j)
        {
            for (int i = 1; i < mResolution.x(); ++i)
            { // skip first and last row: those are determined by the boundary condition
                // TODO: Compute the velocity
                float last_x_velocity = (mVelocity_u->getVertexDataAt({ i, j })).x();                                                                                                                                                    // ... set correct value
                float last_y_velocity = ((mVelocity_v->getVertexDataAt({ i - 1, j + 1 }) + mVelocity_v->getVertexDataAt({ i, j + 1 }) + (mVelocity_v->getVertexDataAt({ i, j }) + mVelocity_v->getVertexDataAt({ i - 1, j }))).x()) / 4; // ... set correct value

                // TODO: Find the last position of the particle (in grid coordinates) using an Euler step
                float last_x = i - last_x_velocity * mResolution.x() * stepSize; // ... set correct value
                float last_y = j - last_y_velocity * mResolution.y() * stepSize; // ... set correct value

                // Make sure the coordinates are inside the boundaries
                if (last_x < 0)
                    last_x = 0;
                if (last_y < 0)
                    last_y = 0;
                if (last_x > mResolution.x() - 0)
                    last_x = mResolution.x() - 0;
                if (last_y > mResolution.y() - 1)
                    last_y = mResolution.y() - 1;

                // Determine corners for bilinear interpolation
                int x_low  = (int)last_x;
                int y_low  = (int)last_y;
                int x_high = std::min(x_low + 1, mResolution.x());
                int y_high = std::min(y_low + 1, mResolution.y() - 1);

                // Compute the interpolation weights
                float x_weight = last_x - x_low;
                float y_weight = last_y - y_low;

                // TODO: Bilinear interpolation
                float new_u = (mVelocity_u->getVertexDataAt({ x_high, y_high }).x() * y_weight * x_weight) + (mVelocity_u->getVertexDataAt({ x_high, y_low }).x() * (1 - y_weight) * x_weight) + (mVelocity_u->getVertexDataAt({ x_low, y_low }).x() * (1 - y_weight) * (1 - x_weight)) + (mVelocity_u->getVertexDataAt({ x_low, y_high }).x() * (1 - x_weight) * y_weight);
                mVelocity_u_tmp->setVertexDataAt({ i, j }, new_u);
            }
        }

        //     the staggered grid indices look like this (with half offsets added already)
        //     |              |
        //     O uij  X pij   | u_i+1,j
        //     |              |
        //     -------O--------
        //     |    v_i,j     |
        //     O uij-1        | u_i+1,j-1
        //     |              |

        // Velocities (v), MAC grid
        for (int j = 1; j < mResolution.y(); ++j)
        { // skip first and last column: those are determined by the boundary condition
            for (int i = 0; i < mResolution.x(); ++i)
            {
                // TODO: Compute the velocity
                float last_y_velocity = (mVelocity_v->getVertexDataAt({ i, j })).x();                                                                                                                                                    // ... set correct value
                float last_x_velocity = ((mVelocity_u->getVertexDataAt({ i, j }) + mVelocity_u->getVertexDataAt({ i + 1, j }) + (mVelocity_u->getVertexDataAt({ i + 1, j - 1 }) + mVelocity_u->getVertexDataAt({ i, j - 1 }))).x()) / 4; // ... set correct value                                                                                                                                                                                                                        // ... set correct value

                // TODO: Find the last position of the particle (in grid coordinates) using an Euler step
                float last_x = i - last_x_velocity * mResolution.x() * stepSize; // ... set correct value
                float last_y = j - last_y_velocity * mResolution.y() * stepSize; // ... set correct value

                // Make sure the coordinates are inside the boundaries
                if (last_x < 0)
                    last_x = 0;
                if (last_y < 0)
                    last_y = 0;
                if (last_x > mResolution.x() - 1)
                    last_x = mResolution.x() - 1;
                if (last_y > mResolution.y() - 0)
                    last_y = mResolution.y() - 0;

                // Determine corners for bilinear interpolation
                int x_low  = (int)last_x;
                int y_low  = (int)last_y;
                int x_high = std::min(x_low + 1, mResolution.x() - 1);
                int y_high = std::min(y_low + 1, mResolution.y() - 0);

                // Compute the interpolation weights
                float x_weight = last_x - x_low;
                float y_weight = last_y - y_low;

                // TODO: Bilinear interpolation
                float new_v = (mVelocity_v->getVertexDataAt({ x_high, y_high }).x() * y_weight * x_weight) + (mVelocity_v->getVertexDataAt({ x_high, y_low }).x() * (1 - y_weight) * x_weight) + (mVelocity_v->getVertexDataAt({ x_low, y_low }).x() * (1 - y_weight) * (1 - x_weight)) + (mVelocity_v->getVertexDataAt({ x_low, y_high }).x() * (1 - x_weight) * y_weight);
                mVelocity_v_tmp->setVertexDataAt({ i, j }, new_v);
            }
        }

        // Copy the values in temp to the original buffers
        for (int j = 0; j < mResolution.y(); ++j)
            for (int i = 1; i < mResolution.x(); ++i)
                mVelocity_u->setVertexDataAt({ i, j }, mVelocity_u_tmp->getVertexDataAt({ i, j }));
        for (int j = 1; j < mResolution.y(); ++j)
            for (int i = 0; i < mResolution.x(); ++i)
                mVelocity_v->setVertexDataAt({ i, j }, mVelocity_v_tmp->getVertexDataAt({ i, j }));
    }

    void FluidSolver::advectDensity()
    {
        // Densities live on the grid centers, the velocities on the MAC grid
        // Separate their computation to avoid confusion

        //     the staggered grid indices look like this (with half offsets added already)
        //     ------------------------------
        //     |                            |
        //     O u_i,j       X p_i,j        O u_i+1,j
        //     |                            |
        //     ------------------------------

        // Densities, grid centers
        for (int j = 0; j < mResolution.y(); ++j)
        {
            for (int i = 0; i < mResolution.x(); ++i)
            {
                // TODO: Compute the velocity
                float last_x_velocity = ((mVelocity_u->getVertexDataAt({ i, j }) + mVelocity_u->getVertexDataAt({ i + 1, j })).x()) / 2; // ... set correct value
                float last_y_velocity = ((mVelocity_v->getVertexDataAt({ i, j }) + mVelocity_v->getVertexDataAt({ i, j + 1 })).x()) / 2; // ... set correct value

                // TODO: Find the last position of the particle (in grid coordinates) using an Euler step
                float last_x = i - last_x_velocity * mResolution.x() * stepSize; // ... set correct value
                float last_y = j - last_y_velocity * mResolution.y() * stepSize; // ... set correct value

                // Make sure the coordinates are inside the boundaries
                const float offset = 0.0; // a trick to fight the dissipation through boundaries is to sample with a small offset
                if (last_x < offset)
                    last_x = offset;
                if (last_y < offset)
                    last_y = offset;
                if (last_x > mResolution.x() - 1 - offset)
                    last_x = mResolution.x() - 1 - offset;
                if (last_y > mResolution.y() - 1 - offset)
                    last_y = mResolution.y() - 1 - offset;

                // Determine corners for bilinear interpolation
                int x_low  = (int)last_x;
                int y_low  = (int)last_y;
                int x_high = std::min(x_low + 1, mResolution.x() - 1);
                int y_high = std::min(y_low + 1, mResolution.y() - 1);

                // Compute the interpolation weights
                float x_weight = last_x - x_low;
                float y_weight = last_y - y_low;

                // TODO: Bilinear interpolation
                float new_dens = (mDensity->getVertexDataAt({ x_high, y_high }).x() * y_weight * x_weight) + (mDensity->getVertexDataAt({ x_high, y_low }).x() * (1 - y_weight) * x_weight) + (mDensity->getVertexDataAt({ x_low, y_low }).x() * (1 - y_weight) * (1 - x_weight)) + (mDensity->getVertexDataAt({ x_low, y_high }).x() * (1 - x_weight) * y_weight);
                mDensity_tmp->setVertexDataAt({ i, j }, new_dens);
            }
        }

        // Copy the values in temp to the original buffers
        for (int j = 0; j < mResolution.y(); ++j)
            for (int i = 0; i < mResolution.x(); ++i)
                mDensity->setVertexDataAt({ i, j }, mDensity_tmp->getVertexDataAt({ i, j }));
    }

 void FluidSolver::buildLaplace2d(Eigen::SparseMatrix<double>& A)
    {
        Eigen::Vector2i resolution = mPressure->getGrid()->getResolution();
        int nx                     = resolution.x(); // Number of grid points along x-axis
        int ny                     = resolution.y(); // Number of grid points along y-axis 
        std::list<Eigen::Triplet<double>> coeff;
        for (int j = 0; j < resolution.y(); ++j)
            for (int i = 0; i < resolution.x(); ++i)
            {
                
                coeff.push_back(Eigen::Triplet<double>(j * resolution.x() + i, j * resolution.x() + i, -4)); 
                
                switch (boundary)
                {

                case EBoundary::Open: // smooth continuation
                                      // example for 4x4 matrix
                                      //   -4  1  0  0  1  0  0  0  0  0  0  0  0  0  0  0
                                      //    1 -4  1  0  0  1  0  0  0  0  0  0  0  0  0  0
                                      // 	0  1 -4  1  0  0  1  0  0  0  0  0  0  0  0  0
                                      // 	0  0  1 -4  0  0  0  1  0  0  0  0  0  0  0  0
                                      //    1  0  0  0 -4  1  0  0  1  0  0  0  0  0  0  0
                                      // 	0  1  0  0  1 -4  1  0  0  1  0  0  0  0  0  0
                                      // 	0  0  1  0  0  1 -4  1  0  0  1  0  0  0  0  0
                                      // 	0  0  0  1  0  0  1 -4  0  0  0  1  0  0  0  0
                                      // 	0  0  0  0  1  0  0  0 -4  1  0  0  1  0  0  0
                                      // 	0  0  0  0  0  1  0  0  1 -4  1  0  0  1  0  0
                                      // 	0  0  0  0  0  0  1  0  0  1 -4  1  0  0  1  0
                                      // 	0  0  0  0  0  0  0  1  0  0  1 -4  0  0  0  1
                                      // 	0  0  0  0  0  0  0  0  1  0  0  0 -4  1  0  0
                                      // 	0  0  0  0  0  0  0  0  0  1  0  0  1 -4  1  0
                                      // 	0  0  0  0  0  0  0  0  0  0  1  0  0  1 -4  1
                                      // 	0  0  0  0  0  0  0  0  0  0  0  1  0  0  1 -4

                    if (i + 1 < nx)
                    {
                        coeff.push_back(Eigen::Triplet<double>(j * resolution.x() + i, j * resolution.x() + i + 1, 1));
                    }

                    // Left neighbor
                    if (i - 1 >= 0)
                    {
                        coeff.push_back(Eigen::Triplet<double>(j * resolution.x() + i, j * resolution.x() + i - 1, 1));
                    }

                    // Bottom neighbor
                    if (j + 1 < ny)
                    {
                        coeff.push_back(Eigen::Triplet<double>(j * resolution.x() + i, j * resolution.x() + i + nx, 1));
                    }

                    // Top neighbor
                    if (j - 1 >= 0)
                    {
                        coeff.push_back(Eigen::Triplet<double>(j * resolution.x() + i, j * resolution.x() + i - nx, 1));
                    }
                    break;
    

                case EBoundary::Closed: // forward/backward difference on boundary   (same pattern of -1's, but the main diagonal contains the number of -1's per row)
                                        // -2  1  0  0  1  0  0  0  0  0  0  0  0  0  0  0
                                        //  1 -3  1  0  0  1  0  0  0  0  0  0  0  0  0  0
                                        //	0  1 -3  1  0  0  1  0  0  0  0  0  0  0  0  0
                                        //	0  0  1 -2  0  0  0  1  0  0  0  0  0  0  0  0
                                        //  1  0  0  0 -3  1  0  0  1  0  0  0  0  0  0  0
                                        //	0  1  0  0  1 -4  1  0  0  1  0  0  0  0  0  0
                                        //	0  0  1  0  0  1 -4  1  0  0  1  0  0  0  0  0
                                        //	0  0  0  1  0  0  1 -3  0  0  0  1  0  0  0  0
                                        //	0  0  0  0  1  0  0  0 -3  1  0  0  1  0  0  0
                                        //	0  0  0  0  0  1  0  0  1 -4  1  0  0  1  0  0
                                        //	0  0  0  0  0  0  1  0  0  1 -4  1  0  0  1  0
                                        //	0  0  0  0  0  0  0  1  0  0  1 -3  0  0  0  1
                                        //	0  0  0  0  0  0  0  0  1  0  0  0 -2  1  0  0
                                        //	0  0  0  0  0  0  0  0  0  1  0  0  1 -3  1  0
                                        //	0  0  0  0  0  0  0  0  0  0  1  0  0  1 -3  1
                                        //	0  0  0  0  0  0  0  0  0  0  0  1  0  0  1 -2

                        // For the 1's off diagonals


                        coeff.push_back(Eigen::Triplet<double>(j * resolution.x() + i, j * resolution.x() + std::max(i - 1, 0), 1));
                 
                        coeff.push_back(Eigen::Triplet<double>(j * resolution.x() + i, j * resolution.x() + std::min(i + 1, resolution.x() - 1), 1));

                        // For the values of the main diagonals
                   
                        coeff.push_back(Eigen::Triplet<double>(j * resolution.x() + i, std::max(j - 1, 0) * resolution.x() + i, 1));
                    
                        coeff.push_back(Eigen::Triplet<double>(j * resolution.x() + i, std::min(j + 1, resolution.y() - 1) * resolution.x() + i, 1));
                    
             

                    break;
                 
                }
            }
       // A.resize(resolution.x() * resolution.y(), resolution.x() * resolution.y());  
        A.setFromTriplets(coeff.begin(), coeff.end());
    }
}

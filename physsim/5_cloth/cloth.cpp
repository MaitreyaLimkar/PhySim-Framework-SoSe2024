#include "cloth.hpp"

#include <vislab/core/array.hpp>
#include <vislab/graphics/mesh.hpp>


namespace physsim
{
    Cloth::Cloth(const Eigen::Vector2i& resolution, const ETopology& topology, const Eigen::Vector3d& origin, const Eigen::Vector3d& axis1, const Eigen::Vector3d& axis2)
        : mResolution(resolution)
        , mTopology(topology)
        , mOrigin(origin)
        , mAxis1(axis1)
        , mAxis2(axis2)
    {
        // allocate mesh
        mesh = std::make_shared<vislab::Mesh>();

        // create buffer to hold the spring positions
        mesh->positions = std::make_shared<vislab::Array3f>();
        mesh->positions->setSize(resolution.prod());

        // create buffer to hold the spring velocities
        velocities = std::make_shared<vislab::Array3f>();
        velocities->setSize(resolution.prod());

        // create buffer to accumulate internal spring forces in
        springForces = std::make_shared<vislab::Array3f>();
        springForces->setSize(resolution.prod());

        // create index buffer
        mesh->indices = std::make_shared<vislab::Array3ui>();
        mesh->indices->setSize((resolution.x() - 1) * (resolution.y() - 1) * 2);
        int index = 0;
        for (int iy = 0; iy < resolution.y() - 1; ++iy)
            for (int ix = 0; ix < resolution.x() - 1; ++ix)
            {
                int i00 = (iy + 0) * resolution.x() + (ix + 0);
                int i01 = (iy + 1) * resolution.x() + (ix + 0);
                int i10 = (iy + 0) * resolution.x() + (ix + 1);
                int i11 = (iy + 1) * resolution.x() + (ix + 1);
                mesh->indices->setValue(index++, Eigen::Vector3u(i00, i01, i10));
                mesh->indices->setValue(index++, Eigen::Vector3u(i01, i11, i10));
            }

        reset();
    }

    void Cloth::advance(double stepSize)
    {
        double Lx  = mAxis1.norm();
        double Ly  = mAxis2.norm();
        double Lxy = (mAxis1 + mAxis2).norm();
        springForces->setZero();

#ifndef _DEBUG
#pragma omp parallel for
#endif
        for (int64_t i = 0; i < mResolution.prod(); ++i)
        {
            // early out if the start of the spring is fixed
            if (mFixed.find(i) != mFixed.end())
                continue;

            // vertex grid indices
            int ix = i % mResolution.x();
            int iy = i / mResolution.x();

            // add spring forces according to the spring topology
            if (mTopology == ETopology::Structural)
            {
                if (ix > 0)
                    addSpringForce(Eigen::Vector2i(ix - 1, iy), Eigen::Vector2i(ix, iy), Lx);
                if (ix < mResolution.x() - 1)
                    addSpringForce(Eigen::Vector2i(ix + 1, iy), Eigen::Vector2i(ix, iy), Lx);
                if (iy > 0)
                    addSpringForce(Eigen::Vector2i(ix, iy - 1), Eigen::Vector2i(ix, iy), Ly);
                if (iy < mResolution.y() - 1)
                    addSpringForce(Eigen::Vector2i(ix, iy + 1), Eigen::Vector2i(ix, iy), Ly);
            }
            else if (mTopology == ETopology::Diagonal)
            {
                if (ix > 0)
                    addSpringForce(Eigen::Vector2i(ix - 1, iy), Eigen::Vector2i(ix, iy), Lx);
                if (ix < mResolution.x() - 1)
                    addSpringForce(Eigen::Vector2i(ix + 1, iy), Eigen::Vector2i(ix, iy), Lx);
                if (iy > 0)
                    addSpringForce(Eigen::Vector2i(ix, iy - 1), Eigen::Vector2i(ix, iy), Ly);
                if (iy < mResolution.y() - 1)
                    addSpringForce(Eigen::Vector2i(ix, iy + 1), Eigen::Vector2i(ix, iy), Ly);
                if (ix > 0 && iy > 0)
                    addSpringForce(Eigen::Vector2i(ix - 1, iy - 1), Eigen::Vector2i(ix, iy), Lxy);
                if (ix > 0 && iy < mResolution.y() - 1)
                    addSpringForce(Eigen::Vector2i(ix - 1, iy + 1), Eigen::Vector2i(ix, iy), Lxy);
                if (ix < mResolution.x() - 1 && iy > 0)
                    addSpringForce(Eigen::Vector2i(ix + 1, iy - 1), Eigen::Vector2i(ix, iy), Lxy);
                if (ix < mResolution.x() - 1 && iy < mResolution.y() - 1)
                    addSpringForce(Eigen::Vector2i(ix + 1, iy + 1), Eigen::Vector2i(ix, iy), Lxy);
            }
        }

        // symplectic euler update
#ifndef _DEBUG
#pragma omp parallel for
#endif
        for (int64_t i = 1; i < mResolution.prod(); ++i)
        {
            // skip the corners, since those are fixed
            if (i == mResolution.x() - 1)
                continue;

            // get current position and velocities
            Eigen::Vector3d x = mesh->positions->getValueDouble(i);
            Eigen::Vector3d v = velocities->getValueDouble(i);

            // get the accumulated internal spring forces
            Eigen::Vector3d f_int = springForces->getValueDouble(i);

            // TODO: compute damping force and external force
            Eigen::Vector3d f_damp = -damping * v;
            Eigen::Vector3d f_ext  = mass * gravity;

            // TODO: sum up internal, damping, and external force, and convert this to an acceleration.
            Eigen::Vector3d a = (f_int + f_damp + f_ext)/mass;

            // TODO: symplectic Euler integration of position and velocity.
            v += stepSize * a;
            x += stepSize * v;
            
            // TODO: set the new position and the new velocity
            mesh->positions->setValueDouble(i, x);
            velocities->setValueDouble(i, v);            
           
        }

        // recompute the vertex normals
        mesh->recomputeVertexNormals();

        // notify that mesh positions have changed
        mesh->positionsChanged.notify(mesh.get());
    }

    void Cloth::reset()
    {
        for (int iy = 0; iy < mResolution.y(); ++iy)
            for (int ix = 0; ix < mResolution.x(); ++ix)
            {
                mesh->positions->setValueDouble(iy * mResolution.x() + ix, mOrigin + ix * mAxis1 + iy * mAxis2);
            }
        velocities->setZero();
    }

    const Eigen::Vector2i& Cloth::resolution() const
    {
        return mResolution;
    }

    void Cloth::pin(const Eigen::Vector2i& gridIndex)
    {
        mFixed.insert(gridIndex.y() * mResolution.x() + gridIndex.x());
    }

    void Cloth::unpin(const Eigen::Vector2i& gridIndex)
    {
        mFixed.erase(gridIndex.y() * mResolution.x() + gridIndex.x());
    }

    void Cloth::addSpringForce(const Eigen::Vector2i& startPoint, const Eigen::Vector2i& endPoint, const double& L)
    {
        // linear array indices for the start and end point
        Eigen::Index linearIndex_startPoint = startPoint.y() * mResolution.x() + startPoint.x();
        Eigen::Index linearIndex_endPoint   = endPoint.y() * mResolution.x() + endPoint.x();

        // get start and end position of spring, as well as current sum of spring forces.
        Eigen::Vector3d startPos = mesh->positions->getValueDouble(linearIndex_startPoint);
        Eigen::Vector3d endPos   = mesh->positions->getValueDouble(linearIndex_endPoint);
        Eigen::Vector3d f_int   = springForces->getValueDouble(linearIndex_endPoint);

        // TODO: compute spring direction and the displacement
        Eigen::Vector3d springVect  = endPos - startPos; 
        double L_x   = springVect.norm(); 
        Eigen::Vector3d springForce = -stiffness * (L_x - L) * springVect.normalized(); 
        
        // TODO: add spring force to f_int
        f_int += springForce;
                       
        // TODO: store the new accumulate spring force f_int in the array springForces
        springForces->setValueDouble(linearIndex_endPoint, f_int); 
      
    }
}

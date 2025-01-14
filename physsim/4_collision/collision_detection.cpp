#include "collision_detection.hpp"

#include "gilbert_johnson_keerthi.hpp"
#include "transformed_mesh.hpp"

#include <vislab/core/array.hpp>

namespace physsim
{
    CollisionDetection::CollisionDetection(const std::vector<std::shared_ptr<RigidBody>>& objects)
        : mObjects(objects)
    {
    }

    void CollisionDetection::computeCollisionDetection(EBroadPhaseMethod broadPhaseMethod, ENarrowPhaseMethod narrowPhaseMethod, double eps, double stepSize)
    {
        clearDataStructures();
        computeBroadPhase(broadPhaseMethod);
        computeNarrowPhase(narrowPhaseMethod);
        applyImpulse(eps, stepSize);
    }

    void CollisionDetection::computeBroadPhase(EBroadPhaseMethod broadPhaseMethod)
    {
        // compute possible collisions
        mOverlappingBodys.clear();

        switch (broadPhaseMethod)
        {
        case EBroadPhaseMethod::None:
        {
            for (size_t i = 0; i < mObjects.size(); i++)
            {
                for (size_t j = i + 1; j < mObjects.size(); j++)
                {
                    if (mObjects[i]->type() == RigidBody::EType::Dynamic || mObjects[j]->type() == RigidBody::EType::Dynamic)
                        mOverlappingBodys.push_back(std::make_pair(i, j));
                }
            }
            break;
        }

        case EBroadPhaseMethod::AABB:
        {
            // compute bounding boxes
            std::vector<Eigen::AlignedBox3d> aabbs(mObjects.size());
            for (size_t i = 0; i < aabbs.size(); i++)
            {
                aabbs[i] = mObjects[i]->shape()->worldBounds();
            }
            for (size_t i = 0; i < mObjects.size(); i++)
            {
                for (size_t j = i + 1; j < mObjects.size(); j++)
                {
                    // add pair of objects to possible collision if their bounding boxes overlap
                    if (mObjects[i]->type() == RigidBody::EType::Dynamic || mObjects[j]->type() == RigidBody::EType::Dynamic)
                    {
                        if (aabbs[i].intersects(aabbs[j]))
                        {
                            mOverlappingBodys.push_back(std::make_pair(i, j));
                        }
                    }
                }
            }
            break;
        }

        case EBroadPhaseMethod::SweepAndPrune:
        {
            // TODO: compute bounding boxes and create intervals on the 3 main axes
            std::vector<Eigen::AlignedBox3d> aabbs(mObjects.size());
            for (size_t i = 0; i < aabbs.size(); i++)
            {
                aabbs[i] = mObjects[i]->shape()->worldBounds();
            }
            using interval = std::tuple<double, double, int>;
            std::vector<interval> intervals[3];
            for (int i = 0; i < 3; ++i)
                intervals[i].resize(aabbs.size());

            for (size_t i = 0; i < aabbs.size(); i++)
            {
                intervals[0][i] = std::make_tuple(aabbs[i].min()[0], aabbs[i].max()[0], i);
                intervals[1][i] = std::make_tuple(aabbs[i].min()[1], aabbs[i].max()[1], i); 
                intervals[2][i] = std::make_tuple(aabbs[i].min()[2], aabbs[i].max()[2], i);
            }

            // TODO: sort intervals in ascending order by beginning of interval
            // std::get<0>(intervals[0][i])
            std::sort(intervals[0].begin(), intervals[0].end());
            std::sort(intervals[1].begin(), intervals[1].end());
            std::sort(intervals[2].begin(), intervals[2].end());

            // TODO: iterate and place overlaps in a set
            std::set<std::pair<int, int>> overlap[3];
            for (int n = 0; n < 3; n++)
            {
                for (size_t i = 0; i < aabbs.size(); i++)
                {
                    double end_i = std::get<1>(intervals[n][i]);
                    for (size_t j = i + 1; j < aabbs.size(); j++)
                    {
                        double start_j = std::get<0>(intervals[n][j]);
                        if (start_j < end_i)
                        {
                            int ind_i = std::get<2>(intervals[n][i]);
                            int ind_j = std::get<2>(intervals[n][j]);
                            if (ind_j > ind_i)
                            {
                                overlap[n].insert(std::make_pair(ind_j, ind_i));
                            }
                            else
                            {
                                overlap[n].insert(std::make_pair(ind_i, ind_j));
                            }
                        }
                        else
                        {
                            break;
                        }
                    }
                }
            }

            std::set<std::pair<int, int>> commonOverlaps;
            std::set_intersection(overlap[0].begin(), overlap[0].end(),
                                  overlap[1].begin(), overlap[1].end(),
                                  std::inserter(commonOverlaps, commonOverlaps.begin()));

            std::set<std::pair<int, int>> finalOverlaps;
            std::set_intersection(commonOverlaps.begin(), commonOverlaps.end(),
                                  overlap[2].begin(), overlap[2].end(),
                                  std::inserter(finalOverlaps, finalOverlaps.begin()));

            mOverlappingBodys.assign(finalOverlaps.begin(), finalOverlaps.end());


            // TODO: compute bounding boxes and create intervals on the 3 main axes

            
            // TODO: sort intervals in ascending order by beginning of interval

            // TODO: iterate and place overlaps in a set
            
            // TODO: grab elements that occurred in all containers for the narrow test
            
            // TODO: pass the intersections on to the narrow phase
            break;
        }
        }
    }

    void CollisionDetection::computeNarrowPhase(ENarrowPhaseMethod narrowPhaseMethod)
    {
        // iterate through all pairs of possible collisions
        for (auto overlap : mOverlappingBodys)
        {
            mPenetratingVertices.clear();
            mPenetratingEdges.clear();
            std::vector<Contact> temp_contacts[2];
            // compute intersection of a with b first and intersection
            // of b with a and store results in temp_contacts
            for (int switcher = 0; switcher < 2; switcher++)
            {
                RigidBody* a =
                    mObjects[(!switcher) ? overlap.first
                                         : overlap.second]
                        .get();
                RigidBody* b =
                    mObjects[(!switcher) ? overlap.second
                                         : overlap.first]
                        .get();

                const vislab::Triangle* a_tri = dynamic_cast<const vislab::Triangle*>(a->shape().get());
                const vislab::Triangle* b_tri = dynamic_cast<const vislab::Triangle*>(b->shape().get());

                TransformedMesh a_tmesh(
                    a_tri->transform,
                    *a_tri->mesh->positions.get(),
                    *a_tri->mesh->indices.get());

                TransformedMesh b_tmesh(
                    b_tri->transform,
                    *b_tri->mesh->positions.get(),
                    *b_tri->mesh->indices.get());

                switch (narrowPhaseMethod)
                {
                    // exhaustive
                case ENarrowPhaseMethod::Exhaustive:
                {
                    // iterate through all faces of first object
                    for (Eigen::Index face = 0; face < a_tmesh.indices.getSize(); face++)
                    {
                        // iterate through all edges of given face
                        for (size_t j = 0; j < 3; j++)
                        {
                            int start = a_tmesh.getPrimitive(face)[j];
                            int end   = a_tmesh.getPrimitive(face)[(j + 1) % 3];

                            // check if there is a collision
                            Eigen::Vector3d v_start   = a_tmesh.getVertex(start);
                            Eigen::Vector3d v_end     = a_tmesh.getVertex(end);
                            Eigen::Vector3d origin    = v_start;
                            Eigen::Vector3d direction = v_end - v_start;
                            vislab::Ray3d ray(
                                origin,
                                direction.normalized(),
                                0, direction.norm());
                            ContactType ct = isColliding(ray, b_tri);

                            // find collision and check for duplicates
                            switch (ct)
                            {
                            case ContactType::Vertex_Face:
                            {
                                auto si = b_tri->closestHit(ray);
                                if (si.isValid())
                                {
                                    double penetration_depth_start = si.normal.dot(v_start - si.position);
                                    double penetration_depth_end   = si.normal.dot(v_end - si.position);
                                    int penetrating_vertex =
                                        penetration_depth_start < penetration_depth_end ? start : end;

                                    if (mPenetratingVertices.find(penetrating_vertex) == mPenetratingVertices.end())
                                    {
                                        Contact temp_col;
                                        temp_col.n     = si.normal.stableNormalized();
                                        temp_col.p     = si.position;
                                        temp_col.a     = a;
                                        temp_col.b     = b;
                                        temp_col.depth = std::min(penetration_depth_start, penetration_depth_end);
                                        temp_col.type  = ContactType::Vertex_Face;
                                        temp_contacts[switcher].push_back(temp_col);
                                        mPenetratingVertices.insert(penetrating_vertex);
                                    }
                                }
                                break;
                            }
                            case ContactType::Edge_Edge:
                            {
                                int orderedStart = std::min(start, end);
                                int orderedEnd   = std::max(start, end);
                                auto pair        = std::make_pair(orderedStart,
                                                                  orderedEnd);
                                // if not already in set
                                if (mPenetratingEdges.find(pair) == mPenetratingEdges.end())
                                {
                                    Contact temp_col =
                                        findEdgeEdgeCollision(
                                            a_tmesh.getVertex(orderedStart),
                                            a_tmesh.getVertex(orderedEnd), b_tmesh);
                                    temp_col.a    = a;
                                    temp_col.b    = b;
                                    temp_col.type = ContactType::Edge_Edge;
                                    temp_contacts[switcher].push_back(temp_col);
                                    mPenetratingEdges.insert(pair);
                                }
                                break;
                            }
                            case ContactType::None:
                            {
                                break;
                            }
                            }
                        }
                    }
                    break;
                }
                case ENarrowPhaseMethod::GilbertJohnsonKeerthi: // TODO: implement narrow-phase collision detection method. Here, a starting point for a GJK algorithm is written.
                {
                    Contact contact = GilbertJohnsonKeerthi::findCollisionGJK(a_tmesh, b_tmesh);
                    if (contact.type != ContactType::None)
                    {
                        contact.a = a;
                        contact.b = b;
                        temp_contacts[switcher].push_back(contact);
                    }
                }
                }
            }
            // look for vertexFace
            bool found = false;
            for (int i = 0; i < 2; i++)
            {
                for (const auto& cont : temp_contacts[i])
                {
                    if (cont.type == ContactType::Vertex_Face)
                    {
                        mContacts.push_back(cont);
                        found = true;
                        break;
                    }
                }
                if (found)
                {
                    break;
                }
            }
            if (found)
            {
                continue;
            }

            // take single edgeedge if possible
            if (temp_contacts[0].size() > 0 &&
                temp_contacts[0].size() < temp_contacts[1].size())
            {
                for (int i = 0; i < temp_contacts[0].size(); i++)
                {
                    mContacts.push_back(temp_contacts[0][i]);
                }
            }
            else if (temp_contacts[1].size() > 0 &&
                     temp_contacts[0].size() >
                         temp_contacts[1].size())
            {
                for (int i = 0; i < temp_contacts[1].size(); i++)
                {
                    mContacts.push_back(temp_contacts[1][i]);
                }
            }
            else if (temp_contacts[0].size() > 0)
            {
                for (int i = 0; i < temp_contacts[0].size(); i++)
                {
                    mContacts.push_back(temp_contacts[0][i]);
                }
            }
            else if (temp_contacts[1].size() > 0)
            {
                for (int i = 0; i < temp_contacts[1].size(); i++)
                {
                    mContacts.push_back(temp_contacts[1][i]);
                }
            }
        }
    }

    ContactType CollisionDetection::isColliding(const vislab::Ray3d& ray, const vislab::Shape* shape)
    {
        int32_t hits = shape->countHits(ray);

        // if hits is odd then the ray enters through one face and
        // does not leave again through another, hence the starting point is
        // inside the object, if the number of intersections between start and
        // end of edge are even, then the edge entering through one and leaving
        // through another face, hence the edge is intersecting the object
        ContactType ret = ContactType::None;
        if (hits % 2 == 1)
        {
            ret = ContactType::Vertex_Face;
        }
        if (hits % 2 == 0 && hits > 0)
        {
            ret = ContactType::Edge_Edge;
        }
        return ret;
    }

    Contact CollisionDetection::findEdgeEdgeCollision(const Eigen::Vector3d& start, const Eigen::Vector3d& end, const TransformedMesh& mesh)
    {
        double minDist = std::numeric_limits<double>::infinity();
        Contact ret;
        for (int i = 0; i < mesh.indices.getSize(); i++)
        {
            Eigen::Vector3d a = mesh.getVertex(mesh.getPrimitive(i)[0]);
            Eigen::Vector3d b = mesh.getVertex(mesh.getPrimitive(i)[1]);
            Eigen::Vector3d c = mesh.getVertex(mesh.getPrimitive(i)[2]);

            Eigen::Vector3d n_face = -(b - a).cross(c - a).normalized();

            for (int j = 0; j < 3; j++)
            {
                Eigen::Vector3d s = mesh.getVertex(mesh.getPrimitive(i)[j]);
                Eigen::Vector3d e = mesh.getVertex(mesh.getPrimitive(i)[(j + 1) % 3]);

                Eigen::Vector3d ea = end - start;
                Eigen::Vector3d eb = e - s;
                Eigen::Vector3d n  = (ea).cross(eb); // direction of shortest distance
                double distance    = n.dot(start - s) / n.norm();

                Eigen::Vector3d plane_normal = n.cross(eb).normalized();
                double t                     = (s - start).dot(plane_normal) / (ea.dot(plane_normal));
                if (n_face.dot(n) < 0 && distance < 0 && -distance < minDist &&
                    t >= 0 && t <= 1)
                {
                    ret.ea    = ea;
                    ret.eb    = eb;
                    ret.n     = n.stableNormalized();
                    ret.p     = start + t * ea;
                    minDist   = -distance;
                    ret.depth = minDist;
                }
            }
        }
        return ret;
    }

    void CollisionDetection::applyImpulse(double eps, double stepSize)
    {
        // see here fore more information: https://en.wikipedia.org/wiki/Collision_response
        if (mContacts.empty())
            return;

        // compute impulse for all contacts
        for (auto& contact : mContacts)
        {
            // compute relative velocity and skip if the bodies are already moving apart
            Eigen::Vector3d vrel_vec = contact.b->velocity(contact.p) - contact.a->velocity(contact.p);
            double vrel              = contact.n.dot(vrel_vec);
            if (vrel < 0)
            {
                // bodies are moving apart
                continue;
            }

            // vectors ra and rb from the centers of mass to the contact point (in world space)
            Eigen::Vector3d ra = contact.p - contact.a->position();
            Eigen::Vector3d rb = contact.p - contact.b->position();
            Eigen::Vector3d n  = contact.n;

            // TODO: compute impulse response
            double ma_inv = contact.a->massInverse();
            double mb_inv = contact.b->massInverse();
            Eigen::Matrix3d Ia_inv = contact.a->inertiaBodyInverse();
            Eigen::Matrix3d Ib_inv = contact.b->inertiaBodyInverse();
            double den_1           = n.dot((Ia_inv * (ra.cross(n))).cross(ra));
            double den_2           = n.dot((Ib_inv * (rb.cross(n))).cross(rb));

            double j_num = -(1 + eps) * vrel;
            double j_den = ma_inv + mb_inv + den_1 + den_2;
            double j = j_num / j_den;

            // TODO: apply impulse forces to the bodies at the contact point
            Eigen::Vector3d f_impulse = j * contact.n / stepSize; 
            contact.a->applyForce(-f_impulse, contact.p);
            contact.b->applyForce(f_impulse, contact.p);
        }
    }

    void CollisionDetection::clearDataStructures()
    {
        mPenetratingEdges.clear();
        mPenetratingVertices.clear();
        mOverlappingBodys.clear();
        mContacts.clear();
    }
}

#include "gilbert_johnson_keerthi.hpp"

#include "transformed_mesh.hpp"

#include <vislab/core/array.hpp>

namespace physsim
{
    Contact GilbertJohnsonKeerthi::findCollisionGJK(const TransformedMesh& A, const TransformedMesh& B)
    {
        Contact ret;
        ret.type = ContactType::None;

        SupportVec a, b, c, d;        // set of coordinates of the simplex S
        Eigen::Vector3d dir(1, 0, 0); // initial direction

        // TODO: find first point on Minkoswki difference (store in c)
        c = support(A, B, dir);
        
        // TODO: revert search direction
        dir = -dir;

        // TODO: find second point on Minkoswki difference (store in b)
        b = support(A, B, dir);

        // TODO: early out, if we cannot contain the origin (dot product)
        
        if (b.P.dot(dir) < 0)
        {
            return ret;
        }

        // the next direction is perpendicular to line towards origin
        dir = (c.P - b.P).cross(-b.P).cross(c.P - b.P);
        if (dir.isZero())
        { // if origin is on the line try another direction
            dir = (c.P - b.P).cross(Eigen::Vector3d(1, 0, 0));
            if (dir.isZero())
                dir = (c.P - b.P).cross(Eigen::Vector3d(0, 1, 0)); // if we were unlucky again, pick another direction
        }
        int dim = 2; // current dimensionality of the simplex

        const int GJK_MAX_NUM_ITERATIONS = 64; // threshold depends on complexity of the meshes. this is rather a safety measure to not get stuck in infinite loops
        for (int iterations = 0; iterations < GJK_MAX_NUM_ITERATIONS; iterations++)
        {
            // TODO: find next point (store in a)
            a = support(A, B, dir);

            // TODO: early out, if we cannot contain the origin (dot product)
            if (a.P.dot(dir) < 0)
            {

                return ret;
            }

            // virtually add point a to simplex (technically, we have always allocate the memory for it)
            dim++;

            // find nearest simplex. returns true if a collision was found
            if (nearestSimplex(a, b, c, d, dim, dir))
            {
                // call the EPA algorithm to locate the contact
                return EPA(a, b, c, d, A, B);
            }
        }
        return ret;
    }

    Eigen::Vector3d GilbertJohnsonKeerthi::support(const TransformedMesh& mesh, const Eigen::Vector3d& d)
    {
        // TODO: find the vertex of V that is farthest into direction d.
        int n = mesh.vertices.getSize();
        Eigen::Index ind;
        double largest_dot = -std::numeric_limits<double>::max();

        for (int i = 0; i < n; i++)
        {
            Eigen::Vector3d farthestVertex = mesh.getVertex(i);
            double dot_product             = farthestVertex.dot(d);
            if (dot_product > largest_dot)
            { 
                largest_dot = dot_product;
                ind         = i;
            }
            
        }
        
        return mesh.getVertex(ind);
    }

    GilbertJohnsonKeerthi::SupportVec GilbertJohnsonKeerthi::support(const TransformedMesh& A, const TransformedMesh& B, const Eigen::Vector3d& d)
    {
        Eigen::Vector3d sup_A = support(A, d);
        Eigen::Vector3d sup_B = support(B, -d);
        return SupportVec(sup_A, sup_B, sup_A - sup_B);
    }

    bool GilbertJohnsonKeerthi::nearestSimplex(SupportVec& a, SupportVec& b, SupportVec& c, SupportVec& d, int& dim, Eigen::Vector3d& search_dir)
    {
        if (dim == 3)
        {
            Eigen::Vector3d n  = (b.P - a.P).cross(c.P - a.P); // triangle normal
            Eigen::Vector3d AO = -a.P;                         // direction to origin

            // Determine which feature is closest to origin, make that the new simplex
            dim = 2;
            if ((b.P - a.P).cross(n).dot(AO) > 0)
            { // Closest to edge AB
                c          = a;
                search_dir = (b.P - a.P).cross(AO).cross(b.P - a.P);
                return false;
            }
            if (n.cross(c.P - a.P).dot(AO) > 0)
            { // Closest to edge AC
                b          = a;
                search_dir = (c.P - a.P).cross(AO).cross(c.P - a.P);
                return false;
            }

            dim = 3;
            if (n.dot(AO) > 0)
            { // Above triangle
                d          = c;
                c          = b;
                b          = a;
                search_dir = n;
                return false;
            }
            // else Below triangle
            d          = b;
            b          = a;
            search_dir = -n;
            return false;
        }
        else if (dim == 4)
        {
            // Get normals of three new faces
            Eigen::Vector3d ABC = (b.P - a.P).cross(c.P - a.P);
            Eigen::Vector3d ACD = (c.P - a.P).cross(d.P - a.P);
            Eigen::Vector3d ADB = (d.P - a.P).cross(b.P - a.P);
            Eigen::Vector3d AO  = -a.P; // dir to origin
            dim                 = 3;
            if (ABC.dot(AO) > 0)
            { // In front of ABC
                d          = c;
                c          = b;
                b          = a;
                search_dir = ABC;
                return false;
            }
            if (ACD.dot(AO) > 0)
            { // In front of ACD
                b          = a;
                search_dir = ACD;
                return false;
            }
            if (ADB.dot(AO) > 0)
            { // In front of ADB
                c          = d;
                d          = b;
                b          = a;
                search_dir = ADB;
                return false;
            }
            // else inside tetrahedron; enclosed!
            return true;
        }
        return false;
    }

    Contact GilbertJohnsonKeerthi::EPA(const SupportVec& a, const SupportVec& b, const SupportVec& c, const SupportVec& d, const TransformedMesh& A, const TransformedMesh& B)
    {
        const double EPA_TOLERANCE        = 0.0001;
        const int EPA_MAX_NUM_FACES       = 64;
        const int EPA_MAX_NUM_LOOSE_EDGES = 32;
        const int EPA_MAX_NUM_ITERATIONS  = 64;
        SupportVec faces[EPA_MAX_NUM_FACES][4]; // Array of faces, each with 3 verts and a normal

        // Init with final simplex from GJK
        faces[0][0] = a;
        faces[0][1] = b;
        faces[0][2] = c;
        faces[0][3] = SupportVec(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), (b.P - a.P).cross(c.P - a.P).normalized()); // ABC
        faces[1][0] = a;
        faces[1][1] = c;
        faces[1][2] = d;
        faces[1][3] = SupportVec(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), (c.P - a.P).cross(d.P - a.P).normalized()); // ACD
        faces[2][0] = a;
        faces[2][1] = d;
        faces[2][2] = b;
        faces[2][3] = SupportVec(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), (d.P - a.P).cross(b.P - a.P).normalized()); // ADB
        faces[3][0] = b;
        faces[3][1] = d;
        faces[3][2] = c;
        faces[3][3] = SupportVec(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), (d.P - b.P).cross(c.P - b.P).normalized()); // BDC

        int num_faces = 4;
        int closest_face;
        for (int iterations = 0; iterations < EPA_MAX_NUM_ITERATIONS; iterations++)
        {
            // Find face that's closest to origin
            double min_dist = faces[0][0].P.dot(faces[0][3].P);
            closest_face    = 0;
            for (int i = 1; i < num_faces; i++)
            {
                double dist = faces[i][0].P.dot(faces[i][3].P);
                if (dist < min_dist)
                {
                    min_dist     = dist;
                    closest_face = i;
                }
            }

            // search normal to face that's closest to origin
            Eigen::Vector3d search_dir = faces[closest_face][3].P;
            SupportVec p               = support(A, B, search_dir);

            if (p.P.dot(search_dir) - min_dist < EPA_TOLERANCE)
            {

                double dist_to_resolve = p.P.dot(search_dir);
                Eigen::Vector3d v0     = faces[closest_face][1].P - faces[closest_face][0].P;
                Eigen::Vector3d v1     = faces[closest_face][2].P - faces[closest_face][0].P;
                Eigen::Vector3d v2     = faces[closest_face][3].P * dist_to_resolve - faces[closest_face][0].P;
                double d00(v0.dot(v0)), d01(v0.dot(v1)), d11(v1.dot(v1)), d20(v2.dot(v0)), d21(v2.dot(v1));
                double denom                 = d00 * d11 - d01 * d01;
                double bary_v                = (d11 * d20 - d01 * d21) / denom;
                double bary_w                = (d00 * d21 - d01 * d20) / denom;
                double bary_u                = 1.0 - bary_v - bary_w;
                Eigen::Vector3d contact_on_A = bary_u * faces[closest_face][0].A + bary_v * faces[closest_face][1].A + bary_w * faces[closest_face][2].A;
                // Eigen::Vector3d contact_on_B = bary_u * faces[closest_face][0].B + bary_v * faces[closest_face][1].B + bary_w * faces[closest_face][2].B;
                Contact ret;
                ret.n     = -faces[closest_face][3].P.normalized();
                ret.p     = contact_on_A;
                ret.type  = ContactType::Vertex_Face; // note that we don't distinguish between edge-edge and vertex-face collisions!
                ret.depth = dist_to_resolve;
                return ret;
            }

            SupportVec loose_edges[EPA_MAX_NUM_LOOSE_EDGES][2]; // keep track of edges we need to fix after removing faces
            int num_loose_edges = 0;

            // Find all triangles that are facing p
            for (int i = 0; i < num_faces; i++)
            {
                if (faces[i][3].P.dot(p.P - faces[i][0].P) > 0) // triangle i faces p, remove it
                {
                    // Add removed triangle's edges to loose edge list.
                    // If it's already there, remove it (both triangles it belonged to are gone)
                    for (int j = 0; j < 3; j++) // Three edges per face
                    {
                        SupportVec current_edge[2] = { faces[i][j], faces[i][(j + 1) % 3] };
                        bool found_edge            = false;
                        for (int k = 0; k < num_loose_edges; k++) // Check if current edge is already in list
                        {
                            if (loose_edges[k][1].P == current_edge[0].P && loose_edges[k][0].P == current_edge[1].P)
                            {
                                // Edge is already in the list, remove it
                                loose_edges[k][0] = loose_edges[num_loose_edges - 1][0]; // Overwrite current edge
                                loose_edges[k][1] = loose_edges[num_loose_edges - 1][1]; // with last edge in list
                                num_loose_edges--;
                                found_edge = true;
                                k          = num_loose_edges; // exit loop because edge can only be shared once
                            }
                        } // endfor loose_edges

                        if (!found_edge)
                        { // add current edge to list
                            if (num_loose_edges >= EPA_MAX_NUM_LOOSE_EDGES)
                                break;
                            loose_edges[num_loose_edges][0] = current_edge[0];
                            loose_edges[num_loose_edges][1] = current_edge[1];
                            num_loose_edges++;
                        }
                    }

                    // Remove triangle i from list
                    faces[i][0] = faces[num_faces - 1][0];
                    faces[i][1] = faces[num_faces - 1][1];
                    faces[i][2] = faces[num_faces - 1][2];
                    faces[i][3] = faces[num_faces - 1][3];
                    num_faces--;
                    i--;
                } // endif p can see triangle i
            }     // endfor num_faces

            // Reconstruct polytope with p added
            for (int i = 0; i < num_loose_edges; i++)
            {
                if (num_faces >= EPA_MAX_NUM_FACES)
                    break;
                faces[num_faces][0]   = loose_edges[i][0];
                faces[num_faces][1]   = loose_edges[i][1];
                faces[num_faces][2]   = p;
                faces[num_faces][3].P = (loose_edges[i][0].P - loose_edges[i][1].P).cross(loose_edges[i][0].P - p.P).normalized();

                // Check for wrong normal to maintain CCW winding
                double bias = 0.000001; // in case dot result is only slightly < 0 (because origin is on face)
                if (faces[num_faces][0].P.dot(faces[num_faces][3].P) + bias < 0)
                {
                    SupportVec temp       = faces[num_faces][0];
                    faces[num_faces][0]   = faces[num_faces][1];
                    faces[num_faces][1]   = temp;
                    faces[num_faces][3].P = -faces[num_faces][3].P;
                }
                num_faces++;
            }
        } // End for iterations
        Contact ret;
        return ret;
    }

}

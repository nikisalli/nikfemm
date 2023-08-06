#ifndef NIK_VIRTUAL_MESH_HPP
#define NIK_VIRTUAL_MESH_HPP

#include <cstdint>
#include <map>

#include "../algebra/build_coo.hpp"
#include "../algebra/simple_vector.hpp"
#include "../drawing/drawing.hpp"
#include "../geometry/common.hpp"

namespace nikfemm {
    struct Elem {
        int verts[3];

        int& operator[](int i) {
            return verts[i];
        }
    };

    // this was copied from the triangle library, I changed some types to be able to reinterpret_cast and use the data in a more intuitive way
    struct MeshData {
        Vector *pointlist = nullptr;                                               /* In / out */
        TRI_REAL *pointattributelist = nullptr;                                      /* In / out */
        int *pointmarkerlist = nullptr;                                          /* In / out */
        int numberofpoints;                                            /* In / out */
        int numberofpointattributes;                                   /* In / out */

        Elem *trianglelist = nullptr;                                             /* In / out */
        TRI_REAL *triangleattributelist = nullptr;                                   /* In / out */
        TRI_REAL *trianglearealist = nullptr;                                         /* In only */
        int *neighborlist = nullptr;                                             /* Out only */
        int numberoftriangles;                                         /* In / out */
        int numberofcorners;                                           /* In / out */
        int numberoftriangleattributes;                                /* In / out */

        int *segmentlist = nullptr;                                              /* In / out */
        int *segmentmarkerlist = nullptr;                                        /* In / out */
        int numberofsegments;                                          /* In / out */

        TRI_REAL *holelist = nullptr;                        /* In / pointer to array copied out */
        int numberofholes;                                      /* In / copied out */

        TRI_REAL *regionlist = nullptr;                      /* In / pointer to array copied out */
        int numberofregions;                                    /* In / copied out */

        int *edgelist = nullptr;                                                 /* Out only */
        int *edgemarkerlist = nullptr;            /* Not used with Voronoi diagram; out only */
        TRI_REAL *normlist = nullptr;                /* Used only with Voronoi diagram; out only */
        int numberofedges;                                             /* Out only */

        int hullsize;                                                  /* Out only */

        inline Vector getElemBarycenter(uint32_t elem_id) {
            Elem myelem = trianglelist[elem_id];
            Vector p1 = pointlist[myelem[0]];
            Vector p2 = pointlist[myelem[1]];
            Vector p3 = pointlist[myelem[2]];
            return (p1 + p2 + p3) / 3.0;
        }

        inline double getArea(uint32_t elem_id) {
            return abs(getOrientedArea(elem_id));
        }

        inline double getOrientedArea(uint32_t elem_id) {
            return 0.5 * getDoubleOrientedArea(elem_id);
        }

        inline double getDoubleOrientedArea(uint32_t elem_id) {
            Elem myelem = trianglelist[elem_id];
            return getDoubleOrientedArea(myelem[0], myelem[1], myelem[2]);
        }

        inline double getArea(uint32_t v1, uint32_t v2, uint32_t v3) {
            return abs(getOrientedArea(v1, v2, v3));
        }

        inline double getOrientedArea(uint32_t v1, uint32_t v2, uint32_t v3) {
            return 0.5 * getDoubleOrientedArea(v1, v2, v3);
        }

        inline double getDoubleOrientedArea(uint32_t v1, uint32_t v2, uint32_t v3) {
            Vector p1 = pointlist[v1];
            Vector p2 = pointlist[v2];
            Vector p3 = pointlist[v3];
            return (p2 - p1) ^ (p3 - p1);
        }

        inline bool pointInElem(uint32_t elem_id, Vector point) {
            Elem myelem = trianglelist[elem_id];
            return pointInTriangle(point, pointlist[myelem[0]], pointlist[myelem[1]], pointlist[myelem[2]]);
        }
    };
    template <typename Prop>
    struct Mesh {
        Vector center = Vector(0, 0);
        double radius = 0;
        Prop default_prop;
        double epsilon;
        double max_triangle_area;

        MeshData data;
        Drawing<Prop> drawing;

        Mesh(double max_triangle_area);
        Mesh();
        ~Mesh();

        void plotRend(cv::Mat* image, double width, double height);
        void plotToFile(uint32_t width, uint32_t height, std::string filename);
        void plot(uint32_t width, uint32_t height);
        void mesh();
        void addKelvinBoundaryConditions(uint32_t boundary_points);
        void kelvinTransformCentered();
        void computeEpsilon();
    };

    template <typename Prop>
    Mesh<Prop>::Mesh(double max_triangle_area) {
        this->max_triangle_area = max_triangle_area;
    }

    // templated member functions must be defined in the header file
    template <typename Prop>
    Mesh<Prop>::Mesh() {
        max_triangle_area = 1e-0;
    }

    template <typename Prop>
    Mesh<Prop>::~Mesh() {
        // we don't know how triangle allocated the memory
        if (data.pointlist != nullptr) free(data.pointlist);
        if (data.pointattributelist != nullptr) free(data.pointattributelist);
        if (data.pointmarkerlist != nullptr) free(data.pointmarkerlist);
        if (data.trianglelist != nullptr) free(data.trianglelist);
        if (data.triangleattributelist != nullptr) free(data.triangleattributelist);
        if (data.neighborlist != nullptr) free(data.neighborlist);
        if (data.segmentlist != nullptr) free(data.segmentlist);
        if (data.segmentmarkerlist != nullptr) free(data.segmentmarkerlist);
        if (data.holelist != nullptr) free(data.holelist);
        if (data.edgelist != nullptr) free(data.edgelist);
        if (data.edgemarkerlist != nullptr) free(data.edgemarkerlist);
    }

    template <typename Prop>
    void Mesh<Prop>::plotRend(cv::Mat* image, double width, double height) {
        // get mesh enclosing rectangle
        float min_x = -1.1 * radius;
        float min_y = -1.1 * radius;
        float max_x = 1.1 * radius;
        float max_y = 1.1 * radius;

        // object to window ratio
        float ratio = 0.9;

        // x scale factor to loosely fit mesh in window (equal in x and y)
        float x_scale = ratio * width / std::max(max_x - min_x, max_y - min_y);
        // y scale factor to loosely fit mesh in window
        float y_scale = ratio * height / std::max(max_x - min_x, max_y - min_y);
        // x offset to center mesh in window
        float x_offset = 0.5 * width - 0.5 * (max_x + min_x) * x_scale;
        // y offset to center mesh in window
        float y_offset = 0.5 * height - 0.5 * (max_y + min_y) * y_scale;

        // render
        for (uint32_t i = 0; i < data.numberoftriangles; i++) {
            // get the triangle
            Elem e = data.trianglelist[i];
            // get the vertices
            Vector v1 = data.pointlist[e[0]];
            Vector v2 = data.pointlist[e[1]];
            Vector v3 = data.pointlist[e[2]];

            if (Vector::distance(v1, Vector(0, 0)) > radius + epsilon ||
                Vector::distance(v2, Vector(0, 0)) > radius + epsilon || 
                Vector::distance(v3, Vector(0, 0)) > radius + epsilon) {
                continue;
            }

            // draw mesh edges
            // draw lines from v1 to v2
            cv::line(*image, cv::Point(x_scale * static_cast<float>(v1.x) + x_offset, y_scale * static_cast<float>(v1.y) + y_offset),
                            cv::Point(x_scale * static_cast<float>(v2.x) + x_offset, y_scale * static_cast<float>(v2.y) + y_offset),
                            cv::Scalar(255, 255, 255), 1);
            cv::line(*image, cv::Point(x_scale * static_cast<float>(v2.x) + x_offset, y_scale * static_cast<float>(v2.y) + y_offset),
                            cv::Point(x_scale * static_cast<float>(v3.x) + x_offset, y_scale * static_cast<float>(v3.y) + y_offset),
                            cv::Scalar(255, 255, 255), 1);
            cv::line(*image, cv::Point(x_scale * static_cast<float>(v3.x) + x_offset, y_scale * static_cast<float>(v3.y) + y_offset),
                            cv::Point(x_scale * static_cast<float>(v1.x) + x_offset, y_scale * static_cast<float>(v1.y) + y_offset),
                            cv::Scalar(255, 255, 255), 1);
        }

        // draw the geometry
        // draw the segments
        for (DrawingSegment s : drawing.segments) {
            cv::line(*image, cv::Point(x_scale * drawing.points[s.p1].x + x_offset,
                                       y_scale * drawing.points[s.p1].y + y_offset), 
                             cv::Point(x_scale * drawing.points[s.p2].x + x_offset, 
                                       y_scale * drawing.points[s.p2].y + y_offset),
                             cv::Scalar(255, 255, 255), 1);
        }

        for (Vector p : drawing.points) {
            cv::circle(*image, cv::Point(x_scale * p.x + x_offset, y_scale * p.y + y_offset), 2, cv::Scalar(255, 255, 255), -1);
        }

        // draw the regions
        for (DrawingRegion r : drawing.regions) {
            Vector pos = r.first;
            // draw white cross
            cv::line(*image, cv::Point(x_scale * pos.x + x_offset - 5, y_scale * pos.y + y_offset),
                                cv::Point(x_scale * pos.x + x_offset + 5, y_scale * pos.y + y_offset),
                                cv::Scalar(0, 0, 255), 1);
            cv::line(*image, cv::Point(x_scale * pos.x + x_offset, y_scale * pos.y + y_offset - 5),
                                cv::Point(x_scale * pos.x + x_offset, y_scale * pos.y + y_offset + 5),
                                cv::Scalar(0, 0, 255), 1);
        }
    }

    template <typename Prop>
    void Mesh<Prop>::plot(uint32_t width, uint32_t height) {
        // create the image
        cv::Mat image = cv::Mat::zeros(height, width, CV_8UC3);

        // render the mesh
        plotRend(&image, width, height);

        // show the image
        cv::imshow("drawing", image);
        // continue if image is closed
        cv::waitKey(0);
    }

    template <typename Prop>
    void Mesh<Prop>::plotToFile(uint32_t width, uint32_t height, std::string filename) {
        // create the image
        cv::Mat image = cv::Mat::zeros(height, width, CV_8UC3);

        // render the mesh
        plotRend(&image, width, height);

        // save the image
        cv::imwrite(filename, image);
    }

    template <typename Prop>
    void Mesh<Prop>::computeEpsilon() {
        // compute epsilon to compare point distances
        double epsilon = std::numeric_limits<double>::infinity();
        for (uint32_t i = 0; i < data.numberoftriangles; i++) {
            Elem myelem = data.trianglelist[i];
            double x1 = data.pointlist[myelem[0]].x;
            double y1 = data.pointlist[myelem[0]].y;
            double x2 = data.pointlist[myelem[1]].x;
            double y2 = data.pointlist[myelem[1]].y;
            double x3 = data.pointlist[myelem[2]].x;
            double y3 = data.pointlist[myelem[2]].y;
            epsilon = std::min(epsilon, Vector::distance(Vector(x1, y1), Vector(x2, y2)));
            epsilon = std::min(epsilon, Vector::distance(Vector(x1, y1), Vector(x3, y3)));
            epsilon = std::min(epsilon, Vector::distance(Vector(x2, y2), Vector(x3, y3)));
        }
        epsilon *= 0.5;
        this->epsilon = epsilon;
        nloginfo("epsilon: %.17g", epsilon);
    }

    template <typename Prop>
    void Mesh<Prop>::mesh() {
        triangulateio in;

        /* set up the input */
        struct TempSegment {
            uint32_t id1;
            uint32_t id2;
        };

        in.numberofpoints = drawing.points.size();
        in.pointlist = (TRI_REAL*)malloc(in.numberofpoints * 2 * sizeof(TRI_REAL));
        in.numberofsegments = drawing.segments.size();
        in.segmentlist = (int*)malloc(in.numberofsegments * 2 * sizeof(int));

        for (uint32_t i = 0; i < drawing.points.size(); i++) {
            in.pointlist[2 * i] = drawing.points[i].x;
            in.pointlist[2 * i + 1] = drawing.points[i].y;
        }

        uint32_t i = 0;
        for (auto s : drawing.segments) {
            in.segmentlist[2 * i] = s.p1;
            in.segmentlist[2 * i + 1] = s.p2;
            i++;
        }

        in.pointmarkerlist = NULL;
        in.numberofpointattributes = 0;
        in.pointattributelist = NULL;
        in.segmentmarkerlist = NULL;
        in.numberofholes = 0;
        in.holelist = NULL;

        in.numberofregions = drawing.regions.size();
        in.regionlist = (TRI_REAL*)malloc(in.numberofregions * 4 * sizeof(TRI_REAL));
        i = 0;

        for (auto r : drawing.regions) {
            // nloginfo("adding region %f %f %lu\n", r.first.x, r.first.y, r.second);
            in.regionlist[4 * i] = r.first.x;
            in.regionlist[4 * i + 1] = r.first.y;
            in.regionlist[4 * i + 2] = r.second;
            in.regionlist[4 * i + 3] = 0;
            i++;
        }
        // printf("----------------------\n");
        char switches[30];
        sprintf(switches, "pzq33AQa%.17g", max_triangle_area);
        

        /* Make necessary initializations so that Triangle can return a */
        /*   triangulation in `mid' and a voronoi diagram in `vorout'.  */

        data.pointlist = (Vector *) NULL;            /* Not needed if -N switch used. */
        /* Not needed if -N switch used or number of point attributes is zero: */
        data.pointattributelist = (TRI_REAL *) NULL;
        data.pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
        data.trianglelist = (Elem *) NULL;          /* Not needed if -E switch used. */
        /* Not needed if -E switch used or number of triangle attributes is zero: */
        data.triangleattributelist = (TRI_REAL *) NULL;
        data.neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
        /* Needed only if segments are output (-p or -c) and -P not used: */
        data.segmentlist = (int *) NULL;
        /* Needed only if segments are output (-p or -c) and -P and -B not used: */
        data.segmentmarkerlist = (int *) NULL;
        data.edgelist = (int *) NULL;             /* Needed only if -e switch used. */
        data.edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */

        triangulate(switches, &in, reinterpret_cast<triangulateio*>(&data), NULL);
        free(in.pointlist);
        free(in.segmentlist);
        free(in.regionlist);
    }

    template <typename Prop>
    void Mesh<Prop>::kelvinTransformCentered() {
        // kelvin transform each Vertex<Prop> in the mesh
        /*

        x* = (R^2 / |x|^2) * x

        */
        const double max_radius_coeff = 100;
        double R_squared = radius * radius;
        double max_x = radius / max_radius_coeff;
        // printf("R^2 = %f\n", R_squared);

        for (uint32_t i = 0; i < data.numberofpoints; i++) {
            if (data.pointmarkerlist[i] == 0) {  // if it's not a boundary point
                Vector v = data.pointlist[i];
                double mag_squared = v.x * v.x + v.y * v.y;
                double scale = R_squared / mag_squared;
                // printf("v = (%f, %f) -> (%f, %f), mag = %f, scale = %f", v->p.x, v->p.y, v->p.x * scale, v->p.y * scale, mag_squared, scale);
                double dist = Vector::distance(v, center);
                if (dist < max_x) {
                // if (false) {
                    nloginfo("kelvin transform too large");
                    data.pointlist[i] = v * ((radius * max_radius_coeff) / dist);
                } else {
                    data.pointlist[i] = v * scale;
                }
            }
        }
    }

    template <typename Prop>
    void Mesh<Prop>::addKelvinBoundaryConditions(uint32_t boundary_points) {
        // kelvin mesh = km
        // this mesh = tm
        // kelvin mesh boundary vertices = kmb
        // this mesh boundary vertices = tmb
        // kelvin mesh vertices = kmv
        // this mesh vertices = tmv

        Mesh<Prop> kelvin_mesh;

        Circle boundary_circle = Circle(Vector(0, 0), radius);
        kelvin_mesh.drawing.drawCircle(boundary_circle, boundary_points);
        // add region near the edge of the circle
        kelvin_mesh.drawing.drawRegion(Vector(0, 0), default_prop);

        kelvin_mesh.mesh();
        kelvin_mesh.center = Vector(0, 0);
        kelvin_mesh.radius = radius;

        // merge km into tm
        // transform km
        kelvin_mesh.kelvinTransformCentered();
        // kelvin_mesh.plot();

        // merge km into tm
        // for each vertex in km, find the corresponding vertex in tm and build a map from kmb to tmb

        // 1. create map from kelvin mesh boundary vertices to this mesh boundary vertices
        // 2. calculate new this mesh point number
        // 3. reallocate memory for this mesh
        // 4. add kelvin mesh non boundary vertices to this mesh
        // 5. create map from old kelvin mesh vertices minus boundary to new mesh vertices
        // 6. calculate new this mesh triangle number
        // 7. reallocate memory for this mesh
        // 8. use map to update triangle vertex indices
        // 9. add kelvin mesh triangles to this mesh

        std::map<uint32_t, uint32_t> kelvin_to_this;
        uint32_t kelvin_number_of_boundary_points = 0;
        for (uint32_t i = 0; i < kelvin_mesh.data.numberofpoints; i++) {
            if (kelvin_mesh.data.pointmarkerlist[i]) {
                for (uint32_t j = 0; j < data.numberofpoints; j++) {
                    if (data.pointmarkerlist[j]) {
                        if (data.pointlist[j].x == kelvin_mesh.data.pointlist[i].x && data.pointlist[j].y == kelvin_mesh.data.pointlist[i].y) {
                            kelvin_to_this[i] = j;
                            kelvin_number_of_boundary_points++;
                            break;
                        }
                    }
                }
            }
        }

        // calculate new this mesh point number
        uint32_t new_point_number = data.numberofpoints + kelvin_mesh.data.numberofpoints - kelvin_number_of_boundary_points;

        // reallocate memory for this mesh
        data.pointlist = (Vector *) realloc(data.pointlist, new_point_number * sizeof(Vector));
        // I actually don't need these
        // data.pointattributelist = (TRI_REAL *) realloc(data.pointattributelist, new_point_number * sizeof(TRI_REAL));
        // data.pointmarkerlist = (int *) realloc(data.pointmarkerlist, new_point_number * sizeof(int));

        // add kelvin mesh non boundary vertices to this mesh
        for (uint32_t i = 0; i < kelvin_mesh.data.numberofpoints; i++) {
            if (!kelvin_mesh.data.pointmarkerlist[i]) {
                data.pointlist[data.numberofpoints] = kelvin_mesh.data.pointlist[i];
                // data.pointmarkerlist[data.numberofpoints] = 0;
                kelvin_to_this[i] = data.numberofpoints;
                data.numberofpoints++;
            }
        }

        // calculate new this mesh triangle number
        uint32_t new_triangle_number = data.numberoftriangles + kelvin_mesh.data.numberoftriangles;

        // reallocate memory for this mesh
        data.trianglelist = (Elem *) realloc(data.trianglelist, new_triangle_number * sizeof(Elem));
        data.triangleattributelist = (TRI_REAL *) realloc(data.triangleattributelist, new_triangle_number * sizeof(TRI_REAL));

        // find the id of the default prop in this mesh
        uint32_t default_prop_id = drawing.getRegionId(default_prop);

        // use map to update triangle vertex indices
        for (uint32_t i = 0; i < kelvin_mesh.data.numberoftriangles; i++) {
            // check if point in map else leave it alone
            for (uint8_t j = 0; j < 3; j++) {
                if (kelvin_to_this.count(kelvin_mesh.data.trianglelist[i][j])) {
                    data.trianglelist[data.numberoftriangles][j] = kelvin_to_this[kelvin_mesh.data.trianglelist[i][j]];
                }
            }
            data.triangleattributelist[data.numberoftriangles] = default_prop_id;
            data.numberoftriangles++;
        }

        data.numberoftriangles = new_triangle_number;
    }
}

#endif
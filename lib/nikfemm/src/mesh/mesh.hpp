#ifndef NIK_VIRTUAL_MESH_HPP
#define NIK_VIRTUAL_MESH_HPP

#include <cstdint>

#include "../algebra/coo.hpp"
#include "../algebra/simple_vector.hpp"
#include "../drawing/drawing.hpp"

namespace nikfemm {
    struct Elem {
        int verts[3];
    };

    // this was copied from the triangle library, I changed some types to be able to reinterpret_cast and use the data in a more intuitive way
    struct MeshData {
        Point *pointlist;                                               /* In / out */
        TRI_REAL *pointattributelist;                                      /* In / out */
        int *pointmarkerlist;                                          /* In / out */
        int numberofpoints;                                            /* In / out */
        int numberofpointattributes;                                   /* In / out */

        Elem *trianglelist;                                             /* In / out */
        TRI_REAL *triangleattributelist;                                   /* In / out */
        TRI_REAL *trianglearealist;                                         /* In only */
        int *neighborlist;                                             /* Out only */
        int numberoftriangles;                                         /* In / out */
        int numberofcorners;                                           /* In / out */
        int numberoftriangleattributes;                                /* In / out */

        int *segmentlist;                                              /* In / out */
        int *segmentmarkerlist;                                        /* In / out */
        int numberofsegments;                                          /* In / out */

        TRI_REAL *holelist;                        /* In / pointer to array copied out */
        int numberofholes;                                      /* In / copied out */

        TRI_REAL *regionlist;                      /* In / pointer to array copied out */
        int numberofregions;                                    /* In / copied out */

        int *edgelist;                                                 /* Out only */
        int *edgemarkerlist;            /* Not used with Voronoi diagram; out only */
        TRI_REAL *normlist;                /* Used only with Voronoi diagram; out only */
        int numberofedges;                                             /* Out only */

        int hullsize;                                                  /* Out only */
    };

    template <typename Prop>
    struct Mesh {
        Point center = Point(0, 0);
        double radius = 0;
        Prop default_prop;

        MeshData data;
        Drawing<Prop> drawing;

        Mesh();
        ~Mesh();

        void plot();
        void mesh();
        void addKelvinBoundaryConditions();
        void kelvinTransformCentered();
    };

    // templated member functions must be defined in the header file
    template <typename Prop>
    Mesh<Prop>::Mesh() {

    }

    template <typename Prop>
    Mesh<Prop>::~Mesh() {

    }

    template <typename Prop>
    void Mesh<Prop>::plot() {
        if (SDL_Init(SDL_INIT_EVERYTHING) != 0) {
            nexit("error initializing SDL: %s\n");
        }
        SDL_Window* win = SDL_CreateWindow("GAME", // creates a window
                                        SDL_WINDOWPOS_CENTERED,
                                        SDL_WINDOWPOS_CENTERED,
                                        2000, 2000, 0);
            
        // triggers the program that controls
        // your graphics hardware and sets flags
        Uint32 render_flags = SDL_RENDERER_ACCELERATED;
    
        // creates a renderer to render our images
        SDL_Renderer* rend = SDL_CreateRenderer(win, -1, 0);

        // clears the window
        SDL_RenderClear(rend);

        // get mesh enclosing rectangle
        float min_x = 1000000;
        float min_y = 1000000;
        float max_x = -1000000;
        float max_y = -1000000;
        for (uint32_t i = 0; i < data.numberofpoints; i++) {
            Point p = data.pointlist[i];
            if (p.x < min_x) min_x = p.x;
            if (p.y < min_y) min_y = p.y;
            if (p.x > max_x) max_x = p.x;
            if (p.y > max_y) max_y = p.y;
        }

        // object to window ratio
        float ratio = 0.9;

        // x scale factor to loosely fit mesh in window (equal in x and y)
        float x_scale = ratio * 2000 / std::max(max_x - min_x, max_y - min_y);
        // y scale factor to loosely fit mesh in window
        float y_scale = ratio * 2000 / std::max(max_x - min_x, max_y - min_y);
        // x offset to center mesh in window
        float x_offset = 0.5 * 2000 - 0.5 * (max_x + min_x) * x_scale;
        // y offset to center mesh in window
        float y_offset = 0.5 * 2000 - 0.5 * (max_y + min_y) * y_scale;

        // render

        while(true){
            // clears the window
            SDL_SetRenderDrawColor(rend, 0, 0, 0, 255);
            SDL_RenderClear(rend);
            // draw the points
            for (uint32_t i = 0; i < data.numberoftriangles; i++) {
                for (uint32_t j = 0; j < 3; j++) {
                    Point p1 = data.pointlist[data.trianglelist[i].verts[j]];
                    Point p2 = data.pointlist[data.trianglelist[i].verts[(j + 1) % 3]];
                    SDL_SetRenderDrawColor(rend, 255, 0, 0, 255);
                    SDL_RenderDrawLine(rend, x_scale * p1.x + x_offset, y_scale * p1.y + y_offset, x_scale * p2.x + x_offset, y_scale * p2.y + y_offset);

                    // draw the vertices
                    SDL_SetRenderDrawColor(rend, 255, 255, 255, 255);
                    SDL_RenderDrawPoint(rend, x_scale * p1.x + x_offset, y_scale * p1.y + y_offset);
                }
            }

            SDL_RenderPresent(rend);

            SDL_Event event;
            if (SDL_PollEvent(&event)) {
                if (event.type == SDL_QUIT) {
                    // destroy the window
                    SDL_DestroyWindow(win);
                    // clean up
                    SDL_Quit();
                    break;
                }
            }

            SDL_KeyboardEvent key = event.key;
            if (key.keysym.sym == SDLK_UP) {
                y_offset += 10;
            } else if (key.keysym.sym == SDLK_DOWN) {
                y_offset -= 10;
            } else if (key.keysym.sym == SDLK_LEFT) {
                x_offset += 10;
            } else if (key.keysym.sym == SDLK_RIGHT) {
                x_offset -= 10;
            } else if (key.keysym.sym == SDLK_w) {
                y_offset += 10;
            } else if (key.keysym.sym == SDLK_s) {
                y_offset -= 10;
            } else if (key.keysym.sym == SDLK_a) {
                y_offset += 10;
            } else if (key.keysym.sym == SDLK_d) {
                y_offset -= 10;
            } else if (key.keysym.sym == SDLK_q) {
                x_scale *= 0.9;
                y_scale *= 0.9;
            } else if (key.keysym.sym == SDLK_e) {
                x_scale *= 1.1;
                y_scale *= 1.1;
            }

            // sleep for 1 second
        }
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
            // printf("adding region %f %f %lu\n", r.first.x, r.first.y, r.second);
            in.regionlist[4 * i] = r.first.x;
            in.regionlist[4 * i + 1] = r.first.y;
            in.regionlist[4 * i + 2] = r.second;
            in.regionlist[4 * i + 3] = 0;
            i++;
        }
        // printf("----------------------\n");
        char switches[] = "pzq20AQa0.01";

        /* Make necessary initializations so that Triangle can return a */
        /*   triangulation in `mid' and a voronoi diagram in `vorout'.  */

        data.pointlist = (Point *) NULL;            /* Not needed if -N switch used. */
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
                Point v = data.pointlist[i];
                double mag_squared = v.x * v.x + v.y * v.y;
                double scale = R_squared / mag_squared;
                // printf("v = (%f, %f) -> (%f, %f), mag = %f, scale = %f\n", v->p.x, v->p.y, v->p.x * scale, v->p.y * scale, mag_squared, scale);
                double dist = geomDistance(v, center);
                if (dist < max_x) {
                // if (false) {
                    #ifdef DEBUG_PRINT
                    printf("kelvin transform too large\n");
                    #endif
                    data.pointlist[i] = v * ((radius * max_radius_coeff) / dist);
                } else {
                    data.pointlist[i] = v * scale;
                }
            }
        }
    }

    template <typename Prop>
    void Mesh<Prop>::addKelvinBoundaryConditions() {
        // kelvin mesh = km
        // this mesh = tm
        // kelvin mesh boundary vertices = kmb
        // this mesh boundary vertices = tmb
        // kelvin mesh vertices = kmv
        // this mesh vertices = tmv

        Mesh<Prop> kelvin_mesh;

        Circle boundary_circle = Circle(Point(0, 0), radius);
        kelvin_mesh.drawing.drawCircle(boundary_circle, BOUNDARY_VERTICES);
        // add region near the edge of the circle
        kelvin_mesh.drawing.drawRegion(Point(0, 0), default_prop);

        kelvin_mesh.mesh();
        kelvin_mesh.center = Point(0, 0);
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
        data.pointlist = (Point *) realloc(data.pointlist, new_point_number * sizeof(Point));
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
                if (kelvin_to_this.count(kelvin_mesh.data.trianglelist[i].verts[j])) {
                    data.trianglelist[data.numberoftriangles].verts[j] = kelvin_to_this[kelvin_mesh.data.trianglelist[i].verts[j]];
                }
            }
            data.triangleattributelist[data.numberoftriangles] = default_prop_id;
            data.numberoftriangles++;
        }

        data.numberoftriangles = new_triangle_number;
    }
}

#endif
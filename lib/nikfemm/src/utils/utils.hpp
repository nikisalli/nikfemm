#ifndef NIK_UTILS_HPP
#define NIK_UTILS_HPP

#include <string>   
#include <iostream>

#include "SDL2/SDL.h"

#include "../geometry/point.hpp"

namespace nikfemm {
    class PolygonShape {
        public:
            PolygonShape(std::vector<SDL_Point> vertices);
            ~PolygonShape();
            SDL_Point GetCenter(void);  
            SDL_Point * GetVertices(void);
            int GetNumberOfVertices(void);
        private:
            SDL_Point * vertices;
            SDL_Point center;
            int length;
    };

    void nexit(std::string message);
    void nassert(bool condition, std::string message);
    SDL_Color val2jet(double v, double vmin, double vmax);
    bool DrawFilledPolygon(SDL_Renderer* rend, PolygonShape poly, const SDL_Color color);
}

#endif
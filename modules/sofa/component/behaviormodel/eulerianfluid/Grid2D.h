/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_COMPONENT_BEHAVIORMODEL_EULERIANFLUID_GRID2D_H
#define SOFA_COMPONENT_BEHAVIORMODEL_EULERIANFLUID_GRID2D_H

#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/helper/rmath.h>
#include <sofa/component/component.h>
#include <iostream>


namespace sofa
{

namespace component
{

namespace behaviormodel
{

namespace eulerianfluid
{

using namespace sofa::defaulttype;
using namespace sofa::helper;

#ifndef NDEBUG
#define DEBUGGRID
#endif

class SOFA_EULERIAN_FLUID_API Grid2D
{
public:

    typedef float real;
    typedef Vec<2,real> vec2;

    struct Cell
    {
        vec2 u; ///< velocity (unit = cell width)
        int type; ///< First particle (or 0 for non-fluid cell)
        int dummy; ///< Align the structure to 16 bytes
        void clear()
        {
            u.clear();
            type = 0;
        }
    };

    int nx,ny,ncell;

    enum { PART_EMPTY=0, PART_WALL=-1, PART_FULL=1 };

    Cell* fdata;
    real* pressure;
    real* levelset;
    //char* distance;

    real t;
    real tend;

    real max_pressure;
    Cell bcell;

    static const unsigned long* obstacles;

    Grid2D();
    ~Grid2D();

    void clear(int _nx,int _ny);

    template<class T>
    static inline T lerp(T a, T b, real f) { return a+(b-a)*f; }

    int clamp_x(int x)
    {
        return (x<0)?0:(x>=nx)?nx-1:x;
    }

    int clamp_y(int y)
    {
        return (y<0)?0:(y>=ny)?ny-1:y;
    }

    int clamp_in_x(int x)
    {
        return (x<1)?1:(x>=nx-1)?nx-2:x;
    }

    int clamp_in_y(int y)
    {
        return (y<1)?1:(y>=ny-1)?ny-2:y;
    }

    int index(int x, int y) const
    {
        return x + y*nx;
    }

    int index(const vec2& p) const
    {
        return index(rnear(p[0]), rnear(p[1]));
    }

    Cell* get(int x, int y, Cell* base) const
    {
        return base+index(x,y);
    }

    const Cell* get(int x, int y, const Cell* base) const
    {
        return base+index(x,y);
    }

    Cell* get(int x, int y)
    {
#ifdef DEBUGGRID
        if (((unsigned)x>=(unsigned)nx) || ((unsigned)y>=(unsigned)ny))
        {
            std::cerr<<"INVALID CELL "<<x<<','<<y<<std::endl;
            return &bcell;
        }
#endif
        return get(x,y,fdata);
    }

    const Cell* get(int x, int y) const
    {
#ifdef DEBUGGRID
        if (((unsigned)x>=(unsigned)nx) || ((unsigned)y>=(unsigned)ny))
        {
            std::cerr<<"INVALID CELL "<<x<<','<<y<<std::endl;
            return &bcell;
        }
#endif
        return get(x,y,fdata);
    }

    Cell* get(const vec2& p)
    {
        return get(rnear(p[0]),rnear(p[1]));
    }

    const Cell* get(const vec2& p) const
    {
        return get(rnear(p[0]),rnear(p[1]));
    }

    template<int C> real interp(const Cell* base, real fx, real fy) const
    {
        return lerp( lerp(get(0,0,base)->u[C],get(1,0,base)->u[C],fx),
                lerp(get(0,1,base)->u[C],get(1,1,base)->u[C],fx), fy );
    }

    template<int C> real interp(vec2 p) const
    {
        p[C] += 0.5;
        int ix = rfloor(p[0]);
        int iy = rfloor(p[1]);
        const Cell* base = get(ix,iy);
        return interp<C>(base, p[0]-ix, p[1]-iy);
    }

    vec2 interp(vec2 p) const
    {
        return vec2( interp<0>(p), interp<1>(p) );
    }

    template<int C> void impulse(Cell* base, real fx, real fy, real i)
    {
        get(0,0,base)->u[C] += i*(1-fx)*(1-fy);
        get(1,0,base)->u[C] += i*(  fx)*(1-fy);
        get(0,1,base)->u[C] += i*(1-fx)*(  fy);
        get(1,1,base)->u[C] += i*(  fx)*(  fy);
    }

    template<int C> void impulse(vec2 p, real i)
    {
        p[C] += 0.5;
        int ix = rfloor(p[0]);
        int iy = rfloor(p[1]);
        Cell* base = get(ix,iy);
        impulse<C>(base, p[0]-ix, p[1]-iy, i);
    }

    void impulse(const vec2& p, const vec2& i)
    {
        impulse<0>(p,i[0]);
        impulse<1>(p,i[1]);
    }

    real* getpressure(int x, int y)
    {
        return pressure + index(x,y);
    }

    real getpressure(vec2 p)
    {
        int ix = rfloor(p[0]);
        int iy = rfloor(p[1]);
        real fx = p[0]-ix;
        real fy = p[1]-iy;
        real* base = getpressure(ix,iy);
        return lerp( lerp(base[index(0,0)],base[index(1,0)],fx),
                lerp(base[index(0,1)],base[index(1,1)],fx), fy );
    }

    real* getlevelset(int x, int y)
    {
        return levelset + index(x,y);
    }

    real getlevelset(vec2 p)
    {
        int ix = rfloor(p[0]);
        int iy = rfloor(p[1]);
        real fx = p[0]-ix;
        real fy = p[1]-iy;
        real* base = getlevelset(ix,iy);
        return lerp( lerp(base[index(0,0)],base[index(1,0)],fx),
                lerp(base[index(0,1)],base[index(1,1)],fx), fy );
    }

    void seed(real height);

    void seed(real height, vec2 normal);

    void seed(vec2 p0, vec2 p1, vec2 velocity=vec2(0,0));

    void step(Grid2D* prev, Grid2D* temp, real dt=0.04, real diff=0.00001);
    void step_init(const Grid2D* prev, Grid2D* temp, real dt, real diff);
    //void step_particles(Grid2D* prev, Grid2D* temp, real dt, real diff);
    void step_levelset(Grid2D* prev, Grid2D* temp, real dt, real diff);
    void step_forces(const Grid2D* prev, Grid2D* temp, real dt, real diff, real scale=1.0);
    void step_surface(const Grid2D* prev, Grid2D* temp, real dt, real diff);
    void step_advect(const Grid2D* prev, Grid2D* temp, real dt, real diff);
    void step_diffuse(const Grid2D* prev, Grid2D* temp, real dt, real diff);
    void step_project(const Grid2D* prev, Grid2D* temp, real dt, real diff);
    void step_color(const Grid2D* prev, Grid2D* temp, real dt, real diff);

    // internal helper function
    //  template<int C> inline real find_velocity(int x, int y, int ind, int ind2, const Grid2D* prev, const Grid2D* temp);

    // Fast Marching Method Levelset Update
    enum Status { FMM_FRONT0 = 0, FMM_FAR = -1, FMM_KNOWN = -2, FMM_BORDER = -3 };
    int* fmm_status;
    int* fmm_heap;
    int fmm_heap_size;

    int fmm_pop();
    void fmm_push(int index);
    void fmm_swap(int entry1, int entry2);
};

} // namespace eulerianfluid

} // namespace behaviormodel

} // namespace component

} // namespace sofa

#endif

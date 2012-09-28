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
#ifndef SOFA_GPU_CUDA_CUDAMATRIX_H
#define SOFA_GPU_CUDA_CUDAMATRIX_H

//#include "host_runtime.h" // CUDA
#include "CudaTypes.h"
#include <iostream>

//#define DEBUG_OUT_MATRIX

#ifdef DEBUG_OUT_MATRIX
#define DEBUG_OUT_M(a) a
#define SPACEP std::cout << "(" << hostIsValid << "," << (deviceIsValid) << ") " ;for(int espaceaff=0;espaceaff<spaceDebug;espaceaff++) std::cout << "  ";spaceDebug++; std::cout << ">"
#define SPACEM std::cout << "(" << hostIsValid << "," << (deviceIsValid) << ") " ;spaceDebug--;for(int espaceaff=0;espaceaff<spaceDebug;espaceaff++) std::cout << "  "; std::cout << "<"
#define SPACEN std::cout << "(" << hostIsValid << "," << (deviceIsValid) << ") " ;for(int espaceaff=0;espaceaff<spaceDebug;espaceaff++) std::cout << "  "; std::cout << "|"
#else
#define DEBUG_OUT_M(a)
#endif


namespace sofa
{

namespace gpu
{

namespace cuda
{

template<class T, class MemoryManager = CudaMemoryManager<T> >
class CudaMatrix
{
public:
    typedef CudaMatrix<T> Matrix;
    typedef T      value_type;
    typedef size_t size_type;

private:
    size_type    sizeX;     ///< Current size of the vector
    size_type    sizeY;     ///< Current size of the vector
    size_type    pitch_device;     ///< Row alignment on the GPU
    size_type    pitch_host;     ///< Row alignment on the GPU
    size_type    allocSizeY;  ///< Allocated size
    void*        devicePointer;  ///< Pointer to the data on the GPU side
    T*           hostPointer;    ///< Pointer to the data on the CPU side
    mutable bool deviceIsValid;  ///< True if the data on the GPU is currently valid
    mutable bool hostIsValid;    ///< True if the data on the CPU is currently valid
    DEBUG_OUT_M(mutable int spaceDebug;)
public:

    CudaMatrix()
        : sizeX ( 0 ), sizeY( 0 ), pitch_device(0), pitch_host ( 0 ), allocSizeY ( 0 ), devicePointer ( NULL ), hostPointer ( NULL ), deviceIsValid ( true ), hostIsValid ( true )
    {
        DEBUG_OUT_M(spaceDebug = 0);
    }

    CudaMatrix(size_t x, size_t y, size_t size)
        : sizeX ( 0 ), sizeY ( 0 ), pitch_device(0), pitch_host ( 0 ), allocSizeY ( 0 ), devicePointer ( NULL ), hostPointer ( NULL ), deviceIsValid ( true ), hostIsValid ( true )
    {
        resize (x,y,size);
        DEBUG_OUT_M(spaceDebug = 0);
    }

    CudaMatrix(const CudaMatrix<T>& v )
        : sizeX ( 0 ), sizeY ( 0 ), pitch_device(0), pitch_host ( 0 ), allocSizeY ( 0 ), devicePointer ( NULL ), hostPointer ( NULL ), deviceIsValid ( true ), hostIsValid ( true )
    {
        *this = v;
        DEBUG_OUT_M(spaceDebug = 0);
    }

    void clear()
    {
        DEBUG_OUT_M(SPACEP << "Clear" << std::endl);
        sizeX = 0;
        sizeY = 0;
        deviceIsValid = true;
        hostIsValid = true;
        DEBUG_OUT_M(SPACEM << "Clear" << std::endl);
    }

    ~CudaMatrix()
    {
        if (hostPointer!=NULL) mycudaFreeHost(hostPointer);
        if (devicePointer!=NULL) mycudaFree(devicePointer);
    }

    size_type getSizeX() const
    {
        return sizeX;
    }

    size_type getSizeY() const
    {
        return sizeY;
    }

    size_type getPitchDevice() const
    {
        return pitch_device;
    }

    size_type getPitchHost() const
    {
        return pitch_host;
    }

    bool empty() const
    {
        return sizeX==0 || sizeY==0;
    }

    void memsetHost(int v = 0)
    {
        DEBUG_OUT_M(SPACEP << "memsetHost" << std::endl);
        MemoryManager::memsetHost(hostPointer,v,pitch_host*sizeY);
        hostIsValid = true;
        deviceIsValid = false;
        DEBUG_OUT_M(SPACEM << "memsetHost" << std::endl);
    }

    void memsetDevice(int v = 0)
    {
        DEBUG_OUT_M(SPACEP << "memsetHost" << std::endl);
        MemoryManager::memsetDevice(0,devicePointer, v, pitch_device*sizeY);
        hostIsValid = false;
        deviceIsValid = true;
        DEBUG_OUT_M(SPACEM << "memsetHost" << std::endl);
    }

    void invalidateDevices()
    {
        hostIsValid = true;
        deviceIsValid = false;
    }

    void invalidatehost()
    {
        hostIsValid = false;
        deviceIsValid = true;
    }

    void fastResize(size_type y,size_type x,size_type WARP_SIZE=MemoryManager::BSIZE)
    {
        DEBUG_OUT_M(SPACEP << "fastResize : " << x << " " << y << " WArp_Size=" << WARP_SIZE << " sizeof(T)=" << sizeof(T) << std::endl);

        if ( x==0 || y==0)
        {
            clear();
            DEBUG_OUT_M(SPACEM << std::endl);
            return;
        }
        if ( sizeX==x && sizeY==y)
        {
            DEBUG_OUT_M(SPACEM << std::endl);
            return;
        }

        size_type d_x = x;
        size_type d_y = y;

        if (WARP_SIZE==0)
        {
            d_x = x;
            d_y = y;
        }
        else
        {
            d_x = ((d_x+WARP_SIZE-1)/WARP_SIZE)*WARP_SIZE;
            d_y = ((d_y+WARP_SIZE-1)/WARP_SIZE)*WARP_SIZE;
        }
        size_type allocSize = d_x*d_y*sizeof(T);

        if ( !sizeX && !sizeY)   //special case anly reserve
        {
            DEBUG_OUT_M(SPACEN << "Is in ( !sizeX && !sizeY)" << std::endl);
            if (allocSize > pitch_host*allocSizeY || pitch_host < d_x*sizeof(T))
            {
                T* prevHostPointer = hostPointer;
                MemoryManager::hostAlloc( (void **) &hostPointer, allocSize ); pitch_host = d_x*sizeof(T);
                DEBUG_OUT_M(SPACEN << "Allocate Host : " << ((int) hostPointer) << " HostPitch = " << pitch_host << std::endl);
                if ( prevHostPointer != NULL ) MemoryManager::hostFree( prevHostPointer );

                void* prevDevicePointer = devicePointer;

                mycudaMallocPitch(&devicePointer, &pitch_device, d_x*sizeof(T), d_y);
                DEBUG_OUT_M(SPACEN << "Allocate Device : " << ((int) devicePointer) << " DevicePitch = " << pitch_device << std::endl);
                if ( prevDevicePointer != NULL ) mycudaFree ( prevDevicePointer );

                allocSizeY = d_y;
            }
        }
        else if (x*sizeof(T) <= pitch_host)
        {
            DEBUG_OUT_M(SPACEN << "Is in (x <= pitch_host)" << std::endl);
            if (d_y > allocSizeY)   // allocate
            {
                DEBUG_OUT_M(SPACEN << "Is in (y > allocSizeY)" << std::endl);
                T* prevHostPointer = hostPointer;
                MemoryManager::hostAlloc( (void **) &hostPointer, allocSize);
                DEBUG_OUT_M(SPACEN << "Allocate Host : " << ((int) hostPointer) << " HostPitch = " << pitch_host << std::endl);
                if (hostIsValid)
                {
                    DEBUG_OUT_M(SPACEN << "MemcpyHost from 0 to " << (pitch_host*sizeY) << std::endl);
                    std::copy ( prevHostPointer, ((T*) ((char*)prevHostPointer)+(pitch_host*sizeY)), hostPointer);
                }
                if ( prevHostPointer != NULL ) MemoryManager::hostFree( prevHostPointer );

                void* prevDevicePointer = devicePointer;
                mycudaMallocPitch(&devicePointer, &pitch_device, pitch_device, d_y); //pitch_device should not be modified
                DEBUG_OUT_M(SPACEN << "Allocate Device : " << ((int) devicePointer) << " DevicePitch = " << pitch_device << std::endl);
                if (deviceIsValid)
                {
                    DEBUG_OUT_M(SPACEN << "MemcpyDevice from 0 to " << (pitch_device*sizeY) << std::endl);
                    MemoryManager::memcpyDeviceToDevice (0, devicePointer, prevDevicePointer, pitch_device*sizeY );
                }
                if ( prevDevicePointer != NULL ) mycudaFree ( prevDevicePointer );

                allocSizeY = d_y;
            }
        }
        else     //necessary to change layout....
        {
            std::cerr << "ERROR : resize column is not implemented in CudaMatrix, you must copy data in another matrix" << std::endl;
            DEBUG_OUT_M(SPACEM << std::endl);
            return;
        }

        sizeX = x; sizeY = y;

        DEBUG_OUT_M(SPACEM << "fastResize" << std::endl);
    }

    void recreate(int nbRow,int nbCol)
    {
        clear();
        fastResize(nbRow,nbCol);
    }


    void resize (size_type y,size_type x,size_t WARP_SIZE=MemoryManager::BSIZE)
    {
        DEBUG_OUT_M(SPACEP << "reisze : " << x << " " << y << " WArp_Size=" << WARP_SIZE << " sizeof(T)=" << sizeof(T) << std::endl);

        if ((x==0) || (y==0))
        {
            clear();
            DEBUG_OUT_M(SPACEM << std::endl);
            return;
        }
        if ((sizeX==x) && (sizeY==y))
        {
            DEBUG_OUT_M(SPACEM << std::endl);
            return;
        }

        size_type d_x = x;
        size_type d_y = y;

        if (WARP_SIZE==0)
        {
            d_x = x;
            d_y = y;
        }
        else
        {
            d_x = ((d_x+WARP_SIZE-1)/WARP_SIZE)*WARP_SIZE;
            d_y = ((d_y+WARP_SIZE-1)/WARP_SIZE)*WARP_SIZE;
        }
        size_type allocSize = d_x*d_y*sizeof(T);

        if ( !sizeX && !sizeY)   //special case anly reserve
        {
            DEBUG_OUT_M(SPACEN << "Is in ( !sizeX && !sizeY)" << std::endl);
            if (allocSize > pitch_host*allocSizeY || pitch_host < d_x*sizeof(T))
            {
                T* prevHostPointer = hostPointer;
                MemoryManager::hostAlloc( (void **) &hostPointer, allocSize ); pitch_host = d_x*sizeof(T);
                DEBUG_OUT_M(SPACEN << "Allocate Host : " << ((unsigned long) hostPointer) << " HostPitch = " << pitch_host << std::endl);
                if ( prevHostPointer != NULL ) MemoryManager::hostFree( prevHostPointer );

                void* prevDevicePointer = devicePointer;
                mycudaMallocPitch(&devicePointer, &pitch_device, d_x*sizeof(T), d_y);
                DEBUG_OUT_M(SPACEN << "Allocate Device : " << ((unsigned long) devicePointer) << " DevicePitch = " << pitch_device << std::endl);
                if ( prevDevicePointer != NULL ) mycudaFree ( prevDevicePointer );

                allocSizeY = d_y;
            }
            else if (pitch_host < d_x*sizeof(T)) pitch_host = d_x*sizeof(T);

            if (hostIsValid)
            {
                DEBUG_OUT_M(SPACEN << "MemsetHost from 0 to " << (pitch_host*y) << std::endl);
                MemoryManager::memsetHost(hostPointer,0,pitch_host*y);
            }
            if (deviceIsValid)
            {
                DEBUG_OUT_M(SPACEN << "MemsetDevice from 0 to " << (pitch_device*y) << std::endl);
                MemoryManager::memsetDevice(0,devicePointer, 0, pitch_device*y);
            }
        }
        else if (x*sizeof(T) <= pitch_host)
        {
            DEBUG_OUT_M(SPACEN << "Is in (x <= pitch_host)" << std::endl);
            if (d_y > allocSizeY)   // allocate
            {
                DEBUG_OUT_M(SPACEN << "Is in (y > allocSizeY)" << std::endl);
                T* prevHostPointer = hostPointer;
                MemoryManager::hostAlloc( (void **) &hostPointer, allocSize);
                DEBUG_OUT_M(SPACEN << "Allocate Host : " << ((unsigned long) hostPointer) << " HostPitch = " << pitch_host << std::endl);
                if (hostIsValid)
                {
                    DEBUG_OUT_M(SPACEN << "MemcpyHost from 0 to " << (pitch_host*sizeY) << std::endl);
                    std::copy ( prevHostPointer, ((T*) ((char*)prevHostPointer)+(pitch_host*sizeY)), hostPointer);
                }
                if ( prevHostPointer != NULL ) MemoryManager::hostFree( prevHostPointer );

                void* prevDevicePointer = devicePointer;
                mycudaMallocPitch(&devicePointer, &pitch_device, pitch_device, d_y); //pitch_device should not be modified
                DEBUG_OUT_M(SPACEN << "Allocate Device : " << ((unsigned long) devicePointer) << " DevicePitch = " << pitch_device << std::endl);
                if (deviceIsValid)
                {
                    DEBUG_OUT_M(SPACEN << "MemcpyDevice from 0 to " << (pitch_device*sizeY) << std::endl);
                    MemoryManager::memcpyDeviceToDevice (0, devicePointer, prevDevicePointer, pitch_device*sizeY );
                }
                if ( prevDevicePointer != NULL ) mycudaFree ( prevDevicePointer );

                allocSizeY = d_y;
            }
            if (y > sizeY)
            {
                DEBUG_OUT_M(SPACEN << "Is in (y > sizeY)" << std::endl);
                if (hostIsValid)
                {
                    DEBUG_OUT_M(SPACEN << "MemsetHost from " << pitch_host*sizeY << " to " << (pitch_host*y) << std::endl);
                    MemoryManager::memsetHost((T*) (((char*)hostPointer)+pitch_host*sizeY),0,pitch_host*(y-sizeY));
                }

                if (deviceIsValid)
                {
                    DEBUG_OUT_M(SPACEN << "MemsetDevice from " << pitch_device*sizeY << " to " << (pitch_device*y) << std::endl);
                    MemoryManager::memsetDevice(0,((char*)devicePointer) + (pitch_device*sizeY), 0, pitch_device*(y-sizeY));
                }
            }
            if (x > sizeX)
            {
                DEBUG_OUT_M(SPACEN << "Is in (x > sizeX)" << std::endl);
                if (hostIsValid)
                {
                    DEBUG_OUT_M(SPACEN << "MemsetHost for each line from " << (sizeX*sizeof(T)) << " to " << (sizeof(T)*(x-sizeX)) << std::endl);
                    for (unsigned j=0; j<y; j++) MemoryManager::memsetHost((T*) (((char*)hostPointer)+pitch_host*j+sizeX*sizeof(T)),0,sizeof(T)*(x-sizeX));
                }

                if (deviceIsValid)
                {
                    DEBUG_OUT_M(SPACEN << "MemsetDevice for each line from " << (sizeX*sizeof(T)) << " to " << (sizeof(T)*(x-sizeX)) << std::endl);
                    for (unsigned j=0; j<y; j++) MemoryManager::memsetDevice(0,((char*)devicePointer)+(pitch_device*j+sizeX*sizeof(T)), 0, sizeof(T)*(x-sizeX));
                }
            }
        }
        else     //necessary to change layout....
        {
            std::cerr << "ERROR : resize column is not implemented in CudaMatrix, you must copy data in another matrix" << std::endl;
            DEBUG_OUT_M(SPACEM);
            return;
        }

        sizeX = x; sizeY = y;

        DEBUG_OUT_M(SPACEM << "reisze" << std::endl);
    }

    void swap ( CudaMatrix<T>& v )
    {
#define VSWAP(type, var) { type t = var; var = v.var; v.var = t; }
        VSWAP ( size_type, sizeX );
        VSWAP ( size_type, sizeY );
        VSWAP ( size_type, pitch_device );
        VSWAP ( size_type, pitch_host );
        VSWAP ( size_type, allocSizeY );
        VSWAP ( void*    , devicePointer );
        VSWAP ( T*       , hostPointer );
        VSWAP ( bool     , deviceIsValid );
        VSWAP ( bool     , hostIsValid );
#undef VSWAP
    }

    const void* deviceRead ( int y=0, int x=0 ) const
    {
        copyToDevice();
        return ((T*) ((char*)devicePointer) + pitch_device*y) + x;

    }

    void* deviceWrite ( int y=0, int x=0 )
    {
        copyToDevice();
        hostIsValid = false;
        return ((T*) (((char*)devicePointer) + pitch_device*y)) + x;
    }

    const T* hostRead ( int y=0, int x=0 ) const
    {
        copyToHost();
        return ((const T*) (((char*) hostPointer)+(y*pitch_host))) + x;
    }

    T* hostWrite ( int y=0, int x=0 )
    {
        copyToHost();
        deviceIsValid = false;
        return ((T*) (((char*) hostPointer)+(y*pitch_host))) + x;
    }

    bool isHostValid() const
    {
        return hostIsValid;
    }

    bool isDeviceValid() const
    {
        return deviceIsValid;
    }

    const T& operator() (size_type y,size_type x) const
    {
        checkIndex (y,x);
        return hostRead(y,x);
    }

    T& operator() (size_type y,size_type x)
    {
        checkIndex (y,x);
        return hostWrite(y,x);
    }

    const T* operator[] (size_type y) const
    {
        checkIndex (y,0);
        return hostRead(y,0);
    }

    T* operator[] (size_type y)
    {
        checkIndex (y,0);
        return hostWrite(y,0);
    }

    const T& getCached (size_type y,size_type x) const
    {
        checkIndex (y,x);
        return ((T*) (((char*) hostPointer)+(y*pitch_host))) + x;
    }

    void operator= ( const CudaMatrix<T,MemoryManager >& m )
    {
        if (&m == this) return;

        sizeX = m.sizeX;
        sizeY = m.sizeY;

        if (sizeY*pitch_host<m.sizeY*m.pitch_host)   //simple case, we simply copy data with the same attribute
        {
            T* prevHostPointer = hostPointer;
            MemoryManager::hostAlloc( (void **) &hostPointer, m.pitch_host * sizeY);
            if ( prevHostPointer != NULL ) MemoryManager::hostFree( prevHostPointer );

            void* prevDevicePointer = devicePointer;
            mycudaMallocPitch(&devicePointer, &pitch_device, m.pitch_device, sizeY);
            if ( prevDevicePointer != NULL ) mycudaFree ( prevDevicePointer );

            allocSizeY = sizeY;
        }
        else
        {
            int allocline = (allocSizeY*pitch_host) / m.pitch_host;
            allocSizeY = allocline * m.pitch_host;
            // Here it's possible that the allocSizeY is < the the real memory allocated, but it's not a problem, it will we deleted at the next resize;
        }

        pitch_host = m.pitch_host;
        pitch_device = m.pitch_device;

        if (m.hostIsValid) std::copy ( m.hostPointer, ((T*) (((char*) m.hostPointer)+(m.pitch_host*m.sizeY))), hostPointer);
        if (m.deviceIsValid) MemoryManager::memcpyDeviceToDevice (0, devicePointer, m.devicePointer, m.pitch_device*m.sizeY );

        hostIsValid = m.hostIsValid;
        deviceIsValid = m.deviceIsValid; /// finally we get the correct device valid
    }

    friend std::ostream& operator<< ( std::ostream& os, const Matrix & mat )
    {
        mat.hostRead();
        os << "[\n";
        for (unsigned j=0; j<mat.getSizeY(); j++)
        {
            os << "[ ";
            for (unsigned i=0; i<mat.getSizeX(); i++)
            {
                os << " " << mat[j][i];
            }
            os << "]\n";
        }
        os << "]\n";
        return os;
    }

protected:
    void copyToHost() const
    {
        if ( hostIsValid ) return;
        DEBUG_OUT_M(SPACEN << "copyToHost" << std::endl);

//#ifndef NDEBUG
        if (mycudaVerboseLevel>=LOG_TRACE) std::cout << "CUDA: GPU->CPU copy of "<<sofa::core::objectmodel::BaseClass::decodeTypeName ( typeid ( *this ) ) <<": "<<sizeX*sizeof(T) <<" B"<<std::endl;
//#endif
        DEBUG_OUT_M(SPACEN << "copyToHost host : " << ((unsigned long) hostPointer) << " pitchH : " << pitch_host << " | device : " << ((unsigned long)devicePointer) << " pitchD : " << pitch_device << " | (" << sizeX*sizeof(T) << "," << sizeY << ")" << std::endl);
        mycudaMemcpyDeviceToHost2D ( hostPointer, pitch_host, devicePointer, pitch_device, sizeX*sizeof(T), sizeY);
        hostIsValid = true;
    }
    void copyToDevice() const
    {
        if ( deviceIsValid ) return;

//#ifndef NDEBUG
        if (mycudaVerboseLevel>=LOG_TRACE) std::cout << "CUDA: CPU->GPU copy of "<<sofa::core::objectmodel::BaseClass::decodeTypeName ( typeid ( *this ) ) <<": "<<sizeX*sizeof(T) <<" B"<<std::endl;
//#endif
        DEBUG_OUT_M(SPACEN << "copyToDevice device : " << ((unsigned long)devicePointer) << " pitchD : " << pitch_device << " | host : " << ((unsigned long) hostPointer) << " pitchH : " << pitch_host << " | (" << sizeX*sizeof(T) << "," << sizeY << ")" << std::endl);
        mycudaMemcpyHostToDevice2D ( devicePointer, pitch_device, hostPointer, pitch_host,  sizeX*sizeof(T), sizeY);
        deviceIsValid = true;
    }

#ifdef NDEBUG
    void checkIndex ( size_type,size_type ) const
    {
    }
#else
    void checkIndex ( size_type y,size_type x) const
    {
        assert (x<this->sizeX);
        assert (y<this->sizeY);
    }
#endif

};

#ifdef DEBUG_OUT_MATRIX
#undef DEBUG_OUT_M
#undef SPACEP
#undef SPACEM
#undef SPACEN
#undef DEBUG_OUT_MATRIX
#else
#undef DEBUG_OUT_M
#endif

} // namespace cuda

} // namespace gpu

} // namespace sofa

#endif

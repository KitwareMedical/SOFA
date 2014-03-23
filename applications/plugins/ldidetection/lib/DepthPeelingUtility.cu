//#include <sofa/gpu/cuda/CudaMath.h>
#include <sofa/gpu/cuda/CudaMath.h>

namespace sofa
{

namespace component
{

namespace collision
{
	extern "C"
	{
		void CollisionVolume_count(const void* layers, void* counts, void* count_per_bloc, const int nlayers, const int npixels, const int bsize, bool self);
		void CollisionVolume_write(const void* layers, const void* counts, const void* bloc_input_pos, const void* bloc_output_pos, void* collisions, const int nlayers, const int npixels, const int bsize, const int nblocs, bool self);
	}

__global__ void CollisionVolume_count_kernel(const float4* layers, unsigned int* counts, unsigned int* count_sums, const int nlayers, const int npixels) //, int offset)
{
    //! Dynamically allocated shared memory for gather
    extern  __shared__ unsigned int temp[];
    extern  __shared__ float ftemp[];
    const int index0 = umul24(blockIdx.x,blockDim.x);
    const int index1 = threadIdx.x;
    const int index2 = index1+index1;
    int index = index0+index1;
    //float vol = 0.0f;
    unsigned int count = 0;
    //layers += index;
    layers += index0;
    counts += index;
    //float4 first = *layers;
/*
    int p;
    p = index1;
    ftemp[p] = ((const float*)layers)[p]; p+=blockDim.x;
    ftemp[p] = ((const float*)layers)[p]; p+=blockDim.x;
    ftemp[p] = ((const float*)layers)[p]; p+=blockDim.x;
    ftemp[p] = ((const float*)layers)[p]; //p+=blockDim.x;
    __syncthreads();
    float4 first=make_float4(ftemp[index4+0],ftemp[index4+1],ftemp[index4+2],ftemp[index4+3]);
    __syncthreads();
*/
    float2 f2 = ((const float2*)layers)[index1];
    ftemp[index1] = (index1&1)?f2.x:f2.y;
    f2 = ((const float2*)layers)[index1+blockDim.x];
    ftemp[index1+blockDim.x] = (index1&1)?f2.x:f2.y;
    __syncthreads();
    float2 first=make_float2(ftemp[index2+0],ftemp[index2+1]);
    __syncthreads();
    for(int L=nlayers-2;L>=0;--L)
    {
        layers += npixels;
        //float4 current = *layers;
/*
        p = index1;
        ftemp[p] = ((const float*)layers)[p]; p+=blockDim.x;
        ftemp[p] = ((const float*)layers)[p]; p+=blockDim.x;
        ftemp[p] = ((const float*)layers)[p]; p+=blockDim.x;
        ftemp[p] = ((const float*)layers)[p]; //p+=blockDim.x;
        __syncthreads();
        float4 current=make_float4(ftemp[index4+0],ftemp[index4+1],ftemp[index4+2],ftemp[index4+3]);
        __syncthreads();
*/
        f2 = ((const float2*)layers)[index1];
        ftemp[index1] = (index1&1)?f2.x:f2.y;
        f2 = ((const float2*)layers)[index1+blockDim.x];
        ftemp[index1+blockDim.x] = (index1&1)?f2.x:f2.y;
        __syncthreads();
        float2 current=make_float2(ftemp[index2+0],ftemp[index2+1]);
        __syncthreads();
        //if (current.z!=0.0 && first.z > 0.5 && current.z < 0.5 && floor(first.y*2) != floor(current.y*2))
        if (current.y!=0.0 && first.y > 0.5 && current.y < 0.5 && floor(first.x*2) != floor(current.x*2))
        {
            ++count;
            //vol += current.x - first.x;
        }
        first = current;
    }

    //*counts = count;
    int n = blockDim.x;

    unsigned int v = count;
    for (int i=1;i<n;)
    {
        temp[index1] = v;
        __syncthreads();
        if (index1>=i)
            v += temp[index1 - i];
        i<<=1;
        temp[n + index1] = v;
        __syncthreads();
        if (index1>=i)
            v += temp[n + index1 - i];
        i<<=1;
    }
    *counts = v - count;
    if (index1 == n-1)
        count_sums[blockIdx.x] = v;

    /*
    temp[index1+index1] = (float)count;
    temp[index1+index1+1] = vol;
    __syncthreads();
    vol = temp[index1]+temp[index1+n];
    __syncthreads();
    */
    /*
    while(offset>0)
    {
        if (index1 >= offset && index1 < n)
            temp[index1] = count;
        __syncthreads();
        if (index1+offset < n)
            count += temp[index1+offset];
        n = offset;
        offset >>= 1;
    }
    if (index1 == 0)
        count_sums[blockIdx.x] = count;
    */
}

__global__ void SelfCollisionVolume_count_kernel(const float4* layers, unsigned int* counts, unsigned int* count_sums, const int nlayers, const int npixels) //, int offset)
{
    //! Dynamically allocated shared memory for gather
    extern  __shared__  unsigned int temp[];
    extern  __shared__  float ftemp[];
    const int index0 = umul24(blockIdx.x,blockDim.x);
    const int index1 = threadIdx.x;
    const int index2 = index1+index1;
    int index = index0+index1;
    //float vol = 0.0f;
    unsigned int count = 0;
    //layers += index;
    layers += index0;
    counts += index;
    //float4 first = *layers;
/*
    int p;
    p = index1;
    ftemp[p] = ((const float*)layers)[p]; p+=blockDim.x;
    ftemp[p] = ((const float*)layers)[p]; p+=blockDim.x;
    ftemp[p] = ((const float*)layers)[p]; p+=blockDim.x;
    ftemp[p] = ((const float*)layers)[p]; //p+=blockDim.x;
    __syncthreads();
    float4 first=make_float4(ftemp[index4+0],ftemp[index4+1],ftemp[index4+2],ftemp[index4+3]);
    __syncthreads();
*/
    float2 f2 = ((const float2*)layers)[index1];
    ftemp[index1] = (index1&1)?f2.x:f2.y;
    f2 = ((const float2*)layers)[index1+blockDim.x];
    ftemp[index1+blockDim.x] = (index1&1)?f2.x:f2.y;
    __syncthreads();
    float first=ftemp[index2+1];
    __syncthreads();

    int inside = floor(2*first)*2-1;
    for(int L=nlayers-2;L>=0;--L)
    {
        layers += npixels;
        //float4 current = *layers;
/*
        p = index1;
        ftemp[p] = ((const float*)layers)[p]; p+=blockDim.x;
        ftemp[p] = ((const float*)layers)[p]; p+=blockDim.x;
        ftemp[p] = ((const float*)layers)[p]; p+=blockDim.x;
        ftemp[p] = ((const float*)layers)[p]; //p+=blockDim.x;
        __syncthreads();
        float4 current=make_float4(ftemp[index4+0],ftemp[index4+1],ftemp[index4+2],ftemp[index4+3]);
        __syncthreads();
*/
        f2 = ((const float2*)layers)[index1];
        ftemp[index1] = (index1&1)?f2.x:f2.y;
        f2 = ((const float2*)layers)[index1+blockDim.x];
        ftemp[index1+blockDim.x] = (index1&1)?f2.x:f2.y;
        __syncthreads();
        float current=ftemp[index2+1];
        __syncthreads();
        if (current > 0.5)
            ++inside;
        else
        {
            --inside;
            if (current!=0.0 && first > 0.5 && inside > 0)
            {
                ++count;
            }
        }
        first = current;
    }

    //*counts = count;
    int n = blockDim.x;

    unsigned int v = count;
    for (int i=1;i<n;)
    {
        temp[index1] = v;
        __syncthreads();
        if (index1>=i)
            v += temp[index1 - i];
        i<<=1;
        temp[n + index1] = v;
        __syncthreads();
        if (index1>=i)
            v += temp[n + index1 - i];
        i<<=1;
    }
    *counts = v - count;
    if (index1 == n-1)
        count_sums[blockIdx.x] = v;
    /*
    temp[index1+index1] = (float)count;
    temp[index1+index1+1] = vol;
    __syncthreads();
    vol = temp[index1]+temp[index1+n];
    __syncthreads();
    */
    /*
    while(offset>0)
    {
        if (index1 >= offset && index1 < n)
            temp[index1] = count;
        __syncthreads();
        if (index1+offset < n)
            count += temp[index1+offset];
        n = offset;
        offset >>= 1;
    }
    if (index1 == 0)
        count_sums[blockIdx.x] = count;
    */
}

void CollisionVolume_count(const void* layers, void* counts, void* count_sums, const int nlayers, const int npixels, const int bsize, bool self)
{
    int nblocs = (npixels+bsize-1)/bsize;
    dim3 threads(bsize,1);
    dim3 grid(nblocs,1);
    /*int offset;
    if (bsize==1)
        offset = 0;
    else
    {
        offset = 1;
        while (offset*2 < bsize)
            offset *= 2;
    }*/
    if (self)
        SelfCollisionVolume_count_kernel<<< grid, threads, 2*bsize * sizeof(unsigned int) >>>((const float4*)layers, (unsigned int*)counts, (unsigned int*)count_sums, nlayers, npixels); //, offset);
    else
        CollisionVolume_count_kernel<<< grid, threads, 2*bsize * sizeof(unsigned int) >>>((const float4*)layers, (unsigned int*)counts, (unsigned int*)count_sums, nlayers, npixels); //, offset);
}

__global__ void CollisionVolume_write_kernel(const float4* layers, const unsigned int* counts, const unsigned int* bloc_input_pos, const unsigned int* bloc_output_pos, float4* collisions, const int nlayers, const int npixels)
{
    //! Dynamically allocated shared memory
    extern  __shared__  unsigned int temp[];
    extern  __shared__  float ftemp[];
    const int index1 = threadIdx.x;
    const int index2 = index1+index1;
    if (index1 == 0)
    {
        temp[0] = bloc_input_pos[blockIdx.x];
        temp[1] = bloc_output_pos[blockIdx.x];
    }
    __syncthreads();
    layers += temp[0]; //+index1;
    collisions += (temp[1] + counts[temp[0]+index1])<<1;
    //float4 first = *layers;
    __syncthreads();
    float2 f2 = ((const float2*)layers)[index1];
    float2 f2b = ((const float2*)layers)[index1+blockDim.x];
    ftemp[index1] = f2.x;
    ftemp[index1+blockDim.x] = f2b.x;
    __syncthreads();
    float4 first;
    first.x = ftemp[index2+0];
    first.z = ftemp[index2+1];
    __syncthreads();
    ftemp[index1] = f2.y;
    ftemp[index1+blockDim.x] = f2b.y;
    __syncthreads();
    first.y = ftemp[index2+0];
    first.w = ftemp[index2+1];
    for(int L=nlayers-2;L>=0;--L)
    {
        layers += npixels;
//        float4 current = *layers;
        __syncthreads();
        f2 = ((const float2*)layers)[index1];
        f2b = ((const float2*)layers)[index1+blockDim.x];
        ftemp[index1] = f2.x;
        ftemp[index1+blockDim.x] = f2b.x;
        __syncthreads();
        float4 current;
        current.x = ftemp[index2+0];
        current.z = ftemp[index2+1];
        __syncthreads();
        ftemp[index1] = f2.y;
        ftemp[index1+blockDim.x] = f2b.y;
        __syncthreads();
        current.y = ftemp[index2+0];
        current.w = ftemp[index2+1];
        if (current.z!=0.0 && first.z > 0.5 && current.z < 0.5 && floor(first.y*2) != floor(current.y*2))
        {
            *collisions = first; ++collisions;
            *collisions = current; ++collisions;
        }
        first = current;
    }
}

__global__ void SelfCollisionVolume_write_kernel(const float4* layers, const unsigned int* counts, const unsigned int* bloc_input_pos, const unsigned int* bloc_output_pos, float4* collisions, const int nlayers, const int npixels)
{
    //! Dynamically allocated shared memory
    extern  __shared__  unsigned int temp[];
    extern  __shared__  float ftemp[];
    const int index1 = threadIdx.x;
    const int index2 = index1+index1;
    if (index1 == 0)
    {
        temp[0] = bloc_input_pos[blockIdx.x];
        temp[1] = bloc_output_pos[blockIdx.x];
    }
    __syncthreads();
    layers += temp[0]; //+index1;
    collisions += (temp[1] + counts[temp[0]+index1])<<1;
    //float4 first = *layers;
    __syncthreads();
    float2 f2 = ((const float2*)layers)[index1];
    float2 f2b = ((const float2*)layers)[index1+blockDim.x];
    ftemp[index1] = f2.x;
    ftemp[index1+blockDim.x] = f2b.x;
    __syncthreads();
    float4 first;
    first.x = ftemp[index2+0];
    first.z = ftemp[index2+1];
    __syncthreads();
    ftemp[index1] = f2.y;
    ftemp[index1+blockDim.x] = f2b.y;
    __syncthreads();
    first.y = ftemp[index2+0];
    first.w = ftemp[index2+1];
    int inside = floor(2*first.z)*2-1;
    for(int L=nlayers-2;L>=0;--L)
    {
        layers += npixels;
//        float4 current = *layers;
        __syncthreads();
        f2 = ((const float2*)layers)[index1];
        f2b = ((const float2*)layers)[index1+blockDim.x];
        ftemp[index1] = f2.x;
        ftemp[index1+blockDim.x] = f2b.x;
        __syncthreads();
        float4 current;
        current.x = ftemp[index2+0];
        current.z = ftemp[index2+1];
        __syncthreads();
        ftemp[index1] = f2.y;
        ftemp[index1+blockDim.x] = f2b.y;
        __syncthreads();
        current.y = ftemp[index2+0];
        current.w = ftemp[index2+1];
        if (current.z > 0.5)
            ++inside;
        else
        {
            --inside;
            if (current.z!=0.0 && first.z > 0.5 && inside > 0)
            {
                *collisions = first; ++collisions;
                *collisions = current; ++collisions;
            }
        }
        first = current;
    }
}

//void CollisionVolume_write(const void* layers, const void* scan_counts, void* collisions, const int nlayers, const int npixels, const int bsize, const int nblocs, bool self)
void CollisionVolume_write(const void* layers, const void* counts, const void* bloc_input_pos, const void* bloc_output_pos, void* collisions, const int nlayers, const int npixels, const int bsize, const int nblocs, bool self)
{
    dim3 threads(bsize,1);
    dim3 grid(nblocs,1);
    if (self)
        SelfCollisionVolume_write_kernel<<< grid, threads, 2*bsize*sizeof(unsigned int) >>>((const float4*)layers, (const unsigned int*)counts, (const unsigned int*)bloc_input_pos, (const unsigned int*)bloc_output_pos, (float4*)collisions, nlayers, npixels);
    else
        CollisionVolume_write_kernel<<< grid, threads, 2*bsize*sizeof(unsigned int) >>>((const float4*)layers, (const unsigned int*)counts, (const unsigned int*)bloc_input_pos, (const unsigned int*)bloc_output_pos, (float4*)collisions, nlayers, npixels);
}

}
}
}

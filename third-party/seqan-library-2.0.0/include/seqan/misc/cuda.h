// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2013 NVIDIA Corporation
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Jacopo Pantaleoni <jpantaleoni@nvidia.com>
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_MISC_CUDA_MISC_H_
#define SEQAN_MISC_CUDA_MISC_H_

//#include <cuda_runtime.h>
//#include <thrust/version.h>
//#if THRUST_VERSION < 100600
//#include <thrust/detail/backend/cuda/arch.h>
//#else
//#include <thrust/system/cuda/detail/arch.h>
//#endif

namespace seqan {

// ============================================================================
// Classes
// ============================================================================

/*
struct Arch
{
    static const __uint32 LOG_WARP_SIZE = 5;
    static const __uint32 WARP_SIZE     = 1u << LOG_WARP_SIZE;
};
*/

// --------------------------------------------------------------------------
// Type FunctionPointer
// --------------------------------------------------------------------------

typedef void (*FunctionPointer)();


// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function nvccShowType()
// --------------------------------------------------------------------------

template <typename TObject>
inline void
nvccShowType(TObject object)
{
    ignoreUnusedVariableWarning(static_cast<Nothing>(object));
}

// --------------------------------------------------------------------------
// Function divideRI()
// --------------------------------------------------------------------------

// Round x/y towards +infinity for integers, used to determine # of blocks/warps etc.
template<typename L, typename R>
inline L divideRI(const L x, const R y)
{
    return L((x + (y - 1)) / y);
}

// --------------------------------------------------------------------------
// Function roundI()
// --------------------------------------------------------------------------

// Round x towards infinity to the next multiple of y.
template<typename L, typename R>
inline L roundI(L x, R y)
{
    return L(y * divideRI(x, y));
}

// --------------------------------------------------------------------------
// Function cudaPrintFreeMemory()
// --------------------------------------------------------------------------

inline void cudaPrintFreeMemory()
{
    size_t free, total;
    cudaMemGetInfo(&free, &total);

    std::cout << "Free " << free / 1024 / 1024 <<  " of " << total / 1024 / 1024 << " MB\n";
}

// --------------------------------------------------------------------------
// Function cudaOccupancy()
// --------------------------------------------------------------------------

inline SEQAN_DEVICE
float cudaOccupancy()
{
    return __popc(__ballot(true)) / 32.0f;
}

// --------------------------------------------------------------------------
// Function cudaSmemAllocationUnit()
// --------------------------------------------------------------------------
// Granularity of shared memory allocation.

inline size_t cudaSmemAllocationUnit(cudaDeviceProp const & /* properties */)
{
    return 512;
}

// --------------------------------------------------------------------------
// Function cudaRegAllocationUnit()
// --------------------------------------------------------------------------
// Granularity of register allocation.

inline size_t cudaRegAllocationUnit(cudaDeviceProp const & properties)
{
    return (properties.major <= 1) ? (properties.minor <= 1 ? 256 : 512) : 64;
}

// --------------------------------------------------------------------------
// Function cudaWarpAllocationMultiple()
// --------------------------------------------------------------------------
// Granularity of warp allocation.

inline size_t cudaWarpAllocationMultiple(cudaDeviceProp const & /* properties */)
{
    return 2;
}

// --------------------------------------------------------------------------
// Function cudaMaxBlocksPerMultiprocessor()
// --------------------------------------------------------------------------

inline size_t cudaMaxBlocksPerMultiprocessor(cudaDeviceProp const & properties)
{
    return properties.major <= 2 ? 8 : 16;
}

// --------------------------------------------------------------------------
// Function cudaKernelGetAttributes()
// --------------------------------------------------------------------------

template <typename TKernel>
inline cudaFuncAttributes cudaKernelGetAttributes(TKernel kernel)
{
    FunctionPointer kernelPointer = reinterpret_cast<FunctionPointer>(kernel);

    cudaFuncAttributes attributes;
    cudaFuncGetAttributes(&attributes, kernelPointer);

    return attributes;
}

// --------------------------------------------------------------------------
// Function cudaRegistersUsed()
// --------------------------------------------------------------------------

template <typename TKernel>
inline size_t cudaRegistersUsed(TKernel kernel)
{
    return cudaKernelGetAttributes(kernel).numRegs;
}

// --------------------------------------------------------------------------
// Function cudaMaxActiveBlocks()
// --------------------------------------------------------------------------

template <typename TKernel, typename TCTASize, typename TSmemSize>
inline size_t cudaMaxActiveBlocks(TKernel kernel, TCTASize ctaSize, TSmemSize dynamicSmemBytes)
{
    int device;
    cudaGetDevice(&device);

    cudaDeviceProp properties;
    cudaGetDeviceProperties(&properties, device);

    return properties.multiProcessorCount * cudaMaxActiveBlocksPerSM(kernel, ctaSize, dynamicSmemBytes, properties);
}

// --------------------------------------------------------------------------
// Function cudaMaxActiveBlocksPerSM()
// --------------------------------------------------------------------------

template <typename TKernel, typename TCTASize, typename TSmemSize>
inline size_t cudaMaxActiveBlocksPerSM(TKernel kernel, TCTASize ctaSize, TSmemSize dynamicSmemBytes,
                                       cudaDeviceProp & properties)
{
    cudaFuncAttributes attributes = cudaKernelGetAttributes(kernel);

    // Determine the maximum number of CTAs that can be run simultaneously per SM.
    // This is equivalent to the calculation done in the CUDA Occupancy Calculator spreadsheet.
    size_t regAllocationUnit      = cudaRegAllocationUnit(properties);
    size_t warpAllocationMultiple = cudaWarpAllocationMultiple(properties);
    size_t smemAllocationUnit     = cudaSmemAllocationUnit(properties);
    size_t maxThreadsPerSM        = properties.maxThreadsPerMultiProcessor;  // 768, 1024, 1536, etc.
    size_t maxBlocksPerSM         = cudaMaxBlocksPerMultiprocessor(properties);

    // Number of warps (round up to nearest whole multiple of warp size & warp allocation multiple).
    size_t numWarps = roundI(divideRI(ctaSize, properties.warpSize), warpAllocationMultiple);

    // Number of regs is regs per thread times number of warps times warp size.
    size_t regsPerCTA = properties.major < 2 ?
            roundI(attributes.numRegs * properties.warpSize * numWarps, regAllocationUnit) :
            roundI(attributes.numRegs * properties.warpSize, regAllocationUnit) * numWarps;

    size_t smemBytes  = attributes.sharedSizeBytes + dynamicSmemBytes;
    size_t smemPerCTA = roundI(smemBytes, smemAllocationUnit);

    size_t ctaLimitRegs    = regsPerCTA > 0 ? properties.regsPerBlock      / regsPerCTA : maxBlocksPerSM;
    size_t ctaLimitSMem    = smemPerCTA > 0 ? properties.sharedMemPerBlock / smemPerCTA : maxBlocksPerSM;
    size_t ctaLimitThreads =                  maxThreadsPerSM              / ctaSize;

  return _min(ctaLimitRegs, _min(ctaLimitSMem, _min(ctaLimitThreads, maxBlocksPerSM)));
}

// --------------------------------------------------------------------------
// Function checkCudaError()
// --------------------------------------------------------------------------

//inline void checkCudaError(const char *message)
//{
//    cudaError_t error = cudaGetLastError();
//    if(error!=cudaSuccess) {
//        fprintf(stderr,"%s: %s\n", message, cudaGetErrorString(error) );
//        exit(1);
//    }
//}

// --------------------------------------------------------------------------
// Functions not yet used.
// --------------------------------------------------------------------------

/*
template <typename KernelFunction>
size_t auto_blocksize(KernelFunction kernel, size_t dynamic_smem_bytes_per_thread = 0)
{
#if THRUST_VERSION < 100600
    return thrust::detail::backend::cuda::arch::max_blocksize_with_highest_occupancy(kernel, dynamic_smem_bytes_per_thread);
#else
    return thrust::system::cuda::detail::arch::max_blocksize_with_highest_occupancy(kernel, dynamic_smem_bytes_per_thread);
#endif
}

inline bool is_tcc_enabled()
{
    int            device;
    cudaDeviceProp device_properties;
    cudaGetDevice(&device);
    cudaGetDeviceProperties( &device_properties, device );
    return device_properties.tccDriver ? true : false;
}

/// a generic syncthreads() implementation to synchronize contiguous
/// blocks of N threads at a time
///
template <__uint32 N>
SEQAN_HOST_DEVICE inline
void syncThreads()
{
    #ifdef __CUDA_ARCH__
    if ((N > cuda::Arch::WARP_SIZE) || (is_pow2<N>() == false))
        __syncthreads();
    #endif
}
*/

}

#endif  // SEQAN_MISC_CUDA_MISC_H_

#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable

KEY_TYPE getValue(DATA_TYPE value) {
    return SORT_KEY;
}

/**
 * Sort a list that is short enough to entirely fit in local memory.  This is executed as
 * a single thread block.
 */
__kernel void sortShortList(__global DATA_TYPE* __restrict__ data, uint length, __local DATA_TYPE* dataBuffer) {
    // Load the data into local memory.
    
    for (int index = get_local_id(0); index < length; index += get_local_size(0))
        dataBuffer[index] = data[index];
    barrier(CLK_LOCAL_MEM_FENCE);

    // Perform a bitonic sort in local memory.

    for (unsigned int k = 2; k < 2*length; k *= 2) {
        for (unsigned int j = k/2; j > 0; j /= 2) {
            for (unsigned int i = get_local_id(0); i < length; i += get_local_size(0)) {
                int ixj = i^j;
                if (ixj > i && ixj < length) {
                    DATA_TYPE value1 = dataBuffer[i];
                    DATA_TYPE value2 = dataBuffer[ixj];
                    bool ascending = ((i&k) == 0);
                    for (unsigned int mask = k*2; mask < 2*length; mask *= 2)
                        ascending = ((i&mask) == 0 ? !ascending : ascending);
                    KEY_TYPE lowKey  = (ascending ? getValue(value1) : getValue(value2));
                    KEY_TYPE highKey = (ascending ? getValue(value2) : getValue(value1));
                    if (lowKey > highKey) {
                        dataBuffer[i] = value2;
                        dataBuffer[ixj] = value1;
                    }
                }
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }
    }

    // Write the data back to global memory.

    for (int index = get_local_id(0); index < length; index += get_local_size(0))
        data[index] = dataBuffer[index];
}

/**
 * Calculate the minimum and maximum value in the array to be sorted.  This kernel
 * is executed as a single work group.
 */
__kernel void computeRange(__global const DATA_TYPE* restrict data, uint length, __global KEY_TYPE* restrict range, __local KEY_TYPE* restrict buffer) {
    KEY_TYPE minimum = MAX_KEY;
    KEY_TYPE maximum = MIN_KEY;

    // Each thread calculates the range of a subset of values.

    for (uint index = get_local_id(0); index < length; index += get_local_size(0)) {
        KEY_TYPE value = getValue(data[index]);
        minimum = min(minimum, value);
        maximum = max(maximum, value);
    }

    // Now reduce them.

    buffer[get_local_id(0)] = minimum;
    barrier(CLK_LOCAL_MEM_FENCE);
    for (uint step = 1; step < get_local_size(0); step *= 2) {
        if (get_local_id(0)+step < get_local_size(0) && get_local_id(0)%(2*step) == 0)
            buffer[get_local_id(0)] = min(buffer[get_local_id(0)], buffer[get_local_id(0)+step]);
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    minimum = buffer[0];
    barrier(CLK_LOCAL_MEM_FENCE);
    buffer[get_local_id(0)] = maximum;
    barrier(CLK_LOCAL_MEM_FENCE);
    for (uint step = 1; step < get_local_size(0); step *= 2) {
        if (get_local_id(0)+step < get_local_size(0) && get_local_id(0)%(2*step) == 0)
            buffer[get_local_id(0)] = max(buffer[get_local_id(0)], buffer[get_local_id(0)+step]);
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    maximum = buffer[0];
    if (get_local_id(0) == 0) {
        range[0] = minimum;
        range[1] = maximum;
    }
}

/**
 * Assign elements to buckets.
 */
__kernel void assignElementsToBuckets(__global const DATA_TYPE* restrict data, uint length, uint numBuckets, __global const KEY_TYPE* restrict range,
        __global uint* bucketOffset, __global uint* restrict bucketOfElement, __global uint* restrict offsetInBucket) {
#ifdef AMD_ATOMIC_WORK_AROUND
    // Do a byte write to force all memory accesses to interactionCount to use the complete path.
    // This avoids the atomic access from causing all word accesses to other buffers from using the slow complete path.
    // The IF actually causes the write to never be executed, its presence is all that is needed.
    // AMD APP SDK 2.4 has this problem.
    if (get_global_id(0) == get_local_id(0)+1)
        ((__global char*)bucketOffset)[sizeof(int)*numBuckets+1] = 0;
#endif
    float minValue = (float) (range[0]);
    float maxValue = (float) (range[1]);
    float bucketWidth = (maxValue-minValue)/numBuckets;
    for (uint index = get_global_id(0); index < length; index += get_global_size(0)) {
#if defined(MAC_AMD_WORKAROUND) && VALUE_IS_INT2
        __global int* d = (__global int*) data;
        int2 element = (int2) (d[2*index], d[2*index+1]);
        float key = (float) getValue(element);
#else
        float key = (float) getValue(data[index]);
#endif
        uint bucketIndex = min((uint) ((key-minValue)/bucketWidth), numBuckets-1);
        offsetInBucket[index] = atom_inc(&bucketOffset[bucketIndex]);
        bucketOfElement[index] = bucketIndex;
    }
}

/**
 * Sum the bucket sizes to compute the start position of each bucket.  This kernel
 * is executed as a single work group.
 */
__kernel void computeBucketPositions(uint numBuckets, __global uint* restrict bucketOffset, __local uint* restrict buffer) {
    uint globalOffset = 0;
    for (uint startBucket = 0; startBucket < numBuckets; startBucket += get_local_size(0)) {
        // Load the bucket sizes into local memory.

        uint globalIndex = startBucket+get_local_id(0);
        buffer[get_local_id(0)] = (globalIndex < numBuckets ? bucketOffset[globalIndex] : 0);
        barrier(CLK_LOCAL_MEM_FENCE);

        // Perform a parallel prefix sum.

        for (uint step = 1; step < get_local_size(0); step *= 2) {
            uint add = (get_local_id(0) >= step ? buffer[get_local_id(0)-step] : 0);
            barrier(CLK_LOCAL_MEM_FENCE);
            buffer[get_local_id(0)] += add;
            barrier(CLK_LOCAL_MEM_FENCE);
        }

        // Write the results back to global memory.

        if (globalIndex < numBuckets)
            bucketOffset[globalIndex] = buffer[get_local_id(0)]+globalOffset;
        globalOffset += buffer[get_local_size(0)-1];
    }
}

/**
 * Copy the input data into the buckets for sorting.
 */
__kernel void copyDataToBuckets(__global const DATA_TYPE* restrict data, __global DATA_TYPE* restrict buckets, uint length, __global const uint* restrict bucketOffset, __global const uint* restrict bucketOfElement, __global const uint* restrict offsetInBucket) {
    for (uint index = get_global_id(0); index < length; index += get_global_size(0)) {
        DATA_TYPE element = data[index];
        uint bucketIndex = bucketOfElement[index];
        uint offset = (bucketIndex == 0 ? 0 : bucketOffset[bucketIndex-1]);
        buckets[offset+offsetInBucket[index]] = element;
    }
}

/**
 * Sort the data in each bucket.
 */
__kernel void sortBuckets(__global DATA_TYPE* restrict data, __global const DATA_TYPE* restrict buckets, uint numBuckets, __global const uint* restrict bucketOffset, __local DATA_TYPE* restrict buffer) {
    for (uint index = get_group_id(0); index < numBuckets; index += get_num_groups(0)) {
        uint startIndex = (index == 0 ? 0 : bucketOffset[index-1]);
        uint endIndex = bucketOffset[index];
        uint length = endIndex-startIndex;
        if (length <= get_local_size(0)) {
            // Load the data into local memory.

            if (get_local_id(0) < length)
                buffer[get_local_id(0)] = buckets[startIndex+get_local_id(0)];
            else
                buffer[get_local_id(0)] = MAX_VALUE;
            barrier(CLK_LOCAL_MEM_FENCE);

            // Perform a bitonic sort in local memory.

            for (uint k = 2; k <= get_local_size(0); k *= 2) {
                for (uint j = k/2; j > 0; j /= 2) {
                    int ixj = get_local_id(0)^j;
                    if (ixj > get_local_id(0)) {
                        DATA_TYPE value1 = buffer[get_local_id(0)];
                        DATA_TYPE value2 = buffer[ixj];
                        bool ascending = (get_local_id(0)&k) == 0;
                        KEY_TYPE lowKey = (ascending ? getValue(value1) : getValue(value2));
                        KEY_TYPE highKey = (ascending ? getValue(value2) : getValue(value1));
                        if (lowKey > highKey) {
                            buffer[get_local_id(0)] = value2;
                            buffer[ixj] = value1;
                        }
                    }
                    barrier(CLK_LOCAL_MEM_FENCE);
                }
            }

            // Write the data to the sorted array.

            if (get_local_id(0) < length)
                data[startIndex+get_local_id(0)] = buffer[get_local_id(0)];
        }
        else {
            // Copy the bucket data over to the output array.

            for (uint i = get_local_id(0); i < length; i += get_local_size(0))
                data[startIndex+i] = buckets[startIndex+i];
            barrier(CLK_GLOBAL_MEM_FENCE);

            // Perform a bitonic sort in global memory.

            for (uint k = 2; k < 2*length; k *= 2) {
                for (uint j = k/2; j > 0; j /= 2) {
                    for (uint i = get_local_id(0); i < length; i += get_local_size(0)) {
                        int ixj = i^j;
                        if (ixj > i && ixj < length) {
                            DATA_TYPE value1 = data[startIndex+i];
                            DATA_TYPE value2 = data[startIndex+ixj];
                            bool ascending = ((i&k) == 0);
                            for (uint mask = k*2; mask < 2*length; mask *= 2)
                                ascending = ((i&mask) == 0 ? !ascending : ascending);
                            KEY_TYPE lowKey  = (ascending ? getValue(value1) : getValue(value2));
                            KEY_TYPE highKey = (ascending ? getValue(value2) : getValue(value1));
                            if (lowKey > highKey) {
                                data[startIndex+i] = value2;
                                data[startIndex+ixj] = value1;
                            }
                        }
                    }
                    barrier(CLK_GLOBAL_MEM_FENCE);
                }
            }
        }
    }
}

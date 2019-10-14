#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <OpenCL/opencl.h>

//const char *KernelSource = "\n" \
//"inline float add(float a,float b)										\n" \
//"{																		\n" \
//"   return a+b;                                                         \n" \
//"}                                                                      \n" \
//"__kernel void square(                                                  \n" \
//"   __global float* input,                                              \n" \
//"   __global float* output,                                             \n" \
//"   const unsigned int count)                                           \n" \
//"{                                                                      \n" \
//"   int i = get_global_id(0);                                           \n" \
//"   if(i < count)                                                       \n" \
//"       //output[i] = input[i] * input[i];                              \n" \
//"       output[i] = add(input[i],input[i]);                             \n" \
//"}                                                                      \n" \
//"\n";
//
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
    int err;                            // error code returned from api calls

	int DATA_SIZE = 1024*1024;
    float *data = (float *)malloc(sizeof(float)*DATA_SIZE);              // original data set given to device
    float *results = (float *)malloc(sizeof(float)*DATA_SIZE);              // original data set given to device
    unsigned int correct;               // number of correct results returned

    FILE* programHandle;
    size_t programSize, KernelSourceSize;
    char *programBuffer, *KernelSource;

    size_t global;                      // global domain size for our calculation
    size_t local;                       // local domain size for our calculation

    cl_device_id device_id;             // compute device id
    cl_context context;                 // compute context
    cl_command_queue commands;          // compute command queue
    cl_program program;                 // compute program
    cl_kernel kernel;                   // compute kernel

    cl_mem input;                       // device memory used for the input array
    cl_mem output;                      // device memory used for the output array

    int i = 0;
    unsigned int count = DATA_SIZE;
    for(i = 0; i < count; i++) data[i] = rand() / (float)RAND_MAX;

    int gpu = 1;
    err = clGetDeviceIDs(NULL, gpu ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU, 1, &device_id, NULL);
    context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
    commands = clCreateCommandQueue(context, device_id, 0, &err);

	//----------------------------------------------------------------------------
    // get size of kernel source
    programHandle = fopen("./simple_add.cl", "r");
    fseek(programHandle, 0, SEEK_END);
    programSize = ftell(programHandle);
    rewind(programHandle);

    // read kernel source into buffer
    programBuffer = (char*) malloc(programSize + 1);
    programBuffer[programSize] = '\0';
    fread(programBuffer, sizeof(char), programSize, programHandle);
    fclose(programHandle);

    // create program from buffer
    program = clCreateProgramWithSource(context, 1,
            (const char**) &programBuffer, &programSize, NULL);
    free(programBuffer);

    // read kernel source back in from program to check
    clGetProgramInfo(program, CL_PROGRAM_SOURCE, 0, NULL, &KernelSourceSize);
    KernelSource = (char*) malloc(KernelSourceSize);
    clGetProgramInfo(program, CL_PROGRAM_SOURCE, KernelSourceSize, KernelSource, NULL);

    program = clCreateProgramWithSource(context, 1, (const char **) & KernelSource, NULL, &err);
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    kernel = clCreateKernel(program, "square", &err);
	//----------------------------------------------------------------------------

    input = clCreateBuffer(context,  CL_MEM_READ_ONLY,  sizeof(float) * count, NULL, NULL);
    output = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * count, NULL, NULL);

    err = clEnqueueWriteBuffer(commands, input, CL_TRUE, 0, sizeof(float) * count, data, 0, NULL, NULL);
    err = 0;
    err  = clSetKernelArg(kernel, 0, sizeof(cl_mem), &input);
    err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &output);
    err |= clSetKernelArg(kernel, 2, sizeof(unsigned int), &count);

    err = clGetKernelWorkGroupInfo(kernel, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(local), &local, NULL);
    global = count;
    err = clEnqueueNDRangeKernel(commands, kernel, 1, NULL, &global, &local, 0, NULL, NULL);
    clFinish(commands);
    err = clEnqueueReadBuffer( commands, output, CL_TRUE, 0, sizeof(float) * count, results, 0, NULL, NULL );

    correct = 0;
    for(i = 0; i < count; i++)
    {
        //if(results[i] == data[i] * data[i])
        if(results[i] == data[i]+data[i])
            correct++;
    }
    printf("Computed '%d/%d' correct values!\n", correct, count);

    clReleaseMemObject(input);
    clReleaseMemObject(output);
    clReleaseProgram(program);
    clReleaseKernel(kernel);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);
    free(KernelSource);

    return 0;
}

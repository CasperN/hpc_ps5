#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <curand.h>
#include <sys/time.h>

#define CENTER (vec_t) {0, 12, 0}
#define LIGHT  (vec_t) {4, 4, -1}
#define RADIUS 6
#define WINDOW_Y 10
#define WINDOW_MAX 10
#define BLOCKS_PER_DIM 22
#define THREADS_PER_BLOCK_DIM 8

typedef struct { float x,y,z; } vec_t;

__device__ vec_t add(vec_t v, vec_t w){
    vec_t res;
    res.x = v.x + w.x;
    res.y = v.y + w.y;
    res.z = v.z + w.z;
    return res;
}

__device__ vec_t scale(float c, vec_t v){
    v.x *= c;
    v.y *= c;
    v.z *= c;
    return v;
}


__device__ float dot(vec_t u, vec_t v){
    return (u.x * v.x) + (u.y * v.y) + (u.z * v.z);
}

__device__ vec_t direction(vec_t end, vec_t start){
    vec_t dir = add(end, scale(-1, start));
    float norm = sqrt(dot(dir, dir));
    return scale(1.0 / norm, dir);
}

__device__ int getGlobalIdx(){
    // https://cs.calvin.edu/courses/cs/374/CUDA/CUDA-Thread-Indexing-Cheatsheet.pdf
    int threadId, blockId;

    blockId = blockIdx.x
            + blockIdx.y * gridDim.x
            + blockIdx.z * gridDim.x * gridDim.y;

    threadId = blockId * blockDim.x * blockDim.y * blockDim.z
             + threadIdx.x
             + threadIdx.y * blockDim.x
             + threadIdx.z * blockDim.x * blockDim.y;

    return threadId;
}

__device__ vec_t sample_direction(curandState_t *states){
    float phi, cos_th, sin_th, r1, r2;
    vec_t v;
    int idx;

    idx = getGlobalIdx();
    r1 = curand_uniform(&states[idx]);
    r2 = curand_uniform(&states[idx]);

    phi    = r1 * 2 * M_PI;
    cos_th = r2 * 2 - 1;
    sin_th = sqrt(1 - cos_th * cos_th);

    v.x = sin_th * cos(phi);
    v.y = sin_th * sin(phi);
    v.z = cos_th;

    return v;
}


typedef struct {
    float *surface, max, y;    // Windows are always parallel to y axis since
    int size;                   // the viewer is assumed to be at the origin
} window_t;                     // facing in y direction


void init_window(window_t *w, int size, float y_coordinate, float w_max){
    w->y = y_coordinate;
    w->max = w_max;
    w->size = size;
    w->surface = NULL;
    // surface is allocated in gpu then copied over here
}

void destroy_window(window_t *w){
    free(w->surface);
}

__device__ int on_surface(window_t *w, vec_t co_ord) {
    return (abs(co_ord.x) < w->max) && (abs(co_ord.z) < w->max);
}

void save_surface(window_t *w, const char *file_name){
    printf("Saving to `%s`\n", file_name);
    FILE *f = fopen(file_name,"w");
    fwrite(w->surface, sizeof(float), w->size * w->size, f);
    fclose(f);
    printf("Done saving\n");
}

__device__ void add_brightness(window_t *w, vec_t co_ord, float b){
    int x,y;
    x = (co_ord.x + w->max) / 2 / w->max * w->size;
    y = (co_ord.z + w->max) / 2 / w->max * w->size;
    atomicAdd(w->surface + x * w->size + y, b);
}


// Ray tracing simulation of a sphere at coordinates `center` and radius `radius`.
// Assumes observer is at origin facing in y direction dowards `window`
__global__ void kernel(window_t window, float* surface, vec_t center,
                       vec_t light, float radius, curandState_t *states,
                       long n_rays, unsigned int seed){

    float radius_sq, center_sq;
    int idx = getGlobalIdx();
    curand_init(seed + idx, idx, 0, &states[idx]);

    window.surface = surface;
    radius_sq = radius * radius;
    center_sq = dot(center, center);

    for(long i=0; i<n_rays; i++){
        vec_t ray, ray_window, ray_sphere, sph_dir, light_dir;
        float brightness, t, x = -1;
        do{
            ray = sample_direction(states);
            ray_window = scale(window.y / ray.y, ray);
            if(!on_surface(&window, ray_window)) continue;
            t = dot(ray, center);
            x = t * t + radius_sq - center_sq;

        } while(x < 0);

        ray_sphere = scale(t - sqrt(x), ray);

        sph_dir = direction(ray_sphere, center);
        light_dir = direction(light, ray_sphere);

        brightness = dot(sph_dir, light_dir);
        brightness = brightness > 0 ? brightness : 0;
        add_brightness(&window, ray_window, brightness);
    }
}


void gpu_ray_tracing(window_t* window, vec_t center, vec_t light, float radius, long n_rays){

    curandState_t* d_states;
    struct timeval start, end;
    size_t surface_size;
    float *d_surface;
    dim3 block, grid;
    long elapsed;
    int n_par;

    block.x = block.y = block.z = THREADS_PER_BLOCK_DIM;
    grid.x = grid.y = grid.z = BLOCKS_PER_DIM;
    n_par = block.x * block.y * block.z * grid.x * grid.y * grid.z;

    printf("Running with %d threads/block, %d blocks/grid\n", block.x, grid.x);
    gettimeofday(&start, 0);

    // Initialize random number seeds and window->surface
    cudaMalloc((void**) &d_states, n_par * sizeof(curandState_t));
    surface_size = window->size * window->size * sizeof(float);
    cudaMalloc((void**) &d_surface, surface_size);
    cudaMemset(&d_surface, 0, surface_size);

    kernel<<<grid, block>>>(
        *window, d_surface, center, light, radius, d_states,
        ceil( (float)n_rays / n_par), time(NULL));

    // Collect surface
    window->surface = (float*) malloc(surface_size);
    cudaMemcpy(window->surface, d_surface, surface_size, cudaMemcpyDeviceToHost);
    cudaFree(d_surface);
    cudaFree(d_states);

    // Log performance and cuda errors
    gettimeofday(&end, 0);
    elapsed = (end.tv_sec - start.tv_sec)*1e6 + end.tv_usec - start.tv_usec;
    printf("Done in %ld usec. Last cuda error: %s\n",
        elapsed, cudaGetErrorName(cudaGetLastError()));

    for(int i=0; i < window->size * window->size; i++)
        if (window->surface[i]) return;
    printf("Warning: Something happened, its all zeros\n");
}


int main(int argc, char const *argv[]) {

    long n_rays;
    int pixels_per_side;
    window_t w;

    if (argc != 3){
        printf("Usage: ./ray_trace [n_rays] [pixels_per_side]\n");
        exit(0);
    }
    n_rays = atol(argv[1]);
    pixels_per_side = atoi(argv[2]);

    init_window(&w, pixels_per_side, WINDOW_Y, WINDOW_MAX);

    gpu_ray_tracing(&w, CENTER, LIGHT, RADIUS, n_rays);

    save_surface(&w, (const char *)"result.bin");
    destroy_window(&w);

    return 0;
}

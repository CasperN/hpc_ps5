#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <cuda.h>

#define CENTER (vec_t) {0, 12, 0}
#define LIGHT (vec_t) {4, 4, -1}
#define RADIUS 6
#define WINDOW_Y 10
#define WINDOW_MAX 10

typedef struct { double x,y,z; } vec_t;

__device__ vec_t Z_Vec = (vec_t) {0,0,0};

__device__ vec_t add(vec_t v, vec_t w){
    vec_t res;
    res.x = v.x + w.x;
    res.y = v.y + w.y;
    res.z = v.z + w.z;
    return res;
}

__device__ vec_t scale(double c, vec_t v){
    v.x *= c;
    v.y *= c;
    v.z *= c;
    return v;
}


__device__ double dot(vec_t u, vec_t v){
    return (u.x * v.x) + (u.y * v.y) + (u.z * v.z);
}

__device__ vec_t direction(vec_t end, vec_t start){
    vec_t dir = add(end, scale(-1, start));
    double norm = sqrt(dot(dir, dir));
    return scale(1/norm, dir);
}

__device__ vec_t sample_direction(){
    double phi, cos_th, sin_th;
    vec_t v;
    phi    = (double) rand() / (double) RAND_MAX * 2 * M_PI;
    cos_th = (double) rand() / (double) RAND_MAX * 2 - 1;
    sin_th = sqrt(1 - cos_th * cos_th);

    v.x = sin_th * cos(phi);
    v.y = sin_th * sin(phi);
    v.z = cos_th;
    return v;
}


typedef struct {
    double *surface, max, y;    // Windows are always parallel to y axis since
    int size;                   // the viewer is assumed to be at the origin
} window_t;                     // facing in y direction


void init_window(window_t *w, int size, double y_coordinate, double w_max){
    w->y = y_coordinate;
    w->max = w_max;
    w->size = size;
    w->surface = NULL;
    // surface is allocated in gpu then copied over here
}

void destroy_window(window_t *w){
    free(w->surface);
}

int on_surface(window_t *w, vec_t co_ord) {
    return (abs(co_ord.x) < w->max) && (abs(co_ord.z) < w->max);
}

void save_surface(window_t *w,  char *file_name){
    FILE *f = fopen(file_name,"w");
    fwrite(w->surface, sizeof(double), w->size * w->size, f);
    fclose(f);
}

__device__ void add_brightness(window_t *w, vec_t co_ord, double b){
    int x,y;
    x = (co_ord.x + w->max) / 2 / w->max * w->size;
    y = (co_ord.z + w->max) / 2 / w->max * w->size;
    w->surface[x * w->size + y] += b;
}

// Ray tracing simulation of a sphere at coordinates `center` and radius `radius`.
// Assumes observer is at origin facing in y direction dowards `window`
__global__ void trace_ray(window_t window, double* surface, vec_t center, vec_t light, double radius){

    vec_t ray, ray_window, ray_sphere, sph_dir, light_dir;
    double radius_sq, center_sq,brightness, t, x = -1;

    window.surface = surface;

    radius_sq = radius * radius;
    center_sq = dot(center, center);

    do{
        ray = sample_direction();                   // Ray originates from the origin
        ray_window = scale(window->y / ray.y, ray); // Ray's intersection with window

        if(!on_surface(&window, ray_window)) continue;

        t = dot(ray, center);                       // t and x are algebra steps needed
        x = t * t + radius_sq - center_sq;          // to solve for ray-sphere intersection

    } while(x < 0);

    ray_sphere = scale(t - sqrt(x), ray);

    sph_dir = direction(ray_sphere, center);
    light_dir = direction(light, ray_sphere);

    brightness = dot(sph_dir, light_dir);
    brightness = brightness > 0 ? brightness : 0;

    add_brightness(&window, ray_window, brightness);
}

void gpu_ray_tracing(window_t* window, vec_t center, vec_t light, double radius, int n_rays){

    int nblocks, th_per_bl;
    double *d_surface;

    nblocks = n_rays;
    th_per_bl = 128;
    size_t surface_size = window->size * window->size, sizeof(double);

    // Allocate memory for surface on gpu
    cudaMalloc(&d_surface, surface_size);
    cudaMemset(&d_surface, 0, surface_size);
    // Launch ray tracing
    trace_ray<<<nblocks, th_per_bl>>>(window, d_surface, center, light, radius);
    // Collect surface
    window->surface = malloc(surface_size);
    cudaMemcpy(&window->surface, d_surface, surface_size, cudaMemcpyDeviceToHost);
    cudaFree(d_surface);
}



int main(int argc, char const *argv[]) {

    int n_rays, pixels_per_side;
    window_t w;

    if (argc != 3){
        printf("Usage: ./ray_trace [n_rays] [pixels_per_side]\n");
        exit(0);
    }
    n_rays = atoi(argv[1]);
    pixels_per_side = atoi(argv[2]);

    init_window(&w, pixels_per_side, WINDOW_Y, WINDOW_MAX);

    gpu_ray_tracing(&w, CENTER, LIGHT, RADIUS, n_rays);

    save_surface(&w, "result.bin");
    destroy_window(&w);

    return 0;
}

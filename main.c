#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <assert.h>

#define CENTER (vec_t) {0, 12, 0}
#define LIGHT (vec_t) {4, 4, -1}
#define RADIUS 6
#define WINDOW_Y 10
#define WINDOW_MAX 10

typedef struct { double x,y,z; } vec_t;

vec_t Z_Vec = (vec_t) {0,0,0};

vec_t add(double a, vec_t v, double b, vec_t w){
    vec_t res;
    res.x = a * v.x + b * w.x;
    res.y = a * v.y + b * w.y;
    res.z = a * v.z + b * w.z;
    return res;
}


double dot(vec_t u, vec_t v){
    return (u.x * v.x) + (u.y * v.y) + (u.z * v.z);
}

vec_t direction(vec_t end, vec_t start){
    vec_t dir = add(1, end, -1, start);
    double norm = sqrt(dot(dir, dir));
    return add(1/norm, dir, 0, Z_Vec);
}

vec_t sample_direction(){
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


// Windows are always parallel to y axis since the viewer is assumed to be
// at the origin facing in y direction
typedef struct { double **surface, max, y; int size; } window_t;

void init_window(window_t *w, int size, double y_coordinate, double w_max){
    w->y = y_coordinate;
    w->max = w_max;
    w->size = size;
    w->surface = malloc(size * sizeof(double*));
    w->surface[0] = calloc(size * size, sizeof(double));
    assert(w->surface);
    assert(w->surface[0]);
    for(int i=0; i<size; i++)
        w->surface[i] = w->surface[0] + i * size;
}

void destroy_window(window_t *w){
    free(w->surface[0]);
    free(w->surface);
}

int on_surface(window_t *w, vec_t co_ord) {
    return (abs(co_ord.x) < w->max) && (abs(co_ord.z) < w->max);
}

void save_surface(window_t *w,  char *file_name){
    FILE *f = fopen(file_name,"w");
    fwrite(w->surface[0], sizeof(double), w->size * w->size, f);
    fclose(f);
}

void add_brightness(window_t *w, vec_t co_ord, double b){
    int x,y;
    x = (co_ord.x + w->max) / 2 / w->max * w->size;
    y = (co_ord.z + w->max) / 2 / w->max * w->size;
    w->surface[x][y] += b;
}

// Ray tracing simulation of a sphere at coordinates `center` and radius `radius`.
// Assumes observer is at origin facing in y direction dowards `window`
void serial_ray_tracing(window_t* window, vec_t center, vec_t light, double radius, int n_rays){

    double radius_sq, center_sq;
    radius_sq = radius * radius;
    center_sq = dot(center, center);

    for(int n=0; n<n_rays; n++){
        vec_t ray, ray_window, ray_sphere, sph_dir, light_dir;
        double brightness, t, x = -1;

        do{
            ray = sample_direction();                           // Ray originates from the origin
            ray_window = add(window->y / ray.y, ray, 0, Z_Vec); // Ray's intersection with window

            if(!on_surface(window, ray_window)) continue;

            t = dot(ray, center);                   // t and x are algebra steps needed
            x = t * t + radius_sq - center_sq;      // to solve for ray-sphere intersection

        } while(x < 0);

        // `t` has to be initialized to exit the loop as `x` is initialized to be negative
        #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
        ray_sphere = add(t - sqrt(x), ray, 0, Z_Vec);
        #pragma GCC diagnostic pop

        sph_dir = direction(ray_sphere, center);
        light_dir = direction(light, ray_sphere);

        brightness = dot(sph_dir, light_dir);
        brightness = brightness > 0 ? brightness : 0;

        add_brightness(window, ray_window, brightness);
    }
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

    serial_ray_tracing(&w, CENTER, LIGHT, RADIUS, n_rays);

    save_surface(&w, "result.bin");
    destroy_window(&w);

    return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>

#define CENTER (vec_t) {0, 12, 0}
#define LIGHT (vec_t) {4, 4, -1}
#define RADIUS 6
#define WINDOW_Y 10
#define WINDOW_MAX 10
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


typedef struct { float x,y,z; } vec_t;

vec_t add(vec_t v, vec_t w){
    vec_t res;
    res.x = v.x + w.x;
    res.y = v.y + w.y;
    res.z = v.z + w.z;
    return res;
}

vec_t scale(float c, vec_t v){
    v.x *= c;
    v.y *= c;
    v.z *= c;
    return v;
}

float dot(vec_t u, vec_t v){
    return (u.x * v.x) + (u.y * v.y) + (u.z * v.z);
}

vec_t direction(vec_t end, vec_t start){
    vec_t dir = add(end, scale(-1, start));
    float norm = sqrt(dot(dir, dir));
    return scale(1/norm, dir);
}

vec_t sample_direction(){
    float phi, cos_th, sin_th;
    vec_t v;
    phi    = (float) rand() / (float) RAND_MAX * 2 * M_PI;
    cos_th = (float) rand() / (float) RAND_MAX * 2 - 1;
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
    w->surface = calloc(size * size, sizeof(float));
    assert(w->surface);
}

void destroy_window(window_t *w){
    free(w->surface);
}

int on_surface(window_t *w, vec_t co_ord) {
    return (fabs(co_ord.x) < w->max) && (fabs(co_ord.z) < w->max);
}

void save_surface(window_t *w,  char *file_name){
    FILE *f = fopen(file_name,"w");
    fwrite(w->surface, sizeof(float), w->size * w->size, f);
    fclose(f);
}

void add_brightness(window_t *w, vec_t co_ord, float b){
    int x,y;
    x = (co_ord.x + w->max) / 2 / w->max * w->size;
    y = (co_ord.z + w->max) / 2 / w->max * w->size;
    w->surface[x * w->size + y] += b;
}

// Ray tracing simulation of a sphere at coordinates `center` and radius `radius`.
// Assumes observer is at origin facing in y direction dowards `window`
void serial_ray_tracing(window_t* window, vec_t center, vec_t light, float radius, int n_rays){

    float radius_sq, center_sq;
    radius_sq = radius * radius;
    center_sq = dot(center, center);

    for(int n=0; n<n_rays; n++){
        vec_t ray, ray_window, ray_sphere, sph_dir, light_dir;
        float brightness, t, x = -1;

        do{
            ray = sample_direction();                     // Ray originates from the origin
            ray_window = scale(window->y / ray.y, ray); // Ray's intersection with window
            t = dot(ray, center);                   // t and x are algebra steps needed
            x = t * t + radius_sq - center_sq;      // to solve for ray-sphere intersection
        } while(!on_surface(window, ray_window) || x < 0);

        ray_sphere = scale(t - sqrt(x), ray);

        sph_dir = direction(ray_sphere, center);
        light_dir = direction(light, ray_sphere);

        brightness = dot(sph_dir, light_dir);
        brightness = brightness > 0 ? brightness : 0;

        add_brightness(window, ray_window, brightness);
    }
}


int main(int argc, char const *argv[]) {

    int n_rays, pixels_per_side;
    struct timeval start, end;
    long elapsed;
    window_t w;

    if (argc != 3){
        printf("Usage: ./ray_trace [n_rays] [pixels_per_side]\n");
        exit(0);
    }
    n_rays = atoi(argv[1]);
    pixels_per_side = atoi(argv[2]);

    init_window(&w, pixels_per_side, WINDOW_Y, WINDOW_MAX);

    gettimeofday(&start, 0);
    serial_ray_tracing(&w, CENTER, LIGHT, RADIUS, n_rays);
    gettimeofday(&end, 0);

    elapsed = (end.tv_sec - start.tv_sec)*1e6 + end.tv_usec - start.tv_usec;
    printf("Done in %ld usec\n", elapsed);

    save_surface(&w, "result.bin");
    destroy_window(&w);

    return 0;
}

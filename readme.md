---
geometry: margin=2cm
output: pdf-latex
---
# High Performance Computing Problem Set 5

## Usage
See makefile


### Block configurations

Since I used a 3 dimension configuration, there are actually $t^3$ threads per
block and $b^3$ blocks in the grid. In the table "threads" and "blocks"
refer to the number of threads / blocks in each dimension. The actual amount of
parallelism is "$(tb)^3$". The following tests are show with 100,000,000 rays
onto a 1000 by 1000 pixel surface.


threads      blocks       seconds
--------    -------     ---------
8                 4      0.799883
8                 8      0.644736
8                10      2.674156
8                16     16.087933
4                 4      1.316118
4                 8      0.809381
4                16      0.841489
4                32     16.189896
2                 8      1.940364
2                16      1.927837
2                32      2.492310
1                32      5.022130
1                64      7.924116

Based on this results, it seems 8^3 blocks of 8^3 threads yields optimal
results. I think the reason more blocks seem to reduce performance is because
the initialization time of `curand_init` dominates the computation of the rays.
I've written my program such that each thread has its own random seed but reuses
it in its computation of $\frac{N}{b^3t^3}$ rays where $N$ is the total number
of rays.

### Performance Comparison
The gpu here is configured to use 8 blocks and 8 threads for each of 3
dimensions (ie 262,144 way parallelism). Observe the GPU seemingly scales better
for small numbers of rays, this is because of the high initialization time
associated with `curand_init`, and copying data back and forth. This
initialization time is insignificant at the $10^9$ scale.
The CPU quickly takes too long for me to wait for it, but it is clearly much
slower than the gpu.

rays         GPU (Secs)   CPU (Secs)
---------   -----------  -----------
$10^{10}$     25.222563         ?
$10^{9 }$      2.871905         ?
$10^{8 }$      0.718318    71.292471
$10^{7 }$      0.495271     7.282018
$10^{6 }$      0.467835     0.733851


### GPU Image Output

![gpu_fig.png]($10{^10}$ rays, 2500^2 pixels)

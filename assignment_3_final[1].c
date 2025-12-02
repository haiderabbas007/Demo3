/*
PX390 Assignment 3 – Bug Fix Notes
==================================

-----------------------------------------------------
A. Header / compilation issues
-----------------------------------------------------

1) Missing stdio header
   Original code only included:
       #include <math.h>
       #include <stdlib.h>
   but used printf, FILE, fopen, fscanf, fclose.

   Problem:
   - Without <stdio.h>, these functions and FILE are implicitly declared.

   Fix:
       #include <stdio.h>
       #include <math.h>
       #include <stdlib.h>

   This makes printf, FILE, fopen, fscanf, fclose visible and removes the
   "implicit declaration" / "unknown type name 'FILE'" errors.

-----------------------------------------------------
B. Parameters, dx, and stability checks
-----------------------------------------------------

2) Unused variable 'nsteps'
   Original:
       long nsteps, nx, ndt;
   Only nx and ndt are used; nsteps is never used.

   Problem:
   - With -Wall -Werror, an unused variable is treated as an error.
   Fix:
   - Removed 'nsteps' entirely and only kept:
       long nx, ndt;

3) Wrong spatial step dx
   Original:
       double dx = lx/nx;
   Specification: there are Nx grid points with the first at x = 0 and the
   final at x = L. That means:
       dx = L / (Nx - 1).

   If we use L/Nx, the last grid point is at L - dx instead of exactly L.

   Fix:
       double dx = lx / (nx - 1.0);

   This matches the spec’s geometry: x_j = j * dx, with j = 0,...,nx-1 gives
   x_0 = 0 and x_{nx-1} = lx.

4) Incorrect stability condition (misuse of ||)
   Original:
       if((c*dt/dx||gamma*dt)>1.0) {
           printf("Timestep too large\n");
           return 1;
       }

   Problem:
   - The expression (c*dt/dx || gamma*dt) uses the logical OR operator on
     two doubles. In C, this converts each term to 0 or 1 BEFORE comparison.
   - So the code does NOT implement the intended test
       C*dt/dx > 1  OR  gamma*dt > 1
     as required by the specification.

   Fix:
       if ( (c * dt / dx > 1.0) || (gamma * dt > 1.0) ) {
           printf("Timestep too large\n");
           return 1;
       }

   This correctly enforces the condition:
   - If C dt/dx > 1 or gamma dt > 1, exit with a warning. :contentReference[oaicite:1]{index=1}  

-----------------------------------------------------
C. Dynamic memory allocation
-----------------------------------------------------

5) malloc with wrong element size
   Original:
       y      = malloc(nx*sizeof(nx));
       y_next = malloc(nx*sizeof(nx));
       z      = malloc(nx*sizeof(nx));
       z_next = malloc(nx*sizeof(nx));

   Problem:
   - sizeof(nx) is the size of a long, but the arrays store doubles.
   - This allocates the wrong number of bytes and can cause out-of-bounds
     memory access or wasted memory.

   Fix:
       y      = malloc(nx * sizeof(double));
       y_next = malloc(nx * sizeof(double));
       z      = malloc(nx * sizeof(double));
       z_next = malloc(nx * sizeof(double));

   This matches the stored type (double) and allocates nx doubles for each.

6) Not exiting on allocation failure
   Original:
       if((y==NULL)||(y_next==NULL)||(z==NULL)||(z_next==NULL)) {
           printf("Allocation error\n");
       }

   Problem:
   - If allocation fails, the code prints a message but continues to use
     NULL pointers, which will definitely crash later.

   Fix:
       if ((y == NULL) || (y_next == NULL) || (z == NULL) || (z_next == NULL)) {
           printf("Allocation error\n");
           free(y);
           free(y_next);
           free(z);
           free(z_next);
           return 1;
       }

   Now the program stops safely if memory allocation fails.

-----------------------------------------------------
D. Initialisation / time variable
-----------------------------------------------------

7) Uninitialised time variable 'ctime'
   Original:
       double ctime;
       ...
       while (ctime<tf) {
       ...
           ctime + dt;

   Problems:
   - ctime was never given an initial value before being used in the while
     condition. 
   - The update line 'ctime + dt;' computes a value and throws it away.
     It does NOT update ctime.

   Fix:
       double ctime = 0.0;   // start at t = 0
       ...
       while (ctime < tf) {
           ...
           ctime += dt;       // advances time

   This makes the time loop correct and removes "used uninitialised" warnings.

8) Incorrect initialisation of Y and never initialising Z
   Original:
       double x;
       ...
       // **************
       // initialisation 
       // **************
       y[x] = exp(-x);

   Problems:
   - x is a double and uninitialised at this point.
   - Using y[x] uses a double as an array index (compiler error).
   - Only one element of y was set; the rest of y and all of z were never
     initialised.
   - Spec requires:
       Y(x,0) = exp(-x)
       Z(x,0) = 0   for all grid points. :contentReference[oaicite:2]{index=2}  

   Fix:
       double ctime = 0.0;
       long j;
       for (j = 0; j < nx; j++) {
           double x = j * dx;
           y[j] = exp(-x);    // initial condition for Y
           z[j] = 0.0;        // initial condition for Z
       }

       // Apply boundary conditions at t = 0
       y[0]      = 1.0;      // Y(0,0) = 1
       z[nx - 1] = 0.0;      // Z(L,0) = 0

   This correctly sets the entire initial profile of Y and Z and matches the
   given IC/BC.

9) Missing output at t = 0  

   Original code only printed inside the while loop, so t=0 was never printed.

   Fix:
       // After initialising Y and Z, output the initial state
       for (j = 0; j < nx; j++) {
           double x = j * dx;
           printf("%g %g %g %g\n", ctime, x, y[j], z[j]);
       }

   Then the time loop outputs every ndt timesteps.

-----------------------------------------------------
E. Time stepping and spatial derivatives
-----------------------------------------------------

10) Out-of-bounds indexing in finite differences
    Original:
        for (j=0; j<nx; j++) {
            double yslope = (y[j+1]-y[j-1])/(2*dx);
            ...
            double zslope = (z[j+1]-z[j-1])/(2*dx);
            ...
        }

    Problem:
    - When j = 0: y[j-1] = y[-1], z[j-1] = z[-1] → out of bounds.
    - When j = nx-1: y[j+1] = y[nx], z[j+1] = z[nx] → out of bounds.

    Fix:
    - Only use central differences for interior points, j = 1..nx-2.
    - Handle boundaries separately.

    New approach (see code below):
        for (j = 1; j < nx - 1; j++) {
            double dYdx = (y[j] - y[j-1]) / dx;      // upwind for Y
            double dZdx = (z[j+1] - z[j]) / dx;      // upwind for Z
            ...
        }
        // then explicitly update j=0 and j=nx-1 with BCs / one-sided formulas.

11) Central difference scheme is not stable for pure advection
    Original:
        double yslope = (y[j+1]-y[j-1])/(2*dx);
        double zslope = (z[j+1]-z[j-1])/(2*dx);

    Problem:
    - For an advection equation solved with forward Euler in time, a central
      difference in space is generally unstable, even if C dt/dx <= 1.
    - The spec explicitly mentions “simple difference scheme” and also imposes
      a CFL-like condition, which is consistent with an upwind scheme.

    Fix (upwind scheme):
    - For Y: ∂Y/∂t = -C ∂Y/∂x - γY can be written as:
          ∂Y/∂t + C ∂Y/∂x = -γY
      So Y moves to the right with speed C. We use a backward (upwind) difference:
          dYdx ≈ (Y_j - Y_{j-1}) / dx

    - For Z: ∂Z/∂t = C ∂Z/∂x + γY can be written as:
          ∂Z/∂t - C ∂Z/∂x = γY
      So Z moves to the left with speed C (or velocity -C). We use a forward
      (upwind) difference:
          dZdx ≈ (Z_{j+1} - Z_j) / dx

    Combined with forward Euler:
        Y^{n+1}_j = Y^n_j + dt (-C dYdx - γ Y^n_j)
        Z^{n+1}_j = Z^n_j + dt ( C dZdx + γ Y^n_j)

12) Missing boundary conditions in the time-stepping
    Spec: boundary conditions are
        Y(0,t) = 1
        Z(L,t) = 0   for all t. :contentReference[oaicite:4]{index=4}  

    Original code:
    - Never explicitly enforced these in the time loop.
    - Evolved all points with the same finite difference update.

    Fix:
    - At the new time level we explicitly set:
        y_next[0]      = 1.0;   // Y(0,t) = 1
        z_next[nx - 1] = 0.0;   // Z(L,t) = 0
    - For the other ends (x = L for Y, x = 0 for Z), we use one-sided upwind
      derivatives consistent with the PDE.

-----------------------------------------------------
F. Pointer updates and memory freeing
-----------------------------------------------------

13) Incorrect update of y and z between time steps
    Original:
        // Copy next values at timestep to y array.
        y = y_next;
        z = z_next;

    Problems:
    - This does NOT copy values; it just makes y and y_next point to the same
      memory, and leaks the original y allocation.
    - The comment is misleading: no copy happens.

    Fix:
    - Use pointer swapping to treat y_next as the new "current" array:
        double *tmp;
        tmp    = y;      y = y_next;      y_next = tmp;
        tmp    = z;      z = z_next;      z_next = tmp;

    This is efficient (no element-wise copy) and keeps both allocations valid.

14) free() in the wrong place and missing frees
    Original (brace structure was off, but logically):
        while (ctime<tf){
            ...
        free(y);
        free(y_next);
    }}

    plus z and z_next were never freed.

    Problems:
    - free(y) inside the while loop would deallocate the arrays after the first
      iteration, but the loop still tries to use them.
    - z and z_next leaked.

    Fix:
    - Move all four free() calls to after the time loop:
        free(y);
        free(y_next);
        free(z);
        free(z_next);

    - Ensure they are only called once, at the end of main.

-----------------------------------------------------
G. Other cleanliness / spec-compliance points
-----------------------------------------------------

15) Output format and times
    - Spec requires output at t = 0, dt*nD, 2*dt*nD, ...:
      I added printing of the initial condition at t=0, and in the loop I output
      whenever ntstep % ndt == 0 after the time step is applied.

16) read_input function
    - The assignment explicitly states that the read_input function and its
      prototype are correct and should not be modified, so I left that function
      unchanged and only fixed how it is used (headers, arguments).

---------------------------------------------------------------------------
END OF NOTES
---------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

long read_input(double *Lx_p, long *nx_p, double *c_p,
                double *tf_p, double *dt_p, long *ndt_p,
                double *S_p, char *fname);

int main(void) {
    
    long nx, ndt;            
    double lx, c, tf, dt;    
    double gamma;            

    char *fname = "input.txt";

    // Read parameters from file
    if (read_input(&lx, &nx, &c, &tf, &dt, &ndt, &gamma, fname)) {
        printf("File read error\n");
        return 2;
    }

    // Need at least 2 grid points to have a domain [0, L]
    if (nx < 2) {
        printf("Need at least 2 grid points\n");
        return 1;
    }

    // Spatial step: Nx points from 0 to L inclusive -> dx = L/(Nx-1)
    double dx = lx / (nx - 1.0);

    // Stability checks: C dt/dx <= 1 and gamma dt <= 1
    if ( (c * dt / dx > 1.0) || (gamma * dt > 1.0) ) {
        printf("Timestep too large\n");
        return 1;
    }

    // -----------------
    // Grid storage
    // -----------------
    double *y, *y_next;
    double *z, *z_next;

    y      = malloc(nx * sizeof(double));
    y_next = malloc(nx * sizeof(double));
    z      = malloc(nx * sizeof(double));
    z_next = malloc(nx * sizeof(double));

    if ((y == NULL) || (y_next == NULL) || (z == NULL) || (z_next == NULL)) {
        printf("Allocation error\n");
        free(y);
        free(y_next);
        free(z);
        free(z_next);
        return 1;
    }

    // -----------------
    // Initial conditions
    // -----------------
    double ctime = 0.0;      // current time
    long j;

    for (j = 0; j < nx; j++) {
        double x = j * dx;
        y[j] = exp(-x);      // Y(x,0) = exp(-x)
        z[j] = 0.0;          // Z(x,0) = 0
    }

    // Apply boundary conditions at t = 0
    y[0]        = 1.0;       // Y(0,t) = 1
    z[nx - 1]   = 0.0;       // Z(L,t) = 0

    // Output initial condition (t = 0)
    for (j = 0; j < nx; j++) {
        double x = j * dx;
        printf("%g %g %g %g\n", ctime, x, y[j], z[j]);
    }

    long ntstep = 0;         // number of timesteps taken so far

    // -----------------
    // Time-stepping loop
    // -----------------
    while (ctime < tf) {
        // Enforce boundary conditions at current time
        y[0]      = 1.0;     // left boundary for Y
        z[nx - 1] = 0.0;     // right boundary for Z

        // Update interior points with upwind differences
        for (j = 1; j < nx - 1; j++) {
            // Upwind derivative for Y moving right with speed c
            double dYdx = (y[j] - y[j - 1]) / dx;

            // Upwind derivative for Z moving left (velocity -c)
            double dZdx = (z[j + 1] - z[j]) / dx;

            // PDE right-hand sides:
            // dY/dt = -c dY/dx - gamma Y
            double dYdt = -c * dYdx - gamma * y[j];

            // dZ/dt =  c dZ/dx + gamma Y
            double dZdt =  c * dZdx + gamma * y[j];

            // Forward Euler update
            y_next[j] = y[j] + dt * dYdt;
            z_next[j] = z[j] + dt * dZdt;
        }

        // Boundaries at next time level:

        // Left boundary for Y
        y_next[0] = 1.0;

        // Right boundary for Y
        {
            double dYdx_right = (y[nx - 1] - y[nx - 2]) / dx;
            double dYdt_right = -c * dYdx_right - gamma * y[nx - 1];
            y_next[nx - 1]    = y[nx - 1] + dt * dYdt_right;
        }

        // Left boundary for Z
        {
            double dZdx_left = (z[1] - z[0]) / dx;
            double dZdt_left =  c * dZdx_left + gamma * y[0];
            z_next[0]        = z[0] + dt * dZdt_left;
        }

        // Right boundary for Z
        z_next[nx - 1] = 0.0;

        // Swap current and next arrays
        double *tmp;
        tmp    = y;      y = y_next;      y_next = tmp;
        tmp    = z;      z = z_next;      z_next = tmp;

        // Advance time
        ctime  += dt;
        ntstep += 1;

        // Output every ndt timesteps
        if (ntstep % ndt == 0) {
            for (j = 0; j < nx; j++) {
                double x = j * dx;
                printf("%g %g %g %g\n", ctime, x, y[j], z[j]);
            }
        }
    }


    free(y);
    free(y_next);
    free(z);
    free(z_next);

    return 0;
}


long read_input(double *Lx_p, long *nx_p, double *c_p,
                double *tf_p, double *dt_p, long *ndt_p,
                double *S_p, char *fname) {
    FILE* fptr = fopen(fname, "r");
    if (fptr == NULL) return 1;
    if (7 != fscanf(fptr, "%lf %ld %lf %lf %lf %ld %lf",
                    Lx_p, nx_p, c_p, tf_p, dt_p, ndt_p, S_p)) {
        return 1;
    }
    fclose(fptr);
    return 0;
}

# mcs_rot_angles


Geometrical definitions used by the code. See [van der Marel & Cioni (2001)][1],
and [van der Marel et al. (2002)][2].

![alt tag](/sideview.png)

**Fig 2, van der Marel et al. (2002)**

![alt tag](/systems.png)

**Fig 3, van der Marel et al. (2002)**

* Line of nodes: intersection of the rotated galaxy plane (x', y', z') and
the sky plane (x, y, z). Formally, it is the x' axis.

* Position angle of the line of nodes (\Theta): counterclockwise rotation
around the z axis, measured form the y axis (North).
\Theta = \theta - 90º, ie: \theta is measured from the x axis (West)

* Inclination (i): clockwise rotation around the new x' axis. It is the
angle between the (x, y)-plane of the sky and the
(x', y')-plane of the galaxy disk.


## Structure

### 1. `rot_angles()`

Main function. Reads input data, obtains the angles for each method,
and produces the final plots.

### 2. `readData()`

Read input data for both galaxies.

### 3. `procParams()`

Reads all the parameters needed to process the data.

1. `inc_PA_grid()`

   Define a grid of (i, theta) angles.
2. `angles2Plane()`

   Convert each pair of (i, theta) angles to its "plane form":
   a*x + b*y + c*z + d = 0.
   1. `gal_theta()`

      Convert position angle 'Theta' (N-->E) to 'theta' (W-->E).

### 4. `galax_struct_dist()`

Process all data.

1. `dist2CloudCenter()`

   3D distance between cluster-galaxy center and its error, in parsec.
   1. `MCs_data()`

      Read fixed data for each galaxy (center, distance, error in distance.
   2. `dist_err_mag_2_pc()`

      Convert the error in distance modulus to an error, in parsec.
   3. `dist_err_2_pts()`

      Error for the cluster-galaxy center distances, in parsec. 
2. `dist_filter()`

   Filter clusters based on their projected angular distances ('rho').
   1. `get_rho_phi()`

      Projected angular distance cluster-galaxy center, and position angle
      of cluster.
3. `i_PA_DeprjDist()`

   Deprojected distance values for each (i, PA) point in the defined grid.
   1. `get_deproj_dist()`

      1. `vdm_2001_dep_dist_kpc()`

         Use Eq. (8) from van der Marel & Cioni (2001).
4. `xyz_coords()`

   Coordinates of filtered clusters in the (x,y,z) system.
5. Methods

   1.  m1_ccc_map
       1. ccc
          Concordance correlation coefficient between deprojected distances
          obtained using ASteCA's distances and the van der Marel's equations.
   2. interp_dens_map
   3. m2_fix_plane_perp_dist
      For each (i, theta) pair, calculate the mean of the perpendicular
      distances of all clusters to that plane.
   4. interp_dens_map
   5. m3_min_perp_distance
      1. perp_error
   6. interp_dens_map
   2. get_angles
      Extract the best fit angles from all the (i, theta) values tested by
      the methods 1 and 2. For method 3, call function below. 
      1. angle_betw_planes
         Transform the (a, b, c, d) plane parameters for the best fit found
         by Method 3, into (i, theta) angles.
   3. ccc_sum_d_for_best_fit
      1. get_deproj_dist
      2. ccc
      3. perp_error
   4. monte_carlo_errors
      1. draw_rand_dep_dist
      2. draw_rand_dist_mod
   5. cov_ellipse

## 5. `make_plots()`

Produce final plots.

__________________________________________________

[Subramanian & Subramaniam (2015)][3] (SMC),
[Subramanian & Subramaniam (2010)][4] (LMC)

![alt tag](/S-LMC_i_pa.png)

![alt tag](/inno_2016.png)


__________________________
[1]: http://adsabs.harvard.edu/abs/2001AJ....122.1807V
[2]: http://adsabs.harvard.edu/abs/2002AJ....124.2639V
[3]: http://arxiv.org/abs/1410.7588
[4]: http://adsabs.harvard.edu/abs/2010A%26A...520A..24S
#In this file we define a sinngle ndt cell message
float64 mean_x       #x coordinate of mean
float64 mean_y       #y coordinate of mean
float64 mean_z       #z coordinate of mean
float64[] cov_matrix #covariance matrix 3x3 Parameters in order xx,xy,xz,yx,yy,yz,zx,zy,zz
float64 occupancy    #occupancy
int64 N              #point count
/* ggrd_velinterpol.c */
int ggrd_find_vel_and_der(double *, double, double, struct ggrd_master *, int, unsigned short, unsigned short, double *, double *, double *);
void ggrd_get_velocities(double *, double *, double *, int, struct ggrd_master *, double, double);
void ggrd_weights(double, double *, int, int, double [(5 +1)][(1 +1)]);
/* ggrd_readgrds.c */
void ggrd_init_vstruc(struct ggrd_master *);
int ggrd_read_vel_grids(struct ggrd_master *, double, unsigned short, unsigned short, char *, unsigned char);
void ggrd_resort_and_check(double *, float *, double *, int, int, unsigned short, double, unsigned short, unsigned short, double, unsigned char *);
void ggrd_read_depth_levels(struct ggrd_master *, int **, char *, unsigned short);
/* ggrd_grdtrack_util.c */
void ggrd_init_master(struct ggrd_master *);
void ggrd_grdinfo(char *);
int ggrd_grdtrack_init_general(unsigned char, char *, char *, char *, struct ggrd_gt *, unsigned char, unsigned char, unsigned char);
int ggrd_grdtrack_rescale(struct ggrd_gt *, unsigned char, unsigned char, unsigned char, double);
unsigned char ggrd_grdtrack_interpolate_rtp(double, double, double, struct ggrd_gt *, double *, unsigned char, unsigned char, double);
unsigned char ggrd_grdtrack_interpolate_lonlatz(double, double, double, struct ggrd_gt *, double *, unsigned char);
unsigned char ggrd_grdtrack_interpolate_xyz(double, double, double, struct ggrd_gt *, double *, unsigned char);
unsigned char ggrd_grdtrack_interpolate_tp(double, double, struct ggrd_gt *, double *, unsigned char, unsigned char);
unsigned char ggrd_grdtrack_interpolate_xy(double, double, struct ggrd_gt *, double *, unsigned char);
void ggrd_grdtrack_free_gstruc(struct ggrd_gt *);
void ggrd_find_spherical_vel_from_rigid_cart_rot(double *, double *, double *, double *, double *);
int ggrd_grdtrack_init(double *, double *, double *, double *, float **, int *, char *, struct GRD_HEADER **, struct GMT_EDGEINFO **, char *, unsigned char *, GMT_LONG *, unsigned char, char *, float **, int *, GMT_LONG, unsigned char, unsigned char, struct GMT_BCR *);
void ggrd_print_layer_avg(float *, float *, int, int, int, FILE *, GMT_LONG *);
unsigned char ggrd_grdtrack_interpolate(double *, unsigned char, struct GRD_HEADER *, float *, struct GMT_EDGEINFO *, int, float *, int, double *, unsigned char, struct GMT_BCR *);
int ggrd_init_thist_from_file(struct ggrd_t *, char *, unsigned char, unsigned char);
void ggrd_gt_interpolate_z(double, float *, int, int *, int *, double *, double *, unsigned char, unsigned char *);
int ggrd_interpol_time(double, struct ggrd_t *, int *, int *, double *, double *);
int interpolate_seafloor_ages(double, double, double, struct ggrd_master *, double *);
FILE *ggrd_open(char *, char *, char *);
void ggrd_vecalloc(double **, int, char *);
void ggrd_vecrealloc(double **, int, char *);
void ggrd_calc_mean_and_stddev(double *, double *, int, double *, double *, double *, unsigned char, unsigned char, double *);
void ggrd_indexx(int, double *, int *);
float ggrd_gt_rms(float *, int);
float ggrd_gt_mean(float *, int);
/* sh_exp_ggrd.c */
void sh_read_spatial_data_from_grd(struct sh_lms *, struct ggrd_gt *, unsigned short, int, double *, double *);

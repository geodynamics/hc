/* sh_exp.c */
void sh_allocate_and_init(struct sh_lms **, int, int, int, int, unsigned short, unsigned short);
void sh_init_expansion(struct sh_lms *, int, int, int, unsigned short, unsigned short);
void sh_free_expansion(struct sh_lms *, int);
void sh_clear_alm(struct sh_lms *);
long double sh_total_power(struct sh_lms *);
long double sh_total_rms(struct sh_lms *);
void sh_compute_power_per_degree(struct sh_lms *, long double *);
long double sh_correlation(struct sh_lms *, struct sh_lms *, int);
long double sh_correlation_per_degree(struct sh_lms *, struct sh_lms *, int, int);
void sh_single_par_and_exp_to_file(struct sh_lms *, char *, unsigned short, unsigned short);
void sh_single_par_and_exp_to_stream(struct sh_lms *, FILE *, unsigned short, unsigned short);
void sh_print_parameters_to_stream(struct sh_lms *, int, int, int, long double, FILE *, unsigned short, unsigned short, unsigned short);
unsigned short sh_read_parameters_from_stream(int *, int *, int *, int *, int *, long double *, int *, FILE *, unsigned short, unsigned short, unsigned short);
void sh_print_coefficients_to_stream(struct sh_lms *, int, FILE *, long double *, unsigned short, unsigned short);
void sh_read_coefficients_from_stream(struct sh_lms *, int, int, FILE *, unsigned short, long double *, unsigned short);
void sh_print_nonzero_coeff(struct sh_lms *, FILE *);
void sh_read_spatial_data_from_stream(struct sh_lms *, FILE *, unsigned short, int, long double *, long double *);
void sh_read_spatial_data(struct sh_lms *, FILE *, unsigned short, int, long double *, long double *);
void sh_compute_spatial_basis(struct sh_lms *, FILE *, unsigned short, long double, long double **, int, unsigned short);
void sh_compute_spectral(long double *, int, unsigned short, long double **, struct sh_lms *, unsigned short);
void sh_compute_spatial(struct sh_lms *, int, unsigned short, long double **, long double *, unsigned short);
void sh_compute_spatial_reg(struct sh_lms *, int, unsigned short, long double **, long double *, int, long double *, int, long double *, unsigned short, unsigned short);
void sh_compute_spatial_irreg(struct sh_lms *, int, long double *, long double *, int, long double *, unsigned short);
void sh_exp_type_error(char *, struct sh_lms *);
void sh_print_plm(long double *, int, int, int, FILE *);
void sh_print_spatial_data_to_stream(struct sh_lms *, int, long double *, unsigned short, long double, FILE *);
void sh_get_coordinates(struct sh_lms *, int, long double *, long double *);
void sh_print_reg_spatial_data_to_stream(struct sh_lms *, int, long double *, unsigned short, long double, long double *, int, long double *, int, FILE *);
void sh_print_irreg_spatial_data_to_stream(struct sh_lms *, int, long double *, unsigned short, long double, long double *, long double *, int, FILE *);
void sh_compute_plm(struct sh_lms *, int, long double **, unsigned short);
void sh_compute_plm_reg(struct sh_lms *, int, long double **, unsigned short, long double *, int);
void sh_get_coeff(struct sh_lms *, int, int, int, unsigned short, long double *);
void sh_write_coeff(struct sh_lms *, int, int, int, unsigned short, long double *);
void sh_add_coeff(struct sh_lms *, int, int, int, unsigned short, long double *);
void sh_copy_lms(struct sh_lms *, struct sh_lms *);
void sh_aexp_equals_bexp_coeff(struct sh_lms *, struct sh_lms *);
void sh_c_is_a_plus_b_coeff(struct sh_lms *, struct sh_lms *, struct sh_lms *);
void sh_scale_expansion_l_factor(struct sh_lms *, long double *);
void sh_scale_expansion(struct sh_lms *, long double);
/* sh_model.c */
void sh_init_model(struct sh_lms_model *, int, int, int, int, int, int, unsigned short);
void sh_free_model(struct sh_lms_model *);
void sh_print_model_coefficients(struct sh_lms_model *, FILE *, unsigned short, unsigned short);
void sh_print_model_spatial_basis(struct sh_lms_model *, FILE *, unsigned short);
void sh_read_model_spatial_data(struct sh_lms_model *, long double **, FILE *, unsigned short);
void sh_compute_model_spectral(struct sh_lms_model *, long double *, unsigned short);
void sh_compute_model_spatial(struct sh_lms_model *, long double **, unsigned short);
void sh_print_model_spatial_data(struct sh_lms_model *, long double *, FILE *, unsigned short);
/* hc_init.c */
void hc_init_parameters(struct hc_parameters *);
void hc_struc_init(struct hcs **);
void hc_init_polsol_struct(struct hc_ps *);
void hc_init_main(struct hcs *, int, struct hc_parameters *);
void hc_init_constants(struct hcs *, long double, char *, unsigned short);
void hc_handle_command_line(int, char **, int, struct hc_parameters *);
void hc_assign_viscosity(struct hcs *, int, long double [4], struct hc_parameters *);
void hc_assign_density(struct hcs *, struct hc_parameters *, int, char *, int, unsigned short, unsigned short, unsigned short);
long double hc_find_dens_scale(long double, long double, unsigned short, long double *, long double *, int);
void hc_init_phase_boundaries(struct hcs *, int, unsigned short);
void hc_assign_plate_velocities(struct hcs *, int, char *, unsigned short, int, unsigned short, unsigned short, unsigned short, unsigned short);
void hc_init_single_plate_exp(char *, struct hcs *, unsigned short, struct sh_lms *, unsigned short, unsigned short, unsigned short, unsigned short);
void hc_init_l_factors(struct hcs *, int);
void hc_get_blank_expansions(struct sh_lms **, int, int, char *);
void hc_struc_free(struct hcs **);
void hc_assign_dd_scaling(int, long double [4], struct hc_parameters *, long double);
void hc_read_scalar_shexp(char *, struct sh_lms **, char *, struct hc_parameters *);
void hc_select_pvel(long double, struct pvels *, struct sh_lms *, unsigned short);
FILE *hc_fopen(char *, char *, char *, char *);
/* hc_solve.c */
void hc_solve(struct hcs *, unsigned short, int, struct sh_lms *, unsigned short, unsigned short, unsigned short, unsigned short, unsigned short, struct sh_lms *, struct sh_lms *, struct sh_lms *, unsigned short, unsigned short);
void hc_sum(struct hcs *, int, struct sh_lms *, struct sh_lms *, int, unsigned short, struct sh_lms *, unsigned short);
void hc_compute_sol_spatial(struct hcs *, struct sh_lms *, long double **, unsigned short);
void hc_compute_dynamic_topography(struct hcs *, struct sh_lms *, struct sh_lms **, unsigned short, unsigned short);
void hc_calc_geoid_corr_four_layer(long double *, struct sh_lms *, struct sh_lms *, struct sh_lms *, struct hc_parameters *, struct hcs *, unsigned short *, long double *, long double *);
/* hc_propagator.c */
void hc_evalpa(int, long double, long double, long double, long double *);
void hc_evppot(int, long double, long double *);
/* hc_polsol.c */
void hc_polsol(struct hcs *, int, long double *, int, long double *, unsigned short, struct sh_lms *, unsigned short, int, long double *, long double *, unsigned short, struct sh_lms *, struct sh_lms *, unsigned short, struct sh_lms *, unsigned short, unsigned short, unsigned short);
/* hc_matrix.c */
void hc_ludcmp_3x3(long double [3][3], int, int *);
void hc_lubksb_3x3(long double [3][3], int, int *, long double *);
/* hc_torsol.c */
void hc_torsol(struct hcs *, int, int, int, long double *, long double **, long double **, struct sh_lms *, struct sh_lms *, long double *, unsigned short);
/* hc_output.c */
void hc_print_spectral_solution(struct hcs *, struct sh_lms *, FILE *, int, unsigned short, unsigned short);
void hc_print_sh_scalar_field(struct sh_lms *, FILE *, unsigned short, unsigned short, unsigned short);
void hc_print_spatial_solution(struct hcs *, struct sh_lms *, long double *, char *, char *, int, unsigned short, unsigned short);
void hc_print_depth_layers(struct hcs *, FILE *, unsigned short);
void hc_print_3x3(long double [3][3], FILE *);
void hc_print_sm(long double [6][4], FILE *);
void hc_print_vector(long double *, int, FILE *);
void hc_print_vector_label(long double *, int, FILE *, char *);
void hc_print_matrix_label(long double *, int, int, FILE *, char *);
void hc_print_vector_row(long double *, int, FILE *);
void hc_compute_solution_scaling_factors(struct hcs *, int, long double, long double, long double *);
void hc_print_poloidal_solution(struct sh_lms *, struct hcs *, int, char *, unsigned short, unsigned short);
void hc_print_toroidal_solution(long double *, int, struct hcs *, int, char *, unsigned short);
void hc_print_vtk(FILE *, long double *, long double *, int, int, unsigned short, int, long double *, int, int);
int hc_print_be_float(long double *, int, FILE *, unsigned short);
int hc_print_float(long double *, int, FILE *);
int hc_read_float(long double *, int, FILE *);
void hc_print_be_int(int *, int, FILE *, unsigned short);
unsigned short hc_is_little_endian(void);
void hc_flip_byte_order(void *, size_t);
void hc_flipit(void *, void *, size_t);
void hc_print_dens_anom(struct hcs *, FILE *, unsigned short, unsigned short);
void hc_print_geoid_kernel(struct sh_lms *, long double *, int, FILE *, unsigned short);
/* hc_input.c */
int hc_read_sh_solution(struct hcs *, struct sh_lms **, FILE *, unsigned short, unsigned short);
/* hc_misc.c */
void hc_hvecalloc(long double **, int, char *);
void hc_dvecalloc(double **, int, char *);
void hc_svecalloc(float **, int, char *);
void hc_ivecalloc(int **, int, char *);
void hc_vecalloc(long double **, int, char *);
void hc_scmplx_vecalloc(struct hc_scmplx **, int, char *);
void hc_svecrealloc(float **, int, char *);
void hc_dvecrealloc(long double **, int, char *);
void hc_vecrealloc(long double **, int, char *);
float hc_vec_rms_diff(long double *, long double *, int);
float hc_vec_rms(long double *, int);
void hc_a_equals_b_svector(float *, float *, int);
void hc_a_equals_b_vector(long double *, long double *, int);
float hc_mean_svec(float *, int);
long double hc_mean_vec(long double *, int);
void hc_zero_dvector(long double *, int);
void hc_zero_lvector(unsigned short *, int);
void hc_get_flt_frmt_string(char *, int, unsigned short);
char *hc_name_boolean(unsigned short);
unsigned short hc_toggle_boolean(unsigned short *);
void hc_advance_argument(int *, int, char **);
void hc_compute_correlation(struct sh_lms *, struct sh_lms *, long double *, int, unsigned short);
void lonlatpv2cv(long double, float, long double *, long double *);
void thetaphipv2cv(long double, float, long double *, long double *);
void lonlatpv2cv_with_base(long double *, long double *, long double *);
void calc_polar_base_at_theta_phi(long double, long double, long double *);
void hc_linear_interpolate(long double *, int, long double, int *, int *, long double *, long double *);
/* hc_extract_sh_layer.c */
/* hc_extract_spatial.c */
/* hc_visc_scan.c */
void visc_scan_out(long double *, struct sh_lms *, struct sh_lms *, struct sh_lms *, struct hc_parameters *, struct hcs *, unsigned short *, unsigned short);
/* rick_sh_c.c */
void rick_compute_allplm(int, int, long double *, long double *, struct rick_module *);
void rick_compute_allplm_reg(int, int, long double *, long double *, struct rick_module *, long double *, int);
void rick_pix2ang(int, int, long double *, long double *, struct rick_module *);
void rick_shc2d(long double *, long double *, int, int, long double *, long double *, struct rick_module *);
void rick_shc2d_reg(long double *, long double *, int, int, long double *, long double *, struct rick_module *, long double *, int, long double *, int, unsigned short);
void rick_shc2d_pre(long double *, long double *, int, long double *, long double *, int, long double *, long double *, struct rick_module *);
void rick_shc2d_pre_reg(long double *, long double *, int, long double *, long double *, int, long double *, long double *, struct rick_module *, long double *, int, long double *, int, unsigned short);
void rick_shc2d_irreg(long double *, long double *, int, int, long double *, long double *, struct rick_module *, long double *, long double *, int);
void rick_shd2c(long double *, long double *, int, int, long double *, long double *, struct rick_module *);
void rick_shd2c_pre(long double *, long double *, int, long double *, long double *, int, long double *, long double *, struct rick_module *);
void rick_init(int, int, int *, int *, int *, struct rick_module *, unsigned short);
void rick_free_module(struct rick_module *, int);
void rick_plmbar1(long double *, long double *, int, int, long double, struct rick_module *);
void rick_gauleg(long double, long double, long double *, long double *, int);
/* rick_fft_c.c */
void rick_cs2ab(long double *, int);
void rick_ab2cs(long double *, int);
void rick_realft_nr(long double *, int, int);
void rick_four1_nr(long double *, int, int);
/* gaussp.c */
/* print_gauss_lat.c */
/* rotvec2vel.c */
FILE *rv_myopen(const char *, const char *);
/* sh_corr.c */
/* simple_test.c */
/* spherepack_sh.c */
/* test_fft.c */
/* sh_test.c */
/* sh_power.c */
/* sh_model.c */
void sh_init_model(struct sh_lms_model *, int, int, int, int, int, int, unsigned short);
void sh_free_model(struct sh_lms_model *);
void sh_print_model_coefficients(struct sh_lms_model *, FILE *, unsigned short, unsigned short);
void sh_print_model_spatial_basis(struct sh_lms_model *, FILE *, unsigned short);
void sh_read_model_spatial_data(struct sh_lms_model *, long double **, FILE *, unsigned short);
void sh_compute_model_spectral(struct sh_lms_model *, long double *, unsigned short);
void sh_compute_model_spatial(struct sh_lms_model *, long double **, unsigned short);
void sh_print_model_spatial_data(struct sh_lms_model *, long double *, FILE *, unsigned short);
/* shana_sh.c */
void shana_compute_allplm(int, int, double *, double *, struct shana_module *);
void shana_pix2ang(int, int, double *, double *, struct shana_module *);
void shana_shc2d(long double *, long double *, int, int, long double *, long double *, struct shana_module *);
void shana_shc2d_pre(long double *, long double *, int, double *, double *, int, float *, float *, struct shana_module *);
void shana_shd2c(long double *, long double *, int, int, long double *, long double *, struct shana_module *);
void shana_shd2c_pre(long double *, long double *, int, double *, double *, int, long double *, long double *, struct shana_module *);
void shana_init(int, int, int *, int *, int *, struct shana_module *);
void shana_free_module(struct shana_module *, int);
void shana_plmbar1(double *, double *, int, int, long double, struct shana_module *);

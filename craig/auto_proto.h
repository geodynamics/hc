/* rick_sh2.c */
void rick_compute_allplm(int, int, double *, double *, struct rick_module *);
void rick_PIx2ang(int, int, double *, double *, struct rick_module *);
void rick_shc2d(double *, double *, int, int, double *, double *, struct rick_module *);
void rick_shc2d_pre(double, double, int, double *, double *, int, float *, float *, struct rick_module *);
void rick_shd2c(double *, double *, int, int, double *, double *, struct rick_module *);
int rick_shd2c_pre(double *, double *, int, double *, double *, int, double *, double *, struct rick_module *);
void rick_init(int, int, int *, int *, int *, struct rick_module *);
void rick_free_module(struct rick_module *, int);
void rick_compute_allplm(int, int, double *, double *, struct rick_module *);
void rick_plmbar1(double *, double *, int, int, double, struct rick_module *);

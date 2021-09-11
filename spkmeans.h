#ifndef SPKMEANS_H_
#define SPKMEANS_H_
#include <stdio.h>

typedef struct COOR
{
    double value;
    struct COOR *next;
} COOR;

typedef struct VEC
{
    COOR *first_coor;
    struct VEC *next;
} VEC;

typedef struct EIGEN
{
    double *vector;
    double value;
} EIGEN;

/*########## inputs treatment ##########*/
VEC *initialize(FILE *fp, int *n, int *d);
double **initialize_data(VEC *head, int n, int d);
/*########## arguments treatment ##########*/
int is_non_negative_integer(char number[]);
int enum_goal(char *goal);
/*########## SPK ##########*/
int determine_k(EIGEN *eigens, int n);
int my_comparator(const void *a, const void *b);
double **gen_T(EIGEN *eigens, int n, int k);
double **initialize_centroids(double **data, int k, int d);
int *initialize_sizes(int k);
double **initialize_sums(int k, int d);
void assign_vec_to_cluster(double **data, int *sizes, double **sums, double **centroids, int n, int k, int d);
int update_clusters(int *sizes, double **sums, double **centroids, int k, int d);
double square_vec(double *vec, int d);
double *subtract_vectors(double *vec1, double *vec2, int d);
void add_vectors(double *vec_sum, double *vec_data, int d);
/*########## WAM ##########*/
double **gen_wam(double **points, int n, int d);
double calc_weight(double *a, double *b, int d);
/*########## DDG ##########*/
double **gen_ddg(double **wam, int n);
double calc_row_sum(double *wam_row, int n);
/*########## LNORM ##########*/
double **gen_lnorm(double **wam, double **ddg, int n);
void mat_sub(double **A, double **B, int n);
double **mat_mul(double **A, double **B, int n);
/*########## JACOBI ##########*/
EIGEN *jacobi(double **A, int n);
void find_pivot(double **A, int *i_tmp, int *j_tmp, int n);
void find_s_c(double **A, int i_tmp, int j_tmp, double *s, double *c);
void update_A(double **A, int i_tmp, int j_tmp, double s, double c, int n);
void update_V(double **V, int i_tmp, int j_tmp, double s, double c, int n);
void V_to_eigens(double **V, double **A, EIGEN *eigens, int n);
int my_comparator(const void *a, const void *b);
int converge(double **A, double *old, int n); /* 0 -> converged (stop) 1 -> didn't converged (continue) */
double calc_off(double **A, int n);
/*########## output ##########*/
void print_matrix(double **mat, int lines, int columns);

#endif
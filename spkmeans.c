#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "spkmeans.h"

/*########## initialize data ##########*/
VEC *initialize(FILE *fp, int *n, int *d)
{
    int fsf, num_of_vec = 0, dimension = 0;
    double value;
    char c;
    VEC *head = (VEC *)malloc(sizeof(VEC));
    VEC *last_vec = head;
    COOR *last_coor = (COOR *)malloc(sizeof(COOR));
    assert(head != NULL);
    assert(last_coor != NULL);
    last_vec->first_coor = last_coor;
    last_vec->next = NULL;

    /* run first iteration to capture the dimension d */
    fsf = fscanf(fp, "%lf%c", &value, &c);
    if (fsf == 2)
    {
        while (c != '\n' && fsf != 1)
        {

            last_coor->value = value;
            last_coor->next = NULL;
            last_coor->next = (COOR *)malloc(sizeof(COOR));
            assert(last_coor->next != NULL);
            last_coor = last_coor->next;

            dimension += 1;
            fsf = fscanf(fp, "%lf%c", &value, &c);
        }
    }
    else
    {
        printf("Invalid Input!\n");
        exit(0);
    }
    num_of_vec += 1;
    last_coor->value = value;
    last_coor->next = NULL;
    last_coor = last_coor->next;
    *d = dimension + 1;

    last_vec->next = (VEC *)malloc(sizeof(VEC));
    assert(last_vec->next != NULL);
    last_vec = last_vec->next;
    last_vec->first_coor = NULL;
    last_vec->next = NULL;

    last_vec->first_coor = (COOR *)malloc(sizeof(COOR));
    assert(last_vec->first_coor != NULL);
    last_coor = last_vec->first_coor;

    /* run from iteration 2 */
    fsf = fscanf(fp, "%lf%c", &value, &c);
    while (fsf != EOF)
    {
        last_coor->value = value;
        last_coor->next = NULL;
        last_coor->next = (COOR *)malloc(sizeof(COOR));
        assert(last_coor->next != NULL);
        last_coor = last_coor->next;
        if (c == '\n' || fsf == 1)
        {
            num_of_vec += 1;

            last_vec->next = (VEC *)malloc(sizeof(VEC));
            assert(last_vec->next != NULL);
            last_vec = last_vec->next;
            last_vec->first_coor = (COOR *)malloc(sizeof(COOR));
            assert(last_vec->first_coor != NULL);
            last_vec->next = NULL;
            last_coor = last_vec->first_coor;
        }
        fsf = fscanf(fp, "%lf%c", &value, &c);
    }
    *n = num_of_vec;
    return head;
}

double **initialize_data(VEC *head, int n, int d)
{
    VEC *last_vec = head;
    COOR *last_coor = last_vec->first_coor;
    VEC *tmp_vec;
    COOR *tmp_coor;
    int i = 0;
    int j = 0;
    double **data = (double **)malloc(n*sizeof(double *));
    double *data_inner = (double *)malloc(n * d * sizeof(double));
    assert(data != NULL);
    assert(data_inner != NULL);

    for (i = 0; i < n; ++i)
    {
        data[i] = data_inner + i * d;
    }
    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < d; ++j)
        {
            data[i][j] = last_coor->value;
            tmp_coor = last_coor;
            last_coor = last_coor->next;
            free(tmp_coor);
        }
        tmp_vec = last_vec;
        last_vec = last_vec->next;
        free(tmp_vec);
        last_coor = last_vec->first_coor;
    }
    return data;
}

/*########## arguments treatment ##########*/
int is_non_negative_integer(char number[])
{
    int i;
    if (number[0] == '-')
    {
        return 0;
    }
    for (i = 0; number[i] != '\0'; i++)
    {
        if (number[i] > '9' || number[i] < '0')
        {
            return 0;
        }
    }
    return 1;
}

int enum_goal(char *arg)
{
    int i, goal = 0;
    char *goals[] = {"spk", "wam", "ddg", "lnorm", "jacobi"};
    for (i = 1; i < 6; i++)
    {
        if (strcmp(goals[i - 1], arg) == 0)
        {
            goal = i;
        }
    }
    return goal;
}

/*########## SPK ##########*/
int determine_k(EIGEN *eigens, int n)
{
    int i, m, arg_max = 1;
    double tmp, delta_max = eigens[1].value - eigens[0].value;
    m = (n / 2) - 1;
    for (i = 1; i <= m; i++)
    {
        tmp = eigens[i + 1].value - eigens[i].value;
        if (tmp > delta_max)
        {
            delta_max = tmp;
            arg_max = i;
        }
    }
    return arg_max + 1;
}

int my_comparator(const void *a, const void *b)
{
    const EIGEN *pa = (const EIGEN *)a;
    const EIGEN *pb = (const EIGEN *)b;
    double da = pa->value;
    double db = pb->value;

    return (da > db) ? 1 : ((da < db) ? -1 : 0);
}

double **gen_T(EIGEN *eigens, int n, int k)
{
    int i, j;
    double denominator = 0;
    double **U = (double **)malloc(n * sizeof(double *));
    double *U_inner = (double *)malloc(n * k * sizeof(double));
    double **T = (double **)malloc(n * sizeof(double *));
    double *T_inner = (double *)malloc(n * k * sizeof(double));
    assert(U != NULL);
    assert(U_inner != NULL);
    assert(T != NULL);
    assert(T_inner != NULL);
    for (i = 0; i < n; i++)
    {
        U[i] = U_inner + i * k;
        T[i] = T_inner + i * k;
    }
    /* generating U */
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < n; j++)
        {
            U[j][i] = eigens[i].vector[j];
        }
    }
    /* generating T */
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < k; j++)
        {
            denominator += U[i][j] * U[i][j];
        }
        denominator = sqrt(denominator);

        for (j = 0; j < k; j++)
        {
            T[i][j] = U[i][j] / denominator;
        }
        denominator = 0;
    }
    free(U[0]);
    free(U);
    return T;
}

double **initialize_centroids(double **data, int k, int d)
{
    int i = 0;
    int j = 0;
    double **centroids = (double **)malloc(k * sizeof(double *));
    double *centroids_inner = (double *)malloc(k * d * sizeof(double));
    assert(centroids != NULL);
    assert(centroids_inner != NULL);

    for (i = 0; i < k; ++i)
    {
        centroids[i] = centroids_inner + i * d;
    }
    for (i = 0; i < k; ++i)
    {
        for (j = 0; j < d; ++j)
        {
            centroids[i][j] = data[i][j];
        }
    }
    return centroids;
}

int *initialize_sizes(int k)
{
    int *sizes = (int *)calloc(k, sizeof(int));
    assert(sizes != NULL);
    /* intializes automatically to array of zeros */
    return sizes;
}

double **initialize_sums(int k, int d)
{
    int i;
    double **sums = (double **)malloc(k * sizeof(double *));
    double *sums_inner = (double *)calloc(k * d, sizeof(double));
    assert(sums != NULL);
    assert(sums_inner != NULL);
    for (i = 0; i < k; ++i)
    {
        sums[i] = sums_inner + i * d;
    }
    return sums;
}

void assign_vec_to_cluster(double **data, int *sizes, double **sums, double **centroids, int n, int k, int d)
{
    double argmin;
    int cluster_index;
    double a;
    double *diff;
    int i;
    int j;
    for (i = 0; i < n; ++i)
    {
        argmin = -1.0;
        cluster_index = -1;
        for (j = 0; j < k; ++j)
        {
            diff = subtract_vectors(data[i], centroids[j], d);
            a = square_vec(diff, d);
            free(diff);
            if (a < argmin || argmin == -1.0)
            {
                argmin = a;
                cluster_index = j;
            }
        }
        sizes[cluster_index] += 1;
        add_vectors(sums[cluster_index], data[i], d);
    }
}

int update_clusters(int *sizes, double **sums, double **centroids, int k, int d)
{
    int changed = 0;
    double new_coor;
    int i;
    int j;
    for (i = 0; i < k; ++i)
    {
        for (j = 0; j < d; ++j)
        {
            new_coor = (sums[i][j] / sizes[i]);
            if (new_coor != centroids[i][j])
            {
                changed = 1;
                centroids[i][j] = new_coor;
            }
            sums[i][j] = 0;
        }
        sizes[i] = 0;
    }
    return changed;
}

double square_vec(double *vec, int d)
{
    double a = 0;
    int i;
    for (i = 0; i < d; ++i)
    {
        a += (vec[i] * vec[i]);
    }
    return a;
}

double *subtract_vectors(double *vec1, double *vec2, int d)
{
    int i;
    double *res = (double *)malloc(d * sizeof(double));
    assert(res != NULL);
    for (i = 0; i < d; ++i)
    {
        res[i] = vec1[i] - vec2[i];
    }
    return res;
}

void add_vectors(double *vec_sum, double *vec_data, int d)
{
    int i;
    for (i = 0; i < d; ++i)
    {
        vec_sum[i] += vec_data[i];
    }
}

/*########## WAM ##########*/
double **gen_wam(double **points, int n, int d)
{
    int i, j;
    double weight;
    double **wam = (double **)malloc(n * sizeof(double *));
    double *wam_inner = (double *)calloc(n * n, sizeof(double));
    assert(wam != NULL);
    assert(wam_inner != NULL);
    for (i = 0; i < n; i++)
    {
        wam[i] = wam_inner + i * n;
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < i; j++)
        {
            weight = calc_weight(points[i], points[j], d);
            wam[i][j] = weight;
            wam[j][i] = weight;
        }
    }
    return wam;
}

double calc_weight(double *a, double *b, int d)
{
    int i;
    double weight, sum = 0;
    for (i = 0; i < d; i++)
    {
        sum += pow((a[i] - b[i]), 2);
    }
    weight = exp(-((sqrt(sum)) / 2));
    return weight;
}

/*########## DDG ##########*/
double **gen_ddg(double **wam, int n)
{
    int i;
    double **ddg = (double **)malloc(n * sizeof(double *));
    double *ddg_inner = (double *)calloc(n * n, sizeof(double));
    assert(ddg != NULL);
    assert(ddg_inner != NULL);
    for (i = 0; i < n; i++)
    {
        ddg[i] = ddg_inner + i * n;
    }
    for (i = 0; i < n; i++)
    {
        ddg[i][i] = calc_row_sum(wam[i], n);
    }
    return ddg;
}

double calc_row_sum(double *wam_row, int n)
{
    int i;
    double ret = 0;
    for (i = 0; i < n; i++)
    {
        ret += wam_row[i];
    }
    return ret;
}

/*########## LNORM ##########*/
double **gen_lnorm(double **wam, double **ddg, int n)
{
    int i;
    double **A, **B;
    double **lnorm = (double **)malloc(n * sizeof(double *));
    double *lnorm_inner = (double *)calloc(n * n, sizeof(double));
    assert(lnorm != NULL);
    assert(lnorm_inner != NULL);
    for (i = 0; i < n; i++)
    {
        lnorm[i] = lnorm_inner + i * n;
        lnorm[i][i] = 1;
        ddg[i][i] = 1/sqrt(ddg[i][i]);
    }
    /*
    In the first step, lnorm plays the role of the identity matrix.
    Notice that "mat_sub" changes lnorm in-place.
    */
    A = mat_mul(ddg, wam, n);
    B = mat_mul(A, ddg, n);
    mat_sub(lnorm, B, n);
    free(A);
    free(B);
    return lnorm;
}

void mat_sub(double **A, double **B, int n)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            A[i][j] -= B[i][j];
        }
    }
}

double **mat_mul(double **A, double **B, int n)
{
    int i, j, k;
    double **ret = (double **)malloc(n * sizeof(double *));
    double *ret_inner = (double *)malloc(n * n * sizeof(double));
    assert(ret != NULL);
    assert(ret_inner != NULL);
    for (i = 0; i < n; i++)
    {
        ret[i] = ret_inner + i * n;
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            ret[i][j] = 0;
            for (k = 0; k < n; k++)
            {
                ret[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return ret;
}

/*########## JACOBI ##########*/
EIGEN *jacobi(double **A, int n)
{
    int i, i_tmp, j_tmp, first = 1, cnt = 0;
    double s, c, old;
    EIGEN *eigens = (EIGEN *)malloc(n * sizeof(EIGEN));
    double **V = (double **)malloc(n * sizeof(double *));
    double *V_inner = (double *)calloc(n * n, sizeof(double));
    assert(eigens != NULL);
    assert(V != NULL);
    assert(V_inner != NULL);
    for (i = 0; i < n; i++)
    {
        V[i] = V_inner + i * n;
        V[i][i] = 1;
    }
    do
    {
        cnt += 1;
        find_pivot(A, &i_tmp, &j_tmp, n);
        find_s_c(A, i_tmp, j_tmp, &s, &c);
        if (first)
        {
            first = 0;
            V[i_tmp][i_tmp] = c;
            V[i_tmp][j_tmp] = s;
            V[j_tmp][i_tmp] = -s;
            V[j_tmp][j_tmp] = c;
            if ((old = calc_off(A, n)) < 1.0e-15)
            {
                break;
            }
            update_A(A, i_tmp, j_tmp, s, c, n); /*in place*/
            continue;
        }
        update_A(A, i_tmp, j_tmp, s, c, n); /*in place*/
        update_V(V, i_tmp, j_tmp, s, c, n); /*in place*/
    } while (cnt != 100 && converge(A, &old, n));
    V_to_eigens(V, A, eigens, n);

    free(V[0]);
    free(V);

    return eigens;
}

void find_pivot(double **A, int *i_tmp, int *j_tmp, int n)
{
    int i, j;
    double pivot_value = 0;
    for (i = 0; i < n; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            if (fabs(A[i][j]) > pivot_value)
            {
                pivot_value = fabs(A[i][j]);
                *i_tmp = i;
                *j_tmp = j;
            }
        }
    }
}

void find_s_c(double **A, int i_tmp, int j_tmp, double *s, double *c)
{
    int sign_theta;
    double t;
    double theta = (A[j_tmp][j_tmp] - A[i_tmp][i_tmp]) / (2 * A[i_tmp][j_tmp]);
    if (theta < 0)
    {
        sign_theta = -1;
    }
    else
    {
        sign_theta = 1;
    }
    t = sign_theta / (fabs(theta) + sqrt(pow(theta, 2) + 1));
    *c = 1 / sqrt(pow(t, 2) + 1);
    *s = t * (*c);
}

void update_A(double **A, int i_tmp, int j_tmp, double s, double c, int n)
{
    int r;
    double A_ii = A[i_tmp][i_tmp], A_jj = A[j_tmp][j_tmp];
    for (r = 0; r < n; r++)
    {
        if (r == i_tmp || r == j_tmp)
        {
            continue;
        }
        A[i_tmp][r] = c * A[r][i_tmp] - s * A[r][j_tmp];
        A[j_tmp][r] = c * A[r][j_tmp] + s * A[r][i_tmp];
        A[r][i_tmp] = A[i_tmp][r];
        A[r][j_tmp] = A[j_tmp][r];
    }

    A[i_tmp][i_tmp] = c * c * A_ii + s * s * A_jj - 2 * s * c * A[i_tmp][j_tmp];
    A[j_tmp][j_tmp] = s * s * A_ii + c * c * A_jj + 2 * s * c * A[i_tmp][j_tmp];
    A[i_tmp][j_tmp] = 0;
    A[j_tmp][i_tmp] = 0;
}

void update_V(double **V, int i_tmp, int j_tmp, double s, double c, int n)
{
    int i;
    double v1, v2;
    for (i = 0; i < n; i++)
    {
        v1 = V[i][i_tmp];
        v2 = V[i][j_tmp];
        V[i][i_tmp] = v1 * c - v2 * s;
        V[i][j_tmp] = v1 * s + v2 * c;
    }
}

void V_to_eigens(double **V, double **A, EIGEN *eigens, int n)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        eigens[i].vector = (double *)malloc(n * sizeof(double));
        assert(eigens[i].vector != NULL);
        eigens[i].value = A[i][i];
        for (j = 0; j < n; j++)
        {
            eigens[i].vector[j] = V[j][i];
        }
    }
}

int converge(double **A, double *old, int n)
{
    double new = calc_off(A, n);
    int converged = ((*old) - new <= 1.0e-15 ? 0 : 1);
    *old = new;
    return converged;
}

double calc_off(double **A, int n)
{
    int i, j;
    double new = 0;
    for (i = 0; i < n; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            new += 2 * A[i][j] * A[i][j];
        }
    }
    return new;
}

/*########## output ##########*/
void print_matrix(double **mat, int lines, int columns)
{
    int i, j;
    for (i = 0; i < lines - 1; i++)
    {
        for (j = 0; j < columns - 1; j++)
        {
            printf("%.4f,", -0.00005 < mat[i][j] && mat[i][j] < 0.0000 ? 0.0000 : mat[i][j]);
        }
        printf("%.4f\n", -0.00005 < mat[i][columns - 1] && mat[i][columns - 1] < 0.0000 ? 0.0000 : mat[i][columns - 1]);
    }
    for (j = 0; j < columns - 1; j++)
    {
        printf("%.4f,", -0.00005 < mat[lines - 1][j] && mat[lines - 1][j] < 0.0000 ? 0.0000 : mat[lines - 1][j]);
    }
    printf("%.4f", -0.00005 < mat[lines - 1][columns - 1] && mat[lines - 1][columns - 1] < 0.0000 ? 0.0000 : mat[lines - 1][columns - 1]);
}

int main(int argc, char *argv[])
{
    int n, k, d, max_iter = 300, goal, changed, i, j;
    FILE *fp;
    VEC *head;
    double **points, **wam, **ddg, **lnorm, **T, **centroids, **sums;
    EIGEN *eigens;
    int *sizes;

    if (argc != 4)
    {
        printf("Invalid Input!\n");
        exit(0);
    }
    else
    {
        if (is_non_negative_integer(argv[1]) == 0)
        {
            printf("Invalid Input!\n");
            exit(0);
        }
        k = atoi(argv[1]);

        goal = enum_goal(argv[2]);

        if (goal == 0)
        {
            printf("Invalid Input!\n");
            exit(0);
        }

        fp = fopen(argv[3], "r");
        assert(fp != NULL);
    }
    head = initialize(fp, &n, &d);
    fclose(fp);
    points = initialize_data(head, n, d);
    if (goal == 1)
    { /* spk */
        if (k >= n)
        {
            printf("Invalid Input!\n");
            exit(0);
        }
        wam = gen_wam(points, n, d);
        ddg = gen_ddg(wam, n);
        lnorm = gen_lnorm(wam, ddg, n);
        eigens = jacobi(lnorm, n);
        qsort(eigens, n, sizeof(EIGEN), my_comparator);
        if (k == 0)
        {
            k = determine_k(eigens, n);
        }
        d = k;
        T = gen_T(eigens, n, k);
        centroids = initialize_centroids(T, k, d);
        sizes = initialize_sizes(k);
        sums = initialize_sums(k, d);

        changed = 1;
        while (changed == 1 && max_iter > 0)
        {
            assign_vec_to_cluster(T, sizes, sums, centroids, n, k, d);
            changed = update_clusters(sizes, sums, centroids, k, d);
            max_iter -= 1;
        }
        free(wam[0]);
        free(wam);
        free(ddg[0]);
        free(ddg);
        free(lnorm[0]);
        free(lnorm);
        for (i = 0; i < n; i++)
        {
            free(eigens[i].vector);
        }
        free(eigens);
        free(points[0]);
        free(points);
        free(T[0]);
        free(T);
        free(sums[0]);
        free(sums);
        free(sizes);
        
        print_matrix(centroids, k, d);

        free(centroids[0]);
        free(centroids);
    }
    else if (goal == 2)
    { /* wam */
        wam = gen_wam(points, n, d);
        print_matrix(wam, n, n);

        free(wam[0]);
        free(wam);
    }
    else if (goal == 3)
    { /* ddg */
        wam = gen_wam(points, n, d);
        ddg = gen_ddg(wam, n);
        print_matrix(ddg, n, n);

        free(wam[0]);
        free(wam);
        free(ddg[0]);
        free(ddg);
    }
    else if (goal == 4)
    { /* lnorm */
        wam = gen_wam(points, n, d);
        ddg = gen_ddg(wam, n);
        lnorm = gen_lnorm(wam, ddg, n);
        print_matrix(lnorm, n, n);

        free(wam[0]);
        free(wam);
        free(ddg[0]);
        free(ddg);
        free(lnorm[0]);
        free(lnorm);
    }
    else
    { /* jacobi */
        wam = gen_wam(points, n, d);
        ddg = gen_ddg(wam, n);
        lnorm = gen_lnorm(wam, ddg, n);
        eigens = jacobi(lnorm, n);

        for (i = 0; i < n - 1; i++)
        {
            printf("%.4f,", -0.00005 < eigens[i].value && eigens[i].value < 0.0000 ? 0.0000 : eigens[i].value);
        }
        printf("%.4f\n", -0.00005 < eigens[n-1].value && eigens[n-1].value < 0.0000 ? 0.0000 : eigens[n-1].value);

        for (i = 0; i < n - 1; i++)
        {
            for (j = 0; j < n - 1; j++)
            {
                printf("%.4f,", -0.00005 < eigens[i].vector[j] && eigens[i].vector[j] < 0.0000 ? 0.0000 : eigens[i].vector[j]);
            }
            printf("%.4f\n", -0.00005 < eigens[i].vector[n-1] && eigens[i].vector[n-1] < 0.0000 ? 0.0000 : eigens[i].vector[n-1]);
        }
        for (j = 0; j < n - 1; j++)
        {
            printf("%.4f,", -0.00005 < eigens[n-1].vector[j] && eigens[n-1].vector[j] < 0.0000 ? 0.0000 : eigens[n-1].vector[j]);
        }
        printf("%.4f", -0.00005 < eigens[n-1].vector[n-1] && eigens[n-1].vector[n-1] < 0.0000 ? 0.0000 : eigens[n-1].vector[n-1]);
        free(wam[0]);
        free(wam);
        free(ddg[0]);
        free(ddg);
        free(lnorm[0]);
        free(lnorm);
        for (i = 0; i < n; i++)
        {
            free(eigens[i].vector);
        }
        free(eigens);
    }
    return 0;
}
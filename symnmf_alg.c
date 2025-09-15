#include "symnmf.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


static int is_delim(int c) {
    return c==' ' || c=='\t' || c==',' || c=='\n' || c=='\r';
}


Matrix mat_alloc(int rows, int cols) {
    Matrix M;
    M.rows = rows; M.cols = cols;
    M.data = (double*)calloc((size_t)rows * (size_t)cols, sizeof(double));
    return M;
}


void mat_free(Matrix *m) {
    if (m && m->data) { free(m->data); m->data = NULL; }
    if (m) { m->rows = 0; m->cols = 0; }
}


void mat_fill(Matrix *m, double v) {
    int i, n;
    n = m->rows * m->cols;
    for (i = 0; i < n; ++i) m->data[i] = v;
}


/* count rows & cols by first pass */
static int count_rows_cols(FILE *fp, int *rows, int *cols) {
    long pos;
    int r, c, in_tok, ch;
    *rows = 0; *cols = 0;
    pos = ftell(fp);
    if (pos < 0) return -1;
    r = 0; c = 0; in_tok = 0;
    while ((ch = fgetc(fp)) != EOF) {
        if (is_delim(ch)) {
            if (in_tok) { c++; in_tok = 0; }
            if (ch == '\n') {
                if (c > 0) {
                    if (*cols == 0) {
                        *cols = c;
                    }
                    r++;
                    c = 0;
                }
            }
        } else {
            in_tok = 1;
        }
    }
    if (in_tok) { 
        c++; 
        in_tok = 0; 
    }
    if (c > 0) { 
        if (*cols == 0) *cols = c;
        r++; 
    }
    *rows = r;
    if (fseek(fp, pos, SEEK_SET) != 0) return -1;
    return (r > 0 && *cols > 0) ? 0 : -1;
}


int read_points_file(const char *path, Matrix *X_out) {
    FILE *fp;
    int n, d, i, j;
    double val;
    fp = fopen(path, "r");
    if (!fp) return -1;
    if (count_rows_cols(fp, &n, &d) != 0) { fclose(fp); return -1; }
    *X_out = mat_alloc(n, d);
    for (i = 0; i < n; ++i) {
        for (j = 0; j < d; ++j) {
            if (fscanf(fp, "%lf", &val) != 1) { 
                fclose(fp); 
                mat_free(X_out); 
                return -1; 
            }
            MAT_AT(X_out, i, j) = val;
        }
    }
    return 0;
}
int print_matrix_csv(const Matrix *M) {
    int i, j;
    for (i = 0; i < M->rows; ++i) {
        for (j = 0; j < M->cols; ++j) {
            printf("%.4f", MAT_AT((Matrix*)M, i, j));
            if (j + 1 < M->cols) printf(",");
        }
        printf("\n");
    }
    return 0;
}
static double sq_euclid(const double *a, const double *b, int d) {
    int t; 
    double s = 0.0, diff;
    for (t = 0; t < d; ++t) { diff = a[t] - b[t]; s += diff * diff; }
    return s;
}
int compute_A(const Matrix *X, Matrix *A_out) {
    int n, d, i, j;
    double dist2;
    n = X->rows; d = X->cols;
    *A_out = mat_alloc(n, n);
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            if (i == j) {
                MAT_AT(A_out, i, j) = 0.0;
            } else {
                dist2 = sq_euclid(&X->data[i*d], &X->data[j*d], d);
                MAT_AT(A_out, i, j) = exp(-dist2 / 2.0);
            }
        }
    }
    return 0;
}

int compute_D(const Matrix *A, Matrix *D_out) {
    int n, i, j;
    double sum;
    n = A->rows;
    *D_out = mat_alloc(n, n);
    mat_fill(D_out, 0.0);
    for (i = 0; i < n; ++i) {
        sum = 0.0;
        for (j = 0; j < n; ++j) {
            sum += MAT_AT((Matrix*)A, i, j);
        }
        MAT_AT(D_out, i, i) = sum;
    }
    return 0;
}
int compute_W(const Matrix *A, const Matrix *D, Matrix *W_out) {
    int n, i, j;
    double dj, di, sdj, sdi;
    n = A->rows;
    *W_out = mat_alloc(n, n);
    for (i = 0; i < n; ++i) {
        di = MAT_AT((Matrix*)D, i, i);
        sdi = di > 0.0 ? sqrt(di) : 0.0;
        for (j = 0; j < n; ++j) {
            dj = MAT_AT((Matrix*)D, j, j);
            sdj = dj > 0.0 ? sqrt(dj) : 0.0;
            if (sdi == 0.0 || sdj == 0.0) MAT_AT(W_out, i, j) = 0.0;
            else MAT_AT(W_out, i, j) = MAT_AT((Matrix*)A, i, j) / (sdi * sdj);
        }
    }
    return 0;
}

int mat_mul(const Matrix *A, const Matrix *B, Matrix *C) {
    int i, j, k;
    double s;
    if (A->cols != B->rows) return -1;
    *C = mat_alloc(A->rows, B->cols);
    for (i = 0; i < A->rows; ++i) {
        for (j = 0; j < B->cols; ++j) {
            s = 0.0;
            for (k = 0; k < A->cols; ++k) {
                s += MAT_AT((Matrix*)A, i, k) * MAT_AT((Matrix*)B, k, j);
            }
            MAT_AT(C, i, j) = s;
        }
    }
    return 0;
}

int mat_transpose(const Matrix *A, Matrix *AT) {
    int i, j;
    *AT = mat_alloc(A->cols, A->rows);
    for (i = 0; i < A->rows; ++i){
        for (j = 0; j < A->cols; ++j){
            MAT_AT(AT, j, i) = MAT_AT((Matrix*)A, i, j);
        }
    }
    return 0;
}

/* beta fixed to 0.5 */
int symnmf_optimize(const Matrix *W, Matrix *H, int max_iter, double eps) {
    int n, k, it, i, j;
    double hij, num, den, factor, newv, d, diff2, beta = 0.5;
    Matrix WH, HT, HTH, HTHH;
    n = H->rows; k = H->cols;
    for (it = 0; it < max_iter; ++it) {
        /* WH = W * H */
        mat_mul(W, H, &WH);
        /* HTH = (H^T * H) */
        mat_transpose(H, &HT);
        mat_mul(&HT, H, &HTH);
        /* HTHH = H * HTH */
        mat_mul(H, &HTH, &HTHH);
        /* Update */
        diff2 = 0.0;
        for (i = 0; i < n; ++i) {
            for (j = 0; j < k; ++j) {
                hij = MAT_AT(H, i, j);
                num = MAT_AT(&WH, i, j);
                den = MAT_AT(&HTHH, i, j);
                factor = 0.0;
                if (den > 0.0) factor = (1.0 - beta) + beta * (num / den);
                else factor = (1.0 - beta);
                {
                    newv = hij * factor;
                    d = newv - hij;
                    MAT_AT(H, i, j) = newv < 0.0 ? 0.0 : newv; /* keep nonnegative */
                    diff2 += d * d;
                }
            }
        }
        if (diff2 < eps) {
            mat_free(&WH); mat_free(&HT); mat_free(&HTH); mat_free(&HTHH);
            break;
        }
        mat_free(&WH); mat_free(&HT); mat_free(&HTH); mat_free(&HTHH);
    }
    return 0;
}
#define _GNU_SOURCE
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kmeans.h"
#include "struct.h"

/* declaring functions:
double euclidean_distance(double *p1, double *p2, int dim);
int kmeans_algorithm(Point* all_points, Centroid* centroids_lst, int K, int d, int iter, int num_points, double eps); */

double euclidean_distance(double *p1, double *p2, int dim){
  int i;
  double distance = 0;
  for(i = 0; i < dim; i++){
    distance += pow((p1[i] - p2[i]), 2);
  }
  return sqrt(distance);
}

int kmeans_algorithm(Point* all_points, Centroid* centroids_lst, int K, int d, int iter, int num_points, double eps) {
    int cnt_iter = 0;
    int i, j, m;
    int min_idx = 0;
    double dis, min_dis;
    double flag = 1.0; /*indicates about epsilon*/
    double *new_coordinates = NULL;
    double *delta = (double*)calloc(K,sizeof(double));
    double max_delta;


    while (cnt_iter < iter && flag != 0.0) {
        for (i=0; i < K; i++) {
            centroids_lst[i].numPoints = 0;
            free(centroids_lst[i].Points);
            centroids_lst[i].Points = NULL;
        }

        for (i=0; i < num_points; i++) {
            min_dis = euclidean_distance(all_points[i].coordinates, centroids_lst[0].coordinates, d);
            min_idx = 0;
            for (j = 0; j < K; j++) {
                 /*computing minimal distance to one of the K clusters*/
                dis = euclidean_distance(all_points[i].coordinates, centroids_lst[j].coordinates, d);
                if (dis < min_dis) {
                    min_idx = j;
                    min_dis = dis;
                }
            }
        /*assign point to it's cluster*/
        centroids_lst[min_idx].Points = (Point*)realloc(centroids_lst[min_idx].Points, (centroids_lst[min_idx].numPoints+1)*sizeof(Point));
        if(centroids_lst[min_idx].Points == NULL) {
            for (i=0; i <K; i++) {
                free(centroids_lst[i].Points);
                }
            free(centroids_lst);
            free(new_coordinates);
            free(delta);
            return 0;
        }
        centroids_lst[min_idx].Points[centroids_lst[min_idx].numPoints]=all_points[i];
        centroids_lst[min_idx].numPoints++;
        }
        /*calculating delta & new centroids*/
        for ( m= 0; m < K; m++) {
            new_coordinates = (double *)calloc(d, sizeof(double));
            if (new_coordinates == NULL) {
                    /*free(new_coordinates);*/
                    free(delta);
                    return 0;
                }
            if (centroids_lst[m].numPoints > 0) {
                if (new_coordinates == NULL) {
                    /*free(new_coordinates);*/
                    free(delta);
                    return 0;
                }
                /* calculating average point of a cluster*/
                for (i = 0; i < centroids_lst[m].numPoints; i++) {
                    Point p = centroids_lst[m].Points[i];
                    for (j = 0; j < d; j++) {
                        new_coordinates[j] += p.coordinates[j];
                    }
                }

                for (i = 0; i < d; i++) {
                    new_coordinates[i] = new_coordinates[i]/(double)centroids_lst[m].numPoints;
                }
 
                /*for each centroid, calculating delta*/
                delta[m] = euclidean_distance(new_coordinates, centroids_lst[m].coordinates ,d);

                /* updating centroids*/
                for (i = 0; i < d; i++) {
                    centroids_lst[m].coordinates[i] = new_coordinates[i];
                }
                free(new_coordinates);
            }
            else {
                free(new_coordinates);
                for (i = 0; i < K; i++) {
                    free(centroids_lst[i].Points);
                }
                free(delta);
                return 0;
            }
        }

        /*checking coveragnce*/
        max_delta = delta[0];
        for (j = 0; j < K; j++) {
            if (delta[j] > max_delta) {
                max_delta = delta[j];
            }
        }
        
        if (max_delta < eps) {
            flag = 0.0; 
        }

    cnt_iter++;
    }

    free(delta);
    return 1;
}
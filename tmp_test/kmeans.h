#ifndef KMEANS_H
#define KMEANS_H

#include "struct.h"   /* for Point and Centroid */

double euclidean_distance(double *p1, double *p2, int dim);

int kmeans_algorithm(Point *all_points,
                     Centroid *centroids_lst,
                     int K,
                     int d,
                     int iter,
                     int num_points,
                     double eps);

#endif /* KMEANS_H */
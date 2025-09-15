import math
import sys

def validate_args():
    if len(sys.argv) not in [2, 3]:
        # invalid input
        print("An Error Has Occured")
        sys.exit(1)

    else:  #len(sys.argv) == 2 or 3
        try:
            K = float(sys.argv[1])
            if K > 1 and K.is_integer():
                centroids = int(K)
            else:
                raise ValueError
        except ValueError:
            print("Incorrect number of clusters!")
            sys.exit(1)
        if len(sys.argv) == 3:
            try:
                iter = float(sys.argv[2])
                if 1000 > iter > 1 and iter.is_integer():
                    iter = int(iter)
                else:
                    raise ValueError
            except ValueError:
                print("Incorrect maximum iteration!")
                sys.exit(1)
        else:
            iter = 400
        # len is 3 => default case: 400 iterations

        return centroids, iter


def Kmeans_algorithm(K, iter):
    def euclidean_distance(p1, p2):
        dis = 0
        for i in range(len(p1)):
            dis += (p1[i] - p2[i]) ** 2
        return math.sqrt(dis)

    epsilon = 0.001
    data_points = []
    for line in sys.stdin:
        if line.strip():
            coordinates = tuple(map(float, line.strip().split(',')))
            data_points.append(coordinates)
    #with open(input_data, 'r') as file:
    #    for line in file:
    #        coordinates = line.strip().split(',')
    #        coordinates = tuple(map(float, coordinates))
    #        data_points.append(coordinates)
    N = len(data_points)
    if K >= N:
        print("Incorrect number of clusters!")
        sys.exit(1)

    centroids_lst = []
    for i in range(K):
        centroids_lst.append(data_points[i])
    dim = len(data_points[0])

    for _ in range(iter):
        clusters = [[] for _ in range(K)]
        for point in data_points:
            min_centroid_dist = euclidean_distance(point, centroids_lst[0])
            min_idx = 0
            for i, centroid in enumerate(centroids_lst):
                dis = euclidean_distance(point, centroid)
                if dis < min_centroid_dist:
                    min_centroid_dist = dis
                    min_idx = i
            clusters[min_idx].append(point)

        delta = []
        for cluster_idx in range(len(clusters)):
            sum_points = [0] * dim
            for point in clusters[cluster_idx]:
                for j in range(dim):
                    sum_points[j] += point[j]
            # check epsilon condition & update centroids
            average_point = tuple(sum_point / len(clusters[cluster_idx]) for sum_point in sum_points)
            delta.append(euclidean_distance(average_point, centroids_lst[cluster_idx]))
            centroids_lst[cluster_idx] = average_point

        if max(delta) < epsilon:
            break

    #correct the output format
    centroids = []
    for point in centroids_lst:
        new_point = tuple("%.4f" % num for num in point)
        centroids.append(new_point)

    #output result
    for centroid in centroids:
        print(",".join(coord for coord in centroid))


def main():
    args = validate_args()
    Kmeans_algorithm(args[0], args[1])

if __name__ == "__main__":
    main()

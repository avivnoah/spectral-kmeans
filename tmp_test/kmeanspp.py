import random
import math
import numpy as np
import pandas as pd
import sys
import mykmeanspp


def validate_args():
    """
    1. We can assume filename is valid & provided.

    """
    if len(sys.argv) <= 4 or len(sys.argv) >= 7:
        # invalid input
        print("An Error Has Occurred")
        sys.exit(1)

    else:  #len(sys.argv) == 5 or 6
        try:
            K = float(sys.argv[1])
        except ValueError:
            print("Invalid number of clusters!")
            sys.exit(1)

        if (not K.is_integer()) or K <= 1:
            print("Invalid number of clusters!")
            sys.exit(1)
        if len(sys.argv) == 6:
            try:
                iter = float(sys.argv[2])
            except ValueError:
                print("Invalid maximum iteration!")
                sys.exit(1)

            if (not iter.is_integer()) or iter <= 1 or iter >= 1000:
                print("Invalid maximum iteration!")
                sys.exit(1)
            try:
                eps = float(sys.argv[3])
            except ValueError:
                print("Invalid epsilon!")
                sys.exit(1)

            if eps < 0:
                print("Invalid epsilon!")
                sys.exit(1)
        else:  # len(sys.argv) == 5, iter = 300
            iter = 300
            try:
                eps = float(sys.argv[2])
            except ValueError:
                print("Invalid epsilon!")
                sys.exit(1)

            if eps < 0:
                print("Invalid epsilon!")
                sys.exit(1)
    K = int(K)
    iter = int(iter)
    return K, iter, eps


def merge_and_read_files(file_name_1, file_name_2):  #works
    # Read the files
    df_file1 = pd.read_csv(file_name_1, header=None)
    df_file2 = pd.read_csv(file_name_2, header=None)

    # Merge the files
    merged_df = pd.merge(
        left=df_file1,
        right=df_file2,
        how='inner',
        left_on=df_file1.columns[0],
        right_on=df_file2.columns[0]
    )
    merged_df.sort_values(by=merged_df.columns[0], inplace=True)

    # Extract the sorted key column as a list
    key_lst = merged_df[merged_df.columns[0]].tolist()

    # Remove the key column and keep only the merged data
    merged_df = merged_df.iloc[:, 1:]

    # convert to numpy array
    points_array = merged_df.to_numpy()

    # Then the values to list
    points_list = merged_df.values.tolist()
    return points_array, key_lst, points_list


def init_centroids(points_array, key_lst, K):
    centroids = []  # K centroids
    selected_centroid_indices = []  # Keys of the centroids
    np.random.seed(1234)

    # randomly choose the first centroid
    first_centroid_idx = np.random.randint(len(key_lst))
    # Making sure the coordinates are float
    first_coordinates = [float(coord) for coord in points_array[first_centroid_idx]]
    centroids.append(first_coordinates)
    selected_centroid_indices.append(first_centroid_idx)

    # choose K-1 centroids using kmeans++
    for _ in range(K-1):

        # Compute the minimum distance from each point to any existing centroid
        distances = np.array([
            np.min([np.linalg.norm(point - centroid) for centroid in centroids])
            for point in points_array
        ])

        # Compute probability distribution proportional to distance
        probabilities = distances / distances.sum()

        # Choose next centroid index based on the computed probabilities
        next_idx = np.random.choice(len(points_array), p=probabilities)
        # Making sure the coordinates are float
        next_coordinates = [float(coord) for coord in points_array[next_idx]]
        centroids.append(next_coordinates)
        selected_centroid_indices.append(next_idx)

    return centroids, selected_centroid_indices


def K_means_algorithm(file_name1, file_name2, K, iter, eps):
    # Load and merge input data
    points_array, key_lst, points_list = merge_and_read_files(file_name1, file_name2)
    N = len(points_array)
    d = len(points_array[0])  # The dimensions of a point

    # Validate cluster count
    if K <= 1 or K >= N:
        print("Invalid number of clusters!")
        sys.exit(1)
    else:
        # Initialize centroids using k-means++
        centroids_lst, chosen_key_lst = init_centroids(points_array, key_lst, K)

        # Run the k-means clustering algorithm using C extension
        final_centroids_lst = mykmeanspp.fit(
            points_list,
            centroids_lst,
            K,
            iter,
            d,
            N,
            eps
        )

        # Format centroids to 4 decimal places
        formatted_centroids = [
            tuple(f"{coord:.4f}" for coord in centroid)
            for centroid in final_centroids_lst
        ]

        # Print final results
        chosen_key_lst = [int(num) for num in chosen_key_lst]
        joined_string = ','.join(map(str, chosen_key_lst))
        print(joined_string)
        for centroid in formatted_centroids:
            print(",".join(coord for coord in centroid))


def main():
    K, iter, eps = validate_args()
    if len(sys.argv) == 5:
        file_name_1 = sys.argv[3]
        file_name_2 = sys.argv[4]
    else:
        file_name_1 = sys.argv[4]
        file_name_2 = sys.argv[5]
    K_means_algorithm(file_name_1, file_name_2, K, iter, eps)


if __name__ == "__main__":
    main()

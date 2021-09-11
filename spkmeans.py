import sys
import numpy as np
import pandas as pd
import spkmeansmodule as spkm


def kmeans_pp_alg(data, n, k, d):
    np.random.seed(0)
    z = 1
    dist_prob = [0.0 for i in range(n)]
    centroids = []
    indexes = []
    index_first_centroid = np.random.choice(n)
    first_centroid = data[index_first_centroid]
    centroids.append(first_centroid)
    indexes.append(index_first_centroid)
    while (z != k):
        d_sum = 0
        for i in range(n):
            d_min = float('inf')
            for j in range(z):
                sub = np.subtract(data[i], centroids[j])
                dist = 0
                for l in range(d):
                    dist += sub[l] * sub[l]
                if (dist < d_min):
                    d_min = dist
                if (d_min == 0):
                    break
            dist_prob[i] = d_min
            d_sum += d_min
        for i in range(n):
            dist_prob[i] = dist_prob[i] / d_sum
        chosen_index = np.random.choice(n, replace=False, p=dist_prob)
        centroids.append(data[chosen_index])
        indexes.append(chosen_index)
        z += 1
    return (centroids, indexes)


def print_matrix(mat, lines, columns):
    for i in range(lines - 1):
        for j in range(columns - 1):
            print(f'{mat[i][j]:.4f}', end=',')
        print(f'{mat[i][columns - 1]:.4f}')
    for j in range(columns - 1):
        print(f'{mat[lines - 1][j]:.4f}', end=',')
    print(f'{mat[lines - 1][columns - 1]:.4f}', end='')


enum_goal = {"spk": 1, "wam": 2, "ddg": 3, "lnorm": 4, "jacobi": 5}

# main#
if (len(sys.argv) == 4):
    if not (sys.argv[1].isdigit()):
        sys.exit("Invalid Input!")
    k = int(sys.argv[1])
    goal = enum_goal.get(sys.argv[2], 0)
    if (goal == 0):
        sys.exit("Invalid Input!")
    file = pd.read_csv(sys.argv[3], header=None)
    #data = file.sort_values(by=0, ignore_index=True).iloc[:, 1:]
    data = pd.DataFrame.to_numpy(file)
    n = len(data)
    d = data.shape[1]
    if (n == 0 or d == 0):
        sys.exit("Invalid input!")
    max_iter = 300

    if (goal == 1):  # spk
        if (k >= n):
            sys.exit("Invalid input!")
        T = spkm.goal_switch(data, 1, n, k, d)
        k = len(T[0])
        d = k
        centroids_indexes_tuple = kmeans_pp_alg(T, n, k, d)
        clusters = spkm.fit(T, centroids_indexes_tuple[0], n, k, d, max_iter)
        for i in range(k - 1):
            print(centroids_indexes_tuple[1][i], end=",")
        print(centroids_indexes_tuple[1][k - 1])
        print_matrix(clusters, k, d)

    elif (goal == 2):
        wam = spkm.goal_switch(data.tolist(), 2, n, k, d)
        print_matrix(wam, n, n)
    elif (goal == 3):
        ddg = spkm.goal_switch(data.tolist(), 3, n, k, d)
        print_matrix(ddg, n, n)
    elif (goal == 4):
        lnorm = spkm.goal_switch(data.tolist(), 4, n, k, d)
        print_matrix(lnorm, n, n)
    else:  # goal == 5
        jacobi = spkm.goal_switch(data.tolist(), 5, n, k, d)
        print_matrix(jacobi, n + 1, n)
else:
    sys.exit("Invalid Input!")

import numpy


def distance_to_similarity(x):
    return 1.0 - (x / (1 + x))


# https://github.com/rubenqba/java-ml/blob/master/src/main/java/net/sf/javaml/distance/dtw/DTWSimilarity.java#L34
# https://blog.csdn.net/liyuefeilong/article/details/45748399
def measure(x: list, y: list):
    dP2P = numpy.zeros((len(x), len(y)))
    for i, x_value in enumerate(x):
        for j, y_value in enumerate(y):
            dP2P[i][j] = (x_value - y_value) ** 2

    if len(x) == 0 or len(y) == 0:
        return None

    if len(x) == 1 or len(y) == 0:
        return distance_to_similarity(dP2P[0][0] ** 0.5)

    D = numpy.zeros((len(x), len(y)))
    D[0][0] = dP2P[0][0]

    for i in range(1, len(x)):
        D[i][0] = dP2P[i][0] + D[i - 1][0]

    if len(y) == 1:
        sum = 0
        for i in range(len(x)):
            sum += D[i][0]
        return distance_to_similarity(sum ** 0.5 / len(x))

    for j in range(1, len(y)):
        D[0][j] = dP2P[0][j] + D[0][j - 1]

    if len(x) == 1:
        sum = 0
        for j in range(len(y)):
            sum += D[0][j]
        return distance_to_similarity(sum ** 0.5 / len(y))

    for i in range(1, len(x)):
        for j in range(1, len(y)):
            steps = [D[i - 1][j - 1], D[i - 1][j], D[i][j - 1]]
            step_min = min(steps)
            D[i][j] = dP2P[i][j] + step_min

    i = len(x) - 1
    j = len(y) - 1
    k = 1
    dist = D[i][j]

    while i + j > 2:
        if i == 0:
            j -= 1
        elif j == 0:
            i -= 1
        else:
            steps = [D[i - 1][j - 1], D[i - 1][j], D[i][j - 1]]
            step_min = min(steps)

            if step_min == steps[0]:
                i -= 1
                j -= 1
            elif step_min == steps[1]:
                i -= 1
            elif step_min == steps[2]:
                j -= 1
        k += 1
        dist += D[i][j]

    return distance_to_similarity(dist ** 0.5 / k)

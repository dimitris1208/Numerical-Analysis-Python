"""
THIS IS ANOTHER VERSION OF THE NEWTON-RAPHSON , BISECTION AND INTERSECTION METHODS AND THEIR ANALYSIS
"""



import math
import random
import numpy as np


def f(x):
    return 94 * (math.cos(x) ** 3) - 24 * math.cos(x) + 177 * (math.sin(x) ** 2) - 108 * (math.sin(x) ** 4) - 72 * (
                math.cos(x) ** 3) * (math.sin(x) ** 2) - 65


def df(x):
    return (216 * (math.cos(x) ** 2)) * (math.sin(x) ** 3) + (
                -144 * (math.cos(x) ** 4) - 282 * (math.cos(x) ** 2) + 354 * math.cos(x) + 24) * math.sin(x)


def ddf(x):
    return (432 - 432 * math.cos(x)) * (math.sin(x) ** 4) + (
                1224 * (math.cos(x) ** 3) - 1296 * (math.cos(x) ** 2) + 564 * math.cos(x) - 354) * (
                math.sin(x) ** 2) - 144 * (math.cos(x) ** 5) - 282 * (math.cos(x) ** 3) + 354 * (
                math.cos(x) ** 2) + 24 * (math.cos(x))


def g(x):
    return np.exp(np.sin(x) ** 3) + (x ** 6) - (2 * x ** 4) - (x ** 3) - 1


def dg(x):
    return np.exp(np.sin(x) ** 3) * 3 * (np.sin(x)) ** 2 * np.cos(x) + (6 * x ** 5) - (8 * x ** 3) - (3 * x ** 2)


def ddg(x):
    return 3 * ((math.e ** (math.sin(x) ** 3)) * math.sin(x) * (
                (math.sin(x) ** 2) * (3 * (math.cos(x) ** 2) * math.sin(x) - 1) + 2 * (math.cos(x) ** 2)) + 2 * (
                            x - 1) * x * (5 * x ** 2 + 5 * x + 1))


# --------------------------------------------------------------------------------------------------

def newton_method(x, acc_dig):
    print("Newton-Raphson Method")
    x_pr = x  # x previous
    x_ne = x_pr - 1 / ((df(x_pr) / f(x_pr)) - (ddf(x_pr) / 2 * df(x_pr)))  # x next
    while (abs(f(x_ne)) > (1 / 10 ** acc_dig)):
        x_pr = x_ne
        x_ne = x_pr - 1 / ((df(x_pr) / f(x_pr)) - (ddf(x_pr) / (2 * df(x_pr))))

    return x_ne


def bisection_method(pos_pos, pos_neg, acc_dig):
    print("Bisection Method")
    loops = 0

    if (pos_pos > pos_neg):
        # bis is a random number that divides the space
        bis = random.uniform(pos_neg, pos_pos)
    else:
        bis = random.uniform(pos_pos, pos_neg)

    while (abs(f(bis)) > (1 / 10 ** acc_dig)):
        loops += 1
        if (f(bis) > 0):
            pos_pos = bis
        else:
            pos_neg = bis
        if (pos_pos > pos_neg):
            bis = random.uniform(pos_neg, pos_pos)
        else:
            bis = random.uniform(pos_pos, pos_neg)
    return bis, loops


def intersection_method(x0, x1, x2, acc_dig):
    print("Intersection Method")
    xn = x0
    xn_1 = x1
    xn_2 = x2

    xn_3 = xn_2 - (((f(xn_2) * (xn_2 - xn_1) * f(xn)) / ((f(xn) - f(xn_1)) * (f(xn_2) - f(xn_1)))) - (
                (f(xn_2) * (xn_2 - xn) * f(xn_1)) / ((f(xn) - f(xn_1)) * (f(xn_2) - f(xn)))))
    while (abs(f(xn_3)) > (1 / 10 ** acc_dig)):
        xn = xn_1
        xn_1 = xn_2
        xn_2 = xn_3
        xn_3 = xn_2 - (((f(xn_2) * (xn_2 - xn_1) * f(xn)) / ((f(xn) - f(xn_1)) * (f(xn_2) - f(xn_1)))) - (
                    (f(xn_2) * (xn_2 - xn) * f(xn_1)) / ((f(xn) - f(xn_1)) * (f(xn_2) - f(xn)))))

    return xn_3


# --------------------------------------------------------------------------------------------------



print(newton_method(1, 5))
print("1st root\n")
print(newton_method(3, 5))
print("2nd root\n")
print(bisection_method(0, 1, 5))
print("1st root\n")
print(bisection_method(1.5, 1, 5))
print("2nd root\n")
print(bisection_method(1.5, 3, 5))
print("3rd root\n")
print(intersection_method(-0.5, 0, 0.6, 5))
print("1st root\n")
print(intersection_method(0, 1, 2, 5))
print("2nd root\n")
print(intersection_method(2, 2.5, 3, 5))
print("3rd root\n")

# --------------------------------------------------------------------------------------------------


def another_newton_method(x, acc_dig):
    print("Another Newton Raphson")
    loops = 0
    x_pr = x  # x previous
    x_ne = x_pr - 1 / ((dg(x_pr) / f(x_pr)) - (ddg(x_pr) / 2 * dg(x_pr)))  # x next
    while (abs(g(x_ne)) > (1 / 10 ** acc_dig)):
        x_pr = x_ne
        x_ne = x_pr - 1 / ((dg(x_pr) / g(x_pr)) - (ddg(x_pr) / (2 * dg(x_pr))))
        loops += 1
    return x_ne, loops


def another_intersection_method(x0, x1, x2, acc_dig):
    print("Another Intersection Method")
    xn = x0
    xn_1 = x1
    xn_2 = x2
    loops = 0
    xn_3 = xn_2 - (((g(xn_2) * (xn_2 - xn_1) * g(xn)) / (
            (g(xn) - g(xn_1)) * (g(xn_2) - g(xn_1)))) - (
                           (g(xn_2) * (xn_2 - xn) * g(xn_1)) / (
                           (g(xn) - g(xn_1)) * (g(xn_2) - g(xn)))))
    while (abs(g(xn_3)) > (1 / 10 ** acc_dig)):
        xn = xn_1
        xn_1 = xn_2
        xn_2 = xn_3
        xn_3 = xn_2 - (((g(xn_2) * (xn_2 - xn_1) * g(xn)) / (
                (g(xn) - g(xn_1)) * (g(xn_2) - g(xn_1)))) - (
                               (g(xn_2) * (xn_2 - xn) * g(xn_1)) / (
                               (g(xn) - g(xn_1)) * (g(xn_2) - g(xn)))))
        loops += 1
    return xn_3, loops


def another_bisection_method(pos_pos, pos_neg, acc_dig):
    print("Another Bisection Method")
    loops = 0
    if (pos_pos > pos_neg):
        # bis is a random number that divides the space
        bis = random.uniform(pos_neg, pos_pos)
    else:
        bis = random.uniform(pos_pos, pos_neg)

    while (abs(g(bis)) > (1 / 10 ** acc_dig)):
        loops += 1
        if (g(bis) > 0):
            pos_pos = bis
        else:
            pos_neg = bis
        if (pos_pos > pos_neg):
            bis = random.uniform(pos_neg, pos_pos)
        else:
            bis = random.uniform(pos_pos, pos_neg)
    return bis, loops


# The original bisection method needed 11 iterations to get to the root x = -1.1976
print(another_bisection_method(-2, 1, 5))


# The original Newton Raphson method needed 6 iterations to get to the root x = -1.1976
print(another_newton_method(-2, 5))


# The original intersection method needed 12 iterations to get to the root x = -1.1976
print(another_intersection_method(-2, -1.5, -1, 5))

#--------------------------------------------------------------------------------------------------

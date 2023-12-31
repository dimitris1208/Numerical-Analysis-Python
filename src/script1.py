import math
import numpy as np
import matplotlib.pyplot as plt


# Defining the function f(x): This function represents the equation for which we want to find the root.
def f(x):
    return np.exp(np.sin(x) ** 3) + (x ** 6) - (2 * x ** 4) - (x ** 3) - 1


""" 
 Defining the function df(x): This function represents the derivative of the equation.
 It is used in the Newton-Raphson method to calculate the next guess for the root. 
 """

def df(x):
    return np.exp(np.sin(x) ** 3) * 3 * (np.sin(x)) ** 2 * np.cos(x) + (6 * x ** 5) - (8 * x ** 3) - (3 * x ** 2)


"""
 Generating the x-values for plotting: The code generates an array of x-values using the numpy.linspace function.
 These values will be used to plot the function 
 """

x_list = np.linspace(-2, 2, 10000)

y_list = f(x_list)


# Plotting the function: The code uses matplotlib.pyplot to plot the function defined by f(x).
plt.figure(num=1, dpi=120)
plt.plot(x_list, y_list)
plt.title("Graph Exercise 1")

plt.savefig("Plot.png")

"""
 Implementing the root finding methods: The code defines three methods:
 bisection_method, newton_raphson_method, and intersection_method. Each method takes two initial
 guesses for the root and iteratively refines the guess until the root is found. 
 """
def intersection_method(pos_pos, pos_neg):
    print("Intersection Method")
    loops = 0
    if pos_neg < pos_pos:
        x_previous = pos_neg
        x_now = pos_pos
    else:
        x_previous = pos_pos
        x_now = pos_neg

    x_next = x_now - (f(x_now) * (x_now - x_previous) / (f(x_now) - f(x_previous)))

    while abs(f(x_next)) > 0.00001:
        x_previous = x_now
        x_now = x_next
        x_next = x_now - f(x_now) * (x_now - x_previous) / (f(x_now) - f(x_previous))
        loops += 1
    print("Root in position :", x_next)
    return loops


def newton_raphson_method(x):
    print("Newton Raphson Method")
    loops = 0
    x_0 = x
    x_1 = x_0 - f(x_0) / df(x_0)
    while abs(f(x_1)) > 0.00001:
        x_0 = x_1
        x_1 = x_0 - f(x_0) / df(x_0)
        loops += 1

    print("Root in position :", x_1)

    return loops


def bisection_method(pos_pos, pos_neg):
    print("Bisection Method ")
    loops = int(math.log(abs(pos_pos - pos_neg) / 0.00001) - math.log(2)) + 1
    print(loops)
    half = (pos_pos + pos_neg) / 2
    i = 0

    while i <= loops:

        if f(half) > 0:
            pos_pos = half
        else:
            pos_neg = half
        half = (pos_neg + pos_pos) / 2
        i += 1

    print("Root in position :", half)
    return loops


"""
 Calling the root finding methods: The code calls the root finding methods with
 different initial guesses and prints the number of iterations required to find the root. 
 """
print("It needed", bisection_method(-2, -1), "loops to find this root")

print("\n")

print("It needed", bisection_method(2, -1), "loops to find this root")

print("\n")

# f(-2) > 0 and f"(-2) > 0 , so f(-2)f"(-2)>0 , so x_0 = -2
print("It needed", newton_raphson_method(-2), "loops to find this root")


print("\n")

# f(2) > 0 and f"(2) > 0 , so f(2)f"(2)>0 , so x_0 = 2
print("It needed", newton_raphson_method(2), "loops to find this root")


print("\n")

print("It needed", intersection_method(-2, -1), "loops to find this root")

print("\n")

print("It needed", intersection_method(2, -1), "loops to find this root")

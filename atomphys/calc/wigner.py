from math import factorial as factorial_int
from math import floor, sqrt


def factorial(x: float) -> int:
    return factorial_int(int(x))


def ishalfint(x: float) -> bool:
    """check if value is half integer

    Arguments:
        x: is this a half integer?

    Returns:
        True if x is a half integer and False if it is not.
    """
    return 2 * x == floor(2 * x)


def isint(x: float) -> bool:
    """checks if value is an integer

    Arguments:
        x: is this an integer?

    Returns:
        True if x is an integer and False if it is not.
    """
    return x == floor(x)


def istriangle(a: float, b: float, c: float) -> bool:
    """checks if triad (a, b, c) obeys the triangle inequality

    Arguments:
        a:
        b:
        c:

    Returns:
        True if the triangle inequality is satisfied and False if it is not.
    """
    return abs(a - b) <= c and c <= a + b

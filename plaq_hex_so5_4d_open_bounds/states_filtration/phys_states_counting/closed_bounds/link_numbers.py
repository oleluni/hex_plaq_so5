import numpy as np


def nl1(ns: int, b: int) -> float or None:
    if ((ns//b) < 1) and ((ns % b) != b-1):
        return ns
    elif (ns % b) == b-1:
        return None
    else:
        return 2*ns - (b+(ns//b))


def nl2(ns: int, b: int) -> float or None:
    if (ns//b) < 1:
        return None
    elif ((ns//b) > 0) and ((ns % b) == b-1):
        return (2*ns+1) - (b+(ns//b)+1)
    else:
        return (2*ns+1) - (b+(ns//b))


def nl3(ns: int, b: int) -> float or None:
    if (ns % b) == 0:
        return None
    elif (ns//b)==0:
        return ns - 1
    else:
        return (2*ns - 2) - (b+(ns//b))


def nl4(ns: int, a: int, b: int) -> float or None:
    if (ns//b) + 1 == a:
        return None
    elif ((ns % b) == b-1) and ((ns//b) + 1 < a):
        return (2*ns + 1) - (b + (ns//b) + 1) + (2*b - 1)
    else:
        return (2 * ns + 1) - (b + (ns // b)) + (2*b - 1)


def nl(ns: int, a: int, b: int) -> tuple:
    return (nl1(ns=ns, b=b), nl2(ns=ns, b=b), nl3(ns=ns, b=b), nl4(ns=ns, a=a, b=b))

# [print(nl(ns=ns, a=5, b=4)) for ns in range(int(20))]

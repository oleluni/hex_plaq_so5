import numpy as np


def get_ravel_from_irrep_so5(_tuple: tuple) -> int:
    def get_ordinal_num(mL:float, mR:float) -> int:
        if (mL, mR) == (0, 0.5):
            return 0
        elif (mL, mR) == (0, -0.5):
            return 1
        elif (mL, mR) == (0.5, 0):
            return 2
        elif (mL, mR) == (-0.5, 0):
            return 3
        else:
            raise ValueError("Incorrect input, `mL` and `mR` may only take values 0, -0.5, 0.5 .")
    ordinal_nums: list = [get_ordinal_num(_tuple[ind], _tuple[ind+1]) for ind in range(0, len(_tuple), 2)]
    ordinal_nums_inv: list = ordinal_nums[::-1]

    mult_by4powers_nums: list = [ordinal_nums_inv[i] * (4**i) for i in range(len(ordinal_nums_inv))]
    return sum(mult_by4powers_nums)
# res1 = get_ravel_from_irrep_so5((0,0.5,0.5, 0, 0.5, 0, 0, -0.5))
# res2 = get_ravel_from_irrep_so5((0,0.5,-0.5, 0, 0.5, 0, 0, 0.5))
#
# print(res1)
# print(res2)

def get_irrep_from_ravel_so5():
    pass
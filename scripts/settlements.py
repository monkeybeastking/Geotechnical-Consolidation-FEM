import numpy as np


def get_settlement(Mv:float, spacing:float, initial_stress:float, Load:float):
    settlement = 0
    for i in initial_stress:
        if i > Load:
            pass # we are not looking at these types of laods
        elif i<= Load:
            settlement_i = Mv * spacing * i
            settlement += settlement_i 
        else:
            break
    return settlement


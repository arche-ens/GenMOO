import math

# QED parameters from Bickerton G., et al. Nature Chem. 2012
QED_PARAMS = {
    "MW":     (2.817, 392.575, 290.749, 2.420, 49.223, 65.371, 104.981),
    "ALOGP":  (3.173, 137.862, 2.535, 4.581, 0.823, 0.576, 131.319),
    "HBD":    (1.619, 1010.051, 0.985, 1e-12, 0.714, 0.921, 258.163),
    "HBA":    (2.949, 160.461, 3.615, 4.436, 0.290, 1.301, 148.776),
    "PSA":    (1.877, 125.223, 62.908, 87.834, 12.020, 28.513, 104.569),
    "ROTB":   (0.010, 272.412, 2.558, 1.566, 1.272, 2.758, 105.442),
    "AROM":   (3.218, 957.737, 2.275, 1e-12, 1.318, 0.376, 312.337),
    "ALERTS": (0.010, 1199.094, -0.090, 1e-12, 0.186, 0.875, 417.725),
}

DESCRIPTOR_ORDER = ["MW", "ALOGP", "HBD", "HBA", "PSA", "ROTB", "AROM", "ALERTS"]


def desirability(x, params):
    """Compute desirability d(x) using ADS function.
    
    Parameters
    ---
        x : float
            Descriptor value
        params : list
            a, b, c, d, e, f, dmax
    """
    a, b, c, d, e, f, dmax = params
    sig1 = 1.0 / (1.0 + math.exp(-(x - c + d/2.0) / e))
    sig2 = 1.0 - 1.0 / (1.0 + math.exp(-(x - c - d/2.0) / f))
    desire = (a + b * sig1 * sig2) / dmax

    return desire


def calculate_qed(values: list[float], weights: list[float]):
    """Compute weighted QED from 8 descriptor values and 8 weights.
    
    Parameters
    ---
        values : list
            List of descriptor values (MW, ALOGP, HBD, HBA, PSA, ROTB, AROM, ALERTS).
        weights : list
    """
    desires = []
    for x, name in zip(values, DESCRIPTOR_ORDER):
        params = QED_PARAMS[name]
        desires.append(desirability(x, params))

    numerator = sum(w * math.log(d) for w, d in zip(weights, desires))
    denominator = sum(weights)
    qed = math.exp(numerator / denominator)

    return qed


# Example usage:
if __name__ == "__main__":
    x_values = [263.315, 0.96, 4, 5, 106.86, 3, 1, 0]  # 8 descriptor values
    # w_values = [0.1417, 0.1108, 0.0709, 0.0639, 0.1604, 0.0757, 0.3288, 0.0479]  # 8 weights
    w_values = [0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95]  # MEAN 8 weights

    qed = round(calculate_qed(x_values, w_values), 3)
    print("QED_w =", qed)

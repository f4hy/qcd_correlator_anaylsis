#!/usr/bin/env python
import numpy as np
import scipy.optimize
import logging
import warnings


def newton_cosh_for_m(i, j, ave_cor, guess, T):

    def f(m):
        numerator = np.exp(-m*i) + np.exp(-1.0*m*(T-i))
        denominator = np.exp(-m*j) + np.exp(-1.0*m*(T-j))
        return numerator/denominator - ave_cor[i]/ave_cor[j]

    logging.debug("newtons method to find cosh (T={}) emass starting guess"
                  " {}, {},{}, fguess={}".format(T, guess, i, j, f(guess)))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            result = scipy.optimize.newton(f, guess, tol=1.48e-6)
        except RuntimeError:
            logging.error("Newtons failed to converge using standard for {},{}".format(i, j))
            logging.error("Newtons failed, fallback is {}".format(guess))
            return guess
        logging.debug("newtons method converged to {} from {}".format(result, guess))
    return result


def newton_sinh_for_m(i, j, ave_cor, guess, T):

    def f(m):
        numerator = np.exp(-m*i) - np.exp(-1.0*m*(T-i))
        denominator = np.exp(-m*j) - np.exp(-1.0*m*(T-j))
        return numerator/denominator - ave_cor[i]/ave_cor[j]

    logging.debug("newtons method to find sinh (T={}) emass starting guess"
                  " {}, {},{}, fguess={}".format(T, guess, i, j, f(guess)))
    try:
        result = scipy.optimize.newton(f, guess)
    except RuntimeError:
        logging.error("Newtons failed to converge using standard for {},{}".format(i, j))
        logging.error("Newtons failed, fallback is {}".format(guess))
        return guess
    logging.debug("newtons method converged to {}".format(result))
    return result

#!/usr/bin/python

from scipy.special import gammaincc;
from math import floor;

def poisson_cdf(k,lam):
	return gammaincc(floor(k+1),lam);

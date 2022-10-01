# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 11:16:31 2022

@author: Jungle
"""

import netCDF4 as nc
import numpy as np
import bisect
import math
import matplotlib.pyplot as plt
import re
import os 


df = nc.Dataset("D:/学科竞赛/管悦程序/zzytest.nc")
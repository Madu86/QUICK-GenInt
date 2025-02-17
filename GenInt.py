#!---------------------------------------------------------------------!
#! Written by Madu Manathunga on 07/01/2021                            !
#!                                                                     !
#! Copyright (C) 2020-2021 Merz lab                                    !
#! Copyright (C) 2020-2021 Götz lab                                    !
#!                                                                     !
#! This source file is a part of QUICK-GenInt code generator and       !
#! is subjected to the terms of the Mozilla Public License, v. 2.0.    !
#! If a copy of the MPL was not distributed with this file, you can    !
#! obtain one at http://mozilla.org/MPL/2.0/.                          !
#!_____________________________________________________________________!

#!---------------------------------------------------------------------!
#! This is the main source file of QUICK integral code generator.      !
#!---------------------------------------------------------------------!

import sys
import os

# Get the absolute path of oei directory and set module search path
GenIntHome=os.path.abspath(os.getcwd())
outdir_name="output"
outdir=os.path.join(GenIntHome,outdir_name)

try:
    os.mkdir(outdir)
except OSError as error:
    print(error)

#common_abs_path=os.path.abspath(os.getcwd())+"/src/common"
#oei_abs_path=os.path.abspath(os.getcwd())+"/src/oei"
#sys.path.append(common_abs_path)
#sys.path.append(oei_abs_path)

import src.oei.one_electron_integral as one_electron_integral

# set function qualifiers. If you want to generate host code, leave an empty string.
# The default is for generating cuda code. 
func_qualifier='__device__ __inline__'

# generate one electron integral source code
one_electron_integral.write_oei(outdir, func_qualifier)

def print_store():

    labels = ("S", "Px", "Py", "Pz", "Dxy", "Dyz", "Dxz", "Dxx", "Dyy", "Dzz", "Fxyz", "Fxxy", "Fxyy", "Fxxz", "Fxzz", "Fyyz", "Fyzz", "Fxxx", "Fyyy", "Fzzz",\
        "Gxxyy", "Gxxzz", "Gyyzz", "Gxxyz", "Gxyyz", "Gxyzz", "Gxxxz", "Gxzzz", "Gxxxy", "Gxyyy", "Gyyyz", "Gyzzz", "Gxxxx", "Gyyyy", "Gzzzz")

    for i in range(0, 35):
        str=''
        for j in range(0,35):
            str += "[" + labels[i] + "|" + labels[j] + "]"
        print(str)

print_store()
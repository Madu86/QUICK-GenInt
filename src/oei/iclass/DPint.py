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
#! This source file contains classes necessary for generating one      !
#! electron integrals. Note that we use vertical recurrence relations  !
#! algorithm developed by Obara and Saika. See J. Chem. Phys. 1986, 84,!
#! 3963−3974 paper for theoretical details.                            !
#!                                                                     !
#!---------------------------------------------------------------------!

import src.common.params as params
import src.common.file_handler as file_handler
from src.oei.iclass.OEint import OEint

# [d|p] class, subclass of OEint
class DPint(OEint):
    def gen_int(self):
        # write code paths for integrals. Note that we use C++ classes here.
        for m in range(0,self.max_m+1):
            if m == 0:
                self.fhc.write("\n/* DP true integral, m=%d */ \n" % (m)) 
                self.fhd.write("\n/* DP true integral, m=%d */ \n" % (m)) 
            else:
                self.fhc.write("\n/* DP auxilary integral, m=%d */ \n" % (m))
                self.fhd.write("\n/* DP auxilary integral, m=%d */ \n" % (m))              

            self.fhc.write("class DPint_%d{ \n" % (m))
            self.fhc.write("public: \n")

            # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
            for i in range(0,6):
                for j in range(0,3):
                    self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (i+4, j+1, self.d_lbl[i], self.p_lbl[j]))

            # write class functions
            self.fhc.write("  %s DPint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz,\n\
                QUICKDouble Zeta, QUICKDouble* YVerticalTemp); \n" % (self.func_qualifier, m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("%s DPint_%d::DPint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz,\n\
                QUICKDouble Zeta, QUICKDouble* YVerticalTemp){ \n\n" % (self.func_qualifier, m, m))
            self.fhd.write("  PSint_%d ps_%d(PAx, PAy, PAz, PCx, PCy, PCz, YVerticalTemp); // construct [p|s] for m=%d \n" % (m, m, m))
            self.fhd.write("  PSint_%d ps_%d(PAx, PAy, PAz, PCx, PCy, PCz, YVerticalTemp); // construct [p|s] for m=%d \n" % (m+1, m+1, m+1))
            self.fhd.write("  DSint_%d ds_%d(PAx, PAy, PAz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [d|s] for m=%d \n" % (m, m, m))            
            self.fhd.write("  DSint_%d ds_%d(PAx, PAy, PAz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [d|s] for m=%d \n\n" % (m+1, m+1, m+1))

            for i in range(0,6):
                for j in range(0,3):
                    tmp_mcal=[params.Mcal[i+4][0], params.Mcal[i+4][1], params.Mcal[i+4][2]]
                    for k in range(0,3):
                        if params.Mcal[j+1][k] != 0:
                            self.fhd.write("  x_%d_%d = %s * ds_%d.x_%d_%d - %s * ds_%d.x_%d_%d; \n" % (i+4, j+1, self.PB[k], m, i+4, 0,\
                            self.PC[k], m+1, i+4, 0))


                            if tmp_mcal[k] != 0:
                                tmp_mcal[k] -= 1
                                tmp_i=params.trans[tmp_mcal[0]][tmp_mcal[1]][tmp_mcal[2]]
                                self.fhd.write("  x_%d_%d += 0.5/Zeta * %f * (ps_%d.x_%d_%d - ps_%d.x_%d_%d); \n" % (i+4, j+1, params.Mcal[i+4][k], m, tmp_i-1, 0, m+1, tmp_i-1, 0))

                            break
            self.fhd.write("\n } \n")


    # generate code to save computed [d|p] integral
    def save_int(self):
        self.fha.write("\n  /* DP integral, m=%d */ \n" % (0))
        self.fha.write("  if(I == 2 && J == 1){ \n")
        self.fha.write("    DPint_0 dp(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n")
        for i in range(0,6):
            for j in range(0,3):
                self.fha.write("    LOCSTORE(store, %d, %d, STOREDIM, STOREDIM) += dp.x_%d_%d;\n" % (i+4, j+1, i+4, j+1))

        # include print statements if debug option is on    
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            for i in range(0,6):
                for j in range(0,3):
                    self.fha.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOCSTORE(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "DP", i+4, j+1, i+4, j+1))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")

    # generate code to save [d|p] integral gradients
    def save_int_grad(self):
        self.fhga.write("\n  /* DP integral gradient, m=%d */ \n" % (0))
        self.fhga.write("  if(I == 2 && J == 1){ \n")
        self.fhga.write("    DSint_0 ds(PAx, PAy, PAz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n")
        self.fhga.write("    PPint_0 pp(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n")
        self.fhga.write("    DDint_0 dd(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, store, YVerticalTemp); \n")
        self.fhga.write("    FPint_0 fp(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n\n")

        for i in range(0,6):
            self.fhga.write("    LOCSTORE(store, %d, %d, STOREDIM, STOREDIM) = ds.x_%d_%d;\n" % (i+4, 0, i+4, 0))

        for i in range(0,3):
            for j in range(0,3):
                self.fhga.write("    LOCSTORE(store, %d, %d, STOREDIM, STOREDIM) = pp.x_%d_%d;\n" % (i+1, j+1, i+1, j+1))

        for i in range(0,6):
            for j in range(0,6):
                self.fhga.write("    LOCSTORE(store, %d, %d, STOREDIM, STOREDIM) = dd.x_%d_%d;\n" % (i+4, j+4, i+4, j+4))

        for i in range(0,10):
            for j in range(0,3):
                self.fhga.write("    LOCSTORE(store, %d, %d, STOREDIM, STOREDIM) = fp.x_%d_%d;\n" % (i+10, j+1, i+10, j+1))

        if OEint.debug == 1:
            self.fhga.write("\n#ifdef DEBUG_OEI \n")

            for i in range(0,6):
                self.fhga.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOCSTORE(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "DS", i+4, 0, i+4, 0))

            for i in range(0,3):
                for j in range(0,3):
                    self.fhga.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOCSTORE(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "PP", i+1, j+1, i+1, j+1))

            for i in range(0,6):
                for j in range(0,6):
                    self.fhga.write("    printf(\"II %%d JJ %%d %s store%d,%d] = %%f \\n\", II, JJ, LOCSTORE(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "DD", i+4, j+4, i+4, j+4))

            for i in range(0,10):
                for j in range(0,3):
                    self.fhga.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOCSTORE(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "FP", i+10, j+1, i+10, j+1))
            self.fhga.write("#endif \n\n")

        self.fhga.write("  } \n")
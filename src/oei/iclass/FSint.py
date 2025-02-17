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

# [f|s] class, subclass of OEint
class FSint(OEint):
    def gen_int(self):
        # write code paths for integrals. Note that we use C++ classes here.
        for m in range(0,self.max_m+1):
            if m == 0:
                self.fhc.write("\n/* FS true integral, m=%d */ \n" % (m)) 
                self.fhd.write("\n/* FS true integral, m=%d */ \n" % (m)) 
            else:
                self.fhc.write("\n/* FS auxilary integral, m=%d */ \n" % (m))
                self.fhd.write("\n/* FS auxilary integral, m=%d */ \n" % (m))              

            self.fhc.write("class FSint_%d{ \n" % (m))
            self.fhc.write("public: \n")

            # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
            self.fhc.write("#ifdef REG_FS \n")
            for i in range(0,10):
                self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (i+10, 0, self.f_lbl[i], "S"))
            self.fhc.write("#endif \n") 

            # write class functions
            self.fhc.write("  %s FSint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz, QUICKDouble TwoZetaInv, QUICKDouble* store, QUICKDouble* YVerticalTemp); \n" % (self.func_qualifier, m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("%s FSint_%d::FSint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz, QUICKDouble TwoZetaInv, QUICKDouble* store, QUICKDouble* YVerticalTemp){ \n\n" % (self.func_qualifier, m, m))
            self.fhd.write("  PSint_%d ps_%d(PAx, PAy, PAz, PCx, PCy, PCz, store, YVerticalTemp); // construct [p|s] for m=%d \n" % (m, m, m))
            self.fhd.write("  PSint_%d ps_%d(PAx, PAy, PAz, PCx, PCy, PCz, store, YVerticalTemp); // construct [p|s] for m=%d \n" % (m+1, m+1, m+1))
            self.fhd.write("  DSint_%d ds_%d(PAx, PAy, PAz, PCx, PCy, PCz, TwoZetaInv, store, YVerticalTemp); // construct [d|s] for m=%d \n" % (m, m, m))
            self.fhd.write("  DSint_%d ds_%d(PAx, PAy, PAz, PCx, PCy, PCz, TwoZetaInv, store, YVerticalTemp); // construct [d|s] for m=%d \n\n" % (m+1, m+1, m+1))

            # save all computed values into class variables that will reside in register/lmem space
            self.fhd.write("#ifdef REG_FS \n")
            for i in range(0,10):
                tmp_mcal=[params.Mcal[i+10][0], params.Mcal[i+10][1], params.Mcal[i+10][2]]
                for j in range(0,3):
                    if params.Mcal[i+10][j] != 0:
                        tmp_mcal[j] -= 1
                        tmp_i=params.trans[tmp_mcal[0]][tmp_mcal[1]][tmp_mcal[2]]
                        self.fhd.write("  x_%d_%d = %s * ds_%d.x_%d_%d - %s * ds_%d.x_%d_%d; \n" % (i+10, 0, self.PA[j], m, tmp_i-1, 0,\
                        self.PC[j], m+1, tmp_i-1, 0))

                        if params.Mcal[i+10][j] > 1:
                            tmp_mcal[j] = params.Mcal[i+10][j] - 2
                            tmp_i=params.trans[tmp_mcal[0]][tmp_mcal[1]][tmp_mcal[2]]
                            self.fhd.write("  x_%d_%d += TwoZetaInv * %f * (ps_%d.x_%d_%d - ps_%d.x_%d_%d); \n" % (i+10, 0, params.Mcal[i+10][j] - 1,\
                                m, tmp_i-1, 0, m+1, tmp_i-1, 0))

                        break
            self.fhd.write("#else \n")

            # save all computed values into store array in global memory
            self.fhd.write("  QUICKDouble val; \n")
            for i in range(0,10):
                tmp_mcal=[params.Mcal[i+10][0], params.Mcal[i+10][1], params.Mcal[i+10][2]]
                for j in range(0,3):
                    if params.Mcal[i+10][j] != 0:
                        tmp_mcal[j] -= 1
                        tmp_i=params.trans[tmp_mcal[0]][tmp_mcal[1]][tmp_mcal[2]]
                        self.fhd.write("  val = %s * ds_%d.x_%d_%d - %s * ds_%d.x_%d_%d; \n" % (self.PA[j], m, tmp_i-1, 0,\
                        self.PC[j], m+1, tmp_i-1, 0))

                        if params.Mcal[i+10][j] > 1:
                            tmp_mcal[j] = params.Mcal[i+10][j] - 2
                            tmp_i=params.trans[tmp_mcal[0]][tmp_mcal[1]][tmp_mcal[2]]
                            self.fhd.write("  val += TwoZetaInv * %f * (ps_%d.x_%d_%d - ps_%d.x_%d_%d); \n" % ( params.Mcal[i+10][j] - 1,\
                                m, tmp_i-1, 0, m+1, tmp_i-1, 0))
                        
                        self.fhd.write("  LOCSTOREFULL(store, %d, %d, STOREDIM, STOREDIM, %d) = val; \n" % (i+10, 0, m))
                        break
            self.fhd.write("#endif \n") 
            self.fhd.write("\n } \n")

    # generate code to save computed [d|s] integral
    def save_int(self):
        self.fha.write("\n  // FS integral, m=%d  \n" % (0))
        self.fha.write("  if(I == 3 && J == 0){ \n")
        self.fha.write("    FSint_0 fs(PAx, PAy, PAz, PCx, PCy, PCz, TwoZetaInv, store, YVerticalTemp); \n")

        self.fha.write("#ifdef REG_FS \n")
        for i in range(0,10):
            self.fha.write("    LOCSTORE(store, %d, %d, STOREDIM, STOREDIM) = fs.x_%d_%d;\n" % (i+10, 0, i+10, 0))
        self.fha.write("#endif \n") 

        # include print statements if debug option is on    
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            for i in range(0,10):
                self.fha.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOCSTORE(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "FS", i+10, 0, i+10, 0))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")

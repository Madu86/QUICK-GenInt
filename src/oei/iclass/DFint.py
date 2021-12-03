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

# [d|f] class, subclass of OEint
class DFint(OEint):
    def gen_int(self):
        # write code paths for integrals. Note that we use C++ classes here.
        for m in range(0,self.max_m+1):
            if m == 0:
                self.fhc.write("\n/* DF true integral, m=%d */ \n" % (m)) 
                self.fhd.write("\n/* DF true integral, m=%d */ \n" % (m)) 
            else:
                self.fhc.write("\n/* DF auxilary integral, m=%d */ \n" % (m))
                self.fhd.write("\n/* DF auxilary integral, m=%d */ \n" % (m))              

            self.fhc.write("class DFint_%d{ \n" % (m))
            self.fhc.write("public: \n")

            # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
            self.fhc.write("#ifdef REG_DF \n")
            for i in range(0,10):
                for j in range(0,6):
                    self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (j+4, i+10, self.d_lbl[j], self.f_lbl[i]))
            self.fhc.write("#endif \n") 

            # write class functions
            self.fhc.write("  %s DFint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz,\n\
                QUICKDouble TwoZetaInv, QUICKDouble* store, QUICKDouble* YVerticalTemp); \n" % (self.func_qualifier, m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("%s DFint_%d::DFint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz,\n\
                QUICKDouble TwoZetaInv, QUICKDouble* store, QUICKDouble* YVerticalTemp){ \n\n" % (self.func_qualifier, m, m))
            self.fhd.write("  PDint_%d pd_%d(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, TwoZetaInv, store, YVerticalTemp); // construct [p|d] for m=%d \n" % (m, m, m))
            self.fhd.write("  SFint_%d sf_%d(PBx, PBy, PBz, PCx, PCy, PCz, TwoZetaInv, store, YVerticalTemp); // construct [s|f] for m=%d \n" % (m, m, m))
            self.fhd.write("  PFint_%d pf_%d(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, TwoZetaInv, store, YVerticalTemp); // construct [p|f] for m=%d \n" % (m, m, m))            
            self.fhd.write("  PDint_%d pd_%d(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, TwoZetaInv, store, YVerticalTemp); // construct [p|d] for m=%d \n" % (m+1, m+1, m+1))
            self.fhd.write("  SFint_%d sf_%d(PBx, PBy, PBz, PCx, PCy, PCz, TwoZetaInv, store, YVerticalTemp); // construct [s|f] for m=%d \n" % (m+1, m+1, m+1))
            self.fhd.write("  PFint_%d pf_%d(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, TwoZetaInv, store, YVerticalTemp); // construct [p|f] for m=%d \n\n" % (m+1, m+1, m+1))             

            # save all computed values into class variables that will reside in register/lmem space
            self.fhd.write("#ifdef REG_DF \n")
            for i in range(0,10):
                for j in range(0,6):
                    tmp_mcal1=[params.Mcal[i+10][0], params.Mcal[i+10][1], params.Mcal[i+10][2]]
                    tmp_mcal2=[params.Mcal[j+4][0], params.Mcal[j+4][1], params.Mcal[j+4][2]]

                    for k in range(0,3):
                        if params.Mcal[j+4][k] != 0:
                            tmp_mcal2[k] -= 1
                            tmp_j=params.trans[tmp_mcal2[0]][tmp_mcal2[1]][tmp_mcal2[2]]
                            iclass_obj="pf"
                            self.fhd.write("  x_%d_%d = %s * %s_%d.x_%d_%d - %s * %s_%d.x_%d_%d; \n" % (j+4, i+10, self.PA[k], iclass_obj, m, tmp_j-1, i+10,\
                            self.PC[k], iclass_obj, m+1, tmp_j-1, i+10))

                            if params.Mcal[j+4][k] > 1:
                                iclass_obj="sf"
                                self.fhd.write("  x_%d_%d += TwoZetaInv * (%s_%d.x_%d_%d - %s_%d.x_%d_%d); \n" % (j+4, i+10, iclass_obj, m, 0, i+10,\
                                iclass_obj, m+1, 0, i+10))

                            if tmp_mcal1[k] > 0:
                                tmp_mcal1[k] -= 1
                                tmp_i=params.trans[tmp_mcal1[0]][tmp_mcal1[1]][tmp_mcal1[2]]
                                iclass_obj="pd"
                                self.fhd.write("  x_%d_%d += TwoZetaInv * %f * (%s_%d.x_%d_%d - %s_%d.x_%d_%d); \n" % (j+4, i+10, params.Mcal[i+10][k], iclass_obj, m, tmp_j-1, tmp_i-1,\
                                iclass_obj, m+1, tmp_j-1, tmp_i-1))
                            break
            self.fhd.write("#else \n")

            # save all computed values into store array in global memory
            self.fhd.write("  QUICKDouble val; \n")
            for i in range(0,10):
                for j in range(0,6):
                    tmp_mcal1=[params.Mcal[i+10][0], params.Mcal[i+10][1], params.Mcal[i+10][2]]
                    tmp_mcal2=[params.Mcal[j+4][0], params.Mcal[j+4][1], params.Mcal[j+4][2]]

                    for k in range(0,3):
                        if params.Mcal[j+4][k] != 0:
                            tmp_mcal2[k] -= 1
                            tmp_j=params.trans[tmp_mcal2[0]][tmp_mcal2[1]][tmp_mcal2[2]]
                            iclass_obj="pf"
                            self.fhd.write("#ifdef REG_PF \n")
                            self.fhd.write("  val = %s * %s_%d.x_%d_%d - %s * %s_%d.x_%d_%d; \n" % (self.PA[k], iclass_obj, m, tmp_j-1, i+10,\
                            self.PC[k], iclass_obj, m+1, tmp_j-1, i+10))
                            self.fhd.write("#else \n")
                            self.fhd.write("  val = %s * LOCSTOREFULL(store, %d, %d, STOREDIM, STOREDIM, %d) - %s * LOCSTOREFULL(store, %d, %d, STOREDIM, STOREDIM, %d); \n" % (self.PA[k], tmp_j-1, i+10, m,\
                            self.PC[k], tmp_j-1, i+10, m+1))
                            self.fhd.write("#endif \n")

                            if params.Mcal[j+4][k] > 1:
                                iclass_obj="sf"

                                self.fhd.write("#ifdef REG_SF \n")
                                self.fhd.write("  val += TwoZetaInv * (%s_%d.x_%d_%d - %s_%d.x_%d_%d); \n" % (iclass_obj, m, 0, i+10,\
                                iclass_obj, m+1, 0, i+10))
                                self.fhd.write("#else \n")
                                self.fhd.write("  val += TwoZetaInv * (LOCSTOREFULL(store, %d, %d, STOREDIM, STOREDIM, %d) - LOCSTOREFULL(store, %d, %d, STOREDIM, STOREDIM, %d)); \n" % (0, i+10, m,\
                                0, i+10, m+1))
                                self.fhd.write("#endif \n")

                            if tmp_mcal1[k] > 0:
                                tmp_mcal1[k] -= 1
                                tmp_i=params.trans[tmp_mcal1[0]][tmp_mcal1[1]][tmp_mcal1[2]]
                                iclass_obj="pd"
                                self.fhd.write("  val += TwoZetaInv * %f * (%s_%d.x_%d_%d - %s_%d.x_%d_%d); \n" % (params.Mcal[i+10][k], iclass_obj, m, tmp_j-1, tmp_i-1,\
                                iclass_obj, m+1, tmp_j-1, tmp_i-1))

                            self.fhd.write("  LOCSTOREFULL(store, %d, %d, STOREDIM, STOREDIM, %d) = val; \n" % (j+4, i+10, m))
                            break
            self.fhd.write("#endif \n") 
            self.fhd.write("\n } \n")

    # generate code to save computed [d|f] integral
    def save_int(self):
        self.fha.write("\n  /* DF integral, m=%d */ \n" % (0))
        self.fha.write("  if(I == 2 && J == 3){ \n")
        self.fha.write("    DFint_0 df(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, TwoZetaInv, store, YVerticalTemp); \n")

        self.fha.write("#ifdef REG_DF \n")
        for i in range(0,10):
            for j in range(0,6):
                self.fha.write("    LOCSTORE(store, %d, %d, STOREDIM, STOREDIM) = df.x_%d_%d;\n" % (j+4, i+10, j+4, i+10))
        self.fha.write("#endif \n") 

        # include print statements if debug option is on    
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            for i in range(0,10):
                for j in range(0,6):
                    self.fha.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOCSTORE(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "DF", j+4, i+10, j+4, i+10))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")


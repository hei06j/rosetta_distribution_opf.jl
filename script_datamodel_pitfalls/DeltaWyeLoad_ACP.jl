#!/usr/bin/env julia
###### AC-OPF using JuMP ######
#
# implementation reference: https://github.com/lanl-ansi/PowerModelsAnnex.jl/blob/master/src/model/ac-opf.jl
# only the built-in AD library is supported

import PowerModelsDistribution
import Ipopt
import JuMP
using rosetta_distribution_opf
import OpenDSSDirect
const PMD = PowerModelsDistribution
const ODSS = OpenDSSDirect


##
ODSS.dss("Redirect ./data/case2_load_3ph_delta_cz.dss")

vm = ODSS.Circuit.AllBusVolts()
cm = ODSS.PDElements.AllCurrentsAllCurrents()
vm_load_delta = vm[4:7]
cd_delta = -cm[5:8]

##
ODSS.dss("Redirect ./data/case2_load_3ph_wye_neutral_cz.dss")

vm = ODSS.Circuit.AllBusVolts()
cm = ODSS.PDElements.AllCurrentsAllCurrents()
vm_load_wye_neutral = vm[4:7]
cd_wye_neutral = -cm[5:8]

##
ODSS.dss("Redirect ./data/case2_load_3ph_wye_floating_cz.dss")

vm = ODSS.Circuit.AllBusVolts()
cm = ODSS.PDElements.AllCurrentsAllCurrents()
vm_load_wye_floating = vm[4:end]
cd_wye_floating = -cm[5:8]


##
ODSS.dss("Redirect ./data/case2_load_3ph_wye_ground_cz.dss")

vm = ODSS.Circuit.AllBusVolts()
cm = ODSS.PDElements.AllCurrentsAllCurrents()
vm_load_wye_grounded = vm[4:7]
cd_wye_grounded = -cm[5:8]

##

abs.([cd_delta cd_wye_neutral cd_wye_grounded cd_wye_floating])
angle.([cd_delta cd_wye_neutral cd_wye_grounded cd_wye_floating])'.*180/pi

abs.([[vm_load_delta ; 0] [vm_load_wye_neutral ; 0] [vm_load_wye_grounded ; 0] vm_load_wye_floating])
angle.([[vm_load_delta ; 0] [vm_load_wye_neutral ; 0] [vm_load_wye_grounded ; 0] vm_load_wye_floating]).*180/pi

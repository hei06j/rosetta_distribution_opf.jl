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
ODSS.dss("Redirect ./data/case3_4wire_kron.dss")

vm = ODSS.Circuit.AllBusVolts()
cm = ODSS.PDElements.AllCurrentsAllCurrents()
vm = vm[4:6]
cm = cm[4:6]

Z = 
[0.227218   0.0592176  0.0592176  0.0592176
 0.0592176  0.227218   0.0592176  0.0592176
 0.0592176  0.0592176  0.227218   0.0592176
 0.0592176  0.0592176  0.0592176  0.227218].+ im *
[0.879326  0.54107   0.475546  0.449141
 0.54107   0.879326  0.516533  0.475546
 0.475546  0.516533  0.879326  0.54107
 0.449141  0.475546  0.54107   0.879326]

In = -1/Z[4,4] * Z[4,1:3]' * cm
[abs.(cm) ; abs(In)]
[angle.(cm) ; angle(In)] * 180/pi


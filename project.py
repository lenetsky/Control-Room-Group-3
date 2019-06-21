#==================================================================
# Perform a MEBT quad (you choose) field gradient
# scan measuring the RMS transverse beam size with
# one Wire Scanner
# 
# Calculate the transverse RMS emittance and errors
# from this scan by using RMS vs. gradient curve and
# the general LSQ method
#==================================================================

import sys
import math
import types
import time
import random

from jarray import *
from java.lang import *
from java.util import *
from java.io import *

from Jama import Matrix

from xal.extension.widgets.plot import BasicGraphData, FunctionGraphsJPanel

from xal.smf.impl import Marker, Quadrupole, RfGap, BPM, RfCavity, ProfileMonitor
from xal.sim.scenario import Scenario, AlgorithmFactory, ProbeFactory

from xal.tools.beam.calc import CalculationsOnBeams

from xal.smf.data import XMLDataManager

#===============================================================
#              MAIN PROGRAM
#===============================================================

accl = XMLDataManager.loadDefaultAccelerator()
accSeq = accl.findSequence("MEBT")

wss = accSeq.getAllNodesOfType(ProfileMonitor.s_strType)

#--- Some of Wire Scanners can be removed from analysis
#wss = wss[0:3]

for ws in wss:
	print('quad: {:16s}  pos={:8.3f}m'.format(ws.getId(), accSeq.getPosition(ws)))

print "==========================="

quads = accSeq.getAllNodesOfType(Quadrupole.s_strType)
for quad in quads:
	print('quad: {:16s}  B={:8.3f}  pos={:8.3f}m'.format(quad.getId(), quad.getDfltField(), accSeq.getPosition(quad)))

print "==========================="

scenario = Scenario.newScenarioFor(accSeq)
scenario.setSynchronizationMode(Scenario.SYNC_MODE_DESIGN)

tracker = AlgorithmFactory.createEnvTrackerAdapt(accSeq)
tracker.setProbeUpdatePolicy(tracker.UPDATE_ALWAYS)
		
#---- initial kinetic energy in [MeV]

def reRun(quad, field_b, ws):
	'''
	quad/ws: selected quad/ws object
	field_b: quad field setpoint
	'''
	quad.setDfltField(field_b)
	print("Setting Quad0 to {}".format(field_b))

	probe = ProbeFactory.getEnvelopeProbe(accSeq,tracker)
#---- peak current in [A]
	peak_current = 0.000
	probe.setBeamCurrent(peak_current)

	eKin_init = probe.getKineticEnergy()/1.0e+6

	#scenario.setStopElementId(ws.getId())
	#scenario.setIncludeStopElement(False)

	scenario.setProbe(probe)
	scenario.resync()
	scenario.run()
	traj = scenario.getProbe().getTrajectory()

	#for quad in quads:
	state = traj.stateForElement(quad.getId())
	pos = state.getPosition()
	length = quad.getLength()
	print "quad=",quad.getId()," pos[m]= %8.3f "%pos, "len[m]=%8.3f"%length
	quad_mtrx = state.getResponseMatrix()  # Q

	#for ws in wss:
	state = traj.stateForElement(ws.getId())
	ws_mtrx = state.getResponseMatrix()  # R
	mtrx_s7 = ws_mtrx.times(quad_mtrx.inverse())  # S 7x7
	mtrx_s = [
		mtrx_s7.getElem(0,0), mtrx_s7.getElem(0,1),
		mtrx_s7.getElem(1,0), mtrx_s7.getElem(1,1)
	]
	# print(mtrx_s)
	pos = state.getPosition()
	print "ws=",ws.getId()," pos[m]= %8.3f "%pos

    #--------- sizes for each WS and transport matrices
	xRMS_Size = state.twissParameters()[0].getEnvelopeRadius()
	print(state.twissParameters()[0])
	#--------elements of the transport matrix
	a11 = ws_mtrx.getElem(0,0)
	a12 = ws_mtrx.getElem(0,1)
	m_row = [a11**2, 2*a11*a12,a12**2]
	
	return xRMS_Size, m_row, mtrx_s

matrx_all_S = []
matrx_all_X = []
sizes_X_arr = []
results = []
field_center = -85  # 1st
field0 = quads[0].getDfltField()
for pos in range(-15, 15):
	print('*'*40)
	field = field_center * (1 + pos/100.)
	xsize, mtrx_row, mtrx_s = reRun(quads[0], field, wss[0])
	matrx_all_X.append(mtrx_row)
	matrx_all_S.append(mtrx_s)
	sizes_X_arr.append(xsize)
	results.append((field, xsize, mtrx_s))

with open('data.txt', 'w') as f:
	for ix, (field, xsize, mtrx_s) in enumerate(results):
		print('Run: {:2d} Field: {:8f} xRMS size: {:12f}'.format(ix, field, xsize))
		f.write('{}    {}   {}\n'.format(field, xsize, mtrx_s[1]))  # S12

field_center = -28  # last
quads[0].setDfltField(field0)
matrx_all_S = []
matrx_all_X = []
sizes_X_arr = []
results = []
for pos in range(-15, 15):
	print('*'*40)
	field = field_center * (1 + pos/100.)
	xsize, mtrx_row, mtrx_s = reRun(quads[-1], field, wss[-1])
	matrx_all_X.append(mtrx_row)
	matrx_all_S.append(mtrx_s)
	sizes_X_arr.append(xsize)
	results.append((field, xsize, mtrx_s))

with open('data1.txt', 'w') as f:
	for ix, (field, xsize, mtrx_s) in enumerate(results):
		print('Run: {:2d} Field: {:8f} xRMS size: {:12f}'.format(ix, field, xsize))
		f.write('{}    {}   {}\n'.format(field, xsize, mtrx_s[1]))  # S12
#=========================================
#  Only x-axis analysis
#=========================================
sigma_rel_err = 0.05

n_ws = len(matrx_all_X)
mMatrX = Matrix(matrx_all_X,n_ws,3)

sigma2Vector = Matrix(n_ws,1)
weightM = Matrix.identity(n_ws,n_ws)

for ind in range(n_ws):
	sigma = sizes_X_arr[ind]*(1.0 +	random.gauss(0.,sigma_rel_err))
	sigma2Vector.set(ind,0,sigma**2)
	err2 = (2*sigma*sigma*sigma_rel_err)**2
	weightM.set(ind,ind,1.0/err2)
	
#=== mwmMatr = (M^T*W*M) =======
mwmMatr = ((mMatrX.transpose()).times(weightM)).times(mMatrX)

#=== corr2ValMatr = [(M^T*W*M)^-1] * M^T * W * Vector(sigma**2)
corr2ValMatr = (((mwmMatr.inverse()).times(mMatrX.transpose())).times(weightM)).times(sigma2Vector)

corr2ErrMatr = mwmMatr.inverse()

print "========================================"
print "<x^2>            = ","%12.5e"%corr2ValMatr.get(0,0)," +- ","%12.5e"%math.sqrt(abs(corr2ErrMatr.get(0,0)))
print "<x*x'>           = ","%12.5e"%corr2ValMatr.get(1,0)," +- ","%12.5e"%math.sqrt(abs(corr2ErrMatr.get(1,1)))
print "<x'^2>           = ","%12.5e"%corr2ValMatr.get(2,0)," +- ","%12.5e"%math.sqrt(abs(corr2ErrMatr.get(2,2)))
print "========================================"

x2 = corr2ValMatr.get(0,0)
xxp = corr2ValMatr.get(1,0)
xp2 = corr2ValMatr.get(2,0)

print "========================================"
print "x2, xxp, xp2 = ",   x2, xxp, xp2
print "========================================"

# Calculate some other parameters at the start of HEBT2
emitX = math.sqrt(x2*xp2-xxp*xxp)
alphX = -xxp/emitX
betaX = x2/emitX
gammaX = xp2/emitX
print "emittance =", emitX, ", apha = ", alphX, ", beta = ", betaX, ", gamma = ", gammaX

#----- emittance error calculation ----------------
sig2_x2 = corr2ErrMatr.get(0,0)
sig2_x_xp = corr2ErrMatr.get(1,1)
sig2_xp2 = corr2ErrMatr.get(2,2)

sig_x2_xp2 = corr2ErrMatr.get(0,2)
sig_x2_x_xp = corr2ErrMatr.get(0,1)
sig_xp2_x_xp = corr2ErrMatr.get(2,1)

sig2_emit2 = xp2**2*sig2_x2 + x2**2*sig2_xp2 + 4*xxp**2*sig2_x_xp + 2*x2*xp2*sig_x2_xp2 - 4*xxp*(xp2*sig_x2_x_xp + x2*sig_xp2_x_xp)
#print "sig2_emit2=",sig2_emit2
sig_emitt = 0.5*math.sqrt(abs(sig2_emit2))/emitX 
print "emitt [pi*mm*mrad] = %8.4f "%(emitX*1.0e+6)," +-  %8.4f "%(sig_emitt*1.0e+6)
print "========================================"
print "Curve fitting method:"
emitX_fit = 2.839e-6
sig_emitt_fit = 0.5*math.sqrt(abs(sig2_emit2))/emitX_fit
print "emitt [pi*mm*mrad] = %8.4f "%(emitX_fit*1e6)," +-  %8.4f "%(sig_emitt_fit*1.0e+6)

#------------------------------------------
#state_init = traj.initialState()
#x_rms_ini = state_init.twissParameters()[0].getEnvelopeRadius()
#
#twissX = state_init.twissParameters()[0]
#print "Initial model Twiss X=",twissX

print "Done."

import time
import java.awt.Color
import ij
from ij import IJ
from ij.plugin.frame import RoiManager
import ij.measure
import ij.gui.Plot
import ij.process.FloatProcessor
import sys
import math
import os.path
from org.apache.commons.math.analysis import MultivariateRealFunction
from jarray import *
from org.apache.commons.math.optimization import SimpleScalarValueChecker, GoalType
from org.apache.commons.math.optimization.direct import NelderMead
from org.apache.commons.math.optimization.fitting import CurveFitter

class Gaussian2D():
#	params[0]=center(x)
#	params[1]=center(y)
#	params[2]=sigma(x)
#	params[3]=sigma(y)
#	params[4]=amp
#	params[5]=yoffset
	__nx=None;
	__ny=None;
	__data=[];
	__params=None;

	def __init__(self, width, height, par):
		self.__nx=width
		self.__ny=height
		self.__params=par
		self.__data=[]

	def compute(self):
		c_x=float(self.__params[0])
		c_y=float(self.__params[1])
		s_x=float(self.__params[2])
		s_y=float(self.__params[3])
		amp=float(self.__params[4])
		y_off=float(self.__params[5])
		nx=self.__nx
		ny=self.__ny
		d=self.__data.append		#points to append method to save time on object lookup inside loop
		for j in range(self.__ny):
			for i in range(self.__nx):
				d(y_off+amp*math.exp(-((i-c_x)/s_x)**2 - ((j-c_y)/s_y)**2))
		return self.__data

	def plot2D(self):
		pix=array(self.__data, 'd')
		imp=ImagePlus("Gaussian2D", FloatProcessor(self.__nx, self.__ny, pix))
		imp.show()	    

	def getParams(self):
		return self.__params

	def getData(self):
		return self.__data

	def setParams(self, params):
		self.__params=params

class Gaussian2DMinimizer(MultivariateRealFunction):
	# for fitting, capable of dealing with constrained optimization problems
	__data=0.0
	__nx=0
	__ny=0
	__lbounds=[]
	__ubounds=[]
	__isConstrained=False

	def __init__(self):
		pass

	def setImage(self, data, width, height, lb=[], ub=[], isConstrained=False):
		self.__data=data
		self.__nx=int(width)
		self.__ny=int(height)
		self.__lbounds=lb
		self.__ubounds=ub
		self.__isConstrained=isConstrained

	def setBounds(self, lb=[],ub=[], isConst=False) :
		self.__lbounds=lb
		self.__ubounds=ub
		self.__isConstrainedFit=isConst

	def origToUnconstrainedParams(self, params) :
		#convert real-params to unconstrained parameters for fitting
		#assuming ub and lb are both specified
		lb=self.__lbounds
		ub=self.__ubounds
		out=[]
		for i in range(len(params)):
			if params[i]<=lb[i]:
				out.append(-math.pi/2)
			elif params[i]>=ub[i]:
				out.append(math.pi/2)
			else :
				 bla=2*(params[i]-lb[i])/float((ub[i]-lb[i]))-1
				 out.append(2*math.pi+math.asin(max(-1,min(1,bla))))
			#print "p[%f]: %f and lb,ub are (%f,%f) --> %f"  %  (i, params[i],lb[i],ub[i],out[i] )
		return out

  	def uncToOriginalParams(self, params) :
  		#convert unconstrained variables to real-parameter space
  		#for now, assume both lb and ub are specified
  		out=[]
  		for i in range(len(params)):
  			bla=float((math.sin(params[i])+1))/2
  			out.append(self.__lbounds[i]+(self.__ubounds[i]-self.__lbounds[i])*bla)
  		return out

	def value(self, params):	 #return residual value
		if self.__isConstrained :
			params=self.uncToOriginalParams(params)

		bla=Gaussian2D(self.__nx, self.__ny, params)
		calc=bla.compute()	    
		residual=0.0
		d=self.__data
		if len(calc)==len(self.__data) : #subtract them
			for i in range(len(calc)) :
				residual+=(calc[i]-d[i])**2
		else:
			#somehow calc and data vectors are not the same length
			raise

		return residual   


class GaussianFit:
	__width=0
	__height=0
	__realData=0.0
	__calcData=0.0
	__steps=[]
	__firstGuess=[]
	__bestGuess=[]
	__maxIter=1000
	__lbounds=[]
	__ubounds=[]
	__ip=[]
	__objSigma=3		#assume 3 pixel object width for first guess
	__fitResult=[]
	__convergenceVal=1e-4

	def reshapeArray(self,arr):
    #reshape 2DJava arrays from ImageJ for fitting in my code
		px=[]
		col =len(arr)
		row=len(arr[0])			#reshape list of pixels to match formatting for 2D Gaussian Fitting routines
		for i in range(row):
			for j in range(col):
				px.append(arr[j][i])
		return px


	def __init__(self, data=0.0, width=0, height=0):
		#constructor needs to check whether it got an ImagePlus, ImageProcessor or a data array
		self.__ip=[]
		pixels=[]
		if isinstance(data, ij.ImagePlus) : #check to see type of object we were sent
			data=data.getProcessor()  #if ImagePlus, just keep the processor

		if isinstance(data, ij.process.FloatProcessor):
			self.__ip=data
			pixels=self.__ip.getFloatArray()
		elif isinstance(data, ij.process.ByteProcessor) or isinstance(data, ij.process.ShortProcessor):
			self.__ip=data
			pixels=self.__ip.getIntArray()
		elif isinstance(data, list):  #just got a list of values
			self.__realData=data
			self.__width=int(width)
			self.__height=int(height)
		else:
			sys.exit('Gaussian Fitter could not initialize because data was not recognized')

		#assign width, height and pixel values based on processor, if one was specified
		if self.__ip!=[]:
			self.__width=self.__ip.getWidth()
			self.__height=self.__ip.getHeight()
			self.__realData=self.reshapeArray(pixels)

	def getBounds(self):
		return (self.__lbounds, self.__ubounds)

	def getData(self):
		return self.__realData

	def getFitParams(self):
		pass

	def getWidth(self):
		return self.__width

	def getHeight(self):
		return self.__height

	def getSteps(self):
		return self.__steps

	def getFirstGuess(self):
		return self.__firstGuess

	def getBestGuess(self):
		return self.__bestGuess

	def getResidual(self, p=[]):
		if p==[]:
			p=self.__bestGuess
		calcMin=Gaussian2DMinimizer()
		calcMin.setImage(self.__realData, self.__width, self.__height)
		return calcMin.value(p)


  	def getConvergenceValue(self):
  		return self.__convergenceVal

	def estimateGaussianParametersCOM(self):
	#try to estimate first guess for fitting a Gaussian based on center of mass
		if self.__ip!=[] :
			bla=self.__ip.getStatistics()

			return [bla.xCenterOfMass, bla.yCenterOfMass, self.__objSigma, self.__objSigma, bla.max, bla.min]		
		else:
			sys.error('Could not estimate centroid, ImageProcessor not defined!')

	def estimateGaussianParametersBinary(self):
	#try to estimate first guess for fitting a Gaussian based on largest binary particle
		if self.__ip!=[] :
			ip_binary=self.__ip.duplicate()
			ip_binary.setAutoThreshold(ij.process.AutoThresholder.Method.IsoData, True) 
			#ip_binary.setAutoThreshold(ImageProcessor.ISODATA2, True) 
			ip_binary.threshold(int(ip_binary.getMinThreshold()))
			imp=ImagePlus('temp',ip_binary)
			#not finished yet
			print 'Threshold: ', int(ip_binary.getMinThreshold())

			#return [bla.xCenterOfMass, bla.yCenterOfMass, self.__objSigma, self.__objSigma, bla.max, bla.min]		
		else:
			sys.error('Could not estimate centroid, ImageProcessor not defined!')

	def setFirstGuess(self, p=[]) :	
		if p!=[] :
			self.__firstGuess=p
		else:
			pass #eventually have code to assemble first guess

	def setBestGuess(self, p=[]) :	
		if p!=[] :
			self.__bestGuess=p
		else:
			pass #eventually have code to assemble first guess

	def setSteps(self, mutiplier=.05, offset=0.05):
		self.__steps=[]
		p=self.__firstGuess
		for i in range(len(p)) :	    #set up array of initial steps for simplex
			if p[i] !=0 :
				self.__steps.append(p[i]*mutiplier)
			else :
				self.__steps.append(offset)

	def setBounds(self, lb, ub):
  		self.__lbounds=lb
  		self.__ubounds=ub

  	def setConvergenceValue(self,val):
  		self.__convergenceVal=val

	def doFit(self, maxIter=10000):
		self.__maxIter=abs(maxIter)
		#set up Gaussain2DMinimizer function, include options for constrained fitting
		#constrained optimization is done as in fminsearchbnd.m MATLAB routine (google it)
		calcMin=Gaussian2DMinimizer()

		calcMin.setBounds(self.__lbounds,self.__ubounds)


		#check to see if there are constraints
		if len(self.__lbounds)==len(self.__firstGuess) and len(self.__ubounds)==len(self.__firstGuess):
			calcMin.setImage(self.__realData, self.__width, self.__height, self.__lbounds,self.__ubounds ,True)
			#printGuessBounds(self.__firstGuess, self.__lbounds, self.__ubounds)

			self.__firstGuess=calcMin.origToUnconstrainedParams(self.__firstGuess)
			#print "Constrained Fit"
			isConstrained=True
		else:
			calcMin.setImage(self.__realData, self.__width, self.__height, [],[] ,False)
			#print "UN Constrained Fit"
			isConstrained=False

		self.setSteps()
		n=array(self.__firstGuess, 'd')
		#set up NelderMead optimizer settings
		nMeadOptimizer=NelderMead()
		convergeChecker=SimpleScalarValueChecker(self.__convergenceVal,-1)	 #play with this later to speed up fitting?
		nMeadOptimizer.setStartConfiguration(array(self.__steps, 'd'))  #define step size in simplex
		nMeadOptimizer.setConvergenceChecker(convergeChecker)
		nMeadOptimizer.setMaxIterations(self.__maxIter);

		#t=time.time()
		#begin fitting at this point
		try:  #catch error if fit fails
			out=nMeadOptimizer.optimize(calcMin, GoalType.MINIMIZE, n)
		except:
			self.__fitResult=None
		else:
			#print "Number of iterations: %f" % nMeadOptimizer.getIterations()
			#print "Time: %.3f" % (time.time()-t)
		#print out.getPoint().tolist()
			if isConstrained:
				self.__fitResult= calcMin.uncToOriginalParams(out.getPoint().tolist())
			else:
				self.__fitResult= out.getPoint().tolist()

		return self.__fitResult  #return None if fit failed


	def setDataFromImg(self, img): #sets fitting data from ImagePlus
		ip=img
		if isinstance(img, ij.ImagePlus) : #check to see type of object we were sent
			ip=img.getProcessor()
		else :
			pass	 #error of some sort

		if isinstance(ip, ij.process.FloatProcessor):
			bla=ip.getFloatArray()
		elif isinstance(ip, ij.process.ByteProcessor) or isinstance(ip, ij.process.ShortProcessor):
			bla=ip.getIntArray()
		else:
			print 'weird format'
			raise

		px=[]
		col=len(bla)
		row=len(bla[0])			#reshape list of pixels to match formatting for 2D Gaussian Fitting routines
		for i in range(row):
			for j in range(col):
				px.append(bla[j][i])
		self.__realData=px


#	params[0]=center(x)
#	params[1]=center(y)
#	params[2]=sigma(x)
#	params[3]=sigma(y)
#	params[4]=amp
#	params[5]=yoffset

class StackFitter:
	#stack fitting methods
	__imp=[]
	__nSlices=[]
	__width=[]
	__height=[]
	__prevFitParams=[]
	__initialGuess=[]
	__fitResults={}
	__sliceFitter=[]
	__keyFrames={}
	__objSigma=3		#assume 3 pixel object width for first guess
	headings=["x_c","y_c","x_width","y_width","Amplitude","offset"]

	def __init__(self, imp=[]):
		self.__imp=imp
		if isinstance(imp, ij.ImagePlus):
			self.__width=imp.getWidth()
			self.__height=imp.getHeight()
			self.__nSlices=imp.getImageStackSize()
			self.__sliceFitter=GaussianFit(imp.getStack().getProcessor(1))
		else:
			sys.exit('Cant fit stack because ImagePlus was not specified')

	def getInitialGuess(self):
		return self.__initialGuess

	def getFitResults(self):
		return self.__fitResults

	def printFitResults(self):
		for fr in sorted(self.__fitResults):
			print 'Frame: ', fr, 'Fit: ', self.__fitResults[fr]

	def showFitResultsTable(self):
		bla=ij.measure.ResultsTable() #new results table
		#self.__fitResults
		#self.printFitResults()

		for fr in sorted(self.__fitResults):
			bla.incrementCounter()
			bla.addValue('frame',fr)
			for i in range(6):
				bla.addValue(''.join(['LB ',self.headings[i]]), self.__fitResults[fr][1][i])
				bla.addValue(self.headings[i], self.__fitResults[fr][0][i])
				bla.addValue(''.join(['UB ',self.headings[i]]), self.__fitResults[fr][2][i])

			bla.addValue("Residual", self.__fitResults[fr][3])
		bla.show("Fit Results")

	def showFitResultsPlotWindow(self, fitColumn=0):
		#make lists of dictionary
		a=[]
		frames=sorted(self.__fitResults)
		for i in frames:
  			a.append(self.__fitResults[i][0][fitColumn])
		bla=ij.gui.Plot(''.join(['Fit Result for: ',self.headings[fitColumn]]), "Frame #", self.headings[fitColumn], array(frames,'f'),array(a,'f')) #new results table
		bla.show()

	def showFitResultsRGBImg(self):
		#create a new RGB stack with the fit oval marked in one color and original data in another
		impResult=ij.ImagePlus("Results",self.__imp.getStack())
		IJ.run(impResult, "RGB Color","")
		frames=sorted(self.__fitResults)
		for fr in frames:
			fitout=self.__fitResults[fr][0]
			ip=	impResult.getStack().getProcessor(fr)
			ip.setColor(java.awt.Color.MAGENTA)
			ip.drawOval(int(round(fitout[0]-fitout[2],0)),int(round(fitout[1]-fitout[3],0)),int(round(fitout[2]*2,0)),int(round(fitout[3]*2,0)))
		impResult.show()
		return impResult

	def getKeyFrames(self):
		return self.__keyFrames

	def setKeyFramesMultiROI(self):
		bla=RoiManager(True)
		rm=bla.getInstance()
		if rm!=None and rm.getCount()>0:
			roiArr=rm.getRoisAsArray()
			for i in range(len(roiArr)):
				name=rm.getName("%s"%i) 
				fr=rm.getSliceNumber(name)
				if fr>0:
					self.__keyFrames[fr]=[roiArr[i].getPolygon().xpoints[0], roiArr[i].getPolygon().ypoints[0]]

	def fitAllSlices(self):
		self.__fitResults={}
		lb=[]
		ub=[]
		params=[]
		moveSpeed=2 # max pixels allowed to move from frame to frame
		maxIntensityFluct=.75 #75% percentage intensity fluctuations allowed
		maxSigmaFluct=0.75 #75% percentage width fluctuations allowed
		self.__initialGuess=[]
		for i in range(1,self.__nSlices+1) :
			if i in self.__keyFrames:  #check to see if keyFrame is defined for this slice, if so, constrained fit
				print "Key frame defined in slice %i -- x: %i, y: %i" % (i, self.__keyFrames[i][0], self.__keyFrames[i][1])
				bla=self.__imp.getImageStack().getProcessor(i).getStatistics()
				#print bla
				self.__initialGuess=[int(self.__keyFrames[i][0]), int(self.__keyFrames[i][1]), self.__objSigma, self.__objSigma, int(bla.max), int(bla.min)]	

			params,residual=self.fitSlice(i,lb,ub)  #do the fit
			#if fit is reasonable, make these parameters first guess for next frame, also adjust lb,ub
			print "Frame %i of %i (%f done)" % (i, self.__nSlices+1,float(i)/(self.__nSlices+1))
			if params[2]<10 and params[2]>1 and params[3]<10 and params[3]>1 and params[4]>3*params[5]:
				self.__initialGuess=params #set intial guess for next frame
				#center position constrains
				lb=[0,0,0,0,0,0]
				ub=[0,0,0,0,0,0]
				lb[0]=max(params[0]-moveSpeed,0)
				lb[1]=max(params[1]-moveSpeed,0)
				ub[0]=min(params[0]+moveSpeed,self.__width)
				ub[1]=min(params[1]+moveSpeed,self.__height)

				#sigma constrains
				lb[2]=params[2]*maxSigmaFluct
				lb[3]=params[2]*maxSigmaFluct
				ub[2]=params[2]*(maxSigmaFluct+1)
				ub[3]=params[2]*(maxSigmaFluct+1)

				#amplitude
				lb[4]=params[4]*maxIntensityFluct
				ub[4]=params[4]*(maxIntensityFluct+1)

				#offset
				lb[5]=0
				ub[5]=params[4]
				self.__fitResults[i]=[params, lb,ub,residual]
			else: #fit was unreasobable, maybe object blinked, for now do unconstrained fit
				lb=[]
				ub=[]
				#self.__initialGuess=[]
				print "Unconstrained fit at frame: %i" % i
			#end if	


	def fitSlice(self, sl=1, lb=[], ub=[], numIter=10000):
	#fits 2D gaussian to slice in stack
		gFit=GaussianFit(self.__imp.getImageStack().getProcessor(sl))
		if self.__initialGuess==[]:  #check to see if initial guess was defined
			#print gFit.estimateGaussianParametersBinary()
			self.__initialGuess=gFit.estimateGaussianParametersCOM()

		if len(lb)==len(self.__initialGuess) and len(ub)==len(self.__initialGuess): #bounded fitting 
			#printGuessBounds(self.__initialGuess,lb,ub)
			gFit.setBounds(lb,ub)

		gFit.setFirstGuess(self.__initialGuess) 
		#lb=[1,1,1,1,0.25*guess[4],0.5*guess[5]]
		#ub=[w-1,h-1,w,h,2*guess[4],0.5*guess[4]]
		#test.setBounds(lb,ub)
		bestVal=gFit.doFit(numIter)
		return (bestVal,gFit.getResidual(bestVal))

	def setInitialGuess(self, params):
		self.__initialGuess=params

	def writeFittingResults(self, filename=""):
		if filename==None or filename=="":  #try to place filename in same folder as imp
			bla=self.__imp.getOriginalFileInfo()  #determine info on original file
			try:	#see if the file has been saved already
				dirpath= bla.directory #determine where it was saved
				fn,dummy=os.path.splitext(bla.fileName)
				fn+="_fit.txt"
				filename=os.path.join(dirpath, fn)
				print "Data saved to: "+filename
			except:
				IJ.showMessage("Save stack to file or define data export file name!")
				raise RuntimeException("Save stack to file")

		#now write to file		
		try:
			f= open(filename, "w" )
			#write header row

			for fr in sorted(self.__fitResults):
				st="%s" % fr
				for i in range(6):
					st+="\t%s"%self.__fitResults[fr][0][i] #saved just the fit values, not the bounds
				f.write(st+'\n')
		finally:
			f.close()

	def writeFitResultsRGB(self, impResult, filename=""):
		if filename==None or filename=="":  #try to place filename in same folder as imp
			bla=self.__imp.getOriginalFileInfo()  #determine info on original file
			try:	#see if the file has been saved already
				dirpath= bla.directory #determine where it was saved
				fn,dummy=os.path.splitext(bla.fileName)
				fn+="_Results.tif"
				filename=os.path.join(dirpath, fn)

				print "Results saved to: "+filename
			except:
				IJ.showMessage("Save stack to file or define data export file name!")
				raise RuntimeException("Save stack to file")

		#now write to file		
		IJ.save(impResult, filename)

def printGuessBounds(guess, lb, ub):
	for i in range(len(guess)):
		print "Parameter [%i]: %f < %f < %f" % (i, lb[i],guess[i], ub[i])

if __name__ == "__main__":

	w=10
	h=20
	p=[5, 8, 2.5,3.5,.7,0]
	p2=[5,5,5,5,10,0]
	maxIter=10000

	#bla=Gaussian2D(w,h,p)
	#bla.compute()
	#bla.plot2D()
	imp=IJ.getImage()
	testingStack=True
	if testingStack:
		test=StackFitter(imp)
		test.setKeyFramesMultiROI()
		print test.getKeyFrames()
		test.fitAllSlices()
		test.writeFittingResults()
		#out=test.getFitResults()
		impResult=test.showFitResultsRGBImg()
		test.writeFitResultsRGB(impResult)
		#test.showFitResultsTable()
		test.showFitResultsPlotWindow(1)
	else:
		test=GaussianFit(imp.getProcessor())
		guess=test.estimateGaussianParametersCOM()
		test.setFirstGuess(guess) 
		lb=[1,1,1,1,0.25*guess[4],0.5*guess[5]]
		ub=[w-1,h-1,w,h,2*guess[4],0.5*guess[4]]

		guess[1]=10
		printGuessBounds(guess, lb, ub)
		print "Residual from first guess: %f" % test.getResidual(guess)
	#lb=[1,1,1,1,1,0]
	#ub=[6,6,6,6,1,6]
		test.setBounds(lb,ub)
		bestVal=test.doFit(10000)
		print "Residual from fit result: %f" % test.getResidual(bestVal)
#	print test.uncToOriginalParams([7.268, 7.0129, 1.57, 1.57, 1.57,1.56])
#	print test.doFit()
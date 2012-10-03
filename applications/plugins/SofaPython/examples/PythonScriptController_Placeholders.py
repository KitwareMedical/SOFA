import Sofa

############################################################################################
# this is a PythonScriptController empty script, with entry points placeholders
# following defs are optionnal entry points, called by the PythonScriptController component;
############################################################################################

class ExampleController(Sofa.PythonScriptController):

	# called once the script is loaded
	def onLoaded(self,node):
		print 'Controller script loaded from node %s'%node.findData('name').value
		return 0

	# optionnally, script can create a graph...
	def createGraph(self,node):
		print 'createGraph called (python side)'
		return 0



	# called once graph is created, to init some stuff...
	def initGraph(self,node):
		print 'initGraph called (python side)'
		return 0



	# called on each animation step
	total_time = 0
	def onBeginAnimationStep(self,dt):
		self.total_time += dt
		#print 'onBeginAnimatinStep (python) dt=%f total time=%f'%(dt,self.total_time)
		return 0
	 
	def onEndAnimationStep(self,dt):
		return 0
	 
	# called when necessary by Sofa framework... 
	def storeResetState(self):
		print 'storeResetState called (python side)'
		return 0

	def reset(self):
		print 'reset called (python side)'
		return 0

	def cleanup(self):
		print 'cleanup called (python side)'
		return 0
	 
	 
	# called when a GUIEvent is received
	def onGUIEvent(self,controlID,valueName,value):
		print 'GUIEvent received: controldID='+controlID+' valueName='+valueName+' value='+value
		return 0 
	 
	 
	 
	# key and mouse events; use this to add some user interaction to your scripts 
	def onKeyPressed(self,k):
		print 'onKeyPressed '+k
		return 0 
	 
	def onKeyReleased(self,k):
		print 'onKeyReleased '+k
		return 0 

	def onMouseButtonLeft(self,x,y,pressed):
		print 'onMouseButtonLeft x='+str(x)+' y='+str(y)+' pressed='+str(pressed)
		return 0
	 
	def onMouseButtonRight(self,x,y,pressed):
		print 'onMouseButtonRight x='+str(x)+' y='+str(y)+' pressed='+str(pressed)
		return 0
	 
	def onMouseButtonMiddle(self,x,y,pressed):
		print 'onMouseButtonMiddle x='+str(x)+' y='+str(y)+' pressed='+str(pressed)
		return 0
	 
	def onMouseWheel(self,x,y,delta):
		print 'onMouseButtonWheel x='+str(x)+' y='+str(y)+' delta='+str(delta)
		return 0
	 
	def onScriptEvent(self,senderNode,eventName,data)
		print 'onScriptEvent eventName='+eventName+' data='+str(data)
		return 0
	 
 

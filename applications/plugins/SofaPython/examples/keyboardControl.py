import Sofa

class KeyboardControl(Sofa.PythonScriptController):
	# called once graph is created, to init some stuff...
	def initGraph(self,node):
		print 'initGraph called (python side)'
		self.MechanicalState = node.getObject('DOFs')
		return 0
	 
	 
	# key and mouse events; use this to add some user interaction to your scripts 
	def onKeyPressed(self,k):
		
		# free_position is a scalar array : [tx,ty,tz,rx,ry,rz,rw]
		free_position=self.MechanicalState.findData('free_position').value
		
		# translation speed
		speed = 1 

		# UP key : front
		if ord(k)==19:
			free_position[2]-=speed
		# DOWN key : rear
		if ord(k)==21:
			free_position[2]+=speed
		# LEFT key : left
		if ord(k)==18:
			free_position[1]-=speed
		# RIGHT key : right
		if ord(k)==20:
			free_position[1]+=speed
		# PAGEUP key : up
		if ord(k)==22:
			free_position[0]-=speed
		# PAGEDN key : down
		if ord(k)==23:
			free_position[0]+=speed
			
		self.MechanicalState.findData('free_position').value=free_position
		return 0 
	 
	 

import os

# concatenate lists for use with data
def cat(x):
    return ' '.join(map(str, x))


# absolute path from filename

def path( name ):
    return os.path.dirname( os.path.abspath( name ) )


# reasonable standard scene
def scene(node):
    
    node.createObject('RequiredPlugin', pluginName = "Compliant" )
    node.createObject('RequiredPlugin', pluginName = "CompliantDev" )

    node.dt = 0.01
    node.gravity = '0 -9.81 0'

    node.createObject('DefaultPipeline', name = 'pipeline')
    node.createObject('BruteForceDetection', name = 'detection')
    
    proximity = node.createObject('NewProximityIntersection', name = 'proximity' )
    
    manager = node.createObject('DefaultContactManager',
                                name = 'manager',
                                response = "FrictionCompliantContact",
                                responseParams = "mu=1" )
    
    style = node.createObject('VisualStyle', 
                              name = 'style',
                              displayFlags='hideBehaviorModels hideCollisionModels hideMappings hideForceFields')

    ode = node.createObject('AssembledSolver',
                            name='ode' )
    
    return node.createChild('scene') 


import numpy as np

import h5py


from xopt.tools import fingerprint

VOCS = {
    'name':'TNK_test',
    'description':'Constrainted test function TNK. See Table V in https://www.iitk.ac.in/kangal/Deb_NSGA-II.pdf',
    'simulation':'test_TNK',
    'variables': {
        'x1':[0, 3.14159],
        'x2':[0, 3.14159]
    },
    'objectives':{
        'y1':'MINIMIZE',
        'y2':'MINIMIZE'
        
    },
    'constraints':{
        'c1': ['GREATER_THAN', 0],
        'c2': ['GREATER_THAN', 0]
        
    },
    'constants':{'a':'dummy_constant'},
    'linked_variables':{'x9':'x1'}
    
    
}


NAME = 'TNK'
BOUND_LOW, BOUND_UP = [0.0, 0.0], [3.14159, 3.14159]  


X_RANGE = [0, 1.4]
Y_RANGE = [0, 1.4]

# Pure number version
def TNK(individual):    
    x1=individual[0]
    x2=individual[1]
    objectives =  (x1, x2)
    constraints = (x1**2+x2**2-1.0 - 0.1*np.cos(16*np.arctan2(x1, x2)), 0.5-(x1-0.5)**2-(x2-0.5)**2 )
    return objectives, constraints

# labeled version
def evaluate_TNK(inputs, extra_option='abc', **params):
    
    info = {'some':'info', 'about':['the', 'run']}
    ind = [inputs['x1'], inputs['x2']]
    
    if ind[0]>BOUND_UP[0]:
        raise ValueError(f'Input greater than {BOUND_UP[0]} ')
        
    name = fingerprint(inputs)
    fname = name +'.h5'
    
    objectives, constraints = TNK(ind)    
    outputs = {'y1':objectives[0], 'y2':objectives[1], 'c1':constraints[0], 'c2':constraints[1],'archive':fname}
    archive_tnk_outputs(inputs,outputs,name=fname)
    
    return outputs


def archive_tnk_outputs(inputs,outputs,name=None):
    dat={'inputs':inputs, 'outputs':outputs}
    
    #if not name:
    #    name = fingerprint(outputs)
    #fname=name+'.h5'
        
    with h5py.File(name, 'w') as h5:
        
        g = h5.create_group('inputs')
        for k,v in inputs.items():
            g.attrs[k]=v
            
        g = h5.create_group('outputs')
        for k,v in outputs.items():
            g.attrs[k]=v
        
    return
        
    
    





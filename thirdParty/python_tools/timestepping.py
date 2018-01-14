#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Contains different timestepping schemes

AUTHOR  - Jarrad Kapor (jarrad.kapor@postgrad.curtin.edu.au)
MODIFIED- 2009.11

"""

'''def ode_calc_dx(step_methd, dt, x_zero , y_zero, velocity_function, vel_fun_dict, step_option_dict=None)
    """#Generic timestepping module...
    
    #velocity_function : takes dictionary and ex,ey, returns velocity field at those points
    """
    (dx_dt, dy_dt) = velocity_function(vel_fun_dict, x_zero, y_zero)
    
    if step_methd == 'euler':
        dx = [i * dt for i in dx_dt]
        dy = [i * dt for i in dy_dt]

    elif step_methd == 'adams-bashforth-one':
        dx = [0.0] * len(x_zero)
        dy = [0.0] * len(x_zero)
        if step_option_dict['old_dx_dt'] is not None and step_option_dict['old_dy_dt'] is not None:
            old_dx_dt = step_option_dict['old_dx_dt'] 
            old_dy_dt = step_option_dict['old_dy_dt'] 
            for i in xrange(len(x_zero)):
                if old_dx_dt[i] is not None and old_dy_dt[i] is not None:
                    thisdx = dt * (1.5 * dx_dt[i] - 0.5 * old_dx_dt[i])
                    thisdy = dt * (1.5 * dy_dt[i] - 0.5 * old_dy_dt[i])
                else:
                    # not sufficient old field data, so use euler instead

                
                dx[i] = thisdx
                dy[i] = thisdy
        else:
            # fall back to euler for the first timestep...
            dx = [i * dt for i in dx_dt]
            dy = [i * dt for i in dy_dt]

    else:
        raise ValueError("Given DVM time stepping method \"" + str(TIME_STEPPING) + "\" does not exist!")
'''





def euler_dx(dt, x_zero , y_zero, velocity_function, vel_fun_dict):
    """Generic euler timestepping module...
    
    velocity_function : takes dictionary and ex,ey, returns velocity field at those points
    """
    (dx_dt, dy_dt) = velocity_function(x_zero, y_zero, vel_fun_dict)
    dx = [i * dt for i in dx_dt]
    dy = [i * dt for i in dy_dt]
    return dx, dy

def adams_bashforth_dx(dt, x_zero , y_zero, velocity_function, vel_fun_dict, old_dx_dt=None, old_dy_dt=None, backup=None, mod_dict_fun=None):
    """Generic timestepping module...
    
    velocity_function : takes dictionary and ex,ey, returns velocity field at those points
    """
    
    (dx_dt, dy_dt) = velocity_function(x_zero, y_zero, vel_fun_dict)
    
    dx = [0.0] * len(x_zero)
    dy = [0.0] * len(x_zero)
    leftover_inds = []
    if old_dx_dt is not None and old_dy_dt is not None:
        for i in xrange(len(x_zero)):
            if old_dx_dt[i] is not None and old_dy_dt[i] is not None:
                dx[i] = dt * (1.5 * dx_dt[i] - 0.5 * old_dx_dt[i])
                dy[i] = dt * (1.5 * dy_dt[i] - 0.5 * old_dy_dt[i])
            else:
                # not sufficient old field data, so use backup instead
                leftover_inds.append(i)
    else:
        # keep track of points not determined through adams-b!
        leftover_inds = range(len(x_zero))
    
    # dont forget about leftovers
    if len(leftover_inds) > 0:
        # print "len leftovers", len(leftover_inds)
        leftover_x = [x_zero[i] for i in leftover_inds]
        leftover_y = [y_zero[i] for i in leftover_inds] 
        leftover_dx_dt = [dx_dt[i] for i in leftover_inds]
        leftover_dy_dt = [dy_dt[i] for i in leftover_inds] 

        # fall back to euler for the first timestep...
        if backup == 'euler' or backup is None:
            left_real_dx = [i * dt for i in leftover_dx_dt]
            left_real_dy = [i * dt for i in leftover_dy_dt]
            left_real_dx_dt = leftover_dx_dt
            left_real_dy_dt = leftover_dy_dt
        elif backup == 'improved-euler' or backup == 'improved-euler1' :
            
            if backup == 'impro-euler1':
                pcmethod = 'hybrid'
                print 'using hybrid p-c scheme'
            else:
                pcmethod = 'standard'
            if mod_dict_fun is None:
                raise TypeError("You asked for a method which does intermediate steps with the dictionary however never provided a function to perform those intermediate modifications of the intermediate positions!")
            
            numtargets = len(leftover_inds)
            # predictor
            # already calculated first velocity
             
            left_new_x = [leftover_x[i] + leftover_dx_dt[i] * dt for i in xrange(numtargets)]
            left_new_y = [leftover_y[i] + leftover_dy_dt[i] * dt for i in xrange(numtargets)]
            # reconstruct new_x...then evaluate at the new_x
            if pcmethod == 'hybrid':
                # uses bashforth to predict some particles motion...not really consistant?
                new_x = [x_zero[i] + dx[i] for i in xrange(len(x_zero))]
                new_y = [y_zero[i] + dy[i] for i in xrange(len(x_zero))]
                for i in xrange(len(leftover_inds)):
                    new_x[leftover_inds[i]] = left_new_x[i]
                    new_y[leftover_inds[i]] = left_new_y[i]
            else:
                new_x = [x_zero[i] + dt * dx_dt[i] for i in xrange(len(x_zero))]
                new_y = [y_zero[i] + dt * dy_dt[i] for i in xrange(len(x_zero))]
            
            # need to modify dictionary to new aproximate locations
            mod_dict_fun(vel_fun_dict, new_x, new_y)
            # corrector
            (new_dx_dt, new_dy_dt) = velocity_function(left_new_x, left_new_y, vel_fun_dict )
            
            left_real_dx_dt = [(leftover_dx_dt[i] + new_dx_dt[i]) / 2.0  for i in xrange(numtargets)]
            left_real_dy_dt = [(leftover_dy_dt[i] + new_dy_dt[i]) / 2.0  for i in xrange(numtargets)]
            left_real_dx = [dt * i for i in left_real_dx_dt]
            left_real_dy = [dt * i for i in left_real_dy_dt]
                
            # need to fix up damage done
            mod_dict_fun(vel_fun_dict, x_zero, y_zero)

        else:
            raise ValueError("Given Adams-bashforth backup stepping method \"" + str(backup) + "\" does not exist!")

        notdodgeo = False
        # reconstruct leftovers!
        for i in xrange(len(leftover_inds)):
            dx[leftover_inds[i]] = left_real_dx[i]
            dy[leftover_inds[i]] = left_real_dy[i]
            if notdodgeo:
                dx_dt[leftover_inds[i]] = left_real_dx_dt[i]
                dy_dt[leftover_inds[i]] = left_real_dy_dt[i]


    return dx, dy, dx_dt, dy_dt


def runge_kutta_dx(dt, x_zero , y_zero, velocity_function, vel_fun_dict, mod_dict_fun=None):
    """mod is the field that gets changed in the space of a substep.
    
    modify function takes in x,y,disctionary
    """
    raise ValueError("Runge-kutta is not implemented PROPERLY yet!!!")
    order = 2
    subdt = dt / order
    numtargets = len(x_zero)
    dxlog = []
    dylog = []
    new_x = x_zero
    new_y = y_zero
    newdx = [0.0] * numtargets
    newdy = [0.0] * numtargets

    for i in xrange(order):
        newdx, newdy = euler_dx(subdt, new_x, new_y, velocity_function, vel_fun_dict)
        new_x = [new_x[i] + newdx[i] for i in xrange(numtargets)]
        new_y = [new_y[i] + newdy[i] for i in xrange(numtargets)]
        if mod_dict_fun is not None: 
            # need to fix up damage done
            mod_dict_fun(vel_fun_dict, new_x, new_y)
        
        dxlog.append(newdx)
        dylog.append(newdy)
    if order == 2:
        dxlog1 = [0.0] * numtargets
        dxlog2 = [0.0] * numtargets
        dxlog = [dxlog[0], dxlog1, dxlog2, dxlog[1]]
        dylog = [dylog[0], dxlog1, dxlog2, dylog[1]]
        order = 6

    dx = [order * (dxlog[0][i] + dxlog[3][i] + 2 * (dxlog[1][i]  + dxlog[2][i])) / 6.0 for i in xrange(numtargets)]
    dy = [order * (dylog[0][i] + dylog[3][i] + 2 * (dylog[1][i]  + dylog[2][i])) / 6.0 for i in xrange(numtargets)]
    
    if mod_dict_fun is not None: 
        # need to fix up damage done
        mod_dict_fun(vel_fun_dict, x_zero, y_zero)

    return dx, dy

def improved_euler_dx(dt, x_zero , y_zero, velocity_function, vel_fun_dict, mod_dict_fun):
    """mod is the field that gets changed in the space of a substep.
    
    modify function takes in x,y,disctionary
    
    explaination on huanns method
    http://calculuslab.deltacollege.edu/ODE/7-C-2/7-C-2-h.html
    """
    numtargets = len(x_zero)
    # predictor
    (dx_dt, dy_dt) = velocity_function(x_zero, y_zero, vel_fun_dict)
    new_x = [x_zero[i] + dx_dt[i] * dt for i in xrange(numtargets)]
    new_y = [y_zero[i] + dy_dt[i] * dt for i in xrange(numtargets)]
    
    # need to modify dictionary to new aproximate locations
    mod_dict_fun(vel_fun_dict, new_x, new_y)
    # corrector
    (new_dx_dt, new_dy_dt) = velocity_function(new_x, new_y, vel_fun_dict)
    
    dx = [dt * (dx_dt[i] + new_dx_dt[i]) / 2.0  for i in xrange(numtargets)]
    dy = [dt * (dy_dt[i] + new_dy_dt[i]) / 2.0  for i in xrange(numtargets)]
        
    # need to fix up damage done
    mod_dict_fun(vel_fun_dict, x_zero, y_zero)

    return dx, dy

def runge_kutta_dx_mod(dt, x_zero , y_zero, velocity_function, vel_fun_dict):
    raise ValueError("Runge-kutta not implemented yet!!!")
    return dx, dy

def test_main():
    # tests for euler
    velfun = lambda dct,x,y: ([dict['test'] /2 for i in range(len(x))], [dict['test'] * 0.25 for i in range(len(x))])
    modfun = lambda dct,x,y: None
    x=range(5)
    y=range(5)
    dict={'test':1.0}
    dt = 0.1
    edx,edy = euler_dx(dt, x, y, velfun, dict)
    print "\n---------------------------------"
    print "euler dx,dy"
    print edx
    print edy
    print "EULER VALIDATED OK ON 2009/11/24"
    """ rdx,rdy = runge_kutta_dx(dt, x, y, velfun, dict)
    print "\n---------------------------------"
    print "runge dx,dy"
    print rdx
    print rdy
    print "runge not fully VALIDATED OK ON 2009/11/30"
    """
    hdx,hdy = improved_euler_dx(dt, x, y, velfun, dict, modfun)
    print "\n---------------------------------"
    print "runge dx,dy"
    print hdx
    print hdy
    print "improved euler not fully VALIDATED OK ON 2009/11/30"
    dx,dy,u,v = adams_bashforth_dx(dt, x, y, velfun, dict)
    print "\n---------------------------------"
    print "no preprovided u,v should equal same as euler "
    print dx
    print dy
    print "\n---------------------------------"
    print u
    print v
    print "with all preprovided u/v"
    ou = [1.0] * len(dx)
    ov = [1.0] * len(dx)
    dx,dy,u,v = adams_bashforth_dx(dt, x, y, velfun, dict, ou, ov)
    print dx
    print dy
    print "\n---------------------------------"
    print u
    print v
    ou[1] = None
    ov[1] = None
    ou[3] = None
    ov[3] = None
    print "providing some nones at 1,3"
    dx,dy,u,v = adams_bashforth_dx(dt, x, y, velfun, dict, ou, ov)
    print dx
    print dy
    print "\n---------------------------------"
    print u
    print v
    print "Adams-bash-eul VALIDATED OK ON 2009/11/24"




if __name__ == '__main__':
    test_main()

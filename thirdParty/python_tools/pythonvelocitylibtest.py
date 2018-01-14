import velocitylib
print '\n\n\n\t\t-----------------elementvel now'
"""ccc=[0.0]*4
x=[0.0]
y=[0.0]
ccc=[0.0]
c=[1.0]
cc = [0.0]
s=[1.0]
ex=[0.5]
ey=[0.5]
u,v = velocitylib.vortex_element_vel_python(x,y,ccc,c,cc,s,ex,ey); uo,vo = velocitylib.vortex_element_vel(x,y,ccc,c,cc,s,ex,ey)
#u,v = velocitylib.vortex_element_vel_python(0,0,0,1,0,1,0.5,0); uo,vo = velocitylib.vortex_element_vel(x,y,ccc,c,cc,s,ex,ey)
print u
print uo
print v
print vo
print [u[i]-uo[i] for i in xrange(len(u))]
print [v[i]-vo[i] for i in xrange(len(u))]
"""



x=[0.69867 ,  0.54589 ,  0.75348   ,0.76529]
y=[0.46133 ,  0.56701,   0.82499  , 0.37496]
s=[0.72289,   0.84098 ,  0.72157   ,0.89466]
c=[0.547485,   0.032465,   0.108338 ,  0.987475]
ex=[0.013244  , 0.472345,   0.138415  , 0.792851]
ey=[0.22945  , 0.92952,   0.22198 ,  0.69124]
cc = [0.37582 ,  0.14753,   0.18299   ,0.35446]
ccc=[0.86905 ,  0.68951 ,  0.44606,   0.85764]

x= [i + 0 for i in x]
c= [i + 0 for i in c]
y= [i + 0 for i in y]
cc= [i + 0 for i in cc]

print '\n\n\n\t\t-----------------lambvel now'

u,v = velocitylib.lamb_vortex_vel_python(x,y,s,c,ex,ey); uo,vo = velocitylib.lamb_vortex_vel(x,y,s,c,ex,ey)

print [u[i]-uo[i] for i in xrange(len(u))]
print [v[i]-vo[i] for i in xrange(len(u))]

print '\n\n\n\t\t-----------------sheetvel now'
u,v = velocitylib.vortex_sheet_vel_python(x,y,c,cc,s,ex,ey); uo,vo = velocitylib.vortex_sheet_vel(x,y,c,cc,s,ex,ey)
print [u[i]-uo[i] for i in xrange(len(u))]
print [v[i]-vo[i] for i in xrange(len(u))]

print '\n\n\n\t\t-----------------elementvel now'
u,v = velocitylib.vortex_element_vel_python(x,y,s,c,cc,s,ex,ey); uo,vo = velocitylib.vortex_element_vel(x,y,s,c,cc,s,ex,ey)
print [u[i]-uo[i] for i in xrange(len(u))]
print [v[i]-vo[i] for i in xrange(len(u))]

u,v = velocitylib.vortex_element_vel_python(x,y,ccc,c,cc,s,ex,ey); uo,vo = velocitylib.vortex_element_vel(x,y,ccc,c,cc,s,ex,ey)

print '\n\n\n\t\t-----------------elementvellinerar now'
print [u[i]-uo[i] for i in xrange(len(u))]
print [v[i]-vo[i] for i in xrange(len(u))]

from vpython import *
#Web VPython 3.2

scene.background = vec(.5,0.5,.55)
shield1 = extrusion(path=[vec(-15.2,0,0), vec(-14.8,0,0)],
    color=color.black, opacity = 0.2,
    shape=[ shapes.rectangle(pos=[0,0],
            length=0.4, height=30, width=0.5), shapes.circle(radius=0.23) ]) 
shield2 = extrusion(path=[vec(-5.2,0,0), vec(-4.8,0,0)],
    color=color.black, opacity = 0.2,
    shape=[ shapes.rectangle(pos=[0,0],
            length=0.4, height=30, width=0.5), shapes.circle(radius=0.23) ])

running = False
bbutton = button(text="Run", pos=scene.title_anchor, bind=Run)

#constants  for particles 
q = -1.602e-19            # in C
m = 1.673e-27           # in Kg

#stage 1 plates
L1 = 15.00               # in meters length of rod 
high_voltage1 = 0       # in Volts
low_voltage1 = 15         # in Volts

#stage 2 plates
L2 = 10.00               # in meters length of rod 
high_voltage2 = 1.00        # in Volts
low_voltage2 = 0
B_in_EPlates = vec(0,0,0.7e-5)

#stage 3 coils
mu = 4*pi*1e-7
i = 1e1    # current of loops (CCW)
n = 40     # number  of loops 
R=  2      # radius  of loops
N = 25     # number  of loops 

#charge
if q < 0:
    c1 = vec(0, 0.8, 1)
else:
    c1 = vec(1, 0, 0)
charge = sphere(pos=vec(-20, 0, 0), color=c1, radius = 0.2, make_trail = True, trail_type="points",
                interval=3, retain=50, trail_radius = 0.2)
    
charge.v = vec(0, 0, 0)
charge.a = vec(0, 0, 0)
POI = charge.pos
charge.q = q
charge.m = m

#stage 1 plates
plate1_left_x = -20.5
plate1_right_x = -18.5

Stage1_L = label(pos = vec(-19.5, 10, 0), text = 'Stage 1')

y_min = -L1/2
y_max = L1/2

plate1_left = cylinder(pos=vec(plate1_left_x,y_min,0), axis=vec(0,L1,0), radius = 0.3, color=vec(1,0,0))
plate1_left.V = high_voltage1

plate1_right = cylinder(pos=vec(plate1_right_x,-L1/2,0), axis=vec(0,L1,0), radius = 0.3, color=vec(0,1,1))
plate1_right.V = low_voltage1

separation1 = plate1_right.pos.x - plate1_left.pos.x 
E1_tot = vec(((plate1_left.V - plate1_right.V)/separation1), 0, 0)


#stage 2 plates
x_min = -15
x_max = x_min + L2 
Stage2_L = label(pos = vec(x_min-x_max, 5, 0), text = 'Stage 2')

plate2_top = cylinder(pos=vec(x_min,1.5,0), axis=vec(L2,0,0), radius = 0.3, color=vec(0,1,1))
plate2_top.V = low_voltage2

plate2_bottom = cylinder(pos=vec(x_min,-1.5,0), axis=vec(L2,0,0), radius = 0.3, color=vec(1,0,0))
plate2_bottom.V = high_voltage2
        
separation2 = plate2_top.pos.y - plate2_bottom.pos.y
E2_tot = vec(0, (plate2_bottom.V - plate2_top.V)/separation2, 0)

#coil arrays and covers for the coils
coil1 = []
coil2= []
coil1_ring = ring(pos = vec(0, 0, 1), axis=vec(0,0,1), radius =R, thickness = 0.2, color = vec(1,0,0))
coil2_ring = ring(pos = vec(0, 0, -1), axis=vec(0,0,1), radius =R, thickness = 0.2, color = vec(1,0,0))

Stage3_L = label(pos = vec(-0, 5, 0), text = 'Stage 3')
              
theta_min = radians(0)
theta_max = radians(360)
angle_tot = theta_max-theta_min
dtheta = angle_tot / (N-1)
ds = R * dtheta

#build ring 1
for theta in arange (theta_min + dtheta/2, theta_max, dtheta):                                          # for ring 1 
    ball = sphere(pos = vec(R*cos(theta), R*sin(theta), 1), radius = 0.05 * ds, color = color.red)
    theta_hat = vec(-sin(theta), cos(theta), 0)
    ball.ds = R*theta_hat
    coil1.append(ball)
        
#build ring 2
for theta in arange (theta_min + dtheta/2, theta_max, dtheta):                                          # for ring 2 
    ball = sphere(pos = vec(R*cos(theta), R*sin(theta), -1), radius = 0.05 * ds, color = color.red)
    theta_hat = vec(-sin(theta), cos(theta), 0)
    ball.ds = R*theta_hat
    coil2.append(ball)
    
# getting the B net total
B_tot_net = get_B(coil1, POI) + get_B(coil2, POI)
B_tot_arrow = arrow(pos=POI, axis=vec(0,0,0), color = color.blue, opacity = 0.5)
# the velocity arrow 
charge_arrow = arrow(pos=POI, axis=vec(0,0,0), color = vec(1,1,0),opacity = 0.5 ) 
# the force arrow
charge_f_arrow = arrow(pos=POI, axis=vec(0,0,0), color = vec(0.4,0.8,0),opacity = 0.45 ) 
#the E field arrow
charge_E_arrow = arrow(pos=POI, axis=vec(0,0,0), color = color.red,opacity = 0.45 )

print("Velocity Yellow arrow")
print("Force Green arrow")
print("B-field Blue arrow")
print("E-field Red arrow")

t=0
dt = 0.5e-5
sim_speed = 1e-4 
f1=graph(width=600,height=225,ymin=-1,
                title="<b><i>Deflection</i> vs <i>t</i></b>" ,
                xtitle="<i>t</i> (s)", ytitle="<i>Deflection</i> (units)",
                foreground = color.black, background = color.white)
f1=gcurve(color=vec(0.05,6,0.84))

#scene.waitfor('click')
while True:
    rate(sim_speed/dt)
    POI = charge.pos
    
    #STAGE 1
    if(charge.pos.x >= plate1_left_x & charge.pos.x <= plate1_right_x):
        F = charge.q * E1_tot
        updatePos(F, charge, dt)
        charge_arrow.pos = POI
        charge_arrow.axis = charge.v*1e-4
        charge_f_arrow.pos = POI
        charge_f_arrow.axis = F*1.3e18
        charge_E_arrow.pos = POI
        charge_E_arrow.axis = E1_tot
    #STAGE 2
    elif ((charge.pos.x) >= x_min & (charge.pos.x) <= x_max):
        F = vec(0,0,0)
        F = charge.q * E2_tot
        F += charge.q * charge.v.cross(B_in_EPlates)
        updatePos(F, charge, dt)
        charge_f_arrow.pos = POI
        charge_f_arrow.axis = F*0.3e21
        charge_arrow.pos = POI
        charge_arrow.axis = charge.v*0.5e-4
        charge_E_arrow.pos = POI
        charge_E_arrow.axis = E2_tot*0.9e1
        B_tot_arrow.pos = POI
        B_tot_arrow.axis = B_in_EPlates*0.5e6
        if ((charge.pos.y >= plate2_top.pos.y) | (charge.pos.y <= plate2_bottom.pos.y)) | (((charge.pos.y >= 0.23/2) | (charge.pos.y <= -0.23/2)) & (charge.pos.x >= (shield2.pos.x - shield2.width/2))):    
            charge_arrow.axis = vec(0,0,0)
            B_tot_arrow.axis = vec(0,0,0)
            charge_f_arrow.axis = vec(0,0,0)
            charge_E_arrow.axis = vec(0,0,0)
            break
        f1.plot( t*1e4,((charge.pos.y)) )
    #STAGE 3
    elif(charge.pos.x >= coil1_ring.pos.x-R & charge.pos.x <= coil1_ring.pos.x+R  & charge.pos.y >= coil1_ring.pos.y-R & charge.pos.y <= coil1_ring.pos.y+R):
        F = vec(0,0,0)
        B_tot_net = get_B(coil1, POI) + get_B(coil2, POI) #          particle moving through the coils 
        F = charge.q * charge.v.cross(B_tot_net)
        updatePos(F, charge, dt)
        charge_arrow.pos = POI
        charge_arrow.axis = charge.v*0.5e-4
        charge_f_arrow.pos = POI
        charge_f_arrow.axis = F*1e18
        B_tot_arrow.pos = POI
        B_tot_arrow.axis = B_tot_net*1e4
        
    else:
        charge.pos += charge.v * dt 
    
        charge_arrow.pos = POI
        charge_arrow.axis = charge.v * 0 
    
        charge_f_arrow.pos = POI
        charge_f_arrow.axis = F * 0 
    
        B_tot_arrow.pos = POI
        B_tot_arrow.axis = B_tot_net*0
    
        charge_E_arrow.pos = POI
        charge_E_arrow.axis = E2_tot*0
        
    if (charge.v.x < 0) & (charge.pos.x <= (shield2.pos.x + shield2.width/2)):
        charge_arrow.axis = vec(0,0,0)
        B_tot_arrow.axis = vec(0,0,0)
        charge_f_arrow.axis = vec(0,0,0)
        charge_E_arrow.axis = vec(0,0,0)
        break
    
    scene.camera.follow(charge)
    t += dt
    
#-------------------------------------------------------------------------------
#Euler-Cromer Method to update position of charge
def updatePos(Force, charge, dt):
    charge.a = Force/charge.m 
    charge.v += charge.a * dt
    charge.pos += charge.v * dt
    
#calculate Net B from ring
def get_B (shape, my_POI):
    B_tot = vec(0,0,0)
    for ball in shape:
        r = my_POI - ball.pos
        B_tot += (n*mu*i / (4*pi)) * (ball.ds).cross(r)/mag(r)**3
    return B_tot
    
#def get_E_approx_net (charge, plate1, plate2):
#    return get_E(plate1, charge.pos) + get_E(plate2, charge.pos)
#    
#def get_E(shape, my_POI):
#    E = vec(0,0,0)
#    for my_ball in shape:
#        r = my_POI - my_ball.pos
#        E += k*my_ball.dq*r/mag(r)**3
#    return E
    
    
#calculate Net B from ring 2
#def get_B2 (shape, my_POI):
#    B_tot2 = vec(0,0,0)
#    for ball in shape:
#        r = my_POI - ball.pos
#        B_tot2 += (n*mu*i / (4*pi)) * (ball.ds).cross(r)/mag(r)**3 
#    return B_tot2

def Run(b):
    global running, remember_dt, dt
    running = not running
    if running: 
        b.text = "Pause"
        dt = remember_dt
    else: 
        b.text = "Run"
        remember_dt = dt
        dt=0

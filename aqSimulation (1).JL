using Printf
using Statistics
using LinearAlgebra
using Plots
using Distances
using DifferentialEquations
using Suppressor
using ProgressBars
using DataFrames
using CSV
using DelimitedFiles
#########################################
########    CONSTANTS           #########
#########################################
eps = 8.854 * 10^-12   #C^2/(N*m^2)
k = 1 / (4 * pi * eps)
mass = 9.109 * 10^-31
elecCharge = -1.602 * 10^-19
num_seconds = 10.0
num_steps = 1000.0
timestep = num_seconds / num_steps
timespan = (0.0, timestep)

#########################################
########    ODE Model           #########
#########################################
function ode!(du, u::Array{Float64,1}, p::Array{Float64,1}, t)
        println(p)
        du[1] = p[1] / mass
        du[2] = p[2] / mass
        du[3] = u[1]
        du[4] = u[2]
end
#########################################
########    FUNCTIONS           #########
#########################################
#solves for kinetic energy
function Kinetic(vel::Array{Float64,1})
        norm_vel = norm(vel)
        return (0.5) * mass * norm_vel^2
end
#solves for potential energy
function PotE(dist::Float64)
        return (k * elecCharge^2) / dist
end
#THIS IS CHEAT VERSION, NEED ACTUAL HEAVISIDE
function heaviside(bound,xpos,ypos)
        if (xpos<bound && xpos>-bound && ypos >-bound && ypos<bound)
                return 1
        else
                return 0
        end
end

aq_comb = []
for a = 0:5.0
        for q in range(45, 80, step=5)
                push!(aq_comb, [a/100,q/100])
        end
end
#aq_comb = [[0.1,0.8], [0.0, 0.1], [0.0,0.2], [0.0, 0.3], [0.0, 0.4], [0.0, 0.6], [0.0, 0.7], [0.0, 0.8], [0.0, 0.9], [-0.1, 0.5], [-0.1, 0.6], [-0.1,0.7], [-0.1, 0.8], [-0.1, 0.9], [-0.1, 1.0], [-0.1, 1.1], [-0.2, 0.7], [-0.2, 0.8], [-0.2, 0.9], [-0.2, 1.0], [-0.2, 1.1], [-0.3, 0.9], [-0.3, 1.0], [-0.3, 1.1], [-0.3, 1.2], [-0.4, 1.0], [-0.4, 1.1], [-0.4, 1.2], [-0.5, 1.1], [-0.5, 1.2], [-0.5, 1.3], [-0.6, 1.2]]



numTrials = 1
r =10.0
omega=10.0
num_e = 5

df = DataFrame(a = Any[], q = Any[], initpos = Any[], r= Any[], omega=Any[],percent_e_confined = Any[])
percentConfinedList = []
bound = r-1.0
upper_x = bound
lower_x = -bound
upper_y = bound
lower_y = -bound
iontrap = true

initpos=Any[]
while(length(initpos)<num_e)
        for i=1:num_e
                pos=[rand(lower_x:upper_x),rand(lower_y:upper_y)]
                if !(pos in initpos)
                        push!(initpos,pos)
                end
        end
end
println(initpos)
#initpos=Any[[7.0, -8.0], [7.0, 11.0], [12.0, -8.0], [-1.0, -9.0], [9.0, -14.0], [-3.0, -2.0], [9.0, -4.0], [-9.0, -12.0], [8.0, 3.0], [-6.0, -4.0]]
global yes=true
for comb = 1:length(aq_comb)


        global electrons = Any[]       #holds current pos,vel,force (x and y comp) for each electron
        global velocities = Any[]      #holds all previous positions for each electrons
        global positions = Any[]       #holds all previous velocities for each electrons

        local pe_now = 0.0
        local counter = 1

        for i = 1:num_e
                push!(positions, [initpos[i]])
                push!(velocities, [[0.0, 0.0]])
                push!(electrons, [initpos[i], [0.0, 0.0], [0.0, 0.0],i])       # 0 velocities, 0 force
        end
        global force_electrons=electrons
        global nonforce_electrons=[]
        for x = 1:Int(num_steps)-1
                local ke_now = 0.0
                local pe_now = 0.0
                local counter = 1
                local pos_now = []
                local force_pos_now = []
                local sub_ke = []


                for e = 1:length(force_electrons)
                        if (heaviside(r,force_electrons[e][1][1],force_electrons[e][1][2]) == 0)
                                push!(nonforce_electrons,force_electrons[e])
                                force_electrons[e]=0
                        end
                end

                global force_electrons=filter(x->x≠0,force_electrons)



                for e = 1:num_e
                        #adds all current electron positions to one array
                        push!(pos_now, electrons[e][1])


                end

                for e =1:length(force_electrons)
                        push!(force_pos_now,force_electrons[e][1])
                end



                #finds total force on each electron
                for e = 1:length(force_electrons)

                        fnet = [0.0, 0.0]
                        #solves for force
                        for j = 1:length(force_electrons)

                                #distance between electrons
                                distance = euclidean(force_electrons[e][1], force_pos_now[j])

                                #adds nothing to force when the loop considers an electron with itself
                                if (distance != 0.0)
                                        #finds the position vector between 2 electrons
                                        pos_vectors = force_electrons[e][1] - force_pos_now[j]

                                        #finds the unit vectors between 2 electrons
                                        unit_vectors = pos_vectors/distance

                                        #find Coulomb force, adds to net force on electron
                                        fcoul= ((k*elecCharge^2)/(distance^2))*unit_vectors

                                        fnet= fnet + fcoul

                                end
                        end

                        #implements ion trap
                        if (iontrap)

                                global ax= aq_comb[comb][1]
                                global qx= aq_comb[comb][2]
                                U= (ax*mass*r^2*omega^2)/(8*elecCharge)
                                V= (qx*mass*r^2*omega^2)/(-4*elecCharge)
                                fx= -2*elecCharge*(U+V*cos(omega*x))*force_electrons[e][1][1]/omega^2
                                fy= -2*elecCharge*(U+V*cos(omega*x))*force_electrons[e][1][2]/omega^2

                                heavi=heaviside(r,force_electrons[e][1][1],force_electrons[e][1][2])
                                fnet=(fnet+[fx,fy])*heavi


                        end
                        #updates the current force on the electron
                        force_electrons[e][3]=fnet

                        #initial values for the ODE
                        u0 = [force_electrons[e][2][1],force_electrons[e][2][2],force_electrons[e][1][1],force_electrons[e][1][2]]

                        #parameters for the ODE (force)
                        param=force_electrons[e][3]

                        @suppress begin

                                #define and solve ODE
                                prob = ODEProblem(ode!,u0,timespan,param)
                                sol = solve(prob,saveat=timestep,save_everystep=false)


                        #update pos and vel for electron using ODE solution. Adds new pos and vel to respective arrays
                                push!(positions[force_electrons[e][4]],[sol[2][3],sol[2][4]])
                                push!(velocities[force_electrons[e][4]],[sqrt(sol[2][1]^2+sol[2][2]^2)])
                                force_electrons[e][2]=[sol[2][1],sol[2][2]]
                                force_electrons[e][1]= [sol[2][3],sol[2][4]]

                        end
                        



                end
                for e=1:length(nonforce_electrons)
                        push!(positions[nonforce_electrons[e][4]],[10^3,10^3])
                        push!(velocities[nonforce_electrons[e][4]],nonforce_electrons[e][2])
                end


        end

        electronsConfined = 0
        for i = 1:num_e
                if (abs(electrons[i][1][1]) < r) && (abs(electrons[i][1][2]) < r)
                        electronsConfined += 1
                end
        end
        global percentConfined = (electronsConfined / num_e)*100

        push!(df, (aq_comb[comb][1], aq_comb[comb][2], initpos,r,omega, percentConfined))
        te=string(comb)*"/"*string(length(aq_comb))
        println("Test: ",te,". Combination of ", [ax,qx], " has average percentConfined ", percentConfined)
end

show(df, summary = false)
CSV.write("/Users/drewj/Documents/atom/uropjulia/MassConfinementSimTest.csv", df)

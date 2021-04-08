using Printf
using Statistics
using LinearAlgebra
using Plots
using Distances
using DifferentialEquations
pyplot()

#########################################
########    CONSTANTS           #########
#########################################

eps = 8.854*10^-12   #C^2/(N*m^2)
k = 1/(4*pi*eps)
mass = 9.109*10^-31
elecCharge = -1.602*10^-19

#########################################
########    SETUP               #########
#########################################

#timepoints

#arrays
electrons=Any[]       #holds current pos,vel,force (x and y comp) for each electron
velocities=Any[]      #holds all previous positions for each electrons
positions=Any[]       #holds all previous velocities for each electrons
ke=Any[]              #holds all kinetic energy values for each timestep
pe=Any[]              #holds all potential energy values for each timestep
totale=Any[]          #holds all total energy values for each timestep


#number of electrons
num_e=2

#timepoints
num_seconds=10
num_steps=20
timestep=10.0/20.0

#initial conditions
initpos=Any[[1.0,1.0],[-1.0,-1.0]]

#add initial conditions to arrays
for i=1:num_e
    push!(positions,initpos[i])
    push!(velocities,[0,0])
    push!(electrons,[initpos[i],[0,0],[0,0]])       # 0 velocities, 0 force
end

println(electrons[1][3])
#########################################
########    ODE Model           #########
#########################################


function ode!(du,u,t)
    du[1] = u[1]/mass
    du[2] = u[2]
end

#########################################
########    FUNCTIONS           #########
#########################################

function changePos(pos,index)
        electrons[index][1]=pos
        push!(positions[index],pos)
end

function changePos(vel,index)
        electrons[index][2]=vel
        push!(velocities[index],vel)
end

function changeFor(force,index)
        electrons[index][3]=force
end

function Kinetic(vel::Array{Float64,1})
        norm_vel = norm(vel)
        return (1/2)*mass*norm_vel^2
end

function PotE(dist::Float64)
        return (k*elecCharge^2)/dist
end

#########################################
########    MAIN LOOP           #########
#########################################

for x = 1:num_e
        global ke_now = 0.0
        global pe_now = 0.0
        global counter = 0
        pos_now = []

        for e = 1:num_e
                push!(pos_now, electrons[e][1])
                println("KE from electron ", e, " is: ", Kinetic(electrons[e][2]))
                global ke_now = ke_now + Kinetic(electrons[e][2])
                println("total KE is ", ke_now)
        end


        for e = 1:num_e
                for i = counter+1:num_e
                        pedist = euclidean(electrons[e][1], pos_now[i])
                        if pedist != 0
                                global pe_now += PotE(abs(pedist))
                        end
                end
                global counter = counter + 1
                #global pe_now = pe_now + PotE(abs(getDistance(getPos(electrons[e]), getPos(staticelectrons[1])))) + PotE(abs(getDistance(getPos(electrons[e]), getPos(staticelectrons[2]))))
        end


        #work with vectors for this part

        for e = 1:num_e
                global fnet = [0.0, 0.0]
                for j in 1:num_e
                        distance = euclidean(electrons[e][1], pos_now[j])
                        if (distance == 0)
                                fnet = fnet + [0.0,0.0]

                        else
                        #println("Distance between electron ", e, " and ",j, " is ", distance )
                                pos_vectors = electrons[e][1] - pos_now[j] #position vector
                        #println("position vector is ", pos_vectors)
                                unit_vectors = pos_vectors/distance
                        #println("unit vector is ", unit_vectors)
                                fnet = fnet + ((k*elecCharge^2)/(distance^2))*unit_vectors
                        #println("fnet is", fnet)

                        #println("loop ", x, " dist b/e electron ", e, " and ", i, " is ", distance )
                        end
                end

                println("net force in loop ", x, " and electron ", e, " is ", fnet)
                electrons[e][3]=fnet
                timespan=(0.0,timestep)

                #for x component
                local u0x = [electrons[e][3][1],electrons[e][2][1]]
                local probx = ODEProblem(ode!,u0x,timespan)
                local solx = solve(probx,save_everystep=false)
                electrons[e][1][1]=sol[1]
                electrons[e][2][1]=sol[2]


                #for y component
                local u0y = [electrons[e][3][2],electrons[e][2][2]]
                local proby = ODEProblem(ode!,u0y,timespan)
                local soly = solve(proby,save_everystep=false)
                electrons[e][1][2]=sol[1]
                electrons[e][2][2]=sol[2]

                push!(positions[e],electrons[e][1])
                push!(velocities[e],electrons[e[2]])
        end


end


for e = 1:num_e

        println("the final net force on electron ", e, " is: ", electrons[e][3])
end

println("the final ke of the system is: ", ke_now)
println("the final pe of the system is: ", pe_now)

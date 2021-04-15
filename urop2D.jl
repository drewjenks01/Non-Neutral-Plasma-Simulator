using Printf
using Statistics
using LinearAlgebra
using Plots
using Distances
using DifferentialEquations
using Suppressor
using ProgressBars


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

#arrays
electrons=Any[]       #holds current pos,vel,force (x and y comp) for each electron
velocities=Any[]      #holds all previous positions for each electrons
positions=Any[]       #holds all previous velocities for each electrons
ke=Any[]              #holds all kinetic energy values for each timestep
pe=Any[]              #holds all potential energy values for each timestep
totale=Any[]          #holds all total energy values for each timestep


#number of electrons
num_e=10

#magnetic field  (if too big then discontinuity error)
B=[0,0,10^-10]

#timepoints
num_seconds=10.0
num_steps=1000.0
timestep=num_seconds/num_steps
timespan=(0.0,timestep)

#initial position conditions
upper_x = 30.0
lower_x = -30.0
upper_y = 30.0
lower_y = -30.0

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

#add initial conditions to arrays
for i=1:num_e
    push!(positions,[initpos[i]])
    push!(velocities,[[0.0,0.0]])
    push!(electrons,[initpos[i],[0.0,0.0],[0.0,0.0]])       # 0 velocities, 0 force
end

#settings
animation=false
confined=true


#########################################
########    ODE Model           #########
#########################################


function ode!(du,u::Array{Float64,1},p::Array{Float64,1},t)
    println(p)
    du[1] = p[1]/mass
    du[2] = p[2]/mass
    du[3] = u[1]
    du[4] = u[2]

end

#########################################
########    FUNCTIONS           #########
#########################################

#solves for kinetic energy
function Kinetic(vel::Array{Float64,1})
        norm_vel = norm(vel)
        return (0.5)*mass*norm_vel^2
end

#solves for potential energy
function PotE(dist::Float64)
        return (k*elecCharge^2)/dist
end



#########################################
########    MAIN LOOP           #########
#########################################

# f=open("positions.txt", "w")
println("-------------------------------------------------------------------------------")
for x = tqdm(1:Int(num_steps)-1)
        local ke_now = 0.0
        local pe_now = 0.0
        local counter = 1
        local pos_now = []


        for e = 1:num_e
                #adds all current electron positions to one array
                push!(pos_now, electrons[e][1])

                #adds KE from each electrons together
                ke_now = ke_now + Kinetic(electrons[e][2])

        end


        for e = 1:num_e

                #solves for PE
                for i = counter+1:num_e
                        pedist = euclidean(electrons[e][1], pos_now[i])
                        pe_now += PotE(abs(pedist))
                end
                counter += 1
        end


        #finds total force on each electron
        for e = 1:num_e
                fnet = [0.0, 0.0]
                #solves for force
                for j = 1:num_e

                        #distance between electrons
                        distance = euclidean(electrons[e][1], pos_now[j])

                        #adds nothing to force when the loop considers an electron with itself
                        if (distance == 0.0)
                                fnet = fnet + [0.0,0.0]

                        else
                                #finds the position vector between 2 electrons
                                pos_vectors = electrons[e][1] - pos_now[j]

                                #finds the unit vectors between 2 electrons
                                unit_vectors = pos_vectors/distance

                                #find Coulomb force, adds to net force on electron
                                fcoul= ((k*elecCharge^2)/(distance^2))*unit_vectors

                                fnet= fnet + fcoul

                        end
                end

                #finds mag force, adds to net force
                if (confined)
                        vel=[electrons[e][2][1],electrons[e][2][2],0]
                        fmag=elecCharge*cross(vel,B)
                        fnet[1]=fnet[1]+fmag[1]
                        fnet[2]=fnet[2]+fmag[2]
                end
                #updates the current force on the electron
                electrons[e][3]=fnet

                #initial values for the ODE
                u0 = [electrons[e][2][1],electrons[e][2][2],electrons[e][1][1],electrons[e][1][2]]

                #parameters for the ODE (force)
                param=electrons[e][3]

                @suppress begin

                        #define and solve ODE
                        prob = ODEProblem(ode!,u0,timespan,param)
                        sol = solve(prob,saveat=timestep,save_everystep=false)


                #update pos and vel for electron using ODE solution. Adds new pos and vel to respective arrays
                        push!(positions[e],[sol[2][3],sol[2][4]])
                        push!(velocities[e],[sol[2][1],sol[2][2]])
                        electrons[e][2]=[sol[2][1],sol[2][2]]
                        electrons[e][1]= [sol[2][3],sol[2][4]]

                end
                if(x==1)
                        global beginningpe=pe_now
                end

                # for e=1:num_e
                #         println(f,positions[e][x+1])
                # end

        end

        #adds energies to respective arrays
        push!(ke,ke_now)
        push!(pe,pe_now)
        push!(totale,ke_now+pe_now)

end

finale=ke[Int(num_steps)-1]+pe[Int(num_steps)-1]

println("the original pe of the system was: ", beginningpe)
println("the final tote of the system is: ", finale)
println("The % error is: ", 100*((finale-beginningpe)/beginningpe))
println("the final ke of the system is: ", ke[Int(num_steps)-1])
println("the final pe of the system is: ", pe[Int(num_steps)-1])
println("the final pos of e1 is: ", electrons[1][1])
println("the final pos of e2 is: ", electrons[2][1])
println("the final vel of e1 is: ", electrons[1][2])
println("the final vel of e2 is: ", electrons[2][2])


#########################################
########    PLOTTING            #########
#########################################

#save data to .txt files
f=open("positions.txt", "w")
for e=1:num_e
        for i=1:(Int(num_steps))
                println(f, positions[e][i])
        end
end
close(f)

v=open("velocities.txt", "w")
for e=1:num_e
        for i=1:Int(num_steps)
                println(v, (velocities[e][i]))
        end
end
close(v)

#animation
global plt=plot(xlims=(-100,100),ylims=(-100,100),legend=false, widen=true)

if (animation==true)
        xdata=Any[]
        ydata=Any[]

        for e=1:num_e
                global xlist=Any[]
                global ylist=Any[]
                for i=1:length(positions[1])
                        push!(xlist, positions[e][i][1])
                        push!(ylist, positions[e][i][2])
                end
                push!(xdata,xlist)
                push!(ydata,ylist)
        end

        global plt=plot(xlims=(-100,100),ylims=(-100,100),legend=false,widen=true)
        global anim = @animate for i = tqdm(1:Int(num_steps))
                for e=1:num_e
                    plot!(plt,xdata[e][i:i+1], ydata[e][i:i+1],linewidth=2,linealpha=1)
                end
        end every 500

        gif(anim, "anim.gif", fps = 5)
end

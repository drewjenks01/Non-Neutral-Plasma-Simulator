

#need to import these if first time using them
using Printf
using Statistics
using LinearAlgebra
using Plots
using Distances
using DifferentialEquations
using Suppressor
using ProgressBars
using CSV
using DataFrames
using PlotlyBase
#import Distributions: Uniform
using Distributions
#########################################
########    CONSTANTS           #########
#########################################
mz = 10.0       #mass:charge ratio
eps = 8.854 * 10^-12   #vacuum permitivity
k = 1 / (4 * pi * eps)
elecCharge = -1.602*10^-19      #charge of electron (Coulombs)
mass =  9.109*10^-31 * mz       #mass of electron (kg)


#########################################
########    SETUP               #########
#########################################
#arrays
electrons = Any[]       #holds current pos,vel,force (x and y comp) for each electron
velocities = Any[]      #holds all previous positions for each electrons
positions = Any[]       #holds all previous velocities for each electrons
ke = Float64[]              #holds all kinetic energy values for each timestep
pe = Float64[]              #holds all potential energy values for each timestep
totale = Any[]          #holds all total energy values for each timestep

#number of electrons
global num_e = 100

#timepoints
num_seconds = 10.0
num_steps = 4000.0
timestep = num_seconds / num_steps
timespan = (0.0, timestep)

#initial position conditions
r=10.0       #1/2 length of each side of trap
bound=r-0.01
upper_x = bound         #these are bounds for creating random starting positions for electrons (making sure they are within box)
lower_x = -bound
upper_y = bound
lower_y = -bound


#creates random initial positions ensuring within trap and no duplicates
while(length(initpos)<num_e)

        pos=[round(rand(Uniform(lower_x,upper_x)),digits=3),round(rand(Uniform(lower_y,upper_y)),digits=3)]
        if norm(pos)<bound && !(pos in initpos)
                push!(initpos,pos)
        end

end

#displaces one electron to test stability
# displace_percent=0.01
# start=initpos[50][1]
# start=start+start*displace_percent
# initpos[50][1]=start


#add initial conditions to arrays
for i=1:num_e
    push!(positions,[initpos[i]])
    push!(velocities,[[0.0]])
    push!(electrons,[initpos[i],[0.0,0.0],[0.0,0.0]])       # 0 velocities, 0 force to start
end




#settings
iontrap = true
dissipation = false


coulomb=Any[]           #used for printing avg coulomb at end
#damp=Any[]

#########################################
########    ODE Model           #########
#########################################
#creates the model in order for the ode to be solved
function ode!(du, u::Array{Float64,1}, p::Array{Float64,1}, t)
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
#checks if electron is still in trap or not, return 1 if true, 0 if not
function heaviside(bound, xpos, ypos)
        if (xpos < bound && xpos > -bound && ypos > -bound && ypos < bound)
                return 1
        else
                return 0
        end
end



electronCount = []              #used for plotting electron count over time
push!(electronCount, num_e)


forceCount = [[] for i = 1:num_e]       #keeps track of forces on electrons left in trap

vels=Any[]      #used for plotting acceleration

#distList = Any[]
displace=Any[]
overdamped=false
#########################################
########    MAIN LOOP           #########
#########################################
# f=open("positions.txt", "w")
println("-------------------------------------------------------------------------------")
anim = @animate for x in tqdm(1:Int(num_steps)-1)
        #= creates an animation of the electrons in the trap. saves every 15th frame (determined in the 'end' statemen of loop)
           tqdm is used for the progress bar showing in the REPL   =#

        #used for plotting accel
        #beforevel=sqrt(electrons[1][1][1]^2+electrons[1][1][2]^2)

        #keeping track of x and y pos for each electron in current step. used for plotting
        xval = []
        yval = []
        for e = 1:length(electrons)
                push!(xval, electrons[e][1][1])
                push!(yval, electrons[e][1][2])
        end

        #plots current positions of electrons
        Plots.scatter(xval, yval, aspect_ratio = 1,label = "Ion")
        plot!([(-r,-r),(-r,r),(r,r),(r,-r)],  seriestype = [:shape], c=:blue, fillalpha = 0.2,linecolor=:blue,linewidth=3,label = "Trap")

        #bounds for plot axes
        ylims!(-(r+5),(r+5))
        xlims!(-(r+5),(r+5))

        #variables used for finding current ke and pe of system
        local ke_now = 0.0
        local pe_now = 0.0
        local counter = 1
        local pos_now = []

        for e = 1:length(electrons)
                #adds all current electron positions to one array
                push!(pos_now, electrons[e][1])
                #adds KE from each electrons together
                ke_now = ke_now + Kinetic(electrons[e][2])
        end

        for e = 1:length(electrons)
                #solves for PE
                for i = counter+1:length(electrons)
                        pedist = euclidean(electrons[e][1], pos_now[i])
                        pe_now += PotE(abs(pedist))
                end
                counter += 1
        end

        #used for finding avg coulombic force on electrons during first step...used for dissipation constant throughout the sim
        if(x==1)
                fnet = [0.0, 0.0]
                for j = 1:length(electrons)
                        #distance between electrons
                        distance = euclidean(electrons[1][1], pos_now[j])
                        #adds nothing to force when the loop considers an electron with itself
                        if (distance != 0.0)
                                #finds the position vector between 2 electrons
                                pos_vectors = electrons[1][1] - pos_now[j]
                                #finds the unit vectors between 2 electrons
                                unit_vectors = pos_vectors / distance
                                #find Coulomb force, adds to net force on electron
                                fcoul =
                                        ((k * elecCharge^2) / (distance^2)) * unit_vectors
                                fnet = fnet + fcoul

                        end
                end
                global avgcoul=sqrt(fnet[1]^2+fnet[2]^2)
        end

        #finds total force on each electron
        for e = 1:length(electrons)
                fnet = [0.0, 0.0]
                #solves for force
                for j = 1:length(electrons)
                        #distance between electrons
                        distance = euclidean(electrons[e][1], pos_now[j])
                        #adds nothing to force when the loop considers an electron with itself
                        if (distance != 0.0)
                                #finds the position vector between 2 electrons
                                pos_vectors = electrons[e][1] - pos_now[j]
                                #finds the unit vectors between 2 electrons
                                unit_vectors = pos_vectors / distance
                                #find Coulomb force, adds to net force on electron
                                fcoul =
                                        ((k * elecCharge^2) / (distance^2)) * unit_vectors
                                fnet = fnet + fcoul
                        end
                end
                #push!(coulomb,sqrt(fnet[1]^2+fnet[2]^2))

                #implements ion trap
                if (iontrap)
                        #trap variables
                        omega = 10                      #self chosen
                        ax = 0.7                      #self chosen
                        qx = 0.7                     #self chosen

                        #from papers
                        U = (ax * mass * r^2 * omega^2) / (8 * elecCharge)
                        V = (qx * mass * r^2 * omega^2) / (-4 * elecCharge)

                        #force components from trap
                        fx = ( -2 * elecCharge * (U + V * cos(omega * x)) * electrons[e][1][1] / r^2 )
                        fy = ( -2 * elecCharge * (U + V * cos(omega * x)) * electrons[e][1][2] / r^2 )

                        #potential from trap
                        pot_0 = 2 * (U + V * cos(omega * x))
                        potential = pot_0 / (2 * r^2) * ((electrons[e][1][1])^2 - (electrons[e][1][2])^2 )
                        pe_trap = abs(elecCharge * potential)
                        pe_now += pe_trap

                        #do we still need this?
                        heavi = heaviside(
                                r,
                                electrons[e][1][1],
                                electrons[e][1][2],
                        )
                        fnet = (fnet + [fx, fy]) * heavi



                end

                constan=18.5
                lowb=0.0
                conspan=[constan:-((constan-lowb)/num_steps):lowb;]
                #sets disspiation
                if(dissipation)

                        #diss constant. self chosen
                        cons=15*avgcoul/(num_e)

                        #check if system overdamped
                        if(x==1)
                                if(cons/mass>omega)
                                        overdamped=true
                                end
                        end

                        #push!(damp,sqrt((electrons[e][2][1]*cons)^2+(electrons[e][2][2]*cons)^2))
                        fnet-= [electrons[e][2][1]*cons,electrons[e][2][2]*cons]

                end



                #updates the current force on the electron
                electrons[e][3] = fnet
                push!(forceCount[e], norm(fnet))
                #initial values for the ODE
                u0 = [
                        electrons[e][2][1],
                        electrons[e][2][2],
                        electrons[e][1][1],
                        electrons[e][1][2],
                ]
                #parameters for the ODE (force)
                param = electrons[e][3]

                #@suppress is so the solver doesnt print something out every iteration
                @suppress begin
                        #define and solve ODE
                        prob = ODEProblem(ode!, u0, timespan, param)
                        sol = solve(
                                prob,
                                saveat = timestep,
                                save_everystep = false,
                        )
                        #update pos and vel for electron using ODE solution. Adds new pos and vel to respective arrays
                        push!(positions[e], [sol[2][3], sol[2][4]])
                        push!(velocities[e], [sqrt(sol[2][1]^2 + sol[2][2]^2)])

                        electrons[e][2] = [sol[2][1], sol[2][2]]
                        electrons[e][1] = [sol[2][3], sol[2][4]]
                end
                if (x == 1)
                        global beginningpe = pe_now
                end
        end

        #finding index of electrons that are out of trap
        delete_index = []
        for e = 1:length(electrons)
                #if abs(electrons[e][1][1]) > r || abs(electrons[e][1][2]) > r
                if abs(electrons[e][1][1])>r || abs(electrons[e][1][2])>r
                        push!(delete_index, e)
                end
        end

        #deleting electrons that are out of trap
        deleteat!(electrons, delete_index)
        deleteat!(positions, delete_index)
        deleteat!(velocities, delete_index)
        deleteat!(forceCount, delete_index)

        #adds energies to respective arrays
        push!(ke, ke_now)
        push!(pe, pe_now)
        push!(totale, ke_now + pe_now)
        push!(electronCount, length(electrons))
        disp=Any[]
        for e =1:length(electrons)
                d=euclidean(electrons[e][1],initpos[e])
                if electrons[e][1][1]-initpos[e][1]<0
                        if electrons[e][1][2]-initpos[e][2]<0
                                d=-1*d
                        end
                end

                push!(disp, d)
        end
        push!(displace,sum(disp)/length(disp))





#=
        if length(positions[1])>1
                using Distances
                i = length(positions[1])
                push!(distList, euclidean(positions[1][i],positions[1][i-1]))
        end
=#

#used for accel
# aftervel=sqrt(electrons[1][1][1]^2+electrons[1][1][2]^2)
# push!(vels,[beforevel,aftervel])

end every 1000
finalpos=Any[]
for e=1:length(electrons)
        push!(finalpos, electrons[e][1])
end

#saving the animation
gif(anim, "/Users/drewj/Documents/atom/uropjulia/Electrons.gif")

#prints final # electrons left
println(length(electrons))
println("overdamped? : ",overdamped)
#println("final pos: ",finalpos)


#printing final values to make sure everything working correctly

# println("the original pe of the system was: ", beginningpe)
println("the final tote of the system is: ", ke[Int(num_steps)-1] + pe[Int(num_steps)-1])
# println("The % error is: ", 100 * ((finale - beginningpe) / beginningpe))
# println("the final ke of the system is: ", ke[Int(num_steps)-1])
#rintln("the final pe of the system is: ", pe[Int(num_steps)-1])
# println("avg coulomb: ",sum(coulomb)/length(coulomb))
#println("avg damp: ",sum(damp)/length(damp))
#println("ratio: ",(sum(coulomb)/length(coulomb))/(sum(damp)/length(damp)))
#########################################
########    PLOTTING            #########
#########################################

#energy
# time = (0:timestep:num_seconds-2*num_seconds/num_steps)
# plotly()
# plot(time, ke, label = "KE", ticks = :native)
# plot!(time, pe, label = "PE", ticks = :native)
# #energyPlot = plot!(time, ke+pe, label = "TOTAL E", title="$num_e electrons; init position: $initpos; $num_steps steps", showaxis = true)
# energyPlot = plot!(
#         time,
#         ke + pe,
#         label = "TOTAL E",
#         title = "energy for $num_e electrons; $num_steps steps",
#         showaxis = true,
#         legend = :outertopleft,
#         ticks = :native,
# )
#
#
# #accel
# accels=Any[]
# for v = 1:length(vels)
#         push!(accels,(vels[v][2]-vels[v][1])/timestep)
# end
# time = (0:timestep:num_seconds-2*num_seconds/num_steps)
# accelplot=plot(time,accels,title="inst accel")
#
# #electron count over time and force on electrons
# if length(electrons)!= 0
#         #plot electron count over time
#         time = (0:num_seconds/num_steps:num_seconds-num_seconds/num_steps)
#         electronCount = (1/num_e) .* electronCount
#         electronPlot = plot(time, electronCount, title = "electron count over time")#, yaxis=:log, xaxis = ((0,num_seconds), [10^-3, 10^-2, 10^-1, 10^0, 10^1], :log))#, xaxis=((0, 10), :log))
#         #yaxis!("electron count", :log10)
#         #xaxis!("time", :log10)
#         #plot net force on each remaining electron over time
#         forceTime = (0:num_seconds/num_steps:num_seconds-2*num_seconds/num_steps)
#         remElectrons = length(electrons)
#         forcePlot = plot(forceTime, forceCount, title = "force on electrons $remElectrons / $num_e over time", ticks =:native)
#
#
#         #plot all at once
#         plot(energyPlot, electronPlot, forcePlot,accelplot, layout = (4,1), legend = false)
# end


#inter-electron distance
# dists=Any[]
# shortarr=Any[]
# avgdist=Any[]
#
#
# for e = 1:length(electrons)
#
#         #solves for force
#         for j = 1:length(electrons)
#                 #distance between electrons
#                 distance = euclidean(electrons[e][1], electrons[j][1])
#                 #adds nothing to force when the loop considers an electron with itself
#                 if (distance != 0.0)
#                         push!(dists,distance)
#                 end
#         end
#         shortarr = sort(dists)
#         fifth = int(num_e/5)
#         avgdist = shortarr[1:fifth]
#         avg = sum(avgdist)/length(avgdist)
#         push!(avgdist, avg)
#
#
#
# end
#
# print("avg dist: ", sum(avgdist)/length(avgdist))

time = (0:timestep:num_seconds-2*num_seconds/num_steps)
dipsl=plot(time,displace,title="Avg displacement from initial positions",legend=false,ylims=(-2,2))
yaxis!("Displacement")
xaxis!("Time")

export anakgplot, sparsematrix_plot

using GLMakie
using Random
using Distributions
using Colors
using LinearAlgebra
using SparseArrays

using AssociativeScenesSimulation

mutable struct State
    neuronpositions::Dict{UInt, Point3f}
    connections::Dict{UInt, Vector{Point3f}}
    activeconnections::Dict{UInt, Vector{Point3f}}
    scenespositions::Dict{UInt, Point3f}
    scenesconnections::Dict{UInt, Vector{Point3f}}
    activeneurons::Vector{UInt}

    State() = new(Dict(), Dict(), Dict(), Dict(), Dict(), Vector())
end

state = State()

function anakgplot(
    file::String=ANAKG_SAMPLE_MAT_FILE;
    resolution=primary_resolution(), camera3d=true, background=:white,
    randomscale=0.3, sides=6, sparsity=1.0,
    neuron_outercolor=RGBA{Float32}(0.39, 0.58, 0.93, 0.44), 
    neuron_innercolor=RGBA{Float32}(0.0, 0.0, 0.5, 0.88),  
    scene_outercolor=RGBA{Float32}(0.7, 0.7, 0.7, 0.44), 
    scene_innercolor=RGBA{Float32}(0.08, 0.08, 0.08, 0.88), 
    active_outercolor=RGBA{Float32}(0.6, 0.98, 0.6, 0.58), 
    active_innercolor=RGBA{Float32}(0.0, 0.39, 0.0, 0.95), 
    neuron_outersize=0.11 , neuron_innersize=0.035,
    scene_outersize=0.18 , scene_innersize=0.058,
    connectionthickness=0.025, connectioncolor=RGBA{Float32}(0.15, 0.25, 0.55, 0.35),
    active_connectioncolor=RGBA{Float32}(0.15, 0.25, 0.55, 0.85),
    scene_connectionthickness=1.5, scene_connectioncolor=RGBA{Float32}(0.13, 0.55, 0.13, 0.85),
    showallconnections=false
)
    global state = State()

    figure, parentscene, camera = createscenes(
        resolution, camera3d, background
    )
    
    data = anakgload(file)
    data = data[1:Int(floor(size(data, 1) * sparsity)), 1:Int(floor(size(data, 2) * sparsity))]
    neuronsnumber = size(data, 1)
    neuronsnumber_perside = neuronsnumber / sides
    sidelength = √(neuronsnumber_perside)
    sidelength_int = Int(floor(sidelength))
    lengthdiff = neuronsnumber - ceil(sidelength_int^2 * sides)
    
    lengthdiff_int = lengthdiff ÷ (sidelength_int * sides)
    lengthdiff_rem = lengthdiff % sidelength_int
    sideranges = LinRange[]
    for i in 1:2sides
        siderange = if i % 2 == 0
            siderange = LinRange(
                0.0, 
                sidelength_int, 
                sidelength_int + lengthdiff_int
            )
            siderange
        else
            if lengthdiff_rem > 0 && i == 2sides - 1
                LinRange(
                    0.0, 
                    sidelength_int, 
                    sidelength_int + 1
                )
            else
                LinRange(0.0, sidelength_int, sidelength_int)
            end
        end
        push!(sideranges, siderange)
    end

    updatepositions(
        state, sides, neuronsnumber, sideranges, sidelength_int, randomscale
    )
    drawscences(
        state, parentscene, 208,
        sidelength_int, 2.5,
        scene_outercolor, scene_innercolor, scene_outersize, scene_innersize,
        scene_connectionthickness, scene_connectioncolor,
        active_outercolor, active_innercolor
    )
    drawneurons(
        state, data, parentscene,
        neuron_outercolor, neuron_innercolor, neuron_outersize, neuron_innersize,
        active_outercolor, active_innercolor
    )
    drawconnections(
        state, data, parentscene, connectionthickness, 
        connectioncolor, active_connectioncolor, showallconnections
    )

    figure
end

function drawscences(
    state::State, parentscene::LScene, scenesnumber::Int,
    sidelength_int::Int, sideoffset::Float64,
    scene_outercolor::RGBA, scene_innercolor::RGBA,
    scene_outersize::Float64, scene_innersize::Float64,
    linewidth::Float64, color,
    active_outercolor, active_innercolor
)
    side_smalloffset = 0.035sidelength_int
    lastside = 2side_smalloffset + sidelength_int
    scenceslength_int = Int(ceil(√(scenesnumber)))
    siderangey = LinRange(0.7, 0.7 + sidelength_int, scenceslength_int)
    siderangez = LinRange(0.0, sidelength_int, scenceslength_int + 1)
    sceneid = 1
    scenespositions = state.scenespositions
    for y=siderangey, z=siderangez
        sceneid > scenesnumber && break

        ypos = y + side_smalloffset
        zpos = z + side_smalloffset
        scenespositions[sceneid] = Point3f(lastside + sideoffset, ypos, zpos)
        sceneid += 1
    end
    
    # activescenes = [scenesnumber ÷ 2, scenesnumber - 3]
    activescenes = [scenesnumber - 3]
    activeindices = findall(x -> x in activescenes, collect(keys(scenespositions)))

    positions = collect(values(scenespositions))
    innercolors = repeat([scene_innercolor], length(scenespositions))
    innercolors[activeindices] .= active_innercolor
    innersizes = repeat([scene_innersize], length(scenespositions))
    meshscatter!(
        parentscene, positions, 
        markersize=innersizes, color=innercolors
    )
    
    outercolors = repeat([scene_outercolor], length(scenespositions))
    outercolors[activeindices] .= active_outercolor
    outersizes = repeat([scene_outersize], length(scenespositions))
    meshscatter!(
        parentscene, positions, 
        markersize=outersizes, color=outercolors
    )

    for sceneid in activescenes
        sceneconnections = state.scenesconnections
        sceneposition = scenespositions[sceneid]
        sceneconnections[sceneid] = Point3f[sceneposition]
        neuronids = Int.(floor.(1.0 .+ (rand(25) * 100)))
        append!(state.activeneurons, neuronids)
        neuronpositions = collect(values(state.neuronpositions))[neuronids]
        for pos in neuronpositions
            append!(sceneconnections[sceneid], Point3f[pos, sceneposition])
        end

        for (id, points) in sceneconnections
            linewidth = linewidth
            lines!(parentscene, points, linewidth=linewidth, color=color)
        end
    end
end

function drawconnections(
    state::State, data::Matrix{Float64}, parentscene::LScene,
    linewidth::Float64, color, active_connectioncolor, showall
)
    connections = state.connections
    activeconnections = state.activeconnections
    neuronpositions = state.neuronpositions
    neuronpositions = collect(values(neuronpositions))
    neuroncounter = data[collect(diagind(data))][1:length(neuronpositions)]
    neuroncounter_scaled = minmax(neuroncounter)

    for x = axes(data)[1][1:length(neuronpositions)]
        for y in 1:(x - 1)
            if data[x, y] > 0.0 && x != y
                xpos = neuronpositions[x]
                ypos = neuronpositions[y]
                if !(x in state.activeneurons && y in state.activeneurons)
                    if haskey(connections, x)
                        append!(connections[x], Point3f[xpos, ypos])
                    else
                        connections[x] = Point3f[xpos, ypos]
                    end
                end
            end
        end
    end

    for x = state.activeneurons, y = state.activeneurons
        if x != y
            xpos = neuronpositions[x]
            ypos = neuronpositions[y]
            if haskey(activeconnections, x)
                append!(activeconnections[x], Point3f[xpos, ypos])
            else
                activeconnections[x] = Point3f[xpos, ypos]
            end
        end
    end

    for key in keys(connections)
        push!(connections[key], neuronpositions[key])
    end

    if showall
        for (id, points) in connections
            linewidth = linewidth + 0.0058 * neuroncounter_scaled[id]
            lines!(parentscene, points, linewidth=linewidth, color=color)
        end
    end
    for (id, points) in activeconnections
        linewidth = linewidth + 1.0 * neuroncounter_scaled[id]
        lines!(parentscene, points, linewidth=linewidth, color=active_connectioncolor)
    end
end

function drawneurons(
    state::State, data::Matrix{Float64}, parentscene::LScene,
    neuron_outercolor::RGBA, neuron_innercolor::RGBA,
    neuron_outersize::Float64, neuron_innersize::Float64,
    active_outercolor, active_innercolor
)
    neuronpositions = collect(values(state.neuronpositions))
    neuroncounter = data[collect(diagind(data))][1:length(neuronpositions)]
    neuroncounter_scaled = minmax(neuroncounter)

    activeindices = findall(
        x -> x in state.activeneurons, collect(keys(neuronpositions))
    )
    
    innercolors = repeat([neuron_innercolor], length(neuronpositions))
    innercolors[activeindices] .= active_innercolor
    innersizes = repeat([neuron_innersize], length(neuronpositions)) .+ 
        neuroncounter_scaled .* 0.15
    meshscatter!(
        parentscene, neuronpositions, 
        markersize=innersizes, color=innercolors
    )
    
    outercolors = repeat([neuron_outercolor], length(neuronpositions))
    outercolors[activeindices] .= active_outercolor
    outersizes = repeat([neuron_outersize], length(neuronpositions)) .+ 
        neuroncounter_scaled .* 0.15
    meshscatter!(
        parentscene, neuronpositions, 
        markersize=outersizes, color=outercolors
    )
end

function updatepositions(
    state::State,
    sides::Int,
    neuronsnumber::Int,
    sideranges::Vector{LinRange},
    sidelength_int::Int,
    randomscale::Float64
)
    sideoffset = 0.035sidelength_int
    lastside = 2sideoffset + sidelength_int

    neuronid = 1
    for i in 1:sides
        if i == 1
            for x=sideranges[1], y=sideranges[2]
                xpos = x + sideoffset + rand(Uniform(-1, 1)) * randomscale
                ypos = y + sideoffset + rand(Uniform(-1, 1)) * randomscale
                state.neuronpositions[neuronid] = Point3f(xpos, ypos, 0)
                neuronid += 1
            end
        elseif i == 2
            for x=sideranges[3], z=sideranges[4]
                xpos = x + sideoffset + rand(Uniform(-1, 1)) * randomscale
                zpos = z + sideoffset + rand(Uniform(-1, 1)) * randomscale
                state.neuronpositions[neuronid] = Point3f(xpos, 0, zpos)
                neuronid += 1
            end
        elseif i == 3
            for y=sideranges[5], z=sideranges[6]
                ypos = y + sideoffset + rand(Uniform(-1, 1)) * randomscale
                zpos = z + sideoffset + rand(Uniform(-1, 1)) * randomscale
                state.neuronpositions[neuronid] = Point3f(0, ypos, zpos)
                neuronid += 1
            end
        elseif i == 4
            for x=sideranges[7], y=sideranges[8]
                xpos = x + sideoffset + rand(Uniform(-1, 1)) * randomscale
                ypos = y + sideoffset + rand(Uniform(-1, 1)) * randomscale
                state.neuronpositions[neuronid] = Point3f(xpos, ypos, lastside)
                neuronid += 1
            end
        elseif i == 5
            for x=sideranges[9], z=sideranges[10]
                xpos = x + sideoffset + rand(Uniform(-1, 1)) * randomscale
                zpos = z + sideoffset + rand(Uniform(-1, 1)) * randomscale
                state.neuronpositions[neuronid] = Point3f(xpos, lastside, zpos)
                neuronid += 1
            end
        elseif i == 6
            for y=sideranges[11], z=sideranges[12]
                neuronid > neuronsnumber && break

                ypos = y + sideoffset + rand(Uniform(-1, 1)) * randomscale
                zpos = z + sideoffset + rand(Uniform(-1, 1)) * randomscale
                state.neuronpositions[neuronid] = Point3f(lastside, ypos, zpos)
                neuronid += 1
            end
        end
    end
end

function createscenes(resolution, camera3d, background)
    if background == :black
        set_theme!(theme_black(), resolution=resolution)
    else
        set_theme!(theme_light(), resolution=resolution)
    end
    GLMakie.activate!(; framerate=120, fullscreen=true)

    figure = Figure()
    rowsize!(figure.layout, 1, Fixed(resolution[2]))

    al = AmbientLight(RGBf(1.0, 1.0, 1.0))

    parentscene = LScene(
        figure[:, 1],
        show_axis=false,
        scenekw = (
            clear=true,
            lights=[al],
            backgroundcolor=background
        )
    )

    camera = if camera3d 
        camera = cam3d!(parentscene.scene)
        camera.attributes.reset[] = Keyboard.m
        camera
    else 
        cam2d!(parentscene.scene)
    end
    # camc = cameracontrols(parentscene.scene)
    # update_cam!(parentscene.scene, camc, Vec3f(0, 5, 5), Vec3f(0.0, 0, 0))

    figure, parentscene, camera
end

function primary_resolution()
    monitor = GLMakie.GLFW.GetPrimaryMonitor()
    videomode = GLMakie.MonitorProperties(monitor).videomode
    return (videomode.width, videomode.height)
end

function minmax(array::AbstractArray)
    min = minimum(array)
    max = maximum(array)
    range = max - min
    (array .- min) ./ range
end

function sparsematrix_plot(file::String=ANAKG_SAMPLE_MAT_FILE)
    data = anakgload(file)
    for i in axes(data)[1]
        data[i, i] = 0.0
    end

    f, ax, plt = spy(
        sparse(data), 
        markersize = 4, marker = :circle, framecolor = :lightgrey, colormap=:RdBu
    )

    hidedecorations!(ax) # remove axis labeling
    ax.title = "Visualization of a ANAKG adjacency matrix"

    f
end
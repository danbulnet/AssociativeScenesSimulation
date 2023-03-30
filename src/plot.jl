export anakgplot

using GLMakie
using Random
using Distributions
using Colors
using LinearAlgebra

using AssociativeScenesSimulation

mutable struct State
    neuronpositions::Dict{UInt, Point3f}

    State() = new(Dict())
end

state = State()

function anakgplot(
    file::String=ANAKG_SAMPLE_MAT_FILE;
    resolution=primary_resolution(), camera3d=true, background=:azure1,
    randomscale=0.3,
    neuron_outercolor=RGBA{Float32}(0.53, 0.81, 1.0, 0.5), 
    neuron_innercolor=RGBA{Float32}(0.0, 0.0, 0.55, 0.88),
    neuron_outersize=0.2 , neuron_innersize=0.07
)
    global state = State()

    set_theme!(theme_black(), resolution=resolution)
    GLMakie.activate!(; framerate=120, fullscreen=true)

    figure, parentscene, camera = createscenes(
        resolution, camera3d, background
    )
    
    sides = 6
    data = anakgload(file)
    neuronsnumber = size(data, 1)
    neuronsnumber_perside = neuronsnumber / sides
    sidelength = √(neuronsnumber_perside)
    sidelength_int = Int(floor(sidelength))
    println("sidelength_int $sidelength_int")
    lengthdiff = neuronsnumber - ceil(sidelength_int^2 * sides)
    println("lengthdiff $lengthdiff")
    
    lengthdiff_int = lengthdiff ÷ (sidelength_int * sides)
    println("lengthdiff_int $lengthdiff_int")
    lengthdiff_rem = lengthdiff % sidelength_int
    println("lengthdiff_rem $lengthdiff_rem")
    sideranges = []
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

    # siderangex, siderangey, siderangez, lengthdiff_rem = if lengthdiff > 3sides
    #     lengthdiff_int = lengthdiff ÷ 3sides
    #     range = LinRange(
    #         0.0, 
    #         sidelength_int + lengthdiff_int, 
    #         sidelength_int + lengthdiff_int
    #     )
    #     range, range, range, lengthdiff % 3sides
    # else
    #     lengthdiff_int = lengthdiff ÷ sides 
    #     range = LinRange(0.0, sidelength_int, sidelength_int)
    #     zrange = LinRange(
    #         0.0, 
    #         sidelength_int + lengthdiff_int, 
    #         sidelength_int + lengthdiff_int
    #     )
    #     range, range, zrange, lengthdiff % 3sides
    # end
    # siderangex = LinRange(0.0, sidelength_int, sidelength_int)
    # siderangey = LinRange(0.0, sidelength_int, sidelength_int)
    # siderangez = LinRange(0.0, sidelength_int, sidelength_int)

    sideoffset = 0.05sidelength_int
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

    neuronpositions = collect(values(state.neuronpositions))
    neuroncounter = data[collect(diagind(data))]
    neuroncounter_scaled = minmax(neuroncounter)
    
    innercolors = repeat([neuron_innercolor], length(neuronpositions))
    innersizes = repeat([neuron_innersize], length(neuronpositions)) .+ 
        neuroncounter_scaled .* 0.15
    meshscatter!(
        parentscene, neuronpositions, 
        markersize=innersizes, color=innercolors
    )
    
    outercolors = repeat([neuron_outercolor], length(neuronpositions))
    outersizes = repeat([neuron_outersize], length(neuronpositions)) .+ 
        neuroncounter_scaled .* 0.15
    meshscatter!(
        parentscene, neuronpositions, 
        markersize=outersizes, color=outercolors
    )

    figure
end

function createscenes(resolution, camera3d, background)
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
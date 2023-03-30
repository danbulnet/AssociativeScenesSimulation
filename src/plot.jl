export anakgplot

using GLMakie

using AssociativeScenesSimulation

function anakgplot(
    file::String=ANAKG_SAMPLE_MAT_FILE;
    resolution=primary_resolution(), camera3d=true
)
    set_theme!(theme_black(), resolution=resolution)
    GLMakie.activate!(; framerate=120, fullscreen=true)

    figure, parentscene, camera = createscenes(resolution, camera3d)
    
    data = anakgload(file)
    neuronsnumber = size(data, 1)
    sidelength = ceil(√(neuronsnumber / 6))
    siderange = LinRange(0.0, sidelength, Int(sidelength))

    xs = Float64[]
    ys = Float64[]
    zs = Float64[]
    for x=siderange, y=siderange
        push!(xs, x)
        push!(ys, y)
        push!(zs, 0)
    end

    println("typeof(xs), $(typeof(xs))")

    # xs = cos.(1:0.5:20)
    # ys = sin.(1:0.5:20)
    # zs = LinRange(0, 3, length(xs))

    meshscatter!(parentscene, xs, ys, zs, markersize = 0.1, color = zs)

    figure
end

function createscenes(resolution, camera3d)
    figure = Figure()
    rowsize!(figure.layout, 1, Fixed(resolution[2]))

    # pl = PointLight(Point3f(0, 0, 15), RGBf(1.0, 0.98, 0.94))
    al = AmbientLight(RGBf(0.58, 0.58, 0.58))

    parentscene = LScene(
        figure[:, 1],
        show_axis=false,
        scenekw = (
            clear=true,
            lights=[al],
            backgroundcolor=:black,
            # ssao = Makie.SSAO(radius=250.0, blur=2, bias=1),
            # lightposition = Vec3f(0, 0, 15),
            shininess=256f0
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
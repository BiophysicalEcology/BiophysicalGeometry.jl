module BiophysicalGeometryAppExt

# The "Build an animal" app: controls in the browser, the package in Julia. Every change writes the Julia code for
# the animal with `AnimalBuilder.animal_code`, evaluates it, and shows what the package makes of it.
# `using BiophysicalGeometry, Bonito, WGLMakie; BiophysicalGeometry.app()`.

using BiophysicalGeometry
using BiophysicalGeometry: CompositeBody, covered_areas
using BiophysicalGeometry.AnimalBuilder: DEFAULTS, PRESETS, animal_code, build_animal
using Bonito
using WGLMakie
using Unitful
using WGLMakie.Colors: hex
import BiophysicalGeometry: app

# ── Controls ──────────────────────────────────────────────────────────────────────────────────────────────────────
#
# A heading, or (setting, label, values): a range of numbers is a slider, strings a drop-down, a Bool a checkbox.

const MASSES = unique(round.(10 .^ (-2:0.01:4); sigdigits = 3))   # a log scale
const CONTROLS = [
    "Body",
    (:mass, "Mass (kg)", MASSES),
    (:posture, "Posture", ["Horizontal", "Upright"]),
    (:pitch, "Body pitch, head up (°)", -30:1.0:80),
    (:torsoShape, "Torso shape", ["Cylinder", "Ellipsoid"]),
    (:torsoRatio, "Torso length / width", 1.1:0.1:6),
    (:fat, "Fat, fraction of torso mass", 0:0.005:0.5),
    "Coat",
    (:backFur, "Depth on the back (m)", 0:0.0005:0.06),
    (:bellyFur, "Depth on the belly (m)", 0:0.0005:0.06),
    (:limbFur, "Depth on head and legs (m)", 0:0.0005:0.06),
    "Head",
    (:headShape, "Shape", ["None", "Sphere", "Ellipsoid"]),
    (:headRatio, "Length / width", 1.0:0.1:3),
    (:headFraction, "Fraction of mass", 0.01:0.0001:0.25),
    (:neck, "Neck", false),
    (:neckFraction, "Neck, fraction of mass", 0.01:0.005:0.15),
    (:neckRatio, "Neck length / width", 0.5:0.1:8),
    (:neckPosture, "Neck posture", ["Forward", "Up"]),
    (:nose, "Nose or beak", ["None", "Nose", "Beak"]),
    (:noseFraction, "Nose, fraction of mass", 0.001:0.001:0.03),
    (:ears, "Ears", ["None", "Cone", "Plate"]),
    (:earFraction, "Ear, fraction of mass, each", 0.0005:0.0005:0.03),
    (:earPosture, "Plate ears", ["Up", "Flat"]),
    (:earRatio, "Plate ear length / width", 0.5:0.1:4),
    (:earFlatness, "Plate ear length / thickness", 3:1.0:40),
    "Legs",
    (:legs, "Number", ["0", "2", "4"]),
    (:legScaling, "Proportions", ["Manual", "Elastic", "Geometric"]),
    (:legFraction, "Fraction of mass, each (Manual)", 0.0005:0.0005:0.2),
    (:legRatio, "Length / width (Manual)", 1:0.1:12),
    (:legTop, "Taper, foot / top (1: cylinder)", 0.1:0.05:1),
    (:hindLegs, "Hind legs (four legs, Manual)", ["Same", "Different"]),
    (:hindFraction, "Hind leg, fraction of mass, each", 0.005:0.005:0.2),
    (:hindRatio, "Hind leg, length / width", 1:0.1:12),
    "Arms (upright)",
    (:arms, "Arms", false),
    (:armFraction, "Fraction of mass, each", 0.005:0.0001:0.1),
    (:armRatio, "Length / width", 2:0.5:16),
    "Wings and tail (horizontal)",
    (:wings, "Wings", ["None", "Folded", "Spread"]),
    (:wingFraction, "Wing, fraction of mass, each", 0.005:0.005:0.15),
    (:tail, "Tail", false),
    (:tailFraction, "Tail, fraction of mass", 0.001:0.001:0.25),
    (:tailRatio, "Tail length / width", 1:0.5:20),
]

nearest(values, x) = argmin(abs.(values .- x))

# A slider's stops include every preset's value, so that each preset is reproduced exactly.
stops(key, values) = sort(unique(vcat(collect(Float64, values), [Float64(last(p)[key]) for p in PRESETS])))
widget(key, values::AbstractVector{<:Real}, x) = (v = stops(key, values); Bonito.Slider(v; value = v[nearest(v, x)]))
widget(key, options::Vector{String}, x) = Bonito.Dropdown(options; index = something(findfirst(==(string(x)), options), 1))
widget(key, ::Bool, x) = Bonito.Checkbox(Bool(x))

value(w::Bonito.Slider) = w.value[]
value(w::Bonito.Dropdown) = w.value[]
value(w::Bonito.Checkbox) = w.value[]

function set!(w::Bonito.Slider, x)
    i = nearest(w.values[], x)
    w.index[] = i
    w.value[] = w.values[][i]
end
set!(w::Bonito.Dropdown, x) = (w.option_index[] = something(findfirst(==(string(x)), w.options[]), 1))
set!(w::Bonito.Checkbox, x) = (w.value[] = Bool(x))

shown(w::Bonito.Slider, key) = key == :mass ? weight(w.value[]) : string(round(w.value[]; sigdigits = 4))
shown(w, key) = ""

# The settings: the preset last chosen, for what has no control (such as density), changed by the controls.
function settings(base, widgets)
    s = merge(base, (; (key => value(w) for (key, w) in widgets)...))
    return merge(s, (; legs = parse(Int, s.legs)))
end
apply!(widgets, preset) = foreach(((key, w),) -> set!(w, preset[key]), widgets)

# ── What the package makes of the animal ──────────────────────────────────────────────────────────────────────────

fmt(x) = x >= 100 ? string(round(Int, x)) : x >= 1 ? string(round(x; sigdigits = 3)) : string(round(x; sigdigits = 3))
area(a) = (x = ustrip(u"m^2", a); x >= 0.1 ? "$(fmt(x)) m²" : "$(fmt(x * 1e4)) cm²")
weight(x) = x >= 1 ? "$(fmt(x)) kg" : "$(fmt(x * 1000)) g"

# The group a part belongs to in the table: its name, with the numbered pairs taken together.
group(name) = (s = string(name); for g in ("leg", "ear", "wing", "arm"); startswith(s, g) && return g * "s"; end; s)

function numbers(animal::CompositeBody, s)
    cov = covered_areas(animal.parts, animal.joins)
    groups = Dict{String,Vector{Any}}()
    order = String[]
    for (name, part) in pairs(animal.parts)
        g = group(name)
        haskey(groups, g) || (push!(order, g); groups[g] = [0.0u"m^2", 0.0u"kg"])
        groups[g][1] += total_area(part) - cov[name]
        groups[g][2] += BiophysicalGeometry.mass(part.shape)
    end
    total = total_area(animal)
    summary = [
        "Outer area" => area(total),
        "Skin area" => area(skin_area(animal)),
        "Hidden by joins" => area(sum(values(cov))),
        "Volume" => "$(fmt(1000 * s.mass / s.density)) L",
        "Meeh coefficient" => string(round(ustrip(u"m^2", total) / s.mass^(2 / 3); digits = 3)),
    ]
    cell(c::Tuple) = DOM.td(c...)
    cell(c) = DOM.td(c)
    row(cells...; head = false) = DOM.tr((head ? DOM.th(c) : cell(c) for c in cells)...)
    DOM.div(
        DOM.table((row(k, v) for (k, v) in summary)...),
        DOM.table(row("Part", "Exposed area", "Mass"; head = true),
            (row((DOM.span("■ "; style = "color: #$(hex(get(PART_COLOURS, g, RGBf(0.6, 0.6, 0.6))));"), g),
                 area(groups[g][1]), weight(ustrip(u"kg", groups[g][2]))) for g in order)...);
        style = "display: flex; gap: 24px; font-size: 13px;")
end

# The two figures are made once per session and redrawn in place, so the view keeps its rotation.
function body_view()
    fig = Figure(size = (620, 440), backgroundcolor = :white)
    ax = Axis3(fig[1, 1]; perspectiveness = 0.3, viewmode = :fitzoom, aspect = :data,
        elevation = π / 7, azimuth = 5π / 4, xlabel = "x (cm)", ylabel = "y (cm)", zlabel = "z (cm)")
    return fig, ax
end

# A colour for each part, the pairs (legs, ears, wings, arms) sharing one.
const PART_COLOURS = Dict(
    "dorsal" => RGBf(0.30, 0.55, 0.75), "ventral" => RGBf(0.90, 0.62, 0.20), "head" => RGBf(0.35, 0.68, 0.45),
    "neck" => RGBf(0.58, 0.45, 0.72), "nose" => RGBf(0.55, 0.42, 0.35), "beak" => RGBf(0.85, 0.65, 0.25),
    "tail" => RGBf(0.50, 0.50, 0.50), "ears" => RGBf(0.85, 0.55, 0.75), "wings" => RGBf(0.39, 0.67, 0.71),
    "arms" => RGBf(0.63, 0.75, 0.35), "legs" => RGBf(0.80, 0.40, 0.45))
part_colour(name) = get(PART_COLOURS, group(name), RGBf(0.6, 0.6, 0.6))

# The outer surface of every part, posed, as one mesh coloured by part: the meshes the package rasterises for the
# silhouette, in cm.
function body_mesh(animal)
    vertices, faces, colours = Point3f[], Int[], RGBf[]
    for (name, part) in pairs(animal.parts)
        pose = getfield(animal.poses, name)
        for grid in BiophysicalGeometry.part_outer_meshes(part.shape, part, 100.0)
            add_grid!(vertices, faces, colours, BiophysicalGeometry.transform_mesh(grid..., pose, 100.0)..., part_colour(name))
        end
    end
    return vertices, permutedims(reshape(faces, 3, :)), colours
end

function add_grid!(vertices, faces, colours, X::Matrix{Float64}, Y::Matrix{Float64}, Z::Matrix{Float64}, colour)
    n0 = length(vertices)
    nr, nc = size(X)
    for i in eachindex(X)
        push!(vertices, Point3f(X[i], Y[i], Z[i]))
        push!(colours, colour)
    end
    index(i, j) = n0 + (j - 1) * nr + i
    for j in 1:nc-1, i in 1:nr-1
        append!(faces, (index(i, j), index(i + 1, j), index(i, j + 1), index(i + 1, j), index(i + 1, j + 1), index(i, j + 1)))
    end
end

# The animal in its part colours, or with every part cut open where it faces the viewer, to show its layers.
function draw_body!(ax, animal, layers::Bool = false)
    empty!(ax)
    if layers
        draw_cutaway!(ax, animal)
    else
        vertices, faces, colours = body_mesh(animal)
        mesh!(ax, vertices, faces; color = colours)
    end
    reset_limits!(ax)
end

function sun_direction(zenith, azimuth)
    z, a = deg2rad(zenith), deg2rad(azimuth)
    return (sin(z) * cos(a), sin(z) * sin(a), cos(z))
end

const SHADOW_RESOLUTION = 120

function shadow_view()
    fig = Figure(size = (300, 260), backgroundcolor = :white)
    ax = Axis(fig[1, 1]; aspect = DataAspect(), titlesize = 13, xlabel = "cm", ylabel = "cm")
    n = SHADOW_RESOLUTION
    image = (x = Observable(LinRange(0.0, 1.0, n)), y = Observable(LinRange(0.0, 1.0, n)),
             z = Observable(zeros(Float32, n, n)))
    heatmap!(ax, image.x, image.y, image.z; colormap = [:white, RGBAf(0.2, 0.2, 0.25, 1.0)], colorrange = (0.0f0, 1.0f0))
    return fig, ax, image
end

function draw_shadow!(ax, image, animal, direction)
    r = silhouette_rasterized(animal, direction; resolution = SHADOW_RESOLUTION, return_image = true)
    # The image keeps its size, so each of its three parts can be updated on its own.
    image.x[] = LinRange(100 .* r.x_range..., size(r.bitmap, 1))
    image.y[] = LinRange(100 .* r.y_range..., size(r.bitmap, 2))
    image.z[] = Float32.(r.bitmap)
    ax.title = "Silhouette to the sun: $(area(r.area))"
    limits!(ax, 100 .* r.x_range..., 100 .* r.y_range...)
end

# Run `f` for the latest request only: requests made while it runs are folded into one more run when it ends, so
# that a dragged slider never queues up a backlog of stale animals.
function latest(f)
    busy, again = Ref(false), Ref(false)
    return function ()
        busy[] && (again[] = true; return)
        busy[] = true
        @async try
            while true
                again[] = false
                f()
                again[] || break
            end
        catch e
            @error "Build an animal" exception = (e, catch_backtrace())
        finally
            busy[] = false
        end
        return
    end
end

# ── The app ───────────────────────────────────────────────────────────────────────────────────────────────────────

const NOTE = "font-size: 12px; color: #555;"
labelled(text, widget, readout = "") = DOM.div(DOM.label(text; style = "display: block; font-weight: 600; margin-top: 6px; font-size: 13px;"),
    DOM.div(widget, " ", DOM.span(readout; style = NOTE)))

function make_app()
    App(; title = "Build an animal · BiophysicalGeometry.jl") do session
        start = last(first(PRESETS))
        widgets = [c[1] => widget(c[1], c[3], start[c[1]]) for c in CONTROLS if c isa Tuple]
        base = Ref(start)
        lookup = Dict(widgets)
        preset = Bonito.Dropdown(collect(first.(PRESETS)); index = 1)
        zenith = Bonito.Slider(collect(0.0:1:90); value = 30.0)
        azimuth = Bonito.Slider(collect(0.0:5:360); value = 90.0)
        layers = Bonito.Checkbox(false)

        status = Observable("")
        code = Observable("")
        body_fig, body_ax = body_view()
        shadow_fig, shadow_ax, image = shadow_view()
        table_slot = Observable{Any}(DOM.div())
        animal = Ref{Any}(nothing)
        applying = Ref(false)

        function draw_shadow()
            animal[] === nothing && return
            draw_shadow!(shadow_ax, image, animal[], sun_direction(zenith.value[], azimuth.value[]))
        end
        function rebuild()
            applying[] && return
            s = settings(base[], widgets)
            try
                text = animal_code(s)
                a = build_animal(text)
                animal[] = a
                code[] = text
                draw_body!(body_ax, a, layers.value[])
                table_slot[] = numbers(a, s)
                draw_shadow()
                status[] = ""
            catch e
                status[] = "This animal cannot be built: " * first(sprint(showerror, e), 300)
            end
        end

        request_rebuild = latest(rebuild)
        request_shadow = latest(draw_shadow)
        for (_, w) in widgets
            on(_ -> request_rebuild(), w.value)
        end
        onany((_...) -> request_shadow(), zenith.value, azimuth.value)
        # Cut open, the cut follows the view as it turns.
        request_redraw = latest(() -> animal[] === nothing || draw_body!(body_ax, animal[], layers.value[]))
        on(_ -> request_redraw(), layers.value)
        onany((_...) -> layers.value[] && request_redraw(), body_ax.azimuth, body_ax.elevation)
        on(preset.option_index) do i
            applying[] = true
            try
                base[] = last(PRESETS[i])
                apply!(widgets, base[])
            finally
                applying[] = false
            end
            request_rebuild()
        end
        rebuild()

        controls = Any[DOM.h3("Build an animal"), labelled("Start from", preset)]
        for c in CONTROLS
            if c isa String
                push!(controls, DOM.h4(c; style = "margin: 12px 0 0 0;"))
            else
                key, label = c
                w = lookup[key]
                push!(controls, labelled(label, w, map(_ -> shown(w, key), w.value)))
            end
        end
        push!(controls, DOM.h4("Sun"; style = "margin: 12px 0 0 0;"),
            labelled("Zenith angle (°)", zenith, zenith.value), labelled("Azimuth, 0° is head-on (°)", azimuth, azimuth.value))
        left = DOM.div(controls...; style = "width: 320px; min-width: 320px; padding: 10px; height: 92vh; overflow-y: auto; font-family: sans-serif;")

        copy = DOM.button("Copy code"; onclick = js"""() => navigator.clipboard.writeText($(code).value)""")
        right = DOM.div(
            DOM.p(status; style = "color: #b00; font-weight: 600;"),
            DOM.div(layers, " Show fat and fur: cut each part open"; style = "font-size: 13px; margin-bottom: 4px;"),
            DOM.div(body_fig, shadow_fig; style = "display: flex; gap: 12px; flex-wrap: wrap; align-items: flex-start;"),
            table_slot,
            DOM.h4("The Julia code for this animal"), copy,
            DOM.p("The animal above was built by running exactly this code with BiophysicalGeometry.jl."; style = NOTE),
            DOM.pre(code; style = "font-size: 11px; background: #f6f6f6; padding: 8px; white-space: pre-wrap;");
            style = "flex: 1; min-width: 0; padding: 10px; font-family: sans-serif;")
        return DOM.div(left, right; style = "display: flex;")
    end
end

function app(; port::Integer = 8080, open::Bool = true, host::AbstractString = "127.0.0.1", proxy_url = nothing)
    # Build every preset once, so that the first change in the browser does not wait on compilation.
    for (_, preset) in PRESETS
        a = build_animal(animal_code(preset))
        numbers(a, preset)
    end
    # and draw one, so that the plotting is compiled too
    a = build_animal(animal_code(last(first(PRESETS))))
    draw_body!(last(body_view()), a)
    draw_body!(last(body_view()), a, true)
    _, ax, image = shadow_view()
    draw_shadow!(ax, image, a, sun_direction(30, 90))
    server = proxy_url === nothing ? Bonito.Server(make_app(), host, port) :
        Bonito.Server(make_app(), host, port; proxy_url)
    url = Bonito.online_url(server, "/")
    @info "Build an animal: running at $url. Close the server with `close(server)`."
    open && Bonito.HTTPServer.openurl(url)
    return server
end

end

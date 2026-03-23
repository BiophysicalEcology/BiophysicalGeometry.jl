# dog_widget.jl — interactive GLMakie widget for the dog organism
# Run: julia --project=examples examples/dog_widget.jl

using GLMakie, BiophysicalGeometry, Unitful, Printf
import BiophysicalGeometry: Cylinder, Ellipsoid, Plate

const DENSITY = 990.0u"kg/m^3"
const FAT_ρ   = 900.0u"kg/m^3"

function build_dog(; mass_kg,
                    head_w, torso_w, legs_w, tail_w, ears_w,
                    fur_mm, fat_pct,
                    head_b, torso_b, legs_b, tail_b, ears_b, ears_c)
    mass = mass_kg * u"kg"
    total_w     = head_w + torso_w + legs_w + tail_w + ears_w
    head_f      = head_w  / total_w
    torso_f     = torso_w / total_w
    leg_f_each  = legs_w  / total_w / 4
    tail_f      = tail_w  / total_w
    ears_f_each = ears_w  / total_w / 2

    fur = Fur(fur_mm * 1e-3 * u"m", 5e-5u"m", 1e8u"m^-2")
    fat = Fat(fat_pct / 100, FAT_ρ)

    Organism(
        parts = (
            head      = BodyPart(Body(Ellipsoid(mass * head_f,     DENSITY, head_b,  1.0),     CompositeInsulation(fur, fat)), 1),
            torso     = BodyPart(Body(Cylinder( mass * torso_f,    DENSITY, torso_b),           CompositeInsulation(fur, fat)), 1),
            legs_fore = BodyPart(Body(Cylinder( mass * leg_f_each, DENSITY, legs_b),            CompositeInsulation(fur, fat)), 2),
            legs_hind = BodyPart(Body(Cylinder( mass * leg_f_each, DENSITY, legs_b),            CompositeInsulation(fur, fat)), 2),
            tail      = BodyPart(Body(Cylinder( mass * tail_f,     DENSITY, tail_b),            CompositeInsulation(fur, fat)), 1),
            ears      = BodyPart(Body(Plate(    mass * ears_f_each,DENSITY, ears_b,  ears_c),   CompositeInsulation(fur, fat)), 2),
        ),
        joins = (
            Join(:torso, :head,      :forward, 0.5),
            Join(:torso, :legs_fore, :lateral, 0.0, π/2,  0.7),
            Join(:torso, :legs_hind, :lateral, 0.0, π/2, -0.7),
            Join(:torso, :tail,      :back,    0.0),
            Join(:head,  :ears,      :top,     0.0),
        )
    )
end

fig = Figure(size=(1500, 700), backgroundcolor=:white)

ax = Axis(fig[1:2, 1]; xlabel="metres", ylabel="metres", aspect=DataAspect(),
          title="Dog (25 kg)  —  Prone")

# sl[1]     = total mass
# sl[2..6]  = part weights
# sl[7]     = fur depth
# sl[8]     = body fat
# sl[9..14] = shape ratios
sg = SliderGrid(fig[1:2, 2],
    (label="Total mass (kg)",          range=1.0:0.5:80.0,     startvalue=25.0,  format=x -> @sprintf("%.1f kg", x)),
    (label="Head weight",              range=0.01:0.01:0.5,    startvalue=0.10,  format=x -> @sprintf("%.2f", x)),
    (label="Torso weight",             range=0.1:0.01:2.0,     startvalue=0.52,  format=x -> @sprintf("%.2f", x)),
    (label="Legs weight (×4)",         range=0.01:0.01:1.0,    startvalue=0.36,  format=x -> @sprintf("%.2f", x)),
    (label="Tail weight",              range=0.005:0.005:0.1,  startvalue=0.01,  format=x -> @sprintf("%.3f", x)),
    (label="Ears weight (×2)",         range=0.001:0.001:0.02, startvalue=0.003, format=x -> @sprintf("%.3f", x)),
    (label="Fur depth (mm)",           range=1.0:0.5:50.0,     startvalue=20.0,  format=x -> @sprintf("%.1f mm", x)),
    (label="Body fat (%)",             range=0.0:0.5:40.0,     startvalue=10.0,  format=x -> @sprintf("%.1f %%", x)),
    (label="Head a/b ratio",           range=0.5:0.1:5.0,      startvalue=2.0,   format=x -> @sprintf("%.1f", x)),
    (label="Torso L/D ratio",          range=0.5:0.1:10.0,     startvalue=3.5,   format=x -> @sprintf("%.1f", x)),
    (label="Legs L/D ratio",           range=1.0:0.1:15.0,     startvalue=6.0,   format=x -> @sprintf("%.1f", x)),
    (label="Tail L/D ratio",           range=1.0:0.5:20.0,     startvalue=8.0,   format=x -> @sprintf("%.1f", x)),
    (label="Ear length/width",         range=0.1:0.1:2.0,      startvalue=0.6,   format=x -> @sprintf("%.1f", x)),
    (label="Ear width (rel. depth)",   range=1.0:0.5:10.0,     startvalue=2.5,   format=x -> @sprintf("%.1f", x)),
    width=400, tellheight=false)

sl = sg.sliders

function redraw!()
    empty!(ax)
    try
        mass_kg = sl[1].value[]
        dog = build_dog(
            mass_kg = mass_kg,
            head_w  = sl[2].value[],
            torso_w = sl[3].value[],
            legs_w  = sl[4].value[],
            tail_w  = sl[5].value[],
            ears_w  = sl[6].value[],
            fur_mm  = sl[7].value[],
            fat_pct = sl[8].value[],
            head_b  = sl[9].value[],
            torso_b = sl[10].value[],
            legs_b  = sl[11].value[],
            tail_b  = sl[12].value[],
            ears_b  = sl[13].value[],
            ears_c  = sl[14].value[],
        )
        draw_organism_prone!(ax, dog)
        autolimits!(ax)
        ax.title = @sprintf("Dog (%.1f kg)  —  Prone", mass_kg)
        geom = geometry(dog)
        text!(ax, @sprintf("SA = %.0f cm²", ustrip(u"cm^2", geom.total_area));
              position=(0.02, 0.02), space=:relative, align=(:left, :bottom),
              fontsize=14, font=:bold)
    catch e
        @warn "Could not draw: $e"
    end
end

for s in sl
    on(s.value) do _; redraw!(); end
end

redraw!()
display(fig)
println("Close the window to exit.")
wait(fig.scene)

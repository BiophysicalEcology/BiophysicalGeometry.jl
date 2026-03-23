# human_widget.jl — interactive GLMakie widget for the human organism
# Run: julia --project=examples examples/human_widget.jl

using GLMakie, BiophysicalGeometry, Unitful, Printf
import BiophysicalGeometry: Cylinder, Ellipsoid

const DENSITY = 1050.0u"kg/m^3"
const FAT_ρ   = 900.0u"kg/m^3"

function build_human(; mass_kg, head_w, torso_w, arms_w, legs_w,
                       hair_mm, fur_body_mm, fat_pct,
                       head_b, torso_b, arms_b, legs_b)
    mass = mass_kg * u"kg"
    total_w     = head_w + torso_w + arms_w + legs_w
    head_f      = head_w  / total_w
    torso_f     = torso_w / total_w
    arms_f_each = arms_w  / total_w / 2
    legs_f_each = legs_w  / total_w / 2

    fur_head = Fur(hair_mm     * 1e-3 * u"m", 7.5e-5u"m", 3e8u"m^-2")
    fur_body = Fur(fur_body_mm * 1e-3 * u"m",   1e-6u"m", 3e8u"m^-2")
    fat      = Fat(fat_pct / 100, FAT_ρ)

    Organism(
        parts = (
            head  = BodyPart(Body(Ellipsoid(mass * head_f,       DENSITY, head_b,  1.0),
                                  CompositeInsulation(fur_head, fat)), 1),
            torso = BodyPart(Body(Cylinder( mass * torso_f,      DENSITY, torso_b),
                                  CompositeInsulation(fur_body, fat)), 1),
            arms  = BodyPart(Body(Cylinder( mass * arms_f_each,  DENSITY, arms_b),
                                  CompositeInsulation(fur_body, fat)), 2),
            legs  = BodyPart(Body(Cylinder( mass * legs_f_each,  DENSITY, legs_b),
                                  CompositeInsulation(fur_body, fat)), 2),
        ),
        joins = (
            Join(:torso, :head, :top,     0.15),
            Join(:torso, :arms, :lateral, 0.0),
            Join(:torso, :legs, :bottom,  0.0),
        )
    )
end

fig = Figure(size=(1500, 800), backgroundcolor=:white)

ax = Axis(fig[1:2, 1]; xlabel="metres", ylabel="metres", aspect=DataAspect(),
          title="Human (70 kg)  —  Upright")

# sl[1]     = total mass
# sl[2..5]  = part weights
# sl[6..8]  = hair/beard/body fur
# sl[9]     = body fat
# sl[10..13]= shape ratios
sg = SliderGrid(fig[1:2, 2],
    (label="Total mass (kg)",       range=20.0:1.0:200.0,   startvalue=70.0,   format=x -> @sprintf("%.0f kg", x)),
    (label="Head weight",           range=0.01:0.005:0.3,   startvalue=0.076,  format=x -> @sprintf("%.3f", x)),
    (label="Torso weight",          range=0.1:0.01:2.0,     startvalue=0.501,  format=x -> @sprintf("%.3f", x)),
    (label="Arms weight (×2)",      range=0.01:0.005:0.4,   startvalue=0.098,  format=x -> @sprintf("%.3f", x)),
    (label="Legs weight (×2)",      range=0.02:0.01:1.0,    startvalue=0.324,  format=x -> @sprintf("%.3f", x)),
    (label="Hair depth (mm)",       range=0.0:0.5:30.0,     startvalue=10.0,   format=x -> @sprintf("%.1f mm", x)),
    (label="Beard depth (mm)",      range=0.0:0.5:30.0,     startvalue=0.0,    format=x -> @sprintf("%.1f mm", x)),
    (label="Body fur depth (mm)",   range=0.5:0.5:20.0,     startvalue=6.0,    format=x -> @sprintf("%.1f mm", x)),
    (label="Body fat (%)",          range=0.0:0.5:50.0,     startvalue=18.0,   format=x -> @sprintf("%.1f %%", x)),
    (label="Head a/b ratio",        range=0.5:0.1:3.0,      startvalue=1.6,    format=x -> @sprintf("%.1f", x)),
    (label="Torso L/D ratio",       range=0.5:0.1:6.0,      startvalue=1.9,    format=x -> @sprintf("%.1f", x)),
    (label="Arms L/D ratio",        range=2.0:0.5:25.0,     startvalue=12.0,   format=x -> @sprintf("%.1f", x)),
    (label="Legs L/D ratio",        range=2.0:0.5:20.0,     startvalue=7.0,    format=x -> @sprintf("%.1f", x)),
    width=400, tellheight=false)

sl = sg.sliders

function redraw!()
    empty!(ax)
    try
        mass_kg = sl[1].value[]
        human = build_human(
            mass_kg     = mass_kg,
            head_w      = sl[2].value[],
            torso_w     = sl[3].value[],
            arms_w      = sl[4].value[],
            legs_w      = sl[5].value[],
            hair_mm     = sl[6].value[],
            fur_body_mm = sl[8].value[],
            fat_pct     = sl[9].value[],
            head_b      = sl[10].value[],
            torso_b     = sl[11].value[],
            arms_b      = sl[12].value[],
            legs_b      = sl[13].value[],
        )
        draw_organism_upright!(ax, human;
            beard_depth = sl[7].value[] * 1e-3 * u"m")
        autolimits!(ax)
        ax.title = @sprintf("Human (%.0f kg)  —  Upright", mass_kg)
        geom = geometry(human)
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

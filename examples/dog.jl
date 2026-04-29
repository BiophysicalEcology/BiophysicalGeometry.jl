# dog.jl
#
# Build a dog from primitives and plot it.
#
# Body: dorsal + ventral half-cylinders (different insulation) joined on
# their flat faces with FullCover.
# Legs: four cylinders attached to the ventral lateral surface as Discs.
# Head: an ellipsoid attached to one end-cap of the dorsal half.
#
# Run from the repository root:
#   julia --project=examples examples/dog.jl

using BiophysicalGeometry
using Unitful
using GLMakie

import BiophysicalGeometry: Cylinder, Sphere, Ellipsoid, HalfCylinder

GLMakie.activate!()

# ── Parameters ────────────────────────────────────────────────────────────────
density = 1000.0u"kg/m^3"
b_body  = 3.0    # body length-to-diameter ratio
b_leg   = 5.0
b_head  = 1.5

m_total_body = 20.0u"kg"      # full body mass; each half is half of this
m_leg_cyl    = 1.0u"kg"       # the cylinder leg this frustum replaces
leg_top_ratio = 0.4           # foot radius / shoulder radius
# A frustum with the same base radius and length as that cylinder has volume
# (1 + t + t²)/3 of the cylinder, so scale mass by the same factor.
m_leg        = m_leg_cyl * (1 + leg_top_ratio + leg_top_ratio^2) / 3
m_head       = 2.0u"kg"

dorsal_fur  = Fur(15.0u"mm", 30.0u"μm", 3000u"cm^-2")
ventral_fur = Fur( 5.0u"mm", 30.0u"μm", 3000u"cm^-2")  # thinner belly fur
limb_fur    = Fur( 8.0u"mm", 30.0u"μm", 3000u"cm^-2")  # legs + head

# ── Parts ─────────────────────────────────────────────────────────────────────
dorsal  = Body(HalfCylinder(m_total_body / 2, density, b_body), dorsal_fur)
ventral = Body(HalfCylinder(m_total_body / 2, density, b_body), ventral_fur)
leg     = Body(Cone(m_leg, density, b_leg, leg_top_ratio), limb_fur)

# Truncate the head's pole_a end so it joins the body flush rather than
# touching at a single point. We pick the truncation so the cut-disc radius
# equals the join Disc radius (head_disc = 0.7 · b_minor → cut disc has
# y/z extents 0.7 · b_minor; that means pole_a_truncation = 1 - sqrt(1 - 0.49)).
head_truncate = 1 - sqrt(1 - 0.7^2)
head    = Body(Ellipsoid(m_head, density, b_head, 1.0, head_truncate), limb_fur)

L_body = dorsal.geometry.length.length_skin
r_body = dorsal.geometry.length.radius_skin
r_leg  = skin_radius(leg)

# Lay the body horizontal. The HalfCylinder local frame has axis along +z and
# bulge along +y, so we want:
#   - dorsal local +z (length axis)  → world +x  (head/tail direction)
#   - dorsal local +y (dorsal bulge) → world +z  (dorsal side up)
# Columns of the rotation matrix are the world-frame images of the local axes.
R_horizontal = [0.0  0.0  1.0;
                1.0  0.0  0.0;
                0.0  1.0  0.0]
root_pose = Pose((0.0u"m", 0.0u"m", 0.0u"m"), R_horizontal)

# ── Composite ─────────────────────────────────────────────────────────────────
# Place legs on the ventral half (they hang from the underside).
# Ventral's :lateral has φ ∈ [0, π]; φ = π/2 is the centre of its bulge,
# which after the dorsal/ventral join points *down* in world space.
# Splay the legs slightly: φ = π/2 ± 0.4 rad.
leg_z_front = 0.20 * L_body
leg_z_back  = 0.80 * L_body
leg_φ_left  = π/2 + 0.35
leg_φ_right = π/2 - 0.35

# The head sits on the dorsal :end_b (the +x world end after rotation).
# The disc must fit on the dorsal half-disc (radius r_body) AND on the head's
# pole, where the notional surface area is π·b² (head minor radius).
head_disc = 0.7 * min(r_body, head.geometry.length.b_semi_minor_skin)

dog = CompositeBody(;
    parts = (; dorsal, ventral, leg_fl=leg, leg_fr=leg, leg_bl=leg, leg_br=leg, head),
    joins = (
        # Dorsal/ventral split. The twist around the joint axis is the 6th DOF
        # left free by the surface-normal alignment; -π/2 keeps ventral's
        # length-axis parallel to dorsal's instead of rotated by 90°.
        Join(:dorsal,  Attachment(:flat, (;), FullCover()),
             :ventral, Attachment(:flat, (;), FullCover());
             twist=-π/2),
        # Head onto dorsal's far end-cap, centred.
        Join(:dorsal, Attachment(:end_b,  (r=0.0u"m", φ=0.0), Disc(head_disc)),
             :head,   Attachment(:pole_a, (;),                Disc(head_disc))),
        # Four legs onto ventral's lateral.
        Join(:ventral, Attachment(:lateral, (z=leg_z_front, φ=leg_φ_left),  Disc(r_leg)),
             :leg_fl,  Attachment(:base,    (r=0.0u"m", φ=0.0),             Disc(r_leg))),
        Join(:ventral, Attachment(:lateral, (z=leg_z_front, φ=leg_φ_right), Disc(r_leg)),
             :leg_fr,  Attachment(:base,    (r=0.0u"m", φ=0.0),             Disc(r_leg))),
        Join(:ventral, Attachment(:lateral, (z=leg_z_back,  φ=leg_φ_left),  Disc(r_leg)),
             :leg_bl,  Attachment(:base,    (r=0.0u"m", φ=0.0),             Disc(r_leg))),
        Join(:ventral, Attachment(:lateral, (z=leg_z_back,  φ=leg_φ_right), Disc(r_leg)),
             :leg_br,  Attachment(:base,    (r=0.0u"m", φ=0.0),             Disc(r_leg))),
    ),
    root = :dorsal,
    root_pose = root_pose,
)

# ── Report ────────────────────────────────────────────────────────────────────
# Rasterised silhouette accounts for part-on-part shadowing; the per-part
# summed `silhouette_area(dog)` is shown alongside as the upper bound.
silhouette_top  = silhouette_rasterized(dog, (0.0, 0.0, 1.0))
silhouette_side = silhouette_rasterized(dog, (1.0, 0.0, 0.0))
println("Dog built.")
println("  total surface area      : ", total_area(dog))
println("  flesh volume            : ", flesh_volume(dog))
println("  silhouette top  (sun↓z) : ", silhouette_top)
println("  silhouette side (sun→x) : ", silhouette_side)
println("  per-part summed (upper) : ", silhouette_area(dog).normal)

# ── Plot ──────────────────────────────────────────────────────────────────────
# Interactive silhouette: 3-D body on the left, rasterised projection on the
# right, with sliders for sun zenith θ and azimuth φ.
fig = plot_body_silhouette(dog; resolution=128)
screen = display(fig)
println()
println("Drag the θ / φ sliders to rotate the sun direction.")
println("Rotate the 3-D view with left-click drag, zoom with scroll.")
println("Close the window to exit.")
wait(screen)

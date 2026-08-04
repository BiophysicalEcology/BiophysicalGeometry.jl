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

dorsal_fibrous  = FibrousLayer(15.0u"mm", 30.0u"μm", 3000u"cm^-2")
ventral_fibrous = FibrousLayer( 5.0u"mm", 30.0u"μm", 3000u"cm^-2")  # thinner belly fur
limb_fibrous    = FibrousLayer( 8.0u"mm", 30.0u"μm", 3000u"cm^-2")  # legs + head

# ── Parts ─────────────────────────────────────────────────────────────────────
dorsal_body  = Body(HalfCylinder(m_total_body / 2, density, b_body), dorsal_fibrous)
ventral_body = Body(HalfCylinder(m_total_body / 2, density, b_body), ventral_fibrous)
leg_body     = Body(Cone(m_leg, density, b_leg, leg_top_ratio), limb_fibrous)

# Truncate the head's pole_a end so it joins the body flush rather than
# touching at a single point. We pick the truncation so the cut-disc radius
# equals the join Disc radius (head_disc = 0.7 · b_minor → cut disc has
# y/z extents 0.7 · b_minor; that means pole_a_truncation = 1 - sqrt(1 - 0.49)).
head_truncate = 1 - sqrt(1 - 0.7^2)
head_body = Body(Ellipsoid(m_head, density, b_head, 1.0, head_truncate), limb_fibrous)

L_body = dorsal_body.geometry.length.length_skin
r_body = dorsal_body.geometry.length.radius_skin
r_leg  = skin_radius(leg_body)

# Lay the body horizontal. The HalfCylinder local frame has axis along +z and
# bulge along +y, so we want:
# - dorsal local +z (length axis)  → world +x  (head/tail direction)
# - dorsal local +y (dorsal bulge) → world +z  (dorsal side up)
# Columns of the rotation matrix are the world-frame images of the local axes.
R_horizontal = [0.0  0.0  1.0;
                1.0  0.0  0.0;
                0.0  1.0  0.0]
root_pose = Pose((0.0u"m", 0.0u"m", 0.0u"m"), R_horizontal)

# ── Composite ─────────────────────────────────────────────────────────────────
# Place legs on the ventral half (they hang from the underside).
# Ventral's Lateral surface has φ ∈ [0, π]; φ = π/2 is the centre of its bulge,
# which after the dorsal/ventral join points *down* in world space.
# Splay the legs slightly: φ = π/2 ± 0.4 rad.
leg_z_front = 0.20 * L_body
leg_z_back  = 0.80 * L_body
leg_φ_left  = π/2 + 0.35
leg_φ_right = π/2 - 0.35

# The head sits on the dorsal EndB (the +x world end after rotation).
# The disc must fit on the dorsal half-disc (radius r_body) AND on the head's
# pole, where the notional surface area is π·b² (head minor radius).
head_disc = 0.7 * min(r_body, head_body.geometry.length.b_semi_minor_skin)

dog = CompositeBody(;
    parts = (;
        dorsal = dorsal_body,
        ventral = ventral_body,
        head = head_body,
        leg_fl = leg_body,
        leg_fr = leg_body,
        leg_bl = leg_body,
        leg_br = leg_body,
    ),
    joins = (
        # Dorsal/ventral split. The twist around the joint axis is the 6th DOF
        # left free by the surface-normal alignment; -π/2 keeps ventral's
        # length-axis parallel to dorsal's instead of rotated by 90°.
        Join(dorsal = Attachment(Flat(), FullCover()),
             ventral = Attachment(Flat(), FullCover()); twist=-π/2),
        # Head onto dorsal's far end-cap, centred.
        Join(dorsal = Attachment(EndB(0.0u"m", 0.0), Disc(head_disc)),
             head = Attachment(PoleA(), Disc(head_disc))),
        # Four legs onto ventral's lateral. On the Cone, EndA is the base.
        Join(ventral = Attachment(Lateral(leg_z_front, leg_φ_left), Disc(r_leg)),
             leg_fl = Attachment(EndA(0.0u"m", 0.0), Disc(r_leg))),
        Join(ventral = Attachment(Lateral(leg_z_front, leg_φ_right), Disc(r_leg)),
             leg_fr = Attachment(EndA(0.0u"m", 0.0), Disc(r_leg))),
        Join(ventral = Attachment(Lateral(leg_z_back, leg_φ_left), Disc(r_leg)),
             leg_bl = Attachment(EndA(0.0u"m", 0.0), Disc(r_leg))),
        Join(ventral = Attachment(Lateral(leg_z_back, leg_φ_right), Disc(r_leg)),
             leg_br = Attachment(EndA(0.0u"m", 0.0), Disc(r_leg))),
    ),
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

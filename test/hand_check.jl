#
# Hand-derived reference values for every test in test/geometry.jl.
#
# Every formula here is written from first principles (no reference to the package
# implementation) and the values are then compared with the numbers pinned in the
# test file. Anything flagged "MISMATCH" means either (a) my derivation is wrong
# or (b) the package implementation diverges from the standard formula and the
# test values would lock in that divergence.
#
# Run from the project root:
#     julia --project=. test/hand_check.jl
#

using Unitful

# === Inputs (mirror test/geometry.jl) ===

mass = 65.0u"kg"
density = 1000.0u"kg/m^3"
axis_ratio_b = 5.0     # for plate: L/W; for cylinder: L/(2R); for ellipsoid: a/b
axis_ratio_c = 5.0     # for plate: L/H
fibrous_thickness = 10.0u"mm"
fibre_diameter = 30.0u"μm"
fibre_density_per_area = 3000u"cm^-2"
fat_fraction = 0.1
fat_density = 901.0u"kg/m^3"

T_fibrous = uconvert(u"m", fibrous_thickness)
D_fibrous = uconvert(u"m", fibre_diameter)
N_fibrous = uconvert(u"m^-2", fibre_density_per_area)

V       = uconvert(u"m^3", mass / density)                       # body volume
V_fat   = uconvert(u"m^3", mass * fat_fraction / fat_density)    # fat volume
V_flesh = V - V_fat

hair_area(skin) = π * (D_fibrous/2)^2 * N_fibrous * skin                 # fibrous cross-section over skin

θ = uconvert(u"rad", 30.0u"°")

# === Sphere: V = 4/3 π R³ ===

R       = uconvert(u"m", cbrt(3V       / (4π)))
R_fibrous   = R + T_fibrous
R_flesh = uconvert(u"m", cbrt(3V_flesh / (4π)))
fat_sphere = R - R_flesh

S_sphere(r)   = 4π * r^2          # closed sphere surface
sil_sphere(r) = π * r^2           # disc projection (orientation-independent)

# === Cylinder (closed, axis_ratio_b = L / (2R), so L = 2 b R) ===
# V = π R² L = 2π b R³

Rc       = uconvert(u"m", cbrt(V       / (2π * axis_ratio_b)))
Lc       = 2 * axis_ratio_b * Rc
Rc_fibrous   = Rc + T_fibrous
Lc_fibrous   = Lc + 2 * T_fibrous
Rc_flesh = uconvert(u"m", cbrt(V_flesh / (2π * axis_ratio_b)))
fat_cyl  = Rc - Rc_flesh

S_cyl(r, l)              = 2π*r*l + 2π*r^2
sil_cyl_normal(r, l)     = 2 * r * l            # body broadside to sun (rectangle face)
sil_cyl_parallel(r)      = π * r^2              # body end-on to sun (end cap)
sil_cyl_zenith(r, l, θ)  = 2 * r * l * sin(θ) + π * r^2 * cos(θ)
                                                # body axis horizontal, θ from vertical:
                                                # θ=0 sun overhead → end cap; θ=π/2 → side

# === Plate (axis_ratio_b = L/W, axis_ratio_c = L/H) ===
# V = L × W × H = L × L/b × L/c = L³/(b c)

Lp       = uconvert(u"m", cbrt(V * axis_ratio_b * axis_ratio_c))
Wp       = Lp / axis_ratio_b
Hp       = Lp / axis_ratio_c
Lp_fibrous   = Lp + 2 * T_fibrous
Wp_fibrous   = Wp + 2 * T_fibrous
Hp_fibrous   = Hp + 2 * T_fibrous
Wp_flesh = uconvert(u"m", cbrt(V_flesh * axis_ratio_b * axis_ratio_c)) / axis_ratio_b
fat_plate = (Wp - Wp_flesh) / 2

S_plate(L, W, H) = 2 * (L*W + L*H + W*H)

# Plate silhouette: largest face = NormalToSun, smallest = ParallelToSun
plate_faces(L, W, H) = (L*W, L*H, W*H)

# === Ellipsoid (prolate spheroid: c = b, a = b × axis_ratio_b) ===
# V = 4/3 π a b c = 4/3 π axis_ratio_b b³

b_e     = uconvert(u"m", cbrt(3V / (4π * axis_ratio_b)))
c_e     = b_e
a_e     = axis_ratio_b * b_e
b_efibrous  = b_e + T_fibrous
c_efibrous  = c_e + T_fibrous
a_efibrous  = a_e + T_fibrous

# Prolate spheroid (c = b < a) surface area: S = 2πb² + 2π (a b / e) asin(e),  e = √(a² - b²)/a
function S_prolate(a, b, c)
    aₘ = ustrip(u"m", a); bₘ = ustrip(u"m", b); cₘ = ustrip(u"m", c)
    e = sqrt(aₘ^2 - cₘ^2) / aₘ
    return (2π * bₘ^2 + 2π * (aₘ * bₘ / e) * asin(e)) * u"m^2"
end

# Ellipsoid silhouette area (shadow on plane perpendicular to view direction n̂):
#   A = π × abc × √( (n_x/a)² + (n_y/b)² + (n_z/c)² )
# Body's long axis (a) along x; sun zenith θ from vertical → n̂ = (sin θ, 0, cos θ).
sil_ell_normal(a, b)       = π * a * b                        # n̂ = (0,0,1) → πab
sil_ell_parallel(b, c)     = π * b * c                        # n̂ = (1,0,0) → πbc
sil_ell_zenith(a, b, c, θ) = π * a * b * c * sqrt((sin(θ)/a)^2 + (cos(θ)/c)^2)

# === Desert iguana / Leopard frog: allometric in mass (g) ===

mg = ustrip(u"g", uconvert(u"g", mass))

di_total    = uconvert(u"m^2", 10.4713 * mg^0.688 * u"cm^2")
di_ventral  = uconvert(u"m^2",  0.425  * mg^0.85  * u"cm^2")
di_normal   = uconvert(u"m^2",  3.798  * mg^0.683 * u"cm^2")
di_parallel = uconvert(u"m^2",  0.694  * mg^0.743 * u"cm^2")
di_intermed = (di_normal + di_parallel) / 2
# Body radius from cylinder approx with L = 4R: V = 4π R³
di_r = uconvert(u"m", cbrt(V / (4π)))

lf_total = uconvert(u"m^2", 12.79 * mg^0.606 * u"cm^2")
lf_r     = uconvert(u"m", cbrt(V / (4π)))

# === Comparison harness ===

const RESULTS = Tuple{String, Any, Any, Float64, Symbol}[]

function check(name, computed, expected)
    diff = abs(computed - expected)
    rel = ustrip(diff) == 0 ? 0.0 : ustrip(diff / abs(expected))
    status = rel < 1e-10 ? :exact : (rel < 1e-6 ? :match : (rel < 1e-3 ? :close : :MISMATCH))
    push!(RESULTS, (name, computed, expected, rel, status))
end

# --- Plate / Naked ---
faces_p_skin = plate_faces(Lp, Wp, Hp)
check("Plate/Naked total_area",       S_plate(Lp, Wp, Hp),    1.2163304590093516u"m^2")
check("Plate/Naked skin_radius",      Wp/2,                   0.11756673438603786u"m")
check("Plate/Naked sil_normal",       max(faces_p_skin...),   0.27643874068394353u"m^2")
check("Plate/Naked sil_parallel",     min(faces_p_skin...),   0.05528774813678871u"m^2")
check("Plate/Naked sil_intermediate", (max(faces_p_skin...) + min(faces_p_skin...))/2,
                                                              0.16586324441036612u"m^2")
check("Plate/Naked volume",           V,                       0.065u"m^3")

# --- Plate / FibrousLayer ---
S_plate_skin = S_plate(Lp, Wp, Hp)
S_plate_fibrous  = S_plate(Lp_fibrous, Wp_fibrous, Hp_fibrous)
faces_p_fibrous  = plate_faces(Lp_fibrous, Wp_fibrous, Hp_fibrous)
check("Plate/Fibrous total_area",         S_plate_fibrous,            1.3504052015217138u"m^2")
check("Plate/Fibrous skin_area",          S_plate_skin,           1.2163304590093516u"m^2")
check("Plate/Fibrous evaporation_area",   S_plate_skin - hair_area(S_plate_skin),
                                                              1.190537258877413u"m^2")
check("Plate/Fibrous insulation_radius",  Wp_fibrous/2,               0.12756673438603786u"m")
check("Plate/Fibrous sil_normal",         max(faces_p_fibrous...),    0.3050547569365926u"m^2")
check("Plate/Fibrous sil_parallel",       min(faces_p_fibrous...),    0.06509308688767174u"m^2")

# --- Plate / FatLayer ---
check("Plate/Fat flesh_radius",       Wp/2 - fat_plate,       0.11304560873468857u"m")
check("Plate/Fat flesh_volume",       V_flesh,                0.057785793562708104u"m^3")

# --- Cylinder / Naked ---
check("Cyl/Naked total_area",         S_cyl(Rc, Lc),          1.1222291434482228u"m^2")
check("Cyl/Naked skin_radius",        Rc,                     0.12742495668987025u"m")
check("Cyl/Naked sil_normal",         sil_cyl_normal(Rc, Lc), 0.32474239174830616u"m^2")
check("Cyl/Naked sil_parallel",       sil_cyl_parallel(Rc),   0.05101041561128287u"m^2")
check("Cyl/Naked sil_intermediate",   (sil_cyl_normal(Rc,Lc) + sil_cyl_parallel(Rc))/2,
                                                              0.18787640367979452u"m^2")
check("Cyl/Naked sil_zenith(30°)",    sil_cyl_zenith(Rc, Lc, θ), 0.20654751165112634u"m^2")

# --- Cylinder / FibrousLayer ---
S_cyl_skin = S_cyl(Rc, Lc)
S_cyl_fibrous  = S_cyl(Rc_fibrous, Lc_fibrous)
check("Cyl/Fibrous total_area",           S_cyl_fibrous,              1.2362029452302272u"m^2")
check("Cyl/Fibrous skin_area",            S_cyl_skin,             1.1222291434482228u"m^2")
check("Cyl/Fibrous evaporation_area",     S_cyl_skin - hair_area(S_cyl_skin),
                                                              1.098431432327489u"m^2")
check("Cyl/Fibrous insulation_radius",    Rc_fibrous,                 0.13742495668987026u"m")
check("Cyl/Fibrous sil_normal",           sil_cyl_normal(Rc_fibrous, Lc_fibrous), 0.355724381353875u"m^2")
check("Cyl/Fibrous sil_zenith(30°)",      sil_cyl_zenith(Rc_fibrous, Lc_fibrous, θ),
                                                              0.22924427552149568u"m^2")

# --- Cylinder / FatLayer ---
check("Cyl/Fat flesh_radius",         Rc_flesh,               0.12252472497618691u"m")

# --- Sphere / Naked ---
check("Sphere/Naked total_area",      S_sphere(R),            0.7817952526648283u"m^2")
check("Sphere/Naked skin_radius",     R,                      0.24942591981125847u"m")
check("Sphere/Naked sil_normal",      sil_sphere(R),          0.19544881316620707u"m^2")
check("Sphere/Naked sil_zenith(30°)", sil_sphere(R),          0.19544881316620707u"m^2")

# --- Sphere / FibrousLayer ---
check("Sphere/Fibrous total_area",        S_sphere(R_fibrous),        0.8457394607097782u"m^2")
check("Sphere/Fibrous skin_area",         S_sphere(R),            0.7817952526648283u"m^2")
check("Sphere/Fibrous evaporation_area",  S_sphere(R) - hair_area(S_sphere(R)),
                                                              0.7652166976637417u"m^2")
check("Sphere/Fibrous insulation_radius", R_fibrous,                  0.25942591981125845u"m")
check("Sphere/Fibrous sil_normal",        sil_sphere(R_fibrous),      0.21143486517744456u"m^2")

# --- Sphere / FatLayer ---
check("Sphere/Fat flesh_radius",      R_flesh,                0.2398340405261943u"m")

# --- Ellipsoid / Naked ---
check("Ell/Naked total_area",         S_prolate(a_e, b_e, c_e), 1.0679282565991794u"m^2")
check("Ell/Naked skin_radius",        b_e,                    0.14586516277963593u"m")
check("Ell/Naked sil_normal",         sil_ell_normal(a_e, b_e),    0.3342127693207218u"m^2")
check("Ell/Naked sil_parallel",       sil_ell_parallel(b_e, c_e),  0.06684255386414435u"m^2")
# Ellipsoid zenith silhouette is INTENTIONALLY not checked here: the package
# implementation uses a 2-angle Euler projection algorithm with the body-
# orientation angle hardcoded to 90°, and the resulting convention does not
# correspond to "horizontal prolate ellipsoid, sun zenith θ from vertical".
# For reference:
#   - textbook horizontal-ellipsoid projection at θ=30°  ≈ 0.2914 m² (Naked)
#                                                        ≈ 0.3158 m² (Fibrous)
#   - package implementation at θ=30°                    ≈ 0.0767 m² (Naked)
#                                                        ≈ 0.0875 m² (Fibrous)
# See the comment block on `silhouette_area(::Ellipsoid, a, b, c, θ)` in
# src/shapes/ellipsoid.jl. Once the orientation convention is clarified
# (e.g. against the original Fortran/NicheMapR source), reinstate one of
# these as a real first-principles check.

# --- Ellipsoid / FibrousLayer ---
S_ell_skin = S_prolate(a_e, b_e, c_e)
S_ell_fibrous  = S_prolate(a_efibrous, b_efibrous, c_efibrous)
check("Ell/Fibrous total_area",           S_ell_fibrous,              1.1587845516526847u"m^2")
check("Ell/Fibrous skin_area",            S_ell_skin,             1.0679282565991794u"m^2")
check("Ell/Fibrous evaporation_area",     S_ell_skin - hair_area(S_ell_skin),
                                                              1.045282036532102u"m^2")
check("Ell/Fibrous insulation_radius",    b_efibrous,                 0.15586516277963594u"m")
check("Ell/Fibrous sil_normal",           sil_ell_normal(a_efibrous, b_efibrous), 0.3620218640142718u"m^2")
# Ell/Fibrous zenith silhouette: same divergence — see note above.

# --- DesertIguana / Naked ---
check("DI/Naked total_area",          di_total,               2.144297264544543u"m^2")
check("DI/Naked skin_radius",         di_r,                   0.1729422736164134u"m")
check("DI/Naked sil_normal",          di_normal,              0.7358254132114781u"m^2")
check("DI/Naked sil_parallel",        di_parallel,            0.26142920031318184u"m^2")
check("DI/Naked sil_intermediate",    di_intermed,            0.49862730676232997u"m^2")

# --- LeopardFrog / Naked ---
check("LF/Naked total_area",          lf_total,               1.055591874642603u"m^2")
check("LF/Naked skin_radius",         lf_r,                   0.1729422736164134u"m")

# === Print results ===

println("="^110)
println(rpad("name", 38), " | ", rpad("computed", 24), " | ", rpad("expected", 24), " | rel error  | status")
println("-"^110)
for (name, computed, expected, rel, status) in RESULTS
    sym = status === :exact    ? "✓ exact"  :
          status === :match    ? "✓ match"  :
          status === :close    ? "≈ close"  :
          "✗ MISMATCH"
    println(rpad(name, 38), " | ", rpad(string(computed), 24), " | ",
            rpad(string(expected), 24), " | ", rpad(string(rel), 10), " | ", sym)
end
println("="^110)

n_exact = count(r -> r[5] === :exact, RESULTS)
n_match = count(r -> r[5] === :match, RESULTS)
n_close = count(r -> r[5] === :close, RESULTS)
n_bad   = count(r -> r[5] === :MISMATCH, RESULTS)
println("Total: $(length(RESULTS)) — exact: $n_exact, match: $n_match, close: $n_close, mismatch: $n_bad")

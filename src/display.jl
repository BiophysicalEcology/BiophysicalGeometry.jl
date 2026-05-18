const _W = 22  # field-label column width

# ── Marker types ──────────────────────────────────────────────────────────────

Base.show(io::IO, ::MIME"text/plain", ::Naked)              = print(io, "Naked()")
Base.show(io::IO, ::MIME"text/plain", ::NormalToSun)        = print(io, "NormalToSun()")
Base.show(io::IO, ::MIME"text/plain", ::ParallelToSun)      = print(io, "ParallelToSun()")
Base.show(io::IO, ::MIME"text/plain", ::Intermediate)       = print(io, "Intermediate()")
Base.show(io::IO, ::MIME"text/plain", ::ZenithAngleVarying) = print(io, "ZenithAngleVarying()")

# ── Shape types ───────────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", s::Sphere)
    println(io, "Sphere")
    println(io, rpad("  mass:", _W), s.mass)
    print(io,   rpad("  density:", _W), s.density)
end

function Base.show(io::IO, ::MIME"text/plain", s::Cylinder)
    println(io, "Cylinder")
    println(io, rpad("  mass:", _W), s.mass)
    println(io, rpad("  density:", _W), s.density)
    print(io,   rpad("  axis_ratio_b:", _W), s.axis_ratio_b)
end

function Base.show(io::IO, ::MIME"text/plain", s::Ellipsoid)
    println(io, "Ellipsoid")
    println(io, rpad("  mass:", _W), s.mass)
    println(io, rpad("  density:", _W), s.density)
    println(io, rpad("  axis_ratio_b:", _W), s.axis_ratio_b)
    print(io,   rpad("  axis_ratio_c:", _W), s.axis_ratio_c)
end

function Base.show(io::IO, ::MIME"text/plain", s::Plate)
    println(io, "Plate")
    println(io, rpad("  mass:", _W), s.mass)
    println(io, rpad("  density:", _W), s.density)
    println(io, rpad("  axis_ratio_b:", _W), s.axis_ratio_b)
    print(io,   rpad("  axis_ratio_c:", _W), s.axis_ratio_c)
end

function Base.show(io::IO, ::MIME"text/plain", s::LeopardFrog)
    println(io, "LeopardFrog")
    println(io, rpad("  mass:", _W), s.mass)
    print(io,   rpad("  density:", _W), s.density)
end

function Base.show(io::IO, ::MIME"text/plain", s::DesertIguana)
    println(io, "DesertIguana")
    println(io, rpad("  mass:", _W), s.mass)
    print(io,   rpad("  density:", _W), s.density)
end

# ── Insulation types ──────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", f::FibrousLayer)
    println(io, "FibrousLayer")
    println(io, rpad("  thickness:", _W), f.thickness)
    println(io, rpad("  fibre_diameter:", _W), f.fibre_diameter)
    print(io,   rpad("  fibre_density:", _W), f.fibre_density)
end

function Base.show(io::IO, ::MIME"text/plain", f::FatLayer)
    println(io, "FatLayer")
    println(io, rpad("  fraction:", _W), f.fraction)
    print(io,   rpad("  density:", _W), f.density)
end

function Base.show(io::IO, ::MIME"text/plain", c::CompositeInsulation)
    n = length(c.layers)
    println(io, "CompositeInsulation ($n layer$(n == 1 ? "" : "s"))")
    for (i, layer) in enumerate(c.layers)
        print(io, "  [$i] ")
        show(io, MIME"text/plain"(), layer)
        i < n && println(io)
    end
end

# ── Internal helper ───────────────────────────────────────────────────────────

function _show_area_lines(io, a::SurfaceAreas, indent, w)
    items = Pair{String,Any}["total" => a.total]
    a.skin !== a.total       && push!(items, "skin"       => a.skin)
    a.convection !== a.total && push!(items, "convection" => a.convection)
    !isnothing(a.ventral)    && push!(items, "ventral"    => a.ventral)
    for (i, (k, v)) in enumerate(items)
        i < length(items) ? println(io, rpad("$indent$k:", w), v) : print(io, rpad("$indent$k:", w), v)
    end
end

# ── SurfaceAreas ──────────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", a::SurfaceAreas)
    println(io, "Surface Areas")
    _show_area_lines(io, a, "  ", _W)
end

# ── Geometry ──────────────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", g::Geometry)
    println(io, "Geometry")
    println(io, rpad("  volume:", _W), g.volume)
    println(io, "  Dimensions:")
    for (k, v) in pairs(g.length)
        println(io, rpad("    $k:", _W + 2), v)
    end
    println(io, "  Surface Areas:")
    _show_area_lines(io, g.area, "    ", _W + 2)
end

# ── Body ──────────────────────────────────────────────────────────────────────

function Base.show(io::IO, ::MIME"text/plain", b::Body)
    header = "Body{$(nameof(typeof(b.shape))), $(nameof(typeof(b.insulation)))}"
    println(io, header)
    println(io, "─" ^ length(header))

    print(io, "Shape:      ")
    show(io, MIME"text/plain"(), b.shape)
    println(io)
    println(io)

    print(io, "Insulation: ")
    show(io, MIME"text/plain"(), b.insulation)
    println(io)
    println(io)

    println(io, "Geometry:")
    println(io, rpad("  volume:", _W), b.geometry.volume)
    println(io, "  Dimensions:")
    for (k, v) in pairs(b.geometry.length)
        println(io, rpad("    $k:", _W + 2), v)
    end
    println(io, "  Surface Areas:")
    _show_area_lines(io, b.geometry.area, "    ", _W + 2)
end

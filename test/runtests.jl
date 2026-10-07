using BiophysicalGeometry
using Aqua
using SafeTestsets
using Test

@testset "Quality assurance" begin
    Aqua.test_unbound_args(BiophysicalGeometry)
    Aqua.test_stale_deps(BiophysicalGeometry)
    Aqua.test_undefined_exports(BiophysicalGeometry)
    Aqua.test_project_extras(BiophysicalGeometry)
    Aqua.test_deps_compat(BiophysicalGeometry)
end

@safetestset "geometry" begin include("geometry.jl") end
@safetestset "composition" begin include("composition.jl") end
@safetestset "shapes" begin include("shapes.jl") end
# Whisk needs Julia 1.12, and running the module needs node.
if VERSION >= v"1.12" && Sys.which("node") !== nothing
    @safetestset "wasm" begin include("wasm.jl") end
else
    @info "Skipping the wasm tests: they need Julia 1.12 and node"
end

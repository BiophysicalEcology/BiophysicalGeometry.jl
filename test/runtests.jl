using BiophysicalGeometry
using SafeTestsets
using Aqua

Aqua.test_unbound_args(BiophysicalGeometry)
Aqua.test_stale_deps(BiophysicalGeometry)
Aqua.test_undefined_exports(BiophysicalGeometry)
Aqua.test_project_extras(BiophysicalGeometry)
Aqua.test_deps_compat(BiophysicalGeometry)
Aqua.test_project_toml_formatting(BiophysicalGeometry)

@safetestset "geometry" begin include("geometry.jl") end
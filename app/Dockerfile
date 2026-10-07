# The "Build an animal" app. Build from the repository root:
#
#   docker build -f app/Dockerfile -t animal-builder .
#   docker run -p 8080:8080 animal-builder

FROM julia:1.11

WORKDIR /BiophysicalGeometry.jl
COPY Project.toml ./
COPY src ./src
COPY ext ./ext
COPY app/Project.toml ./app/

ENV JULIA_PROJECT=/BiophysicalGeometry.jl/app
RUN julia -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

COPY app/run.jl ./app/
ENV HOST=0.0.0.0 PORT=8080 OPEN=false
EXPOSE 8080
CMD ["julia", "app/run.jl"]

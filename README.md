# PARSECFoils.jl

Julia package to handle and convert geometries in the PARSEC airfoil parameterization.

## Installation

```
]add PARSECFoils
```

## Usage

To use PARSECFoils, you may resort to the following functions, along with their usage examples in `examples/test.ipynb`:

### Convert PARSEC parameters to `(x, t, y)` coordinates (`y` standing for camber @ chordwise control points)

```
    function parsec2xty( # pre-defined chordwise point grid
        x::AbstractVector;
        leading_edge_radius::Real = 0.01,
        x_trailing_edge::Real = 1.0,
        y_trailing_edge::Real = 0.0,
        trailing_edge_thickness::Real = 0.0,
        α::Real = 0.0,
        β::Real = deg2rad(18.0),
        x_upper::Real = 0.35,
        y_upper::Real = 0.06,
        upper_curvature::Real = -0.5,
        x_lower::Real = 0.35,
        y_lower::Real = 0.06,
        lower_curvature::Real = -0.5,
    )
```

### Convert airfoil from `(x, t, y)` format to Selig format

```
    function xty2selig(
        x::AbstractVector, t::AbstractVector, y::AbstractVector
    )
```

### Convert Selig-format airfoil to xty format (using a pre-defined chordwise point grid `xs`)

```
    function selig2xty(
        xs::AbstractVector, pts::AbstractMatrix
    )
```

### Convert PARSEC definitions to Selig format airfoil with pre-defined chordwise point grid `xs`

```
    parsec2selig(xs::AbstractVector; kwargs...) = xty2selig(
        parsec2xty(xs; kwargs...)...
    )
```

### Use the BFGS optimization algorithm to find the PARSEC parameters that best fit an airfoil in `(x, t, y)` format.

`kwargs` are passed to Optim.jl.

Returns a dictionary of parameters.

```
    function xty2parsec(
        x::AbstractVector, t::AbstractVector, y::AbstractVector;
        kwargs...
    )
```

### Use the BFGS optimization algorithm to find the PARSEC parameters that best fit an airfoil in Selig format.

`kwargs` are passed to Optim.jl

```
    selig2parsec(
        x::AbstractVector, pts::AbstractMatrix;
        kwargs...
    ) = xty2parsec(
        selig2xty(x, pts)...;
        kwargs...
    )
```

module PARSECFoils

    export parsec2xty, xty2selig, selig2xty, parsec2selig, selig2parsec, xty2parsec

    using LinearAlgebra
    using Interpolations

    using Optim

    """
    Obtain coefficients for a given airfoil surface (intra or extra)
    """
    function pcoeff(
        rle::Real,
        xte::Real, 
        yte::Real,
        x_cre::Real,
        y_cre::Real,
        ypp::Real,
        theta::Real;
        intra::Bool = false,
    )

        c1 = sqrt(2 * rle)

        # if intra, invert signals
        if intra
            yte = - yte
            y_cre = - y_cre
            ypp = - ypp
            theta = - theta
            c1 = - c1
        end

        A = [
            xte ^ 1.5 xte ^ 2.5 xte ^ 3.5 xte ^ 4.5 xte ^ 5.5;
            x_cre ^ 1.5 x_cre ^ 2.5 x_cre ^ 3.5 x_cre ^ 4.5 x_cre ^ 5.5;
            1.5 * sqrt(xte) 2.5 * xte ^ 1.5 3.5 * xte ^ 2.5 4.5 * xte ^ 3.5 5.5 * xte ^ 4.5;
            1.5 * sqrt(x_cre) 2.5 * x_cre ^ 1.5 3.5 * x_cre ^ 2.5 4.5 * x_cre ^ 3.5 5.5 * x_cre ^ 4.5;
            0.75 / sqrt(x_cre) 3.75 * sqrt(x_cre) 8.75 * x_cre ^ 1.5 15.75 * x_cre ^ 2.5 24.75 * x_cre ^ 3.5
        ]

        b = [
            yte - c1 * sqrt(xte),
            y_cre - c1 * sqrt(x_cre),
            tan(theta) - 0.5 * c1 / sqrt(xte),
            - 0.5 * c1 / sqrt(x_cre),
            ypp + 0.25 * c1 * x_cre ^ (-1.5)
        ]

        [
            c1;
            (A \ b)
        ]

    end

    """
    Obtain PARSEC coefficients at a given set of x axis positions
    """
    function polyval(
        x,
        c::AbstractVector
    )

        sum(
            t -> let (i, cc) = t 
                (
                    @. x ^ (i - 0.5) * cc
                )
            end,
            enumerate(c)
        )

    end

    const parsec_parameters = [
        :leading_edge_radius,
        :y_trailing_edge,
        :α,
        :β,
        :x_upper,
        :y_upper,
        :upper_curvature,
        :x_lower,
        :y_lower,
        :lower_curvature,
        :trailing_edge_thickness,
    ]

    const opt_x0 = [
        0.01, 0.0, 0.0, deg2rad(18.0), 0.35, 0.06, -0.5, 0.35, 0.06, -0.5, 0.01
    ]

    const lower_bounds = [
        0.001,
        -0.2,
        -deg2rad(45.0),
        -deg2rad(45.0),
        0.1,
        -0.3,
        -1.0,
        0.1,
        -0.3,
        -1.0,
        0.0,
    ]
    const upper_bounds = [
        0.5,
        0.2,
        deg2rad(45.0),
        deg2rad(45.0),
        0.9,
        0.3,
        0.0,
        0.9,
        0.3,
        0.0,
        0.05,
    ]

    """
    ```
        function parsec2xty(
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

    Convert PARSEC parameters to `(x, t, y)` coordinates (`y` standing for camber @ chordwise control points)
    """
    function parsec2xty(
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

        cupp = pcoeff(
            leading_edge_radius, x_trailing_edge, y_trailing_edge + trailing_edge_thickness / 2,
            x_upper, y_upper, upper_curvature, -β / 2 - α
        )
        clow = pcoeff(
            leading_edge_radius, x_trailing_edge, - y_trailing_edge + trailing_edge_thickness / 2,
            x_lower, y_lower, lower_curvature, -β / 2 + α; intra = true,
        )

        yupp = polyval(
            x, cupp
        )
        ylow = polyval(
            x, clow
        )

        y = @. (yupp + ylow) / 2
        t = (yupp - ylow)

        (x, t, y)

    end

    """
    ```
        function xty2selig(
            x::AbstractVector, t::AbstractVector, y::AbstractVector
        )
    ```

    Convert airfoil from `(x, t, y)` format to Selig format
    """
    function xty2selig(
        x::AbstractVector, t::AbstractVector, y::AbstractVector
    )

        flipcat = (e, i) -> [
            e[end:-1:2];
            i
        ]

        [
            flipcat(x, x) flipcat(y .+ t ./ 2, y .- t ./ 2)
        ]

    end

    """
    ```
        function selig2xty(
            xs::AbstractVector, pts::AbstractMatrix
        )
    ```

    Convert Selig-format airfoil to xty format (using a pre-defined chordwise point grid `xs`)
    """
    function selig2xty(
        xs::AbstractVector, pts::AbstractMatrix
    )

        x, y = eachcol(pts)

        x = x .- minimum(x)

        init_dx = maximum(x)

        x = x .* (1.0 / init_dx) 

        ile = argmin(x)

        xupp, yupp = x[ile:-1:1], y[ile:-1:1]
        xlow, ylow = x[ile:end], y[ile:end]

        upper = linear_interpolation(
            xupp, yupp
        )
        lower = linear_interpolation(
            xlow, ylow
        )

        yupp = upper.(xs)
        ylow = lower.(xs)

        t = yupp .- ylow
        y = (yupp .+ ylow) ./ 2

        (xs, t, y)

    end

    """
    ```
        parsec2selig(xs::AbstractVector; kwargs...) = xty2selig(
            parsec2xty(xs; kwargs...)...
        )
    ```

    Convert PARSEC definitions to Selig format airfoil with pre-defined chordwise point grid `xs`
    """
    parsec2selig(xs::AbstractVector; kwargs...) = xty2selig(
        parsec2xty(xs; kwargs...)...
    )

    _trapz(y::AbstractVector, x::AbstractVector) = sum(
        (y[1:(end - 1)] .+ y[2:end]) .* (x[2:end] .- x[1:(end - 1)]) ./ 2
    )

    """
    ```
        function xty2parsec(
            x::AbstractVector, t::AbstractVector, y::AbstractVector;
            kwargs...
        )
    ```

    Use the BFGS optimization algorithm to find the PARSEC parameters that best fit an airfoil in xty format.

    `kwargs` are passed to Optim.jl
    """
    function xty2parsec(
        x::AbstractVector, t::AbstractVector, y::AbstractVector;
        kwargs...
    )

        p2kwargs = p -> Dict(
            [
                prop => v for (prop, v) in zip(
                    parsec_parameters, p
                )
            ]...
        )

        penalty = p -> let (_, tt, yy) = parsec2xty(
            x; p2kwargs(p)...
        )

            _trapz(
                (tt .- t) .^ 2, x
            ) + _trapz(
                (yy .- y) .^ 2, x
            )

        end

        x0 = copy(opt_x0)

        res = optimize(
            penalty,
            lower_bounds, upper_bounds,
            x0,
            Fminbox(BFGS());
            kwargs...
        )

        if !Optim.converged(res)
            @warn "PARSEC airfoil parameter optimization did not converge"
        end

        p2kwargs(
            Optim.minimizer(res)
        )

    end

    """
    ```
        selig2parsec(
            x::AbstractVector, pts::AbstractMatrix;
            kwargs...
        ) = xty2parsec(
            selig2xty(x, pts)...;
            kwargs...
        )
    ```

    Use the BFGS optimization algorithm to find the PARSEC parameters that best fit an airfoil in Selig format.

    `kwargs` are passed to Optim.jl
    """
    selig2parsec(
        x::AbstractVector, pts::AbstractMatrix;
        kwargs...
    ) = xty2parsec(
        selig2xty(x, pts)...;
        kwargs...
    )

end

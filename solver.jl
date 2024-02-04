import SciMLBase
using Unitful: Quantity, ustrip
using DiffEqBase: InternalITP, IntervalNonlinearProblem, ReturnCode, nextfloat_tdir

function SciMLBase.solve(prob::IntervalNonlinearProblem{IP, Tuple{T, T}}, alg::InternalITP,
    args...;
    maxiters = 1000, kwargs...) where {IP, T<:Quantity}
    f = Base.Fix2(prob.f, prob.p)
    left, right = prob.tspan # a and b
    
    fl, fr = f(left), f(right)
    ϵ = T(eps(T))
    if iszero(fl)
        return SciMLBase.build_solution(prob, alg, left, fl;
            retcode = ReturnCode.ExactSolutionLeft, left = left,
            right = right)
    elseif iszero(fr)
        return SciMLBase.build_solution(prob, alg, right, fr;
            retcode = ReturnCode.ExactSolutionRight, left = left,
            right = right)
    end

    #defining variables/cache
    k1 = T(alg.k1)
    k2 = T(alg.k2)
    n0 = T(alg.n0)
    n_h = ceil(log2(abs(right - left) / (2 * ϵ)))
    mid = (left + right) / 2
    x_f = (fr * left - fl * right) / (fr - fl)
    xt = left
    xp = left
    r = zero(left) #minmax radius
    δ = zero(left) # truncation error
    σ = 1.0
    ϵ_s = ϵ * 2^(ustrip(n_h) + ustrip(n0))
    i = 0 #iteration
    while i <= maxiters
        #mid = (left + right) / 2
        span = abs(right - left)
        r = ϵ_s - (span / 2)
        δ = k1 * (ustrip(span)^ustrip(k2))

        ## Interpolation step ##
        x_f = left + (right - left) * (fl / (fl - fr))

        ## Truncation step ##
        σ = sign(mid - x_f)
        if δ <= abs(mid - x_f)
            xt = x_f + (σ * δ)
        else
            xt = mid
        end

        ## Projection step ##
        if abs(xt - mid) <= r
            xp = xt
        else
            xp = mid - (σ * r)
        end

        ## Update ##
        tmin, tmax = minmax(left, right)
        xp >= tmax && (xp = prevfloat(tmax))
        xp <= tmin && (xp = nextfloat(tmin))
        yp = f(xp)
        yps = yp * sign(fr)
        if yps > 0
            right = xp
            fr = yp
        elseif yps < 0
            left = xp
            fl = yp
        else
            left = prevfloat_tdir(xp, prob.tspan...)
            right = xp
            return SciMLBase.build_solution(prob, alg, left, f(left);
                retcode = ReturnCode.Success, left = left,
                right = right)
        end
        i += 1
        mid = (left + right) / 2
        ϵ_s /= 2

        if nextfloat_tdir(left, prob.tspan...) == right
            return SciMLBase.build_solution(prob, alg, left, fl;
                retcode = ReturnCode.FloatingPointLimit, left = left,
                right = right)
        end
    end
    return SciMLBase.build_solution(prob, alg, left, fl; retcode = ReturnCode.MaxIters,
        left = left, right = right)
end
include("Graph.jl")
using LinearAlgebra
import Base: *

struct GaussianState{T<:Number}
    B::Matrix{T}
    Gsize::Vector{Int64}
    Gpos::Vector{Vector{Int64}}
end

struct Evo{T<:AbstractMatrix}
    M::T
end

function GaussianState(G::LatGraph,pos::Vector{Vector})
    # L is dim of single particles Hilbert space. N is num of filling states.
    L = length(G.pos);
    N = length(pos);
    B = zeros(ComplexF64,L,N);
    for i in 1:N
        j = findfirst(x->x==pos[i],G.pos)
        B[j,i] = 1
    end
    GaussianState(B,G.size,G.pos)
end

function Evo(
    H::AbstractMatrix,
    dt::Real
)
    M = exp(im*dt*H);
    Evo(M)
end

function density(G::GaussianState)
    corr = conj(G.B)*transpose(G.B);
    dv = [Real(corr[i,i]) for i in 1:length(G.Gpos)];
    den = zeros(G.Gsize...)
    for i in 1:length(G.Gpos)
        den[G.Gpos[i]...] = dv[i]
    end
    den
end

function entropy(s::GaussianState, i::AbstractVector{<:Integer})
    B = eltype(s) <: Real ? Float64.(s.B[i,:]) : ComplexF64.(s.B[i,:])
    vals = svdvals(B).^2
    EE = 0.0
    for x in vals
        if x < 1
            x < 1e-14 && continue
            EE -= x*log(x)+(1-x)*log(1-x)
        else
            @assert x-1 < 1e-6 "Got a Schmidt value λ = $x."
        end
    end
    EE
end

function entropy(f,s::GaussianState)
    i = findall(f,s.Gpos);
    entropy(s,i)
end


function expm(A::AbstractMatrix, order::Integer=30)
    mat = I + A / order
    order -= 1
    while order > 0
        mat = A * mat
        mat ./= order
        mat += I
        order -= 1
    end
    mat
end

function *(E::Evo, G::GaussianState)
    M_new = E.M*G.B;
    GaussianState(Matrix(qr(M_new).Q),G.Gsize,G.Gpos)
end
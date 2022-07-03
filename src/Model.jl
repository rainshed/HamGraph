include("GraphGuass.jl")
using DelimitedFiles,Plots

function Model1(
    t1::Real,
    t2::Real,
    γ::Real,
    N::Integer
    )
    """
    H(k) = 
    
    0        t1+γ*sin(ky)+t2*exp(ikx)
    t1-γ*sin(ky)+t2*exp(-ikx)       0

    OBC on x direction, PBC on y direction
    """ 
    G = LatGraph(x->true,[N,N,2])
    # onsite t hopping
    #AddEdge!(G,[0,0,1],t1,BC="PBC")
    AddEdge!(G,[0,0,-1],t1,BC="PBC")
    # x direction hopping
    AddEdge!(G,[1,0,-1],t2,BC="OBC")
    AddEdge!(G,[-1,0,1],t2,BC="OBC")
    # y direction hopping
    AddEdge!(G,[0,1,1],im*γ,BC="OBC")
    AddEdge!(G,[0,-1,-1],im*γ,BC="OBC")
    AddEdge!(G,[0,1,-1],-im*γ,BC="OBC")
    AddEdge!(G,[0,-1,1],-im*γ,BC="OBC")
    G
end



function initPos(N::Int,n::Int)
    pos = Vector[];
    for i in 1:n
        for j in 1:n
            x = N÷2-n÷2+i
            y = N÷2-n÷2+j
            push!(pos,[x,y,1])
            push!(pos,[x,y,2])
        end
    end
    pos
end

function test()
    N=5;
    n=2;
    G = Model1(1,1,0.2,N);
    H = ConMat(G);
    evo = Evo(H,0.5);
    pos = initPos(N,n);
    stat = GaussianState(G,pos);
    T = 100;
    den = Vector(undef,T);
    ent = Vector(undef,T);
    for i in 1:T
        stat = evo*stat;
        den[i] = density(stat);
        ent[i] = entropy(x->x[1]<N÷2,stat);
    end
    
    d = zeros(T,N,N);
    for t in 1:T
        for i in 1:N
            for j in 1:N
                d[t,j,i] = den[t][i,j,1]+den[t][i,j,2];
            end
        end
    end
    ent
    
end


#------------------------------------------------------------------

function Model2(L::Int,shape::String)
    """
        H = cos(kx)+i*sin(ky)
    """
    if shape=="Square"
        G = SquareGraph(L)
    elseif shape=="LX"
        G = LatGraph([[L÷2+1,1],[1,L÷2+1],[L÷2+1,L],[L,L÷2+1]],[L,L])
    end

    AddEdge!(G,[1,0],1);
    AddEdge!(G,[-1,0],1);
    AddEdge!(G,[0,1],1);
    AddEdge!(G,[0,-1],-1);

    G
end

function initState2(G::LatGraph,n::Int64)
    N = G.size[1];
    pos = Vector[];
    for i in 1:n
        for j in 1:n
            x = N÷2-n÷2+i
            y = N÷2-n÷2+j
            push!(pos,[x,y])
        end
    end
    GaussianState(G,pos)
end

function main2(L,n,T)
    shape = "Square";

    G = Model2(L,shape);
    state = initState2(G,n);
    H = ConMat(G);
    evo = Evo(H,0.5);


    den = Vector(undef,T+1);
    den[1] = density(state)
    #ent = Vector(undef,T);
    for i in 1:T
        state = evo*state;
        den[i+1] = density(state);
        #ent[i] = entropy(x->x[1]<N÷2,stat);
    end
    den
end

L = 51
n = 10
T = 10
den = main2(L,n,T);
@gif for i in 1:T
    heatmap(1:L,1:L,den[i])
end

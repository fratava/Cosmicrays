include("JF2012_simple.jl")
include("RG4_D.jl")
include("Propa_JF2012_simple.jl")

using MPI
using SkyCoords
#Pkg.add("SkyCoords")


MPI.Init()
comm = MPI.COMM_WORLD
const rank = MPI.Comm_rank(comm)
const size = MPI.Comm_size(comm)
for k = 1:size
    if  k-1==rank
        print("Hello World, I am process $rank of $size running on name. \n" )
    end
end
#const name = MPI.Comm_name(comm)
############################### Constantes
const Coord_EQ = 0 # 1 para mapa en EQUATORIAL; 0 para GALACTIC
const granularity=1
const Charge=1 # 1 para p, 8 para O, 26 para Fe
const Mass =1 # 1 para p, 16 para O, 56 para Fe
const Energy =60
############################### Constantes
const N_lon = 360
const N_lat = 180
Desviacion=zeros(N_lat*granularity,N_lon*granularity)
 
   
for i = 1:N_lat*granularity
    if  (i-1)%size == rank
        print(" i=$i ")
        for j = 1:N_lon*granularity
            if Coord_EQ == 1
                eq_dec=90.-(i-.5)/granularity
                eq_ra= (j-.5)/granularity
                c1 = FK5Coords{2000}(deg2rad(eq_ra), deg2rad(eq_dec))  # inputs are ra, dec in radians
                c2 = convert(GalCoords, c1)
                galactic_l=rad2deg(c2.l)
                galactic_b=rad2deg(c2.b)
                alpha=Propa_JF2012_simple(Charge,Mass,Energy, galactic_l, galactic_b,0);    
                Desviacion[i,j]=alpha
            else
                galactic_b=90.-(i-.5)/granularity
                galactic_l= (j-.5)/granularity-180.
                alpha=round(float(Propa_JF2012_simple(Charge,Mass,Energy, galactic_l, galactic_b,0)),12);
                Desviacion[i,j]=alpha
            end
        end
    
    end    
end
data2=MPI.Reduce(Desviacion, MPI.SUM, 0, comm)

if rank==0
    println("Granularity= $granularity  N_lon=$N_lon N_lat=$N_lat \n")
    data21 = reshape (data2, N_lat*granularity,N_lon*granularity)
    if Coord_EQ == 1
        writedlm("Desviacion_EQ_p60g5.dat", data21) 
    else
        writedlm("Desviacion_GAL_p60g5.dat", data21) 

    end
end

MPI.Finalize()  

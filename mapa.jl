#Using libraries
using PyPlot

#Include some code
include("JF2012_simple.jl")
include("RG4_D.jl")
include("Propa_JF2012_simple.jl")

#Features constantes
Charge=26 # 1 for p, 8 for O, 26 for Fe
Mass=56 # 1 for p, 16 for O, 56 for Fe
Energy=60 #Eev
N_lon=360
N_lat=180
granularity=1
latitude=N_lat*granularity
longitude=N_lon*granularity
Deviation=zeros(latitude,longitude)
print("Theta=")

#Propagation
for i= 1:latitude
    if i%granularity==0
        print(i/granularity,"-")
    end
    
    for j= 1:longitude
        lat=90-(i-0.5)/granularity
        lon=(j-0.5)/granularity-180.0
alpha = round(float(Propa_JF2012_simple(Charge,Mass,Energy,lon,lat,0)),12);
        
        Deviation[i,j]=alpha
    end
end


#Plot
PyPlot.figure(figsize=(12,7))
PyPlot.subplot(111,projection="hammer")
x=linspace(-pi,pi,longitude)
y=linspace(pi/2,-pi/2,latitude)
PyPlot.pcolormesh(x,y,Deviation,cmap="jet",shading="interp")
PyPlot.colorbar()
PyPlot.grid(true)
PyPlot.show()

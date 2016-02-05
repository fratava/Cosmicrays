function RG4_D(f,XX0)
    N=50; # Cambiar este no. de iteraciones si se requiere
    c=2.998*10^(8.); # Usamos MKS
    kpc=3.0857*10^19.;
    kyear=3.154*10^10.;
    t=0;niter=0;
    h=20*kpc/c/N/2.;
    #Esta es la función Runge-Kutta con sus coeficientes
    D=0;
    while (D<20.0*kpc )
        D=sqrt(XX0[1]*XX0[1]+XX0[2]*XX0[2]+XX0[3]*XX0[3]);
        k1=f(XX0,t)
        k2=f(XX0+(0.5*k1*h),t+.5*h)
        k3=f(XX0+(0.5*k2*h),t+.5*h)
        k4=f(XX0+(k3*h),t+h)
        XX0=XX0 + ((h/6) * (k1 + 2 * (k2 + k3) + k4)) 
        t += h;niter += 1;
    end 
    return XX0,t,niter #Regresa un arreglo donde están 
end

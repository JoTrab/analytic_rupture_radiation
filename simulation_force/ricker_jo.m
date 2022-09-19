function s=ricker_jo(f, t)


ricker = (1-2*pi^2*f.^2.*t.^2).*exp(-pi^2*f^2.*t.^2)

% (f, length=0.512, dt=0.001):
%     t = np.linspace(-length/2, (length-dt)/2, length/dt)
%     y = (1.-2.*(np.pi**2)*(f**2)*(t**2))*np.exp(-(np.pi**2)*(f**2)*(t**2))

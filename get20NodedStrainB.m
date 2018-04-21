function strain_mat = get20NodedStrainB(in1,in2,zeta,eta,nu)
%GET20NODEDSTRAINB
%    STRAIN_MAT = GET20NODEDSTRAINB(IN1,IN2,ZETA,ETA,NU)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    11-Apr-2018 16:38:23

height = in2(:,3);
thickness = in2(:,1);
width = in2(:,2);
strain_mat = reshape([(zeta.*5.0e-1-eta.*nu.*2.5e-1-eta.*zeta.*5.0e-1-nu.*zeta.*5.0e-1+eta.*nu.^2.*2.5e-1+eta.^2.*nu.*2.5e-1-eta.^2.*2.5e-1-nu.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1+2.5e-1)./thickness,0.0,0.0,0.0,(nu.*5.0e-1-eta.*nu.*5.0e-1-eta.*zeta.*2.5e-1-nu.*zeta.*5.0e-1+eta.*zeta.^2.*2.5e-1+eta.^2.*zeta.*2.5e-1-eta.^2.*2.5e-1-zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1+2.5e-1)./height,(eta.*5.0e-1-eta.*nu.*5.0e-1-eta.*zeta.*5.0e-1-nu.*zeta.*2.5e-1+nu.*zeta.^2.*2.5e-1+nu.^2.*zeta.*2.5e-1-nu.^2.*2.5e-1-zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1+2.5e-1)./width,0.0,(eta.*5.0e-1-eta.*nu.*5.0e-1-eta.*zeta.*5.0e-1-nu.*zeta.*2.5e-1+nu.*zeta.^2.*2.5e-1+nu.^2.*zeta.*2.5e-1-nu.^2.*2.5e-1-zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1+2.5e-1)./width,0.0,(nu.*5.0e-1-eta.*nu.*5.0e-1-eta.*zeta.*2.5e-1-nu.*zeta.*5.0e-1+eta.*zeta.^2.*2.5e-1+eta.^2.*zeta.*2.5e-1-eta.^2.*2.5e-1-zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1+2.5e-1)./height,0.0,(zeta.*5.0e-1-eta.*nu.*2.5e-1-eta.*zeta.*5.0e-1-nu.*zeta.*5.0e-1+eta.*nu.^2.*2.5e-1+eta.^2.*nu.*2.5e-1-eta.^2.*2.5e-1-nu.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1+2.5e-1)./thickness,0.0,0.0,(nu.*5.0e-1-eta.*nu.*5.0e-1-eta.*zeta.*2.5e-1-nu.*zeta.*5.0e-1+eta.*zeta.^2.*2.5e-1+eta.^2.*zeta.*2.5e-1-eta.^2.*2.5e-1-zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1+2.5e-1)./height,(eta.*5.0e-1-eta.*nu.*5.0e-1-eta.*zeta.*5.0e-1-nu.*zeta.*2.5e-1+nu.*zeta.^2.*2.5e-1+nu.^2.*zeta.*2.5e-1-nu.^2.*2.5e-1-zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1+2.5e-1)./width,(zeta.*5.0e-1-eta.*nu.*2.5e-1-eta.*zeta.*5.0e-1-nu.*zeta.*5.0e-1+eta.*nu.^2.*2.5e-1+eta.^2.*nu.*2.5e-1-eta.^2.*2.5e-1-nu.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1+2.5e-1)./thickness,0.0,(zeta.*5.0e-1+eta.*nu.*2.5e-1-eta.*zeta.*5.0e-1-nu.*zeta.*5.0e-1-eta.*nu.^2.*2.5e-1-eta.^2.*nu.*2.5e-1+eta.^2.*2.5e-1+nu.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1-2.5e-1)./thickness,0.0,0.0,0.0,(nu.*5.0e-1-eta.*nu.*5.0e-1+eta.*zeta.*2.5e-1+nu.*zeta.*5.0e-1+eta.*zeta.^2.*2.5e-1-eta.^2.*zeta.*2.5e-1-eta.^2.*2.5e-1-zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1+2.5e-1)./height,(eta.*5.0e-1-eta.*nu.*5.0e-1+eta.*zeta.*5.0e-1+nu.*zeta.*2.5e-1+nu.*zeta.^2.*2.5e-1-nu.^2.*zeta.*2.5e-1-nu.^2.*2.5e-1-zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1+2.5e-1)./width,0.0,(eta.*5.0e-1-eta.*nu.*5.0e-1+eta.*zeta.*5.0e-1+nu.*zeta.*2.5e-1+nu.*zeta.^2.*2.5e-1-nu.^2.*zeta.*2.5e-1-nu.^2.*2.5e-1-zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1+2.5e-1)./width,0.0,(nu.*5.0e-1-eta.*nu.*5.0e-1+eta.*zeta.*2.5e-1+nu.*zeta.*5.0e-1+eta.*zeta.^2.*2.5e-1-eta.^2.*zeta.*2.5e-1-eta.^2.*2.5e-1-zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1+2.5e-1)./height,0.0,(zeta.*5.0e-1+eta.*nu.*2.5e-1-eta.*zeta.*5.0e-1-nu.*zeta.*5.0e-1-eta.*nu.^2.*2.5e-1-eta.^2.*nu.*2.5e-1+eta.^2.*2.5e-1+nu.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1-2.5e-1)./thickness,0.0,0.0,(nu.*5.0e-1-eta.*nu.*5.0e-1+eta.*zeta.*2.5e-1+nu.*zeta.*5.0e-1+eta.*zeta.^2.*2.5e-1-eta.^2.*zeta.*2.5e-1-eta.^2.*2.5e-1-zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1+2.5e-1)./height,(eta.*5.0e-1-eta.*nu.*5.0e-1+eta.*zeta.*5.0e-1+nu.*zeta.*2.5e-1+nu.*zeta.^2.*2.5e-1-nu.^2.*zeta.*2.5e-1-nu.^2.*2.5e-1-zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1+2.5e-1)./width,(zeta.*5.0e-1+eta.*nu.*2.5e-1-eta.*zeta.*5.0e-1-nu.*zeta.*5.0e-1-eta.*nu.^2.*2.5e-1-eta.^2.*nu.*2.5e-1+eta.^2.*2.5e-1+nu.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1-2.5e-1)./thickness,0.0,(zeta.*5.0e-1-eta.*nu.*2.5e-1+eta.*zeta.*5.0e-1-nu.*zeta.*5.0e-1+eta.*nu.^2.*2.5e-1-eta.^2.*nu.*2.5e-1+eta.^2.*2.5e-1+nu.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1-2.5e-1)./thickness,0.0,0.0,0.0,(nu.*5.0e-1+eta.*nu.*5.0e-1-eta.*zeta.*2.5e-1+nu.*zeta.*5.0e-1-eta.*zeta.^2.*2.5e-1-eta.^2.*zeta.*2.5e-1-eta.^2.*2.5e-1-zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1+2.5e-1)./height,(eta.*5.0e-1-eta.*nu.*5.0e-1+eta.*zeta.*5.0e-1-nu.*zeta.*2.5e-1-nu.*zeta.^2.*2.5e-1+nu.^2.*zeta.*2.5e-1+nu.^2.*2.5e-1+zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1-2.5e-1)./width,0.0,(eta.*5.0e-1-eta.*nu.*5.0e-1+eta.*zeta.*5.0e-1-nu.*zeta.*2.5e-1-nu.*zeta.^2.*2.5e-1+nu.^2.*zeta.*2.5e-1+nu.^2.*2.5e-1+zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1-2.5e-1)./width,0.0,(nu.*5.0e-1+eta.*nu.*5.0e-1-eta.*zeta.*2.5e-1+nu.*zeta.*5.0e-1-eta.*zeta.^2.*2.5e-1-eta.^2.*zeta.*2.5e-1-eta.^2.*2.5e-1-zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1+2.5e-1)./height,0.0,(zeta.*5.0e-1-eta.*nu.*2.5e-1+eta.*zeta.*5.0e-1-nu.*zeta.*5.0e-1+eta.*nu.^2.*2.5e-1-eta.^2.*nu.*2.5e-1+eta.^2.*2.5e-1+nu.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1-2.5e-1)./thickness,0.0,0.0,(nu.*5.0e-1+eta.*nu.*5.0e-1-eta.*zeta.*2.5e-1+nu.*zeta.*5.0e-1-eta.*zeta.^2.*2.5e-1-eta.^2.*zeta.*2.5e-1-eta.^2.*2.5e-1-zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1+2.5e-1)./height,(eta.*5.0e-1-eta.*nu.*5.0e-1+eta.*zeta.*5.0e-1-nu.*zeta.*2.5e-1-nu.*zeta.^2.*2.5e-1+nu.^2.*zeta.*2.5e-1+nu.^2.*2.5e-1+zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1-2.5e-1)./width,(zeta.*5.0e-1-eta.*nu.*2.5e-1+eta.*zeta.*5.0e-1-nu.*zeta.*5.0e-1+eta.*nu.^2.*2.5e-1-eta.^2.*nu.*2.5e-1+eta.^2.*2.5e-1+nu.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1-2.5e-1)./thickness,0.0,(zeta.*5.0e-1+eta.*nu.*2.5e-1+eta.*zeta.*5.0e-1-nu.*zeta.*5.0e-1-eta.*nu.^2.*2.5e-1+eta.^2.*nu.*2.5e-1-eta.^2.*2.5e-1-nu.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1+2.5e-1)./thickness,0.0,0.0,0.0,(nu.*5.0e-1+eta.*nu.*5.0e-1+eta.*zeta.*2.5e-1-nu.*zeta.*5.0e-1-eta.*zeta.^2.*2.5e-1+eta.^2.*zeta.*2.5e-1-eta.^2.*2.5e-1-zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1+2.5e-1)./height,(eta.*5.0e-1-eta.*nu.*5.0e-1-eta.*zeta.*5.0e-1+nu.*zeta.*2.5e-1-nu.*zeta.^2.*2.5e-1-nu.^2.*zeta.*2.5e-1+nu.^2.*2.5e-1+zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1-2.5e-1)./width,0.0,(eta.*5.0e-1-eta.*nu.*5.0e-1-eta.*zeta.*5.0e-1+nu.*zeta.*2.5e-1-nu.*zeta.^2.*2.5e-1-nu.^2.*zeta.*2.5e-1+nu.^2.*2.5e-1+zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1-2.5e-1)./width,0.0,(nu.*5.0e-1+eta.*nu.*5.0e-1+eta.*zeta.*2.5e-1-nu.*zeta.*5.0e-1-eta.*zeta.^2.*2.5e-1+eta.^2.*zeta.*2.5e-1-eta.^2.*2.5e-1-zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1+2.5e-1)./height,0.0,(zeta.*5.0e-1+eta.*nu.*2.5e-1+eta.*zeta.*5.0e-1-nu.*zeta.*5.0e-1-eta.*nu.^2.*2.5e-1+eta.^2.*nu.*2.5e-1-eta.^2.*2.5e-1-nu.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1+2.5e-1)./thickness,0.0,0.0,(nu.*5.0e-1+eta.*nu.*5.0e-1+eta.*zeta.*2.5e-1-nu.*zeta.*5.0e-1-eta.*zeta.^2.*2.5e-1+eta.^2.*zeta.*2.5e-1-eta.^2.*2.5e-1-zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1+2.5e-1)./height,(eta.*5.0e-1-eta.*nu.*5.0e-1-eta.*zeta.*5.0e-1+nu.*zeta.*2.5e-1-nu.*zeta.^2.*2.5e-1-nu.^2.*zeta.*2.5e-1+nu.^2.*2.5e-1+zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1-2.5e-1)./width,(zeta.*5.0e-1+eta.*nu.*2.5e-1+eta.*zeta.*5.0e-1-nu.*zeta.*5.0e-1-eta.*nu.^2.*2.5e-1+eta.^2.*nu.*2.5e-1-eta.^2.*2.5e-1-nu.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1+2.5e-1)./thickness,0.0,(zeta.*5.0e-1+eta.*nu.*2.5e-1-eta.*zeta.*5.0e-1+nu.*zeta.*5.0e-1+eta.*nu.^2.*2.5e-1-eta.^2.*nu.*2.5e-1-eta.^2.*2.5e-1-nu.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1+2.5e-1)./thickness,0.0,0.0,0.0,(nu.*5.0e-1-eta.*nu.*5.0e-1+eta.*zeta.*2.5e-1-nu.*zeta.*5.0e-1-eta.*zeta.^2.*2.5e-1-eta.^2.*zeta.*2.5e-1+eta.^2.*2.5e-1+zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1-2.5e-1)./height,(eta.*5.0e-1+eta.*nu.*5.0e-1-eta.*zeta.*5.0e-1+nu.*zeta.*2.5e-1-nu.*zeta.^2.*2.5e-1+nu.^2.*zeta.*2.5e-1-nu.^2.*2.5e-1-zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1+2.5e-1)./width,0.0,(eta.*5.0e-1+eta.*nu.*5.0e-1-eta.*zeta.*5.0e-1+nu.*zeta.*2.5e-1-nu.*zeta.^2.*2.5e-1+nu.^2.*zeta.*2.5e-1-nu.^2.*2.5e-1-zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1+2.5e-1)./width,0.0,(nu.*5.0e-1-eta.*nu.*5.0e-1+eta.*zeta.*2.5e-1-nu.*zeta.*5.0e-1-eta.*zeta.^2.*2.5e-1-eta.^2.*zeta.*2.5e-1+eta.^2.*2.5e-1+zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1-2.5e-1)./height,0.0,(zeta.*5.0e-1+eta.*nu.*2.5e-1-eta.*zeta.*5.0e-1+nu.*zeta.*5.0e-1+eta.*nu.^2.*2.5e-1-eta.^2.*nu.*2.5e-1-eta.^2.*2.5e-1-nu.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1+2.5e-1)./thickness,0.0,0.0,(nu.*5.0e-1-eta.*nu.*5.0e-1+eta.*zeta.*2.5e-1-nu.*zeta.*5.0e-1-eta.*zeta.^2.*2.5e-1-eta.^2.*zeta.*2.5e-1+eta.^2.*2.5e-1+zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1-2.5e-1)./height,(eta.*5.0e-1+eta.*nu.*5.0e-1-eta.*zeta.*5.0e-1+nu.*zeta.*2.5e-1-nu.*zeta.^2.*2.5e-1+nu.^2.*zeta.*2.5e-1-nu.^2.*2.5e-1-zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1+2.5e-1)./width,(zeta.*5.0e-1+eta.*nu.*2.5e-1-eta.*zeta.*5.0e-1+nu.*zeta.*5.0e-1+eta.*nu.^2.*2.5e-1-eta.^2.*nu.*2.5e-1-eta.^2.*2.5e-1-nu.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1+2.5e-1)./thickness,0.0,(zeta.*5.0e-1-eta.*nu.*2.5e-1-eta.*zeta.*5.0e-1+nu.*zeta.*5.0e-1-eta.*nu.^2.*2.5e-1+eta.^2.*nu.*2.5e-1+eta.^2.*2.5e-1+nu.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1-2.5e-1)./thickness,0.0,0.0,0.0,(nu.*5.0e-1-eta.*nu.*5.0e-1-eta.*zeta.*2.5e-1+nu.*zeta.*5.0e-1-eta.*zeta.^2.*2.5e-1+eta.^2.*zeta.*2.5e-1+eta.^2.*2.5e-1+zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1-2.5e-1)./height,(eta.*5.0e-1+eta.*nu.*5.0e-1+eta.*zeta.*5.0e-1-nu.*zeta.*2.5e-1-nu.*zeta.^2.*2.5e-1-nu.^2.*zeta.*2.5e-1-nu.^2.*2.5e-1-zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1+2.5e-1)./width,0.0,(eta.*5.0e-1+eta.*nu.*5.0e-1+eta.*zeta.*5.0e-1-nu.*zeta.*2.5e-1-nu.*zeta.^2.*2.5e-1-nu.^2.*zeta.*2.5e-1-nu.^2.*2.5e-1-zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1+2.5e-1)./width,0.0,(nu.*5.0e-1-eta.*nu.*5.0e-1-eta.*zeta.*2.5e-1+nu.*zeta.*5.0e-1-eta.*zeta.^2.*2.5e-1+eta.^2.*zeta.*2.5e-1+eta.^2.*2.5e-1+zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1-2.5e-1)./height,0.0,(zeta.*5.0e-1-eta.*nu.*2.5e-1-eta.*zeta.*5.0e-1+nu.*zeta.*5.0e-1-eta.*nu.^2.*2.5e-1+eta.^2.*nu.*2.5e-1+eta.^2.*2.5e-1+nu.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1-2.5e-1)./thickness,0.0,0.0,(nu.*5.0e-1-eta.*nu.*5.0e-1-eta.*zeta.*2.5e-1+nu.*zeta.*5.0e-1-eta.*zeta.^2.*2.5e-1+eta.^2.*zeta.*2.5e-1+eta.^2.*2.5e-1+zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1-2.5e-1)./height,(eta.*5.0e-1+eta.*nu.*5.0e-1+eta.*zeta.*5.0e-1-nu.*zeta.*2.5e-1-nu.*zeta.^2.*2.5e-1-nu.^2.*zeta.*2.5e-1-nu.^2.*2.5e-1-zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1+2.5e-1)./width,(zeta.*5.0e-1-eta.*nu.*2.5e-1-eta.*zeta.*5.0e-1+nu.*zeta.*5.0e-1-eta.*nu.^2.*2.5e-1+eta.^2.*nu.*2.5e-1+eta.^2.*2.5e-1+nu.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1-2.5e-1)./thickness,0.0,(zeta.*5.0e-1+eta.*nu.*2.5e-1+eta.*zeta.*5.0e-1+nu.*zeta.*5.0e-1+eta.*nu.^2.*2.5e-1+eta.^2.*nu.*2.5e-1+eta.^2.*2.5e-1+nu.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1-2.5e-1)./thickness,0.0,0.0,0.0,(nu.*5.0e-1+eta.*nu.*5.0e-1+eta.*zeta.*2.5e-1+nu.*zeta.*5.0e-1+eta.*zeta.^2.*2.5e-1+eta.^2.*zeta.*2.5e-1+eta.^2.*2.5e-1+zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1-2.5e-1)./height,(eta.*5.0e-1+eta.*nu.*5.0e-1+eta.*zeta.*5.0e-1+nu.*zeta.*2.5e-1+nu.*zeta.^2.*2.5e-1+nu.^2.*zeta.*2.5e-1+nu.^2.*2.5e-1+zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1-2.5e-1)./width,0.0,(eta.*5.0e-1+eta.*nu.*5.0e-1+eta.*zeta.*5.0e-1+nu.*zeta.*2.5e-1+nu.*zeta.^2.*2.5e-1+nu.^2.*zeta.*2.5e-1+nu.^2.*2.5e-1+zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1-2.5e-1)./width,0.0,(nu.*5.0e-1+eta.*nu.*5.0e-1+eta.*zeta.*2.5e-1+nu.*zeta.*5.0e-1+eta.*zeta.^2.*2.5e-1+eta.^2.*zeta.*2.5e-1+eta.^2.*2.5e-1+zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1-2.5e-1)./height,0.0,(zeta.*5.0e-1+eta.*nu.*2.5e-1+eta.*zeta.*5.0e-1+nu.*zeta.*5.0e-1+eta.*nu.^2.*2.5e-1+eta.^2.*nu.*2.5e-1+eta.^2.*2.5e-1+nu.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1-2.5e-1)./thickness,0.0,0.0,(nu.*5.0e-1+eta.*nu.*5.0e-1+eta.*zeta.*2.5e-1+nu.*zeta.*5.0e-1+eta.*zeta.^2.*2.5e-1+eta.^2.*zeta.*2.5e-1+eta.^2.*2.5e-1+zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1-2.5e-1)./height,(eta.*5.0e-1+eta.*nu.*5.0e-1+eta.*zeta.*5.0e-1+nu.*zeta.*2.5e-1+nu.*zeta.^2.*2.5e-1+nu.^2.*zeta.*2.5e-1+nu.^2.*2.5e-1+zeta.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1-2.5e-1)./width,(zeta.*5.0e-1+eta.*nu.*2.5e-1+eta.*zeta.*5.0e-1+nu.*zeta.*5.0e-1+eta.*nu.^2.*2.5e-1+eta.^2.*nu.*2.5e-1+eta.^2.*2.5e-1+nu.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1-2.5e-1)./thickness,0.0,(zeta.*5.0e-1-eta.*nu.*2.5e-1+eta.*zeta.*5.0e-1+nu.*zeta.*5.0e-1-eta.*nu.^2.*2.5e-1-eta.^2.*nu.*2.5e-1-eta.^2.*2.5e-1-nu.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1+2.5e-1)./thickness,0.0,0.0,0.0,(nu.*5.0e-1+eta.*nu.*5.0e-1-eta.*zeta.*2.5e-1-nu.*zeta.*5.0e-1+eta.*zeta.^2.*2.5e-1-eta.^2.*zeta.*2.5e-1+eta.^2.*2.5e-1+zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1-2.5e-1)./height,(eta.*5.0e-1+eta.*nu.*5.0e-1-eta.*zeta.*5.0e-1-nu.*zeta.*2.5e-1+nu.*zeta.^2.*2.5e-1-nu.^2.*zeta.*2.5e-1+nu.^2.*2.5e-1+zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1-2.5e-1)./width,0.0,(eta.*5.0e-1+eta.*nu.*5.0e-1-eta.*zeta.*5.0e-1-nu.*zeta.*2.5e-1+nu.*zeta.^2.*2.5e-1-nu.^2.*zeta.*2.5e-1+nu.^2.*2.5e-1+zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1-2.5e-1)./width,0.0,(nu.*5.0e-1+eta.*nu.*5.0e-1-eta.*zeta.*2.5e-1-nu.*zeta.*5.0e-1+eta.*zeta.^2.*2.5e-1-eta.^2.*zeta.*2.5e-1+eta.^2.*2.5e-1+zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1-2.5e-1)./height,0.0,(zeta.*5.0e-1-eta.*nu.*2.5e-1+eta.*zeta.*5.0e-1+nu.*zeta.*5.0e-1-eta.*nu.^2.*2.5e-1-eta.^2.*nu.*2.5e-1-eta.^2.*2.5e-1-nu.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1+2.5e-1)./thickness,0.0,0.0,(nu.*5.0e-1+eta.*nu.*5.0e-1-eta.*zeta.*2.5e-1-nu.*zeta.*5.0e-1+eta.*zeta.^2.*2.5e-1-eta.^2.*zeta.*2.5e-1+eta.^2.*2.5e-1+zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1-2.5e-1)./height,(eta.*5.0e-1+eta.*nu.*5.0e-1-eta.*zeta.*5.0e-1-nu.*zeta.*2.5e-1+nu.*zeta.^2.*2.5e-1-nu.^2.*zeta.*2.5e-1+nu.^2.*2.5e-1+zeta.^2.*2.5e-1-eta.*nu.*zeta.*5.0e-1-2.5e-1)./width,(zeta.*5.0e-1-eta.*nu.*2.5e-1+eta.*zeta.*5.0e-1+nu.*zeta.*5.0e-1-eta.*nu.^2.*2.5e-1-eta.^2.*nu.*2.5e-1-eta.^2.*2.5e-1-nu.^2.*2.5e-1+eta.*nu.*zeta.*5.0e-1+2.5e-1)./thickness,0.0,((nu-1.0).*(eta.^2-1.0).*5.0e-1)./thickness,0.0,0.0,0.0,(zeta.*-5.0e-1+eta.^2.*zeta.*5.0e-1+eta.^2.*5.0e-1-5.0e-1)./height,(eta.*(nu-1.0).*(zeta+1.0).*1.0)./width,0.0,(eta.*(nu-1.0).*(zeta+1.0).*1.0)./width,0.0,(zeta.*-5.0e-1+eta.^2.*zeta.*5.0e-1+eta.^2.*5.0e-1-5.0e-1)./height,0.0,((nu-1.0).*(eta.^2-1.0).*5.0e-1)./thickness,0.0,0.0,(zeta.*-5.0e-1+eta.^2.*zeta.*5.0e-1+eta.^2.*5.0e-1-5.0e-1)./height,(eta.*(nu-1.0).*(zeta+1.0).*1.0)./width,((nu-1.0).*(eta.^2-1.0).*5.0e-1)./thickness,0.0,(zeta.*(eta+1.0).*(nu-1.0).*1.0)./thickness,0.0,0.0,0.0,(eta.*-5.0e-1+eta.*zeta.^2.*5.0e-1+zeta.^2.*5.0e-1-5.0e-1)./height,((nu-1.0).*(zeta.^2-1.0).*5.0e-1)./width,0.0,((nu-1.0).*(zeta.^2-1.0).*5.0e-1)./width,0.0,(eta.*-5.0e-1+eta.*zeta.^2.*5.0e-1+zeta.^2.*5.0e-1-5.0e-1)./height,0.0,(zeta.*(eta+1.0).*(nu-1.0).*1.0)./thickness,0.0,0.0,(eta.*-5.0e-1+eta.*zeta.^2.*5.0e-1+zeta.^2.*5.0e-1-5.0e-1)./height,((nu-1.0).*(zeta.^2-1.0).*5.0e-1)./width,(zeta.*(eta+1.0).*(nu-1.0).*1.0)./thickness,0.0,((nu-1.0).*(eta.^2-1.0).*-5.0e-1)./thickness,0.0,0.0,0.0,(zeta.*5.0e-1-eta.^2.*zeta.*5.0e-1+eta.^2.*5.0e-1-5.0e-1)./height,(eta.*(nu-1.0).*(zeta-1.0).*-1.0)./width,0.0,(eta.*(nu-1.0).*(zeta-1.0).*-1.0)./width,0.0,(zeta.*5.0e-1-eta.^2.*zeta.*5.0e-1+eta.^2.*5.0e-1-5.0e-1)./height,0.0,((nu-1.0).*(eta.^2-1.0).*-5.0e-1)./thickness,0.0,0.0,(zeta.*5.0e-1-eta.^2.*zeta.*5.0e-1+eta.^2.*5.0e-1-5.0e-1)./height,(eta.*(nu-1.0).*(zeta-1.0).*-1.0)./width,((nu-1.0).*(eta.^2-1.0).*-5.0e-1)./thickness,0.0,(zeta.*(eta-1.0).*(nu-1.0).*-1.0)./thickness,0.0,0.0,0.0,(eta.*5.0e-1-eta.*zeta.^2.*5.0e-1+zeta.^2.*5.0e-1-5.0e-1)./height,((nu-1.0).*(zeta.^2-1.0).*-5.0e-1)./width,0.0,((nu-1.0).*(zeta.^2-1.0).*-5.0e-1)./width,0.0,(eta.*5.0e-1-eta.*zeta.^2.*5.0e-1+zeta.^2.*5.0e-1-5.0e-1)./height,0.0,(zeta.*(eta-1.0).*(nu-1.0).*-1.0)./thickness,0.0,0.0,(eta.*5.0e-1-eta.*zeta.^2.*5.0e-1+zeta.^2.*5.0e-1-5.0e-1)./height,((nu-1.0).*(zeta.^2-1.0).*-5.0e-1)./width,(zeta.*(eta-1.0).*(nu-1.0).*-1.0)./thickness,0.0,((nu+1.0).*(eta.^2-1.0).*-5.0e-1)./thickness,0.0,0.0,0.0,(zeta.*5.0e-1-eta.^2.*zeta.*5.0e-1-eta.^2.*5.0e-1+5.0e-1)./height,(eta.*(nu+1.0).*(zeta+1.0).*-1.0)./width,0.0,(eta.*(nu+1.0).*(zeta+1.0).*-1.0)./width,0.0,(zeta.*5.0e-1-eta.^2.*zeta.*5.0e-1-eta.^2.*5.0e-1+5.0e-1)./height,0.0,((nu+1.0).*(eta.^2-1.0).*-5.0e-1)./thickness,0.0,0.0,(zeta.*5.0e-1-eta.^2.*zeta.*5.0e-1-eta.^2.*5.0e-1+5.0e-1)./height,(eta.*(nu+1.0).*(zeta+1.0).*-1.0)./width,((nu+1.0).*(eta.^2-1.0).*-5.0e-1)./thickness,0.0,(zeta.*(eta+1.0).*(nu+1.0).*-1.0)./thickness,0.0,0.0,0.0,(eta.*5.0e-1-eta.*zeta.^2.*5.0e-1-zeta.^2.*5.0e-1+5.0e-1)./height,((nu+1.0).*(zeta.^2-1.0).*-5.0e-1)./width,0.0,((nu+1.0).*(zeta.^2-1.0).*-5.0e-1)./width,0.0,(eta.*5.0e-1-eta.*zeta.^2.*5.0e-1-zeta.^2.*5.0e-1+5.0e-1)./height,0.0,(zeta.*(eta+1.0).*(nu+1.0).*-1.0)./thickness,0.0,0.0,(eta.*5.0e-1-eta.*zeta.^2.*5.0e-1-zeta.^2.*5.0e-1+5.0e-1)./height,((nu+1.0).*(zeta.^2-1.0).*-5.0e-1)./width,(zeta.*(eta+1.0).*(nu+1.0).*-1.0)./thickness,0.0,((nu+1.0).*(eta.^2-1.0).*5.0e-1)./thickness,0.0,0.0,0.0,(zeta.*-5.0e-1+eta.^2.*zeta.*5.0e-1-eta.^2.*5.0e-1+5.0e-1)./height,(eta.*(nu+1.0).*(zeta-1.0).*1.0)./width,0.0,(eta.*(nu+1.0).*(zeta-1.0).*1.0)./width,0.0,(zeta.*-5.0e-1+eta.^2.*zeta.*5.0e-1-eta.^2.*5.0e-1+5.0e-1)./height,0.0,((nu+1.0).*(eta.^2-1.0).*5.0e-1)./thickness,0.0,0.0,(zeta.*-5.0e-1+eta.^2.*zeta.*5.0e-1-eta.^2.*5.0e-1+5.0e-1)./height,(eta.*(nu+1.0).*(zeta-1.0).*1.0)./width,((nu+1.0).*(eta.^2-1.0).*5.0e-1)./thickness,0.0,(zeta.*(eta-1.0).*(nu+1.0).*1.0)./thickness,0.0,0.0,0.0,(eta.*-5.0e-1+eta.*zeta.^2.*5.0e-1-zeta.^2.*5.0e-1+5.0e-1)./height,((nu+1.0).*(zeta.^2-1.0).*5.0e-1)./width,0.0,((nu+1.0).*(zeta.^2-1.0).*5.0e-1)./width,0.0,(eta.*-5.0e-1+eta.*zeta.^2.*5.0e-1-zeta.^2.*5.0e-1+5.0e-1)./height,0.0,(zeta.*(eta-1.0).*(nu+1.0).*1.0)./thickness,0.0,0.0,(eta.*-5.0e-1+eta.*zeta.^2.*5.0e-1-zeta.^2.*5.0e-1+5.0e-1)./height,((nu+1.0).*(zeta.^2-1.0).*5.0e-1)./width,(zeta.*(eta-1.0).*(nu+1.0).*1.0)./thickness,0.0,((eta-1.0).*(nu.^2-1.0).*5.0e-1)./thickness,0.0,0.0,0.0,(nu.*-1.0+eta.*nu.*1.0-nu.*zeta.*1.0+eta.*nu.*zeta.*1.0)./height,((zeta+1.0).*(nu.^2-1.0).*5.0e-1)./width,0.0,((zeta+1.0).*(nu.^2-1.0).*5.0e-1)./width,0.0,(nu.*-1.0+eta.*nu.*1.0-nu.*zeta.*1.0+eta.*nu.*zeta.*1.0)./height,0.0,((eta-1.0).*(nu.^2-1.0).*5.0e-1)./thickness,0.0,0.0,(nu.*-1.0+eta.*nu.*1.0-nu.*zeta.*1.0+eta.*nu.*zeta.*1.0)./height,((zeta+1.0).*(nu.^2-1.0).*5.0e-1)./width,((eta-1.0).*(nu.^2-1.0).*5.0e-1)./thickness,0.0,((eta+1.0).*(nu.^2-1.0).*-5.0e-1)./thickness,0.0,0.0,0.0,(nu.*-1.0-eta.*nu.*1.0-nu.*zeta.*1.0-eta.*nu.*zeta.*1.0)./height,((zeta+1.0).*(nu.^2-1.0).*-5.0e-1)./width,0.0,((zeta+1.0).*(nu.^2-1.0).*-5.0e-1)./width,0.0,(nu.*-1.0-eta.*nu.*1.0-nu.*zeta.*1.0-eta.*nu.*zeta.*1.0)./height,0.0,((eta+1.0).*(nu.^2-1.0).*-5.0e-1)./thickness,0.0,0.0,(nu.*-1.0-eta.*nu.*1.0-nu.*zeta.*1.0-eta.*nu.*zeta.*1.0)./height,((zeta+1.0).*(nu.^2-1.0).*-5.0e-1)./width,((eta+1.0).*(nu.^2-1.0).*-5.0e-1)./thickness,0.0,((eta+1.0).*(nu.^2-1.0).*5.0e-1)./thickness,0.0,0.0,0.0,(nu.*-1.0-eta.*nu.*1.0+nu.*zeta.*1.0+eta.*nu.*zeta.*1.0)./height,((zeta-1.0).*(nu.^2-1.0).*5.0e-1)./width,0.0,((zeta-1.0).*(nu.^2-1.0).*5.0e-1)./width,0.0,(nu.*-1.0-eta.*nu.*1.0+nu.*zeta.*1.0+eta.*nu.*zeta.*1.0)./height,0.0,((eta+1.0).*(nu.^2-1.0).*5.0e-1)./thickness,0.0,0.0,(nu.*-1.0-eta.*nu.*1.0+nu.*zeta.*1.0+eta.*nu.*zeta.*1.0)./height,((zeta-1.0).*(nu.^2-1.0).*5.0e-1)./width,((eta+1.0).*(nu.^2-1.0).*5.0e-1)./thickness,0.0,((eta-1.0).*(nu.^2-1.0).*-5.0e-1)./thickness,0.0,0.0,0.0,(nu.*-1.0+eta.*nu.*1.0+nu.*zeta.*1.0-eta.*nu.*zeta.*1.0)./height,((zeta-1.0).*(nu.^2-1.0).*-5.0e-1)./width,0.0,((zeta-1.0).*(nu.^2-1.0).*-5.0e-1)./width,0.0,(nu.*-1.0+eta.*nu.*1.0+nu.*zeta.*1.0-eta.*nu.*zeta.*1.0)./height,0.0,((eta-1.0).*(nu.^2-1.0).*-5.0e-1)./thickness,0.0,0.0,(nu.*-1.0+eta.*nu.*1.0+nu.*zeta.*1.0-eta.*nu.*zeta.*1.0)./height,((zeta-1.0).*(nu.^2-1.0).*-5.0e-1)./width,((eta-1.0).*(nu.^2-1.0).*-5.0e-1)./thickness,0.0],[6,60]);

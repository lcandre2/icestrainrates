%example script for stake advection
x=0:0.01:1;y=0:0.01:1;
[X,Y]=meshgrid(x,y);
%just making up a velocity field here to create velocity grid
u=0.5+(X-2)./((X-2).^2+Y.^2);v=Y./((X-2).^2+Y.^2);
%velocity vector plot
% quiver(X,Y,u,v)
% xlim([0 1]);ylim([0 1])
%streamline plot (starting the streamlines at x=1 and y=y)
%this is just for visualization of the velocity field
streamline(X,Y,u,v,ones(101,1),y)
hold on
%positioning and advecting an array of particles/stakes
%initialize stake grid (5*5=25 stakes)
xp=0.2:0.1:0.6;yp=0.1:0.1:0.5;
[Xpi,Ypi]=meshgrid(xp,yp);
Xp=Xpi;Yp=Ypi;
plot(Xp,Yp, 'r*')
hold on
dt=0.1;nt=20;%time step for advection and number of steps
for it=1:nt%advecting for nt time-steps
    plot(Xp,Yp,'c.')%just for tracking motion to visualize - would not use in actual code
    pause(0.2)
    %integrating dx/dt=u,dy/dt=v over a time step dt using Improved Euler Integration
    %look up interp2 - interpolating from velocity grid to stake position
    %up, vp - velocities at positions at beginning of time-step
    up=interp2(X,Y,u,Xp,Yp);vp=interp2(X,Y,v,Xp,Yp);
    %first step of Improved Euler is Forward Euler 
    Xpstar=Xp+up*dt;Ypstar=Yp+vp*dt;
    %second step of Improved Euler
    upstar=interp2(X,Y,u,Xpstar,Ypstar);vpstar=interp2(X,Y,v,Xpstar,Ypstar);
    %final positions after each Improved Euler Time Step
    Xp=Xp+(up+upstar)/2*dt;Yp=Yp+(vp+vpstar)/2*dt;
end
%plot final position after nt timesteps
plot(Xp,Yp,'r*')
    %now, can use [Xpi,Ypi] and [Xp,Yp] to make strain rate calculations by
    %choosing different diamonds from the stake array
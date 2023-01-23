%example script for stake advection and strain rate calculation
%illustrating on a single stake diamond
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%making a hypotherical grid and hypothetical velocity dataset
%hypothetical x-y grid
x=0:0.01:1;y=0:0.01:1;
[X,Y]=meshgrid(x,y);
%just making up a hypothetical velocity field here to create velocity grid
u=0.05+0.1*(X-2)./((X-2).^2+Y.^2);v=0.1*Y./((X-2).^2+Y.^2);
%velocity vector plot
 %quiver(X,Y,u,v,0)
% xlim([0 1]);ylim([0 1])
%streamline plot (starting the streamlines at x=1 and y=y)
%this is just for visualization of the velocity field - we don't need this
%in the final code
%streamline(X,Y,u,v,ones(101,1),y)
%hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for illustration take the coordinates of the reference point associated
%with the ROI as (0.5,0.5).  The diamond is centered at the same point
xroi=0.5;yroi=0.5;
%now specify the half-length l of the diamond's legs and the orientation
%alpha in degrees
%made a loop here - taking alpha = 9,18,27......90 degrees (11 entries)
alphas=0:9:90;
eresults=zeros(length(alphas),3);%initialize matrix into which to write out strain rate components
for i=1:length(alphas)
l=0.1;alpha=alphas(i);%specifying l and alpha is from the array
lof1=l*cosd(alpha-90);lof2=l*sind(alpha-90);
%creating a single diamond (5 points) with axes along x and y
%the 5 points of the diamond are at (0.5,0.5),(0.4,0.5),(0.6,0.5),(0.5,0.4),(0.5,0.6)
%Xpi=[0.5,0.4,0.6,0.5,0.5];Ypi=[0.5,0.5,0.5,0.4,0.6];%previous code
Xpi=[xroi,xroi-lof1,xroi+lof1,xroi+lof2,xroi-lof2];
Ypi=[yroi,yroi-lof2,yroi+lof2,yroi-lof1,yroi+lof1];
%plotting the points of the diamond to visualize
Xp=Xpi;Yp=Ypi;
% plot(Xp,Yp, 'r*')
% hold on
%%%NO CHANGE FROM HERE TO.....
%now advecting the 5 points of the diamond
dt=0.001;nt=50;%hypothetical time step for advection and number of steps
%advection loop
for it=1:nt%advecting for nt time-steps
    %integrating dx/dt=u,dy/dt=v over a time step dt using Improved Euler Integration
    %using interp2 - interpolating from velocity grid to stake position
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
% plot(Xp,Yp,'ro')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now on to the strain rate calculations - calculating the lengths of
%different segments in the strain diamond (a,b,c,d)
%initial lengths based on initial coordinates
%(we know a10, a20, c10, c20 should equal l, others equal l*sqrt(2))
%so, these may not be necessary, but can use as a sanity check
a10=sqrt((Xpi(5)-Xpi(1))^2+(Ypi(5)-Ypi(1))^2);
a20=sqrt((Xpi(4)-Xpi(1))^2+(Ypi(4)-Ypi(1))^2);
b10=sqrt((Xpi(5)-Xpi(2))^2+(Ypi(5)-Ypi(2))^2);
b20=sqrt((Xpi(4)-Xpi(3))^2+(Ypi(4)-Ypi(3))^2);
c10=sqrt((Xpi(2)-Xpi(1))^2+(Ypi(2)-Ypi(1))^2);
c20=sqrt((Xpi(3)-Xpi(1))^2+(Ypi(3)-Ypi(1))^2);
d10=sqrt((Xpi(4)-Xpi(2))^2+(Ypi(4)-Ypi(2))^2);
d20=sqrt((Xpi(5)-Xpi(3))^2+(Ypi(5)-Ypi(3))^2);
%final lengths based on final coordinates
a1f=sqrt((Xp(5)-Xp(1))^2+(Yp(5)-Yp(1))^2);
a2f=sqrt((Xp(4)-Xp(1))^2+(Yp(4)-Yp(1))^2);
b1f=sqrt((Xp(5)-Xp(2))^2+(Yp(5)-Yp(2))^2);
b2f=sqrt((Xp(4)-Xp(3))^2+(Yp(4)-Yp(3))^2);
c1f=sqrt((Xp(2)-Xp(1))^2+(Yp(2)-Yp(1))^2);
c2f=sqrt((Xp(3)-Xp(1))^2+(Yp(3)-Yp(1))^2);
d1f=sqrt((Xp(4)-Xp(2))^2+(Yp(4)-Yp(2))^2);
d2f=sqrt((Xp(5)-Xp(3))^2+(Yp(5)-Yp(3))^2);
%next calculating logarithmic strain rates
tott=nt*dt;%total time for virtual advection
edot0=0.5/tott*(log(a1f/a10)+log(a2f/a10));
edot45=0.5/tott*(log(b1f/b10)+log(b2f/b10));
edot90=0.5/tott*(log(c1f/c10)+log(c2f/c10));
edot135=0.5/tott*(log(d1f/d10)+log(d2f/d10));
%%%%%%NO CHANGE TO HERE.....
%finally calculate strain rate components
%general A matrix used here in place of old matrix
%a legs
ca=cosd(alpha);sa=sind(alpha);
%b legs
cam45=cosd(alpha)/sqrt(2)+sind(alpha)/sqrt(2);%cos(alpha-45)
sam45=-cosd(alpha)/sqrt(2)+sind(alpha)/sqrt(2);%sin(alpha-45)
%c legs
cam90=sind(alpha);%cos(alpha-90)
sam90=-cosd(alpha);%sin(alpha-90)
%d legs
cap45=cosd(alpha)/sqrt(2)-sind(alpha)/sqrt(2);%cos(alpha+45)
sap45=cosd(alpha)/sqrt(2)+sind(alpha)/sqrt(2);%sin(alpha+45)
%
A=[ca^2,2*ca*sa,sa^2;cam45^2,2*cam45*sam45,sam45^2;...
    cam90^2,2*cam90*sam90,sam90^2;cap45^2,2*cap45*sap45,sap45^2];
%lsqmatrix=[-1/4,1/4,3/4,1/4;0,1/2,0,1/2;3/4,1/4,-1/4,1/4];%previous
%version - there was a typo, in the second row, needed -1/2 in the last
%column
lsqmatrix=inv(A'*A)*A';
edotvector=lsqmatrix*[edot0;edot45;edot90;edot135];
%components of the strain rate tensor
exxdot=edotvector(1);
exydot=edotvector(2);
eyydot=edotvector(3);
eresults(i,:)=[exxdot,exydot,eyydot];%write estimated components into row of eresults
% hold off
end
eresults


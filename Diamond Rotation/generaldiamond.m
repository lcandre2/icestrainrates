 %alpha denotes the orientation of the "a" legs of the strain diamond with respect to the x-axis
%For the standard case used in Karen's paper, alpha is 90 degrees
%construting the A matrix following the explanation in the supplementary material in
%Karen's paper
%input alpha
alpha=27;%cah change alpha here - with alpha = 90 degrees, recover the A matrix in Karen's paper
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
lsqmatrix=inv(A'*A)*A';
%PART 2
%checking formulas for coordinates of the diamond based on alpha
%for illustration take the coordinates of the reference point associated
%with the ROI as (0.5,0.5).  The diamond is centered at the same point
xroi=0.5;yroi=0.5;
%now specify the half-length l of the diamond's legs and the orientation alpha
l=0.1;alpha=10;
%calculate offsets
lof1=l*cosd(alpha-90);lof2=l*sind(alpha-90);
%create coordinates of diamond - same convention as in the stake_example code
Xpi=[xroi,xroi-lof1,xroi+lof1,xroi+lof2,xroi-lof2];
Ypi=[yroi,yroi-lof2,yroi+lof2,yroi-lof1,yroi+lof1];

function [intXVel,intYVel] = locInterp2(rowCoord,colCoord,sqVx,sqVy)
% This function carries out a bilinear interpolation at the points
% (rowCoord, colCoord) within a local square. In this case, it is written
% to give back the interpolated values withing a local square of
% x-direction and y-direction velocity values.

ULxVel = sqVx(floor(rowCoord),floor(colCoord));
URxVel = sqVx(floor(rowCoord),ceil(colCoord));
LLxVel = sqVx(ceil(rowCoord),floor(colCoord));
LRxVel = sqVx(ceil(rowCoord),ceil(colCoord));

ULyVel = sqVy(floor(rowCoord),floor(colCoord));
URyVel = sqVy(floor(rowCoord),ceil(colCoord));
LLyVel = sqVy(ceil(rowCoord),floor(colCoord));
LRyVel = sqVy(ceil(rowCoord),ceil(colCoord));

topInterp = colCoord - floor(colCoord);
sideInterp = rowCoord - floor(rowCoord);

topxVel = ULxVel + (URxVel-ULxVel)*topInterp;
botxVel = LLxVel + (LRxVel-LLxVel)*topInterp;

topyVel = ULyVel + (URyVel-ULyVel)*topInterp;
botyVel = LLyVel + (LRyVel-LLyVel)*topInterp;

intXVel = topxVel + (botxVel-topxVel)*sideInterp;
intYVel = topyVel + (botyVel-topyVel)*sideInterp;

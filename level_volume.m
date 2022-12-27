function l=level_volume(Volume,level,V)

%function that linearly interpolates the level values versus volume values within the reservoir

%INPUT PARAMETERS
%Volume:....Volume data points at which the volume rating curve is given [m^3]
%level:.....level data points at which the volume rating curve is given [m]
%V:.........volume at which the function evaluate the corresponding level within the reservoir  [m^3]


%OUTPUT
%l=level correspondent to the "V" value

index=find(Volume<=V,1,'last');
a = level(index);
b = level(index+1)-level(index);
c = V-Volume(index);
d = Volume(index+1)-Volume(index);
l = a + b*c/d;
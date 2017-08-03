function varargout = polyfill(varargin)
% fill polygon with line or color, plot border?
%
% EXAMPLES
%     opts = struct('angdis', [60 1],...      : deg dis
%                   'lcolor', [.7 .6 .5],...  : vector or char
%                   'fcolor', 'b',..          : vector or char
%                   'bcolor', 'r');           : vector or char
%     polyfill([px;py],opts);
%
%     polyfill([px;py],'b',[10 7.5]);
% points; data points must be anticlockwise
% angdis: filling line angle and distance
% lcolor: filling line color
% fcolor: filling color
% bcolor: border line color

if (nargin < 1)
  error('Too few input arguments');
end

% polyfill(points, [options,] ...)
points = varargin{1};           % data points [x;y]

paramPairs = {varargin{2:end}}; % opts = paramPairs
if nargin > 1
  if isstruct(paramPairs{1})    % is a struct
    pcell = LocalToCell(paramPairs{1}); % LocalToCell convert a struct to {field1,val1,field2,val2,...}
    paramPairs = {pcell{:}, paramPairs{2:end}};
  end
end

% do some validity checking on param-value pairs
if (rem(length(paramPairs),2) ~= 0)
  error(['Invalid input syntax. Optional parameters and values' ...
	 ' must be in pairs.']);
end

angdis = -1;
lcolor = -1;
fcolor = -1;
bcolor = -1;

% Process param-value pairs
args = {};
for k = 1:2:length(paramPairs)
    param = lower(paramPairs{k});
    if ~ischar(param)        % parameter name is string
        error('Optional parameter names must be strings');
    end
    value = paramPairs{k+1}; % value of the parameter
    
    switch (param)
        case 'angdis'        % numeric vector 2
            angdis = value;
        case 'lcolor'        % numeric vector 3 or char
            lcolor = colr2vect(value);
        case 'fcolor'        % numeric vector 3 or char
            fcolor = colr2vect(value);
        case 'bcolor'        % numeric vector 3 or char
            bcolor = colr2vect(value);            
        otherwise
            error(['Unrecognized option ' param '.']);
    end
end

lnwdth = 0.4;

lthp = size(points,2);
newp = [points points(:,1)];
hold on
if fcolor(1) ~= -1
    fill(newp(1,:),newp(2,:),'color',fcolor);
end
if bcolor(1) ~= -1 
    plot(newp(1,:),newp(2,:),'color',bcolor);
end

if length(angdis) > 1
% (1) check crossed, anticlockwise
for i = 1:lthp-1
    for j = i+1:lthp
        if j==i+1 % in line check
            abc1 = abcl(newp(:,i:i+1));
            abc2 = abcl(newp(:,j:j+1));
            if abc1(1)*abc2(2)==abc1(2)*abc2(1)
                error('The polygon is CROSSED!');
            end
        elseif j==lthp & i==1 % in line check
            abc1 = abcl(newp(:,i:i+1));
            abc2 = abcl(newp(:,j:j+1));
            if abc1(1)*abc2(2)==abc1(2)*abc2(1)
                error('The polygon IS crossed!');
            end
        else
            if p_lline(newp(:,i:i+1),newp(:,j:j+1))
                error('The polygon IS crossed!');
            end
        end
    end
end
% check anticlockwise
newp = [points(:,end) newp];
anglep = zeros(1,lthp);
for i = 1:lthp
    anglep(i) = anglell(newp(:,i:i+2));
end
anglea = sum(anglep);
if anglea==180*(lthp-2) % not anticlockwise
    points = [points(:,1) fliplr(points(:,2:end))];
    anglep = [anglep(1) fliplr(anglep(2:end))];
    anglep = 360-anglep;
end

% (2) add triangle to be convex polygon
% mark the angle
markp = zeros(1,lthp);
j = 1;
for i = 1:lthp
    if anglep(i) < 180
        markp(j) = i;
        j = j+1;
        i = i+1;
    end
end
% translate to convex polygon
if j==1
    sgnpg = 0;
else
    ltri = j-1;
    trip = zeros(ltri,3);
    lisp = [lthp 1:lthp 1];
    lasp = newp;
    k = -1;
    for i = 1:ltri
        trip(i,:) = lisp(markp(i):markp(i)+2);
        k = k+1;
        lasp = [lasp(:,1:markp(i)-k) lasp(:,markp(i)-k+2:lthp-k+2)];
    end
    sgnpg = 1;
end

% (3) plot line
xmin = min(points(1,:));
xmax = max(points(1,:));
ymin = min(points(2,:));
ymax = max(points(2,:));
if angdis(1) < 90
    a0 = -tan(angdis(1)/180*pi);
    b0 = 1;
    c0 = -(a0*xmax+ymin+angdis(2)/2);
    cm = -(a0*xmin+ymax);
elseif angdis(1)==90
    a0 = 1;
    b0 = 0;
    c0 = -(xmin+angdis(2)/2);
    cm = -xmax;
else
    a0 = -tan(angdis(1)/180*pi);
    b0 = 1;
    c0 = -(a0*xmin+ymin+angdis(2)/2);
    cm = -(a0*xmax+ymax);
end
newp = [points points(:,1)];
if lcolor == -1
    lcolor = [0.7 0.6 0.5];
end
while c0 > cm
    abc0 = [a0 b0 c0];
    crsp = zeros(2,lthp);
    j = 1;
    for i = 1:lthp
        two = newp(:,i:i+1);
        pln = p_line(abc0,two);
        if pln(1)==1
            crsp(1,j) = pln(2);
            crsp(2,j) = pln(3);
            j = j+1;
        end
    end
    if j~=1 % have cross points, plot line
        if sgnpg==0 % no added triangel
            plot(crsp(1,1:2),crsp(2,1:2),'color',lcolor,'linewidth',lnwdth)
        else
            cros = zeros(2,j-1);
            xm = max(crsp(1,:));
            ym = max(crsp(2,:));
            if crsp(1,1)==crsp(1,2)
                xy = 2;
            else
                xy = 1;
            end
            for i = 1:j-1 % sort
                [cros(xy,i) dn] = min(crsp(xy,1:j-1));
                cros(3-xy,i) = crsp(3-xy,dn);
                crsp(1,dn) = xm+1;
                crsp(2,dn) = ym+1;
            end
            for i=1:j-2
                sgnp = 0; % cross with triangle or not
                two = cros(:,i:i+1);
                tri = zeros(2,3);
                for k = 1:ltri
                    tri(:,1) = points(:,trip(k,1));
                    tri(:,2) = points(:,trip(k,2));
                    tri(:,3) = points(:,trip(k,3));
                    if line_tri(two,tri)==1
                        sgnp = 1;
                        break;
                    end
                end
                if sgnp~=1
                    plot(two(1,1:2),two(2,1:2),'color',lcolor,'linewidth',lnwdth)
                end
            end
        end
    end
    c0 = c0-angdis(2);
end
end  % end line plot
hold off

% subfunctions ------------------------ 
% sort two points x & y from min to max
function newp = sortxy(oldp)
newp = oldp;
if oldp(1,1)>oldp(1,2)
    newp(1,1) = oldp(1,2);
    newp(1,2) = oldp(1,1);
end
if oldp(2,1)>oldp(2,2)
    newp(2,1) = oldp(2,2);
    newp(2,2) = oldp(2,1);
end

% get line parameters ax+by+c=0 from two points
function  abc = abcl(two)     
abc = zeros(1,3);
if two(1,1) == two(1,2)
    abc = [1 0 -two(1,1)];
else 
    abc(1) = (two(2,2)-two(2,1))/(two(1,1)-two(1,2));
    abc(2) = 1;
    abc(3) = (two(1,2)*two(2,1)-two(1,1)*two(2,2))/(two(1,1)-two(1,2));
end

% get the angle of line
function angle = anglel(two)     
abc = abcl(two);
if abc(1) == 0
    angle = 0;
elseif abc(2) == 0
    angle = 90;
else 
    vcos =( -sign(abc(1))*abc(2))/(abc(1)^2+abc(2)^2)^0.5;
    angle = 180/pi*acos(vcos);
end

% get the angle of two radials, anticlockwise
function angle = anglell(tri)     
angle1 = anglel(tri(:,1:2));
angle2 = anglel(tri(:,2:3));
trn = tri;
trn(1,:) = tri(1,:)-tri(1,2);
trn(2,:) = tri(2,:)-tri(2,2);
if (trn(1,1)<0 & trn(2,1)<=0)|(trn(1,1)>=0 & trn(2,1)<0)
    angle1 = angle1+180;
end
if (trn(1,3)<0 & trn(2,3)<=0)|(trn(1,3)>=0 & trn(2,3)<0)
    angle2 = angle2+180;
end
angle = angle2-angle1;
if angle < 0
    angle = angle+360;
end

% get the intersect point of line and radial, not include superposition
function sgn_p = p_line(abc,two)
sgn_p = zeros(1,3);
abc2  = abcl(two);
a1 = abc(1);
a2 = abc2(1);
b1 = abc(2);
b2 = abc2(2);
c1 = abc(3);
c2 = abc2(3);
if a1*b2~=a2*b1
    x = (b1*c2-b2*c1)/(a1*b2-a2*b1);
    y = (a2*c1-a1*c2)/(a1*b2-a2*b1);
    two = sortxy(two);
%   if (two(1,1)<=x & x<=two(1,2)) & (two(2,1)<=y & y<=two(2,2))  
%   here; matlab has a error less than 1e-15!!!
    if (x-two(1,1)>=-1e-15 & x-two(1,2)<=1e-15) & (y-two(2,1)>=-1e-15 & y-two(2,2)<=1e-15)
        sgn_p(1) = 1;
        sgn_p(2) = x;
        sgn_p(3) = y;
    end
end
     
% is radial and radial intersect? not include end points
function sgn = p_lline(tw1,tw2)
sgn  = 0;
abc1 = abcl(tw1);
abc2 = abcl(tw2);
a1 = abc1(1); a2 = abc2(1);
b1 = abc1(2); b2 = abc2(2);
c1 = abc1(3); c2 = abc2(3);
if a1*b2~=a2*b1
    x = (b1*c2-b2*c1)/(a1*b2-a2*b1);
    y = (a2*c1-a1*c2)/(a1*b2-a2*b1);
    sgn = 1;
    tw1 = sortxy(tw1);
    tw2 = sortxy(tw2);
    if tw1(1,1)~=tw1(1,2) & (x<=tw1(1,1) | tw1(1,2)<=x)     % tw1 x
        sgn = 0;
    elseif tw1(2,1)~=tw1(2,2) & (y<=tw1(2,1) | tw1(2,2)<=y) % tw1 y
        sgn = 0;
    elseif tw2(1,1)~=tw2(1,2) & (x<=tw2(1,1) | tw2(1,2)<=x) % tw2 x
        sgn = 0;
    elseif tw2(2,1)~=tw2(2,2) & (y<=tw2(2,1) | tw2(2,2)<=y) % tw2 y
        sgn = 0;
    end
end

% is radial and triangle intersect? include end points
function sgn = line_tri(two,tri)
sgn = 0;
cen = sum(two,2)/2;
pl1 = tri(:,1:2);
pl2 = tri(:,2:3);
pl3 = [tri(:,3) tri(:,1)];
abc1= abcl([cen tri(:,3)]);
abc2= abcl([cen tri(:,1)]);
abc3= abcl([cen tri(:,2)]);
sg  = zeros(1,3);
tmpva = p_line(abc1,pl1);
sg(1) = tmpva(1);
tmpva = p_line(abc2,pl2);
sg(2) = tmpva(1);
tmpva = p_line(abc3,pl3);
sg(3) = tmpva(1);
if sum(sg)>=3
    sgn  = 1;
end

% change color to three elements vector
function color = colr2vect(value)
color = [0.3,0.3,0.3];
if ischar(value) 
    switch (value)
        case 'y' 
            color = [1 1 0]; % yellow
        case 'm' 
            color = [1 0 1]; % magenta 
        case 'c' 
            color = [0 1 1]; % cyan
        case 'r' 
            color = [1 0 0]; % red
        case 'g' 
            color = [0 1 0]; % green
        case 'b' 
            color = [0 0 1]; % blue
        case 'w' 
            color = [1 1 1]; % white 
        case 'k' 
            color = [0 0 0]; % black
    end
else    
    color = value;
end

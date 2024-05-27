function error_area(X,Y,barlength,color,alpha,varargin)

X = X(:);
Y = Y(:);
barlength = barlength(:);

if nargin == 6
    linestyle = varargin{1};
    linewidth = 1;
    
elseif nargin == 7
    linestyle = varargin{1};
    linewidth = varargin{2};
else
    linestyle = '-';
    linewidth = 1;
end

Yabove = Y+barlength;
Ylow = Y-barlength;
X = X';
Yabove = Yabove';
Ylow = Ylow';

fill([X(1) X(1:end) fliplr([X(1:end) X(end)])],[Yabove(1) Yabove fliplr([Ylow Ylow(end)])],color);
hold on
plot(X,Y,'color',color,'linewidth',linewidth); 
h = get(gca,'children');

if nargin>6
    set(h(2),'facealpha',alpha,'linestyle','none');
    set(h(1),'linestyle',linestyle,'linewidth',linewidth);

else
    set(h(2),'facealpha',alpha,'linestyle','none');
    set(h(1),'linestyle','none','linewidth',linewidth);%here put the linestyle to none instead of '-' so we dont see the line

end
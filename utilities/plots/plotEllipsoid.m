
function plotEllipsoid(x0,y0, horizontal_radius, vertical_radius,color)
    hold on
    a=horizontal_radius; % horizontal radius
    b=vertical_radius; % vertical radius
    t=-pi:0.01:pi;
    x=x0+a*cos(t);
    y=y0+b*sin(t);
    fill(x,y,color,'EdgeColor','none','FaceColor',color);
end
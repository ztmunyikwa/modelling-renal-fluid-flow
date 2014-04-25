function dxy = diffxy(t, xy)

%split xy into variables in our equations

x = xy(1);
xdot = xy(2);
y = xy(3);


% define the derivatives of these variables from equations
xdot = xdot;
ydot = 3*x + 2*y + 5;
xdoubledot = 3 - ydot + 2*xdot;
%



%return the derivatives in dxy in the right order



dxy = [xdot; xdoubledot; ydot]



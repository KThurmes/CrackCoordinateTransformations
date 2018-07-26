function [Grid] = ContourMe_I_nint(xfrom, xto, Nx, yfrom, yto, Ny, func,nint)
%==========================================================================
% ContourMe_I(xfrom, xto, Nx, yfrom, yto, Ny, func)                (01.23.09)
%
%   Contour the imaginary part of the specified complex function.
%
% Arguments:
%
%   xfrom   starting x-value for the domain
%   xto     ending x-value for the domain
%   Nx      number of grid columns
%
%   yfrom   starting y-value for the domain
%   yto     ending y-value for the domain
%   Ny      number of grid rows
%
%   func    function to contour;  must take one complex argument.
%
% Returns:
%
%   Grid    Ny x Nx matrix of values of func at the rid nodes.
%
% Example Usage:
%
%   G = ContourMe_I(1,2,11,1,2,11,@(z)Omega(1,-1,z));
%==========================================================================
Grid = zeros(Ny,Nx);

X = linspace(xfrom, xto, Nx);
Y = linspace(yfrom, yto, Ny);

for row = 1:Ny
    for col = 1:Nx
        Grid(row,col) = func( complex( X(col), Y(row) ) );
    end
end

figure;
contour(X, Y, imag(Grid), nint, 'b');
axis equal
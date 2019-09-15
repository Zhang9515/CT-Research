function [p,ellipse]=Diskphantom(varargin)

%PHANTOM3D Three-dimensional analogue of MATLAB Shepp-Logan phantom
%   P = PHANTOM3D(DEF,N) generates a 3D head phantom that can   
%   be used to test 3-D reconstruction algorithms.
%
%   DEF is a string that specifies the type of head phantom to generate.
%   Valid values are: 
%         
%      'Shepp-Logan'            A test image used widely by researchers in
%                               tomography
%      'Modified Shepp-Logan'   (default) A variant of the Shepp-Logan phantom
%                               in which the contrast is improved for better  
%                               visual perception.
%
%   N is a scalar that specifies the grid size of P.
%   If you omit the argument, N defaults to 64.
% 
%   P = PHANTOM3D(E,N) generates a user-defined phantom, where each row
%   of the matrix E specifies an ellipsoid in the image.  E has ten columns,
%   with each column containing a different parameter for the ellipsoids:
%   
%     Column 1:  A      the additive intensity value of the ellipsoid
%     Column 2:  a      the length of the x semi-axis of the ellipsoid 
%     Column 3:  b      the length of the y semi-axis of the ellipsoid
%     Column 4:  c      the length of the z semi-axis of the ellipsoid
%     Column 5:  x0     the x-coordinate of the center of the ellipsoid
%     Column 6:  y0     the y-coordinate of the center of the ellipsoid
%     Column 7:  z0     the z-coordinate of the center of the ellipsoid
%     Column 8:  phi    phi Euler angle (in degrees) (rotation about z-axis)

%
%   For purposes of generating the phantom, the domains for the x-, y-, and 
%   z-axes span [-1,1].  Columns 2 through 7 must be specified in terms
%   of this range.
%
%   [P,E] = PHANTOM3D(...) returns the matrix E used to generate the phantom.
%
%   Class Support
%   -------------
%   All inputs must be of class double.  All outputs are of class double.
%
%   Remarks
%   -------
%   For any given voxel in the output image, the voxel's value is equal to the
%   sum of the additive intensity values of all ellipsoids that the voxel is a 
%   part of.  If a voxel is not part of any ellipsoid, its value is 0.  
%
%   The additive intensity value A for an ellipsoid can be positive or negative;
%   if it is negative, the ellipsoid will be darker than the surrounding pixels.
%   Note that, depending on the values of A, some voxels may have values outside
%   the range [0,1].
%    
%   Example
%   -------
%        ph = phantom3d(128);
%        figure, imshow(squeeze(ph(64,:,:)))
%
%   Copyright 2005 Matthias Christian Schabel (matthias @ stanfordalumni . org)
%   University of Utah Department of Radiology
%   Utah Center for Advanced Imaging Research
%   729 Arapeen Drive
%   Salt Lake City, UT 84108-1218
%   
%   This code is released under the Gnu Public License (GPL). For more information, 
%   see : http://www.gnu.org/copyleft/gpl.html
%
%   Portions of this code are based on phantom.m, copyrighted by the Mathworks
%

[ellipse,n] = parse_inputs(varargin{:});

p = zeros([n n n]);

rng =  ( (0:n-1)-(n-1)/2 ) / ((n-1)/2); 

[x,y,z] = meshgrid(rng,rng,rng);

coord = [flatten(x); flatten(y); flatten(z)];

p = flatten(p);

for k = 1 : size(ellipse,1)    
   A = ellipse(k,1);            % Amplitude change for this ellipsoid
   asq = ellipse(k,2)^2;        % a^2
   bsq = ellipse(k,3)^2;        % b^2
   csq = ellipse(k,4)^2;        % c^2
   x0 = ellipse(k,5);           % x offset
   y0 = ellipse(k,6);           % y offset
   z0 = ellipse(k,7);           % z offset
   phi = ellipse(k,8)*pi/180;   % first Euler angle in radians
   
   cphi = cos(phi);
   sphi = sin(phi);
   
   % Euler rotation matrix
   alpha = [ cphi   sphi  0 ;
            -sphi  cphi 0 ;
            0     0    1 ];        
   
   % rotated ellipsoid coordinates
   coordp = alpha*coord;
   shift =  ( alpha * [ x0 ; y0 ; z0 ] )';
   idx = find((coordp(1,:)-shift(1)).^2./asq + (coordp(2,:)-shift(2)).^2./bsq + (coordp(3,:)-shift(3)).^2./csq <= 1);
   p(idx) = p(idx) + A;
end

p = reshape(p,[n n n]);

return;


function out = flatten(in)

out = reshape(in,[1 numel(in)]);

return;
   
   
function [e,n] = parse_inputs(varargin)
%  e is the m-by-10 array which defines ellipsoids
%  n is the size of the phantom brain image

n = 128;     % The default size
e = [];
defaults = {'shepp-logan', 'modified shepp-logan', 'yu-ye-wang', 'elimateonepair', 'translatezonepair'...
    , 'translatexonepair'};

for i=1:nargin
   if ischar(varargin{i})         % Look for a default phantom
      def = lower(varargin{i});
      idx = strmatch(def, defaults);
      if isempty(idx)
         eid = sprintf('Images:%s:unknownPhantom',mfilename);
         msg = 'Unknown default phantom selected.';
         error(eid,'%s',msg);
      end
      switch defaults{idx}
      case 'shepp-logan'
         e = shepp_logan;
      case 'modified shepp-logan'
         e = modified_shepp_logan;
      case 'yu-ye-wang'
         e = yu_ye_wang;
      case 'elimateonepair'
         e = ElimateOnePair;
      case 'translatezonepair'
         e = translatezonepair;
      case 'translatexonepair'
         e = translatexonepair;
      end
   elseif numel(varargin{i})==1 
      n = varargin{i};            % a scalar is the image size
   elseif ndims(varargin{i})==2 && size(varargin{i},2)==10 
      e = varargin{i};            % user specified phantom
   else
      eid = sprintf('Images:%s:invalidInputArgs',mfilename);
      msg = 'Invalid input arguments.';
      error(eid,'%s',msg);
   end
end

% ellipse is not yet defined
if isempty(e)                    
   e = Disk ;
end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Default head phantoms:   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
function e = Disk
%
%   This head phantom is the same as the Shepp-Logan except 
%   the intensities are changed to yield higher contrast in
%   the image.  Taken from Toft, 199-200.
%      
%         A      a     b     c     x0      y0      z0    phi 
%        -----------------------------------------------------
e =    [  1    .7  .7   .06      0       0       0        0      
            1    .7  .7   .06      0       0      .24        0
            1    .7  .7   .06      0       0       -.24        0
            1    .7  .7   .06      0       0       .48        0
            1    .7  .7   .06      0       0       -.48        0
            1    .7  .7   .06      0       0       .72        0
            1    .7  .7   .06      0       0       -.72      0 ] ;
       
return;
          
function e = ElimateOnePair
%
%   This head phantom is the same as the Shepp-Logan except 
%   the intensities are changed to yield higher contrast in
%   the image.  Taken from Toft, 199-200.
%      
%         A      a     b     c     x0      y0      z0    phi 
%        -----------------------------------------------------
e =    [  1    .7  .7   .06      0       0       0        0      
            1    .7  .7   .06      0       0      .24        0
            1    .7  .7   .06      0       0       -.24        0
            1    .7  .7   .06      0       0       .72        0
            1    .7  .7   .06      0       0       -.72      0 ] ;
       
return;

function e = translatezonepair
%
%   This head phantom is the same as the Shepp-Logan except 
%   the intensities are changed to yield higher contrast in
%   the image.  Taken from Toft, 199-200.
%      
%         A      a     b     c     x0      y0      z0    phi 
%        -----------------------------------------------------
e =    [  1    .7  .7   .06      0       0       0        0      
            1    .7  .7   .06      0       0      .24        0
            1    .7  .7   .06      0       0       -.24        0
            1    .7  .7   .06      0       0       .48        0
            1    .7  .7   .06      0       0       -.48        0
            1    .7  .7   .06      0       0       .72        0
            1    .7  .7   .06      0       0       -.72      0 ] ;
       
return;

function e = translatexonepair
%
%   This head phantom is the same as the Shepp-Logan except 
%   the intensities are changed to yield higher contrast in
%   the image.  Taken from Toft, 199-200.
%      
%         A      a     b     c     x0      y0      z0    phi 
%        -----------------------------------------------------
e =    [  1    .7  .7   .06      0       0       0        0      
            1    .7  .7   .06      0       0      .24        0
            1    .7  .7   .06      0       0       -.24        0
            1    .7  .7   .06      0.2       0       .48        0
            1    .7  .7   .06      0.2       0       -.48        0
            1    .7  .7   .06      0       0       .72        0
            1    .7  .7   .06      0       0       -.72      0 ] ;
       
return;

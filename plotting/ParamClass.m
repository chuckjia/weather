classdef ParamClass
    %PARAMCLASS Summary of this class goes here
    %   Detailed explanation goes here
    
   properties
      x0
      xf
      zf
      Nx
      Nz
      tf
      Nt
      num_csv
   end
   
   methods
      function param = ParamClass(x0, xf, zf, Nx, Nz, tf, Nt, num_csv)
         if nargin > 0
            param.x0 = x0;
            param.xf = xf;
            param.zf = zf;
            param.Nx = Nx;
            param.Nz = Nz;
            param.tf = tf;
            param.Nt = Nt;
            param.num_csv = num_csv;
         end
      end
   end

end


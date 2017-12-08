function [ ibound ] = findBoundary( A )
   nnu = size(A,1); 
   ibound = zeros(nnu,1);

   for i=1:nnu
      if ( A(i,i) == 1 ) 
         ibound(i) = 1;
      end 
   end
end


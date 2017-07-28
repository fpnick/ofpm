function [ measure ] = measureDiagDom( A )
%MEASUREDIAGDOM Summary of this function goes here
%   Detailed explanation goes here

    absDiag = abs(diag(A));
    absOffDiag = sum(abs(A),2) - absDiag;
    measure = min(absDiag - absOffDiag);

end


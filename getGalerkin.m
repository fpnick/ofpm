%> Computes Galerkin operator.
%>
%> @param[in]  A  Fine level matrix.
%> @param[in]  P  Prolongation.
%> 
%> @param[out] B  Galerkin operator.
function B = getGalerkin(A,P)
    B = P' * A * P;
end

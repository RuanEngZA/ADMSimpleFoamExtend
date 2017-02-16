{

	  tmp <fvVectorMatrix> UarbForEqn
	  (
	    fvm::div(phiFor, UarbForward) 
	    ==
	    fvm::laplacian(DS, UarbForward) + 
	    fvm::Sp(Large_for*1e30, UarbForward) -
	    Large_for*U*1e30 - 
	    fvm::Sp(DatStop*1e30, UarbForward)
	  );
	  
	  solve(UarbForEqn());
	  
	  tmp <fvVectorMatrix> UarbBackEqn
	  (
	    fvm::div(phiBack, UarbBackward) 
	    ==
	    fvm::laplacian(DS, UarbBackward) + 
	    fvm::Sp(Large_back*1e10, UarbBackward) -
	    Large_back*U*1e10 - 
	    fvm::Sp(DatStop*1e10, UarbBackward)
	  );
	 	  
	  solve(UarbBackEqn());

}
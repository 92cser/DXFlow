{
	Info << "Handle the cell whose alpha is slightly larger or smaller than 1.0" << endl;
	scalar c0 = 0, c1 = 0;
	forAll(alpha1,i)
	{	
		scalar error_abs = mag(alpha1[i] - 1.0);
		if((error_abs < 1e-5) || MassSource[i])   
		{		//2 times is enough
			for(int j = 0; j < 2; ++j)
			{
				c0 = mag(alpha1[i] - 1.0)/alpha1[i];
				c1 = 1.0 - c0;
				alpha1[i] = c0*alpha1[i] + c1;
			}
		}
		
		/*if(alpha1[i] > 1.0)   
		{		
			for(int j = 0; j < 2; ++j)
			{
				c0 = mag(alpha1[i] - 1.0)/alpha1[i];
				c1 = 1.0 - c0;
				alpha1[i] = c0*alpha1[i] + c1;
			}
		}

		if(alpha1[i] < 0.0)   
		{		
			for(int j = 0; j < 2; ++j)
			{
				alpha1[i] *= 1e-3;
			}
		}*/
	}
}


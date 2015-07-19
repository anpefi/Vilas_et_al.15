
/* selection.c */

# include "libhdr"
# include "ranlib.h"
# define NmaxI 101 // max ind
# define NmaxL 1001 // max loci 
# define NmaxA 257 // max alleles

/* ************************************************************************* */

int	a, i, j, l, h, m, z, t, n, w, k, gen, ind, loci, all, lociQ, alleles, allelesQ, cifras, cifrasQ;
int genotype[NmaxI][NmaxL][3], genotypeQ[NmaxI][NmaxL][3], parent[NmaxI][NmaxL][3], parentQ[NmaxI][NmaxL][3];
int RM[NmaxL], ind_sel, p1, p2, h1, h2, GEN, renumber[NmaxL][1001], code[NmaxL][1001], all_S[NmaxL], all_N[NmaxL];

double x;
double genval[NmaxI], efecto[NmaxL][NmaxA], genval[NmaxI], alelos, alelos_S;
double VE, AA, Aa, aa, HS, HS_S, q_s, q_s_S, pm_a[NmaxI], PM_A, GEN_VAL; 
double fQ[NmaxI][NmaxL][NmaxA], fsQ[NmaxL][NmaxA], min_efectQ[NmaxL], max_efectQ[NmaxL], rangoQ[NmaxL];
double f[NmaxI][NmaxL][NmaxA], fs[NmaxL][NmaxA];

FILE *fqtl, *fntr, *fphe, *fefe, *fout, *fdat, *floc;

main()
{
	getseed();	
	getinputs();
	headings();
	input_file();
	for (gen=0; gen<=GEN; gen++)
	{
		selected_genes();
		neutral_genes();
		fprintf(fout, "REP  %d	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f	%6.4f\n", gen, PM_A, GEN_VAL, HS_S,	HS, alelos_S, alelos);
		order();
		mating();	
	}
	fclose(fout);
//	fclose(floc);
	writeseed();
	return(0);
}

/* ************************************************************************* */

getinputs()
{
	getintandskip("Number of individuals (max 5000) for population:",&ind,1,5000);
	getintandskip("Number of individuals (max 5000) selected for mating:",&ind_sel,1,5000);
	VE=1.0;
	
	getintandskip("Number of generations:",&GEN,1,infinity);
}

/* ************************************************************************* */

headings()
{
	if ((fout = fopen ("data", "a"))==NULL)
	{
		printf ("No se pudo abrir el fichero data\n");
		exit (1);
	}
	
//	fprintf (fout, "gen	PM_A	GEN_VAL	HS_S	HS_N	all_S	all_N\n");
	
}

/* ************************************************************************* */

input_file()
{
	
	if ((fntr = fopen ("ntrl_genotype", "r"))==NULL)
	{
		printf ("No se pudo abrir el fichero ntrl_genotype\n");
		exit (1);
	}
	
	fscanf (fntr, "%d %d %d", &loci, &alleles, &cifras);
	
	for (i=1; i<=ind; i++)
	for (l=1; l<=loci; l++)
	for (h=1; h<=2; h++)
	{
		fscanf(fntr,"%03d", &a); 
		genotype[i][l][h] = a; 
	}
	fclose(fntr);
		
	if ((fqtl = fopen ("quanti_genotype", "r"))==NULL)
	{
		printf ("No se pudo abrir el fichero quanti_genotype\n");
		exit (1);
	}
	
	fscanf (fqtl, "%d %d %d", &lociQ, &allelesQ, &cifrasQ);
	
	for (i=1; i<=ind; i++)
	for (l=1; l<=lociQ; l++)
	for (h=1; h<=2; h++)
	{
		fscanf(fqtl,"%03d", &a); 
		genotypeQ[i][l][h] = a; 
	}
	fclose(fqtl);
	/* No hemos guardado el fenotipo
	if ((fphe = fopen ("quanti_phenotype", "r"))==NULL)
	{
		printf ("No se pudo abrir el fichero quanti_phenotype\n");
		exit (1);
	}
	
	for (i=1; i<=ind; i++)
	{
		fscanf (fphe, "%lf", &x);
		pm_a[i] = x;
	}
	fclose(fphe);
	*/	
	if ((fefe = fopen ("efecto", "r"))==NULL)
	{
		printf ("No se pudo abrir el fichero efecto\n");
		exit (1);
	}
			
	for (l=1; l<=lociQ; l++)
	for (all=1; all<=allelesQ; all++)
	{
		fscanf (fefe, "%lf", &x);
		efecto[l][all] = x;
	}
	fclose(fefe);	
		
}

/* ************************************************************************* */

selected_genes()
{
	PM_A = 0.0; GEN_VAL = 0.0;
	
	/*
	if (gen == 0)
	{
		for (i=1; i<=ind; i++)
		{
			PM_A += pm_a[i] / (double) ind;
		}
		
		for (i=1; i<=ind; i++)
		{
			for (l=1; l<=lociQ;l++)
			{
				for (h=1; h<=2; h++)		
				{
					genval[i] += efecto[l][genotypeQ[i][l][h]];
			 	}
			}

			GEN_VAL += genval[i] /(double) ind;
		}
	}
	
	else
	{*/
		for (i=1; i<=ind; i++)	genval[i] = 0.0;
		
		for (i=1; i<=ind; i++)
		{
			for (l=1; l<=lociQ;l++)
			{
				for (h=1; h<=2; h++)		
				{
					genval[i] += efecto[l][genotypeQ[i][l][h]];
			 	}
			}
			pm_a[i] = genval[i] + normal(0.0, sqrt(VE));
			PM_A += pm_a[i] / (double) ind;
			GEN_VAL += genval[i] /(double) ind;
		}
	//}
	
	HS_S = 0.0;  
	
	for (l=1; l<=lociQ; l++)
	{

		AA=0.0; Aa=0.0; 
    	
		for (i=1; i<=ind; i++)
		{
			if (genotypeQ[i][l][1]==genotypeQ[i][l][2])	AA+=1.0;
			else 										Aa+=1.0;			
		}

		q_s_S = (AA/(double)ind)+(Aa/(2.0*(double)ind));

		HS_S += (2.0 * q_s_S * (1.0 - q_s_S)) / (double)lociQ;			        	
	} 


	
	/* ***** order alleles ***** */
	for (l=1; l<=lociQ; l++)
	{
		z=0;
		for (i=1; i<=ind; i++)
		for (h=1;h<=2; h++)
		{
			z++;
			renumber[l][z] = genotypeQ[i][l][h];
		}
	}

	/* ***** new code ***** */
	
	for (l=1; l<=lociQ; l++)
	{
		for (t=1; t<(2*ind); t++)
		{
			for (n=(t+1); n<=(2*ind); n++)
			{
				if (renumber[l][t] > renumber[l][n])
				{
					w=renumber[l][t];
					renumber[l][t]=renumber[l][n];
					renumber[l][n]=w;
				}
			}
		}
		k=1;
		code[l][k]=renumber[l][1];
		for (t=2; t<=(2*ind); t++)
		{
			if (renumber[l][t]>renumber[l][t-1])
			{
				k++;
				code[l][k]=renumber[l][t];
			}
		}
		all_S[l]=k;
	}
	
	/* ***** frequencies for each individual ***** */
	
	for (l=1; l<=lociQ; l++)
	for (i=1; i<=ind; i++)
	for (k=1; k<=all_S[l]; k++)
	{
		fQ[i][l][k] = 0.0;
		
		for (h=1; h<=2; h++)
		{
			if (genotypeQ[i][l][h] == code[l][k])
			{
				fQ[i][l][k] += 0.5;
			}
		}
	}	
	
	/* Average frequency of populations */

	for (l=1; l<=lociQ; l++)
	for (k=1; k<=all_S[l]; k++)
	{
		fsQ[l][k] = 0.0;
		
		for(i=1; i<=ind; i++)
		{
			fsQ[l][k] += fQ[i][l][k]/(double)ind; 
		}
	}
	
	for (l=1; l<=lociQ; l++)
	{
		for (k=1; k<=all_S[l]; k++)
		{
			if (k == 1)
			{
				min_efectQ[l] = efecto[l][code[l][k]];
				max_efectQ[l] = efecto[l][code[l][k]];
			}
			else
			{
				if (max_efectQ[l] < efecto[l][code[l][k]])		max_efectQ[l] = efecto[l][code[l][k]];
				if (min_efectQ[l] > efecto[l][code[l][k]])		min_efectQ[l] = efecto[l][code[l][k]];
			}
		}
		rangoQ[l] = max_efectQ[l] - min_efectQ[l];	
	}

	
	alelos_S = 0.0;
	
	for (l=1; l<=lociQ; l++)
	{
		alelos_S += (double)all_S[l]/(double)lociQ;
	}
		
}

/* ************************************************************************* */

neutral_genes()
{
	HS = 0.0;  
	
	for (l=1; l<=loci; l++)
	{

		AA=0.0; Aa=0.0; 
    	
		for (i=1; i<=ind; i++)
		{
			if (genotype[i][l][1]==genotype[i][l][2])	AA+=1.0;
			else 										Aa+=1.0;			
		}

		q_s = (AA/(double)ind)+(Aa/(2.0*(double)ind));

		HS += (2.0 * q_s * (1.0 - q_s)) / (double)loci;			        	
	} 

	/* ***** order alleles ***** */
	
	for (l=1; l<=loci; l++)
	{
		z=0;
		for (i=1; i<=ind; i++)
		for (h=1;h<=2; h++)
		{
			z++;
			renumber[l][z] = genotype[i][l][h];
		}
	}

	/* ***** new code ***** */
	
	for (l=1; l<=loci; l++)
	{
		for (t=1; t<(2*ind); t++)
		{
			for (n=(t+1); n<=(2*ind); n++)
			{
				if (renumber[l][t] > renumber[l][n])
				{
					w=renumber[l][t];
					renumber[l][t]=renumber[l][n];
					renumber[l][n]=w;
				}
			}
		}
		k=1;
		code[l][k]=renumber[l][1];
		for (t=2; t<=(2*ind); t++)
		{
			if (renumber[l][t]>renumber[l][t-1])
			{
				k++;
				code[l][k]=renumber[l][t];
			}
		}
		all_N[l]=k;
	}
	
	/* ***** frequencies for each individual ***** */
	
	for (l=1; l<=loci; l++)
	for (i=1; i<=ind; i++)
	for (k=1; k<=all_N[l]; k++)
	{
		f[i][l][k] = 0.0;
		
		for (h=1; h<=2; h++)
		{
			if (genotype[i][l][h] == code[l][k])
			{
				f[i][l][k] += 0.5;
			}
		}
	}	
	
	/* Average frequency of populations */

	for (l=1; l<=loci; l++)
	for (k=1; k<=all_N[l]; k++)
	{
		fs[l][k] = 0.0;
		
		for(i=1; i<=ind; i++)
		{
			fs[l][k] += f[i][l][k]/(double)ind; 
		}
	}
	
	
	alelos = 0.0;
	
	for (l=1; l<=loci; l++)
	{
		alelos += (double)all_N[l]/(double)loci;
	}

}

/* ************************************************************************* */

order()
{
	for (i=1; i<ind; i++)
	for (j=i+1; j<=ind; j++)
	{
		if (pm_a[j] > pm_a[i])
		{
			x=pm_a[j];	pm_a[j]=pm_a[i];	pm_a[i]=x;
			
			for (l=1; l<=loci; l++)
			for (h=1; h<=2; h++) 
			{
				a=genotype[j][l][h];	genotype[j][l][h]=genotype[i][l][h];	genotype[i][l][h]=a;
			}
			for (l=1; l<=lociQ; l++)
			for (h=1; h<=2; h++)
			{
				a=genotypeQ[j][l][h];	genotypeQ[j][l][h]=genotypeQ[i][l][h];	genotypeQ[i][l][h]=a;	
			}
		}
	}	
}


/* ************************************************************************* */

mating()
{
	// MATING NEUTRAL GENES
	
	for (i=1; i<=ind_sel; i++)
	{
		for (l=1; l<=loci; l++)
		{
			for (h=1; h<=2; h++)		
			{
				parent[i][l][h] = genotype[i][l][h];
			}
		}
		
		for (l=1; l<=lociQ; l++)
		{
			for (h=1; h<=2; h++)		
			{
				parentQ[i][l][h] = genotypeQ[i][l][h];
			}
		}
	}
	
	for	(j=1; j<=ind; j++)
	{
		p1 = (int)(uniform()*ind_sel)+1;
		do
		{
			p2 = (int)(uniform()*ind_sel)+1;
		}
		while (p2 == p1);
		
		
		for (l=1; l<=loci; l++)
		{
			h1 = (int)(uniform()*2)+1;
			h2 = (int)(uniform()*2)+1;

			genotype[j][l][1] = parent[p1][l][h1];
			genotype[j][l][2] = parent[p2][l][h2];
		}
			
		for (l=1; l<=lociQ; l++)
		{
			h1 = (int)(uniform()*2)+1;
			h2 = (int)(uniform()*2)+1;

			genotypeQ[j][l][1] = parentQ[p1][l][h1];
			genotypeQ[j][l][2] = parentQ[p2][l][h2];
			
		}
	}
	
}

/* ************************************************************************* */







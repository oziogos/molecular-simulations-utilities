
#include"topo_diff.h"

void topo(int atoms,char **species,int bonds,int *init_B1,int *init_B2,char **init_B_type,
          char ***global_species,char ****global_bond_type_array,char ****global_angle_type_array,char ****global_dihedral_type_array,char ****global_improper_type_array,
          int **global_neighbors,int ***global_bonds_array,int ***global_angles_array,int ***global_dihedrals_array,int ***global_impropers_array,
          int *angles_core,int *dihedrals_core,int *impropers_core,
          int *atom_types,int *bond_types,int *angle_types,int *dihedral_types,int *improper_types,
          int output_flag,char *current_folder,char *file_name);

char **general_species_registry, ***general_bond_types_registry, ***general_angle_types_registry, ***general_dihedral_types_registry, ***general_improper_types_registry;
int *general_neighbors_registry, **bonds_registry, **general_angles_registry, **general_dihedrals_registry, **general_impropers_registry;
int general_angles, general_dihedrals, general_impropers;
int general_atom_types, general_bond_types, general_angle_types, general_dihedral_types, general_improper_types;

int main(int argc,char **argv)
{
	FILE *fp;
	char current_folder[cmax_length],file_path[cmax_length],buffer[cmax_length],word[cmax_length];
	int atoms,i,j,k,int_buffer,bonds;
	double *x1,*y1,*z1;
	char **species;
	int *B1,*B2;
	char **B_type;
	double *Bval_1,*Aval_1,*Dval_1,r2;
	int ii,jj,kk;
	double xi,yi,zi,xj,yj,zj,xk,yk,zk,rji2,rjk2,rji,rjk,dot,arg,theta;
	int ll;
	double xl,yl,zl,cross,norm1,norm2,sign,coeff,phi,S;
	double *x2,*y2,*z2;
	double *Bval_2,*Aval_2,*Dval_2;
	double val;
	double *Bscore,*Ascore,*Dscore;
	double *B_based_scale,Bmax,rvdw;
	double *A_based_scale,Amax;
	double *D_based_scale,Dmax;
	
	//
	
	getcwd(current_folder,cmax_length);
	
	// read reference xyz file
	sprintf(file_path,"%s/%s",current_folder,argv[1]);
	fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
	fgets(buffer,cmax_length,fp);sscanf(buffer,"%d",&atoms);
	fgets(buffer,cmax_length,fp);
	x1=(double*)malloc(atoms*sizeof(double));
	y1=(double*)malloc(atoms*sizeof(double));
	z1=(double*)malloc(atoms*sizeof(double));
	species=(char**)malloc(atoms*sizeof(char*));for(i=0;i<atoms;++i)species[i]=(char*)malloc(sub_length*sizeof(char));
	for(i=0;i<atoms;++i){fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf\t%lf\t%lf",species[i],&x1[i],&y1[i],&z1[i]);}
	fclose(fp);
	
	// read reference mol2 file for bond topology only!
	sprintf(file_path,"%s/%s",current_folder,argv[2]);
	fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
	fgets(buffer,cmax_length,fp);
	fgets(buffer,cmax_length,fp);
	fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d",&int_buffer,&bonds);
	while(fgets(buffer,cmax_length,fp)!=NULL)if(strcmp(buffer,"@<TRIPOS>BOND\n")==0)break;
	B1=(int*)malloc(bonds*sizeof(int));
	B2=(int*)malloc(bonds*sizeof(int));
	B_type=(char**)malloc(bonds*sizeof(char*));for(i=0;i<bonds;++i)B_type[i]=(char*)malloc(sub_length*sizeof(char));
	for(i=0;i<bonds;++i){fgets(buffer,cmax_length,fp);sscanf(buffer,"%d\t%d\t%d\%s",&int_buffer,&B1[i],&B2[i],B_type[i]);}
	fclose(fp);
	
	// topo call
	topo(atoms,species,bonds,B1,B2,B_type,
         &general_species_registry,&general_bond_types_registry,&general_angle_types_registry,&general_dihedral_types_registry,&general_improper_types_registry,
         &general_neighbors_registry,&bonds_registry,&general_angles_registry,&general_dihedrals_registry,&general_impropers_registry,
         &general_angles,&general_dihedrals,&general_impropers,
         &general_atom_types,&general_bond_types,&general_angle_types,&general_dihedral_types,&general_improper_types,
         1,current_folder,"topo_0_0.txt");
	
	// geometric calculations
	Bval_1=(double*)malloc(bonds*sizeof(double));
	Aval_1=(double*)malloc(general_angles*sizeof(double));
	Dval_1=(double*)malloc(general_dihedrals*sizeof(double));
	// bonds
	for(i=0;i<bonds;++i)
	{
		r2=(x1[B2[i]-1]-x1[B1[i]-1])*(x1[B2[i]-1]-x1[B1[i]-1])+(y1[B2[i]-1]-y1[B1[i]-1])*(y1[B2[i]-1]-y1[B1[i]-1])+(z1[B2[i]-1]-z1[B1[i]-1])*(z1[B2[i]-1]-z1[B1[i]-1]);
		Bval_1[i]=sqrt(r2);
		//printf("%lf\n",Bval_1[i]);
	}
	// angles
	for(i=0;i<general_angles;++i)
	{
		ii=general_angles_registry[i][1]-1;jj=general_angles_registry[i][2]-1;kk=general_angles_registry[i][3]-1;
        xi=x1[ii];yi=y1[ii];zi=z1[ii];
        xj=x1[jj];yj=y1[jj];zj=z1[jj];
        xk=x1[kk];yk=y1[kk];zk=z1[kk];
        rji2=(xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj);
        rjk2=(xk-xj)*(xk-xj)+(yk-yj)*(yk-yj)+(zk-zj)*(zk-zj);
        rji=sqrt(rji2);
        rjk=sqrt(rjk2);
        dot=(xi-xj)*(xk-xj)+(yi-yj)*(yk-yj)+(zi-zj)*(zk-zj);
        arg=dot/(rji*rjk);
        if(arg>1.0)arg=1.0;
        if(arg<-1.0)arg=-1.0;
        theta=acos(arg)*180.0/pi;
        Aval_1[i]=theta;
        //printf("%lf\n",Aval_1[i]);
	}
	// dihedrals
	for(i=0;i<general_dihedrals;++i)
	{
		ii=general_dihedrals_registry[i][1]-1;jj=general_dihedrals_registry[i][2]-1;kk=general_dihedrals_registry[i][3]-1;ll=general_dihedrals_registry[i][4]-1;
        xi=x1[ii];yi=y1[ii];zi=z1[ii];
        xj=x1[jj];yj=y1[jj];zj=z1[jj];
        xk=x1[kk];yk=y1[kk];zk=z1[kk];
        xl=x1[ll];yl=y1[ll];zl=z1[ll];
        cross=(-xj*yi + xk*yi + xi*yj - xk*yj - xi*yk + xj*yk)*(-xk*yj + xl*yj + xj*yk - xl*yk - xj*yl + xk*yl) + (xj*zi - xk*zi - xi*zj + xk*zj + xi*zk - xj*zk)*(xk*zj - xl*zj - xj*zk + xl*zk + xj*zl - xk*zl) + (-yj*zi + yk*zi + yi*zj - yk*zj - yi*zk + yj*zk)*(-yk*zj + yl*zj + yj*zk - yl*zk - yj*zl + yk*zl);
        norm1=sqrt((-xj*yi + xk*yi + xi*yj - xk*yj - xi*yk + xj*yk)*(-xj*yi + xk*yi + xi*yj - xk*yj - xi*yk + xj*yk)+(xj*zi - xk*zi - xi*zj + xk*zj + xi*zk - xj*zk)*(xj*zi - xk*zi - xi*zj + xk*zj + xi*zk - xj*zk)+(-yj*zi + yk*zi + yi*zj - yk*zj - yi*zk + yj*zk)*(-yj*zi + yk*zi + yi*zj - yk*zj - yi*zk + yj*zk));
        norm2=sqrt((-xk*yj + xl*yj + xj*yk - xl*yk - xj*yl + xk*yl)*(-xk*yj + xl*yj + xj*yk - xl*yk - xj*yl + xk*yl)+(xk*zj - xl*zj - xj*zk + xl*zk + xj*zl - xk*zl)*(xk*zj - xl*zj - xj*zk + xl*zk + xj*zl - xk*zl)+(-yk*zj + yl*zj + yj*zk - yl*zk - yj*zl + yk*zl)*(-yk*zj + yl*zj + yj*zk - yl*zk - yj*zl + yk*zl));
        sign=(-xk*yj + xl*yj + xj*yk - xl*yk - xj*yl + xk*yl)*(-zi + zj) + (-yi + yj)*(xk*zj - xl*zj - xj*zk + xl*zk + xj*zl - xk*zl) + (-xi + xj)*(-yk*zj + yl*zj + yj*zk - yl*zk - yj*zl + yk*zl);
        coeff=cross/(norm1*norm2);
        if(coeff>1.0)coeff=1.0;
        if(coeff<-1.0)coeff=-1.0;
        phi=0.0;
        S=1.0;
        if(sign>0.0){
            phi=acos(coeff)*180.0/pi;S=1.0;}
        else if (sign<0.0){phi=-acos(coeff)*180.0/pi;S=-1.0;}
        Dval_1[i]=phi;
        //printf("%lf\n",Dval_1[i]);
	}
	
	// read second xyz file
	sprintf(file_path,"%s/%s",current_folder,argv[3]);
	fp=fopen(file_path,"r");if(fp==NULL){printf("Could not locate %s\n",file_path);exit(-1);}
	fgets(buffer,cmax_length,fp);
	fgets(buffer,cmax_length,fp);
	x2=(double*)malloc(atoms*sizeof(double));
	y2=(double*)malloc(atoms*sizeof(double));
	z2=(double*)malloc(atoms*sizeof(double));
	for(i=0;i<atoms;++i){fgets(buffer,cmax_length,fp);sscanf(buffer,"%s\t%lf\t%lf\t%lf",word,&x2[i],&y2[i],&z2[i]);}
	fclose(fp);
	
	// geometric calculations
	Bval_2=(double*)malloc(bonds*sizeof(double));
	Aval_2=(double*)malloc(general_angles*sizeof(double));
	Dval_2=(double*)malloc(general_dihedrals*sizeof(double));
	// bonds
	for(i=0;i<bonds;++i)
	{
		r2=(x2[B2[i]-1]-x2[B1[i]-1])*(x2[B2[i]-1]-x2[B1[i]-1])+(y2[B2[i]-1]-y2[B1[i]-1])*(y2[B2[i]-1]-y2[B1[i]-1])+(z2[B2[i]-1]-z2[B1[i]-1])*(z2[B2[i]-1]-z2[B1[i]-1]);
		Bval_2[i]=sqrt(r2);
	}
	// angles
	for(i=0;i<general_angles;++i)
	{
		ii=general_angles_registry[i][1]-1;jj=general_angles_registry[i][2]-1;kk=general_angles_registry[i][3]-1;
        xi=x2[ii];yi=y2[ii];zi=z2[ii];
        xj=x2[jj];yj=y2[jj];zj=z2[jj];
        xk=x2[kk];yk=y2[kk];zk=z2[kk];
        rji2=(xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj);
        rjk2=(xk-xj)*(xk-xj)+(yk-yj)*(yk-yj)+(zk-zj)*(zk-zj);
        rji=sqrt(rji2);
        rjk=sqrt(rjk2);
        dot=(xi-xj)*(xk-xj)+(yi-yj)*(yk-yj)+(zi-zj)*(zk-zj);
        arg=dot/(rji*rjk);
        if(arg>1.0)arg=1.0;
        if(arg<-1.0)arg=-1.0;
        theta=acos(arg)*180.0/pi;
        Aval_2[i]=theta;
	}
	// dihedrals
	for(i=0;i<general_dihedrals;++i)
	{
		ii=general_dihedrals_registry[i][1]-1;jj=general_dihedrals_registry[i][2]-1;kk=general_dihedrals_registry[i][3]-1;ll=general_dihedrals_registry[i][4]-1;
        xi=x2[ii];yi=y2[ii];zi=z2[ii];
        xj=x2[jj];yj=y2[jj];zj=z2[jj];
        xk=x2[kk];yk=y2[kk];zk=z2[kk];
        xl=x2[ll];yl=y2[ll];zl=z2[ll];
        cross=(-xj*yi + xk*yi + xi*yj - xk*yj - xi*yk + xj*yk)*(-xk*yj + xl*yj + xj*yk - xl*yk - xj*yl + xk*yl) + (xj*zi - xk*zi - xi*zj + xk*zj + xi*zk - xj*zk)*(xk*zj - xl*zj - xj*zk + xl*zk + xj*zl - xk*zl) + (-yj*zi + yk*zi + yi*zj - yk*zj - yi*zk + yj*zk)*(-yk*zj + yl*zj + yj*zk - yl*zk - yj*zl + yk*zl);
        norm1=sqrt((-xj*yi + xk*yi + xi*yj - xk*yj - xi*yk + xj*yk)*(-xj*yi + xk*yi + xi*yj - xk*yj - xi*yk + xj*yk)+(xj*zi - xk*zi - xi*zj + xk*zj + xi*zk - xj*zk)*(xj*zi - xk*zi - xi*zj + xk*zj + xi*zk - xj*zk)+(-yj*zi + yk*zi + yi*zj - yk*zj - yi*zk + yj*zk)*(-yj*zi + yk*zi + yi*zj - yk*zj - yi*zk + yj*zk));
        norm2=sqrt((-xk*yj + xl*yj + xj*yk - xl*yk - xj*yl + xk*yl)*(-xk*yj + xl*yj + xj*yk - xl*yk - xj*yl + xk*yl)+(xk*zj - xl*zj - xj*zk + xl*zk + xj*zl - xk*zl)*(xk*zj - xl*zj - xj*zk + xl*zk + xj*zl - xk*zl)+(-yk*zj + yl*zj + yj*zk - yl*zk - yj*zl + yk*zl)*(-yk*zj + yl*zj + yj*zk - yl*zk - yj*zl + yk*zl));
        sign=(-xk*yj + xl*yj + xj*yk - xl*yk - xj*yl + xk*yl)*(-zi + zj) + (-yi + yj)*(xk*zj - xl*zj - xj*zk + xl*zk + xj*zl - xk*zl) + (-xi + xj)*(-yk*zj + yl*zj + yj*zk - yl*zk - yj*zl + yk*zl);
        coeff=cross/(norm1*norm2);
        if(coeff>1.0)coeff=1.0;
        if(coeff<-1.0)coeff=-1.0;
        phi=0.0;
        S=1.0;
        if(sign>0.0){
            phi=acos(coeff)*180.0/pi;S=1.0;}
        else if (sign<0.0){phi=-acos(coeff)*180.0/pi;S=-1.0;}
        Dval_2[i]=phi;
	}
	
	// assign scores
	sprintf(file_path,"%s/log.dat",current_folder);
	fp=fopen(file_path,"w+");
	Bscore=(double*)malloc(bonds*sizeof(double));
	Ascore=(double*)malloc(general_angles*sizeof(double));
	Dscore=(double*)malloc(general_dihedrals*sizeof(double));
	for(i=0;i<bonds;++i)
	{
		//Bscore[i]=fabs((Bval_2[i]-Bval_1[i])/Bval_1[i]*100.0);
		Bscore[i]=fabs((Bval_2[i]-Bval_1[i]));
		fprintf(fp,"%lf\t\t%lf\t\t%lf\t\t%d\t%d\n",Bval_1[i],Bval_2[i],Bscore[i],bonds_registry[i][1],bonds_registry[i][2]);
	}
	for(i=0;i<general_angles;++i)
	{
		//Ascore[i]=fabs((Aval_2[i]-Aval_1[i])/Aval_1[i]*100.0);
		Ascore[i]=fabs((Aval_2[i]-Aval_1[i]));
		fprintf(fp,"%lf\t\t%lf\t\t%lf\t\t%d\t%d\t%d\n",Aval_1[i],Aval_2[i],Ascore[i],general_angles_registry[i][1],general_angles_registry[i][2],general_angles_registry[i][3]);
	}
	for(i=0;i<general_dihedrals;++i)
	{
		val=Dval_2[i];
		if(Dval_1[i]*Dval_2[i]<0.0)val=-val;
		//Dscore[i]=fabs((val-Dval_1[i])/Dval_1[i]*100.0);
		Dscore[i]=fabs((val-Dval_1[i]));
		fprintf(fp,"%lf\t\t%lf\t\tactual:\t\t%lf\t\t%lf\t\t%d\t%d\t%d\t%d\n",Dval_1[i],val,Dval_2[i],Dscore[i],general_dihedrals_registry[i][1],general_dihedrals_registry[i][2],general_dihedrals_registry[i][3],general_dihedrals_registry[i][4]);
	}
	fclose(fp);
	
	// bonds based scaling
	B_based_scale=(double*)malloc(atoms*sizeof(double));
	for(i=0;i<atoms;++i)B_based_scale[i]=0.0;
	for(i=0;i<bonds;++i)
	{
		for(j=1;j<=2;++j)if(B_based_scale[bonds_registry[i][j]-1]<Bscore[i])B_based_scale[bonds_registry[i][j]-1]=Bscore[i];
	}
	Bmax=0.0;
	for(i=0;i<atoms;++i)if(B_based_scale[i]>Bmax)Bmax=B_based_scale[i];
	//for(i=0;i<atoms;++i)B_based_scale[i]=B_based_scale[i]/Bmax;
	sprintf(file_path,"%s/bond.csv",current_folder);
	fp=fopen(file_path,"w+");
	fprintf(fp,"X,Y,Z,radius,scale\n");
	for(i=0;i<atoms;++i)
	{
		rvdw=1.0;
		if(strcmp(species[i],"H")==0)rvdw=Hvdw;
		if(strcmp(species[i],"C")==0)rvdw=Cvdw;
		if(strcmp(species[i],"N")==0)rvdw=Nvdw;
		fprintf(fp,"%lf,%lf,%lf,%lf,%lf\n",x1[i],y1[i],z1[i],rvdw,B_based_scale[i]);
	}
	fclose(fp);
	
	// angles based scaling
	A_based_scale=(double*)malloc(atoms*sizeof(double));
	for(i=0;i<atoms;++i)A_based_scale[i]=0.0;
	for(i=0;i<general_angles;++i)
	{
		for(j=1;j<=3;++j)if(A_based_scale[general_angles_registry[i][j]-1]<Ascore[i])A_based_scale[general_angles_registry[i][j]-1]=Ascore[i];
	}
	Amax=0.0;
	for(i=0;i<atoms;++i)if(A_based_scale[i]>Amax)Amax=A_based_scale[i];
	//for(i=0;i<atoms;++i)A_based_scale[i]=A_based_scale[i]/Amax;
	sprintf(file_path,"%s/ang.csv",current_folder);
	fp=fopen(file_path,"w+");
	fprintf(fp,"X,Y,Z,radius,scale\n");
	for(i=0;i<atoms;++i)
	{
		rvdw=1.0;
		if(strcmp(species[i],"H")==0)rvdw=Hvdw;
		if(strcmp(species[i],"C")==0)rvdw=Cvdw;
		if(strcmp(species[i],"N")==0)rvdw=Nvdw;
		fprintf(fp,"%lf,%lf,%lf,%lf,%lf\n",x1[i],y1[i],z1[i],rvdw,A_based_scale[i]);
	}
	fclose(fp);
	
	// dihedral based scaling
	D_based_scale=(double*)malloc(atoms*sizeof(double));
	for(i=0;i<atoms;++i)D_based_scale[i]=0.0;
	for(i=0;i<general_dihedrals;++i)
	{
		for(j=1;j<=4;++j)if(D_based_scale[general_dihedrals_registry[i][j]-1]<Dscore[i])D_based_scale[general_dihedrals_registry[i][j]-1]=Dscore[i];
	}
	Dmax=0.0;
	for(i=0;i<atoms;++i)if(D_based_scale[i]>Dmax)Dmax=D_based_scale[i];
	//for(i=0;i<atoms;++i)D_based_scale[i]=D_based_scale[i]/Dmax;
	//for(i=0;i<atoms;++i)printf("%d\t%lf\n",i+1,D_based_scale[i]);
	sprintf(file_path,"%s/dih.csv",current_folder);
	fp=fopen(file_path,"w+");
	fprintf(fp,"X,Y,Z,radius,scale\n");
	for(i=0;i<atoms;++i)
	{
		rvdw=1.0;
		if(strcmp(species[i],"H")==0)rvdw=Hvdw;
		if(strcmp(species[i],"C")==0)rvdw=Cvdw;
		if(strcmp(species[i],"N")==0)rvdw=Nvdw;
		fprintf(fp,"%lf,%lf,%lf,%lf,%lf\n",x1[i],y1[i],z1[i],rvdw,D_based_scale[i]);
	}
	fclose(fp);
	
	// output
	printf("bonds max score:\t%lf\n",Bmax);
	printf("angles max score:\t%lf\n",Amax);
	printf("dihedrals max score:\t%lf\n",Dmax);
	
	// free
	
	for(i=0;i<atoms;++i)free(species[i]);free(species);free(x1);free(y1);free(z1);
	for(i=0;i<bonds;++i)free(B_type[i]);free(B_type);free(B1);free(B2);
	
	for(j=0;j<bonds;++j)free(bonds_registry[j]);free(bonds_registry);
    for(j=0;j<general_angles;++j)free(general_angles_registry[j]);free(general_angles_registry);
    for(j=0;j<general_dihedrals;++j)free(general_dihedrals_registry[j]);free(general_dihedrals_registry);
    if(general_impropers>0){for(j=0;j<general_impropers;++j)free(general_impropers_registry[j]);free(general_impropers_registry);}
    for(k=0;k<general_atom_types;++k)free(general_species_registry[k]);free(general_species_registry);
    for(k=0;k<general_bond_types;++k)
        for(j=0;j<3;++j)
            free(general_bond_types_registry[k][j]);
    for(k=0;k<general_bond_types;++k)free(general_bond_types_registry[k]);
    free(general_bond_types_registry);
    for(k=0;k<general_angle_types;++k)
        for(j=0;j<5;++j)
            free(general_angle_types_registry[k][j]);
    for(k=0;k<general_angle_types;++k)free(general_angle_types_registry[k]);
    free(general_angle_types_registry);
    for(k=0;k<general_dihedral_types;++k)
        for(j=0;j<7;++j)
            free(general_dihedral_types_registry[k][j]);
    for(k=0;k<general_dihedral_types;++k)free(general_dihedral_types_registry[k]);
    free(general_dihedral_types_registry);
    free(general_neighbors_registry);
    if(general_impropers>0)
    {
        for(k=0;k<general_improper_types;++k)
            for(j=0;j<4;++j)
                free(general_improper_types_registry[k][j]);
        for(k=0;k<general_improper_types;++k)free(general_improper_types_registry[k]);
        free(general_improper_types_registry);
    }
    
    free(Bval_1);
    free(Aval_1);
    free(Dval_1);
    
    free(x2);
    free(y2);
    free(z2);
    
    free(Bval_2);
    free(Aval_2);
    free(Dval_2);
    
    free(Bscore);
    free(Ascore);
    free(Dscore);
    
    free(B_based_scale);
    free(A_based_scale);
    free(D_based_scale);
    
    //
	
	return 0;
}

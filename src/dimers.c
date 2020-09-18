#define cmax_length 1000
#define sub_length 10
#define Hmass 1.00794
#define Cmass 12.0107
#define Nmass 14.0067
#define Omass 15.999
#define Smass 32.065
#define Fmass 18.9984
#define Clmass 35.453
#define Brmass 79.904
#define Imass 126.90447
#define debug 1
#define equal_diff 1e-5
#define CoM_jmol "Xx"
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>
#include<string.h>
int equal(double A,double B);
void resolve_species(int atomic_Z,char *species);
int main(int argc,char **argv)
{
    FILE *fp,*fp_dimer_only;
    char current_folder[cmax_length],file_path[cmax_length],buffer[cmax_length],mol_name[cmax_length],species[sub_length];
    char **atom_name,**atom_type,**subst_name,**bond_type;
    int i,j,k,found,success,num_atoms,num_bonds,num_subst,Rentries,dobreak;
    int *atom_id,*subst_id,*bond_id,*origin_atom_id,*target_atom_id,*atomic_Z,*dimer_i,*dimer_j;
    double cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma,R,cutoff;
    double *x,*y,*z,*mass,*X,*Y,*Z,*M,*Rarray;

    // read arguments
    cutoff=atof(argv[2]);
    
    // open mol2
    getcwd(current_folder,cmax_length);
    sprintf(file_path,"%s/%s",current_folder,argv[1]);
    fp=fopen(file_path,"r");if(fp==NULL){printf("! ... Could not locate %s ... Exiting ... !\n",file_path);exit(-1);}

    // locate and read @<TRIPOS>MOLECULE info
    success=0;
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        if(strcmp(buffer,"@<TRIPOS>MOLECULE\n")==0){success=1;printf("$ ... Located @<TRIPOS>MOLECULE ...\n");break;}
    }
    if(success==0){printf("! ... Unable to locate @<TRIPOS>MOLECULE ... Exiting ... !\n");exit(-1);}
    fgets(buffer,cmax_length,fp);success=sscanf(buffer,"%s",mol_name);printf("$ ... Molecule name = %s ...\n",mol_name);if(success!=1){printf("! ... Could not resolve system name ... Exiting ... !\n");exit(-1);}
    fgets(buffer,cmax_length,fp);success=sscanf(buffer,"%d\t%d\t%d",&num_atoms,&num_bonds,&num_subst);if(success!=3){printf("! ... Could not resolve atoms/bonds/molecules information ... Exiting ... !\n");exit(-1);}printf("$ ... num_atoms/num_bonds/num_subst = %d/%d/%d ...\n",num_atoms,num_bonds,num_subst);
    // preallocate arrays
    atom_name=(char**)malloc(num_atoms*sizeof(char*));for(i=0;i<num_atoms;++i)atom_name[i]=(char*)malloc(sub_length*sizeof(char));
    atom_type=(char**)malloc(num_atoms*sizeof(char*));for(i=0;i<num_atoms;++i)atom_type[i]=(char*)malloc(sub_length*sizeof(char));
    subst_name=(char**)malloc(num_atoms*sizeof(char*));for(i=0;i<num_atoms;++i)subst_name[i]=(char*)malloc(sub_length*sizeof(char));
    bond_type=(char**)malloc(num_bonds*sizeof(char*));for(i=0;i<num_bonds;++i)bond_type[i]=(char*)malloc(sub_length*sizeof(char));
    atom_id=(int*)malloc(num_atoms*sizeof(int));
    subst_id=(int*)malloc(num_atoms*sizeof(int));
    bond_id=(int*)malloc(num_bonds*sizeof(int));
    origin_atom_id=(int*)malloc(num_bonds*sizeof(int));
    target_atom_id=(int*)malloc(num_bonds*sizeof(int));
    atomic_Z=(int*)malloc(num_atoms*sizeof(int));
    x=(double*)malloc(num_atoms*sizeof(double));
    y=(double*)malloc(num_atoms*sizeof(double));
    z=(double*)malloc(num_atoms*sizeof(double));
    mass=(double*)malloc(num_atoms*sizeof(double));
    rewind(fp);
    
    // locate and read @<TRIPOS>ATOM info
    success=0;
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        if(strcmp(buffer,"@<TRIPOS>ATOM\n")==0){success=1;printf("$ ... Located @<TRIPOS>ATOM ...\n");break;}
    }
    if(success==0){printf("! ... Unable to locate @<TRIPOS>ATOM ... Exiting ... !\n");exit(-1);}
    // read atomic info
    for(i=0;i<num_atoms;++i)
    {
        fgets(buffer,cmax_length,fp);success=sscanf(buffer,"%d\t%s\t%lf\t%lf\t%lf\%s\t%d\t%s",&atom_id[i],atom_name[i],&x[i],&y[i],&z[i],atom_type[i],&subst_id[i],subst_name[i]);if(success!=8){printf("! ... Error reading atomic entry %d: %s! Exiting ... !\n",i+1,buffer);exit(-1);}
    }
    rewind(fp);
    
    // locate and read @<TRIPOS>BOND info
    success=0;
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        if(strcmp(buffer,"@<TRIPOS>BOND\n")==0){success=1;printf("$ ... Located @<TRIPOS>BOND ...\n");break;}
    }
    if(success==0){printf("! ... Unable to locate @<TRIPOS>BOND ... Exiting ... !\n");exit(-1);}
    // read bonding info
    for(i=0;i<num_bonds;++i)
    {
        fgets(buffer,cmax_length,fp);success=sscanf(buffer,"%d\t%d\t%d\t%s",&bond_id[i],&origin_atom_id[i],&target_atom_id[i],bond_type[i]);
        if(success!=4){printf("! ... Error reading bond entry %d: %s! Exiting ... !\n",i+1,buffer);exit(-1);}
    }
    rewind(fp);

    // locate and read @<TRIPOS>CRYSIN info
    success=0;
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        if(strcmp(buffer,"@<TRIPOS>CRYSIN\n")==0){success=1;printf("$ ... Located @<TRIPOS>CRYSIN ...\n");break;}
    }
    if(success==0){printf("! ... Unable to locate @<TRIPOS>CRYSIN ... Exiting ... !\n");exit(-1);}
    // read CRYSIN
    fgets(buffer,cmax_length,fp);sscanf(buffer,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf",&cell_a,&cell_b,&cell_c,&cell_alpha,&cell_beta,&cell_gamma);
    printf("$ ... Cell info ...\n$ ... [ a b c ] = [ %lf %lf %lf ] ...\n$ ... [ alpha beta gamma ] = [ %lf %lf %lf ] ...\n",cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma);
    
    fclose(fp);
    
    // resolve species
    for(i=0;i<num_atoms;++i)
    {
        atomic_Z[i]=0;
        // H
        if(strcmp(atom_type[i],"H")==0 || strcmp(atom_type[i],"H.spc")==0 || strcmp(atom_type[i],"H.t3p")==0){atomic_Z[i]=1;mass[i]=Hmass;}
        // C
        if(strcmp(atom_type[i],"C.3")==0 || strcmp(atom_type[i],"C.2")==0 || strcmp(atom_type[i],"C.1")==0 || strcmp(atom_type[i],"C.ar")==0 || strcmp(atom_type[i],"C.cat")==0){atomic_Z[i]=6;mass[i]=Cmass;}
        // N
        if(strcmp(atom_type[i],"N.3")==0 || strcmp(atom_type[i],"N.2")==0 || strcmp(atom_type[i],"N.1")==0 || strcmp(atom_type[i],"N.ar")==0 || strcmp(atom_type[i],"N.am")==0 || strcmp(atom_type[i],"N.pl3")==0 || strcmp(atom_type[i],"N.4")==0){atomic_Z[i]=7;mass[i]=Nmass;}
        // O
        if(strcmp(atom_type[i],"O.3")==0 || strcmp(atom_type[i],"O.2")==0 || strcmp(atom_type[i],"O.co2")==0 || strcmp(atom_type[i],"O.spc")==0 || strcmp(atom_type[i],"O.t3p")==0){atomic_Z[i]=8;mass[i]=Omass;}
        // S
        if(strcmp(atom_type[i],"S.3")==0 || strcmp(atom_type[i],"S.2")==0 || strcmp(atom_type[i],"S.O")==0 || strcmp(atom_type[i],"S.O2")==0){atomic_Z[i]=16;mass[i]=Smass;}
        // F
        if(strcmp(atom_type[i],"F")==0){atomic_Z[i]=9;mass[i]=Fmass;}
        // Cl
        if(strcmp(atom_type[i],"Cl")==0){atomic_Z[i]=17;mass[i]=Clmass;}
        // Br
        if(strcmp(atom_type[i],"Br")==0){atomic_Z[i]=35;mass[i]=Brmass;}
        // I
        if(strcmp(atom_type[i],"I")==0){atomic_Z[i]=53;mass[i]=Imass;}
        if(atomic_Z[i]==0){printf("! ... Unsupported atom type %s ... Exiting ... !\n",atom_type[i]);exit(-1);}
    }
    printf("$ ... Species and masses (amu) resolved ...\n");

    // calculate CoMs
    X=(double*)malloc(num_subst*sizeof(double));
    Y=(double*)malloc(num_subst*sizeof(double));
    Z=(double*)malloc(num_subst*sizeof(double));
    M=(double*)malloc(num_subst*sizeof(double));
    // subst_id sanity check
    for(i=0;i<num_atoms;++i){success=0;for(j=0;j<num_subst;++j){if(subst_id[i]==j+1){success=1;break;}}if(success==0){printf("! ... subst_id sanity check failed ... Exiting ... !\n");exit(-1);}}
    for(i=0;i<num_subst;++i){X[i]=0.0;Y[i]=0.0;Z[i]=0.0;M[i]=0.0;}
    for(i=0;i<num_atoms;++i)
    {
        X[subst_id[i]-1]=X[subst_id[i]-1]+x[i]*mass[i];
        Y[subst_id[i]-1]=Y[subst_id[i]-1]+y[i]*mass[i];
        Z[subst_id[i]-1]=Z[subst_id[i]-1]+z[i]*mass[i];
        M[subst_id[i]-1]=M[subst_id[i]-1]+mass[i];
    }
    for(i=0;i<num_subst;++i){X[i]=X[i]/M[i];Y[i]=Y[i]/M[i];Z[i]=Z[i]/M[i];}
    // conditional preview
    if(debug==1)
    {
        sprintf(file_path,"%s/%s_PREVIEW.xyz",current_folder,mol_name);
        fp=fopen(file_path,"w+");
        fprintf(fp,"%d\nALL:%s\n",num_atoms,mol_name);
        for(i=0;i<num_atoms;++i)
        {
            resolve_species(atomic_Z[i],species);
            fprintf(fp,"%s\t%lf\t%lf\t%lf\n",species,x[i],y[i],z[i]);
        }
        fprintf(fp,"%d\nCoMs:%s\n",num_subst,mol_name);
        for(i=0;i<num_subst;++i)fprintf(fp,"%s\t%lf\t%lf\t%lf\n",CoM_jmol,X[i],Y[i],Z[i]);
        fclose(fp);
    }
    printf("$ ... Molecular CoMs: Done ...\n");
    
    // LOCATE DIMERS
    // populate: conditional bootstrap
    Rentries=1;
    Rarray=(double*)malloc(Rentries*sizeof(double));
    dimer_i=(int*)malloc(Rentries*sizeof(int));
    dimer_j=(int*)malloc(Rentries*sizeof(int));
    success=0;
    for(i=0;i<num_subst-1;++i)
    {
        for(j=i+1;j<num_subst;++j)
        {
            R=(X[j]-X[i])*(X[j]-X[i])+(Y[j]-Y[i])*(Y[j]-Y[i])+(Z[j]-Z[i])*(Z[j]-Z[i]);R=sqrt(R);
            if(R<cutoff)
            {
                Rarray[0]=R;
                dimer_i[0]=i;
                dimer_j[0]=j;
                success=1;dobreak=1;break;
            }
        }
        if(dobreak==1)break;
    }
    if(success==0){printf("! ... cutoff %lf INADEQUATE to find a dimer - revise ... Exiting ... !\n",cutoff);exit(-1);}
    // populate rest
    for(i=0;i<num_subst-1;++i)
    {
        for(j=i+1;j<num_subst;++j)
        {
            //
            R=(X[j]-X[i])*(X[j]-X[i])+(Y[j]-Y[i])*(Y[j]-Y[i])+(Z[j]-Z[i])*(Z[j]-Z[i]);R=sqrt(R);
            if(R<cutoff)
            {
                found=0;
                for(k=0;k<Rentries;++k)
                {
                    if(equal(Rarray[k],R)==1){found=1;break;}
                }
                if(found==0)
                {
                    ++Rentries;
                    Rarray=(double*)realloc(Rarray,Rentries*sizeof(double));
                    dimer_i=(int*)realloc(dimer_i,Rentries*sizeof(int));
                    dimer_j=(int*)realloc(dimer_j,Rentries*sizeof(int));
                    Rarray[Rentries-1]=R;
                    dimer_i[Rentries-1]=i;
                    dimer_j[Rentries-1]=j;
                }
            }
            //
        }
    }
    //
    printf("$ ... Located %d dimers - CoM distance based ONLY! ...\n",Rentries);
    for(i=0;i<Rentries;++i)printf("$ ... Dimer %d {%d-%d} -- R=%lf ...\n",i+1,dimer_i[i]+1,dimer_j[i]+1,Rarray[i]);
    
    // print xyz preview file for dimers
    sprintf(file_path,"%s/%s_DIMERS_cutoff_%lf.xyz",current_folder,mol_name,cutoff);
    fp=fopen(file_path,"w+");
    for(i=0;i<Rentries;++i)
    {
        sprintf(file_path,"%s/%s_DIMER_%d_cutoff_%lf.xyz",current_folder,mol_name,i+1,cutoff);
        fp_dimer_only=fopen(file_path,"w+");
        // resolve dimer atoms; currently not using information from @<TRIPOS>SUBSTRUCTURE (note-to-self: revise?)
        k=0;
        for(j=0;j<num_atoms;++j)if(subst_id[j]==dimer_i[i]+1 || subst_id[j]==dimer_j[i]+1)++k;
        //-- print dimer only --------------------------------------------------
        fprintf(fp,"%d\n%s DIMER %d(%d,%d) - R=%lf\n",k,mol_name,i+1,dimer_i[i]+1,dimer_j[i]+1,Rarray[i]);
        fprintf(fp_dimer_only,"%d\n%s DIMER %d(%d,%d) - R=%lf\n",k,mol_name,i+1,dimer_i[i]+1,dimer_j[i]+1,Rarray[i]);
        // molecule "i"
        for(j=0;j<num_atoms;++j)
        {
            if(subst_id[j]==dimer_i[i]+1)
            {
                resolve_species(atomic_Z[j],species);
                fprintf(fp,"%s\t%lf\t%lf\t%lf\n",species,x[j],y[j],z[j]);
                fprintf(fp_dimer_only,"%s\t%lf\t%lf\t%lf\n",species,x[j],y[j],z[j]);
            }
        }
        // molecule "j"
        for(j=0;j<num_atoms;++j)
        {
            if(subst_id[j]==dimer_j[i]+1)
            {
                resolve_species(atomic_Z[j],species);
                fprintf(fp,"%s\t%lf\t%lf\t%lf\n",species,x[j],y[j],z[j]);
                fprintf(fp_dimer_only,"%s\t%lf\t%lf\t%lf\n",species,x[j],y[j],z[j]);
            }
        }
        fclose(fp_dimer_only);
        //-- print dimer and all CoMs ------------------------------------------
        fprintf(fp,"%d\n%s DIMER %d(%d,%d) - R=%lf\n",k+num_subst,mol_name,i+1,dimer_i[i]+1,dimer_j[i]+1,Rarray[i]);
        // molecule "i"
        for(j=0;j<num_atoms;++j)
        {
            if(subst_id[j]==dimer_i[i]+1)
            {
                resolve_species(atomic_Z[j],species);
                fprintf(fp,"%s\t%lf\t%lf\t%lf\n",species,x[j],y[j],z[j]);
            }
        }
        // molecule "j"
        for(j=0;j<num_atoms;++j)
        {
            if(subst_id[j]==dimer_j[i]+1)
            {
                resolve_species(atomic_Z[j],species);
                fprintf(fp,"%s\t%lf\t%lf\t%lf\n",species,x[j],y[j],z[j]);
            }
        }
        // CoMs
        for(j=0;j<num_subst;++j)fprintf(fp,"%s\t%lf\t%lf\t%lf\n",CoM_jmol,X[j],Y[j],Z[j]);
        //-- print dimer, all CoMs and rest ------------------------------------
        fprintf(fp,"%d\n%s DIMER %d(%d,%d) - R=%lf\n",num_atoms+num_subst,mol_name,i+1,dimer_i[i]+1,dimer_j[i]+1,Rarray[i]);
        // molecule "i"
        for(j=0;j<num_atoms;++j)
        {
            if(subst_id[j]==dimer_i[i]+1)
            {
                resolve_species(atomic_Z[j],species);
                fprintf(fp,"%s\t%lf\t%lf\t%lf\n",species,x[j],y[j],z[j]);
            }
        }
        // molecule "j"
        for(j=0;j<num_atoms;++j)
        {
            if(subst_id[j]==dimer_j[i]+1)
            {
                resolve_species(atomic_Z[j],species);
                fprintf(fp,"%s\t%lf\t%lf\t%lf\n",species,x[j],y[j],z[j]);
            }
        }
        // rest
        for(j=0;j<num_atoms;++j)
        {
            if(subst_id[j]!=dimer_i[i]+1 && subst_id[j]!=dimer_j[i]+1)
            {
                sprintf(species,"%s","Li");
                fprintf(fp,"%s\t%lf\t%lf\t%lf\n",species,x[j],y[j],z[j]);
            }
        }
        // CoMs
        for(j=0;j<num_subst;++j)fprintf(fp,"%s\t%lf\t%lf\t%lf\n",CoM_jmol,X[j],Y[j],Z[j]);
    
    }
    fclose(fp);
    
    
    // free
    for(i=0;i<num_atoms;++i){free(atom_name[i]);free(atom_type[i]);free(subst_name[i]);}
    for(i=0;i<num_bonds;++i)free(bond_type[i]);
    free(atom_name);free(atom_type);free(subst_name);free(bond_type);
    free(atom_id);free(subst_id);free(bond_id);free(origin_atom_id);free(target_atom_id);free(atomic_Z);free(dimer_i);free(dimer_j);
    free(x);free(y);free(z);free(mass);free(X);free(Y);free(Z);free(M);free(Rarray);

    return 0;
}

int equal(double A,double B);
int equal(double A,double B)
{
    int value;
    if(fabs(A-B)<equal_diff){value=1;}else{value=0;}
    return value;
}

void resolve_species(int atomic_Z,char *species);
void resolve_species(int atomic_Z,char *species)
{
    if(atomic_Z==1)sprintf(species,"%s","H");
    if(atomic_Z==6)sprintf(species,"%s","C");
    if(atomic_Z==7)sprintf(species,"%s","N");
    if(atomic_Z==8)sprintf(species,"%s","O");
    if(atomic_Z==9)sprintf(species,"%s","F");
    if(atomic_Z==16)sprintf(species,"%s","S");
    if(atomic_Z==17)sprintf(species,"%s","Cl");
    if(atomic_Z==35)sprintf(species,"%s","Br");
    if(atomic_Z==53)sprintf(species,"%s","I");
}

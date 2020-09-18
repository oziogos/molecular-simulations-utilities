#define cmax_length 1000
#define sub_length 10
#define pi 3.14159265359
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>
#include<string.h>
void resolve_atomic_number(int *atomic_Z,char *atom_type);
void resolve_species(int atomic_Z,char *species);
double proj(double Ax,double Ay,double Az,double Bx,double By,double Bz);
void tri_wrap(int atoms,double ax,double ay,double az,double bx,double by,double bz,double cx,double cy,double cz,double *x,double *y,double *z,int *nx,int *ny,int *nz);
int main(int argc,char **argv)
{
    FILE *fp;
    char current_folder[cmax_length],file_path[cmax_length],buffer[cmax_length],mol_name[cmax_length],species[sub_length],word[cmax_length];
    char **atom_name,**atom_type,**subst_name,**bond_type;
    int i,success,num_atoms,num_bonds,num_subst;
    int *atom_id,*subst_id,*bond_id,*origin_atom_id,*target_atom_id,*atomic_Z;
    double cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma;
    double *x,*y,*z;
    
    char **species_final,c1[cmax_length],c2[cmax_length];
    int j,k,l,repeat_a,repeat_b,repeat_c,ai,bi,ci,*nx,*ny,*nz,atoms_final,*root_atom,int_buffer,sub_start,sub_stop;
    double ax,ay,az,bx,by,bz,cx,cy,cz,*x_final,*y_final,*z_final;

    // read arguments: define repetition along the unit cell vectors a, b, and c
    repeat_a=atoi(argv[2]);
    repeat_b=atoi(argv[3]);
    repeat_c=atoi(argv[4]);
    
    //-- Start reading the mol2 input file -------------------------------------
    //--------------------------------------------------------------------------
    
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
    // TRIPOS vectors
    atom_name=(char**)malloc(num_atoms*sizeof(char*));for(i=0;i<num_atoms;++i)atom_name[i]=(char*)malloc(sub_length*sizeof(char));
    atom_type=(char**)malloc(num_atoms*sizeof(char*));for(i=0;i<num_atoms;++i)atom_type[i]=(char*)malloc(sub_length*sizeof(char));
    subst_name=(char**)malloc(num_atoms*sizeof(char*));for(i=0;i<num_atoms;++i)subst_name[i]=(char*)malloc(sub_length*sizeof(char));
    bond_type=(char**)malloc(num_bonds*sizeof(char*));for(i=0;i<num_bonds;++i)bond_type[i]=(char*)malloc(sub_length*sizeof(char));
    atom_id=(int*)malloc(num_atoms*sizeof(int));
    subst_id=(int*)malloc(num_atoms*sizeof(int));
    bond_id=(int*)malloc(num_bonds*sizeof(int));
    origin_atom_id=(int*)malloc(num_bonds*sizeof(int));
    target_atom_id=(int*)malloc(num_bonds*sizeof(int));
    x=(double*)malloc(num_atoms*sizeof(double));
    y=(double*)malloc(num_atoms*sizeof(double));
    z=(double*)malloc(num_atoms*sizeof(double));
    // atomic number vector for species identification
    atomic_Z=(int*)malloc(num_atoms*sizeof(int));
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

    // locate and read @<TRIPOS>SUBSTRUCTURE
    success=0;
    while(fgets(buffer,cmax_length,fp)!=NULL)
    {
        if(strcmp(buffer,"@<TRIPOS>SUBSTRUCTURE\n")==0){success=1;printf("$ ... Located @<TRIPOS>SUBSTRUCTURE ...\n");break;}
    }
    if(success==0){printf("! ... Unable to locate @<TRIPOS>SUBSTRUCTURE ... Exiting ... !\n");exit(-1);}
    // read SUBSTRUCTURE
    root_atom=(int*)malloc(num_subst*sizeof(int));
    for(i=0;i<num_subst;++i)
    {
        fgets(buffer,cmax_length,fp);printf("%s",buffer);sscanf(buffer,"%d\t%s\t%d",&int_buffer,word,&root_atom[i]);
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
    fgets(buffer,cmax_length,fp);sscanf(buffer,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%s\t%s",&cell_a,&cell_b,&cell_c,&cell_alpha,&cell_beta,&cell_gamma,c1,c2);
    printf("$ ... Cell info ...\n$ ... [ a b c ] = [ %lf %lf %lf ] ...\n$ ... [ alpha beta gamma ] = [ %lf %lf %lf ] ...\n",cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma);
    
    fclose(fp);
    
    //-- End of mol2 file reading ----------------------------------------------
    //--------------------------------------------------------------------------
    
    // resolve atomic number
    for(i=0;i<num_atoms;++i)resolve_atomic_number(&atomic_Z[i],atom_type[i]);
    printf("$ ... Species resolved ...\n");
    
    // supercell
    // cast given unit cell according to the geometric construction found at:
    // http://lammps.sandia.gov/doc/Section_howto.html#howto-12 (or just search lammps triclinic cell online...)
    // from {cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma} to {ax,bx,by,cx,cy,cz}
    ax=cell_a;
    ay=0.0;
    az=0.0;
    bx=cell_b*cos(cell_gamma*pi/180.0);
    by=cell_b*sin(cell_gamma*pi/180.0);
    bz=0.0;
    cx=cell_c*cos(cell_beta*pi/180.0);
    cy=(cell_b*cell_c*cos(cell_alpha*pi/180.0)-bx*cx)/by;
    cz=cell_c*cell_c-cx*cx-cy*cy;cz=sqrt(cz);
    
    // atoms after periodic repetition
    atoms_final=num_atoms*repeat_a*repeat_b*repeat_c;
    
    // preallocations
    species_final=(char**)malloc(atoms_final*sizeof(char*));for(i=0;i<atoms_final;++i)species_final[i]=(char*)malloc(sub_length*sizeof(char));
    x_final=(double*)malloc(atoms_final*sizeof(double));
    y_final=(double*)malloc(atoms_final*sizeof(double));
    z_final=(double*)malloc(atoms_final*sizeof(double));
    nx=(int*)malloc(atoms_final*sizeof(int));
    ny=(int*)malloc(atoms_final*sizeof(int));
    nz=(int*)malloc(atoms_final*sizeof(int));
    // populate (unwrapped version); species identification is based on atomic numbers
    j=-1;
    for(ai=0;ai<repeat_a;++ai)
        for(bi=0;bi<repeat_b;++bi)
            for(ci=0;ci<repeat_c;++ci)
                for(i=0;i<num_atoms;++i){++j;
                    resolve_species(atomic_Z[i],species_final[j]);
                    x_final[j]=x[i]+ai*ax+bi*bx+ci*cx;
                    y_final[j]=y[i]+ai*ay+bi*by+ci*cy;
                    z_final[j]=z[i]+ai*az+bi*bz+ci*cz;nx[j]=0;ny[j]=0;nz[j]=0;}
    // update unit cell according to repetitions
    ax=repeat_a*ax;
    bx=repeat_b*bx;by=repeat_b*by;
    cx=repeat_c*cx;cy=repeat_c*cy;cz=repeat_c*cz;
    
    // preview: unwrapped
    sprintf(file_path,"%s/%s_PREVIEW.xyz",current_folder,argv[1]);
    fp=fopen(file_path,"w+");
    fprintf(fp,"%d\n\n",atoms_final+1+3+4);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",0.0,0.0,0.0);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",ax,ay,az);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",bx,by,bz);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",cx,cy,cz);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",ax+bx,ay+by,az+bz);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",bx+cx,by+cy,bz+cz);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",cx+ax,cy+ay,cz+az);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",cx+ax+bx,cy+ay+by,cz+az+bz);
    for(i=0;i<atoms_final;++i)fprintf(fp,"%s\t%lf\t%lf\t%lf\n",species_final[i],x_final[i],y_final[i],z_final[i]);
    
    //-- PERIODIC WRAPPING -----------------------------------------------------
    tri_wrap(atoms_final,ax,ay,az,bx,by,bz,cx,cy,cz,x_final,y_final,z_final,nx,ny,nz);
    
    // preview: wrapped
    fprintf(fp,"%d\n\n",atoms_final+1+3+4);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",0.0,0.0,0.0);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",ax,ay,az);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",bx,by,bz);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",cx,cy,cz);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",ax+bx,ay+by,az+bz);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",bx+cx,by+cy,bz+cz);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",cx+ax,cy+ay,cz+az);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",cx+ax+bx,cy+ay+by,cz+az+bz);
    for(i=0;i<atoms_final;++i)fprintf(fp,"%s\t%lf\t%lf\t%lf\n",species_final[i],x_final[i],y_final[i],z_final[i]);
    
    // preview: unwrapped
    fprintf(fp,"%d\n\n",atoms_final+1+3+4);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",0.0,0.0,0.0);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",ax,ay,az);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",bx,by,bz);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",cx,cy,cz);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",ax+bx,ay+by,az+bz);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",bx+cx,by+cy,bz+cz);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",cx+ax,cy+ay,cz+az);
    fprintf(fp,"Xx\t%lf\t%lf\t%lf\n",cx+ax+bx,cy+ay+by,cz+az+bz);
    for(i=0;i<atoms_final;++i)fprintf(fp,"%s\t%lf\t%lf\t%lf\n",species_final[i],x_final[i]+nx[i]*ax+ny[i]*bx+nz[i]*cx,y_final[i]+nx[i]*ay+ny[i]*by+nz[i]*cy,z_final[i]+nx[i]*az+ny[i]*bz+nz[i]*cz);
    
    fclose(fp);

    //--------------------------------------------------------------------------
    
    // write extended crystal mol2
    sprintf(file_path,"%s/%s_CRYSTAL_%d_%d_%d.mol2",current_folder,argv[1],repeat_a,repeat_b,repeat_c);
    fp=fopen(file_path,"w+");
    // header
    fprintf(fp,"@<TRIPOS>MOLECULE\n");
    fprintf(fp,"%s\n",mol_name);
    fprintf(fp," %d %d %d %d %d\n",atoms_final,num_bonds*repeat_a*repeat_b*repeat_c,num_subst*repeat_a*repeat_b*repeat_c,0,0);
    fprintf(fp,"SMALL\nNO_CHARGES\n****\nGenerated from tri_crystal util\n\n");
    // write atoms
    fprintf(fp,"@<TRIPOS>ATOM\n");
    k=0;    // counter for subst
    j=-1;   // counter for total atoms
    for(ai=0;ai<repeat_a;++ai){
        for(bi=0;bi<repeat_b;++bi){
            for(ci=0;ci<repeat_c;++ci){
                // loop on subst
                for(i=0;i<num_subst-1;++i)
                {
                    ++k;
                    sub_start=root_atom[i]-1;sub_stop=root_atom[i+1]-2;
                    for(l=sub_start;l<=sub_stop;++l)
                    {
                        ++j;
                        fprintf(fp,"%d\t%s\t%lf\t%lf\t%lf\t%s\t%d\tRES%d\t%lf\n",j+1,atom_name[l],
                                x_final[j]+nx[j]*ax+ny[j]*bx+nz[j]*cx,
                                y_final[j]+nx[j]*ay+ny[j]*by+nz[j]*cy,
                                z_final[j]+nx[j]*az+ny[j]*bz+nz[j]*cz,
                                atom_type[l],k,k,0.0);
                    }
                    //printf("%d-%d\n",sub_start,sub_stop);
                }
                ++k;
                sub_start=root_atom[num_subst-1]-1;sub_stop=num_atoms-1;
                for(l=sub_start;l<=sub_stop;++l)
                {
                    ++j;
                    fprintf(fp,"%d\t%s\t%lf\t%lf\t%lf\t%s\t%d\tRES%d\t%lf\n",j+1,atom_name[l],
                            x_final[j]+nx[j]*ax+ny[j]*bx+nz[j]*cx,
                            y_final[j]+nx[j]*ay+ny[j]*by+nz[j]*cy,
                            z_final[j]+nx[j]*az+ny[j]*bz+nz[j]*cz,
                            atom_type[l],k,k,0.0);
                }
                //printf("%d-%d\n",sub_start,sub_stop);
                /*
                for(i=0;i<num_atoms;++i){++j;
                    //      1 ca      2.2420  10.3128   3.8502   C.ar       1 RES1   0.0000
                    fprintf(fp,"%d\t%s\t%lf\t%lf\t%lf\t%s\t%d\tRES%d\t%lf\n",j+1,atom_name[i],
                            x_final[j]+nx[j]*ax+ny[j]*bx+nz[j]*cx,
                            y_final[j]+nx[j]*ay+ny[j]*by+nz[j]*cy,
                            z_final[j]+nx[j]*az+ny[j]*bz+nz[j]*cz,
                            atom_type[i],k,k,0.0);
                }
                */
            }
        }
    }
    // write bonds
    fprintf(fp,"@<TRIPOS>BOND\n");
    k=0;
    for(i=0;i<repeat_a*repeat_b*repeat_c;++i)
        for(j=0;j<num_bonds;++j){
            ++k;
            fprintf(fp,"%d\t%d\t%d\t%s\n",k,origin_atom_id[j]+i*num_atoms,target_atom_id[j]+i*num_atoms,bond_type[j]);
        }
    // write substructure
    fprintf(fp,"@<TRIPOS>SUBSTRUCTURE\n");
    j=0;
    k=-1;
    for(ai=0;ai<repeat_a;++ai){
        for(bi=0;bi<repeat_b;++bi){
            for(ci=0;ci<repeat_c;++ci){
                ++k;
                for(i=0;i<num_subst;++i)
                {
                    ++j;
                    fprintf(fp,"%d\tRES%d\t%d\tGROUP\t0\t****\t****\t0\n",j,j,root_atom[i]+k*num_atoms);
                }
            }}}
    // write cell
    fprintf(fp,"@<TRIPOS>CRYSIN\n");
    fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%s\t%s\n",cell_a*repeat_a,cell_b*repeat_b,cell_c*repeat_c,cell_alpha,cell_beta,cell_gamma,c1,c2);
    fclose(fp);
    
    //--------------------------------------------------------------------------
    
    // free
    for(i=0;i<num_atoms;++i){free(atom_name[i]);free(atom_type[i]);free(subst_name[i]);}
    for(i=0;i<num_bonds;++i)free(bond_type[i]);
    free(atom_name);free(atom_type);free(subst_name);free(bond_type);
    free(atom_id);free(subst_id);free(bond_id);free(origin_atom_id);free(target_atom_id);free(atomic_Z);
    free(x);free(y);free(z);free(root_atom);
    
    for(i=0;i<atoms_final;++i)free(species_final[i]);free(species_final);free(x_final);free(y_final);free(z_final);free(nx);free(ny);free(nz);
    
    return 0;
}

double proj(double Ax,double Ay,double Az,double Bx,double By,double Bz);
double proj(double Ax,double Ay,double Az,double Bx,double By,double Bz)
{
    double dot,norm,res;
    dot=Ax*Bx+Ay*By+Az*Bz;
    norm=Ax*Ax+Ay*Ay+Az*Az;
    res=dot/norm;
    return res;
}

void resolve_atomic_number(int *atomic_Z,char *atom_type);
void resolve_atomic_number(int *atomic_Z,char *atom_type)
{
    *atomic_Z=0;
    // H
    if(strcmp(atom_type,"H")==0 || strcmp(atom_type,"H.spc")==0 || strcmp(atom_type,"H.t3p")==0){*atomic_Z=1;}
    // B
    if(strcmp(atom_type,"B")==0){*atomic_Z=5;}
    // C
    if(strcmp(atom_type,"C.3")==0 || strcmp(atom_type,"C.2")==0 || strcmp(atom_type,"C.1")==0 || strcmp(atom_type,"C.ar")==0 || strcmp(atom_type,"C.cat")==0){*atomic_Z=6;}
    // N
    if(strcmp(atom_type,"N.3")==0 || strcmp(atom_type,"N.2")==0 || strcmp(atom_type,"N.1")==0 || strcmp(atom_type,"N.ar")==0 || strcmp(atom_type,"N.am")==0 || strcmp(atom_type,"N.pl3")==0 || strcmp(atom_type,"N.4")==0){*atomic_Z=7;}
    // O
    if(strcmp(atom_type,"O.3")==0 || strcmp(atom_type,"O.2")==0 || strcmp(atom_type,"O.co2")==0 || strcmp(atom_type,"O.spc")==0 || strcmp(atom_type,"O.t3p")==0){*atomic_Z=8;}
    // S
    if(strcmp(atom_type,"S.3")==0 || strcmp(atom_type,"S.2")==0 || strcmp(atom_type,"S.O")==0 || strcmp(atom_type,"S.O2")==0){*atomic_Z=16;}
    // F
    if(strcmp(atom_type,"F")==0){*atomic_Z=9;}
    // Cl
    if(strcmp(atom_type,"Cl")==0){*atomic_Z=17;}
    // Br
    if(strcmp(atom_type,"Br")==0){*atomic_Z=35;}
    // I
    if(strcmp(atom_type,"I")==0){*atomic_Z=53;}
    if(*atomic_Z==0){printf("! ... Unsupported atom type %s ... Exiting ... !\n",atom_type);exit(-1);}
}

void resolve_species(int atomic_Z,char *species);
void resolve_species(int atomic_Z,char *species)
{
    if(atomic_Z==1)sprintf(species,"%s","H");
    if(atomic_Z==5)sprintf(species,"%s","B");
    if(atomic_Z==6)sprintf(species,"%s","C");
    if(atomic_Z==7)sprintf(species,"%s","N");
    if(atomic_Z==8)sprintf(species,"%s","O");
    if(atomic_Z==9)sprintf(species,"%s","F");
    if(atomic_Z==16)sprintf(species,"%s","S");
    if(atomic_Z==17)sprintf(species,"%s","Cl");
    if(atomic_Z==35)sprintf(species,"%s","Br");
    if(atomic_Z==53)sprintf(species,"%s","I");
}

void tri_wrap(int atoms,double ax,double ay,double az,double bx,double by,double bz,double cx,double cy,double cz,double *x,double *y,double *z,int *nx,int *ny,int *nz);
void tri_wrap(int atoms,double ax,double ay,double az,double bx,double by,double bz,double cx,double cy,double cz,double *x,double *y,double *z,int *nx,int *ny,int *nz)
{
    // Based on:
    // Tuckerman, M. Statistical Mechanics: Theory and Molecular Simulation; Oxford University Press, 2010.
    // Appendix B; p. 655
    
    double V,norm,ax_hat,ay_hat,az_hat,bx_hat,by_hat,bz_hat,cx_hat,cy_hat,cz_hat;
    double a0x_hat=1.0,a0y_hat=0.0,a0z_hat=0.0;
    double b0x_hat=0.0,b0y_hat=1.0,b0z_hat=0.0;
    double c0x_hat=0.0,c0y_hat=0.0,c0z_hat=1.0;
    double projaa0,projab0,projac0,projba0,projbb0,projbc0,projca0,projcb0,projcc0;
    double denom,a_norm_REF,b_norm_REF,c_norm_REF,nom1,nom2,nom3,a,b,c,a_sign,b_sign,c_sign,a_norm,b_norm,c_norm;
    
    int i;
    
    V=ax*by*cz;
    
    // normalize
    norm=ax*ax+ay*ay+az*az;norm=sqrt(norm);
    ax_hat=ax/norm;
    ay_hat=ay/norm;
    az_hat=az/norm;
    norm=bx*bx+by*by+bz*bz;norm=sqrt(norm);
    bx_hat=bx/norm;
    by_hat=by/norm;
    bz_hat=bz/norm;
    norm=cx*cx+cy*cy+cz*cz;norm=sqrt(norm);
    cx_hat=cx/norm;
    cy_hat=cy/norm;
    cz_hat=cz/norm;
    
    // calculate projections
    projaa0=proj(ax_hat,ay_hat,az_hat,a0x_hat,a0y_hat,a0z_hat);
    projab0=proj(ax_hat,ay_hat,az_hat,b0x_hat,b0y_hat,b0z_hat);
    projac0=proj(ax_hat,ay_hat,az_hat,c0x_hat,c0y_hat,c0z_hat);
    projba0=proj(bx_hat,by_hat,bz_hat,a0x_hat,a0y_hat,a0z_hat);
    projbb0=proj(bx_hat,by_hat,bz_hat,b0x_hat,b0y_hat,b0z_hat);
    projbc0=proj(bx_hat,by_hat,bz_hat,c0x_hat,c0y_hat,c0z_hat);
    projca0=proj(cx_hat,cy_hat,cz_hat,a0x_hat,a0y_hat,a0z_hat);
    projcb0=proj(cx_hat,cy_hat,cz_hat,b0x_hat,b0y_hat,b0z_hat);
    projcc0=proj(cx_hat,cy_hat,cz_hat,c0x_hat,c0y_hat,c0z_hat);
    
    // projection denominator
    denom=(-projac0*projbb0*projca0+projab0*projbc0*projca0+projac0*projba0*projcb0-projaa0*projbc0*projcb0-projab0*projba0*projcc0+projaa0*projbb0*projcc0);
    
    // oblique supercell vector lengths
    a_norm_REF=ax*ax;a_norm_REF=sqrt(a_norm_REF);
    b_norm_REF=bx*bx+by*by;b_norm_REF=sqrt(b_norm_REF);
    c_norm_REF=cx*cx+cy*cy+cz*cz;c_norm_REF=sqrt(c_norm_REF);
    
    for(i=0;i<atoms;++i)
    {
        // projections
        nom1=-1.0*(projbc0*projcb0*x[i]-projbb0*projcc0*x[i]-projbc0*projca0*y[i]+projba0*projcc0*y[i]+projbb0*projca0*z[i]-projba0*projcb0*z[i]);
        nom2=(projac0*projcb0*x[i]-projab0*projcc0*x[i]-projac0*projca0*y[i]+projaa0*projcc0*y[i]+projab0*projca0*z[i]-projaa0*projcb0*z[i]);
        nom3=-1.0*(projac0*projbb0*x[i]-projab0*projbc0*x[i]-projac0*projba0*y[i]+projaa0*projbc0*y[i]+projab0*projba0*z[i]-projaa0*projbb0*z[i]);
        
        // coordinates with respect to oblique system
        a=nom1/denom;
        b=nom2/denom;
        c=nom3/denom;
        
        // +/- direction
        if((a*ax_hat-0.0)*(ax-0.0)+(a*ay_hat-0.0)*(0.0-0.0)+(a*az_hat-0.0)*(0.0-0.0) < 0.0 ){a_sign=-1.0;}else{a_sign=+1.0;}
        if((b*bx_hat-0.0)*(bx-0.0)+(b*by_hat-0.0)*(by-0.0)+(b*bz_hat-0.0)*(0.0-0.0) < 0.0 ){b_sign=-1.0;}else{b_sign=+1.0;}
        if((c*cx_hat-0.0)*(cx-0.0)+(c*cy_hat-0.0)*(cy-0.0)+(c*cz_hat-0.0)*(cz-0.0) < 0.0 ){c_sign=-1.0;}else{c_sign=+1.0;}
        
        // image calculation
        a_norm=(a*ax_hat)*(a*ax_hat)+(a*ay_hat)*(a*ay_hat)+(a*az_hat)*(a*az_hat);a_norm=sqrt(a_norm);
        b_norm=(b*bx_hat)*(b*bx_hat)+(b*by_hat)*(b*by_hat)+(b*bz_hat)*(b*bz_hat);b_norm=sqrt(b_norm);
        c_norm=(c*cx_hat)*(c*cx_hat)+(c*cy_hat)*(c*cy_hat)+(c*cz_hat)*(c*cz_hat);c_norm=sqrt(c_norm);
        nx[i]=(int)floor(a_norm/a_norm_REF);if(a_sign<0.0)nx[i]=-nx[i]-1;
        ny[i]=(int)floor(b_norm/b_norm_REF);if(b_sign<0.0)ny[i]=-ny[i]-1;
        nz[i]=(int)floor(c_norm/c_norm_REF);if(c_sign<0.0)nz[i]=-nz[i]-1;
        
        // periodic wrapping
        x[i]=x[i]-nx[i]*ax-ny[i]*bx-nz[i]*cx;
        y[i]=y[i]-nx[i]*(0.0)-ny[i]*by-nz[i]*cy;
        z[i]=z[i]-nx[i]*(0.0)-ny[i]*(0.0)-nz[i]*cz;
        
    }
    
}

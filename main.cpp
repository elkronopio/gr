///Calculo de g(r)

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string.h>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<iomanip>
#include <sstream>

using namespace std;

#define pi 3.14159265358979
#define npartMax 50000
#define MDINLINE static inline
//#define ROTATION
//#define MDLC


int PERIODIC[3]={1,1,1};
char archivo[30];



MDINLINE double dround(double x) { return floor(x+0.5); }


struct particula
{
    double pos[3];
    double vel[3];
    double force[3];
    double diam;
    double mass;
    double dip[3];
    int id;
    
};

struct parametros
{
    double t;
    double lbox[3];
    double lboxh[3];
    double lboxi[3];
    int nacc;
    int npart;
    double vol;
    double gap;
};

particula part[1000];
parametros param;

MDINLINE double distancia(particula p1, particula p2)
{
  int i;
  double dis=0.0;
  double res[3];
  for(i=0;i<3;i++) {
    res[i] = p2.pos[i] - p1.pos[i];
    if (PERIODIC[i])
      res[i] -= dround(res[i]*param.lboxi[i])*param.lbox[i];
    dis+=res[i]*res[i];
  }
  dis=sqrt(dis);
}





int main(int argc, char** argv)
{
    int importBinary (ifstream *filein);
    int nbytes, nbytesTotal;
    double gr[20000];
    int i,j, nl,nlim;
    double dis,dnorm,dr;
    char directorio[200];
    double t1,t2,dt,cada;
    int nconfCiclo,nconf,count;
    int cicloIni, cicloFin, nconfFin;
    double v1,v2;
    double d1,d2;
    
    param.gap=5.0;
    stringstream ss;
    string pp;// =argv[1];
    pp="/home/jcfernandez/simulacion/ROTACIONES/ROT_OFF/DARWAAN/phi0.2/H100/config.bin";
    ss<<pp;
    ss>>directorio;
    
    dr=0.01;
    nlim=0;
//    cada=1;
//    cicloIni=20;
//    cicloFin=25;
    cout<<"dir "<<directorio<<endl;
    ifstream filein(directorio,ios::in|ios::binary);
    nbytes = importBinary (&filein);
    filein.seekg (0,ios::end);
    nbytesTotal =filein.tellg();

    cout<<nbytes<<"  "<<nbytesTotal<<endl;
    nconf=int(nbytesTotal/nbytes);
    count=0;
    int nconfini=nconf-100;
    cout<<"Numerto total de configuraciones: "<<nconf<<endl;
    filein.seekg (nconfini*nbytes,ios::beg);
    while(count<100)
    {
        if(count%10==0) cout<<count<<"  "<<nlim<<endl;
        importBinary (&filein);
        for(i=0;i<param.npart-1;i++)
        {
            for(j=i+1;j<param.npart;j++)
            {
                dis=distancia(part[i],part[j]);

                if (dis <param.lboxh[2])
                {
                     nl=int(floor((dis)/dr));
                     
                     if(nl>nlim) nlim=nl;
                     d1=nl*dr;
                     d2=d1+dr;
                     /////////////////////////
//                     v1=4.0/3.0*pi*d1*d1*d1;
//                     v2=4.0/3.0*pi*d2*d2*d2;
                     ///////////////////////
                     if(d2>part[i].pos[2])
                     {
                        ///Hay corte con la pared inferior
                        v2=2.0/3.0*pi*(1+part[i].pos[2]/d2)*d2*d2*d2;
                       // v2=4.0/3.0*pi*d2*d2*d2;
                         if(d1>part[i].pos[2])
                         {
                             v1=2.0/3.0*pi*(1+part[i].pos[2]/d1)*d1*d1*d1;
                             //v1=4.0/3.0*pi*d1*d1*d1;
                         }
                         else
                         {
                             v1=4.0/3.0*pi*d1*d1*d1;
                         }

                     }
                     else
                     {
                         if(d2>param.lbox[2]-part[i].pos[2])
                         {
                             ///Hay corte con la pared superior
                              v2=2.0/3.0*pi*(1+(param.lbox[2]-part[i].pos[2])/d2)*d2*d2*d2;
                              //cout<<v2/(4.0/3.0*pi*d2*d2*d2)<<endl;
                              //v2=4.0/3.0*pi*d2*d2*d2;
                              if(d1>param.lbox[2]-part[i].pos[2])
                              {
                                  v1=2.0/3.0*pi*(1+(param.lbox[2]-part[i].pos[2])/d1)*d1*d1*d1;
                              //    v1=4.0/3.0*pi*d1*d1*d1;
                              }
                              else
                              {
                                  v1=4.0/3.0*pi*d1*d1*d1;
                              }
                         }
                         else
                         {
                             ///No hay corte con las paredes superior ni inferior
                             //gr[nl]+=2.0/(4.0/3.0*pi*dr*dr*dr*(3.0*nl*nl+3.0*nl+1));
                             v2=4.0/3.0*pi*d2*d2*d2;
                              v1=4.0/3.0*pi*d1*d1*d1;
                         }

                     }
                     gr[nl]+=2.0/(v2-v1);

                    //gr[nl]+=2;
                }
            }
        }
        count++;
    }

    cout<<nlim<<endl;

    ofstream fout("gr.dat",ios::out);
    fout<<dr<<"  "<<0.<<endl;
    for(i=1;i<nlim-1;i++)
    {
        dis=(i+0.5)*dr;
        //v1=4.0/3.0*pi*(i*dr)*(i*dr)*(i*dr);
        //v2=4.0/3.0*pi*((i+1)*dr)*((i+1)*dr)*((i+1)*dr);
       // dnorm=param.vol/param.npart/(4.0/3.0*pi*dr*dr*dr*(3.0*i*i+3.0*i+1))/param.npart/count;
        //dnorm=param.vol/param.npart/(v2-v1)/param.npart/count;
        //dnorm=param.vol/param.npart/param.npart/count;
        dnorm=param.vol/param.npart/param.npart/count;

        fout<<dis<<"  "<<double(gr[i])*dnorm<<endl;
   //fout<<dis<<"  "<<double(gr[i])/param.npart/count<<endl;
      //  cout<<dis<<"  "<<double(igr[i])<<"  "<<dnorm<<"  "<<count<<endl;
    }


    fout.close();

//	nbytes = importBinary (&filein,0,0);
//	importBinary (&filein,nbytes,1);
//	t1=param.t;
//	importBinary (&filein,nbytes,2);
//	t2=param.t;
//	dt=t2-t1;
//	nconfCiclo=int(floor(2.0*pi/0.1/dt)+1);
//    nlim=int(param.lboxh[0]/dr);
//    nconf=cicloIni*nconfCiclo;
//    nconfFin=cicloFin*nconfCiclo;
//    //nconfFin=nconf+100;
//
//    for(i=0;i<20000;i++) igr[i]=0;
//
//
//    count=0;
//    while(!filein.eof() && nconf<nconfFin)
//    {
//        importBinary (&filein,nbytes,nconf);
//        for(i=0;i<param.npart-1;i++)
//        {
//            for(j=i+1;j<param.npart;j++)
//            {
//                dis=distancia(i,j);
//                if (dis <param.lboxh[0])
//                {
//                     nl=int(floor((dis)/dr));
//                     igr[nl]+=2;
//                }
//            }
//        }
//        nconf+=cada;
//        count++;
//    }
//        ofstream fout("gr.dat",ios::out);
//        fout<<0.5<<"  "<<0.<<endl;
//        for(i=1;i<nlim-1;i++)
//        {
//            dis=(i+0.5)*dr;
//            dnorm=param.vol/param.npart/(4.0/3.0*pi*dr*dr*dr*(3.0*i*i+3.0*i+1))/param.npart;
//            fout<<dis<<"  "<<double(igr[i])*dnorm/count<<endl;
//        }
//
//
//    fout.close();
	filein.close();
}


int importBinary (ifstream *filein)
{
    int id;
    int nbytes=40;
    (*filein).read(reinterpret_cast<char *>(&param.t), sizeof(double));
    (*filein).read(reinterpret_cast<char *>(&param.nacc), sizeof(int));
    (*filein).read(reinterpret_cast<char *>(&param.npart), sizeof(int));
    (*filein).read(reinterpret_cast<char *>(&param.lbox[0]), sizeof(double));
    (*filein).read(reinterpret_cast<char *>(&param.lbox[1]), sizeof(double));
#ifndef MDLC
    (*filein).read(reinterpret_cast<char *>(&param.lbox[2]), sizeof(double));
#else
    (*filein).read(reinterpret_cast<char *>(&param.lbox[2]), sizeof(double));
    param.lbox[2]=param.lbox[0]-param.gap;
#endif
    for (int i=0;i<3;i++) 
    {
        param.lboxi[i]=1.0/param.lbox[i];
        param.lboxh[i]=0.5*param.lbox[i];
    }

  //  cout<<param.lbox[0]<<endl;
    for (int j=0;j<param.npart;j++)
    {
        (*filein).read(reinterpret_cast<char *>(&part[j].pos[0]), sizeof(double));
        (*filein).read(reinterpret_cast<char *>(&part[j].pos[1]), sizeof(double));
        (*filein).read(reinterpret_cast<char *>(&part[j].pos[2]), sizeof(double));
        
        
#ifndef ROTATION
        (*filein).read(reinterpret_cast<char *>(&part[j].vel[0]), sizeof(double));
        (*filein).read(reinterpret_cast<char *>(&part[j].vel[1]), sizeof(double));
        (*filein).read(reinterpret_cast<char *>(&part[j].vel[2]), sizeof(double));
        (*filein).read(reinterpret_cast<char *>(&part[j].force[0]), sizeof(double));
        (*filein).read(reinterpret_cast<char *>(&part[j].force[1]), sizeof(double));
        (*filein).read(reinterpret_cast<char *>(&part[j].force[2]), sizeof(double));
        (*filein).read(reinterpret_cast<char *>(&part[j].diam), sizeof(double));
        (*filein).read(reinterpret_cast<char *>(&id), sizeof(int));
        nbytes+=84;
#else
        (*filein).read(reinterpret_cast<char *>(&part[j].dip[0]), sizeof(double));
        (*filein).read(reinterpret_cast<char *>(&part[j].dip[1]), sizeof(double));
        (*filein).read(reinterpret_cast<char *>(&part[j].dip[2]), sizeof(double));
        part[j].diam=1.0;
        (*filein).read(reinterpret_cast<char *>(&id), sizeof(int));
        nbytes+=52;
#endif
        part[j].mass=pow(part[j].diam,3.0);
    }

    param.vol=(param.lbox[0])*(param.lbox[1])*(param.lbox[2]);
    return nbytes;
}



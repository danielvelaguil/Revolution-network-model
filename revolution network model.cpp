#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <cstddef>


using namespace std;
const int N = 100000; 
int t=1000;
float P=0.5; // probability for creating new connections
float m=0.1; //probability of death of new connections
const int I=int(N*.2);
long seed=1;
int contador[N+I];


int k=10.0;
int n[N];
int s[int(N+I)];
//int s=10;
float matrizp[N+I][100];
//float matriz[12000][120000];
int neig[N+I][100];//neig=np.zeros((N+I,int(N/2)+I))

float matrizp2[N+I][100];
int neig2[N+I][100];//neig=np.zeros((N+I,int(N/2)+I))
int nv2[N+I];

int firings[N+I][1000];//firings=np.zeros((N,t))
int nv[N+I];//nv=np.zeros(N+I)#numero de vecinos
int bandera[N];

int s_update(int s,int n){
 /*int n_bandera=banderas;
 if(int(banderas)==0)
  {*/
  if(s>=1 && s<n-1)
   { s+=1;}
  else if(s==n-1)
    {s=0;
     //n_bandera=0;
     }
  //}
  //else
   //s=0;
  return s;//,n_bandera;
}


int son_vecinos(int h, int h2){
int y;
int bandera2=0;
for(y=0;y<nv2[h];y++)
 if(neig2[h][y]==h2)
 {bandera2=1;
  break;}
 return bandera2;
}


float clustering(){
float conta;
float conta2=0.0;
for(int h=0;h<N;h+=10) // para todos los nodos h
 {
 conta=0.0;
 //conta2=0.0;
 for(int l=0;l<nv2[h];l++)// para todos los vecinos de h
   for(int h2=0;h2<int(nv2[neig2[h][l]]);h2++)// grado del vecino de h
    for(int l2=0;l2<int(nv2[neig2[neig2[h][l]][h2]]);l2++)// para todos los vecinos del vecino de h
     if(int(neig2[neig2[neig2[h][l]][h2]][l2])==int(h)) // ver si el vecino de mi vecino soy yo
     {
      conta+=1.0;
      //if(h==3)
       //cout<<neig2[h][l]<<" "<<neig2[neig2[h][l]][h2]<<" "<<h<<endl;
      }
 if(nv2[h]==0 || nv2[h]==1)
  {
  //cout<<0<<endl;
  conta2+=0;}
 else
  {//cout<<conta/float(nv2[h]*(nv2[h]-1.0))<<endl;
  conta2=conta2+conta/float(nv2[h]*(nv2[h]-1.0));}
 }
 //cout<<conta2/float(N)<<endl;//" "<<float(N)<<" "<<conta/float(nv2[h]*(nv2[h]-1))<<endl;
 return conta2/float(N*.1);//conta2/float(N);
}









int mata_conexion(int h){ /// funcion que elimina todos los vecinos de h
  int g;
  int h2;
  float peso=0.0;
  int vecino=0;
  int y;
  if(nv[h]<nv2[h]) //si tiene conexiones nuevas
   for(g=nv[h];g<nv2[h];g++)
    {
     matrizp2[h][nv2[h]-1]=0.0; //peso del ultimo vecino de h = 0
     y=neig2[h][nv2[h]-1];//ultimo vecino de h = y
     peso=matrizp2[y][nv2[y]-1];  // guardamos el ultimo peso de y
     matrizp2[y][nv2[y]-1]=0; // igualamos a cero el utimo peso de y
     vecino=neig2[y][nv2[y]-1]; // guardamos el ultimo vecino de y
     nv2[h]+=-1; // quito grado a h
     for(h2=nv[y];h2<nv2[y];h2++) //buscamos en los vecinos de y
      {
        if(neig2[y][h2]==h)// si el vecino de y es h
         {
         neig2[y][h2]=vecino; //reemplazamos el vecino
         matrizp2[y][h2]=peso;//reemplazamos el peso
         nv2[y]+=-1; //bajamos el grado
         break; //salimos del for
         }
      }

   }
   return 0;
}


int evolution(float r,float sigma, int ns){
 //
  float xx;
  int ll2=0;
  int i,j,l;
  int h,o;
  if(sigma>-10)
  {
   for(h=0;h<N;h++)
   {
   nv2[h]=nv[h];//esta seria conexion inicial
   for(i=0;i<nv2[h];i++)
    {
     neig2[h][i]=neig[h][i];//solo se actualiza cuando sigma=2
     matrizp2[h][i]=matrizp[h][i];
    }
   }
   }

  for(h=0;h<N;h++)
  {
   contador[h]=0;
   n[h]=ns;
   }
  for(h=0;h<N;h++)
   bandera[h]=0;

  for(h=0;h<N;h++)
   for(o=0;o<t;o++)
    firings[h][o]=0;


 for(h=0;h<N+I;h++)
  s[h]=0;
 float F=0.0;
 //arreglo de vecinos
  float pmax=2*float(sigma)/float(k);  //pmax=2*sigma/k
  int q;
  float rep=1;
  if(sigma<2)
   rep=1;
  int Inh[N];//=np.zeros(N)
  float p[t];//p=[]
  for(h=0;h<t;h++)
   p[h]=0.0;
  float lambd=1.0-exp(-r);//lambd=1-np.exp(-r)
  //cout<<"lambd="<<lambd<<endl;
  //evolucion temporal
  for(i=0;i<t;i++)// i in range(0,t):
  {
   float conteo_b=0.0;
   for(h=0;h<N;h++)
     conteo_b+=float(bandera[h]);
   //cout<<conteo_b<<endl;
   float formula=1.0;//(1.0+1.*float(conteo_b)/float(N));
   //cout<<formula<<endl;
   for(j=0;j<N;j++)// j in range(0,N):
   {
    if(s[j]==0)//:#si las neuronas estan en estado estacionario s=s0 probar:
    {
     if(((double)rand()/((double)RAND_MAX+1.0))<lambd)
      {s[j]=1;
      firings[j][i]=1;}
     for(l=0;l<int(nv2[j]);l++)// actividad por vecinos
      if(s[int(neig2[j][int(l)])]==1)// si el estado de los vecinos es 1 entonces con la prob de la matriz
       if(((double)rand()/((double)RAND_MAX+1.0))<matrizp2[j][l]*pmax*formula)//[int(l)])]+Inh[j])
        {s[j]=1; //se activa el vecinos
        firings[j][i]=1;
         }

     if(s[j]==1)//si las neuronas se activaron checar si va ha ser inhibida (entra en estado estacionario)
      {
      contador[j]=1;
      for(q=0;q<nv2[j];q++)
       if(neig2[j][q]>N) //verificar si tiene una neurona inhibidora conectada
        if(((double)rand()/((double)RAND_MAX+1.0))<matrizp2[j][q]*pmax*rep) //
         {n[j]+=0; // si es asi se incrementa el estado refractario
         //mata_conexion(j);
         bandera[j]=1;//     pmax=2 sigma/k   si sigma=0 and K=10  pmax=0        si sigma=2 and K=10  pmax=0.2
         //cout<<"entro"<<endl;
         }
   //p[i]=sum(firings[:,i])/N)
     }
     else //si las neuronas no se activaron entonces actualizar el valor del contador
      if(contador[j]>0)//el contador es -1 si recientemente se activo
        contador[j]+=-1;
    }
    else //#si no esta en estado estacionario actualizar estado
     {//s[j]=s_update(s[j],n[j]);
      if(contador[j]>0)//el contador es -1 si recientemente se activo
        contador[j]+=-1;
      if(bandera[j]==1)
       {
        s[j]=s_update(s[j],n[j]);
        if(s[j]==0)
         bandera[j]=0;
       }
      else
       s[j]=0;
     //cout<<s[j]<<" "<<bandera[j]<<endl;
     }
    //if(contador[j]>0)//el contador es -1 si recientemente se activo
     // contador[j]+=-1;
    p[i]=p[i]+float(firings[j][i]);///float(N);
     }// fin de for de los nodos
   //for(q=0;q<N;q++)
   int h2;
   int x;
   int y;
   int l2;

   int vecinos_a[100];
   //cout<<i<<","<<i<<endl;
   //cout<<"entro"<<endl; Hebbian rule 2 --triangle
   int va; //numero de vecinos activos
   if(sigma>-10.0)
   {
   for(x=0;x<N;x++) /////////////////// nuevas conexiones /////////////////////////////////////////////
   {
   // if(contador[x]>0)
    if(contador[x]>0) //el vecino de en medio no necesita estar activo
    {
     va=0;
     for(l=0;l<nv2[x];l++)// para todos los vecinos de x, l es el grado de x
      if(contador[neig2[x][l]]>0) //si el vecino de x esta activo
       {
        if(nv2[neig2[x][l]]<90)
         {
         vecinos_a[va]=neig2[x][l];
         va+=1;
         }
       }
    for(y=0;y<round(P*va*(va-1)/2);y++)//probability of creating new connections
     {
     l=0;
     l2=0;
     while(l==l2)
     {
     l=round((double)rand()/((double)RAND_MAX+1.0)*(va-1));
     l2=round((double)rand()/((double)RAND_MAX+1.0)*(va-1));
     }
     if(son_vecinos(vecinos_a[l],vecinos_a[l2])==0 && nv2[vecinos_a[l]]<90 && nv2[vecinos_a[l2]]<90)
      {                          //conexion si previamente no habia
      neig2[vecinos_a[l]][int(nv2[vecinos_a[l]])]=vecinos_a[l2];//establece coneccion entre i y j
      matrizp2[vecinos_a[l]][int(nv2[vecinos_a[l]])]=((double)rand()/((double)RAND_MAX+1.0));//peso entre i y j
      neig2[vecinos_a[l2]][int(nv2[vecinos_a[l2]])]=vecinos_a[l];//establece coneccion entre i y j
      matrizp2[vecinos_a[l2]][int(nv2[vecinos_a[l2]])]=matrizp2[vecinos_a[l]][int(nv2[vecinos_a[l]])];
      nv2[vecinos_a[l]]+=1;
      nv2[vecinos_a[l2]]+=1;
      }
     }
    }
   }
   }

  if(i>499)
  {
  F=float(F)+float(p[i]);///float(t);
  }
  /////////////////    death of new connections   ///////////////////////////

  int sum_=0;
  for(h=0;h<N;h++)
   sum_+=nv2[h]-nv[h];
  int num_con=round(sum_); // total new connections

  float peso=0.0;
  int vecino=0;
  if(num_con>0.1*N)
  for(h=0;h<round(num_con*m*0.5);h++) //probability of removing new connections
  {
   int ban=0;
   while (ban==0)
   {
   x=round((double)rand()/((double)RAND_MAX+1.0)*N);
   if(nv[x]<nv2[x]) //si tiene conexiones nuevas
    ban=1;
   }
    if(ban==ban)
    {
     matrizp2[x][nv2[x]-1]=0.0; //peso del ultimo vecino de x = 0
     y=neig2[x][nv2[x]-1];//ultimo vecino de x = y
     peso=matrizp2[y][nv2[y]-1];  // guardamos el ultimo peso de y
     matrizp2[y][nv2[y]-1]=0; // igualamos a cero el utimo peso de y
     vecino=neig2[y][nv2[y]-1]; // guardamos el ultimo vecino de y
     nv2[x]+=-1; // quito grado a x
     for(h2=nv[y];h2<nv2[y];h2++) //buscamos en los vecinos de y
      {
        if(neig2[y][h2]==x)// si el vecino de y es x
         {
         neig2[y][h2]=vecino; //reemplazamos el vecino
         matrizp2[y][h2]=peso;//reemplazamos el peso
         nv2[y]+=-1; //bajamos el grado
         break; //salimos del for
         }
      }

    }
  } // fin del loop de los nodos
 //////////////////show clustering, degree and activity rho //////////////
  /*
  sum_=0;
  for(h=0;h<N;h++)
   sum_+=nv2[h]-nv[h];
  num_con=round(sum_); // total new connections


  if(ll2==0)
   {
   xx=clustering();
   cout<<i<<","<<xx<<","<<10+float(num_con)/float(N)<<","<<float(p[i])<<endl;
   //cout<<i<<" "<<float(p[i])<<
   }
  else
   {
    cout<<i<<","<<xx<<","<<10+float(num_con)/float(N)<<","<<float(p[i])<<endl;
   }
   ll2=1-ll2;*/
  /////////////////////////////////////////////////////////////////////////////


  } // final time loop 

  /*for(h=0;h<N;h++)
   sum_+=nv2[h];
  cout<<i<<","<<float(sum_)/float(N)-2.0<<","<<float(p[i])<<endl;*/
  return float(F);///(float(t)*float(N));
  //
}



int main(int argc, char ** argv){
  srand(seed);
  int h,o;
  for(h=0;h<N+I;h++)
  {
   for(o=0;o<100;o++)
    matrizp[h][o]=0.0;
   contador[h]=0;
  }
  for(h=0;h<N+I;h++)
   for(o=0;o<100;o++)
    neig[h][o]=0;

 int i,j,l;
 for(h=0;h<N+I;h++)
  nv[h]=0;
 //cout<<"entro"<<endl;
 //conectividad inhibidores
 for(i=0;i<N;i++)//for i in range(0,N):
  for(j=0;j<I;j++)//for j in range(0,I):
   if(((double)rand()/((double)RAND_MAX+1.0))<float(k)/float(N))
    {
    //cout<<"entro"<<endl;
    matrizp[i][int(nv[i])]=((double)rand()/((double)RAND_MAX+1.0));
    matrizp[N+j][nv[N+j]]=matrizp[i][int(N+j)];
    neig[i][int(nv[i])]=N+j;
    neig[N+j][int(nv[N+j])]=i;
    nv[i]+=1;
    nv[N+j]+=1;
    l=0;
   }
 // Erdos-Renyi network
 for(h=0;h<round(N*k/2);h++)
 {
  i=0;
  j=0;
  while(i==j)
  {
   i=round((double)rand()/((double)RAND_MAX+1.0)*N);
   j=round((double)rand()/((double)RAND_MAX+1.0)*N);
  }
    neig[i][int(nv[i])]=j;//establece coneccion entre i y j
    matrizp[i][int(nv[i])]=((double)rand()/((double)RAND_MAX+1.0));//peso entre i y j
    matrizp[j][int(nv[j])]=matrizp[i][int(nv[i])];//peso entre i y j
    neig[j][int(nv[j])]=i;
    nv[i]+=1;
    //origen->destino
    //cout<<0<<","<<j<<","<<i<<endl;
    nv[j]+=1;
 }

  float a=0.0;
  //cout<<a<<endl;
  int array_ns[]={5};
  int deltas_n[]={1,2,5,10};
  float sigmas_[]={.0,.1,.2,.3,.4,.5,.6,.7,.8,.9,.93,.96,1.,1.02,1.05,1.08,1.1,1.13,1.16,1.2,1.25,1.3,1.35,1.4,1.5,1.6,1.7,1.8,1.9,2.0};
  cout<<"r"<<","<<"sigma"<<","<<"F"<<endl;
  float r=0.00001;
  for(i=0;i<30;i++)
    {
    a=evolution(r,float(sigmas_[i]),5);//array_ns[j]);
    cout<<r<<","<<float(sigmas_[i])<<","<<float(a)/float(N*500)<<endl;
   }
  }

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

float closed_interval_rand(float x0, float x1)
{   
        return x0 + (x1 - x0) * rand() / ((float) RAND_MAX);
        
}
float calculusFunction ( float x, float y) /*glavna funkcija,saljemo X i Y kordinatu i dobijamo Z*/
{
    return 3*pow((1-x),2)*exp(-pow(x,2)-pow(y+1,2))-10*(x/5-pow(x,3)-pow(y,5))*exp(-pow(x,2)-pow(y,2))-1/3*exp(-pow(x+1,2)-pow(y,2));
}
int bd (int Gd, int Gg ,int n, float x)
{
 return floor((x-Gd)*(pow(2,n)-1)/(Gg-Gd));
}
float maxfunction( float *z,int populacija) //da seo dredi najveci element u nizu
{
    float max = z[0]; 
    for (int j = 1; j < populacija; j++) 
        if (z[j] > max) 
            max = z[j];
    return max;

}
float minfunction( float *z,int populacija) //da seo dredi najveci element u nizu
{
    float min = z[0]; 
    for (int j = 1; j < populacija; j++) 
        if (z[j] < min) 
            min = z[j];
    return min;

}
float minFitnes(float z,float max)
{
    return max - z;
}
float maxFitnes(float z,float min)
{
    return z-min;
}
float sumFitnes( float *fintes, int populacija)
{
    float sum = 0; 
    for (int i = 0; i < populacija; i++) 
    sum += fintes[i]; 

    return sum;
}
float Pjedinke(float fitnesJedinke, float ocjenaPopulacije)
{
    return fitnesJedinke/ocjenaPopulacije;
}
void iscrtaj (int iteracija,int populacija, float *Xpopulacia, float *Ypopulacia, int *populacia_kodX,int  *populacia_kodY ,float *Zpopulacia, float *fitnesF,float ocjenaPopulacije, float *Pkumulativna)
{
    if(iteracija == 0)
    {
        printf("Pocetna populacija:\n");
        printf("  i          x          y       bdX     bdY            f(x,y)      ff(x,y)      p        q\n");
        printf("----------------------------------------------------------------------------------------------------------------------\n");
    }
    else
    {
        /* code */
    }
     for (int i=0;i<populacija;i++)
    {
        float Pjedinka = Pjedinke(fitnesF[i],ocjenaPopulacije);
        printf("%3d      %7.4f   %7.4f   %5d    %5d     %9.4f       %9.4f   %6.4f   %6.4f\n",i, Xpopulacia[i],Ypopulacia[i], populacia_kodX[i], populacia_kodY[i],Zpopulacia[i],fitnesF[i],Pjedinka,Pkumulativna[i]);
    }
    printf("Ocjena populacije %f", ocjenaPopulacije);

}
void selectionSort(float arr[], int n) ;
void swap(float *xp, float *yp) ;
int findElement (float *Pjedinkii,float Ejedinke, int populacija)
{
    for(int i=0;i<populacija;i++)
    {
        if( Ejedinke== Pjedinkii[i])
        { 
            return i;
        }
    }
    return 0;
}
void ukrstanje( unsigned int*jedinka1, unsigned int *jedinka2, int tackaukstanja)
{
    unsigned int prva,druga;
    prva=*jedinka1;
    druga=*jedinka2;
    unsigned char a1 = (prva >> 24) & 0xFF; // First byte/ high byte
    unsigned char b1 = (prva >> 16) & 0xFF;
    unsigned char c1 = (prva >> 8) & 0xFF;
    unsigned char d1 = (prva) & 0xFF; // Last byte/ smallest byte
    
    unsigned char a2 = (druga >> 24)& 0xFF; // First byte/ high byte
    unsigned char b2 = (druga >> 16)& 0xFF;
    unsigned char c2 = (druga >> 8)& 0xFF;
    unsigned char d2 = (druga)& 0xFF; // Last byte/ smallest byte
    unsigned int result1,result2;
        if(tackaukstanja==0)
         {
             result1 = (a1 << 24) + (b1 << 16) + (c1 << 8) + d2;
             result2 = (a2 << 24) + (b2 << 16) + (c2 << 8) + d1;
         }
         else if (tackaukstanja==1)
         {
             result1 = (a1 << 24) + (b1 << 16) + (c2 << 8) + d2;
             result2 = (a2 << 24) + (b2 << 16) + (c1 << 8) + d1;
         }
         else if (tackaukstanja==2)
         {
             result1 = (a1 << 24) + (b2 << 16) + (c2 << 8) + d2;
             result2 = (a2 << 24) + (b1 << 16) + (c1 << 8) + d1;
         }
         else
         {
             result1 = (a2 << 24) + (b2 << 16) + (c2 << 8) + d2;
             result2 = (a1 << 24) + (b1 << 16) + (c1 << 8) + d1;
         }
        *jedinka1=result1;
        *jedinka2=result2;
        result2 = (a2 << 24) + (b1 << 16) + (c2 << 8) + d1;

}
int main()
{
    int brojac = 0, /*boj uzastopnog ponavljanja iteracija*/
    iteracija = 0, /*broj iteracije*/
    maxBrPonavljanja = 20, /*broj uzastopnog ponavljanja istog resenja, nakon kojeg se prekida program*/
    maxBrIteracija = 1000,  /*maksimalan broj iteracija*/
    minOrmax = 0, /*ako je 0 onda se trazi minimum ako je 1 onda maksimum*/
    Gd = -3, /*donja granica*/
    Gg = 3, /*gornja granica*/
    preciznost = 5; /*decimalna preciznost*/
    const int populacija = 20; /*velicina populacija*/
    float pKombinacija = 0.5, /*vjerovatnoca kombinacije*/
    pMutacija = 0.05;   /*vjerovatnoca mutacije*/ 
    int elite = 2; //koliko direkt prolazi u novu genraciju
    int n = log2f((Gg-Gd)*powf(10,preciznost)+1) + 1 ; /*n broj bita potrebnih za kodovanje*/
    printf("%d\n",n);
   
     srand(time(0)); //Ova funkcija se koristi da svaki put kad pokrenemo kod dobijemo drugacije random brojeve
    float *Xpopulacia, *Ypopulacia; /*populacija X i Y kordinate*/
    Xpopulacia = (float *)calloc(populacija, sizeof(float));
    Ypopulacia = (float *)calloc(populacija, sizeof(float));
    for (int i=0;i<populacija;i++) /*generisanje brojeva iz opsega Gd i Gg za kordinate X i Y*/
    {
        Xpopulacia[i] = closed_interval_rand(Gd,Gg);
        Ypopulacia[i] = closed_interval_rand(Gd,Gg);
    }
    float *Zpopulacia; /*Z populacia*/
    Zpopulacia = (float *)calloc(populacija, sizeof(float));
     for (int i=0;i<populacija;i++)
    {
        Zpopulacia[i] = calculusFunction(Xpopulacia[i],Ypopulacia[i]); /*populacia dobijena od X i Y kordinate*/
    }

    int *populacia_kodX,*populacia_kodY; /*kodovane vrijednosti za X i Y*/
    populacia_kodX = (int *)calloc(populacija, sizeof(int));
    populacia_kodY = (int *)calloc(populacija, sizeof(int));

     for (int i=0;i<populacija;i++) /*Ovde se poziva funkcia  bd za kodovanje*/
    {
       populacia_kodX[i]=bd(Gd,Gg,n,Xpopulacia[i]);
       populacia_kodY[i]=bd(Gd,Gg,n,Ypopulacia[i]);
    }
    float ocjenaPopulacije; //ocjena populacije,suma ff(x)
    float *fitnesF; //fitness funkcija
    fitnesF = (float *)calloc(populacija, sizeof(float));
    float *Pkumultiva;
    Pkumultiva = (float *)calloc(populacija, sizeof(float));
    float *Pjedinkii; //vjetovatnoca jedinki
    Pjedinkii=(float *)calloc(populacija, sizeof(float));
    float *sortPjedinki; //sortirana vjerovatnoca jedinki,kasnije korisna
    sortPjedinki=(float *)calloc(populacija, sizeof(float));
    float *novaPopulacijaX,*novaPopulacijaY; //medjupopulacija X i Y
    novaPopulacijaX=(float *)calloc(populacija, sizeof(float));
    novaPopulacijaY=(float *)calloc(populacija, sizeof(float));



    while (iteracija<100 && brojac<maxBrPonavljanja)
    {
        if(minOrmax == 0)
        {
            float max = maxfunction(Zpopulacia,populacija);  //Da nadjemo najveci element u nizu
             for (int i=0;i<populacija;i++) 
             {
                  fitnesF[i]=minFitnes(Zpopulacia[i],max);
             }
        }
        else
        {
             float min = minfunction(Zpopulacia,populacija);  //Da nadjemo najveci element u nizu
             for (int i=0;i<populacija;i++) 
             {
                  fitnesF[i]=maxFitnes(Zpopulacia[i],min);
             }
        }
        ocjenaPopulacije = sumFitnes(fitnesF,populacija); //dobijamo ocjenu trenutne populacije
        for (int i=0;i<populacija;i++)
        {
            Pjedinkii[i]= Pjedinke(fitnesF[i],ocjenaPopulacije);
            if (i == 0)
            {
                Pkumultiva[i]=Pjedinke(fitnesF[i],ocjenaPopulacije);
            }
            else
            {
                 Pkumultiva[i] =  Pkumultiva[i - 1] + Pjedinke(fitnesF[i],ocjenaPopulacije);
            }
                
        }
        iscrtaj(iteracija,populacija,Xpopulacia,Ypopulacia,populacia_kodX,populacia_kodY,Zpopulacia,fitnesF,ocjenaPopulacije,Pkumultiva);
        printf("Selekcija \n");
        if(elite > 0 && elite < populacija) //elitne jedinke
        {
            for (int i=0;i<populacija;i++) 
             {
                  sortPjedinki[i]=Pjedinkii[i];
             }
             selectionSort(sortPjedinki,populacija);
             for(int i=0;i< elite;i++) //ovde vrsimo dodavanje elitnih jedinke u novu populaciju
             {
                novaPopulacijaX[i]=Xpopulacia[findElement(Pjedinkii,sortPjedinki[i],populacija)];
                novaPopulacijaY[i]=Ypopulacia[findElement(Pjedinkii,sortPjedinki[i],populacija)];
             }
        }
        for(int i=elite; i<populacija;i++) //dodavanje populacije sa dobrom vjerovatnocom u medjupopulaciju    
        {
            float r =closed_interval_rand(0,1);
            int j = 0;
            while (Pkumultiva[j]<r)
            {
                j=j+1;
            
            }
            novaPopulacijaX[i]= Xpopulacia[j];
            novaPopulacijaY[i] =Ypopulacia[j];
        }
        
        
        for(int i=0; i<elite;i++) //dodavanje samo elite
        {
            Xpopulacia[i]=novaPopulacijaX[i];
            Ypopulacia[i] = novaPopulacijaY[i];
        }
        for(int i=elite;i<populacija;) 
        {
            int r1 =floor(closed_interval_rand(0,populacija-1));
            int r2 = floor(closed_interval_rand(0,populacija-1));
            
            float rp =closed_interval_rand(0,1);
            if(rp<pKombinacija)
            {
                int t = (int)floor(closed_interval_rand(0,n)) % 4;
                unsigned int jedinka1X = bd(Gd,Gg,n,novaPopulacijaX[r1]);
                unsigned int jedinka2X = bd(Gd,Gg,n,novaPopulacijaX[r2]);
                unsigned int jedinka1Y = bd(Gd,Gg,n,novaPopulacijaY[r1]);
                unsigned int jedinka2Y = bd(Gd,Gg,n,novaPopulacijaY[r2]);
               
    
                 float jedinka1Xfloat,jedinka2Xfloat,jedinka1Yfloat,jedinka2Yfloat;

               ukrstanje(&jedinka1X,&jedinka2X,t);
                ukrstanje(&jedinka1Y,&jedinka2Y,t);
                 jedinka1Xfloat = Gd+(Gg-Gd)*(int)jedinka1X/(pow(2,n)-1);
                 jedinka2Xfloat = Gd+(Gg-Gd)*(int)jedinka2X/(pow(2,n)-1);
                 jedinka1Yfloat = Gd+(Gg-Gd)*(int)jedinka1Y/(pow(2,n)-1);
                 jedinka2Yfloat = Gd+(Gg-Gd)*(int)jedinka2Y/(pow(2,n)-1);

                novaPopulacijaX[r1]=jedinka1Xfloat;
                novaPopulacijaY[r1]=jedinka1Yfloat;
                novaPopulacijaX[r2]=jedinka2Xfloat;
                novaPopulacijaY[r2]=jedinka2Yfloat;


            }
            Xpopulacia[i]=novaPopulacijaX[r1];
            Ypopulacia[i]= novaPopulacijaY[r1];
            i++;
            if(i<populacija)
            {   
                Xpopulacia[i]=novaPopulacijaX[r2];
                Ypopulacia[i]= novaPopulacijaY[r2];
                i++;
            }

        }
        
         iteracija++; 

    }
    


    return 0;
}
void swap(float *xp, float *yp)  
{  
    float temp = *xp;  
    *xp = *yp;  
    *yp = temp;  
}  
  
void selectionSort(float arr[], int n)  
{  
    int i, j, min_idx;  
  
    // One by one move boundary of unsorted subarray  
    for (i = 0; i < n-1; i++)  
    {  
        // Find the minimum element in unsorted array  
        min_idx = i;  
        for (j = i+1; j < n; j++)  
        if (arr[j] > arr[min_idx])  
            min_idx = j;  
  
        // Swap the found minimum element with the first element  
        swap(&arr[min_idx], &arr[i]);  
    }  
}  
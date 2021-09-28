#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define P 2147483647

typedef struct
{
    long *d;
} fila;

typedef struct
{
    int f; // cantidad de filas
    int c; // cantidad de columnas
    fila *m;
} matriz;

//Funciones del archivo formato lab 1
long SumaP( long a , long b ) ;
long RestaP( long a , long b ) ;
long MultP( long a , long b) ;
long InvP( long a ) ;
//Funciones nuestras
void generar_sistema(int n, int m, matriz *M);
void imprimir_sistema(matriz *M);
void imprimir_solucion(matriz *M, matriz *X);
void triangulacion_superior(matriz *M);
void triangulacion_inferior(matriz *M);
void gaussiana_triangular_superior(matriz *M, matriz *X);
void eliminacion_gaussiana_completa(matriz *M, matriz *X);
void matriz_inversa(matriz *M, matriz *I);
void leer_archivo(matriz *M);

int main()
{
    matriz M, X;
    int i, j = 0;
    double t0, t1, time1;
    //matriz *xd, *xd2;
    //xd = &M;
    //xd2 = &X;
    //t0 = clock();
    //printf("%ld \n", clock());
    srand(time(NULL));
    //long *c;
    t0 = clock();
    /*
    for ( i = 0 ; i < 600 ; i++)
    {
        generar_sistema(300, 301, &M);
        gaussiana_triangular_superior(&M, &X);
    
    t1 = clock();
    time1 = ((t1 - t0) / CLOCKS_PER_SEC);
    printf("TIEMPO TOTAL: %f \n ", time1);
    */
    //generar_sistema(6, 7, &M);
    leer_archivo(&M);
    //printf("--------- \n");
    //gaussiana_triangular_superior(&M, &X);
    //triangulacion_inferior(&M);
    //eliminacion_gaussiana_completa(&M, &X);
    matriz_inversa(&M, &X);
    //printf("--------- \n");
    //printf("--------- \n");
    //t1 = clock();
    //double time = ((t1 - t0) / CLOCKS_PER_SEC);
    //printf("\n %f \n", time);
}


void generar_sistema(int n, int m, matriz *M)
{
    fila F;
    int i, j;
    M->f = n;
    M->c = m;
    M->m = (fila *)malloc(n * sizeof(fila));

    for (i = 0 ; i < n ; i++)
    {
        F.d = (long *)malloc(m * sizeof(long));
        (M->m)[i] = F;
        for (j = 0 ; j < m ; j++)
        {   
            (F.d)[j] = (long)(rand() % P);
        }
    }
    //printf("%ld ", (M->m)[0].d[1]); //EJEMPLO, NO BORRAR: se accede al elemento de la fila 0 columna 1 
}


void imprimir_sistema(matriz *M)
{   
    int i, j;
    for (i = 0 ; i < M->f ; i++)
    {
        printf("[ ");
        for (j = 0 ; j < M->c ; j++)
        {
            printf("%ld ", (M->m)[i].d[j]);
        }
        printf("] \n");
    }
}


void imprimir_solucion(matriz *M, matriz *X)
{
    int i;
    printf("-----------------------------------------------\n");
    printf("MATRIZ TRIANGULADA: \n");
    imprimir_sistema(M);
    printf("-----------------------------------------------\n");
    for (i = 0 ; i < X->f ; i++)
    {
        printf("SOLUCION X%d : %ld \n", i , (X->m)[i].d[0]);
    }
    printf("-----------------------------------------------\n");
}


void triangulacion_superior(matriz *M)
{
    int i, j, k;
    long multiplicador, num_multiplicado;

    for (k = 0 ; k < ((M->f)-1) ; k++) //n pasos (n cantidad de filas)
    {
        for (i = k + 1 ; i < M->f ; i++)
        {
            multiplicador = MultP( (M->m)[i].d[k] , InvP( (M->m)[k].d[k] ) );
            (M->m)[i].d[k] = 0;
            for (j = k + 1 ; j < M->f ; j++)
            {
                num_multiplicado = MultP( multiplicador , (M->m)[k].d[j] );
                (M->m)[i].d[j] = RestaP( (M->m)[i].d[j] , num_multiplicado );
            }
            //la matriz B
            num_multiplicado = MultP( multiplicador, (M->m)[k].d[((M->c)-1)] );
            (M->m)[i].d[((M->c)-1)] = RestaP( (M->m)[i].d[((M->c)-1)] , num_multiplicado );
        }
    }
}


void triangulacion_inferior(matriz *M)
{
    int i, j, k;
    long multiplicador, num_multiplicado;

    for (k = ((M->f)-1) ; k > 0 ; k--)
    {   
        for (i = k - 1 ; i >= 0 ; i--)
        {
            multiplicador = MultP( (M->m)[i].d[k] , InvP( (M->m)[k].d[k] ) );
            (M->m)[i].d[k] = 0;
            for (j = k - 1 ; j >= 0 ; j--)
            {
                num_multiplicado = MultP( multiplicador , (M->m)[k].d[j] );
                (M->m)[i].d[j] = RestaP( (M->m)[i].d[j] , num_multiplicado );
            }
            //la matriz B
            num_multiplicado = MultP( multiplicador, (M->m)[k].d[((M->c)-1)] );
            (M->m)[i].d[((M->c)-1)] = RestaP( (M->m)[i].d[((M->c)-1)] , num_multiplicado );
        }
    }
}


void gaussiana_triangular_superior(matriz *M, matriz *X)
{   
    fila F;
    int i, j, k;
    long multiplicador, num_multiplicado, sumatoria, resta;
    double t0, t1, time;
    X->f = M->f;
    X->c = 1;
    X->m = (fila *)malloc(X->f * sizeof(fila));
    
    if ( M->f <= 6 && M->c <= 7)
    {
        printf("-----------------------------------------------\n");
        printf("MATRIZ INICIAL ( %d x %d ): \n", M->f , M->c);
        imprimir_sistema(M);
    }

    t0 = clock(); //Inicio medición de tiempo

    triangulacion_superior(M); //Triangulamos superiormente la matriz M
    
    //INICIO CREACIÓN VECTOR SOLUCION

    for (i = 0; i < X->f ; i++)
    {   
        F.d = (long *)malloc(X->c * sizeof(long));
        (X->m)[i] = F;
        if (i == ((X->f)-1))    //Ultimo elemento del vector solucion
        {
            (F.d)[0] = MultP( (M->m)[i].d[((M->c)-1)] , InvP( (M->m)[i].d[((M->c)-2)] ) );  
        }
        else
        {
            (F.d)[0] = (M->m)[i].d[((M->c)-1)];
        }
    }

    //FIN CREACIÓN VECTOR SOLUCION
   
    //INICIO SUSTITUCIÓN INVERSA

    for (i = ((X->f)-2) ; i >= 0 ; i--)
    {   
        sumatoria = 0;
        for (j = i + 1 ; j < X->f; j++)
        {
            sumatoria =  SumaP( sumatoria , MultP( (M->m)[i].d[j] , (X->m)[j].d[0] ) );
        }
        resta = RestaP( (M->m)[i].d[((M->c)-1)] , sumatoria );
        (X->m)[i].d[0] = MultP( resta , InvP( (M->m)[i].d[i] ) );
    }
    
    //FIN SUSTITUCIÓN INVERSA

    t1 = clock();   //Fin medición de tiempo

    //PARA IMPRIMIR LAS SOLUCIONES DEL SISTEMA SI LA MATRIZ ES DE HASTA 6X7 (AUMENTADA)
    if ( M->f <= 6 && M->c <= 7)
    {
        imprimir_solucion(M, X);
    }
    
    time = ((t1 - t0) / CLOCKS_PER_SEC);
    printf("FINALIZADO. \n");
    printf("TIEMPO TOTAL EN SEGUNDOS: %f  \n" , time);
    printf("-----------------------------------------------\n");
}


void eliminacion_gaussiana_completa(matriz *M, matriz *X)
{
    fila F;
    int i;
    double t0, t1, time;
    X->f = M->f;
    X->c = 1;
    X->m = (fila *)malloc(X->f * sizeof(fila));

    if ( M->f <= 6 && M->c <= 7)
    {
        printf("-----------------------------------------------\n");
        printf("MATRIZ INICIAL ( %d x %d ): \n", M->f , M->c);
        imprimir_sistema(M);
    }

    t0 = clock(); //inicio medicion de tiempo

    //Proceso de triangulación de la matriz
    triangulacion_superior(M);
    triangulacion_inferior(M);

    //INICIO CREACIÓN VECTOR SOLUCION

    for (i = 0; i < X->f ; i++)
    {   
        F.d = (long *)malloc(X->c * sizeof(long));
        (X->m)[i] = F;
        (F.d)[0] = MultP( (M->m)[i].d[((M->c)-1)]  , InvP ( (M->m)[i].d[i] ) ); //Solucion Xi del sistema
    }

    //FIN CREACIÓN VECTOR SOLUCION

    t1 = clock(); //fin de medicion de tiempo

    if ( M->f <= 6 && M->c <= 7)
    {   
        imprimir_solucion(M, X);
    }

    time = ((t1 - t0) / CLOCKS_PER_SEC);
    printf("FINALIZADO. \n");
    printf("TIEMPO TOTAL EN SEGUNDOS: %f  \n" , time);
    printf("-----------------------------------------------\n");
}
//NO FUNCIONA
void matriz_inversa(matriz *M, matriz *I)
{
    fila F;
    int i, j, k;
    long multiplicador, num_multiplicado;
    I->c = M->c;
    I->f = M->f;
    I->m = (fila *)malloc(I->f * sizeof(fila));
    printf("MATRIZ ORIGINAL \n");
    imprimir_sistema(M);
    //INICIO CREACIÓN MATRIZ IDENTIDAD

    for (i = 0 ; i < I->f ; i++)
    {
        F.d = (long *)malloc(I->c * sizeof(long));
        (I->m)[i] = F;
        for (j = 0 ; j < I->c ; j++)
        {
            if (i == j)
            {
                (F.d)[j] = 1;
            }
            else
            {
                (F.d)[j] = 0;
            }
        }
    }

    printf("MATRIZ IDENTIDAD \n");
    imprimir_sistema(I);

    //FIN CREACIÓN MATRIZ IDENTIDAD
    
    //INICIO TRIANGULACIÓN SUPERIOR 

    for (k = 0 ; k < ((M->f)-1) ; k++) //n pasos (n cantidad de filas)
    {
        for (i = k + 1 ; i < M->f ; i++)
        {
            multiplicador = MultP( (M->m)[i].d[k] , InvP( (M->m)[k].d[k] ) );
            (M->m)[i].d[k] = 0;
            (I->m)[i].d[k] = RestaP( (I->m)[i].d[k] , MultP( multiplicador, (M->m)[k].d[k] ) );
            for (j = k + 1 ; j < M->f ; j++)
            {
                num_multiplicado = MultP( multiplicador , (M->m)[k].d[j] );
                //(M->m)[k].d[j] = MultP( multiplicador , (M->m)[k].d[j] );
                (M->m)[i].d[j] = RestaP( (M->m)[i].d[j] , num_multiplicado );
                num_multiplicado = MultP( multiplicador , (I->m)[k].d[j] );
                //(I->m)[k].d[j] = MultP(multiplicador, (I->m)[k].d[j]); //matriz identidad
                (I->m)[i].d[j] = RestaP( (I->m)[i].d[j] , num_multiplicado ); //matriz identidad
            }
            //(M->m)[k].d[k] = MultP( (M->m)[k].d[k] , multiplicador );
            //(I->m)[k].d[k] = MultP( (I->m)[k].d[k] , multiplicador ); //matriz identidad
        }
    }
    printf("MATRIZ TRIANGULADA SUPERIORMENTE \n");
    imprimir_sistema(M);
    printf("MATRIZ IDENTIDAD NUEVA \n");
    imprimir_sistema(I);
    //FIN TRIANGULACIÓN SUPERIOR

    //INICIO TRIANGULACIÓN INFERIOR
    
    for (k = ((M->f)-1) ; k > 0 ; k--) 
    {   
        for (i = k - 1 ; i >= 0 ; i--)
        {
            multiplicador = MultP( (M->m)[i].d[k] , InvP( (M->m)[k].d[k] ) );
            (M->m)[i].d[k] = 0;
            (I->m)[i].d[k] = RestaP( (I->m)[i].d[k] , MultP( multiplicador, (M->m)[k].d[k] ) );
            for (j = k - 1 ; j >= 0 ; j--)
            {
                num_multiplicado = MultP( multiplicador , (M->m)[k].d[j] );
                //(M->m)[k].d[j] = MultP( multiplicador , (M->m)[k].d[j] );
                (M->m)[i].d[j] = RestaP( (M->m)[i].d[j] , num_multiplicado );
                num_multiplicado = MultP( multiplicador , (I->m)[k].d[j] );
                //(I->m)[k].d[j] = MultP(multiplicador, (I->m)[k].d[j]); //matriz identidad
                (I->m)[i].d[j] = RestaP( (I->m)[i].d[j] , num_multiplicado ); //matriz identidad
            }
            //(M->m)[k].d[k] = MultP( (M->m)[k].d[k] , multiplicador );
            //(I->m)[k].d[k] = MultP( (I->m)[k].d[k] , multiplicador ); //matriz identidad
        }
    }
    
    printf("MATRIZ TRIANGULADA INFERIORMENTE \n");
    imprimir_sistema(M);
    printf("MATRIZ IDENTIDAD NUEVA \n");
    imprimir_sistema(I);

    long inv;
    printf("----------\n");
    imprimir_sistema(M);
    printf("----------\n");
    imprimir_sistema(I);
    for (i = 0 ; i < M->f ; i++)
    {
        for (j = 0 ; j < M->c ; j++)
        {
            if (i == j)
            {
                inv = InvP( (M->m)[i].d[j] );
                (M->m)[i].d[j] = MultP( (M->m)[i].d[j] , inv );
                (I->m)[i].d[j] = MultP( (I->m)[i].d[j] , inv );
            }else
            {
                (I->m)[i].d[j] = MultP( (I->m)[i].d[j] , inv );    
            }
        }
    }
    printf("----------\n");
    imprimir_sistema(M);
    printf("----------\n");
    imprimir_sistema(I);
    //(I->m)[0].d[0] = MultP( (I->m)[0].d[0] , InvP( (M->m)[0].d[0] ) ); // se multiplica el elemento I[0][0] por 1/M[0][0]
    //imprimir_sistema(M);
    //imprimir_sistema(I);
    //FIN TRIANGULACIÓN INFERIOR
/*
    //INICIO CONVERTIR DIAGONAL EN 1

    long inv;

    for (i = 0 ; i < M->f ; i++)
    {
        for (j = i ; j < M->c ; j++)
        {
            inv = InvP((M->m)[i].d[j]);
            (M->m)[i].d[j] = MultP( (M->m)[i].d[j] , inv ); // SE PUEDE CAMBIAR POR 1
            (I->m)[i].d[j] = MultP( (I->m)[i].d[j] , inv ); //SE MULTIPLICA POR 1/Mij
        }
    }


    //FIN CONVERTIR DIAGONAL EN 1
    */
    //imprimir_sistema(I);
}


void leer_archivo(matriz *M)
{
    int n = 0, m = 0;
    long num;
    FILE *archivo;
    char bus;

    fila F;
    int i, j;

    archivo = fopen("matriz2.txt", "r");

    if (archivo == NULL)
    {
        printf("\nError de apertura del archivo. \n\n");
    }
    else
    {
        while ((bus = fgetc(archivo)) != ' ')
        {
            if (bus != ' ')
                n = n * 10 + (long)bus - 48;
        }

        while ((bus = fgetc(archivo)) != '\n')
        {
            //if (bus != '\n')
                m = m * 10 + (long)bus - 48;
        }

        //printf("n= %d \n", n);
        //printf("m= %d \n", m);
    }

    M->f = n;
    M->c = m;
    M->m = (fila *)malloc(n * sizeof(fila));
    for (i = 0; i < n; i++)
    {
        F.d = (long *)malloc(m * sizeof(long));
        (M->m)[i] = F;
        for (j = 0; j < m; j++)
        {
            num=0;
            while ((bus = fgetc(archivo)) != ' ')
            {
                if(bus =='\n'||bus == ' ')
                    break;
                    num = num * 10 + (long)bus - 48;
            }
            (F.d)[j] = num;
            //printf("%ld ", (F.d)[j]);
        }
        //printf("\n");
    }
    fclose(archivo);
}


long SumaP( long a , long b )
{
    long c;
    c = a + b ;
    if ( c < P )
    {
        return(c);
    }
    else
    {
        return(c-P);
    }
}


long RestaP( long a , long b )
{
    if ( a < b )
    {
        return((P+a)-b);
    }
    else
    {
        return(a-b);
    }
}


long MultP( long a , long b)
{
    long c,d;
    c = a * b ;
    d = c >> 31 ;
    c = ( c ^ ( d << 31 ) ) + d ;
    if ( c < P )
    {
        return(c);
    }
    else
    {
        return(c-P);
    }
}


long InvP( long A )
{
    long a,b,s1,s2,r,u;
    a = A ;
    b = P ;
    s1 = 1 ;
    s2 = 0 ;
    while ( ( a % 2 ) == 0 )
    {
        a >>= 1 ;
        if ( ( s1 % 2 ) == 0 )
        {
            s1 >>= 1 ;
        }
        else
        {
            s1 = (s1+P)>>1 ;
        }
    }
    if ( b > a )
    {
        r = b ;
        b = a ;
        a = r ;
        u = s2 ;
        s2 = s1 ;
        s1 = u ;
    }
    while ( ( b != 0 ) && ( b != a ) )
    {
        r = a - b ;
        u = RestaP( s1 , s2 ) ;
        while ( ( r % 2 ) == 0 )
        {
            r >>= 1 ;
            if ( ( u % 2 ) == 0 )
            {
                u >>= 1 ;
            }
            else
            {
                u = (u+P)>>1 ;
            }
        }
        if ( r < b )
        {
            a = b ;
            b = r ;
            s1 = s2 ;
            s2 = u ;
        }
        else
        {
            a = r ;
            s1 = u ;
        }
    }
    return(s1);
}
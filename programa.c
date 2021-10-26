#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define P (long long)3781996138433268199
#define Pa (long long)1761129190
#define Pb (long long)892783079
#define bAJO (long long)2147483647

typedef struct
{
    long long *d;
} fila;

typedef struct
{
    int f; // cantidad de filas
    int c; // cantidad de columnas
    fila *m;
} matriz;

//Funciones del archivo formato lab 1
long long SumaP( long long a , long long b ) ;
long long RestaP( long long a , long long b ) ;
long long MultP( long long a , long long b) ;
long long InvP( long long a ) ;

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
    // FALTA LIMPIAR LA MEMORIA
    matriz M, X;
    int i, j = 0;
    double t0, t1, time1;
    //matriz *xd, *xd2;
    //xd = &M;
    //xd2 = &X;
    //t0 = clock();
    //printf("%ld \n", clock());
    //srand(time(NULL));
    //long long *c;
    //t0 = clock();
    //generar_sistema(2500, 2501, &M);
    //printf("se creó el sistema, se inicia la triangulación");
    //gaussiana_triangular_superior(&M, &X);
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
    //printf("Primo P => %lld\n", P);
    leer_archivo(&M);
    //imprimir_sistema(&M);
    //printf("--------- \n");
    //gaussiana_triangular_superior(&M, &X);
    //triangulacion_inferior(&M);
    //eliminacion_gaussiana_completa(&M, &X);
    //gaussiana_triangular_superior(&M, &X);
    //printf("--------- \n");
    //printf("--------- \n");
    //t1 = clock();
    //double time = ((t1 - t0) / CLOCKS_PER_SEC);
    //printf("\n %f \n", time);
    matriz_inversa(&M, &X);
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
        F.d = (long long *)malloc(m * sizeof(long long));
        (M->m)[i] = F;
        for (j = 0 ; j < m ; j++)
        {   
            (F.d)[j] = (long long)(rand() % P);
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
            printf("%lld ", (M->m)[i].d[j]);
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
        printf("SOLUCION X%d : %lld \n", i , (X->m)[i].d[0]);
    }
}


void triangulacion_superior(matriz *M)
{
    int i, j, k;
    long long multiplicador, num_multiplicado;

    for (k = 0 ; k < ((M->f)-1) ; k++) 
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
    long long multiplicador, num_multiplicado;

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
    long long multiplicador, num_multiplicado, sumatoria, resta;
    double t0, t1, time;
    X->f = M->f;
    X->c = 1;
    X->m = (fila *)malloc(X->f * sizeof(fila));
    
    if ( M->f <= 6 && M->c <= 7)
    {
        printf("-----------------------------------------------\n");
        printf("MATRIZ INICIAL ( %d x %d ): \n", M->f , M->c);
        imprimir_sistema(M);
        printf("-----------------------------------------------\n");
    }

    t0 = clock(); //Inicio medición de tiempo

    triangulacion_superior(M); //Triangulamos superiormente la matriz M
    
    //INICIO CREACIÓN VECTOR SOLUCION
    
    for (i = 0; i < X->f ; i++)
    {   
        F.d = (long long *)malloc(X->c * sizeof(long long));
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
   
    //INICIO SUSTITUCIÓN REVERSA
    
    for (i = ((X->f)-2) ; i > -1 ; i--)
    {   
        sumatoria = 0;
        for (j = i + 1 ; j < X->f; j++)
        {
            sumatoria =  SumaP( sumatoria , MultP( (M->m)[i].d[j] , (X->m)[j].d[0] ) );
        }
        resta = RestaP( (M->m)[i].d[((M->c)-1)] , sumatoria );
        (X->m)[i].d[0] = MultP( resta , InvP( (M->m)[i].d[i] ) );
    }
    
    //FIN SUSTITUCIÓN REVERSA

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
        F.d = (long long *)malloc(X->c * sizeof(long long));
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
    long long multiplicador, num_multiplicado;
    I->c = M->c;
    I->f = M->f;
    I->m = (fila *)malloc(I->f * sizeof(fila));
    printf("MATRIZ ORIGINAL \n");
    imprimir_sistema(M);
    //INICIO CREACIÓN MATRIZ IDENTIDAD

    for (i = 0 ; i < I->f ; i++)
    {
        F.d = (long long *)malloc(I->c * sizeof(long long));
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
    long long prueba;
    for (k = 0 ; k < ((M->f)-1) ; k++) 
    {
        printf("K = %d \n", k);
        for (i = k + 1 ; i < M->f ; i++)
        {
            printf("i = %d \n", i);
            //multiplicador = MultP( (M->m)[i].d[k] , InvP( (M->m)[k].d[k] ) );
            printf("Multiplicador = M[%d][%d] * (1/M[%d][%d]) => %lld * (1/%lld) \n", i, k, k, k, (M->m)[i].d[k], (M->m)[k].d[k]);
            //prueba = (I->m)[i].d[k] * (1/(M->m)[k].d[k]);
            multiplicador = (M->m)[i].d[k] * (1/(M->m)[k].d[k]);
            printf("Resultado Multiplicador = %lld \n", multiplicador);
            (I->m)[i].d[k] = (I->m)[i].d[k] - (multiplicador * (I->m)[k].d[k]);
            //(I->m)[k].d[k] = MultP( (I->m)[k].d[k], multiplicador);
            //(I->m)[i].d[k] = RestaP( (I->m)[i].d[k] , MultP( multiplicador, (M->m)[k].d[k] ) );
            printf("I[%d][%d] = %lld - (%lld * %lld) \n", i, k, (I->m)[i].d[k], multiplicador, (M->m)[k].d[k]);
            (M->m)[i].d[k] = 0;
            printf("M[%d][%d] = 0 \n", i, k);
            for (j = k + 1 ; j < M->f ; j++)
            {
                printf("j = %d \n", j);
                //num_multiplicado = MultP( multiplicador , (M->m)[k].d[j] );
                printf("num_multiplicado = %lld * %lld \n", multiplicador, (M->m)[k].d[j]);
                num_multiplicado = multiplicador * (M->m)[k].d[j];
                //(M->m)[k].d[j] = MultP( multiplicador , (M->m)[k].d[j] );
                //(M->m)[i].d[j] = RestaP( (M->m)[i].d[j] , num_multiplicado );
                printf("M[%d][%d] = %lld - %lld \n", i, j, (M->m)[i].d[j], num_multiplicado);
                (M->m)[i].d[j] = (M->m)[i].d[j] - num_multiplicado;
                //num_multiplicado = MultP( multiplicador , (I->m)[k].d[j] );
                printf("num_multiplicado = %lld * %lld \n", multiplicador, (I->m)[k].d[j]);
                num_multiplicado = multiplicador * (I->m)[k].d[j];
                //(I->m)[k].d[j] = MultP(multiplicador, (I->m)[k].d[j]); //matriz identidad
                //(I->m)[i].d[j] = RestaP( (I->m)[i].d[j] , num_multiplicado ); //matriz identidad
                printf("I[%d][%d] = %lld - %lld \n", i, j, (I->m)[i].d[j], num_multiplicado);
                (I->m)[i].d[j] = (I->m)[i].d[j] - num_multiplicado; //matriz identidad
            }
            //(M->m)[k].d[k] = MultP( (M->m)[k].d[k] , multiplicador );
            //(I->m)[k].d[k] = MultP( (I->m)[k].d[k] , multiplicador ); //matriz identidad
            printf("MATRIZ M \n");
            imprimir_sistema(M);
            printf("MATRIZ M-1 \n");
            imprimir_sistema(I);
        }
    }
    printf("MATRIZ TRIANGULADA SUPERIORMENTE \n");
    //imprimir_sistema(M);
    //printf("MATRIZ IDENTIDAD NUEVA \n");
    //imprimir_sistema(I);
    //FIN TRIANGULACIÓN SUPERIOR

    //INICIO TRIANGULACIÓN INFERIOR
    
    for (k = ((M->f)-1) ; k > 0 ; k--) 
    {   
        printf("K = %d \n", k);
        for (i = k - 1 ; i >= 0 ; i--)
        {
            printf("i = %d \n", i);
            //multiplicador = MultP( (M->m)[i].d[k] , InvP( (M->m)[k].d[k] ) );
            printf("Multiplicador = M[%d][%d] * (1/M[%d][%d]) => %lld * (1/%lld) \n", i, k, k, k, (M->m)[i].d[k], (M->m)[k].d[k]);
            //prueba = (I->m)[i].d[k] * (1/(M->m)[k].d[k]);
            multiplicador = (M->m)[i].d[k] * (1/(M->m)[k].d[k]);
            printf("Resultado Multiplicador = %lld \n", multiplicador);
            (I->m)[i].d[k] = (I->m)[i].d[k] - (multiplicador * (I->m)[k].d[k]);
            //(I->m)[k].d[k] = MultP( (I->m)[k].d[k], multiplicador);
            //(I->m)[i].d[k] = RestaP( (I->m)[i].d[k] , MultP( multiplicador, (M->m)[k].d[k] ) );
            printf("I[%d][%d] = %lld - (%lld * %lld) \n", i, k, (I->m)[i].d[k], multiplicador, (M->m)[k].d[k]);
            (M->m)[i].d[k] = 0;
            printf("M[%d][%d] = 0 \n", i, k);
            for (j = k - 1 ; j > 0 ; j--)
            {
                printf("j = %d \n", j);
                //num_multiplicado = MultP( multiplicador , (M->m)[k].d[j] );
                printf("num_multiplicado = %lld * %lld \n", multiplicador, (M->m)[k].d[j]);
                num_multiplicado = multiplicador * (M->m)[k].d[j];
                //(M->m)[k].d[j] = MultP( multiplicador , (M->m)[k].d[j] );
                //(M->m)[i].d[j] = RestaP( (M->m)[i].d[j] , num_multiplicado );
                printf("M[%d][%d] = %lld - %lld \n", i, j, (M->m)[i].d[j], num_multiplicado);
                (M->m)[i].d[j] = (M->m)[i].d[j] - num_multiplicado;
                //num_multiplicado = MultP( multiplicador , (I->m)[k].d[j] );
                printf("num_multiplicado = %lld * %lld \n", multiplicador, (I->m)[k].d[j]);
                num_multiplicado = multiplicador * (I->m)[k].d[j];
                //(I->m)[k].d[j] = MultP(multiplicador, (I->m)[k].d[j]); //matriz identidad
                //(I->m)[i].d[j] = RestaP( (I->m)[i].d[j] , num_multiplicado ); //matriz identidad
                printf("I[%d][%d] = %lld - %lld \n", i, j, (I->m)[i].d[j], num_multiplicado);
                (I->m)[i].d[j] = (I->m)[i].d[j] - num_multiplicado; //matriz identidad
            }
            //(M->m)[k].d[k] = MultP( (M->m)[k].d[k] , multiplicador );
            //(I->m)[k].d[k] = MultP( (I->m)[k].d[k] , multiplicador ); //matriz identidad
            printf("MATRIZ M \n");
            imprimir_sistema(M);
            printf("MATRIZ M-1 \n");
            imprimir_sistema(I);
        }
    }
    
    printf("MATRIZ TRIANGULADA INFERIORMENTE \n");
    imprimir_sistema(M);
    //printf("MATRIZ IDENTIDAD NUEVA \n");
    //imprimir_sistema(I);

    long long inv;
    //printf("----------\n");
    //imprimir_sistema(M);
    //printf("----------\n");
    //imprimir_sistema(I);
    
    for (i = 0 ; i < M->f ; i++)
    {
        for (j = 0 ; j < M->c ; j++)
        {
            if (i == j)
            {
                inv = 1 / (M->m)[i].d[j];
                (I->m)[i].d[j] = (I->m)[i].d[j] * inv;
                (M->m)[i].d[j] = (M->m)[i].d[j] * inv;
                //inv = InvP( (M->m)[i].d[j] );
                //(I->m)[i].d[j] = MultP( (I->m)[i].d[j] , inv );
                //(M->m)[i].d[j] = MultP( (M->m)[i].d[j] , inv );
            }
            else
            {
                (I->m)[i].d[j] = (I->m)[i].d[j] * inv;    
                //(I->m)[i].d[j] = MultP( (I->m)[i].d[j] , inv );    
            }
        }
    }
    
    printf("MATRIZ ORIGINAL\n");
    imprimir_sistema(M);
    printf("NUEVA MATRIZ IDENTIDAD\n");
    imprimir_sistema(I);
    //(I->m)[0].d[0] = MultP( (I->m)[0].d[0] , InvP( (M->m)[0].d[0] ) ); // se multiplica el elemento I[0][0] por 1/M[0][0]
    //imprimir_sistema(M);
    //imprimir_sistema(I);
    //FIN TRIANGULACIÓN INFERIOR
/*
    //INICIO CONVERTIR DIAGONAL EN 1

    long long inv;

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
    long long num;
    FILE *archivo;
    char bus;

    fila F;
    int i, j;

    archivo = fopen("matriz2021.txt", "r");

    if (archivo == NULL)
    {
        printf("\nError de apertura del archivo. \n\n");
    }
    else
    {
        while ((bus = fgetc(archivo)) != ' ')
        {
            if (bus != ' ')
                n = n * 10 + (long long)bus - 48;
        }

        while ((bus = fgetc(archivo)) != '\n')
        {
            //if (bus != '\n')
            m = m * 10 + (long long)bus - 48;
        }
    }

    M->f = n;
    M->c = m;
    M->m = (fila *)malloc(n * sizeof(fila));

    for (i = 0; i < n; i++)
    {
        F.d = (long long *)malloc(m * sizeof(long long));
        (M->m)[i] = F;
        for (j = 0; j < m; j++)
        {
            num = 0;
            while ((bus = fgetc(archivo)) != ' ')
            {
                if (bus =='\n' || bus == ' ')
                    break;
                    num = num * 10 + (long long)bus - 48;
            }
            (F.d)[j] = num;
        }
    }

    fclose(archivo);
}


long long SumaP( long long a , long long b )
{
    long long c;
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


long long RestaP( long long a , long long b )
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


long long MultP( long long a , long long b)
{
    long long a0,a1,b0,b1,d0,d1,d2,d3 ;
    long long c ;
    
    a0 = a & bAJO ;
    a1 = a>>31 ;
    b0 = b & bAJO ;
    b1 = b>>31 ;
    d0 = a0 * b0 ;
    d3 = a1 * b1 ;
    a1 += a0 ;
    b1 += b0 ;
    d2 = d0 + d3 ;
    d1 = a1 * b1 ;
    d1 -= d2 ;
    d3 += ( d1>>31 ) ;
    d2 = d1 & bAJO ;
    d2 += d0>>31 ;
    d1 = d0 & bAJO ;
    a1 = d3 / Pa ;
    b1 = ( ( d3 % Pa )<<31 ) + d2 - ( a1 * Pb ) ;
    while ( b1 < 0 )
    {
        b1 += P ;
    }
    while ( b1 >= P )
    {
        b1 -= P ;
    }
    a0 = b1 / Pa ;
    b0 = ( ( b1 % Pa )<<31 ) + d1 - ( a0 * Pb ) ;
    while ( b0 < 0 )
    {
        b0 += P ;
    }
    while ( b0 >= P )
    {
        b0 -= P ;
    }
    return(b0);
}


long long InvP( long long A )
{
    long long a,b,s1,s2,r,u;
    a = A ;
    b = P ;
    s1 = 1 ;
    s2 = 0 ;
    if ( A == 0 )
    {
        printf("Error, division entre 0\n");
        return(0);
    }
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
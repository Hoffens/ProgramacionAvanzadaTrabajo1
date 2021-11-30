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
void imprimir_sistema(matriz *M, int n, int m);
void imprimir_solucion(matriz *M, matriz *X);
void triangulacion_superior(matriz *M);
void triangulacion_inferior(matriz *M);
double gaussiana_triangular_superior(matriz *M, matriz *X);
double eliminacion_gaussiana_completa(matriz *M, matriz *X);
double matriz_inversa(matriz *M, matriz *I);
double solucion_matriz_inversa(matriz *M, matriz *X);  // Utiliza la matriz inversa para dar solución al sistema

void leer_archivo(matriz *M);

int main()
{
    // falta limpiar la memoria
    matriz M, X;
    int i, j = 0;
    double t0, t1, time1;
    //matriz *xd, *xd2;
    //xd = &M;
    //xd2 = &X;
    //t0 = clock();
    //printf("%ld \n", clock());
    srand(time(NULL));
    //long long *c;
    //t0 = clock();
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
    //leer_archivo(&M);
    //imprimir_sistema(&M);
    //printf("--------- \n");
    //gaussiana_triangular_superior(&M, &X);
    //triangulacion_inferior(&M);
    generar_sistema(1000, 1001, &M);
    printf("Se genero correctamente el sistema, se inicia la eliminacion gaussiana completa...\n");
    time1 = eliminacion_gaussiana_completa(&M, &X);
    //matriz_inversa(&M, &X);
    //matriz_inversa(&M, &X);
    //time1 = solucion_matriz_inversa(&M, &X);
    printf("TIEMPO TOTAL EN SEGUNDOS: %f  \n" , time1);
    getchar();
    free(M.m);
    free(X.m);
    //printf("--------- \n");
    //printf("--------- \n");
    //t1 = clock();
    //double time = ((t1 - t0) / CLOCKS_PER_SEC);
    //printf("\n %f \n", time);
    return 0;
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


void imprimir_sistema(matriz *M, int n, int m)
{   
    int i, j;
    for (i = 0 ; i < n ; i++)
    {
        printf("[ ");
        for (j = 0 ; j < m ; j++)
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
    imprimir_sistema(M, M->f, M->c);
    printf("-----------------------------------------------\n");
    for (i = 0 ; i < X->f ; i++)
    {
        printf("SOLUCION X%d : %lld \n", i , (X->m)[i].d[0]);
    }
    printf("-----------------------------------------------\n");
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


double gaussiana_triangular_superior(matriz *M, matriz *X)
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
        imprimir_sistema(M, M->f, M->c);
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
    
    return time = ((t1 - t0) / CLOCKS_PER_SEC);
}


double eliminacion_gaussiana_completa(matriz *M, matriz *X)
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
        imprimir_sistema(M, M->f, M->c);
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

    return time = ((t1 - t0) / CLOCKS_PER_SEC);
}


double matriz_inversa(matriz *M, matriz *I) 
{   
    int i, j, k;
    long long temp;
    double t0, t1, time;
    fila F; //aux fila
    I->c = M->c - 1; //sacamos la columna aumentada
    I->f = M->f;
    I->m = (fila *)malloc(I->f * sizeof(fila));

    t0 = clock();   // Inicio medición de tiempo

    //CREACIÓN DE LA MATRIZ M AUMENTADA CON LA IDENTIDAD 
    for (i = 0 ; i < M->f; i++) 
    {
        F.d = (long long *)malloc((I->c * 2) * sizeof(long long));
        (I->m)[i] = F;
        for (j = 0; j < I->c * 2; j++) 
        {
            if (j < I->c) 
            {
                (I->m)[i].d[j] = (M->m)[i].d[j];
            }
            else
            {
                (I->m)[i].d[j] = 0;
            }
        }
    }
     
    // Se agregan los 1's en la diagonal de la matriz aumentada 
    for (i = 0; i < I->f; i++)
    {
        for (j = 0; j < I->f * 2; j++)
        {
            if (j == (i + I->f))    // diagonal de la matriz aumentada
                (I->m)[i].d[j] = 1;
        }
    }


    // Intercambio de filas
    for (i = I->f - 1; i > 0; i--)
    {
        if ((I->m)[i-1].d[0] < (I->m)[i].d[0])
        {
            F.d = (long long *)malloc(I->f * sizeof(long long));
            F = (I->m)[i];  // Guardamos la fila en una variable auxiliar
            (I->m)[i] = (I->m)[i-1];
            (I->m)[i-1] = F;
        }
    }

    // Reemplaza una fila por la suma de si misma y la constante 
    for (i = 0; i < I->f; i++)
    {
        for (j = 0; j < I->f; j++)
        {
            if (j != i)
            {
                temp = MultP( (I->m)[j].d[i], InvP( (I->m)[i].d[i] ) );  
                for (k = 0; k < 2 * I->f; k++) 
                {
                    (I->m)[j].d[k] = RestaP( (I->m)[j].d[k], MultP( (I->m)[i].d[k], temp ) );
                }
            }
        }
    }

    // Multiplica cada fila por un numero distinto de 0 y divide los numeros de la fila por los numeros de la diagonal
    for (i = 0; i < I->f; i++)
    {
        temp = (I->m)[i].d[i];
        for (j = 0; j < 2 * I->f; j++)
        {
            (I->m)[i].d[j] = MultP( (I->m)[i].d[j], InvP( temp ) );
        }
    }

    t1 = clock();   // Fin de medición de tiempo
    
    // Muestra en pantalla solamente si el sistema es pequeño, 
    // la matriz aumentada final. A la izquierda queda la identidad y a la derecha la inversa
    if ( I->f <= 6 && I->c <= 7) 
    {
        printf("Matriz inversa (aumentada): \n");
        imprimir_sistema(I, I->f, I->c * 2); 
    }

    time = ((t1 - t0) / CLOCKS_PER_SEC);
    
    return time;
}


double solucion_matriz_inversa(matriz *M, matriz *I) 
{   
    int i, j, k, p;
    long long resultado;
    double t0, t1, time, timeInv;
    
    // Obtenemos la inversa de M
    timeInv = matriz_inversa(M, I);
    // Ahora I es una matriz aumentada de la forma I|M^-1

    t0 = clock();   // Inicio medición de tiempo

    // Multiplicamos la matriz I (desde la columna I->C) con la matriz B (Ultima columna de la matriz M) para obtener las soluciones
    for (i = 0; i < I->f; i++)
    {
        for (j = 0; j < 1; j++) 
        {   
            p = I->c;   // la matriz inversa empieza en la columna I->C 
            resultado = 0;
            for (k = 0; k < M->f; k++)
            {
                resultado = SumaP( resultado, MultP((I->m)[i].d[p], (M->m)[k].d[M->c - 1]) );
                p++;
            }

            if ( I->f <= 6 && I->c <= 7)    // Muestra las soluciones solamente si el sistema es pequeño
                printf("SOLUCION X%d: %lld\n", i, resultado);
        }
    }
    
    t1 = clock();   // Fin de medición de tiempo

    // Sumatoria del tiempo de la funcion de matriz inversa + tiempo de la multiplicacion para obtener las soluciones
    time = ((t1 - t0) / CLOCKS_PER_SEC) + timeInv;  

    return time;
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
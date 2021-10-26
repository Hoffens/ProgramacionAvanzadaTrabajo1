#include <stdio.h> /* esenciales para la comunicacion con el programa */
#include <stdlib.h> /* generador aleatorio basico, control de memoria, conversiones basicas, control de procesos */
//#include <math.h> /* funciones matematicas basicas como floor, ceil, log2, etc. */
#include <time.h> /* medicion de tiempo */

#define P (long long)3781996138433268199
#define Pa (long long)1761129190
#define Pb (long long)892783079
#define bAJO (long long)2147483647



long long SumaP( long long a , long long b ) ;
long long RestaP( long long a , long long b ) ;
long long MultP( long long a , long long b) ;
long long InvP( long long a ) ;


int main()
{
    int i,j,k,n,m,sel ;
    clock_t tiempo1,tiempo2;
    long long rr;
    long long *L1,*L2;

    printf("Elegir el tipo de operacion a ejecutar:\n");
    printf("0) Repetir el menu\n");
    printf("1) medir la aritmetica mod p\n");
    printf("2) verificar el inverso modular\n");
    printf("9) Salir del programa\n");
    sel = 0 ;
    do
    {
        printf("Seleccion (0 para repetir el menu): ");
        scanf("%d",&sel);
        switch (sel)
        {
            case 0:
                printf("0) Repetir el menu\n");
                printf("1) medir la aritmetica mod p\n");
                printf("2) verificar el inverso modular\n");
                printf("9) Salir del programa\n");
                break ;
            case 1:
                printf("Entrar el numero n de pruebas: ");
                scanf("%d",&n);
                L1 = (long long*) malloc(n*sizeof(long long)) ;
                L2 = (long long*) malloc(n*sizeof(long long)) ;
                for ( i = 0 ; i < n ; i++ )
                {
                    L1[i] = ( ( (long long)(rand() % 13436) )<<48 ) + ( ( ( (long long)rand() ) % 65536 )<<32 ) + ( ( ( (long long)rand() ) % 65536 )<<16 ) + ( ( (long long)rand() ) % 65536 ) ;
                    L2[i] = ( ( (long long)(rand() % 13436) )<<48 ) + ( ( ( (long long)rand() ) % 65536 )<<32 ) + ( ( ( (long long)rand() ) % 65536 )<<16 ) + ( ( (long long)rand() ) % 65536 ) ;
                }
                tiempo1 = clock() ;
                for ( i = 0 ; i < n ; i++ )
                {
                    for ( j = 0 ; j < n ; j++ )
                    {
                        rr = SumaP( L1[i] , L2[j] ) ;
                    }
                }
                tiempo2 = clock() ;
                printf("Tiempo suma (n^2 pruebas): %f\n" , ( (double)tiempo2 - (double)tiempo1 ) / ( (double)CLOCKS_PER_SEC ) ) ;
                tiempo1 = clock() ;
                for ( i = 0 ; i < n ; i++ )
                {
                    for ( j = 0 ; j < n ; j++ )
                    {
                        rr = RestaP( L1[i] , L2[j] ) ;
                    }
                }
                tiempo2 = clock() ;
                printf("Tiempo resta (n^2 pruebas): %f\n" , ( (double)tiempo2 - (double)tiempo1 ) / ( (double)CLOCKS_PER_SEC ) ) ;
                tiempo1 = clock() ;
                for ( i = 0 ; i < n ; i++ )
                {
                    for ( j = 0 ; j < n ; j++ )
                    {
                        rr = MultP( L1[i] , L2[j] ) ;
                    }
                }
                tiempo2 = clock() ;
                printf("Tiempo multiplicacion (n^2 pruebas): %f\n" , ( (double)tiempo2 - (double)tiempo1 ) / ( (double)CLOCKS_PER_SEC ) ) ;
                tiempo1 = clock() ;
                for ( i = 0 ; i < n ; i++ )
                {
                    rr = InvP( L1[i] ) ;
                }
                tiempo2 = clock() ;
                printf("Tiempo de inversos (n pruebas): %f\n" , ( (double)tiempo2 - (double)tiempo1 ) / ( (double)CLOCKS_PER_SEC ) ) ;
                free( L1 ) ;
                free( L2 ) ;
                break ;
            case 2:
                printf("Entrar el numero de pruebas: ");
                scanf("%d",&n);
                L1 = (long long*) malloc(n*sizeof(long long)) ;
                for ( i = 0 ; i < n ; i++ )
                {
                    L1[i] = ( ( (long long)(rand() % 13436) )<<48 ) + ( ( ( (long long)rand() ) % 65536 )<<32 ) + ( ( ( (long long)rand() ) % 65536 )<<16 ) + ( ( (long long)rand() ) % 65536 ) ;
                }
                printf("Verificando los inversos...\n");
                for ( i = 0 ; i < n ; i++ )
                {
                    rr = InvP( L1[i] ) ;
                    if ( MultP( L1[i] , rr ) != 1 )
                    {
                        printf("Error en inverso de %lld:\n",L1[i]);
                    }
                }
                printf("Terminado\n");
                free( L1 ) ;
                break ;
            case 9:
                break ;
            default:
                printf( "Seleccion invalida, reintentar (9 para salir): ");
                break ;
        }

    }
    while ( sel != 9 ) ;

    return(0);
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






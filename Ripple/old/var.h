#define N 9
#define NRK 31
#define Pi 3.141592653589793238462643
#define eps 0.0001
#define Ni 180
#define hmax 0.0004   //must be less than 0.0002
#define hmin 0.00000001
#define TINY 1.0e-10

#define SAFETY 0.9
#define GROW 1.2
#define PGROW -0.2
#define SHRNK 0.15
#define PSHRNK -0.25
#define ERRCON 1.89e-3
#define MAXTRY 50
#define GAM (1.0/2.0)

/* parameters */
//#define A21 2.0
//#define A31 (48.0/25.0)
//#define A32 (6.0/25.0)
//#define C21 -8.0
//#define C31 (372.0/25.0)
//#define C32 (12.0/5.0)
//#define C41 (-112.0/125.0)
//#define C42 (-54.0/125.0)
//#define C43 (-2.0/5.0)
//#define B1 (19.0/9.0)
//#define B2 (1.0/2.0)
//#define B3 (25.0/108.0)
//#define B4 (125.0/108.0)
//#define E1 (17.0/54.0)
//#define E2 (7.0/36.0)
//#define E3 0.0
//#define E4 (125.0/108.0)
//#define C1X (1.0/2.0)
//#define C2X (-3.0/2.0)
//#define C3X (121.0/50.0)
//#define C4X (29.0/250.0)
//#define A2X 1.0
//#define A3X (3.0/5.0)

/* Cash-Karp R-K parameters */
#define A21 0.2
#define A31 (3./40.)
#define A32 (9./40.)
#define A41 0.3
#define A42 -0.9
#define A43 1.2
#define A51 (-11./54.)
#define A52 2.5
#define A53 (-70./27.)
#define A54 (35./27.)
#define A61 (-1631./55296.)
#define A62 (175./512.)
#define A63 (575./13824.)
#define A64 (44275./110592.)
#define A65 (253./4096.)
#define B1 (37./378)
#define B3 (250./621.)
#define B4 (125./594.)
#define B6 (512./1771.)
#define E1 (B1-2825./27648.)
#define E3 (B3-18575./48384.)
#define E4 (B4-13525./55296.)
#define E5 (-277./14336.)
#define E6 (B6-0.25)
#define A2X 0.2
#define A3X 0.3
#define A4X 0.6
#define A5X 1.0
#define A6X 0.875


/* Kaps-Rentrop parameters */
//#define GAM 0.231
//#define A21 2.0
//#define A31 4.52470820736
//#define A32 4.16352878860
//#define C21 -5.70167533877
//#define C31 6.02015272865
//#define C32 0.159750684673
//#define C41 -1.856343618677
//#define C42 -8.50538085819
//#define C43 -2.08407513602
//#define B1 3.95750374663
//#define B2 4.62489238836
//#define B3 0.617477263873
//#define B4 1.282612945268
//#define E1 -2.30215540292
//#define E2 -3.07363448539
//#define E3 0.873280801802
//#define E4 1.282612945268
//#define C1X GAM
//#define C2X -0.396296677520e-01
//#define C3X 0.550778939579
//#define C4X -0.553509845700e-01
//#define A2X 0.462
//#define A3X 0.880208333333
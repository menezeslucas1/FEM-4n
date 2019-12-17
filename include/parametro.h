extern int PRIM_ITER;

void solucaoInicial(double *sol, int *ID, double **COORD, int nnos);

//===================================================================
double funcaoFonte(double x, double y);

//===================================================================
double difusao();

//===================================================================
double reacao();

//===================================================================
double velocx(double x, double y, double t);

//===================================================================
double velocy(double x, double y, double t);

//===================================================================
double parametroCb();

//====================================================================
//kmax gmres
int kmax_gmres();

//====================================================================
//lmax gmres
int lmax_gmres();

//====================================================================
//tolerancia gmres
double tol_gmres();

//====================================================================
//kkmax metodo DV
int kkmax_metodoDV();

//====================================================================
//tolerancia metodoDV
double tol_metodoDV();

//===================================================================
void lerCondFront(double *bv, int *id, int nosx, int nosy);

double gn(double x,double y);

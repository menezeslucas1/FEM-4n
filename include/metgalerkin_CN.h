//==============================================================================
void matrizVetorLocal_CN(int nNosEl, int e, double alpha, double delta_t, double **E, double **EE, double *Fe, int *id, int **ien, double **coordenada, double *bv, int *ND, int **HN, double t);

//==============================================================================

void montaSistemaGlobal_CN(int nNosEl, double ***A, double ***AA, double *F, int *id, int **ien, double **coordenada, double *bv, int neq, int nel, double alpha, double delta_t, int *ND, int **HN, double t);

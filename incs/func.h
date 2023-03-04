void inp_mesh(void);
void inp_gdf(void);
void deriv_lr(void);
void deriv_du(void);
void step(void);
double flow_lr(void);
double flow_du(void);
void new_gdf(void);

void turb_lr(void);
void turb_du(void);


void out(void);
void out_tec(void);
void memo_on(void);

void riemi(double *par,int *alarm);

void deriv_lr_5(void);
void deriv_du_5(void);

double Riem_SC(double sc_l,double sc_r,double u_o);

double rhsuv(int i, int j);

double  deriv_t(int i, int j);

void smoof(void);

void riemi_H2(double* rqp_io, int *alarm);


#define max(a,b)  (((a) > (b)) ? (a) : (b))
#define min(a,b)  (((a) < (b)) ? (a) : (b))

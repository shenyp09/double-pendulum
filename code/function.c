/**************************************************************
*
* 文 件 名： function.c
* 
* 文件描述：
*     该文件包含了SOR迭代方法和共轭梯度法的函数的函数体
*
* 创建人 ： 谌阳平
*
* 创建时间：2012/12/11
* 
* Character Set：UTF-8
*
* 版  本： 1.0  
*
***************************************************************/

#include <stdio.h>
#include <math.h>

#ifndef __FUNCTION__
#define __FUNCTION__ 
#include "function.h"
#endif

//宏定义区
#define tmax 500   //模拟总时长
#define g 9.8       //重力加速度
#define den 1       //文件数据写入密度


//全局变量定义区


double theta1,theta2,dtheta1[4],dtheta2[4],ddtheta1[4],ddtheta2[4],Energy,Energy_bef,dEnergy,Energy0;
double m1=1,m2=1,l1,l2,m3;
double K1,K2,u;
int step=0;
double dt0=0.001;
double dt;
double Di3=1e-10;
int cnt3;



/**************************************************************
*
* 函数名称：E_etest
*
* 函数输入变量：void
*
* 返回值：double
*
* 函数功能：能量误差计算
*
* 版本：1.0
*
* 创建人：谌阳平 
*
* 创建时间：2012/12/11
*
***************************************************************/
double E_etest(void);

/**************************************************************
*
* 函数名称：Q_calc
*
* 函数输入变量：void
*
* 返回值：void
*
* 函数功能：广义阻力Q计算
*
* 版本：1.0
*
* 创建人：谌阳平 
*
* 创建时间：2012/12/18
*
***************************************************************/
void Q_calc(double *Q);

void Initialize(void){
    double cost1,sint1,cost2,sint2,cost2mt1,sint2mt1;

    //初始化系统物理参数
    m1 = 1;
    m2 = 1;
    l1 = 1;
    l2 = 1;
    m3 = 1;
    u=0;
    theta1 = 45;
    theta2 = 0;
    dtheta1[0] = 0;
    dtheta2[0] = 0;


    theta1 = theta1 * M_PI / 180;
    theta2 = theta2 * M_PI / 180;
    dtheta1[0] = dtheta1[0] * M_PI / 180;
    dtheta2[0] = dtheta2[0] * M_PI / 180;
    ddtheta1[0]=0;
    ddtheta2[0]=0;

    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    //计算初始能量
    Energy=1/6.*(-3*g*((m2+2*m3)*l2*cost2+(m1+2*(m2+m3))*l1*cost1)+(m1+3*(m2+m3))*l1*l1*dtheta1[0]*dtheta1[0]+(m2+3*m3)*l2*l2*dtheta2[0]*dtheta2[0]+3*(m2+2*m3)*l1*l2*dtheta1[0]*dtheta2[0]*cost2mt1);
    Energy_bef=Energy;
    Energy0=Energy;
    dEnergy=Energy-Energy_bef;

    dt=dt0;
}

void move(void){
    double cost1,sint1,cost2,sint2,cost2mt1,sint2mt1;
    double A11,A12,A21,A22;
    double ddtheta1_bef,ddtheta2_bef,dtheta1_tmp,dtheta2_tmp;
    double theta1_tmp,theta2_tmp;
    double b1,b2;
    double Q[2];
    int i;

    //平移之前的数据
    for (i = 3; i > 0 ; --i)
    {
        dtheta1[i]=dtheta1[i-1];
        ddtheta1[i]=ddtheta1[i-1];
        dtheta2[i]=dtheta2[i-1];
        ddtheta2[i]=ddtheta2[i-1];
    }
    theta1_tmp=theta1;
    theta2_tmp=theta2;
    dtheta1_tmp=dtheta1[0];
    dtheta2_tmp=dtheta2[0];

    //计算三角函数
    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    //计算系数矩阵和常数矩阵
    A11=(m1+3*m2+3*m3)*l1*l1/3;
    A12=(m2+2*m3)*l1*l2*cost2mt1/2;
    A21=(m2+2*m3)*l1*l2*cost2mt1/2;
    A22=(m2+3*m3)*l2*l2/3;

    Q_calc(Q);

    b1=-l1*(m1/2+m2+m3)*g*sint1+(m2/2+m3)*l1*l2*dtheta2[0]*dtheta2[0]*sint2mt1;
    b2=-l2*(m2/2+m3)*(g*sint2+l1*dtheta1[0]*dtheta1[0]*sint2mt1);

    b1+=Q[0];
    b2+=Q[1];

    //计算角加速度
    ddtheta1_bef=(b1*A22-b2*A12)/(A11*A22-A21*A12);
    ddtheta2_bef=(b2*A11-b1*A21)/(A11*A22-A21*A12);

    //计算预估角度和角速度
    theta1=theta1+dtheta1[0]*dt+ddtheta1_bef*dt*dt/2;
    theta2=theta2+dtheta2[0]*dt+ddtheta2_bef*dt*dt/2;
    dtheta1[0]=dtheta1[0]+ddtheta1_bef*dt;
    dtheta2[0]=dtheta2[0]+ddtheta2_bef*dt;

    //计算三角函数
    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    //计算系数矩阵和常数矩阵
    A11=(m1+3*m2+3*m3)*l1*l1/3;
    A12=(m2+2*m3)*l1*l2*cost2mt1/2;
    A21=(m2+2*m3)*l1*l2*cost2mt1/2;
    A22=(m2+3*m3)*l2*l2/3;
 
    Q_calc(Q);

    b1=-l1*(m1/2+m2+m3)*g*sint1+(m2/2+m3)*l1*l2*dtheta2[0]*dtheta2[0]*sint2mt1;
    b2=-l2*(m2/2+m3)*(g*sint2+l1*dtheta1[0]*dtheta1[0]*sint2mt1);

    b1+=Q[0];
    b2+=Q[1];

    //计算校正后的角加速度
    ddtheta1[0]=(b1*A22-b2*A12)/(A11*A22-A21*A12);
    ddtheta2[0]=(b2*A11-b1*A21)/(A11*A22-A21*A12);

    //计算角速度
    dtheta1[0]=dtheta1_tmp+(ddtheta1[0]+ddtheta1_bef)*dt/2;
    dtheta2[0]=dtheta2_tmp+(ddtheta2[0]+ddtheta2_bef)*dt/2;

    ddtheta1[0]=ddtheta1_bef;
    ddtheta2[0]=ddtheta2_bef;

    //计算能量及与初始能量之差
    Energy=1/6.*(-3*g*((m2+2*m3)*l2*cost2+(m1+2*(m2+m3))*l1*cost1)+(m1+3*(m2+m3))*l1*l1*dtheta1[0]*dtheta1[0]+(m2+3*m3)*l2*l2*dtheta2[0]*dtheta2[0]+3*(m2+2*m3)*l1*l2*dtheta1[0]*dtheta2[0]*cost2mt1);
    dEnergy=Energy-Energy0;
    Energy_bef=Energy;
}

void move_RK(double h){
    double cost1,sint1,cost2,sint2,cost2mt1,sint2mt1;
    double k[5][5];
    double A11,A12,A21,A22;
    double b1,b2;
    double theta1_tmp,theta2_tmp;
    int i;
    double Q[2];

    //平移之前的数据
    for (i = 3; i > 0 ; --i)
    {
        dtheta1[i]=dtheta1[i-1];
        ddtheta1[i]=ddtheta1[i-1];
        dtheta2[i]=dtheta2[i-1];
        ddtheta2[i]=ddtheta2[i-1];
    }
    theta1_tmp=theta1;
    theta2_tmp=theta2;

    //计算三角函数
    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    //计算系数矩阵
    A11=(m1+3*m2+3*m3)*l1*l1/3;
    A12=(m2+2*m3)*l1*l2*cost2mt1/2;
    A21=(m2+2*m3)*l1*l2*cost2mt1/2;
    A22=(m2+3*m3)*l2*l2/3;

    Q_calc(Q);

    b1=-l1*(m1/2+m2+m3)*g*sint1+(m2/2+m3)*l1*l2*dtheta2[0]*dtheta2[0]*sint2mt1;
    b2=-l2*(m2/2+m3)*(g*sint2+l1*dtheta1[0]*dtheta1[0]*sint2mt1);

    b1+=Q[0];
    b2+=Q[1];

    //计算K矩阵
    k[1][1]=(b1*A22-b2*A12)/(A11*A22-A21*A12);
    k[2][1]=(b2*A11-b1*A21)/(A11*A22-A21*A12);
    k[3][1]=dtheta1[0];
    k[4][1]=dtheta2[0];

    //更新角度和角速度
    dtheta1[0]=dtheta1[1]+h*k[1][1]/2;
    dtheta2[0]=dtheta2[1]+h*k[2][1]/2;
    theta1=theta1_tmp+k[3][1]*h/2;
    theta2=theta2_tmp+k[4][1]*h/2;

    //计算三角函数
    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    //计算系数矩阵和常数矩阵
    A11=(m1+3*m2+3*m3)*l1*l1/3;
    A12=(m2+2*m3)*l1*l2*cost2mt1/2;
    A21=(m2+2*m3)*l1*l2*cost2mt1/2;
    A22=(m2+3*m3)*l2*l2/3;

    Q_calc(Q);

    b1=-l1*(m1/2+m2+m3)*g*sint1+(m2/2+m3)*l1*l2*dtheta2[0]*dtheta2[0]*sint2mt1;
    b2=-l2*(m2/2+m3)*(g*sint2+l1*dtheta1[0]*dtheta1[0]*sint2mt1);

    b1+=Q[0];
    b2+=Q[1];


    //计算K矩阵
    k[1][2]=(b1*A22-b2*A12)/(A11*A22-A21*A12);
    k[2][2]=(b2*A11-b1*A21)/(A11*A22-A21*A12);
    k[3][2]=dtheta1[0];
    k[4][2]=dtheta2[0];

    //更新角度和角速度
    dtheta1[0]=dtheta1[1]+h*k[1][2]/2;
    dtheta2[0]=dtheta2[1]+h*k[2][2]/2;
    theta1=theta1_tmp+k[3][2]*h/2;
    theta2=theta2_tmp+k[4][2]*h/2;

    //计算三角函数
    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    //计算系数矩阵和常数矩阵
    A11=(m1+3*m2+3*m3)*l1*l1/3;
    A12=(m2+2*m3)*l1*l2*cost2mt1/2;
    A21=(m2+2*m3)*l1*l2*cost2mt1/2;
    A22=(m2+3*m3)*l2*l2/3;

    Q_calc(Q);

    b1=-l1*(m1/2+m2+m3)*g*sint1+(m2/2+m3)*l1*l2*dtheta2[0]*dtheta2[0]*sint2mt1;
    b2=-l2*(m2/2+m3)*(g*sint2+l1*dtheta1[0]*dtheta1[0]*sint2mt1);

    b1+=Q[0];
    b2+=Q[1];


    //计算K矩阵
    k[1][3]=(b1*A22-b2*A12)/(A11*A22-A21*A12);
    k[2][3]=(b2*A11-b1*A21)/(A11*A22-A21*A12);
    k[3][3]=dtheta1[0];
    k[4][3]=dtheta2[0];

    //更新角度和角速度
    dtheta1[0]=dtheta1[1]+h*k[1][3];
    dtheta2[0]=dtheta2[1]+h*k[2][3];
    theta1=theta1_tmp+k[3][3]*h;
    theta2=theta2_tmp+k[4][3]*h;

    //计算三角函数
    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    //计算系数矩阵和常数矩阵
    A11=(m1+3*m2+3*m3)*l1*l1/3;
    A12=(m2+2*m3)*l1*l2*cost2mt1/2;
    A21=(m2+2*m3)*l1*l2*cost2mt1/2;
    A22=(m2+3*m3)*l2*l2/3;

    Q_calc(Q);

    b1=-l1*(m1/2+m2+m3)*g*sint1+(m2/2+m3)*l1*l2*dtheta2[0]*dtheta2[0]*sint2mt1;
    b2=-l2*(m2/2+m3)*(g*sint2+l1*dtheta1[0]*dtheta1[0]*sint2mt1);

    b1+=Q[0];
    b2+=Q[1];


    //计算K矩阵
    k[1][4]=(b1*A22-b2*A12)/(A11*A22-A21*A12);
    k[2][4]=(b2*A11-b1*A21)/(A11*A22-A21*A12);
    k[3][4]=dtheta1[0];
    k[4][4]=dtheta2[0];

    //计算角度和角速度
    dtheta1[0]=dtheta1[1]+h*(k[1][1]+2*k[1][2]+2*k[1][3]+k[1][4])/6;
    dtheta2[0]=dtheta2[1]+h*(k[2][1]+2*k[2][2]+2*k[2][3]+k[2][4])/6;
    theta1=theta1_tmp+h*(k[3][1]+2*k[3][2]+2*k[3][3]+k[3][4])/6;
    theta2=theta2_tmp+h*(k[4][1]+2*k[4][2]+2*k[4][3]+k[4][4])/6;

    //计算三角函数
    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    //计算能量
    Energy=1/6.*(-3*g*((m2+2*m3)*l2*cost2+(m1+2*(m2+m3))*l1*cost1)+(m1+3*(m2+m3))*l1*l1*dtheta1[0]*dtheta1[0]+(m2+3*m3)*l2*l2*dtheta2[0]*dtheta2[0]+3*(m2+2*m3)*l1*l2*dtheta1[0]*dtheta2[0]*cost2mt1);
    dEnergy=Energy-Energy0;
    Energy_bef=Energy;
}




double move_RKF(double h, int mode){
    double cost1,sint1,cost2,sint2,cost2mt1,sint2mt1;
    double k[5][7];
    double A11,A12,A21,A22;
    double b1,b2;
    double theta1_tmp,theta2_tmp;
    double theta1_t,theta2_t;
    double dt1[4],dt2[4],ddt1[4],ddt2[4];
    double E_t1,E_t2;
    int i;
    double Q[2];
    double E4,E5;

    //保存原始数据
    for (i = 3; i >= 0 ; --i)
    {
        dt1[i]=dtheta1[i];
        ddt1[i]=ddtheta1[i];
        dt2[i]=dtheta2[i];
        ddt2[i]=ddtheta2[i];
    }
    theta1_t=theta1;
    theta2_t=theta2;
    E_t1=Energy;
    E_t2=Energy_bef;

    //平移之前的数据
    for (i = 3; i > 0 ; --i)
    {
        dtheta1[i]=dtheta1[i-1];
        ddtheta1[i]=ddtheta1[i-1];
        dtheta2[i]=dtheta2[i-1];
        ddtheta2[i]=ddtheta2[i-1];
    }
    theta1_tmp=theta1;
    theta2_tmp=theta2;

    //计算三角函数
    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    //计算系数矩阵
    A11=(m1+3*m2+3*m3)*l1*l1/3;
    A12=(m2+2*m3)*l1*l2*cost2mt1/2;
    A21=(m2+2*m3)*l1*l2*cost2mt1/2;
    A22=(m2+3*m3)*l2*l2/3;

    Q_calc(Q);

    b1=-l1*(m1/2+m2+m3)*g*sint1+(m2/2+m3)*l1*l2*dtheta2[0]*dtheta2[0]*sint2mt1;
    b2=-l2*(m2/2+m3)*(g*sint2+l1*dtheta1[0]*dtheta1[0]*sint2mt1);

    b1+=Q[0];
    b2+=Q[1];

    //计算K矩阵
    k[1][1]=(b1*A22-b2*A12)/(A11*A22-A21*A12);
    k[2][1]=(b2*A11-b1*A21)/(A11*A22-A21*A12);
    k[3][1]=dtheta1[0];
    k[4][1]=dtheta2[0];

    //更新角度和角速度
    dtheta1[0]=dtheta1[1]+h*k[1][1]/4;
    dtheta2[0]=dtheta2[1]+h*k[2][1]/4;
    theta1=  theta1_tmp + h*k[3][1]/4;
    theta2=  theta2_tmp + h*k[4][1]/4;

    //计算三角函数
    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    //计算系数矩阵和常数矩阵
    A11=(m1+3*m2+3*m3)*l1*l1/3;
    A12=(m2+2*m3)*l1*l2*cost2mt1/2;
    A21=(m2+2*m3)*l1*l2*cost2mt1/2;
    A22=(m2+3*m3)*l2*l2/3;

    Q_calc(Q);

    b1=-l1*(m1/2+m2+m3)*g*sint1+(m2/2+m3)*l1*l2*dtheta2[0]*dtheta2[0]*sint2mt1;
    b2=-l2*(m2/2+m3)*(g*sint2+l1*dtheta1[0]*dtheta1[0]*sint2mt1);

    b1+=Q[0];
    b2+=Q[1];


    //计算K矩阵
    k[1][2]=(b1*A22-b2*A12)/(A11*A22-A21*A12);
    k[2][2]=(b2*A11-b1*A21)/(A11*A22-A21*A12);
    k[3][2]=dtheta1[0];
    k[4][2]=dtheta2[0];

    //更新角度和角速度
    dtheta1[0]=dtheta1[1]+3./32*h*k[1][1]+9./32*h*k[1][2];
    dtheta2[0]=dtheta2[1]+3./32*h*k[2][1]+9./32*h*k[2][2];
    theta1 = theta1_tmp + 3./32*h*k[3][1]+9./32*h*k[3][2];
    theta2 = theta2_tmp + 3./32*h*k[4][1]+9./32*h*k[4][2];

    //计算三角函数
    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    //计算系数矩阵和常数矩阵
    A11=(m1+3*m2+3*m3)*l1*l1/3;
    A12=(m2+2*m3)*l1*l2*cost2mt1/2;
    A21=(m2+2*m3)*l1*l2*cost2mt1/2;
    A22=(m2+3*m3)*l2*l2/3;

    Q_calc(Q);

    b1=-l1*(m1/2+m2+m3)*g*sint1+(m2/2+m3)*l1*l2*dtheta2[0]*dtheta2[0]*sint2mt1;
    b2=-l2*(m2/2+m3)*(g*sint2+l1*dtheta1[0]*dtheta1[0]*sint2mt1);

    b1+=Q[0];
    b2+=Q[1];


    //计算K矩阵
    k[1][3]=(b1*A22-b2*A12)/(A11*A22-A21*A12);
    k[2][3]=(b2*A11-b1*A21)/(A11*A22-A21*A12);
    k[3][3]=dtheta1[0];
    k[4][3]=dtheta2[0];

    //更新角度和角速度
    dtheta1[0]=dtheta1[1]+1932./2197*h*k[1][1]-7200./2197*h*k[1][2]+7296./2197*h*k[1][3];
    dtheta2[0]=dtheta2[1]+1932./2197*h*k[2][1]-7200./2197*h*k[2][2]+7296./2197*h*k[2][3];
    theta1=  theta1_tmp + 1932./2197*h*k[3][1]-7200./2197*h*k[3][2]+7296./2197*h*k[3][3];
    theta2=  theta2_tmp + 1932./2197*h*k[4][1]-7200./2197*h*k[4][2]+7296./2197*h*k[4][3];

    //计算三角函数
    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    //计算系数矩阵和常数矩阵
    A11=(m1+3*m2+3*m3)*l1*l1/3;
    A12=(m2+2*m3)*l1*l2*cost2mt1/2;
    A21=(m2+2*m3)*l1*l2*cost2mt1/2;
    A22=(m2+3*m3)*l2*l2/3;

    Q_calc(Q);

    b1=-l1*(m1/2+m2+m3)*g*sint1+(m2/2+m3)*l1*l2*dtheta2[0]*dtheta2[0]*sint2mt1;
    b2=-l2*(m2/2+m3)*(g*sint2+l1*dtheta1[0]*dtheta1[0]*sint2mt1);

    b1+=Q[0];
    b2+=Q[1];


    //计算K矩阵
    k[1][4]=(b1*A22-b2*A12)/(A11*A22-A21*A12);
    k[2][4]=(b2*A11-b1*A21)/(A11*A22-A21*A12);
    k[3][4]=dtheta1[0];
    k[4][4]=dtheta2[0];

    dtheta1[0]=dtheta1[1]+439./216*h*k[1][1]-8*h*k[1][2]+3680./513*h*k[1][3]-845./4104*h*k[1][4];
    dtheta2[0]=dtheta2[1]+439./216*h*k[2][1]-8*h*k[2][2]+3680./513*h*k[2][3]-845./4104*h*k[2][4];
    theta1 = theta1_tmp + 439./216*h*k[3][1]-8*h*k[3][2]+3680./513*h*k[3][3]-845./4104*h*k[3][4];
    theta2 = theta2_tmp + 439./216*h*k[4][1]-8*h*k[4][2]+3680./513*h*k[4][3]-845./4104*h*k[4][4];

    //计算三角函数
    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    //计算系数矩阵和常数矩阵
    A11=(m1+3*m2+3*m3)*l1*l1/3;
    A12=(m2+2*m3)*l1*l2*cost2mt1/2;
    A21=(m2+2*m3)*l1*l2*cost2mt1/2;
    A22=(m2+3*m3)*l2*l2/3;

    Q_calc(Q);

    b1=-l1*(m1/2+m2+m3)*g*sint1+(m2/2+m3)*l1*l2*dtheta2[0]*dtheta2[0]*sint2mt1;
    b2=-l2*(m2/2+m3)*(g*sint2+l1*dtheta1[0]*dtheta1[0]*sint2mt1);

    b1+=Q[0];
    b2+=Q[1];


    //计算K矩阵
    k[1][5]=(b1*A22-b2*A12)/(A11*A22-A21*A12);
    k[2][5]=(b2*A11-b1*A21)/(A11*A22-A21*A12);
    k[3][5]=dtheta1[0];
    k[4][5]=dtheta2[0];

    //计算四阶角度和角速度
    dtheta1[0]=dtheta1[1]+h*(25./216*k[1][1]+1408./2565*k[1][3]+2197./4104*k[1][4]-k[1][5]/5);
    dtheta2[0]=dtheta2[1]+h*(25./216*k[2][1]+1408./2565*k[2][3]+2197./4104*k[2][4]-k[2][5]/5);
    theta1=  theta1_tmp + h*(25./216*k[3][1]+1408./2565*k[3][3]+2197./4104*k[3][4]-k[3][5]/5);
    theta2=  theta2_tmp + h*(25./216*k[4][1]+1408./2565*k[4][3]+2197./4104*k[4][4]-k[4][5]/5);
    
    // //计算三角函数
    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    if(mode==0){
        Energy=1/6.*(-3*g*((m2+2*m3)*l2*cost2+(m1+2*(m2+m3))*l1*cost1)+(m1+3*(m2+m3))*l1*l1*dtheta1[0]*dtheta1[0]+(m2+3*m3)*l2*l2*dtheta2[0]*dtheta2[0]+3*(m2+2*m3)*l1*l2*dtheta1[0]*dtheta2[0]*cost2mt1);
        dEnergy=Energy-Energy0;
        return 0;
    }

    E4=1/6.*(-3*g*((m2+2*m3)*l2*cost2+(m1+2*(m2+m3))*l1*cost1)+(m1+3*(m2+m3))*l1*l1*dtheta1[0]*dtheta1[0]+(m2+3*m3)*l2*l2*dtheta2[0]*dtheta2[0]+3*(m2+2*m3)*l1*l2*dtheta1[0]*dtheta2[0]*cost2mt1);


    dtheta1[0]=dtheta1[1]-8./27*h*k[1][1]+2*h*k[1][2]-3544./2565*h*k[1][3]+1859./4104*h*k[1][4]-11./40*h*k[1][5];
    dtheta2[0]=dtheta2[1]-8./27*h*k[2][1]+2*h*k[2][2]-3544./2565*h*k[2][3]+1859./4104*h*k[2][4]-11./40*h*k[2][5];
    theta1=  theta1_tmp - 8./27*h*k[3][1]+2*h*k[3][2]-3544./2565*h*k[3][3]+1859./4104*h*k[3][4]-11./40*h*k[3][5];
    theta2=  theta2_tmp - 8./27*h*k[4][1]+2*h*k[4][2]-3544./2565*h*k[4][3]+1859./4104*h*k[4][4]-11./40*h*k[4][5];

    //计算三角函数
    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;

    //计算系数矩阵和常数矩阵
    A11=(m1+3*m2+3*m3)*l1*l1/3;
    A12=(m2+2*m3)*l1*l2*cost2mt1/2;
    A21=(m2+2*m3)*l1*l2*cost2mt1/2;
    A22=(m2+3*m3)*l2*l2/3;

    Q_calc(Q);

    b1=-l1*(m1/2+m2+m3)*g*sint1+(m2/2+m3)*l1*l2*dtheta2[0]*dtheta2[0]*sint2mt1;
    b2=-l2*(m2/2+m3)*(g*sint2+l1*dtheta1[0]*dtheta1[0]*sint2mt1);

    b1+=Q[0];
    b2+=Q[1];


    //计算K矩阵
    k[1][6]=(b1*A22-b2*A12)/(A11*A22-A21*A12);
    k[2][6]=(b2*A11-b1*A21)/(A11*A22-A21*A12);
    k[3][6]=dtheta1[0];
    k[4][6]=dtheta2[0];

    //计算角度和角速度
    dtheta1[0]=dtheta1[1]+h*(16./135*k[1][1]+6656./12825*k[1][3]+28561./56430*k[1][4]-9./50*k[1][5]+2./55*k[1][6]);
    dtheta2[0]=dtheta2[1]+h*(16./135*k[2][1]+6656./12825*k[2][3]+28561./56430*k[2][4]-9./50*k[2][5]+2./55*k[2][6]);
    theta1 = theta1_tmp + h*(16./135*k[3][1]+6656./12825*k[3][3]+28561./56430*k[3][4]-9./50*k[3][5]+2./55*k[3][6]);
    theta2 = theta2_tmp + h*(16./135*k[4][1]+6656./12825*k[4][3]+28561./56430*k[4][4]-9./50*k[4][5]+2./55*k[4][6]);

    //计算三角函数
    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;


    //计算能量
    E5=1/6.*(-3*g*((m2+2*m3)*l2*cost2+(m1+2*(m2+m3))*l1*cost1)+(m1+3*(m2+m3))*l1*l1*dtheta1[0]*dtheta1[0]+(m2+3*m3)*l2*l2*dtheta2[0]*dtheta2[0]+3*(m2+2*m3)*l1*l2*dtheta1[0]*dtheta2[0]*cost2mt1);

    //恢复初始条件
    for (i = 0; i < 4 ; ++i)
    {
        dtheta1[i]=dt1[i];
        ddtheta1[i]=ddt1[i];
        dtheta2[i]=dt2[i];
        ddtheta2[i]=ddt2[i];
    }
    theta1=theta1_t;
    theta2=theta2_t;
    Energy=E_t1;
    Energy_bef=E_t2;

    return fabs(E4-E5);
}

double t_calc(double h){
    double Dc3,s,hs;

    //更新时间步长
    Dc3=move_RKF(h, 1);

    s=pow(Di3/Dc3, 0.2)*0.8;

    hs=s*h;

    if(hs>3*dt0) 
    {
        hs=3*dt0;
    }else if(hs<0.05*dt0){
        hs=0.05*dt0;
    }

    return hs;
}


void Q_calc(double* Q){
    int sym1=1,sym2=1,sym3=0;
    double x0;
    double Q11,Q12,Q21=0,Q22;
    double cost1,sint1,cost2,sint2,cost2mt1,sint2mt1;

    //计算三角函数
    cost1 = cos(theta1);
    sint1 = sin(theta1);
    cost2 = cos(theta2);
    sint2 = sin(theta2);
    cost2mt1 = cost2 * cost1 + sint2 * sint1;
    sint2mt1 = sint2 * cost1 - sint1 * cost2;
    if (dtheta1[0]>=0)
    {
        sym1=1;
    }else {
        sym1=-1;
    }

    x0=-l1*dtheta1[0]*cost2mt1/dtheta2[0];

    if (x0>=l2||x0<=0)
    {
        if ( (x0+l2)/2*dtheta2[0]+l1*dtheta1[0]*cost2mt1>0 )
        {
            sym2=1;
        }else {
            sym2=-1;
        }
    } else {
        sym2=0;
        if (l1*dtheta1[0]*cost2mt1>0)
        {
            sym3=1;
        } else{
            sym3=-1;
        }
    }

    Q11=-u*dtheta1[0]*dtheta1[0]*l1*l1*l1*l1/4;

    if (sym1==-1)
    {
        Q11=-Q11;
    }

    if (sym2!=0)
    {
        Q12=-u*l1*cost2mt1/dtheta2[0]*(pow(l2*dtheta2[0]+l1*dtheta1[0]*cost2mt1,3)-pow(l1*dtheta1[0]*cost2mt1,3))/3;
        if (sym2==-1)
        {
            Q12=-Q12;
        }
    }else {
        Q12=-u*l1*cost2mt1/dtheta2[0]*(-pow(l2*dtheta2[0]+l1*dtheta1[0]*cost2mt1,3)-pow(l1*dtheta1[0]*cost2mt1,3))/3;
        if (sym3==-1)
        {
            Q12=-Q12;
        }
    }


    Q21=0;

    if (sym2!=0)
    {
        Q22=-u/(dtheta2[0]*dtheta2[0])*((pow(l2*dtheta2[0]+l1*dtheta1[0]*cost2mt1,4)-pow(l1*dtheta1[0]*cost2mt1,4))/4-l1*dtheta1[1]*cost2mt1*(pow(l2*dtheta2[0]+l1*dtheta1[0]*cost2mt1,3)-pow(l1*dtheta1[0]*cost2mt1,3))/3);
        if (sym2==-1)
        {
            Q22=-Q22;
        }
    }else {
        Q22=-u/(dtheta2[0]*dtheta2[0])*((-pow(l2*dtheta2[0]+l1*dtheta1[0]*cost2mt1,4)-pow(l1*dtheta1[0]*cost2mt1,4))/4-l1*dtheta1[1]*cost2mt1*(-pow(l2*dtheta2[0]+l1*dtheta1[0]*cost2mt1,3)-pow(l1*dtheta1[0]*cost2mt1,3))/3);
        if (sym3==-1)
        {
            Q22=-Q22;
        }
    }

    if (dtheta2[0]==0.0)
    {
        Q12=0;
        Q22=0;
    }



    Q[0]=Q11+Q12;
    Q[1]=Q21+Q22;

}
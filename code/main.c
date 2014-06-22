/**************************************************************
*
* 文 件 名： main.c
* 
* 文件描述：
*        计算机模拟物理第三次作业主程序
*
* 创建人 ： 谌阳平
*
* 创建时间：2012/12/2
*
* 修改时间：2012/12/19
*
* 版  本： 1.1   
* 
* Character Set：UTF-8
*
* 编译平台：Mac OS X 10.8.2 with Intel Core i5 and 8GB RAM
*
* 编译器版本：gcc version 4.2.1 (Based on Apple Inc. build 5658) (LLVM build 2336.11.00)
*
***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "function.h"



extern double dt;
extern double dt0;
extern int step;


int main(int argc, char const *argv[])
{
    double t=0; //时间变量
    int mode=0; //计算模式变量
    clock_t start, finish;  //计时变量
    double  duration;   //总花费时间


    //计算模式选取
    printf("Please input the mode of computing(1. Verlet  2. Runge-Kutta  3. Adaptive Stepsize Method):");
    scanf("%d",&mode);

    start = clock();  //计时开始

    //设定初始条件
    Initialize(); 
    
    //在文件init_data.dat中保存模拟参数
    FILE *fp0 = fopen("/Users/Fermi/Desktop/data_plot/origin/data/init_data.dat", "w");
    fprintf(fp0, "%-8d %-12.3e %-12.3e %-12.3e %-12.3e %-12.3e %-12.3e %-12.3e %8d\n",tmax ,dt  ,m1 ,m2 ,m3 ,l1, l2,u,den);
    fclose(fp0);

    //在data.dat中记录模拟数据
    FILE *fp = fopen("/Users/Fermi/Desktop/data_plot/origin/data/data.dat", "w");

    while(t < tmax){
        //写入数据
        fprintf(fp, "%-12.5e %-12.5e %-12.5e %-12.5e %-12.5e %-18.10e %15.8e %15.8e\n",t,theta1*180/M_PI,theta2*180/M_PI,dtheta1[0]*180/M_PI,dtheta2[0]*180/M_PI,Energy,dEnergy,dt);
        
        //根据选取的模式进行计算
        if (mode==1)
        {
            //Verlet算法
            move();
        }else if (mode==2)
        {
            //经典Runge-Kutta算法
            move_RK(dt0);  
        }else{
            //自适应步长的Runge-Kutta算法

            //更新时间步长
            dt=t_calc(dt);
            move_RKF(dt,0);
        }

        t+=dt;
        step++;
    }

    //关闭数据记录文件指针
    fclose(fp);

    //停止计时，计算程序运行时间
    finish = clock(); 
    duration = (double)(finish - start)/CLOCKS_PER_SEC;  
    printf( "Computing cost %.6f seconds\n\n", duration ); 

    //显示平均时间步长
    printf("Average time stepsize is %.6f ms\n",(double)1000*tmax/(step-1) );
    //显示初始化时间步长
    printf("Initial time stepsize is %.6f ms\n", dt0*1000);


    printf("\nProject is Finished\n\n");
	
}
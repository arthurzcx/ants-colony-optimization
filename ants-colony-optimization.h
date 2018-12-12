#ifndef ACO_ANTS_H_
#define ACO_ANTS_H_

#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#include <math.h>

namespace NsAntsColonyOptimization
{

	// Add a small constant to avoid dividing by zero
	#define HEURISTIC(x)     (1.0 / ((double) x + 0.0001))

	#define EPSILON_ACO            1.0E-32

	#define MAX_ANTS       1024    //Maximum number of ants
	#define MAX_NEIGHBOURS 512     // Maximum number of nearest neighbors in the candidate set

	#define INFTY                 LONG_MAX
	#define MAXIMUM_NO_TRIES      100

	#define TRUE  1
	#define FALSE 0

	#define MAX(x,y)        ((x)>=(y)?(x):(y))
	#define MIN(x,y)        ((x)<=(y)?(x):(y))
	
	//Random generator constant
	#define IA 16807
	#define IM 2147483647
	#define AM (1.0/IM)
	#define IQ 127773
	#define IR 2836
	#define MASK 123459876


	enum ACO_Method
	{
		AS = 1,  //as(Ant System)
		EAS = 2,  //eas(Elite Ant System)
		RAS = 3, //ras(rank-based Ant System)
		MMAS = 4,  //mmas(Max-Min Ant System)
		BWAS = 5,  //bwas(Best-Worst Ant System)
		ACS = 6  	//acs(Ant Colony System)
	};

	struct pointACO
	{
		double x;
		double y;
	};

	struct problem
	{
		long int      optimum; 
		long int      n_near;  
		struct pointACO  *nodeptr;
		long int      **distance;
		long int      **nn_list; 
	};

	typedef struct 
	{
	  long int  *tour;
	  char      *visited;
	  long int  tour_length;
	} ant_struct;

	class  global_aco_ants
	{
	public:
		global_aco_ants()
		{
			best_so_far_ant = NULL;
			pheromone = NULL;
			total = NULL;
			prob_of_selection = NULL;
			restart_best_ant = NULL;
			pant = NULL;
			pinstance.distance = NULL;
			pinstance.nodeptr = NULL;
			pinstance.nn_list = NULL;

			iteration = 0;
			n_ants = 0;   
			nn_ants = 0;  

			rho = 0.0;
			alpha = 0.0; 
			beta = 0.0;
			q_0 = 0.0;

			as_flag = 0;  
			eas_flag = 0;
			ras_flag = 0;
			mmas_flag = 0;
			bwas_flag = 0;
			acs_flag = 0;

			elitist_ants = 0;
			ras_ranks = 0;

			trail_max = 0.0;
			trail_min = 0.0;					
			trail_0 = 0.0;
			u_gb = 0;

			optimal = 0;
			branch_fac = 0.0;
			max_tours = 0;
			branching_factor = 0.0;
			lambda = 0.0;
			restart_found_best =0;
			found_branching = 0.0;
			
			restart_time = 0.0;
			restart_iteration = 0;
			start_time=0;
			seed = 1234567;
			max_time = 0.0;
			n_tours = 0;

			ls_flag = 0;
			nn_ls = 0;
			cityNum = 0;

			pinstance.n_near = 0;
			pinstance.optimum = 0;
		}

	public:
		long int iteration;	
		ant_struct *best_so_far_ant;

	private:
		double   **pheromone;
		double   **total;//global pheromone

		double   *prob_of_selection; //selective probability

		long int n_ants;        // ants counter
		long int nn_ants;    //The length of the nearest neighbor list

		double rho;            //信息素蒸发系数
		double alpha;         //信息素重要程度
		double beta;           //启发式因子重要程度
		double q_0;          // 遍历过程中的最佳选择概率

		long int as_flag;     // = 1, 运行Ant System 算法标志
		long int eas_flag;     // = 1, 运行 Elitist Ant System 算法标志
		long int ras_flag;    // = 1, 运行 Rank-Based Ant System 算法标志
		long int mmas_flag;   // = 1, 运行 MAX-MIN Ant System 算法标志
		long int bwas_flag;   // = 1, 运行Best-Worst Ant System 算法标志
		long int acs_flag;    // = 1, 运行 Ant Ccolony System 算法标志

		long int elitist_ants;    //EAS算法的附加参数，定义精英蚂蚁的数量

		long int ras_ranks;     // RBAS算法附加参数

		double   trail_max;      // MMAS算法中的最大信息素值 
		double   trail_min;		// MMAS算法中的最小信息素值 
		long int u_gb;            // 使用目前最优的蚂蚁进行每次 u_gb 迭代更新
								  //同时是MMAS 算法用于调度目前最佳更新的参数

		double   trail_0;        //在 ACS  and BWAS算法中初始化信息素水平

		double   max_time; //最大运行时间
		long int n_tours;//当前旅程号

		long int optimal; //优化长度
		double   branch_fac;//分支比例系数
		long int max_tours;//最大旅程
		double   branching_factor;//分支因子
		double   lambda;//lambda分支因子
		long int restart_found_best;
		double   found_branching;
		double   restart_time;
		long int restart_iteration;
		ant_struct *restart_best_ant;

		clock_t start_time;//计时器起点

		long int seed;	//随机数种子
		long int ls_flag;          //指示是否使用本地搜索
		long int nn_ls;            //在本地搜索中使用的最近邻居列表的最大深度 
		long int cityNum;          //实例中城市的数目

		problem pinstance; //TSP问题实例
		ant_struct* pant; //蚁群


	//函数定义
	public:
		//初始化
		void init_program(int method, int cityNum, pointACO* nodeptr, long int* unClosedNodeptr, long int unClosedNum, bool** sortRelation);
		//终止判断
		long int termination_condition();
		//构建解决方案
		void construct_solutions();
		//本地搜索
		void local_search();
		//更新统计变量
		void update_statistics();
		//更新信息素
		void pheromone_trail_update();
		//更新统计信息及检查算法是否收敛
		void search_control_and_statistics();
		//释放内存
		void release();

	private:
		//===============ACO_ants=====================//
		
		// BEGIN:    信息素操作函数

		//初始化信息素
		void init_pheromone_trails(double initial_trail);
		//信息素蒸发
		void evaporation();
		//模拟素蒸发
		void evaporation_nn_list();
		//全局信息素更新
		void global_update_pheromone(ant_struct *a);
		//根据权重更新信息素
		void global_update_pheromone_weighted(ant_struct *a, long int weight);
		//信息素和启发式因子的综合效应
		void compute_total_information();
		//计算nn_list列表中信息素和启发式因子的综合效应
		void compute_nn_list_total_information();

		// BEGIN:   蚂蚁构建解决方案函数

		//清空禁忌表
		void ant_empty_memory(ant_struct *a);
		//放置蚂蚁
		void place_ant(ant_struct *a, long int phase);
		//选择最好
		void choose_best_next(ant_struct *a, long int phase);
		//领域选择下一个
		void neighbour_choose_best_next(ant_struct *a, long int phase);
		//选择距离最近的
		void choose_closest_next(ant_struct *a, long int phase);
		//领域选择并移动
		void neighbour_choose_and_move_to_next(ant_struct *a, long int phase);

		//BEGIN:  与蚂蚁相关的辅助函数
		
		//找到最好蚂蚁
		long int find_best();
		//招到最差蚂蚁
		long int find_worst();
		//复制蚂蚁解决方案
		void copy_from_to(ant_struct *a1, ant_struct *a2);
		//为蚂蚁分配内存
		void allocate_ants();
		//计算旅程长度并生成最近邻
		long int nn_tour();
		//计算两只蚂蚁之间的距离
		long int distance_between_ants(ant_struct *a1, ant_struct *a2);


		//BEGIN： MAX-MIN Ant System特有的进程

		//模拟MMAS的信息素蒸发
		void mmas_evaporation_nn_list();
		//保持信息素在范围内
		void check_nn_list_pheromone_trail_limits();
		// 保持信息素在限制内
		void check_pheromone_trail_limits();


		//BEGIN:  Ant Colony System特有进程 

		//ACS信息素加强
		void global_acs_pheromone_update(ant_struct *a);
		//更新蚂蚁路径上信息素
		void local_acs_pheromone_update(ant_struct *a, long int phase);

		//BEGIN: Best Worst Ant System 特有进程

		//BWAS中最差蚂蚁更新
		void bwas_worst_ant_update(ant_struct *a1, ant_struct *a2);
		//BWAS最好蚂蚁信息素突变
		void bwas_pheromone_mutation();


		//==============ACO_InOut=================//

		//设置默认参数
		void set_default_parameters(int method);
		//人口统计信息
		void population_statistics();
		//计算平均节点分支因子
		double node_branching(double l);
		//生成信息素矩阵
		void generate_pheromone_matrix();
		//生成全局信息素矩阵
		void generate_total_matrix();

		//==================ACO_acotsp==================//
		
		//试验开始初始化
		void init_try();
		//AS更新
		void as_update();
		//EAS更新
		void eas_update();
		//RAS更新
		void ras_update();
		//MMAS更新
		void mmas_update();
		//BWAS更新
		void bwas_update();
		//ACS更新
		void acs_global_update();


		//=================ACO_ls  本地搜索==================//

		//生成随机数
		long int * generate_random_permutation(long int n);
		//2-opt旅程
		void two_opt_first(long int *tour);
		//2.5-opt旅程
		void two_h_opt_first(long int *tour);
		//3-opt旅程
		void three_opt_first(long int *tour);

		//====================ACO_timer 计时器==============//

		//计数器开始
		void start_timers();
		//读取当前时间
		double elapsed_time();

		//================ACO_TSP TSP问题处理函数========//
		
		//不同的距离计算函数
		long int round_distance(long int i, long int j);
		//不同的距离计算函数
		long int ceil_distance(long int i, long int j);
		//不同的距离计算函数
		long int geo_distance(long int i, long int j);
		//不同的距离计算函数
		inline long int att_distance(long int i, long int j);
		//计算旅程长度
		long int compute_tour_length(long int *t);
		//计算城际距离
		void compute_distances(long int* unClosedNodeptr, long int unClosedNum, bool** sortRelation);
		//计算最近邻列表
		void compute_nn_lists();


		//===============数学计算================/

		//平均值（long int 数组输入）
		inline double mean(long int *values, long int max);
		//平均值（double 数组输入）
		inline double meanr(double *values, long int max);
		//计算标准差
		inline double std_deviation(long int *values, long int i, double mean);
		//排序数组的辅助函数
		inline void swap(long int v[], long int i, long int j);
		//排序数组的递归函数（快速排序）
		inline void sort(long int v[], long int left, long int right);
		//排序数组的辅助函数
		inline void swap2(long int v[], long int v2[], long int i, long int j);
		//排序一个数组的递归函数；第二个数组执行相同的交换序列
		inline void sort2(long int v[], long int v2[], long int left, long int right);
		//产生一个均匀分布在[0,1]内的随机数
		inline double ran01(long *idum);

	};

} //namespace NsAntsColonyOptimization
#endif //ACO_ANTS_H_